library(data.table)
library(xgboost)
library(lubridate)
library(writexl)

# ============================================================
# BLOOD CULTURE REALLOCATION ANALYSIS
# ============================================================
# Estimates the gain in positive blood cultures if the N
# lowest-risk BCs actually sent in 2024 were replaced by the
# N highest-risk emergency admissions that did NOT receive a BC
# (no BC within +-48h of admission).
#
# BC side:  model applied to 2024 BCs → rank ascending by risk
#           → actual outcomes known, so positives lost = exact
# NBC side: predictions from NBC_9 pipeline → rank descending
#           → positives gained = sum(pred_calibrated) [expected value]
#
# Caveat: NBC predictions are from a model trained on BC patients.
# True NBC positivity is unknown and likely lower due to clinician
# selection. Expected gains should be treated as an upper bound.
# ============================================================

out_dir      <- "Output Files/CSVs/NonBC_Admissions/2024"
excel_out <- "Output Files/Excels/NBC"

# ============================================================
# 1. LOAD AND FEATURE-ENGINEER FULL BC COHORT
#    (needed for winsorization bounds matching training)
# ============================================================

bc_data <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))
bc_data[, CollectionDateTime := ymd_hms(CollectionDateTime)]
bc_data[, Outcome := fifelse(CultivabilityCategory %in%
  c("Gram-positive", "Gram-negative", "Polymicrobial",
    "Fungi", "Anaerobes", "Other"), 1L, 0L)]

bc_data[, Sex_num           := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]
bc_data[, SOCIAL_IMDPercentile := fifelse(is.na(SOCIAL_IMDPercentile),
  median(SOCIAL_IMDPercentile, na.rm = TRUE), SOCIAL_IMDPercentile)]
bc_data[, AdmissionType_num := fcase(
  AdmissionMethod == "Emergency", 1L,
  AdmissionMethod == "Elective",  2L,
  default = 0L)]
bc_data[, SourceType_num := fcase(
  CultureSourceCategory == "Peripheral", 1L,
  CultureSourceCategory %in% c("Central Line", "PICC", "Port/Portacath", "Arterial", "Femoral"), 2L,
  default = 0L)]
for (col in c("Source_Urinary", "Source_Respiratory", "Source_Skin_Soft_Tissue", "Source_Other",
              "HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before",
              "Carbapenem_Before", "Any_Prior_Other_Pos", "Prior_Urine_Pos",
              "PositiveCulturesInPrior28Days")) {
  bc_data[is.na(get(col)), (col) := 0L]
}
bc_data[, AnyGramPos_Prior  := fifelse(GramPositiveCulturesInPrior28Days > 0, 1, 0)]
bc_data[, AnyGramNeg_Prior  := fifelse(GramNegativeCulturesInPrior28Days > 0, 1, 0)]
bc_data[, VS_AnyFever       := fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), 1, 0)]
bc_data[, VS_AnyHypotension := fifelse(VS_AnyHypotension == 1 & !is.na(VS_AnyHypotension), 1, 0)]
bc_data[is.na(LOS_BeforeCulture_days), LOS_BeforeCulture_days := 0]

cont_vars <- c("VS_Baseline_HR", "VS_Baseline_Temp", "VS_Baseline_SBP", "VS_Baseline_RR",
               "VS_Baseline_NEWS", "BT_Closest_CRP", "BT_Closest_WCC", "BT_Closest_Neut",
               "BT_Closest_Lymph", "BT_Closest_Alb", "BT_Closest_Plt", "BT_Closest_Creat",
               "BT_Closest_Bili", "BT_Closest_Urea", "BT_Closest_Hb", "BT_Closest_Na",
               "BT_Closest_Eosin", "BT_Closest_Mono", "BT_Closest_ALP", "BT_Closest_HCT",
               "BG_Closest_Lactate", "BG_Closest_pH", "BG_Closest_pCO2", "LOS_BeforeCulture_days")
for (col in cont_vars) {
  p01 <- quantile(bc_data[[col]], 0.01, na.rm = TRUE)
  p99 <- quantile(bc_data[[col]], 0.99, na.rm = TRUE)
  bc_data[get(col) < p01, (col) := p01]
  bc_data[get(col) > p99, (col) := p99]
}

# ============================================================
# 2. FILTER TO 2024 AND APPLY MODEL
# ============================================================

bc_2024 <- bc_data[!is.na(CollectionDateTime) & year(CollectionDateTime) == 2024]
rm(bc_data)

feature_cols <- c(
  "AgeAtCulture", "Sex_num",
  "AdmissionType_num", "LOS_BeforeCulture_days", "SourceType_num",
  "PositiveCulturesInPrior28Days",
  "Source_Respiratory",
  "HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before", "Carbapenem_Before",
  "Any_Prior_Other_Pos", "Prior_Urine_Pos",
  "VS_Baseline_Temp", "VS_Baseline_SBP",
  "VS_Baseline_NEWS", "VS_AnyFever", "VS_Baseline_HR",
  "BT_Closest_CRP", "BT_Closest_WCC", "BT_Closest_Neut", "BT_Closest_Lymph",
  "BT_Closest_Alb", "BT_Closest_Plt", "BT_Closest_Creat", "BT_Closest_Bili",
  "BT_Closest_Na",
  "BT_Closest_Eosin", "BT_Closest_Mono", "BT_Closest_ALP",
  "BG_Closest_Lactate", "BG_Closest_pH"
)
for (col in feature_cols) if (!col %in% names(bc_2024)) bc_2024[, (col) := NA_real_]

model       <- xgb.load("Models/XGBoost/FullCohort/xgb_fullcohort.json")
platt_coefs <- readRDS("Models/XGBoost/FullCohort/xgb_fullcohort_platt_calibration.rds")

pred_raw <- predict(model, xgb.DMatrix(as.matrix(bc_2024[, ..feature_cols])))
pred_logit <- qlogis(pmax(pmin(pred_raw, 0.999), 0.001))
pred_cal   <- plogis(platt_coefs[1] + platt_coefs[2] * pred_logit)
bc_2024[, pred_calibrated := pred_cal]

# ============================================================
# 3. LOAD NBC 2024 PREDICTIONS
# ============================================================

nbc <- fread(file.path(out_dir, "nonbc_admissions_with_predictions.csv"),
             colClasses = c(ClusterID = "character"))

bc_sorted  <- bc_2024[order(pred_calibrated)]     # lowest risk first
nbc_sorted <- nbc[order(-pred_calibrated)]         # highest risk first

n_bc      <- nrow(bc_2024)
n_nbc     <- nrow(nbc_sorted)
base_pos  <- sum(bc_2024$Outcome)
base_rate <- mean(bc_2024$Outcome)

# ============================================================
# 4. REALLOCATION SWEEP
# ============================================================

swap_sizes <- c(500, 1000, 2000, 3000, 5000, 7500, 10000, 15000, 20000, n_bc)
n_comparable <- min(n_bc, n_nbc)
crossover    <- which(nbc_sorted$pred_calibrated[1:n_comparable] <
                      bc_sorted$pred_calibrated[1:n_comparable])
max_k_xgb   <- if (length(crossover) > 0) crossover[1] - 1L else n_comparable
swap_sizes   <- sort(unique(c(swap_sizes[swap_sizes < max_k_xgb], max_k_xgb)))

sweep <- rbindlist(lapply(swap_sizes, function(k) {
  culled <- bc_sorted[1:k]
  added  <- nbc_sorted[1:k]
  lost   <- sum(culled$Outcome)
  gained <- sum(added$pred_calibrated)
  new_total <- base_pos - lost + gained
  data.table(
    Swaps                  = k,
    BC_Culled_MeanRisk_pct = round(100 * mean(culled$pred_calibrated), 2),
    BC_Culled_MaxRisk_pct  = round(100 * max(culled$pred_calibrated),  2),
    Positives_Lost_Actual  = lost,
    NBC_Added_MeanRisk_pct = round(100 * mean(added$pred_calibrated),  2),
    Positives_Gained_Exp   = round(gained),
    Net_Change             = round(gained - lost),
    New_Total_Positives    = round(new_total),
    Abs_Increase_pct_pts   = round(100 * new_total / n_bc - 100 * base_rate, 2),
    Rel_Increase_pct       = round(100 * (new_total - base_pos) / base_pos, 1)
  )
}))

# ============================================================
# 5. DETAILED SNAPSHOT AT 10,000 SWAPS
# ============================================================

k         <- 10000
culled_10k <- bc_sorted[1:k]
added_10k  <- nbc_sorted[1:k]

detail <- data.table(
  Metric = c(
    "Year",
    "Total 2024 BCs sent", "Total 2024 positives", "Baseline positivity rate",
    "NBC pool size (no BC within 48h)",
    "--- CULLED (10k lowest-risk BCs) ---",
    "Mean predicted risk", "Max predicted risk",
    "Actual positives among culled", "Actual positivity among culled",
    "--- ADDED (10k highest-risk NBCs) ---",
    "Mean predicted risk", "Max predicted risk",
    "Expected positives (sum of pred)",
    "--- NET RESULT ---",
    "Positives lost", "Positives gained (expected)",
    "Net change in positives",
    "New total positives (expected)",
    "New positivity rate (expected)",
    "Absolute increase (percentage points)",
    "Relative increase"
  ),
  Value = c(
    "2024",
    format(n_bc, big.mark = ","),
    format(base_pos, big.mark = ","),
    sprintf("%.1f%%", 100 * base_rate),
    format(n_nbc, big.mark = ","),
    "",
    sprintf("%.2f%%", 100 * mean(culled_10k$pred_calibrated)),
    sprintf("%.2f%%", 100 * max(culled_10k$pred_calibrated)),
    as.character(sum(culled_10k$Outcome)),
    sprintf("%.2f%%", 100 * mean(culled_10k$Outcome)),
    "",
    sprintf("%.2f%%", 100 * mean(added_10k$pred_calibrated)),
    sprintf("%.2f%%", 100 * max(added_10k$pred_calibrated)),
    sprintf("%.0f", sum(added_10k$pred_calibrated)),
    "",
    as.character(sum(culled_10k$Outcome)),
    sprintf("%.0f", sum(added_10k$pred_calibrated)),
    sprintf("+%.0f", sum(added_10k$pred_calibrated) - sum(culled_10k$Outcome)),
    sprintf("%.0f", base_pos - sum(culled_10k$Outcome) + sum(added_10k$pred_calibrated)),
    sprintf("%.1f%%", 100 * (base_pos - sum(culled_10k$Outcome) + sum(added_10k$pred_calibrated)) / n_bc),
    sprintf("+%.1f pp", 100 * (base_pos - sum(culled_10k$Outcome) + sum(added_10k$pred_calibrated)) / n_bc - 100 * base_rate),
    sprintf("+%.1f%%", 100 * (sum(added_10k$pred_calibrated) - sum(culled_10k$Outcome)) / base_pos)
  )
)

fwrite(sweep,  file.path(out_dir, "reallocation_sweep.csv"))
fwrite(detail, file.path(out_dir, "reallocation_10k_detail.csv"))

write_xlsx(list(
  Sweep_All_Thresholds = as.data.frame(sweep),
  Detail_10k_Swap      = as.data.frame(detail)
), file.path(excel_out, "NBC_Reallocation_Analysis.xlsx"))
