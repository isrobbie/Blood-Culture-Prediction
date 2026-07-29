################################################################################
# ED-specific deduplication sensitivity analysis
#
# Fixed production model (xgb_fullcohort.json + Platt calibration, no retraining)
# evaluated on the ED test subset (CultureLocation == "ED", >= 2024-01-01) under
# four cohort definitions:
#   1. Original    — all ED cultures (matches the 0.876 ED subgroup figure)
#   2. Same-day    — one culture per (patient, day); any positive = positive
#   3. Per-episode — first culture per 14-day episode per patient
#   4. Per-patient — first-ever culture per patient
#
# Reports AUROC (DeLong 95% CI), AUPRC, n, positives, positivity for each.
################################################################################

library(data.table)
library(xgboost)
library(pROC)
library(PRROC)

set.seed(42)

model_dir <- "Models/XGBoost/FullCohort"
csv_dir   <- "Output Files/CSVs/XGBoost/FullCohort"

# ---------------------------------------------------------------------------
# Load + feature-engineer (identical to Do_XGBoost_FullCohort.r)
# ---------------------------------------------------------------------------
bc_data <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))

bc_data <- bc_data[CultivabilityCategory %in% c(
  "Gram-positive", "Gram-negative", "Polymicrobial", "Fungi", "Anaerobes", "Other",
  "Culture-negative", "Contaminant - single", "Contaminant - recurrent",
  "Polymicrobial - contaminant")]

bc_data[, Outcome := fifelse(CultivabilityCategory %in% c(
  "Gram-positive", "Gram-negative", "Polymicrobial", "Fungi", "Anaerobes", "Other"), 1L, 0L)]

bc_data[, Sex_num := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]
median_imd <- median(bc_data$SOCIAL_IMDPercentile, na.rm = TRUE)
bc_data[is.na(SOCIAL_IMDPercentile), SOCIAL_IMDPercentile := median_imd]

bc_data[, AdmissionType_num := fcase(
  AdmissionMethod == "Emergency", 1,
  AdmissionMethod == "Elective", 2,
  default = 0)]

bc_data[, SourceType_num := fcase(
  CultureSourceCategory == "Peripheral", 1,
  CultureSourceCategory %in% c("Central Line", "PICC", "Port/Portacath", "Arterial", "Femoral"), 2,
  default = 0)]

for (col in c("Source_Urinary", "Source_Respiratory", "Source_Skin_Soft_Tissue", "Source_Other",
              "HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before", "Carbapenem_Before",
              "Any_Prior_Other_Pos", "Prior_Urine_Pos")) {
  bc_data[is.na(get(col)), (col) := 0]
}
bc_data[is.na(PositiveCulturesInPrior28Days), PositiveCulturesInPrior28Days := 0]

bc_data[, VS_AnyFever := fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), 1, 0)]
bc_data[, VS_AnyHypotension := fifelse(VS_AnyHypotension == 1 & !is.na(VS_AnyHypotension), 1, 0)]
bc_data[is.na(LOS_BeforeCulture_days), LOS_BeforeCulture_days := 0]

bc_data[, CollectionDate := as.Date(CollectionDateTime)]

feature_cols <- c(
  "AgeAtCulture", "Sex_num", "AdmissionType_num", "LOS_BeforeCulture_days",
  "SourceType_num", "PositiveCulturesInPrior28Days", "Source_Respiratory",
  "HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before", "Carbapenem_Before",
  "Any_Prior_Other_Pos", "Prior_Urine_Pos",
  "VS_Baseline_Temp", "VS_Baseline_SBP", "VS_Baseline_NEWS", "VS_AnyFever", "VS_Baseline_HR",
  "BT_Closest_CRP", "BT_Closest_WCC", "BT_Closest_Neut", "BT_Closest_Lymph",
  "BT_Closest_Alb", "BT_Closest_Plt", "BT_Closest_Creat", "BT_Closest_Bili", "BT_Closest_Na",
  "BT_Closest_Eosin", "BT_Closest_Mono", "BT_Closest_ALP",
  "BG_Closest_Lactate", "BG_Closest_pH")

# ---------------------------------------------------------------------------
# ED test subset (fixed model, no retraining)
# ---------------------------------------------------------------------------
split_date <- as.Date("2024-01-01")
ed_test    <- bc_data[CollectionDate >= split_date & CultureLocation == "ED"]

model      <- xgb.load(file.path(model_dir, "xgb_fullcohort.json"))
platt_coef <- readRDS(file.path(model_dir, "xgb_fullcohort_platt_calibration.rds"))

predict_cal <- function(dt) {
  raw   <- predict(model, xgb.DMatrix(data = as.matrix(dt[, ..feature_cols])))
  logit <- qlogis(pmax(pmin(raw, 0.999), 0.001))
  plogis(platt_coef[1] + platt_coef[2] * logit)
}

# ---------------------------------------------------------------------------
# Evaluate a cohort definition
# ---------------------------------------------------------------------------
evaluate <- function(dt, label) {
  p  <- predict_cal(dt)
  r  <- roc(dt$Outcome, p, quiet = TRUE)
  ci <- ci.auc(r, method = "delong")
  ap <- pr.curve(scores.class0 = p[dt$Outcome == 1],
                 scores.class1 = p[dt$Outcome == 0])$auc.integral
  data.table(
    Approach   = label,
    N          = nrow(dt),
    Positives  = sum(dt$Outcome),
    Pos_Rate   = round(100 * mean(dt$Outcome), 1),
    AUROC      = round(as.numeric(ci[2]), 3),
    AUROC_low  = round(as.numeric(ci[1]), 3),
    AUROC_high = round(as.numeric(ci[3]), 3),
    AUPRC      = round(ap, 3))
}

# 1. Original — all ED cultures
res <- list(evaluate(ed_test, "Original (all ED)"))

# 2. Same-day — one per (patient, day); any positive = positive
setorder(ed_test, ClusterID, CollectionDate, CollectionDateTime)
ed_sd <- ed_test[, .SD[1L], by = .(ClusterID, CollectionDate)]
ed_sd[, Outcome := NULL]
ed_sd <- merge(ed_sd,
               ed_test[, .(Outcome = max(Outcome)), by = .(ClusterID, CollectionDate)],
               by = c("ClusterID", "CollectionDate"))
res[[2]] <- evaluate(ed_sd, "Same-day dedup")

# 3. Per-episode — first culture per 14-day episode per patient
setorder(ed_test, ClusterID, CollectionDateTime)
ed_ep <- copy(ed_test)
ed_ep[, PrevDate := shift(CollectionDate), by = ClusterID]
ed_ep[, Gap := as.numeric(CollectionDate - PrevDate)]
ed_ep[, NewEp := fifelse(is.na(Gap) | Gap > 14, 1L, 0L)]
ed_ep[NewEp == 1, EpNo := seq_len(.N), by = ClusterID]
ed_ep[, EpNo := nafill(EpNo, type = "locf"), by = ClusterID]
ed_ep <- ed_ep[, .SD[1L], by = .(ClusterID, EpNo)]
res[[3]] <- evaluate(ed_ep, "Per-episode")

# 4. Per-patient — first-ever culture per patient
setorder(ed_test, ClusterID, CollectionDateTime)
ed_pat <- ed_test[, .SD[1L], by = ClusterID]
res[[4]] <- evaluate(ed_pat, "Per-patient")

summary_dt <- rbindlist(res)
fwrite(summary_dt, file.path(csv_dir, "ed_cohort_approach_comparison.csv"))
