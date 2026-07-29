library(data.table)
library(xgboost)
library(lubridate)
library(writexl)
library(ggplot2)

# ============================================================
# APPLY BLOOD CULTURE POSITIVITY MODEL TO NON-BC ADMISSIONS
# ============================================================
# Sections 9-11 from Do_Apply_Model_NonBC_Admissions.r.
# Loads the prepared NBC cohort, applies feature engineering
# matching the XGBoost training pipeline, generates predictions,
# and characterises the high-risk cohort (pred > 5%).
# ============================================================

nbc_df <- fread("Output Files/CSVs/NonBC_Admissions/nbc_with_other_micro.csv",
                colClasses = c(ClusterID = "character"))
nbc_df[, AdmissionDate := as.POSIXct(AdmissionDate, tz = "UTC")]

out_dir <- "Output Files/CSVs/NonBC_Admissions/2024"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 9A. LOAD BC COHORT FOR TRAINING-CONSISTENT PREPROCESSING
# ============================================================

bc_data <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))
bc_data <- bc_data[CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial",
                                                  "Culture-negative", "Contaminant - single", "Contaminant - recurrent",
                                                  "Polymicrobial - contaminant")]

bc_data[, Outcome := fifelse(CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial"), 1L, 0L)]

median_imd_bc <- median(bc_data$SOCIAL_IMDPercentile, na.rm = TRUE)

rm(bc_data)

# ============================================================
# 9B. FEATURE ENGINEERING (MATCHING XGBOOST TRAINING SCRIPT)
# ============================================================

nbc_df[, Sex_num           := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]
nbc_df[, AdmissionType_num := 1L]
nbc_df[, LOS_BeforeCulture_days := 0]
nbc_df[, SourceType_num    := 1L]

nbc_df[is.na(SOCIAL_IMDPercentile), SOCIAL_IMDPercentile := median_imd_bc]

nbc_df[is.na(PositiveCulturesInPrior28Days), PositiveCulturesInPrior28Days := 0L]
nbc_df[, AnyGramPos_Prior := fifelse(GramPositiveCulturesInPrior28Days > 0, 1, 0)]
nbc_df[, AnyGramNeg_Prior := fifelse(GramNegativeCulturesInPrior28Days > 0, 1, 0)]

for (col in c("Source_Urinary", "Source_Respiratory", "Source_Skin_Soft_Tissue", "Source_Other")) {
  if (col %in% names(nbc_df)) nbc_df[is.na(get(col)), (col) := 0L]
}

for (col in c("HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before", "Carbapenem_Before",
              "Any_Prior_Other_Pos", "Prior_Urine_Pos")) {
  if (col %in% names(nbc_df)) nbc_df[is.na(get(col)), (col) := 0L]
}

nbc_df[, VS_AnyFever       := fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), 1, 0)]
nbc_df[, VS_AnyHypotension := fifelse(VS_AnyHypotension == 1 & !is.na(VS_AnyHypotension), 1, 0)]

# ============================================================
# 9C. BUILD FEATURE MATRIX (32 FEATURES, SAME ORDER AS TRAINING)
# ============================================================

feature_cols <- c(
  "AgeAtCulture", "Sex_num",
  "AdmissionType_num", "LOS_BeforeCulture_days",
  "SourceType_num",
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
for (col in feature_cols) {
  if (!col %in% names(nbc_df)) nbc_df[, (col) := NA_real_]
}

# ============================================================
# 9D. LOAD MODEL AND APPLY PLATT CALIBRATION
# ============================================================

model       <- xgb.load("Models/XGBoost/FullCohort/xgb_fullcohort.json")
platt_coefs <- readRDS("Models/XGBoost/FullCohort/xgb_fullcohort_platt_calibration.rds")

pred_matrix <- as.matrix(nbc_df[, ..feature_cols])
dpred       <- xgb.DMatrix(data = pred_matrix)
pred_raw    <- predict(model, dpred)

pred_logit      <- qlogis(pmax(pmin(pred_raw, 0.999), 0.001))
pred_calibrated <- plogis(platt_coefs[1] + platt_coefs[2] * pred_logit)

nbc_df[, pred_raw        := pred_raw]
nbc_df[, pred_calibrated := pred_calibrated]

# ============================================================
# 10. RESULTS SUMMARY
# ============================================================

n_admissions <- nrow(nbc_df)
n_patients   <- uniqueN(nbc_df$ClusterID)

results <- data.table(
  Metric = c(
    "Study Period Start", "Study Period End",
    "Total Emergency Admissions (no BC within 24h)",
    "Unique Patients",
    "Mean Predicted Positive Rate (raw)",
    "Mean Predicted Positive Rate (calibrated)",
    "Median Predicted Positive Rate (calibrated)",
    "SD Predicted Positive Rate (calibrated)",
    "Pct Predicted > 10%", "Pct Predicted > 20%", "Pct Predicted > 30%",
    "Pct with Vital Signs Data", "Pct with Blood Test Data", "Pct with Blood Gas Data",
    "Mean Age at Admission", "Pct Male",
    "Pct with Prior Antibiotics (14d)",
    "Pct with Prior Positive BC (28d)",
    "Pct with Prior Other Positive Micro (28d)"
  ),
  Value = c(
    "2024-01-01", "2024-12-31",
    format(n_admissions, big.mark = ","),
    format(n_patients, big.mark = ","),
    sprintf("%.1f%%", 100 * mean(pred_raw, na.rm = TRUE)),
    sprintf("%.1f%%", 100 * mean(pred_calibrated, na.rm = TRUE)),
    sprintf("%.1f%%", 100 * median(pred_calibrated, na.rm = TRUE)),
    sprintf("%.1f%%", 100 * sd(pred_calibrated, na.rm = TRUE)),
    sprintf("%.1f%%", 100 * mean(pred_calibrated > 0.10)),
    sprintf("%.1f%%", 100 * mean(pred_calibrated > 0.20)),
    sprintf("%.1f%%", 100 * mean(pred_calibrated > 0.30)),
    sprintf("%.1f%%", 100 * mean(!is.na(nbc_df$VS_Baseline_HR))),
    sprintf("%.1f%%", 100 * mean(!is.na(nbc_df$BT_Closest_CRP))),
    sprintf("%.1f%%", 100 * mean(!is.na(nbc_df$BG_Closest_Lactate))),
    sprintf("%.1f",   mean(nbc_df$AgeAtCulture, na.rm = TRUE)),
    sprintf("%.1f%%", 100 * mean(nbc_df$DEMOGRAPHICS_LinkedSex == "M", na.rm = TRUE)),
    sprintf("%.1f%%", 100 * mean(nbc_df$HasAntibiotics_Before == 1)),
    sprintf("%.1f%%", 100 * mean(nbc_df$PositiveCulturesInPrior28Days > 0)),
    sprintf("%.1f%%", 100 * mean(nbc_df$Any_Prior_Other_Pos == 1))
  )
)

pred_dist <- data.table(
  Bin = c("0-5%","5-10%","10-15%","15-20%","20-30%","30-50%",">50%"),
  N   = c(sum(pred_calibrated <= 0.05),
          sum(pred_calibrated > 0.05 & pred_calibrated <= 0.10),
          sum(pred_calibrated > 0.10 & pred_calibrated <= 0.15),
          sum(pred_calibrated > 0.15 & pred_calibrated <= 0.20),
          sum(pred_calibrated > 0.20 & pred_calibrated <= 0.30),
          sum(pred_calibrated > 0.30 & pred_calibrated <= 0.50),
          sum(pred_calibrated > 0.50))
)
pred_dist[, Pct := round(100 * N / sum(N), 1)]

fwrite(results,   file.path(out_dir, "nonbc_admissions_model_results.csv"))
fwrite(pred_dist, file.path(out_dir, "nonbc_admissions_prediction_distribution.csv"))
fwrite(nbc_df,    file.path(out_dir, "nonbc_admissions_with_predictions.csv"))

write_xlsx(list(
  Results              = as.data.frame(results),
  Prediction_Distribution = as.data.frame(pred_dist)
), file.path(out_dir, "NonBC_Admissions_Model_Application.xlsx"))

# ============================================================
# 11. CHARACTERISE HIGH-RISK COHORT (PRED_CALIBRATED > 5%)
# ============================================================

high_risk <- nbc_df[pred_calibrated > 0.05]
fwrite(high_risk, file.path(out_dir, "nonbc_highRisk_gt5pct_full.csv"))

pct_flag <- function(x) sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
med_iqr  <- function(x) {
  if (all(is.na(x))) return("NA")
  sprintf("%.1f [%.1f-%.1f]",
    median(x, na.rm = TRUE),
    quantile(x, 0.25, na.rm = TRUE),
    quantile(x, 0.75, na.rm = TRUE))
}

n_hr    <- nrow(high_risk)
n_pt_hr <- uniqueN(high_risk$ClusterID)

char_rows <- list(
  list(Category = "Overview",          Variable = "Admissions, N",                                Value = format(n_hr, big.mark = ",")),
  list(Category = "Overview",          Variable = "Unique patients, N",                            Value = format(n_pt_hr, big.mark = ",")),
  list(Category = "Overview",          Variable = "Predicted risk (calibrated): median [IQR]",     Value = sprintf("%.1f%% [%.1f%%-%.1f%%]", 100*median(high_risk$pred_calibrated), 100*quantile(high_risk$pred_calibrated,0.25), 100*quantile(high_risk$pred_calibrated,0.75))),
  list(Category = "Demographics",      Variable = "Age at admission: median [IQR]",                Value = med_iqr(high_risk$AgeAtCulture)),
  list(Category = "Demographics",      Variable = "Male sex",                                      Value = pct_flag(high_risk$DEMOGRAPHICS_LinkedSex == "M")),
  list(Category = "Demographics",      Variable = "IMD Percentile: median [IQR]",                  Value = med_iqr(high_risk$SOCIAL_IMDPercentile)),
  list(Category = "Prior BC",          Variable = "Any prior positive BC (28d)",                   Value = pct_flag(high_risk$PositiveCulturesInPrior28Days > 0)),
  list(Category = "Prior BC",          Variable = "Any prior Gram-positive BC (28d)",              Value = pct_flag(high_risk$GramPositiveCulturesInPrior28Days > 0)),
  list(Category = "Prior BC",          Variable = "Any prior Gram-negative BC (28d)",              Value = pct_flag(high_risk$GramNegativeCulturesInPrior28Days > 0)),
  list(Category = "Antibiotics",       Variable = "Prior antibiotics (14d)",                       Value = pct_flag(high_risk$HasAntibiotics_Before == 1)),
  list(Category = "Antibiotics",       Variable = "Broad spectrum (14d)",                          Value = pct_flag(high_risk$BroadSpectrum_Before == 1)),
  list(Category = "Antibiotics",       Variable = "3rd gen+ cephalosporin (14d)",                  Value = pct_flag(high_risk$Ceph3rdGenPlus_Before == 1)),
  list(Category = "Antibiotics",       Variable = "Carbapenem (14d)",                              Value = pct_flag(high_risk$Carbapenem_Before == 1)),
  list(Category = "Antibiotic Source", Variable = "Abx indication: Urinary",                      Value = pct_flag(high_risk$Source_Urinary == 1)),
  list(Category = "Antibiotic Source", Variable = "Abx indication: Respiratory",                   Value = pct_flag(high_risk$Source_Respiratory == 1)),
  list(Category = "Antibiotic Source", Variable = "Abx indication: Skin/Soft tissue",              Value = pct_flag(high_risk$Source_Skin_Soft_Tissue == 1)),
  list(Category = "Antibiotic Source", Variable = "Abx indication: Other",                         Value = pct_flag(high_risk$Source_Other == 1)),
  list(Category = "Prior Other Micro", Variable = "Any prior positive non-BC micro (28d)",         Value = pct_flag(high_risk$Any_Prior_Other_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive urine culture (28d)",            Value = pct_flag(high_risk$Prior_Urine_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive respiratory (28d)",              Value = pct_flag(high_risk$Prior_Resp_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive skin/soft tissue (28d)",         Value = pct_flag(high_risk$Prior_SkinSoftTissue_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive sterile fluid (28d)",            Value = pct_flag(high_risk$Prior_SterileFluid_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive CSF (28d)",                      Value = pct_flag(high_risk$Prior_CSF_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive GI (28d)",                       Value = pct_flag(high_risk$Prior_GI_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive device/line tip (28d)",          Value = pct_flag(high_risk$Prior_Device_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive genital (28d)",                  Value = pct_flag(high_risk$Prior_Genital_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive ear (28d)",                      Value = pct_flag(high_risk$Prior_Ear_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive eye (28d)",                      Value = pct_flag(high_risk$Prior_Eye_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive bone (28d)",                     Value = pct_flag(high_risk$Prior_Bone_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior positive other source (28d)",             Value = pct_flag(high_risk$Prior_Other_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Total prior positive non-BC cultures: median [IQR]", Value = med_iqr(high_risk$Prior_Total_Other_Cultures_Pos)),
  list(Category = "Prior Other Micro", Variable = "Prior Gram-positive non-BC organism (28d)",     Value = pct_flag(high_risk$Prior_GramPos_Pos == 1)),
  list(Category = "Prior Other Micro", Variable = "Prior Gram-negative non-BC organism (28d)",     Value = pct_flag(high_risk$Prior_GramNeg_Pos == 1)),
  list(Category = "Vital Signs",       Variable = "NEWS: median [IQR]",                            Value = med_iqr(high_risk$VS_Baseline_NEWS)),
  list(Category = "Vital Signs",       Variable = "Heart rate: median [IQR]",                      Value = med_iqr(high_risk$VS_Baseline_HR)),
  list(Category = "Vital Signs",       Variable = "Temperature: median [IQR]",                     Value = med_iqr(high_risk$VS_Baseline_Temp)),
  list(Category = "Vital Signs",       Variable = "Systolic BP: median [IQR]",                     Value = med_iqr(high_risk$VS_Baseline_SBP)),
  list(Category = "Vital Signs",       Variable = "Any fever >=38.0 C",                            Value = pct_flag(high_risk$VS_AnyFever_38 == 1)),
  list(Category = "Vital Signs",       Variable = "Any hypotension (SBP <90)",                     Value = pct_flag(high_risk$VS_AnyHypotension == 1)),
  list(Category = "Blood Tests",       Variable = "CRP (mg/L): median [IQR]",                      Value = med_iqr(high_risk$BT_Closest_CRP)),
  list(Category = "Blood Tests",       Variable = "WCC (x10^9/L): median [IQR]",                   Value = med_iqr(high_risk$BT_Closest_WCC)),
  list(Category = "Blood Tests",       Variable = "Neutrophils (x10^9/L): median [IQR]",           Value = med_iqr(high_risk$BT_Closest_Neut)),
  list(Category = "Blood Tests",       Variable = "Lymphocytes (x10^9/L): median [IQR]",           Value = med_iqr(high_risk$BT_Closest_Lymph)),
  list(Category = "Blood Tests",       Variable = "Albumin (g/L): median [IQR]",                   Value = med_iqr(high_risk$BT_Closest_Alb)),
  list(Category = "Blood Tests",       Variable = "Platelets (x10^9/L): median [IQR]",             Value = med_iqr(high_risk$BT_Closest_Plt)),
  list(Category = "Blood Tests",       Variable = "Creatinine (umol/L): median [IQR]",             Value = med_iqr(high_risk$BT_Closest_Creat)),
  list(Category = "Blood Tests",       Variable = "Bilirubin (umol/L): median [IQR]",              Value = med_iqr(high_risk$BT_Closest_Bili)),
  list(Category = "Blood Tests",       Variable = "Urea (mmol/L): median [IQR]",                   Value = med_iqr(high_risk$BT_Closest_Urea)),
  list(Category = "Blood Tests",       Variable = "Haemoglobin (g/L): median [IQR]",               Value = med_iqr(high_risk$BT_Closest_Hb)),
  list(Category = "Blood Tests",       Variable = "Sodium (mmol/L): median [IQR]",                 Value = med_iqr(high_risk$BT_Closest_Na)),
  list(Category = "Blood Tests",       Variable = "Eosinophils (x10^9/L): median [IQR]",           Value = med_iqr(high_risk$BT_Closest_Eosin)),
  list(Category = "Blood Tests",       Variable = "Monocytes (x10^9/L): median [IQR]",             Value = med_iqr(high_risk$BT_Closest_Mono)),
  list(Category = "Blood Tests",       Variable = "ALP (IU/L): median [IQR]",                      Value = med_iqr(high_risk$BT_Closest_ALP)),
  list(Category = "Blood Tests",       Variable = "Haematocrit: median [IQR]",                      Value = med_iqr(high_risk$BT_Closest_HCT)),
  list(Category = "Blood Gas",         Variable = "Lactate (mmol/L): median [IQR]",                Value = med_iqr(high_risk$BG_Closest_Lactate)),
  list(Category = "Blood Gas",         Variable = "pH: median [IQR]",                              Value = med_iqr(high_risk$BG_Closest_pH)),
  list(Category = "Blood Gas",         Variable = "pCO2 (kPa): median [IQR]",                      Value = med_iqr(high_risk$BG_Closest_pCO2))
)

char_dt <- rbindlist(lapply(char_rows, as.data.table))

pred_dist_hr <- data.table(
  Bin = c("5-10%","10-15%","15-20%","20-30%","30-50%",">50%"),
  N   = c(sum(high_risk$pred_calibrated > 0.05 & high_risk$pred_calibrated <= 0.10),
          sum(high_risk$pred_calibrated > 0.10 & high_risk$pred_calibrated <= 0.15),
          sum(high_risk$pred_calibrated > 0.15 & high_risk$pred_calibrated <= 0.20),
          sum(high_risk$pred_calibrated > 0.20 & high_risk$pred_calibrated <= 0.30),
          sum(high_risk$pred_calibrated > 0.30 & high_risk$pred_calibrated <= 0.50),
          sum(high_risk$pred_calibrated > 0.50))
)
pred_dist_hr[, Pct := round(100 * N / sum(N), 1)]

write_xlsx(list(
  Characterisation        = as.data.frame(char_dt),
  Prediction_Distribution = as.data.frame(pred_dist_hr)
), file.path(out_dir, "NonBC_HighRisk_gt5pct_Characterisation.xlsx"))

# ============================================================
# 12. COMBINED RISK DISTRIBUTION — BC Positive / BC Negative / NBC Admission (2024)
# ============================================================

fig_dir_nbc <- "Output Files/Figures/NBC/NonBC_Admissions"
dir.create(fig_dir_nbc, showWarnings = FALSE, recursive = TRUE)

bc_pred_path <- "Output Files/CSVs/XGBoost/FullCohort/test_predictions.csv"

if (file.exists(bc_pred_path)) {
  bc_preds <- fread(bc_pred_path)

  plot_df <- rbind(
    data.table(Prediction = bc_preds[outcome == 1]$pred,  Group = "BC Positive"),
    data.table(Prediction = bc_preds[outcome == 0]$pred,  Group = "BC Negative"),
    data.table(Prediction = pred_calibrated,               Group = "NBC Admission")
  )
  plot_df[, Group := factor(Group, levels = c("BC Positive", "BC Negative", "NBC Admission"))]

  fill_vals <- c("BC Positive" = "#E94F37", "BC Negative" = "#2E86AB", "NBC Admission" = "#5C9E6B")

  density_base <- list(
    geom_density(alpha = 0.4, linewidth = 1.2, bw = 0.01),
    scale_fill_manual(values = fill_vals),
    scale_color_manual(values = fill_vals),
    labs(fill = "Group", color = "Group", x = "Predicted Probability", y = "Density"),
    theme_minimal(),
    theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")
  )

  p_combined <- ggplot(plot_df, aes(x = Prediction, fill = Group, color = Group)) +
    density_base +
    scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.1)) +
    labs(title = "Risk Distribution - Predicted Probability by Group (2024)",
         subtitle = "BC positive cultures, BC negative cultures, and NBC emergency admissions")

  ggsave(file.path(fig_dir_nbc, "01_risk_density_combined.png"),
         p_combined, width = 10, height = 7, dpi = 300)

  # Distribution table (risk bins)
  risk_bins  <- c(0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.001)
  bin_labels <- c("0-5%", "5-10%", "10-15%", "15-20%", "20-30%", "30-50%", ">50%")

  make_risk_bins <- function(preds, grp) {
    b  <- cut(preds, breaks = risk_bins, labels = bin_labels, right = FALSE, include.lowest = TRUE)
    dt <- as.data.table(table(Bin = b))
    dt[, Pct   := round(100 * N / sum(N), 1)]
    dt[, Group := grp]
    dt[, .(Group, Bin, N, Pct)]
  }

  risk_dist_table <- rbindlist(list(
    make_risk_bins(bc_preds[outcome == 1]$pred, "BC Positive"),
    make_risk_bins(bc_preds[outcome == 0]$pred, "BC Negative"),
    make_risk_bins(pred_calibrated,              "NBC Admission")
  ))

  fwrite(risk_dist_table, file.path(fig_dir_nbc, "01_risk_distribution_by_group.csv"))

  write_xlsx(list(Risk_Distribution = as.data.frame(risk_dist_table)),
             file.path(fig_dir_nbc, "01_risk_distribution_by_group.xlsx"))
}

# ============================================================
# 13. WARD-LEVEL BC vs NBC COMPARISON
# Join NBC predictions to WardStays to get admitting ward,
# then compare with pre-saved BC ward performance.
# ============================================================

bc_ward_path   <- "Output Files/CSVs/XGBoost/FullCohort/performance_by_ward_enriched.csv"
ward_stays_path <- "Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_InpatWardStays.csv"

if (file.exists(bc_ward_path) && file.exists(ward_stays_path)) {

  bc_ward <- fread(bc_ward_path)

  ward_stays <- fread(ward_stays_path,
                      colClasses = c(ClusterID = "character"),
                      select = c("SpellID", "WardStartDate", "WardName"))
  ward_stays[, SpellID       := as.character(trimws(SpellID))]
  ward_stays[, WardStartDate := as.POSIXct(WardStartDate, format = "%Y-%m-%d %H:%M:%S")]

  nbc_df[, SpellID := as.character(trimws(SpellID))]
  first_ward <- ward_stays[order(SpellID, WardStartDate)][, .SD[1], by = SpellID]

  nbc_ward_joined <- merge(
    nbc_df[, .(SpellID, pred_calibrated)],
    first_ward[, .(SpellID, WardName)],
    by = "SpellID", all.x = TRUE
  )
  nbc_ward_joined[, WardName := trimws(as.character(WardName))]

  # Load Youden threshold from saved XGBoost thresholds
  thresh_path <- "Output Files/CSVs/XGBoost/FullCohort/sensitivity_locked_thresholds.csv"
  youden_threshold <- if (file.exists(thresh_path)) {
    fread(thresh_path)[Operating_Point == "Youden (optimal)"]$Threshold
  } else { 0.06 }

  nbc_ward_summary <- nbc_ward_joined[
    !is.na(WardName) & WardName != "",
    .(N_NBC              = .N,
      NBC_HighRisk_Pct   = round(100 * mean(pred_calibrated > 0.05), 1),
      NBC_Mean_Risk_Pct  = round(100 * mean(pred_calibrated), 1),
      NBC_Flagged_N      = sum(pred_calibrated >= youden_threshold),
      NBC_Expected_Pos   = round(sum(pred_calibrated[pred_calibrated >= youden_threshold]))),
    by = WardName
  ][N_NBC >= 1000][order(-NBC_Mean_Risk_Pct)]

  fwrite(nbc_ward_summary,
         file.path(out_dir, "nbc_ward_predicted_positivity.csv"))

  # Merge NBC summary with BC ward performance (all NBC wards; BC data where available)
  ward_comp <- merge(
    nbc_ward_summary[, .(WardName, N_NBC, NBC_HighRisk_Pct, NBC_Mean_Risk_Pct,
                         NBC_Flagged_N, NBC_Expected_Pos)],
    bc_ward[, .(WardName = Group, N_BC = N, N_BC_Pos = N_Positive,
                BC_Pos_Pct = Prevalence, BC_Flag_Pct = Flagged_Pct)],
    by = "WardName", all.x = TRUE
  )
  ward_comp[, Has_BC := !is.na(N_BC)]
  ward_comp[, Label  := factor(WardName, levels = ward_comp[order(NBC_Mean_Risk_Pct)]$WardName)]

  fwrite(ward_comp, file.path(out_dir, "ward_bc_nbc_comparison.csv"))

  library(ggplot2)
  library(cowplot)

  panel_theme_mid <- function() {
    theme_minimal() %+replace%
      theme(plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.minor = element_blank())
  }

  col_nbc_n  <- "#5C9E6B"
  col_nbc_hr <- "#F4A261"
  col_bc_n   <- "#2E86AB"
  col_bc_pos <- "#E94F37"
  col_na     <- "gray80"

  # BC panels first — pC is leftmost so it carries the y-axis ward labels
  ward_comp[, N_BC_plot := fifelse(is.na(N_BC), 0L, as.integer(N_BC))]
  pC <- ggplot(ward_comp, aes(x = Label, y = N_BC_plot, fill = Has_BC)) +
    geom_col(width = 0.7) +
    geom_text(data = ward_comp[Has_BC == TRUE],
              aes(label = format(N_BC, big.mark = ",")), hjust = -0.1, size = 2.6) +
    geom_text(data = ward_comp[Has_BC == FALSE],
              aes(x = Label, y = 0, label = "—"), hjust = -0.1,
              size = 3, color = "gray60", inherit.aes = FALSE) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = col_bc_n, "FALSE" = col_na), guide = "none") +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) +
    labs(x = "", y = "N Cultures Sent", title = "BC Volume") +
    theme_minimal() +
    theme(plot.title  = element_text(face = "bold", size = 11, hjust = 0.5),
          axis.text.y = element_text(size = 8),
          panel.grid.minor = element_blank())

  # NBC panels — no y-axis labels
  pA <- ggplot(ward_comp, aes(x = Label, y = N_NBC)) +
    geom_col(fill = col_nbc_n, width = 0.7) +
    geom_text(aes(label = format(N_NBC, big.mark = ",")), hjust = -0.1, size = 2.6) +
    coord_flip() +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) +
    labs(x = "", y = "N Admissions (no BC)", title = "NBC Volume") +
    panel_theme_mid()

  pB <- ggplot(ward_comp, aes(x = Label, y = NBC_Mean_Risk_Pct)) +
    geom_col(fill = col_nbc_hr, width = 0.7, alpha = 0.85) +
    geom_text(aes(label = paste0(NBC_Mean_Risk_Pct, "%")), hjust = -0.1, size = 2.6) +
    coord_flip() +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.25))) +
    labs(x = "", y = "Mean Predicted Risk %", title = "NBC Mean Predicted Risk") +
    panel_theme_mid()

  pD <- ggplot(ward_comp,
               aes(x = Label, y = fifelse(is.na(BC_Pos_Pct), 0, BC_Pos_Pct),
                   fill = Has_BC)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_text(data = ward_comp[Has_BC == TRUE],
              aes(label = paste0(BC_Pos_Pct, "%")), hjust = -0.1, size = 2.6) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = col_bc_pos, "FALSE" = col_na), guide = "none") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.25))) +
    labs(x = "", y = "BC Positivity %", title = "BC Actual % Positive") +
    panel_theme_mid()

  # Panel E: actual N BC positive
  col_bc_pos_n <- "#C0392B"
  pE <- ggplot(ward_comp, aes(x = Label, y = fifelse(is.na(N_BC_Pos), 0L, N_BC_Pos),
                               fill = Has_BC)) +
    geom_col(width = 0.7) +
    geom_text(data = ward_comp[Has_BC == TRUE],
              aes(label = N_BC_Pos), hjust = -0.1, size = 2.6) +
    geom_text(data = ward_comp[Has_BC == FALSE],
              aes(x = Label, y = 0, label = "\u2014"), hjust = -0.1,
              size = 3, color = "gray60", inherit.aes = FALSE) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = col_bc_pos_n, "FALSE" = col_na), guide = "none") +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) +
    labs(x = "", y = "N BC Positive", title = "BC Actual Positives") +
    panel_theme_mid()

  # Panel F: NBC expected positives if above-Youden patients were cultured
  col_nbc_ep <- "#2D6A4F"
  pF <- ggplot(ward_comp, aes(x = Label, y = NBC_Expected_Pos)) +
    geom_col(fill = col_nbc_ep, width = 0.7, alpha = 0.85) +
    geom_text(aes(label = NBC_Expected_Pos), hjust = -0.1, size = 2.6) +
    coord_flip() +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) +
    labs(x = "", y = "Expected Positives",
         title = sprintf("NBC Expected Positives\n(if pred \u2265%.3f cultured)", youden_threshold)) +
    panel_theme_mid() +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5))

  title_grob <- ggdraw() + draw_label(
    sprintf("Ward Comparison: NBC Admissions vs BC Performance (2024)\nNBC wards with \u22651,000 admissions | Ordered by NBC mean predicted risk | Youden threshold = %.3f",
            youden_threshold),
    fontface = "bold", size = 12, x = 0.5, hjust = 0.5)

  body <- plot_grid(pC, pD, pE, pA, pB, pF, nrow = 1,
                    rel_widths = c(1.5, 1, 0.8, 1, 1, 0.8), align = "h", axis = "tb")
  fig_out <- plot_grid(title_grob, body, ncol = 1, rel_heights = c(0.05, 1))

  row_h <- max(7, nrow(ward_comp) * 0.45 + 3)
  ggsave(file.path(fig_dir_nbc, "02_ward_bc_nbc_comparison.png"),
         fig_out, width = 22, height = row_h, dpi = 300)
}
