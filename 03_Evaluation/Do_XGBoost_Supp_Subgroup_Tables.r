################################################################################
# XGBoost Full Cohort — Supplementary tables for Supplementary Figures 2 and 3
#
# Full Table 2 metrics (AUROC with DeLong CI, AUPRC, sensitivity, specificity,
# PPV, NPV with Wilson CIs) and confusion matrix counts at the Youden threshold
# for the data-completeness patterns (Supp Fig 2) and subgroups (Supp Fig 3),
# plus confusion matrices for the full test set across operating thresholds.
################################################################################

library(data.table)
library(xgboost)
library(pROC)
library(PRROC)
library(writexl)

csv_dir   <- "Output Files/CSVs/XGBoost/FullCohort"
excel_dir <- "Output Files/Publication_Tables_Figures/Per_Model_Summaries"
model_dir <- "Models/XGBoost/FullCohort"

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

################################################################################
# DATA + PREDICTIONS (feature engineering identical to Do_XGBoost_FullCohort.r)
################################################################################

bc <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
            colClasses = c(ClusterID = "character"))
bc <- bc[CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial",
                                       "Fungi", "Anaerobes", "Other",
                                       "Culture-negative", "Contaminant - single",
                                       "Contaminant - recurrent", "Polymicrobial - contaminant")]
bc[, Outcome := fifelse(CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial",
                                                      "Fungi", "Anaerobes", "Other"), 1L, 0L)]
bc[, Sex_num := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]
bc[is.na(SOCIAL_IMDPercentile), SOCIAL_IMDPercentile := median(bc$SOCIAL_IMDPercentile, na.rm = TRUE)]
bc[, AdmissionType_num := fcase(AdmissionMethod == "Emergency", 1,
                                AdmissionMethod == "Elective", 2, default = 0)]
bc[, SourceType_num := fcase(
  CultureSourceCategory == "Peripheral", 1,
  CultureSourceCategory %in% c("Central Line", "PICC", "Port/Portacath", "Arterial", "Femoral"), 2,
  default = 0)]
zero_fill <- c("Source_Respiratory", "PositiveCulturesInPrior28Days",
               "HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before",
               "Carbapenem_Before", "Any_Prior_Other_Pos", "Prior_Urine_Pos", "LOS_BeforeCulture_days")
for (col in zero_fill) set(bc, which(is.na(bc[[col]])), col, 0)
bc[, VS_AnyFever := fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), 1, 0)]

test <- bc[as.Date(CollectionDateTime) >= as.Date("2024-01-01")]

model <- xgb.load(file.path(model_dir, "xgb_fullcohort.json"))
platt <- readRDS(file.path(model_dir, "xgb_fullcohort_platt_calibration.rds"))
raw   <- predict(model, xgb.DMatrix(as.matrix(test[, ..feature_cols])))
test[, pred := plogis(platt[1] + platt[2] * qlogis(pmax(pmin(raw, 0.999), 0.001)))]

threshold <- coords(roc(test$Outcome, test$pred, quiet = TRUE), "best", ret = "threshold")$threshold

################################################################################
# SUBSET DEFINITIONS (as in Do_XGBoost_FullCohort_FiguresAndTables.r)
################################################################################

test[, VS := !is.na(VS_Baseline_Temp) & !is.na(VS_Baseline_HR) & !is.na(VS_Baseline_SBP) & !is.na(VS_Baseline_RR)]
test[, RB := !is.na(BT_Closest_WCC) & !is.na(BT_Closest_Neut) & !is.na(BT_Closest_Lymph) &
             !is.na(BT_Closest_Hb) & !is.na(BT_Closest_Plt) &
             !is.na(BT_Closest_Na) & !is.na(BT_Closest_Creat) & !is.na(BT_Closest_Urea) &
             !is.na(BT_Closest_Bili) & !is.na(BT_Closest_ALP) & !is.na(BT_Closest_Alb) &
             !is.na(BT_Closest_CRP)]
test[, BG := !is.na(BG_Closest_Lactate) & !is.na(BG_Closest_pH)]

completeness <- list(
  "All complete"                              = quote(VS & RB & BG),
  "Vital signs + bloods complete, blood gas incomplete" = quote(VS & RB & !BG),
  "Vital signs + blood gas complete, bloods incomplete" = quote(VS & !RB & BG),
  "Bloods + blood gas complete, vital signs incomplete" = quote(!VS & RB & BG),
  "Only vital signs complete"                 = quote(VS & !RB & !BG),
  "Only bloods complete"                      = quote(!VS & RB & !BG),
  "Only blood gas complete"                   = quote(!VS & !RB & BG),
  "None complete"                             = quote(!VS & !RB & !BG)
)

test[, AgeGroup := fcase(AgeAtCulture < 18, "Paediatric (<18)", AgeAtCulture < 65, "Adult (18-64)",
                         AgeAtCulture >= 65, "Elderly (65+)")]
test[, Setting := fcase(CultureLocation == "ED", "Emergency department",
                        CultureLocation == "Inpatient Ward", "Inpatient", default = "Other")]
test[, Acquisition := fifelse(LOS_BeforeCulture_days < 2, "Community-acquired (<48h)", "Hospital-acquired (>=48h)")]
test[, IMD := fcase(SOCIAL_IMDPercentile <= 20, "Q1 (most deprived)", SOCIAL_IMDPercentile <= 40, "Q2",
                    SOCIAL_IMDPercentile <= 60, "Q3", SOCIAL_IMDPercentile <= 80, "Q4",
                    SOCIAL_IMDPercentile <= 100, "Q5 (least deprived)")]
test[, Sex := fifelse(DEMOGRAPHICS_LinkedSex == "M", "Male", "Female")]

subgroups <- list(
  "Age group"    = list(col = "AgeGroup",    vals = c("Paediatric (<18)", "Adult (18-64)", "Elderly (65+)")),
  "Sex"          = list(col = "Sex",         vals = c("Female", "Male")),
  "Care setting" = list(col = "Setting",     vals = c("Emergency department", "Inpatient", "Other")),
  "Acquisition"  = list(col = "Acquisition", vals = c("Community-acquired (<48h)", "Hospital-acquired (>=48h)")),
  "IMD quintile" = list(col = "IMD",         vals = c("Q1 (most deprived)", "Q2", "Q3", "Q4", "Q5 (least deprived)")),
  "Ethnicity"    = list(col = "SOCIAL_Ethnicity", vals = c("White", "Asian", "Black", "Mixed", "Other", "Unknown"))
)

################################################################################
# METRICS
################################################################################

wilson <- function(x, n) {
  if (n == 0) return(c(NA, NA, NA))
  ci <- prop.test(x, n, correct = FALSE)$conf.int
  100 * c(x / n, ci[1], ci[2])
}
fmt_pct <- function(v) if (is.na(v[1])) "—" else sprintf("%.1f%% (%.1f–%.1f)", v[1], v[2], v[3])

metrics <- function(d, group, subgroup, thr = threshold) {
  y <- d$Outcome == 1; p <- d$pred; flag <- p >= thr
  TP <- sum(flag & y); FP <- sum(flag & !y); TN <- sum(!flag & !y); FN <- sum(!flag & y)
  ok <- sum(y) >= 10 && sum(!y) >= 10
  r  <- if (ok) roc(d$Outcome, p, quiet = TRUE, direction = "<") else NULL
  ci <- if (ok) as.numeric(ci.auc(r, method = "delong")) else rep(NA, 3)
  ap <- if (ok) pr.curve(scores.class0 = p[y], scores.class1 = p[!y])$auc.integral else NA
  data.table(Group = group, Subgroup = subgroup, N = nrow(d), Positives = sum(y),
             Prevalence = 100 * mean(y), AUROC = ci[2], AUROC_Lo = ci[1], AUROC_Hi = ci[3], AUPRC = ap,
             Sens = list(wilson(TP, TP + FN)), Spec = list(wilson(TN, TN + FP)),
             PPV = list(wilson(TP, TP + FP)), NPV = list(wilson(TN, TN + FN)),
             TP = TP, FP = FP, TN = TN, FN = FN)
}

overall <- metrics(test, "Overall", "All test cultures")

supp_fig2 <- rbind(overall, rbindlist(lapply(names(completeness), function(nm)
  metrics(test[eval(completeness[[nm]])], "Data completeness", nm))))

supp_fig3 <- rbind(overall, rbindlist(lapply(names(subgroups), function(g) {
  s <- subgroups[[g]]
  rbindlist(lapply(s$vals, function(v) metrics(test[get(s$col) == v], g, v)))
})))

thresholds <- c(0.01, 0.02, threshold, 0.10, 0.20)
thr_tab <- rbindlist(lapply(thresholds, function(t)
  metrics(test, "Overall", sprintf("%.3f", t), thr = t)))

################################################################################
# OUTPUT
################################################################################

raw_cols <- function(dt) dt[, .(Group, Subgroup, N, Positives, Prevalence = round(Prevalence, 2),
  AUROC = round(AUROC, 3), AUROC_Lo = round(AUROC_Lo, 3), AUROC_Hi = round(AUROC_Hi, 3), AUPRC = round(AUPRC, 3),
  Sensitivity = round(sapply(Sens, `[`, 1), 2), Specificity = round(sapply(Spec, `[`, 1), 2),
  PPV = round(sapply(PPV, `[`, 1), 2), NPV = round(sapply(NPV, `[`, 1), 2), TP, FP, TN, FN)]

formatted <- function(dt, thr_label = FALSE) dt[, .(
  Group, Subgroup,
  N           = formatC(N, format = "d", big.mark = ","),
  `Positives (%)` = sprintf("%s (%.1f%%)", formatC(Positives, format = "d", big.mark = ","), Prevalence),
  `AUROC (95% CI)` = fifelse(is.na(AUROC), "—", sprintf("%.3f (%.3f–%.3f)", AUROC, AUROC_Lo, AUROC_Hi)),
  AUPRC       = fifelse(is.na(AUPRC), "—", sprintf("%.3f", AUPRC)),
  `Sensitivity (95% CI)` = sapply(Sens, fmt_pct), `Specificity (95% CI)` = sapply(Spec, fmt_pct),
  `PPV (95% CI)` = sapply(PPV, fmt_pct), `NPV (95% CI)` = sapply(NPV, fmt_pct),
  TP = formatC(TP, format = "d", big.mark = ","), FP = formatC(FP, format = "d", big.mark = ","),
  TN = formatC(TN, format = "d", big.mark = ","), FN = formatC(FN, format = "d", big.mark = ","))]

thr_tab[, `:=`(Group = "Operating threshold", AUROC = NA, AUPRC = NA,
               Subgroup = fifelse(Subgroup == sprintf("%.3f", threshold),
                                  sprintf("%.3f (Youden)", threshold), Subgroup))]
combined <- rbind(supp_fig2, supp_fig3[Group != "Overall"], thr_tab)

fwrite(raw_cols(combined),  file.path(csv_dir, "supp_table_subset_metrics.csv"))
fwrite(formatted(combined), file.path(csv_dir, "supp_table_subset_metrics_formatted.csv"))
write_xlsx(list(Summary_Formatted = formatted(combined), Summary_Raw = raw_cols(combined)),
           file.path(excel_dir, "XGBoost_Subset_Metrics.xlsx"))
