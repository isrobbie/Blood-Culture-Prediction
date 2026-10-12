################################################################################
# XGBoost Full Cohort — SHAP Consistency Across Data Subsets
#
# Mean |SHAP| for the final model (no retraining) on: training set, full test
# set, the 5,000-culture test sample used for Figure 2, test subgroups, and the
# 2024+ front-door (strict timing) cohort. Rankings are compared to the full
# test set via Spearman correlation and top-10 overlap.
################################################################################

library(data.table)
library(xgboost)
library(writexl)

csv_dir   <- "Output Files/CSVs/XGBoost/FullCohort"
excel_dir <- "Output Files/Publication_Tables_Figures/Per_Model_Summaries"
model     <- xgb.load("Models/XGBoost/FullCohort/xgb_fullcohort.json")

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

# Feature engineering identical to Do_XGBoost_FullCohort.r
prep <- function(path) {
  dt <- fread(path, colClasses = c(ClusterID = "character"))
  dt <- dt[CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial",
                                         "Fungi", "Anaerobes", "Other",
                                         "Culture-negative", "Contaminant - single",
                                         "Contaminant - recurrent", "Polymicrobial - contaminant")]
  dt[, Sex_num := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]
  dt[, AdmissionType_num := fcase(AdmissionMethod == "Emergency", 1,
                                  AdmissionMethod == "Elective", 2, default = 0)]
  dt[, SourceType_num := fcase(
    CultureSourceCategory == "Peripheral", 1,
    CultureSourceCategory %in% c("Central Line", "PICC", "Port/Portacath", "Arterial", "Femoral"), 2,
    default = 0)]
  zero_fill <- c("Source_Respiratory", "PositiveCulturesInPrior28Days",
                 "HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before",
                 "Carbapenem_Before", "Any_Prior_Other_Pos", "Prior_Urine_Pos",
                 "LOS_BeforeCulture_days")
  for (col in zero_fill) set(dt, which(is.na(dt[[col]])), col, 0)
  dt[, VS_AnyFever := fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), 1, 0)]
  dt[, CollectionDate := as.Date(CollectionDateTime)]
  dt
}

split_date <- as.Date("2024-01-01")
bc   <- prep("Output Files/CSVs/blood_cultures_with_other_micro.csv")
fd   <- prep("Output Files/CSVs/trial_blood_cultures_with_other_micro.csv")[CollectionDate >= split_date]
test <- bc[CollectionDate >= split_date]

set.seed(42)
fig2_idx <- sample(seq_len(nrow(test)), min(5000, nrow(test)))

subsets <- list(
  "Training set (pre-2024)"        = bc[CollectionDate < split_date],
  "Test set (all)"                 = test,
  "Test set (Figure 2 sample)"     = test[fig2_idx],
  "Test: Emergency Department"     = test[CultureLocation == "ED"],
  "Test: Inpatient ward"           = test[CultureLocation == "Inpatient Ward"],
  "Test: Paediatric (<18)"         = test[AgeAtCulture < 18],
  "Test: Adult (18+)"              = test[AgeAtCulture >= 18],
  "Test: Front-door timing"        = fd
)

mean_abs_shap <- function(dt) {
  contrib <- predict(model, xgb.DMatrix(as.matrix(dt[, ..feature_cols])), predcontrib = TRUE)
  colMeans(abs(contrib[, feature_cols, drop = FALSE]))
}

shap_long <- rbindlist(lapply(names(subsets), function(s) {
  v <- mean_abs_shap(subsets[[s]])
  data.table(Subset = s, N = nrow(subsets[[s]]), Feature = names(v),
             Mean_Abs_SHAP = round(v, 4), Rank = frank(-v, ties.method = "min"))
}))

ref     <- shap_long[Subset == "Test set (all)"]
ref_top <- ref[Rank <= 10, Feature]

consistency <- shap_long[, {
  m <- merge(.SD, ref[, .(Feature, Ref_Rank = Rank)], by = "Feature")
  .(N               = N[1],
    Spearman_vs_Test = round(cor(m$Rank, m$Ref_Rank, method = "spearman"), 3),
    Top10_Overlap    = sum(Feature[Rank <= 10] %in% ref_top),
    Top5             = paste(Feature[order(Rank)][1:5], collapse = ", "))
}, by = Subset]

rank_wide <- dcast(shap_long, Feature ~ Subset, value.var = "Rank")
setcolorder(rank_wide, c("Feature", names(subsets)))
setorderv(rank_wide, "Test set (all)")

fwrite(shap_long,   file.path(csv_dir, "shap_consistency_long.csv"))
fwrite(consistency, file.path(csv_dir, "shap_consistency_summary.csv"))
write_xlsx(list(Summary_Formatted = consistency, Summary_Raw = shap_long, Rank_By_Subset = rank_wide),
           file.path(excel_dir, "XGBoost_SHAP_Consistency.xlsx"))
