library(data.table)
library(xgboost)
library(ggplot2)

# ============================================================
# FIGURE 2 - SHAP summary plot (feature contributions to model
#            predictions, ordered by mean absolute SHAP value).
# Reconstructed from the saved model (the original plotting script
# was not preserved). Per-observation SHAP contributions are not
# persisted anywhere (they would embed patient-level feature values),
# so this script recomputes them from the trained model applied to
# the same temporal hold-out test set used in 02_Models/Do_XGBoost_FullCohort.r.
# Run from the repository root; requires access to the underlying
# IORD data extract used to build blood_cultures_with_other_micro.csv.
# ============================================================

fig_dir <- "Output Files/Figures/Report"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

bc_data <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))

bc_data <- bc_data[CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial",
                                                  "Fungi", "Anaerobes", "Other",
                                                  "Culture-negative", "Contaminant - single", "Contaminant - recurrent",
                                                  "Polymicrobial - contaminant")]

bc_data[, Sex_num := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]
median_imd <- median(bc_data$SOCIAL_IMDPercentile, na.rm = TRUE)
bc_data[is.na(SOCIAL_IMDPercentile), SOCIAL_IMDPercentile := median_imd]
bc_data[, AdmissionType_num := fcase(
  AdmissionMethod == "Emergency", 1, AdmissionMethod == "Elective", 2, default = 0)]
bc_data[, SourceType_num := fcase(
  CultureSourceCategory == "Peripheral", 1,
  CultureSourceCategory %in% c("Central Line", "PICC", "Port/Portacath", "Arterial", "Femoral"), 2,
  default = 0)]
for (col in c("Source_Urinary", "Source_Respiratory", "Source_Skin_Soft_Tissue", "Source_Other")) {
  bc_data[is.na(get(col)), (col) := 0]
}
bc_data[is.na(PositiveCulturesInPrior28Days), PositiveCulturesInPrior28Days := 0]
for (col in c("HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before", "Carbapenem_Before")) {
  bc_data[is.na(get(col)), (col) := 0]
}
for (col in c("Any_Prior_Other_Pos", "Prior_Urine_Pos")) {
  bc_data[is.na(get(col)), (col) := 0]
}
bc_data[, VS_AnyFever := fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), 1, 0)]
bc_data[is.na(LOS_BeforeCulture_days), LOS_BeforeCulture_days := 0]

bc_data[, CollectionDate := as.Date(CollectionDateTime)]
test_data <- bc_data[CollectionDate >= as.Date("2024-01-01")]

feature_cols <- c(
  "AgeAtCulture", "Sex_num", "AdmissionType_num", "LOS_BeforeCulture_days", "SourceType_num",
  "PositiveCulturesInPrior28Days", "Source_Respiratory",
  "HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before", "Carbapenem_Before",
  "Any_Prior_Other_Pos", "Prior_Urine_Pos",
  "VS_Baseline_Temp", "VS_Baseline_SBP", "VS_Baseline_NEWS", "VS_AnyFever", "VS_Baseline_HR",
  "BT_Closest_CRP", "BT_Closest_WCC", "BT_Closest_Neut", "BT_Closest_Lymph",
  "BT_Closest_Alb", "BT_Closest_Plt", "BT_Closest_Creat", "BT_Closest_Bili", "BT_Closest_Na",
  "BT_Closest_Eosin", "BT_Closest_Mono", "BT_Closest_ALP",
  "BG_Closest_Lactate", "BG_Closest_pH"
)

test_matrix <- as.matrix(test_data[, ..feature_cols])
model <- xgb.load("Models/XGBoost/FullCohort/xgb_fullcohort.json")

# ---- SHAP contributions (per observation, per feature) ------

contrib <- predict(model, test_matrix, predcontrib = TRUE)
contrib <- contrib[, feature_cols, drop = FALSE]

top_n <- 20
mean_abs_shap <- sort(colMeans(abs(contrib)), decreasing = TRUE)
top_features  <- names(mean_abs_shap)[seq_len(min(top_n, length(mean_abs_shap)))]

shap_long <- rbindlist(lapply(top_features, function(f) {
  vals <- test_matrix[, f]
  rank_pct <- rank(vals, na.last = "keep", ties.method = "average") / sum(!is.na(vals))
  data.table(Feature = f, SHAP = contrib[, f], FeatureValuePct = rank_pct)
}))
shap_long[, Feature := factor(Feature, levels = rev(top_features))]

figure2 <- ggplot(shap_long, aes(SHAP, Feature, colour = FeatureValuePct)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_jitter(height = 0.3, alpha = 0.55, size = 0.7) +
  scale_colour_gradient(low = "#f0b429", high = "#2166ac", name = "Feature value\n(percentile)",
                         labels = c("Low", "High"), breaks = c(0, 1)) +
  labs(title = "Figure 2 - SHAP summary plot", x = "SHAP value (impact on predicted probability)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "Figure2_SHAP_Summary.png"), figure2, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "Figure2_SHAP_Summary.pdf"), figure2, width = 8, height = 6)

fwrite(data.table(Feature = names(mean_abs_shap), Mean_Abs_SHAP = mean_abs_shap),
       "Output Files/CSVs/XGBoost/FullCohort/shap_summary_top20.csv")
