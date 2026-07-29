################################################################################
# Evaluate Main XGBoost Model on Strict-Timing Cohort
#
# Applies the main model (trained on looser timing data, pre-2024) to the
# 2024-onwards patients who were processed through the stricter blood/vitals
# timing pipeline (Trial_BC_5 / Trial_BC_8).
#
# Main model:  Models/XGBoost/FullCohort/xgb_fullcohort.json
# Platt coefs: Models/XGBoost/FullCohort/xgb_fullcohort_platt_calibration.rds
# Test data:   Output Files/CSVs/trial_blood_cultures_with_other_micro.csv (>= 2024-01-01)
################################################################################

library(data.table)
library(xgboost)
library(pROC)
library(PRROC)
library(ggplot2)
library(mgcv)
library(dplyr)
library(tidyr)
library(writexl)

set.seed(42)

fig_dir <- "Output Files/Figures/Models/XGBoost/MainModel_StrictCohort_Eval"
csv_dir <- "Output Files/CSVs/XGBoost/MainModel_StrictCohort_Eval"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# LOAD DATA — strict-timing pipeline, 2024+ only
################################################################################

bc_data <- fread("Output Files/CSVs/trial_blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))

bc_data <- bc_data[CultivabilityCategory %in% c(
  "Gram-positive", "Gram-negative", "Polymicrobial",
  "Fungi", "Anaerobes", "Other",
  "Culture-negative", "Contaminant - single", "Contaminant - recurrent",
  "Polymicrobial - contaminant"
)]

bc_data[, CollectionDate := as.Date(CollectionDateTime)]
bc_data <- bc_data[CollectionDate >= as.Date("2024-01-01")]

bc_data[, Outcome := fifelse(CultivabilityCategory %in% c(
  "Gram-positive", "Gram-negative", "Polymicrobial",
  "Fungi", "Anaerobes", "Other"), 1L, 0L)]

################################################################################
# FEATURE ENGINEERING — identical to Do_XGBoost_FullCohort.r
################################################################################

bc_data[, Sex_num := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]

median_imd <- median(bc_data$SOCIAL_IMDPercentile, na.rm = TRUE)
bc_data[is.na(SOCIAL_IMDPercentile), SOCIAL_IMDPercentile := median_imd]

bc_data[, AdmissionType_num := fcase(
  AdmissionMethod == "Emergency", 1,
  AdmissionMethod == "Elective",  2,
  default = 0
)]

bc_data[, SourceType_num := fcase(
  CultureSourceCategory == "Peripheral", 1,
  CultureSourceCategory %in% c("Central Line", "PICC", "Port/Portacath", "Arterial", "Femoral"), 2,
  default = 0
)]

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

bc_data[, VS_AnyFever    := fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), 1, 0)]
bc_data[, VS_AnyHypotension := fifelse(VS_AnyHypotension == 1 & !is.na(VS_AnyHypotension), 1, 0)]
bc_data[is.na(LOS_BeforeCulture_days), LOS_BeforeCulture_days := 0]

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

missing_cols <- setdiff(feature_cols, names(bc_data))
if (length(missing_cols) > 0) stop("Missing feature columns: ", paste(missing_cols, collapse = ", "))

test_matrix <- as.matrix(bc_data[, ..feature_cols])
test_label  <- bc_data$Outcome

################################################################################
# LOAD MODEL AND CALIBRATION
################################################################################

model      <- xgb.load("Models/XGBoost/FullCohort/xgb_fullcohort.json")
platt_coef <- readRDS("Models/XGBoost/FullCohort/xgb_fullcohort_platt_calibration.rds")

# Raw XGBoost predictions
dtest   <- xgb.DMatrix(data = test_matrix)
raw_pred <- predict(model, dtest)

# Apply Platt calibration using saved intercept + slope
logit_raw   <- qlogis(pmax(pmin(raw_pred, 0.999), 0.001))
test_pred   <- plogis(platt_coef[1] + platt_coef[2] * logit_raw)

################################################################################
# PERFORMANCE METRICS
################################################################################

test_roc <- roc(test_label, test_pred, quiet = TRUE)

coords_opt  <- coords(test_roc, "best", ret = c("threshold", "sensitivity", "specificity"))
threshold   <- coords_opt$threshold
pred_class  <- ifelse(test_pred >= threshold, 1L, 0L)

TP <- sum(pred_class == 1 & test_label == 1)
TN <- sum(pred_class == 0 & test_label == 0)
FP <- sum(pred_class == 1 & test_label == 0)
FN <- sum(pred_class == 0 & test_label == 1)

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
ppv         <- TP / (TP + FP)
npv         <- TN / (TN + FN)
accuracy    <- (TP + TN) / (TP + TN + FP + FN)
f1          <- 2 * ppv * sensitivity / (ppv + sensitivity)
brier       <- mean((test_pred - test_label)^2)

pr_obj  <- pr.curve(scores.class0 = test_pred[test_label == 1],
                    scores.class1 = test_pred[test_label == 0],
                    curve = TRUE)

auroc_ci <- ci.auc(test_roc)

set.seed(42)
n_boot <- 1000
boot_auprcs <- vapply(seq_len(n_boot), function(b) {
  idx <- sample(length(test_label), replace = TRUE)
  bp  <- test_pred[idx]; bl <- test_label[idx]
  if (sum(bl) > 0 && sum(bl) < length(bl))
    pr.curve(scores.class0 = bp[bl == 1], scores.class1 = bp[bl == 0])$auc.integral
  else NA_real_
}, numeric(1))
auprc_ci <- quantile(boot_auprcs, c(0.025, 0.975), na.rm = TRUE)

model_summary <- data.frame(
  Metric = c("AUROC", "AUROC_CI_Lower", "AUROC_CI_Upper",
             "AUPRC", "AUPRC_CI_Lower", "AUPRC_CI_Upper",
             "Sensitivity", "Specificity", "PPV", "NPV",
             "F1", "Accuracy", "Brier",
             "N_Total", "N_Positive", "Prevalence_Pct",
             "Threshold"),
  Value  = c(round(as.numeric(auc(test_roc)), 3),
             round(auroc_ci[1], 3), round(auroc_ci[3], 3),
             round(pr_obj$auc.integral, 3),
             round(auprc_ci[1], 3), round(auprc_ci[2], 3),
             round(sensitivity, 3), round(specificity, 3),
             round(ppv, 3), round(npv, 3),
             round(f1, 3), round(accuracy, 3), round(brier, 4),
             nrow(bc_data), sum(test_label),
             round(100 * mean(test_label), 1),
             round(threshold, 4))
)
fwrite(model_summary, file.path(csv_dir, "model_summary.csv"))

################################################################################
# FIGURES
################################################################################

# 1. ROC Curve
roc_df <- data.frame(FPR = 1 - test_roc$specificities, TPR = test_roc$sensitivities)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "#2E86AB", linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = 0, ymax = TPR), fill = "#2E86AB", alpha = 0.1) +
  annotate("text", x = 0.6, y = 0.25,
           label = sprintf("AUROC = %.3f\n(95%% CI: %.3f-%.3f)",
                           auc(test_roc), auroc_ci[1], auroc_ci[3]),
           size = 5, fontface = "bold") +
  labs(title = "ROC Curve — Main Model on Strict-Timing Cohort (2024+)",
       subtitle = sprintf("N = %s | Prevalence = %.1f%%",
                          format(nrow(bc_data), big.mark = ","), 100 * mean(test_label)),
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 13)) +
  coord_equal()

ggsave(file.path(fig_dir, "01_roc_curve.png"), p_roc, width = 8, height = 8, dpi = 300)

# 2. Precision-Recall Curve
pr_df <- data.frame(Recall = pr_obj$curve[, 1], Precision = pr_obj$curve[, 2])

p_pr <- ggplot(pr_df, aes(x = Recall, y = Precision)) +
  geom_line(color = "#E94F37", linewidth = 1.2) +
  geom_hline(yintercept = mean(test_label), linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.7, y = 0.3,
           label = sprintf("AUPRC = %.3f\nBaseline = %.3f", pr_obj$auc.integral, mean(test_label)),
           size = 5, fontface = "bold") +
  labs(title = "Precision-Recall Curve — Main Model on Strict-Timing Cohort (2024+)",
       subtitle = sprintf("Positive rate: %.1f%%", 100 * mean(test_label)),
       x = "Recall (Sensitivity)", y = "Precision (PPV)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 13)) +
  ylim(0, 1)

ggsave(file.path(fig_dir, "02_precision_recall_curve.png"), p_pr, width = 8, height = 8, dpi = 300)

# 3. Calibration Plot
test_df <- data.frame(pred = test_pred, actual = test_label)
test_df$bin <- cut(test_df$pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)

cal_df <- test_df %>%
  group_by(bin) %>%
  summarise(predicted = mean(pred), observed = mean(actual), n = n(), .groups = "drop") %>%
  filter(!is.na(bin))

cal_bam <- gam(actual ~ s(pred, k = 10), data = test_df, family = binomial)
pred_grid <- data.frame(pred = seq(0.01, max(test_pred) * 0.95, length.out = 200))
bam_pred  <- predict(cal_bam, newdata = pred_grid, type = "response", se.fit = TRUE)
pred_grid$observed <- bam_pred$fit
pred_grid$lower    <- pmax(0, bam_pred$fit - 1.96 * bam_pred$se.fit)
pred_grid$upper    <- pmin(1, bam_pred$fit + 1.96 * bam_pred$se.fit)

p_cal <- ggplot() +
  geom_ribbon(data = pred_grid, aes(x = pred, ymin = lower, ymax = upper),
              fill = "#2E86AB", alpha = 0.2) +
  geom_line(data = pred_grid, aes(x = pred, y = observed),
            color = "#2E86AB", linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_point(data = cal_df, aes(x = predicted, y = observed, size = n),
             color = "#2E86AB", alpha = 0.6) +
  scale_size_continuous(range = c(2, 10), labels = scales::comma) +
  labs(title = "Calibration Plot — Main Model on Strict-Timing Cohort (2024+)",
       x = "Predicted Probability", y = "Observed Proportion", size = "N") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 13)) +
  coord_equal() + xlim(0, 1) + ylim(0, 1)

ggsave(file.path(fig_dir, "03_calibration_plot.png"), p_cal, width = 8, height = 8, dpi = 300)

# 4. Risk Density
risk_df <- data.frame(
  Prediction = test_pred,
  Outcome = factor(ifelse(test_label == 1, "Positive", "Negative"),
                   levels = c("Positive", "Negative"))
)

p_density <- ggplot(risk_df, aes(x = Prediction, fill = Outcome, color = Outcome)) +
  geom_density(alpha = 0.4, linewidth = 1.2) +
  scale_fill_manual(values  = c("Positive" = "#E94F37", "Negative" = "#2E86AB"),
                    labels  = c("Positive Culture", "Negative Culture")) +
  scale_color_manual(values = c("Positive" = "#E94F37", "Negative" = "#2E86AB"),
                     labels = c("Positive Culture", "Negative Culture")) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.1)) +
  labs(title = "Risk Distribution — Main Model on Strict-Timing Cohort (2024+)",
       subtitle = sprintf("AUROC = %.3f", auc(test_roc)),
       x = "Predicted Probability", y = "Density",
       fill = "Actual Outcome", color = "Actual Outcome") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 13), legend.position = "bottom")

ggsave(file.path(fig_dir, "04_risk_density.png"), p_density, width = 10, height = 7, dpi = 300)

################################################################################
# SAVE EXCEL SUMMARY
################################################################################

excel_dir <- "Output Files/Publication_Tables_Figures/Per_Model_Summaries"
dir.create(excel_dir, showWarnings = FALSE, recursive = TRUE)

cohort_info <- data.frame(
  Item  = c("Model", "Train cohort", "Test cohort",
            "Vital-sign timing (test)", "Blood-test timing (test)",
            "N cultures (test)", "N positive (test)", "Prevalence (%)"),
  Value = c("XGBoost FullCohort (xgb_fullcohort.json)",
            "Standard pipeline, < 2024-01-01",
            "Strict-timing pipeline (Trial_BC_5/BC_8), >= 2024-01-01",
            "-48 h to -5 min before culture",
            "-48 h to -2 h before culture (closest pre-culture value)",
            as.character(nrow(bc_data)),
            as.character(sum(test_label)),
            as.character(round(100 * mean(test_label), 1)))
)

write_xlsx(list(
  Model_Summary  = model_summary,
  Cohort_Info    = cohort_info
), file.path(excel_dir, "MainModel_StrictCohort_Eval.xlsx"))

fwrite(data.table(pred = test_pred, outcome = test_label),
       file.path(csv_dir, "test_predictions.csv"))
