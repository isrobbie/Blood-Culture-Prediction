################################################################################
# BAM Model - Full Cohort (43 Features)
# NO FILTERING - uses all available cultures
# Uses 44 consolidated features with median imputation
################################################################################

library(data.table)
library(mgcv)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(PRROC)
library(gridExtra)

set.seed(42)

# Create output folders
fig_dir <- "Output Files/Figures/Models/BAM/FullCohort"
csv_dir <- "Output Files/CSVs/BAM/FullCohort"
model_dir <- "Models/BAM/FullCohort"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# 1. LOAD AND PREPARE DATA
################################################################################

bc_data <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))

bc_data <- bc_data[CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial",
                                                  "Fungi", "Anaerobes", "Other",
                                                  "Culture-negative", "Contaminant - single", "Contaminant - recurrent",
                                                  "Polymicrobial - contaminant")]

bc_data[, Outcome := fifelse(CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial",
                                                            "Fungi", "Anaerobes", "Other"), 1L, 0L)]

# Feature engineering (CONSOLIDATED - same as XGBoost)

bc_data[, Sex_f := factor(fifelse(is.na(DEMOGRAPHICS_LinkedSex), "Unknown", DEMOGRAPHICS_LinkedSex),
                          levels = c("F", "M", "Unknown"))]
median_imd <- median(bc_data$SOCIAL_IMDPercentile, na.rm = TRUE)
bc_data[is.na(SOCIAL_IMDPercentile), SOCIAL_IMDPercentile := median_imd]

# Admission Type: consolidated factor (Other, Emergency, Elective, Unknown)
bc_data[, AdmissionType_f := factor(fcase(
  is.na(AdmissionMethod), "Unknown",
  AdmissionMethod == "Emergency", "Emergency",
  AdmissionMethod == "Elective", "Elective",
  default = "Other"
), levels = c("Other", "Emergency", "Elective", "Unknown"))]

# Culture Source Type: consolidated factor (Other, Peripheral, Central, Unknown)
bc_data[, SourceType_f := factor(fcase(
  is.na(CultureSourceCategory), "Unknown",
  CultureSourceCategory == "Peripheral", "Peripheral",
  CultureSourceCategory %in% c("Central Line", "PICC", "Port/Portacath", "Arterial", "Femoral"), "Central",
  default = "Other"
), levels = c("Other", "Peripheral", "Central", "Unknown"))]

bc_data[is.na(PositiveCulturesInPrior28Days), PositiveCulturesInPrior28Days := 0]

# Infection source (top 4) - optional, included when available
for (col in c("Source_Urinary", "Source_Respiratory", "Source_Skin_Soft_Tissue", "Source_Other")) {
  bc_data[is.na(get(col)), (col) := 0]
  bc_data[, paste0(col, "_f") := factor(fifelse(get(col) == 1, "Yes", "No"))]
}

bc_data[, HasInfectionSource := (Source_Urinary == 1 | Source_Respiratory == 1 | 
                                  Source_Skin_Soft_Tissue == 1 | Source_Other == 1)]

# Antibiotics (top 4)
for (col in c("HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before", "Carbapenem_Before")) {
  bc_data[is.na(get(col)), (col) := 0]
  bc_data[, paste0(col, "_f") := factor(fifelse(get(col) == 1, "Yes", "No"))]
}

# Prior micro (top 2)
for (col in c("Any_Prior_Other_Pos", "Prior_Urine_Pos")) {
  bc_data[is.na(get(col)), (col) := 0]
  bc_data[, paste0(col, "_f") := factor(fifelse(get(col) == 1, "Yes", "No"))]
}

bc_data[, VS_AnyFever_f := factor(fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), "Yes", "No"))]
bc_data[, VS_AnyHypotension_f := factor(fifelse(VS_AnyHypotension == 1 & !is.na(VS_AnyHypotension), "Yes", "No"))]

bc_data[is.na(LOS_BeforeCulture_days), LOS_BeforeCulture_days := 0]

# Create missingness indicators BEFORE imputation (for deployment readiness analysis)
bc_data[, VS_Complete_orig := !is.na(VS_Baseline_Temp) & !is.na(VS_Baseline_HR) & 
                              !is.na(VS_Baseline_SBP) & !is.na(VS_Baseline_RR)]
bc_data[, FBC_Complete_orig := !is.na(BT_Closest_WCC) & !is.na(BT_Closest_Neut) & 
                               !is.na(BT_Closest_Lymph) & !is.na(BT_Closest_Hb) & !is.na(BT_Closest_Plt)]
bc_data[, UE_Complete_orig := !is.na(BT_Closest_Na) & !is.na(BT_Closest_Creat) & !is.na(BT_Closest_Urea)]
bc_data[, LFT_Complete_orig := !is.na(BT_Closest_Bili) & !is.na(BT_Closest_ALP) & !is.na(BT_Closest_Alb)]
bc_data[, CRP_Complete_orig := !is.na(BT_Closest_CRP)]
bc_data[, BG_Complete_orig := !is.na(BG_Closest_Lactate) & !is.na(BG_Closest_pH)]
bc_data[, RoutineBloods_Complete_orig := FBC_Complete_orig & UE_Complete_orig & LFT_Complete_orig & CRP_Complete_orig]

# Impute continuous variables (median imputation)
cont_vars <- c("VS_Baseline_HR", "VS_Baseline_Temp", "VS_Baseline_SBP", "VS_Baseline_RR", "VS_Baseline_NEWS",
               "BT_Closest_CRP", "BT_Closest_WCC", "BT_Closest_Neut", "BT_Closest_Lymph",
               "BT_Closest_Alb", "BT_Closest_Plt", "BT_Closest_Creat", "BT_Closest_Bili",
               "BT_Closest_Urea", "BT_Closest_Hb", "BT_Closest_Na",
               "BT_Closest_Eosin", "BT_Closest_Mono", "BT_Closest_ALP", "BT_Closest_HCT",
               "BG_Closest_Lactate", "BG_Closest_pH", "BG_Closest_pCO2",
               "LOS_BeforeCulture_days")

for (col in cont_vars) {
  med <- median(bc_data[[col]], na.rm = TRUE)
  bc_data[is.na(get(col)), (col) := med]
}

################################################################################
# 2. TRAIN/TEST SPLIT
################################################################################

bc_data[, CollectionDate := as.Date(CollectionDateTime)]
split_date <- as.Date("2024-01-01")
train_data_dt <- bc_data[CollectionDate <  split_date]
test_data_dt  <- bc_data[CollectionDate >= split_date]

train_data <- as.data.frame(train_data_dt)
test_data  <- as.data.frame(test_data_dt)

# Class weights
n_pos <- sum(train_data$Outcome == 1)
n_neg <- sum(train_data$Outcome == 0)
train_data$class_weight <- ifelse(train_data$Outcome == 1, n_neg / n_pos, 1)

################################################################################
# 3. CROSS-VALIDATION
################################################################################

# Create folds by patient
patient_list <- unique(train_data_dt[, .(ClusterID, Outcome)])
patient_list[, ever_positive := max(Outcome), by = ClusterID]
patient_unique <- unique(patient_list[, .(ClusterID, ever_positive)])
pos_pts <- patient_unique[ever_positive == 1]$ClusterID
neg_pts <- patient_unique[ever_positive == 0]$ClusterID

set.seed(42)
pos_folds <- sample(rep(1:5, length.out = length(pos_pts)))
neg_folds <- sample(rep(1:5, length.out = length(neg_pts)))
patient_fold_map <- data.table(
  ClusterID = c(pos_pts, neg_pts),
  fold = c(pos_folds, neg_folds)
)

cv_aucs  <- numeric(5)
cv_auprcs <- numeric(5)
train_data_dt[, pred_oof := NA_real_]  # Store out-of-fold predictions

for (fold in 1:5) {
  val_pts <- patient_fold_map[fold == !!fold]$ClusterID
  fold_train <- as.data.frame(train_data_dt[!ClusterID %in% val_pts])
  fold_val <- as.data.frame(train_data_dt[ClusterID %in% val_pts])
  
  # Class weights
  fold_train$class_weight <- ifelse(fold_train$Outcome == 1, sum(fold_train$Outcome == 0) / sum(fold_train$Outcome == 1), 1)
  
  # Train model
  fold_model <- bam(
    Outcome ~ 
      s(AgeAtCulture, k = 4) + Sex_f +
      AdmissionType_f + s(LOS_BeforeCulture_days, k = 4) +
      SourceType_f +
      s(PositiveCulturesInPrior28Days, k = 4) +
      Source_Respiratory_f + Source_Other_f +
      HasAntibiotics_Before_f + BroadSpectrum_Before_f + Ceph3rdGenPlus_Before_f + Carbapenem_Before_f +
      Any_Prior_Other_Pos_f + Prior_Urine_Pos_f +
      s(VS_Baseline_Temp, k = 8) + s(VS_Baseline_SBP, k = 4) +
      s(VS_Baseline_NEWS, k = 4) + VS_AnyFever_f +
      s(BT_Closest_CRP, k = 8) + s(BT_Closest_WCC, k = 4) + s(BT_Closest_Neut, k = 4) + s(BT_Closest_Lymph, k = 4) +
      s(BT_Closest_Alb, k = 4) + s(BT_Closest_Plt, k = 4) + s(BT_Closest_Creat, k = 4) + s(BT_Closest_Bili, k = 4) +
      s(BT_Closest_Na, k = 4) +
      s(BT_Closest_Eosin, k = 4) + s(BT_Closest_Mono, k = 4) + s(BT_Closest_ALP, k = 4) +
      s(BG_Closest_Lactate, k = 8) + s(BG_Closest_pH, k = 4),
    data = fold_train, family = binomial(link = "logit"),
    weights = class_weight, discrete = TRUE, select = TRUE
  )
  
  # Evaluate and store out-of-fold predictions
  fold_pred <- predict(fold_model, newdata = fold_val, type = "response")
  fold_auc <- as.numeric(auc(roc(fold_val$Outcome, fold_pred, quiet = TRUE)))
  cv_aucs[fold] <- fold_auc
  valid <- !is.na(fold_pred)
  cv_auprcs[fold] <- tryCatch(
    pr.curve(
      scores.class0 = fold_pred[valid & fold_val$Outcome == 1],
      scores.class1 = fold_pred[valid & fold_val$Outcome == 0]
    )$auc.integral,
    error = function(e) NA_real_
  )
  
  # Store out-of-fold predictions for calibration
  train_data_dt[ClusterID %in% val_pts, pred_oof := fold_pred]
}

cv_results <- data.frame(
  Fold = 1:5, AUROC = cv_aucs, AUPRC = cv_auprcs,
  Mean_AUROC = mean(cv_aucs), SD_AUROC = sd(cv_aucs)
)

################################################################################
# 4. FIT FINAL BAM MODEL (simplified to match XGBoost features)
################################################################################

bam_model <- bam(
  Outcome ~ 
    s(AgeAtCulture, k = 4) + Sex_f +
    AdmissionType_f + s(LOS_BeforeCulture_days, k = 4) +
    SourceType_f +
    s(PositiveCulturesInPrior28Days, k = 4) +
    Source_Respiratory_f + Source_Other_f +
    HasAntibiotics_Before_f + BroadSpectrum_Before_f + Ceph3rdGenPlus_Before_f + Carbapenem_Before_f +
    Any_Prior_Other_Pos_f + Prior_Urine_Pos_f +
    s(VS_Baseline_Temp, k = 8) + s(VS_Baseline_SBP, k = 4) +
    s(VS_Baseline_NEWS, k = 4) + VS_AnyFever_f +
    s(BT_Closest_CRP, k = 8) + s(BT_Closest_WCC, k = 4) + s(BT_Closest_Neut, k = 4) + s(BT_Closest_Lymph, k = 4) +
    s(BT_Closest_Alb, k = 4) + s(BT_Closest_Plt, k = 4) + s(BT_Closest_Creat, k = 4) + s(BT_Closest_Bili, k = 4) +
    s(BT_Closest_Na, k = 4) +
    s(BT_Closest_Eosin, k = 4) + s(BT_Closest_Mono, k = 4) + s(BT_Closest_ALP, k = 4) +
    s(BG_Closest_Lactate, k = 8) + s(BG_Closest_pH, k = 4),
  data = train_data, family = binomial(link = "logit"),
  weights = class_weight, discrete = TRUE, select = TRUE
)

################################################################################
# 5. EVALUATE MODEL
################################################################################

train_pred_raw <- train_data_dt$pred_oof  # Use out-of-fold predictions!
test_pred_raw <- predict(bam_model, newdata = test_data, type = "response")

################################################################################
# PLATT SCALING CALIBRATION (fitted on out-of-fold predictions)
################################################################################

# Fit Platt model on OUT-OF-FOLD predictions (unbiased!)
train_pred_logit <- qlogis(pmax(pmin(train_pred_raw, 0.999), 0.001))
platt_model <- glm(train_data$Outcome ~ train_pred_logit, family = binomial)

test_pred_logit <- qlogis(pmax(pmin(test_pred_raw, 0.999), 0.001))
test_pred <- predict(platt_model, newdata = data.frame(train_pred_logit = test_pred_logit), type = "response")
train_pred <- predict(platt_model, newdata = data.frame(train_pred_logit = train_pred_logit), type = "response")

brier_raw <- mean((test_pred_raw - test_data$Outcome)^2)
brier_calibrated <- mean((test_pred - test_data$Outcome)^2)

train_roc <- roc(train_data$Outcome, train_pred, quiet = TRUE)
test_roc <- roc(test_data$Outcome, test_pred, quiet = TRUE)

coords_opt <- coords(test_roc, "best", ret = c("threshold", "sensitivity", "specificity"))
threshold <- coords_opt$threshold

test_pred_class <- ifelse(test_pred >= threshold, 1, 0)
TP <- sum(test_pred_class == 1 & test_data$Outcome == 1)
TN <- sum(test_pred_class == 0 & test_data$Outcome == 0)
FP <- sum(test_pred_class == 1 & test_data$Outcome == 0)
FN <- sum(test_pred_class == 0 & test_data$Outcome == 1)

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
ppv <- TP / (TP + FP)
npv <- TN / (TN + FN)
accuracy <- (TP + TN) / (TP + TN + FP + FN)
f1 <- 2 * ppv * sensitivity / (ppv + sensitivity)

################################################################################
# 6. GENERATE FIGURES (same as Filtered cohort)
################################################################################

# 1. Parametric Coefficients (Factor variables - Odds Ratios)
param_coef <- summary(bam_model)$p.table
param_df <- data.frame(
  Variable = rownames(param_coef),
  Estimate = param_coef[, "Estimate"],
  SE = param_coef[, "Std. Error"],
  OR = exp(param_coef[, "Estimate"]),
  OR_lower = exp(param_coef[, "Estimate"] - 1.96 * param_coef[, "Std. Error"]),
  OR_upper = exp(param_coef[, "Estimate"] + 1.96 * param_coef[, "Std. Error"]),
  P_value = param_coef[, "Pr(>|z|)"],
  row.names = NULL
)
param_df <- param_df[param_df$Variable != "(Intercept)", ]
param_df$Sig <- ifelse(param_df$P_value < 0.001, "***", 
                       ifelse(param_df$P_value < 0.01, "**", 
                              ifelse(param_df$P_value < 0.05, "*", "")))

p_or <- ggplot(param_df, aes(x = reorder(Variable, OR), y = OR)) +
  geom_point(size = 3, color = "#E94F37") +
  geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.2, color = "#E94F37") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  coord_flip() +
  labs(title = "BAM - Odds Ratios for Categorical Variables",
       subtitle = sprintf("Test AUROC = %.3f", auc(test_roc)),
       x = "", y = "Odds Ratio (95% CI)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "01_odds_ratios.png"), p_or, width = 10, height = 8, dpi = 300)

# 2. Smooth Terms Significance
smooth_coef <- summary(bam_model)$s.table
smooth_df <- data.frame(
  Variable = gsub("s\\(|\\)", "", rownames(smooth_coef)),
  edf = smooth_coef[, "edf"],
  Chi_sq = smooth_coef[, "Chi.sq"],
  P_value = smooth_coef[, "p-value"],
  row.names = NULL
)
smooth_df$Sig <- ifelse(smooth_df$P_value < 0.001, "***", 
                        ifelse(smooth_df$P_value < 0.01, "**", 
                               ifelse(smooth_df$P_value < 0.05, "*", "ns")))

p_smooth_sig <- ggplot(smooth_df, aes(x = reorder(Variable, Chi_sq), y = Chi_sq, fill = Sig)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("***" = "#2E86AB", "**" = "#5DA5C7", "*" = "#8BC4E0", "ns" = "#D3D3D3")) +
  coord_flip() +
  labs(title = "BAM - Smooth Term Significance (Chi-squared)",
       x = "", y = "Chi-squared Statistic") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "none")

ggsave(file.path(fig_dir, "02_smooth_term_significance.png"), p_smooth_sig, width = 10, height = 10, dpi = 300)

# 3. Smooth Term Plots (Top 12)
smooth_df_ordered <- smooth_df[order(-smooth_df$Chi_sq), ]
top_smooths <- head(smooth_df_ordered$Variable, 12)

# Get all smooth term names from the model
all_smooths <- names(bam_model$smooth)

png(file.path(fig_dir, "03_smooth_term_effects.png"), width = 2400, height = 2000, res = 150)
par(mfrow = c(4, 3), mar = c(4, 4, 3, 1))
for (var_name in top_smooths) {
  # Find which smooth index this variable corresponds to
  smooth_idx <- which(sapply(bam_model$smooth, function(s) s$label == paste0("s(", var_name, ")")))
  
  if (length(smooth_idx) > 0) {
    plot(bam_model, select = smooth_idx, shade = TRUE, 
         main = var_name, xlab = var_name, ylab = "Effect on log-odds",
         ylim = c(-2, 2), cex.main = 1.2)
  }
}
dev.off()

# 4. ROC Curve
roc_df <- data.frame(FPR = 1 - test_roc$specificities, TPR = test_roc$sensitivities)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "#E94F37", linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = 0, ymax = TPR), fill = "#E94F37", alpha = 0.1) +
  annotate("text", x = 0.6, y = 0.3, 
           label = sprintf("AUROC = %.3f\n(95%% CI: %.3f-%.3f)", 
                          auc(test_roc), ci.auc(test_roc)[1], ci.auc(test_roc)[3]),
           size = 5, fontface = "bold") +
  labs(title = "ROC Curve - BAM",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14)) +
  coord_equal()

ggsave(file.path(fig_dir, "04_roc_curve.png"), p_roc, width = 8, height = 8, dpi = 300)

# 5. Precision-Recall Curve
pr_obj <- pr.curve(scores.class0 = test_pred[test_data$Outcome == 1], 
                   scores.class1 = test_pred[test_data$Outcome == 0], 
                   curve = TRUE)

pr_df <- data.frame(Recall = pr_obj$curve[, 1], Precision = pr_obj$curve[, 2])

p_pr <- ggplot(pr_df, aes(x = Recall, y = Precision)) +
  geom_line(color = "#E94F37", linewidth = 1.2) +
  geom_hline(yintercept = mean(test_data$Outcome), linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.7, y = 0.3, 
           label = sprintf("AUPRC = %.3f\nBaseline = %.3f", pr_obj$auc.integral, mean(test_data$Outcome)),
           size = 5, fontface = "bold") +
  labs(title = "Precision-Recall Curve - BAM",
       subtitle = sprintf("Positive rate: %.1f%%", 100 * mean(test_data$Outcome)),
       x = "Recall (Sensitivity)", y = "Precision (PPV)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14)) +
  ylim(0, 1)

ggsave(file.path(fig_dir, "05_precision_recall_curve.png"), p_pr, width = 8, height = 8, dpi = 300)

# 6. Calibration Plot with 95% Confidence Band (BAM-based)
test_df_cal <- data.frame(pred = test_pred, actual = test_data$Outcome)
test_df_cal$bin <- cut(test_df_cal$pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)

cal_df <- test_df_cal %>%
  group_by(bin) %>%
  summarise(predicted = mean(pred), observed = mean(actual), n = n()) %>%
  filter(!is.na(bin))

cal_bam <- mgcv::gam(actual ~ s(pred, k = 10), data = test_df_cal, family = binomial)
pred_grid <- data.frame(pred = seq(0.01, max(test_pred) * 0.95, length.out = 200))
bam_pred <- predict(cal_bam, newdata = pred_grid, type = "response", se.fit = TRUE)

pred_grid$observed <- bam_pred$fit
pred_grid$lower <- pmax(0, bam_pred$fit - 1.96 * bam_pred$se.fit)
pred_grid$upper <- pmin(1, bam_pred$fit + 1.96 * bam_pred$se.fit)

p_cal <- ggplot() +
  geom_ribbon(data = pred_grid, aes(x = pred, ymin = lower, ymax = upper), 
              fill = "#E94F37", alpha = 0.2) +
  geom_line(data = pred_grid, aes(x = pred, y = observed), 
            color = "#E94F37", linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_point(data = cal_df, aes(x = predicted, y = observed, size = n), 
             color = "#E94F37", alpha = 0.6) +
  scale_size_continuous(range = c(2, 10), labels = scales::comma) +
  labs(title = "Calibration Plot with 95% Confidence Band - BAM",
       x = "Predicted Probability", y = "Observed Proportion", size = "N") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14)) +
  coord_equal() + xlim(0, 1) + ylim(0, 1)

ggsave(file.path(fig_dir, "06_calibration_plot.png"), p_cal, width = 8, height = 8, dpi = 300)

# 7. Confusion Matrix
cm_df <- data.frame(
  Predicted = factor(c("Negative", "Negative", "Positive", "Positive"), levels = c("Positive", "Negative")),
  Actual = factor(c("Negative", "Positive", "Negative", "Positive"), levels = c("Positive", "Negative")),
  Count = c(TN, FN, FP, TP),
  Label = c(paste0("TN\n", format(TN, big.mark = ",")), 
            paste0("FN\n", format(FN, big.mark = ",")),
            paste0("FP\n", format(FP, big.mark = ",")), 
            paste0("TP\n", format(TP, big.mark = ",")))
)

p_cm <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Count)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = Label), size = 6, fontface = "bold") +
  scale_fill_gradient(low = "#f0f0f0", high = "#E94F37") +
  labs(title = "Confusion Matrix - BAM",
       subtitle = sprintf("Youden threshold = %.3f (one potential operating point) | Accuracy = %.1f%%", threshold, 100 * accuracy),
       x = "Actual", y = "Predicted") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "none",
        axis.text = element_text(size = 12))

ggsave(file.path(fig_dir, "07_confusion_matrix.png"), p_cm, width = 7, height = 6, dpi = 300)

# 8. Performance Metrics Summary
metrics_df <- data.frame(
  Metric = c("AUROC", "AUPRC", "Sensitivity", "Specificity", "PPV", "NPV", "F1 Score", "Accuracy"),
  Value = c(auc(test_roc), pr_obj$auc.integral, sensitivity, specificity, ppv, npv, f1, accuracy),
  Display = c(sprintf("%.3f", auc(test_roc)), sprintf("%.3f", pr_obj$auc.integral),
              sprintf("%.1f%%", 100*sensitivity), sprintf("%.1f%%", 100*specificity),
              sprintf("%.1f%%", 100*ppv), sprintf("%.1f%%", 100*npv),
              sprintf("%.3f", f1), sprintf("%.1f%%", 100*accuracy))
)

p_metrics <- ggplot(metrics_df, aes(x = reorder(Metric, Value), y = Value)) +
  geom_bar(stat = "identity", fill = "#E94F37", width = 0.6) +
  geom_text(aes(label = Display), hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2)) +
  labs(title = "Model Performance Metrics - BAM",
       subtitle = sprintf("Test Set (N = %s) | Threshold-dependent metrics at Youden optimum", format(nrow(test_data), big.mark = ",")),
       x = "", y = "Value") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "08_performance_metrics.png"), p_metrics, width = 9, height = 6, dpi = 300)

# 9. Decision Curve Analysis - Net Benefit & Absolute Cost
test_label <- test_data$Outcome
dca_thresholds <- seq(0.01, 0.50, by = 0.01)
n_test <- length(test_label)
prevalence <- mean(test_label)

dca_results <- lapply(dca_thresholds, function(pt) {
  pred_pos <- test_pred >= pt
  tp <- sum(pred_pos & test_label == 1)
  fp <- sum(pred_pos & test_label == 0)
  fn <- sum(!pred_pos & test_label == 1)
  
  nb_model <- (tp / n_test) - (fp / n_test) * (pt / (1 - pt))
  nb_treat_all <- prevalence - (1 - prevalence) * (pt / (1 - pt))
  nb_treat_none <- 0
  
  cost_model <- (fp / n_test) * (pt / (1 - pt)) + (fn / n_test)
  cost_treat_all <- (1 - prevalence) * (pt / (1 - pt))
  cost_treat_none <- prevalence
  
  data.frame(Threshold = pt, NB_Model = nb_model, NB_TreatAll = nb_treat_all, NB_TreatNone = nb_treat_none,
             Cost_Model = cost_model, Cost_TreatAll = cost_treat_all, Cost_TreatNone = cost_treat_none)
})
dca_df <- do.call(rbind, dca_results)

dca_long <- pivot_longer(dca_df[, c("Threshold", "NB_Model", "NB_TreatAll", "NB_TreatNone")],
                         cols = starts_with("NB_"), names_to = "Strategy", values_to = "NetBenefit")
dca_long$Strategy <- factor(gsub("NB_", "", dca_long$Strategy), levels = c("Model", "TreatAll", "TreatNone"))

p_dca_nb <- ggplot(dca_long, aes(x = Threshold, y = NetBenefit, color = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
  scale_color_manual(values = c("Model" = "#E94F37", "TreatAll" = "#2E86AB", "TreatNone" = "gray50"),
                     labels = c("BAM Model", "Treat All", "Treat None")) +
  scale_linetype_manual(values = c("Model" = "solid", "TreatAll" = "dashed", "TreatNone" = "dotted"),
                        labels = c("BAM Model", "Treat All", "Treat None")) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.5, 0.1)) +
  labs(title = "Decision Curve Analysis - Net Benefit",
       subtitle = "Net benefit across threshold probabilities",
       x = "Threshold Probability", y = "Net Benefit", color = "Strategy", linetype = "Strategy") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom") +
  coord_cartesian(ylim = c(-0.05, max(dca_df$NB_Model) * 1.1))

ggsave(file.path(fig_dir, "09_decision_curve_net_benefit.png"), p_dca_nb, width = 9, height = 7, dpi = 300)

cost_long <- pivot_longer(dca_df[, c("Threshold", "Cost_Model", "Cost_TreatAll", "Cost_TreatNone")],
                          cols = starts_with("Cost_"), names_to = "Strategy", values_to = "Cost")
cost_long$Strategy <- factor(gsub("Cost_", "", cost_long$Strategy), levels = c("Model", "TreatAll", "TreatNone"))

p_dca_cost <- ggplot(cost_long, aes(x = Threshold, y = Cost, color = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Model" = "#E94F37", "TreatAll" = "#2E86AB", "TreatNone" = "gray50"),
                     labels = c("BAM Model", "Treat All", "Treat None")) +
  scale_linetype_manual(values = c("Model" = "solid", "TreatAll" = "dashed", "TreatNone" = "dotted"),
                        labels = c("BAM Model", "Treat All", "Treat None")) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 0.5, 0.1)) +
  labs(title = "Decision Curve Analysis - Absolute Cost",
       subtitle = "Weighted cost of errors across threshold probabilities",
       x = "Threshold Probability", y = "Normalized Cost", color = "Strategy", linetype = "Strategy") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")

ggsave(file.path(fig_dir, "10_decision_curve_cost.png"), p_dca_cost, width = 9, height = 7, dpi = 300)

fwrite(dca_df, file.path(csv_dir, "decision_curve_analysis.csv"))

# 10. Risk Distribution Plot
risk_dist_df <- data.frame(
  Prediction = test_pred,
  Outcome = factor(ifelse(test_label == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))
)

p_risk_dist <- ggplot(risk_dist_df, aes(x = Prediction, fill = Outcome)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 50, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Positive" = "#E94F37", "Negative" = "#2E86AB"),
                    labels = c("Positive Culture", "Negative Culture")) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.1)) +
  labs(title = "Risk Distribution - Predicted Probability by Outcome",
       subtitle = sprintf("Test Set: %s positive, %s negative | AUROC = %.3f",
                         format(sum(test_label == 1), big.mark = ","),
                         format(sum(test_label == 0), big.mark = ","), auc(test_roc)),
       x = "Predicted Probability", y = "Count", fill = "Actual Outcome") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")

ggsave(file.path(fig_dir, "11_risk_distribution.png"), p_risk_dist, width = 10, height = 7, dpi = 300)

p_risk_density <- ggplot(risk_dist_df, aes(x = Prediction, fill = Outcome, color = Outcome)) +
  geom_density(alpha = 0.4, linewidth = 1.2) +
  scale_fill_manual(values = c("Positive" = "#E94F37", "Negative" = "#2E86AB"),
                    labels = c("Positive Culture", "Negative Culture")) +
  scale_color_manual(values = c("Positive" = "#E94F37", "Negative" = "#2E86AB"),
                     labels = c("Positive Culture", "Negative Culture")) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.1)) +
  labs(title = "Risk Distribution - Probability Density by Outcome",
       subtitle = sprintf("Separation between positive and negative cultures | AUROC = %.3f", auc(test_roc)),
       x = "Predicted Probability", y = "Density", fill = "Actual Outcome", color = "Actual Outcome") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")

ggsave(file.path(fig_dir, "12_risk_density.png"), p_risk_density, width = 10, height = 7, dpi = 300)

# Save model summary
model_summary <- data.frame(
  Metric = c("AUROC", "AUPRC", "Sensitivity", "Specificity", "PPV", "NPV", "F1", "Train-Test Gap", "Brier_Raw", "Brier_Calibrated"),
  Value = c(round(auc(test_roc), 3), round(pr_obj$auc.integral, 3),
            round(sensitivity, 3), round(specificity, 3),
            round(ppv, 3), round(npv, 3), round(f1, 3),
            round(auc(train_roc) - auc(test_roc), 3),
            round(brier_raw, 4), round(brier_calibrated, 4))
)

# Save outputs
saveRDS(bam_model, file.path(model_dir, "bam_fullcohort_model.rds"))
saveRDS(platt_model, file.path(model_dir, "bam_fullcohort_platt_calibration.rds"))
fwrite(param_df, file.path(csv_dir, "parametric_coefficients.csv"))
fwrite(smooth_df, file.path(csv_dir, "smooth_term_statistics.csv"))
fwrite(cv_results, file.path(csv_dir, "cross_validation_results.csv"))

# Save model summary as Excel
library(writexl)
excel_dir <- "Output Files/Publication_Tables_Figures/Per_Model_Summaries"
dir.create(excel_dir, showWarnings = FALSE, recursive = TRUE)
write_xlsx(list(
  Model_Summary = model_summary,
  CV_Results = cv_results
), file.path(excel_dir, "BAM_FullCohort_Summary.xlsx"))

# Calibration and cohort summaries
calib_summary <- data.frame(
  Metric = c("Brier_Raw", "Brier_Calibrated", "Brier_Improvement_Pct"),
  Value = c(brier_raw, brier_calibrated, 100 * (brier_raw - brier_calibrated) / brier_raw)
)
fwrite(calib_summary, file.path(csv_dir, "calibration_summary.csv"))

# Cohort summary (with date range for TRIPOD 5b)
cohort_summary <- data.frame(
  Metric = c("Total_Cultures", "Total_Patients", "Train_N", "Test_N", "Positive_Rate", 
             "Infection_Source_Pct", "Data_Start_Date", "Data_End_Date"),
  Value = c(nrow(bc_data), uniqueN(bc_data$ClusterID), nrow(train_data), nrow(test_data),
            100 * mean(bc_data$Outcome), 100 * mean(bc_data$HasInfectionSource),
            as.character(as.Date(min(bc_data$CollectionDateTime, na.rm = TRUE))),
            as.character(as.Date(max(bc_data$CollectionDateTime, na.rm = TRUE))))
)
fwrite(cohort_summary, file.path(csv_dir, "cohort_summary.csv"))

################################################################################
# TRIPOD FAIRNESS: Subgroup Performance Analysis (Items 3c, 14, 23a)
################################################################################

test_data_eval <- as.data.table(test_data)
test_data_eval[, pred := test_pred]
test_data_eval[, pred_class := fifelse(pred >= threshold, 1L, 0L)]

test_data_eval[, AgeGroup := fcase(
  AgeAtCulture < 18, "Pediatric (<18)",
  AgeAtCulture < 65, "Adult (18-64)",
  AgeAtCulture >= 65, "Elderly (65+)",
  default = NA_character_
)]

test_data_eval[, IMD_Quintile := fcase(
  SOCIAL_IMDPercentile <= 20, "Q1 (Most Deprived)",
  SOCIAL_IMDPercentile <= 40, "Q2",
  SOCIAL_IMDPercentile <= 60, "Q3",
  SOCIAL_IMDPercentile <= 80, "Q4",
  SOCIAL_IMDPercentile <= 100, "Q5 (Least Deprived)",
  default = NA_character_
)]

test_data_eval[, PatientSetting := fcase(
  CultureLocation == "ED", "A&E",
  CultureLocation == "Inpatient Ward", "Inpatient",
  default = "Other"
)]

calc_subgroup_metrics <- function(dt, subgroup_col, subgroup_val) {
  sub_dt <- dt[get(subgroup_col) == subgroup_val]
  if (nrow(sub_dt) < 50 || sum(sub_dt$Outcome) < 10) {
    return(data.table(
      Subgroup_Type = subgroup_col, Subgroup = subgroup_val,
      N = nrow(sub_dt), N_Positive = sum(sub_dt$Outcome),
      Prevalence = NA_real_, AUROC = NA_real_, Sensitivity = NA_real_,
      Specificity = NA_real_, PPV = NA_real_, NPV = NA_real_
    ))
  }
  sub_roc <- tryCatch(roc(sub_dt$Outcome, sub_dt$pred, quiet = TRUE), error = function(e) NULL)
  if (is.null(sub_roc)) {
    return(data.table(
      Subgroup_Type = subgroup_col, Subgroup = subgroup_val,
      N = nrow(sub_dt), N_Positive = sum(sub_dt$Outcome),
      Prevalence = 100 * mean(sub_dt$Outcome), AUROC = NA_real_,
      Sensitivity = NA_real_, Specificity = NA_real_, PPV = NA_real_, NPV = NA_real_
    ))
  }
  TP <- sum(sub_dt$pred_class == 1 & sub_dt$Outcome == 1)
  TN <- sum(sub_dt$pred_class == 0 & sub_dt$Outcome == 0)
  FP <- sum(sub_dt$pred_class == 1 & sub_dt$Outcome == 0)
  FN <- sum(sub_dt$pred_class == 0 & sub_dt$Outcome == 1)
  data.table(
    Subgroup_Type = subgroup_col, Subgroup = subgroup_val,
    N = nrow(sub_dt), N_Positive = sum(sub_dt$Outcome),
    Prevalence = round(100 * mean(sub_dt$Outcome), 1),
    AUROC = round(as.numeric(auc(sub_roc)), 3),
    Sensitivity = round(TP / (TP + FN), 3),
    Specificity = round(TN / (TN + FP), 3),
    PPV = round(TP / (TP + FP), 3),
    NPV = round(TN / (TN + FN), 3)
  )
}

subgroup_results <- rbindlist(list(
  data.table(Subgroup_Type = "Overall", Subgroup = "All",
             N = nrow(test_data_eval), N_Positive = sum(test_data_eval$Outcome),
             Prevalence = round(100 * mean(test_data_eval$Outcome), 1),
             AUROC = round(auc(test_roc), 3),
             Sensitivity = round(sensitivity, 3), Specificity = round(specificity, 3),
             PPV = round(ppv, 3), NPV = round(npv, 3)),
  calc_subgroup_metrics(test_data_eval, "AgeGroup", "Pediatric (<18)"),
  calc_subgroup_metrics(test_data_eval, "AgeGroup", "Adult (18-64)"),
  calc_subgroup_metrics(test_data_eval, "AgeGroup", "Elderly (65+)"),
  calc_subgroup_metrics(test_data_eval, "DEMOGRAPHICS_LinkedSex", "M"),
  calc_subgroup_metrics(test_data_eval, "DEMOGRAPHICS_LinkedSex", "F"),
  calc_subgroup_metrics(test_data_eval, "PatientSetting", "A&E"),
  calc_subgroup_metrics(test_data_eval, "PatientSetting", "Inpatient"),
  calc_subgroup_metrics(test_data_eval, "PatientSetting", "Other"),
  calc_subgroup_metrics(test_data_eval, "IMD_Quintile", "Q1 (Most Deprived)"),
  calc_subgroup_metrics(test_data_eval, "IMD_Quintile", "Q2"),
  calc_subgroup_metrics(test_data_eval, "IMD_Quintile", "Q3"),
  calc_subgroup_metrics(test_data_eval, "IMD_Quintile", "Q4"),
  calc_subgroup_metrics(test_data_eval, "IMD_Quintile", "Q5 (Least Deprived)")
))

if ("SOCIAL_Ethnicity" %in% names(test_data_eval)) {
  ethnicity_vals <- unique(test_data_eval$SOCIAL_Ethnicity[!is.na(test_data_eval$SOCIAL_Ethnicity)])
  ethnicity_results <- rbindlist(lapply(ethnicity_vals, function(eth) {
    calc_subgroup_metrics(test_data_eval, "SOCIAL_Ethnicity", eth)
  }))
  subgroup_results <- rbindlist(list(subgroup_results, ethnicity_results))
}

fwrite(subgroup_results, file.path(csv_dir, "subgroup_fairness_analysis.csv"))

subgroup_plot_data <- subgroup_results[!is.na(AUROC) & Subgroup != "All"]
subgroup_plot_data[, Subgroup_Label := paste0(Subgroup, " (N=", format(N, big.mark = ","), ")")]

if (nrow(subgroup_plot_data) > 0) {
  p_subgroup <- ggplot(subgroup_plot_data, aes(x = reorder(Subgroup_Label, AUROC), y = AUROC)) +
    geom_point(aes(size = N), color = "#E94F37") +
    geom_hline(yintercept = auc(test_roc), linetype = "dashed", color = "#2E86AB", linewidth = 1) +
    geom_text(aes(label = sprintf("%.3f", AUROC)), hjust = -0.3, size = 3) +
    coord_flip() +
    facet_wrap(~Subgroup_Type, scales = "free_y", ncol = 1) +
    scale_size_continuous(range = c(2, 8), labels = scales::comma, name = "N") +
    labs(title = "TRIPOD Fairness: BAM Performance by Subgroup",
         subtitle = sprintf("Dashed line = overall AUROC (%.3f) | Youden threshold", auc(test_roc)),
         x = "", y = "AUROC") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom") +
    ylim(0.5, 1)
  ggsave(file.path(fig_dir, "13_subgroup_fairness_auroc.png"), p_subgroup, width = 10, height = 12, dpi = 300)
}

################################################################################
# MISSINGNESS-STRATIFIED PERFORMANCE (Deployment Readiness Analysis)
# Per Tsvetanova et al. - understand performance by data availability pattern
# Note: BAM uses median imputation, but we track ORIGINAL missingness
################################################################################

# Copy original missingness indicators to test data
test_data_eval[, VS_Complete := VS_Complete_orig]
test_data_eval[, FBC_Complete := FBC_Complete_orig]
test_data_eval[, UE_Complete := UE_Complete_orig]
test_data_eval[, LFT_Complete := LFT_Complete_orig]
test_data_eval[, CRP_Complete := CRP_Complete_orig]
test_data_eval[, BG_Complete := BG_Complete_orig]
test_data_eval[, RoutineBloods_Complete := RoutineBloods_Complete_orig]

# Function to calculate metrics for a subset
calc_missingness_metrics <- function(dt, subset_expr, subset_name) {
  sub_dt <- dt[eval(parse(text = subset_expr))]
  if (nrow(sub_dt) < 50 || sum(sub_dt$Outcome) < 10) {
    return(data.table(
      Pattern = subset_name, N = nrow(sub_dt), N_Positive = sum(sub_dt$Outcome),
      Prevalence = NA_real_, AUROC = NA_real_, AUROC_CI_Lower = NA_real_, AUROC_CI_Upper = NA_real_
    ))
  }
  sub_roc <- tryCatch(roc(sub_dt$Outcome, sub_dt$pred, quiet = TRUE), error = function(e) NULL)
  if (is.null(sub_roc)) {
    return(data.table(
      Pattern = subset_name, N = nrow(sub_dt), N_Positive = sum(sub_dt$Outcome),
      Prevalence = 100 * mean(sub_dt$Outcome), AUROC = NA_real_, 
      AUROC_CI_Lower = NA_real_, AUROC_CI_Upper = NA_real_
    ))
  }
  ci <- ci.auc(sub_roc)
  data.table(
    Pattern = subset_name, N = nrow(sub_dt), N_Positive = sum(sub_dt$Outcome),
    Prevalence = round(100 * mean(sub_dt$Outcome), 1),
    AUROC = round(as.numeric(auc(sub_roc)), 3),
    AUROC_CI_Lower = round(ci[1], 3), AUROC_CI_Upper = round(ci[3], 3)
  )
}

# Individual domain analysis (complete vs missing)
individual_patterns <- rbindlist(list(
  calc_missingness_metrics(test_data_eval, "TRUE", "Overall"),
  calc_missingness_metrics(test_data_eval, "VS_Complete == TRUE", "VS Complete"),
  calc_missingness_metrics(test_data_eval, "VS_Complete == FALSE", "VS Missing"),
  calc_missingness_metrics(test_data_eval, "FBC_Complete == TRUE", "FBC Complete"),
  calc_missingness_metrics(test_data_eval, "FBC_Complete == FALSE", "FBC Missing"),
  calc_missingness_metrics(test_data_eval, "UE_Complete == TRUE", "U&E Complete"),
  calc_missingness_metrics(test_data_eval, "UE_Complete == FALSE", "U&E Missing"),
  calc_missingness_metrics(test_data_eval, "LFT_Complete == TRUE", "LFT Complete"),
  calc_missingness_metrics(test_data_eval, "LFT_Complete == FALSE", "LFT Missing"),
  calc_missingness_metrics(test_data_eval, "CRP_Complete == TRUE", "CRP Complete"),
  calc_missingness_metrics(test_data_eval, "CRP_Complete == FALSE", "CRP Missing"),
  calc_missingness_metrics(test_data_eval, "RoutineBloods_Complete == TRUE", "Routine Bloods Complete"),
  calc_missingness_metrics(test_data_eval, "RoutineBloods_Complete == FALSE", "Routine Bloods Missing"),
  calc_missingness_metrics(test_data_eval, "BG_Complete == TRUE", "Blood Gas Complete"),
  calc_missingness_metrics(test_data_eval, "BG_Complete == FALSE", "Blood Gas Missing")
))

# Trio combinations: VS, Routine Bloods, Blood Gas
trio_patterns <- rbindlist(list(
  calc_missingness_metrics(test_data_eval, 
    "VS_Complete & RoutineBloods_Complete & BG_Complete", 
    "All Complete (VS + Bloods + BG)"),
  calc_missingness_metrics(test_data_eval, 
    "!VS_Complete & !RoutineBloods_Complete & !BG_Complete", 
    "All Missing (VS + Bloods + BG)"),
  calc_missingness_metrics(test_data_eval, 
    "VS_Complete & RoutineBloods_Complete & !BG_Complete", 
    "VS + Bloods Complete, BG Missing"),
  calc_missingness_metrics(test_data_eval, 
    "VS_Complete & !RoutineBloods_Complete & BG_Complete", 
    "VS + BG Complete, Bloods Missing"),
  calc_missingness_metrics(test_data_eval, 
    "!VS_Complete & RoutineBloods_Complete & BG_Complete", 
    "Bloods + BG Complete, VS Missing"),
  calc_missingness_metrics(test_data_eval, 
    "VS_Complete & !RoutineBloods_Complete & !BG_Complete", 
    "Only VS Complete"),
  calc_missingness_metrics(test_data_eval, 
    "!VS_Complete & RoutineBloods_Complete & !BG_Complete", 
    "Only Bloods Complete"),
  calc_missingness_metrics(test_data_eval, 
    "!VS_Complete & !RoutineBloods_Complete & BG_Complete", 
    "Only BG Complete")
))

# Pair combinations
pair_patterns <- rbindlist(list(
  calc_missingness_metrics(test_data_eval, "VS_Complete & RoutineBloods_Complete", "VS + Bloods Complete"),
  calc_missingness_metrics(test_data_eval, "!VS_Complete & !RoutineBloods_Complete", "VS + Bloods Missing"),
  calc_missingness_metrics(test_data_eval, "VS_Complete & BG_Complete", "VS + BG Complete"),
  calc_missingness_metrics(test_data_eval, "!VS_Complete & !BG_Complete", "VS + BG Missing"),
  calc_missingness_metrics(test_data_eval, "RoutineBloods_Complete & BG_Complete", "Bloods + BG Complete"),
  calc_missingness_metrics(test_data_eval, "!RoutineBloods_Complete & !BG_Complete", "Bloods + BG Missing")
))

# Combine all patterns
all_missingness_patterns <- rbindlist(list(
  individual_patterns[, Category := "Individual Domains"],
  trio_patterns[, Category := "Trio Combinations"],
  pair_patterns[, Category := "Pair Combinations"]
))

fwrite(all_missingness_patterns, file.path(csv_dir, "missingness_stratified_performance.csv"))

# Summary statistics on data availability
availability_summary <- data.table(
  Domain = c("Vital Signs", "FBC", "U&E", "LFT", "CRP", "Routine Bloods (All)", "Blood Gas"),
  N_Complete = c(sum(test_data_eval$VS_Complete), sum(test_data_eval$FBC_Complete),
                 sum(test_data_eval$UE_Complete), sum(test_data_eval$LFT_Complete),
                 sum(test_data_eval$CRP_Complete), sum(test_data_eval$RoutineBloods_Complete),
                 sum(test_data_eval$BG_Complete)),
  Pct_Complete = c(100*mean(test_data_eval$VS_Complete), 100*mean(test_data_eval$FBC_Complete),
                   100*mean(test_data_eval$UE_Complete), 100*mean(test_data_eval$LFT_Complete),
                   100*mean(test_data_eval$CRP_Complete), 100*mean(test_data_eval$RoutineBloods_Complete),
                   100*mean(test_data_eval$BG_Complete))
)
availability_summary[, Pct_Complete := round(Pct_Complete, 1)]
fwrite(availability_summary, file.path(csv_dir, "data_availability_summary.csv"))

# Visualization: Individual domain impact (matching XGBoost style)
individual_plot_data <- individual_patterns[Pattern != "Overall" & !is.na(AUROC)]
individual_plot_data[, Status := fifelse(grepl("Complete", Pattern), "Complete", "Missing")]
individual_plot_data[, Status := factor(Status, levels = c("Complete", "Missing"))]
individual_plot_data[, Domain := gsub(" Complete| Missing", "", Pattern)]

p_individual <- ggplot(individual_plot_data, aes(x = reorder(Domain, -AUROC), y = AUROC, fill = Status)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = AUROC_CI_Lower, ymax = AUROC_CI_Upper), 
                position = position_dodge(width = 0.8), width = 0.2, color = "gray30") +
  geom_hline(yintercept = individual_patterns[Pattern == "Overall"]$AUROC, 
             linetype = "dashed", color = "red", linewidth = 1) +
  geom_text(aes(label = sprintf("%.3f", AUROC)), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 2.8) +
  scale_fill_manual(values = c("Complete" = "#2E86AB", "Missing" = "#E94F37"),
                    labels = c("Data Complete", "Data Missing")) +
  labs(title = "BAM Performance by Data Availability (Median Imputed)",
       subtitle = sprintf("Dashed line = overall AUROC (%.3f) | Error bars = 95%% CI", 
                         individual_patterns[Pattern == "Overall"]$AUROC),
       x = "", y = "AUROC", fill = "Data Status") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  coord_cartesian(ylim = c(0.5, 1))

ggsave(file.path(fig_dir, "14_missingness_impact_individual.png"), p_individual, width = 12, height = 8, dpi = 300)

# Visualization: Trio combinations (matching XGBoost style)
trio_plot_data <- trio_patterns[!is.na(AUROC) & N >= 50]
trio_plot_data[, Pattern_Short := gsub(" \\(VS \\+ Bloods \\+ BG\\)", "", Pattern)]

p_trio <- ggplot(trio_plot_data, aes(x = reorder(Pattern_Short, AUROC), y = AUROC)) +
  geom_bar(stat = "identity", fill = "#2E86AB", width = 0.7) +
  geom_errorbar(aes(ymin = AUROC_CI_Lower, ymax = AUROC_CI_Upper), width = 0.2) +
  geom_hline(yintercept = individual_patterns[Pattern == "Overall"]$AUROC, 
             linetype = "dashed", color = "red", linewidth = 1) +
  geom_text(aes(label = sprintf("%.3f\n(N=%s)", AUROC, format(N, big.mark = ","))), 
            hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "BAM Performance by VS + Routine Bloods + Blood Gas Combinations",
       subtitle = sprintf("Dashed line = overall AUROC (%.3f) | Median imputation used", 
                         individual_patterns[Pattern == "Overall"]$AUROC),
       x = "", y = "AUROC") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14)) +
  scale_y_continuous(limits = c(0.5, 1.1))

ggsave(file.path(fig_dir, "15_missingness_impact_trio.png"), p_trio, width = 12, height = 8, dpi = 300)

# Heatmap: Trio availability patterns
trio_heatmap_data <- test_data_eval[, .(
  N = .N, 
  Prevalence = round(100 * mean(Outcome), 1),
  AUROC = {
    if (.N >= 50 && sum(Outcome) >= 10) {
      r <- tryCatch(roc(Outcome, pred, quiet = TRUE), error = function(e) NULL)
      if (!is.null(r)) round(as.numeric(auc(r)), 3) else NA_real_
    } else NA_real_
  }
), by = .(VS_Complete, RoutineBloods_Complete, BG_Complete)]

trio_heatmap_data[, VS_Label := fifelse(VS_Complete, "VS ✓", "VS ✗")]
trio_heatmap_data[, Bloods_Label := fifelse(RoutineBloods_Complete, "Bloods ✓", "Bloods ✗")]
trio_heatmap_data[, BG_Label := fifelse(BG_Complete, "BG ✓", "BG ✗")]
trio_heatmap_data[, Cell_Label := sprintf("AUROC: %.3f\nN: %s\nPrev: %.1f%%", 
                                          AUROC, format(N, big.mark = ","), Prevalence)]
trio_heatmap_data[is.na(AUROC), Cell_Label := sprintf("N: %s\n(too few)", format(N, big.mark = ","))]

p_heatmap <- ggplot(trio_heatmap_data, aes(x = Bloods_Label, y = VS_Label, fill = AUROC)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Cell_Label), size = 3) +
  facet_wrap(~BG_Label) +
  scale_fill_gradient2(low = "#2E86AB", mid = "#F5F5F5", high = "#E94F37", 
                       midpoint = individual_patterns[Pattern == "Overall"]$AUROC,
                       na.value = "gray80", limits = c(0.5, 1)) +
  labs(title = "BAM Performance Heatmap: Data Availability Patterns",
       subtitle = "VS × Routine Bloods × Blood Gas combinations (median imputed)",
       x = "", y = "", fill = "AUROC") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, "16_missingness_heatmap.png"), p_heatmap, width = 10, height = 6, dpi = 300)

################################################################################
# TABLE 2 METRICS — saved for compilation script
################################################################################

set.seed(42)
n_boot <- 1000
boot_auprcs_test <- vapply(seq_len(n_boot), function(b) {
  idx <- sample(length(test_label), replace = TRUE)
  bp  <- test_pred[idx]; bl <- test_label[idx]
  valid <- !is.na(bp)
  if (sum(bl[valid] == 1) > 0 && sum(bl[valid] == 0) > 0)
    tryCatch(
      pr.curve(scores.class0 = bp[valid & bl == 1], scores.class1 = bp[valid & bl == 0])$auc.integral,
      error = function(e) NA_real_
    )
  else NA_real_
}, numeric(1))

auroc_ci   <- ci.auc(test_roc)
auprc_ci   <- quantile(boot_auprcs_test, c(0.025, 0.975), na.rm = TRUE)

results_dir <- "Output Files/Publication_Tables_Figures/Metrics"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

table2_metrics <- data.frame(
  Model                = "BAM",
  CV_AUROC_Median      = round(median(cv_aucs), 3),
  CV_AUROC_Q1          = round(quantile(cv_aucs, 0.25), 3),
  CV_AUROC_Q3          = round(quantile(cv_aucs, 0.75), 3),
  CV_AUPRC_Median      = round(median(cv_auprcs), 3),
  CV_AUPRC_Q1          = round(quantile(cv_auprcs, 0.25), 3),
  CV_AUPRC_Q3          = round(quantile(cv_auprcs, 0.75), 3),
  Test_AUROC           = round(as.numeric(auc(test_roc)), 3),
  Test_AUROC_CI_Lower  = round(auroc_ci[1], 3),
  Test_AUROC_CI_Upper  = round(auroc_ci[3], 3),
  Test_AUPRC           = round(pr_obj$auc.integral, 3),
  Test_AUPRC_CI_Lower  = round(auprc_ci[1], 3),
  Test_AUPRC_CI_Upper  = round(auprc_ci[2], 3),
  Test_Sensitivity     = round(sensitivity, 3),
  Test_Specificity     = round(specificity, 3),
  Test_PPV             = round(ppv, 3),
  Test_NPV             = round(npv, 3),
  Test_Brier           = round(brier_calibrated, 3)
)

fwrite(table2_metrics, file.path(results_dir, "table2_metrics_bam.csv"))