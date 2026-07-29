################################################################################
# XGBoost Gram-Stain Prediction — FULL COHORT
# Clinically actionable at point of culture collection (before knowing result)
# 4-class: 0=Negative/Contaminant  1=Gram-positive  2=Gram-negative  3=Other pathogen
# Same 31 features as binary FullCohort model
################################################################################

library(data.table)
library(xgboost)
library(pROC)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(PRROC)
library(mgcv)
library(cowplot)
library(writexl)

set.seed(42)

fig_dir   <- "Output Files/Figures/Models/XGBoost/GramPredictFullCohort"
csv_dir   <- "Output Files/CSVs/XGBoost/GramPredictFullCohort"
model_dir <- "Models/XGBoost/GramPredictFullCohort"
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(csv_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

bc_data <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))

bc_data <- bc_data[CultivabilityCategory %in% c(
  "Gram-positive", "Gram-negative", "Polymicrobial",
  "Fungi", "Anaerobes", "Other",
  "Culture-negative", "Contaminant - single", "Contaminant - recurrent",
  "Polymicrobial - contaminant"
)]

# 4-class label
# 0 = Negative (no significant growth): culture-neg + all contaminants
# 1 = Gram-positive
# 2 = Gram-negative
# 3 = Other pathogen: Polymicrobial, Fungi, Anaerobes, Other
bc_data[, GramLabel := fcase(
  CultivabilityCategory %in% c("Culture-negative",
    "Contaminant - single", "Contaminant - recurrent",
    "Polymicrobial - contaminant"),         0L,
  CultivabilityCategory == "Gram-positive", 1L,
  CultivabilityCategory == "Gram-negative", 2L,
  default = 3L   # Polymicrobial, Fungi, Anaerobes, Other
)]

class_labels <- c("Negative/Contaminant", "Gram-positive", "Gram-negative", "Other pathogen")
num_class    <- 4L

################################################################################
# FEATURE ENGINEERING (identical to binary FullCohort model)
################################################################################

bc_data[, Sex_num := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]

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

bc_data[, VS_AnyFever       := fifelse(VS_AnyFever_38    == 1 & !is.na(VS_AnyFever_38),    1, 0)]
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
  "VS_Baseline_NEWS", "VS_AnyFever",
  "BT_Closest_CRP", "BT_Closest_WCC", "BT_Closest_Neut", "BT_Closest_Lymph",
  "BT_Closest_Alb", "BT_Closest_Plt", "BT_Closest_Creat", "BT_Closest_Bili",
  "BT_Closest_Na",
  "BT_Closest_Eosin", "BT_Closest_Mono", "BT_Closest_ALP",
  "BG_Closest_Lactate", "BG_Closest_pH"
)

################################################################################
# TEMPORAL SPLIT
################################################################################

bc_data[, CollectionDate := as.Date(CollectionDateTime)]
split_date <- as.Date("2024-01-01")
train_data <- bc_data[CollectionDate <  split_date]
test_data  <- bc_data[CollectionDate >= split_date]

################################################################################
# 5-FOLD CV (patient-level, stratified by class)
################################################################################

k_folds <- 5L
fold_assignments <- data.table(ClusterID = character(), Fold = integer())

for (lbl in 0:(num_class - 1L)) {
  pats <- unique(train_data[GramLabel == lbl]$ClusterID)
  pats <- sample(pats)
  for (i in seq_len(k_folds)) {
    fold_pats <- pats[seq(i, length(pats), by = k_folds)]
    fold_assignments <- rbind(fold_assignments,
                              data.table(ClusterID = fold_pats, Fold = i))
  }
}
fold_assignments <- fold_assignments[!duplicated(ClusterID)]
train_data <- merge(train_data, fold_assignments, by = "ClusterID", all.x = TRUE)

################################################################################
# MODEL PARAMETERS
################################################################################

params <- list(
  objective        = "multi:softprob",
  eval_metric      = "mlogloss",
  num_class        = num_class,
  max_depth        = 6,
  eta              = 0.05,
  min_child_weight = 10,
  subsample        = 0.8,
  colsample_bytree = 0.5,
  gamma            = 2,
  lambda            = 30,
  alpha            = 0,
  nthread          = 8
)

################################################################################
# 5-FOLD CV FOR OOF PREDICTIONS
################################################################################

oof_cols <- c("oof_Neg", "oof_GP", "oof_GN", "oof_Other")
train_data[, (oof_cols) := .(NA_real_, NA_real_, NA_real_, NA_real_)]
cv_mlogloss <- numeric(k_folds)

for (fold in seq_len(k_folds)) {
  tr  <- train_data[Fold != fold]
  val <- train_data[Fold == fold]

  dtr  <- xgb.DMatrix(data = as.matrix(tr[,  ..feature_cols]), label = tr$GramLabel)
  dval <- xgb.DMatrix(data = as.matrix(val[, ..feature_cols]), label = val$GramLabel)

  m_fold <- xgb.train(
    params = params, data = dtr, nrounds = 600,
    watchlist = list(train = dtr, val = dval),
    early_stopping_rounds = 50, verbose = 0
  )

  pmat <- matrix(predict(m_fold, dval), ncol = num_class, byrow = TRUE)
  train_data[Fold == fold, (oof_cols) :=
               .(pmat[, 1], pmat[, 2], pmat[, 3], pmat[, 4])]

  cv_mlogloss[fold] <- m_fold$best_score
}

################################################################################
# FINAL MODEL
################################################################################

dtrain <- xgb.DMatrix(data = as.matrix(train_data[, ..feature_cols]),
                      label = train_data$GramLabel)
dtest  <- xgb.DMatrix(data = as.matrix(test_data[, ..feature_cols]),
                      label = test_data$GramLabel)

model <- xgb.train(
  params = params, data = dtrain, nrounds = 600,
  watchlist = list(train = dtrain, test = dtest),
  early_stopping_rounds = 50, verbose = 0
)

xgb.save(model, file.path(model_dir, "xgb_grampredict_fullcohort.json"))

# Test predictions
test_pred_mat <- matrix(predict(model, dtest), ncol = num_class, byrow = TRUE)
colnames(test_pred_mat) <- c("P_Neg", "P_GP", "P_GN", "P_Other")

test_pred_class <- apply(test_pred_mat, 1, which.max) - 1L
test_true       <- test_data$GramLabel

accuracy <- mean(test_pred_class == test_true)

################################################################################
# PER-CLASS METRICS (one-vs-rest)
################################################################################

gram_colors <- c(
  "Negative/Contaminant" = "#6A9E5B",
  "Gram-positive"        = "#E94F37",
  "Gram-negative"        = "#2E86AB",
  "Other pathogen"       = "#F4A261"
)

calc_ovr <- function(true_labels, pred_mat, class_idx, class_name) {
  binary_true <- as.integer(true_labels == class_idx)
  binary_prob <- pred_mat[, class_idx + 1L]
  if (sum(binary_true) < 5L) return(NULL)
  roc_obj <- tryCatch(roc(binary_true, binary_prob, quiet = TRUE), error = function(e) NULL)
  if (is.null(roc_obj)) return(NULL)
  ci_val  <- ci.auc(roc_obj)
  thresh  <- coords(roc_obj, "best", ret = "threshold")$threshold
  pred_b  <- as.integer(binary_prob >= thresh)
  TP <- sum(pred_b == 1L & binary_true == 1L)
  TN <- sum(pred_b == 0L & binary_true == 0L)
  FP <- sum(pred_b == 1L & binary_true == 0L)
  FN <- sum(pred_b == 0L & binary_true == 1L)
  data.table(
    Class          = class_name,
    N              = sum(binary_true),
    Prevalence_Pct = round(100 * mean(binary_true), 1),
    AUROC          = round(as.numeric(auc(roc_obj)), 3),
    AUROC_CI_Lower = round(ci_val[1], 3),
    AUROC_CI_Upper = round(ci_val[3], 3),
    Sensitivity    = round(TP / (TP + FN), 3),
    Specificity    = round(TN / (TN + FP), 3),
    PPV            = round(TP / (TP + FP), 3),
    NPV            = round(TN / (TN + FN), 3)
  )
}

ovr_metrics <- rbindlist(list(
  calc_ovr(test_true, test_pred_mat, 0L, "Negative/Contaminant"),
  calc_ovr(test_true, test_pred_mat, 1L, "Gram-positive"),
  calc_ovr(test_true, test_pred_mat, 2L, "Gram-negative"),
  calc_ovr(test_true, test_pred_mat, 3L, "Other pathogen")
))

macro_auroc <- round(mean(ovr_metrics$AUROC), 3)
fwrite(ovr_metrics, file.path(csv_dir, "per_class_metrics.csv"))

################################################################################
# FIGURE 1: CONFUSION MATRIX
################################################################################

conf_mat  <- table(
  Predicted = factor(class_labels[test_pred_class + 1L], levels = class_labels),
  Actual    = factor(class_labels[test_true      + 1L], levels = class_labels)
)
conf_pct  <- round(100 * prop.table(conf_mat, margin = 2), 1)

conf_df   <- merge(as.data.frame(conf_pct), as.data.frame(conf_mat),
                   by = c("Predicted", "Actual"))
names(conf_df)[3:4] <- c("Pct", "N")
conf_df$Label <- paste0(conf_df$N, "\n(", conf_df$Pct, "%)")

p_conf <- ggplot(conf_df, aes(x = Actual, y = Predicted, fill = Pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Label), size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#2E86AB", limits = c(0, 100),
                      name = "% of Actual") +
  labs(title    = "Confusion Matrix — Gram Prediction (Full Cohort)",
       subtitle = sprintf("Accuracy = %.1f%% | Macro AUROC = %.3f | N test = %s",
                          100 * accuracy, macro_auroc,
                          format(nrow(test_data), big.mark = ",")),
       x = "Actual Category", y = "Predicted Category") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text  = element_text(size = 10),
        axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(fig_dir, "01_confusion_matrix.png"), p_conf, width = 9, height = 7, dpi = 300)

################################################################################
# FIGURE 2: PER-CLASS ROC CURVES
################################################################################

roc_df <- rbindlist(lapply(0:3, function(ci) {
  bt  <- as.integer(test_true == ci)
  bp  <- test_pred_mat[, ci + 1L]
  ro  <- tryCatch(roc(bt, bp, quiet = TRUE), error = function(e) NULL)
  if (is.null(ro)) return(NULL)
  auroc_val <- round(as.numeric(auc(ro)), 3)
  data.table(
    FPR   = 1 - ro$specificities,
    TPR   = ro$sensitivities,
    label = sprintf("%s (AUROC=%.3f)", class_labels[ci + 1L], auroc_val),
    Class = class_labels[ci + 1L]
  )
}))

color_map <- setNames(gram_colors, sprintf("%s (AUROC=%.3f)", class_labels,
                                           ovr_metrics$AUROC))

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = label)) +
  geom_line(linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = color_map) +
  labs(title    = "One-vs-Rest ROC Curves — Gram Prediction (Full Cohort)",
       subtitle = sprintf("Macro AUROC = %.3f", macro_auroc),
       x = "False Positive Rate", y = "True Positive Rate", color = "") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom") +
  coord_equal() +
  guides(color = guide_legend(nrow = 2))

ggsave(file.path(fig_dir, "02_roc_curves_per_class.png"), p_roc, width = 8, height = 9, dpi = 300)

################################################################################
# FIGURE 3: PREDICTED PROBABILITY DISTRIBUTIONS BY TRUE CLASS
################################################################################

prob_long <- as.data.table(test_pred_mat)
prob_long[, True := factor(class_labels[test_true + 1L], levels = class_labels)]

prob_melt <- melt(prob_long, id.vars = "True",
                  variable.name = "Probability_For", value.name = "Probability")
prob_melt[, Probability_For := factor(
  fcase(
    Probability_For == "P_Neg",   "P(Negative/Contaminant)",
    Probability_For == "P_GP",    "P(Gram-positive)",
    Probability_For == "P_GN",    "P(Gram-negative)",
    Probability_For == "P_Other", "P(Other pathogen)"
  ),
  levels = c("P(Negative/Contaminant)", "P(Gram-positive)",
             "P(Gram-negative)", "P(Other pathogen)")
)]

p_prob <- ggplot(prob_melt, aes(x = Probability, fill = True, color = True)) +
  geom_density(alpha = 0.3, linewidth = 0.7) +
  facet_wrap(~Probability_For, ncol = 2, scales = "free_y") +
  scale_fill_manual(values  = gram_colors) +
  scale_color_manual(values = gram_colors) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Predicted Probability Distributions by True Category",
       x = "Predicted Probability", y = "Density",
       fill = "True Category", color = "True Category") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "03_probability_distributions.png"),
       p_prob, width = 12, height = 9, dpi = 300)

################################################################################
# FIGURE 4: OVERALL FEATURE IMPORTANCE
################################################################################

importance <- xgb.importance(model = model)
importance[, Pct_Gain := round(100 * Gain, 2)]

p_imp <- ggplot(importance[seq_len(min(25L, nrow(importance)))],
                aes(x = reorder(Feature, Gain), y = Pct_Gain)) +
  geom_bar(stat = "identity", fill = "#2E86AB", width = 0.7) +
  geom_text(aes(label = paste0(Pct_Gain, "%")), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title    = "XGBoost Feature Importance (Top 25) — Gram Prediction (Full Cohort)",
       subtitle = sprintf("Macro AUROC = %.3f | %d features", macro_auroc, length(feature_cols)),
       x = "", y = "% Gain") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 9)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave(file.path(fig_dir, "04_feature_importance_gain.png"),
       p_imp, width = 10, height = 10, dpi = 300)
fwrite(importance, file.path(csv_dir, "feature_importance.csv"))

################################################################################
# FIGURE 5: PER-CLASS FEATURE IMPORTANCE (top 15 per class via tree filtering)
################################################################################

n_rounds_used <- model$best_iteration

class_imp_list <- rbindlist(lapply(0:(num_class - 1L), function(ci) {
  tree_idx <- seq(ci, (n_rounds_used - 1L) * num_class + ci, by = num_class)
  imp <- tryCatch(xgb.importance(model = model, trees = tree_idx), error = function(e) NULL)
  if (is.null(imp) || nrow(imp) == 0L) return(NULL)
  imp[, Class    := class_labels[ci + 1L]]
  imp[, Pct_Gain := round(100 * Gain, 2)]
  head(imp[order(-Gain)], 15L)
}))

p_class_imp <- ggplot(class_imp_list,
                      aes(x = reorder(Feature, Gain), y = Pct_Gain, fill = Class)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~Class, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = gram_colors) +
  labs(title    = "Per-Class Feature Importance (Top 15) — Gram Prediction (Full Cohort)",
       subtitle = "Gain from trees belonging to each class only",
       x = "", y = "% Gain") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 10))

ggsave(file.path(fig_dir, "05_per_class_importance.png"),
       p_class_imp, width = 20, height = 8, dpi = 300)
fwrite(class_imp_list, file.path(csv_dir, "per_class_importance.csv"))

################################################################################
# FIGURE 6: CALIBRATION PER CLASS (one-vs-rest)
################################################################################

cal_plots <- lapply(0:3, function(ci) {
  bt  <- as.integer(test_true == ci)
  bp  <- test_pred_mat[, ci + 1L]
  cal_df <- data.frame(pred = bp, actual = bt)
  bam_fit <- tryCatch(gam(actual ~ s(pred, k = 8), data = cal_df, family = binomial),
                      error = function(e) NULL)
  if (is.null(bam_fit)) return(NULL)
  pg       <- data.frame(pred = seq(0.01, 0.99, length.out = 200))
  bam_pred <- predict(bam_fit, newdata = pg, type = "response", se.fit = TRUE)
  pg$fit   <- bam_pred$fit
  pg$lower <- pmax(0, bam_pred$fit - 1.96 * bam_pred$se.fit)
  pg$upper <- pmin(1, bam_pred$fit + 1.96 * bam_pred$se.fit)
  bin_df   <- cal_df |>
    mutate(bin = cut(pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) |>
    group_by(bin) |>
    summarise(predicted = mean(pred), observed = mean(actual), n = n(), .groups = "drop") |>
    filter(!is.na(bin))
  cls_auroc <- ovr_metrics[Class == class_labels[ci + 1L]]$AUROC
  cls_col   <- gram_colors[class_labels[ci + 1L]]
  ggplot() +
    geom_ribbon(data = pg, aes(x = pred, ymin = lower, ymax = upper),
                fill = cls_col, alpha = 0.2) +
    geom_line(data = pg, aes(x = pred, y = fit), color = cls_col, linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_point(data = bin_df, aes(x = predicted, y = observed, size = n),
               color = cls_col, alpha = 0.7) +
    scale_size_continuous(range = c(2, 8), labels = scales::comma) +
    labs(title    = sprintf("Calibration: %s", class_labels[ci + 1L]),
         subtitle = sprintf("AUROC = %.3f", cls_auroc),
         x = "Predicted Probability", y = "Observed Proportion", size = "N") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 11)) +
    coord_equal() + xlim(0, 1) + ylim(0, 1)
})

p_cal <- plot_grid(plotlist = Filter(Negate(is.null), cal_plots), nrow = 2, labels = "AUTO")
ggsave(file.path(fig_dir, "06_calibration_per_class.png"),
       p_cal, width = 16, height = 12, dpi = 300)

################################################################################
# FIGURE 7: GRAM-SPECIFIC PERFORMANCE vs BINARY POSITIVITY MODEL
# For GP and GN specifically: what is our positive predictive performance
# from whole cohort (i.e. sensitivity and PPV for each gram class)?
################################################################################

# Among true positives (GP or GN only), what fraction does the model correctly
# assign to the right gram class?
gp_gn_only <- test_data[test_true %in% c(1L, 2L)]
gp_gn_pred <- test_pred_class[test_true %in% c(1L, 2L)]
gp_gn_true <- test_true[test_true %in% c(1L, 2L)]

gp_gn_among_pos <- data.table(
  Class = c("Gram-positive (among GP+GN only)", "Gram-negative (among GP+GN only)"),
  N = c(sum(gp_gn_true == 1L), sum(gp_gn_true == 2L)),
  Correctly_Classified = c(sum(gp_gn_pred == 1L & gp_gn_true == 1L),
                           sum(gp_gn_pred == 2L & gp_gn_true == 2L)),
  Accuracy_Among_GpGn = c(
    round(mean(gp_gn_pred[gp_gn_true == 1L] == 1L), 3),
    round(mean(gp_gn_pred[gp_gn_true == 2L] == 2L), 3)
  )
)

fwrite(gp_gn_among_pos, file.path(csv_dir, "gp_gn_classification_among_positives.csv"))

################################################################################
# SAVE ALL OUTPUTS
################################################################################

conf_save <- as.data.table(conf_mat)
setnames(conf_save, c("Predicted", "Actual", "N"))
fwrite(conf_save, file.path(csv_dir, "confusion_matrix.csv"))

fwrite(data.table(
  P_Neg       = test_pred_mat[, 1],
  P_GP        = test_pred_mat[, 2],
  P_GN        = test_pred_mat[, 3],
  P_Other     = test_pred_mat[, 4],
  PredClass   = class_labels[test_pred_class + 1L],
  TrueClass   = class_labels[test_true + 1L],
  Correct     = test_pred_class == test_true
), file.path(csv_dir, "test_predictions.csv"))

model_summary <- data.table(
  Metric = c("Macro_AUROC",
             "AUROC_Negative", "AUROC_GP", "AUROC_GN", "AUROC_Other",
             "Accuracy",
             "CV_MLogLoss_Median",
             "N_Train", "N_Test",
             "N_Neg_Test", "N_GP_Test", "N_GN_Test", "N_Other_Test"),
  Value  = c(macro_auroc,
             ovr_metrics[Class == "Negative/Contaminant"]$AUROC,
             ovr_metrics[Class == "Gram-positive"]$AUROC,
             ovr_metrics[Class == "Gram-negative"]$AUROC,
             ovr_metrics[Class == "Other pathogen"]$AUROC,
             round(accuracy, 3),
             round(median(cv_mlogloss), 4),
             nrow(train_data), nrow(test_data),
             sum(test_true == 0L), sum(test_true == 1L),
             sum(test_true == 2L), sum(test_true == 3L))
)
fwrite(model_summary, file.path(csv_dir, "model_summary.csv"))

excel_dir <- "Output Files/Publication_Tables_Figures/Per_Model_Summaries"
dir.create(excel_dir, showWarnings = FALSE, recursive = TRUE)
write_xlsx(list(
  Model_Summary     = as.data.frame(model_summary),
  Per_Class_Metrics = as.data.frame(ovr_metrics),
  Confusion_Matrix  = as.data.frame(conf_save),
  GP_GN_Among_Pos   = as.data.frame(gp_gn_among_pos)
), file.path(excel_dir, "XGBoost_GramPredictFullCohort_Summary.xlsx"))
