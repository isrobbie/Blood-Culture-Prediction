################################################################################
# XGBoost Bootstrap Uncertainty Estimation - Full Cohort
# Purpose: Generate prediction intervals via bootstrap aggregation
# Parallelized with low CPU priority for shared server use
################################################################################

library(data.table)
library(xgboost)
library(pROC)
library(parallel)
library(ggplot2)
library(tidyr)
library(scales)

set.seed(42)

# Configuration
N_BOOTSTRAPS <- 1000
THREADS_PER_MODEL <- 2  # Lower threads per model to allow parallelization
N_PARALLEL <- 16        # Number of models to train in parallel (16 * 2 = 32 threads, kind on 36-core server)

# Output directories
fig_dir <- "Output Files/Figures/Models/XGBoost/FullCohort"
csv_dir <- "Output Files/CSVs/XGBoost/FullCohort"
model_dir <- "Models/XGBoost/FullCohort"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)

# Load and prepare data (same as main script)
bc_data <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))

bc_data <- bc_data[CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial",
                                                  "Culture-negative", "Contaminant - single", "Contaminant - recurrent",
                                                  "Polymicrobial - contaminant")]

bc_data[, Outcome := fifelse(CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial"), 1L, 0L)]

# Feature engineering (same as main script)
bc_data[, Sex_num := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]
median_imd <- median(bc_data$SOCIAL_IMDPercentile, na.rm = TRUE)
bc_data[is.na(SOCIAL_IMDPercentile), SOCIAL_IMDPercentile := median_imd]

bc_data[, AdmissionType_num := fcase(
  AdmissionMethod == "Emergency", 1,
  AdmissionMethod == "Elective", 2,
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

bc_data[, VS_AnyFever := fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), 1, 0)]
bc_data[, VS_AnyHypotension := fifelse(VS_AnyHypotension == 1 & !is.na(VS_AnyHypotension), 1, 0)]
bc_data[is.na(LOS_BeforeCulture_days), LOS_BeforeCulture_days := 0]

bc_data[, CollectionDate := as.Date(CollectionDateTime)]
split_date <- as.Date("2024-01-01")
train_data <- bc_data[CollectionDate <  split_date]
test_data  <- bc_data[CollectionDate >= split_date]

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

train_matrix <- as.matrix(train_data[, ..feature_cols])
test_matrix <- as.matrix(test_data[, ..feature_cols])
train_label <- train_data$Outcome
test_label <- test_data$Outcome

scale_pos_weight <- sum(train_label == 0) / sum(train_label == 1)

# Load Platt calibration model from main script
platt_model <- readRDS(file.path(model_dir, "xgb_fullcohort_platt_calibration.rds"))

# Load main model to compute threshold (Youden index)
library(pROC)
main_model <- xgb.load(file.path(model_dir, "xgb_fullcohort.json"))
main_pred_raw <- predict(main_model, xgb.DMatrix(data = test_matrix))
main_pred_logit <- qlogis(pmax(pmin(main_pred_raw, 0.999), 0.001))
main_pred <- predict(platt_model, newdata = data.frame(train_pred_logit = main_pred_logit), type = "response")
main_roc <- roc(test_label, main_pred, quiet = TRUE)
threshold <- coords(main_roc, "best", ret = "threshold")$threshold

# Bootstrap function for a single iteration
train_bootstrap_model <- function(boot_idx, train_matrix, train_label, test_matrix, 
                                   feature_cols, scale_pos_weight, threads) {
  set.seed(boot_idx)
  
  # Bootstrap sample (with replacement) at patient level would be ideal,

  # but for speed we sample rows directly
  boot_idx_rows <- sample(1:nrow(train_matrix), replace = TRUE)
  boot_train_matrix <- train_matrix[boot_idx_rows, ]
  boot_train_label <- train_label[boot_idx_rows]
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 8,
    eta = 0.0553,
    min_child_weight = 10,
    subsample = 0.9,
    colsample_bytree = 0.4,
    gamma = 3.971,
    lambda = 50,
    alpha = 0,
    scale_pos_weight = scale_pos_weight,
    nthread = threads
  )
  
  dtrain <- xgb.DMatrix(data = boot_train_matrix, label = boot_train_label)
  dtest <- xgb.DMatrix(data = test_matrix)
  
  model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 500,
    watchlist = list(train = dtrain),
    early_stopping_rounds = 50,
    verbose = 0
  )
  
  # Return raw predictions (uncalibrated)
  predict(model, dtest)
}

# Run bootstraps in parallel
message(sprintf("Starting %d bootstrap iterations with %d parallel workers...", N_BOOTSTRAPS, N_PARALLEL))
message(sprintf("Each model uses %d threads. Total max threads: %d", THREADS_PER_MODEL, N_PARALLEL * THREADS_PER_MODEL))
start_time <- Sys.time()

cl <- makeCluster(N_PARALLEL)
clusterExport(cl, c("train_matrix", "train_label", "test_matrix", "feature_cols", 
                    "scale_pos_weight", "THREADS_PER_MODEL", "train_bootstrap_model"))
clusterEvalQ(cl, {
  library(xgboost)
})

bootstrap_predictions <- parLapply(cl, 1:N_BOOTSTRAPS, function(i) {
  train_bootstrap_model(i, train_matrix, train_label, test_matrix, 
                        feature_cols, scale_pos_weight, THREADS_PER_MODEL)
})

stopCluster(cl)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")
message(sprintf("Bootstrap complete in %.1f minutes", as.numeric(elapsed)))

# Convert to matrix: rows = test samples, cols = bootstrap iterations
pred_matrix_raw <- do.call(cbind, bootstrap_predictions)

# Apply Platt calibration to each bootstrap prediction
pred_matrix_calibrated <- apply(pred_matrix_raw, 2, function(preds) {
  pred_logit <- qlogis(pmax(pmin(preds, 0.999), 0.001))
  predict(platt_model, newdata = data.frame(train_pred_logit = pred_logit), type = "response")
})

# Save raw bootstrap matrix for later figure regeneration (without re-running bootstrap)
saveRDS(pred_matrix_calibrated, file.path(model_dir, "bootstrap_pred_matrix.rds"))

# Compute prediction intervals
test_results <- data.table(
  idx = 1:nrow(test_matrix),
  Outcome = test_label,
  pred_mean = rowMeans(pred_matrix_calibrated),
  pred_median = apply(pred_matrix_calibrated, 1, median),
  pred_lower_95 = apply(pred_matrix_calibrated, 1, quantile, probs = 0.025),
  pred_upper_95 = apply(pred_matrix_calibrated, 1, quantile, probs = 0.975),
  pred_lower_80 = apply(pred_matrix_calibrated, 1, quantile, probs = 0.10),
  pred_upper_80 = apply(pred_matrix_calibrated, 1, quantile, probs = 0.90),
  pred_sd = apply(pred_matrix_calibrated, 1, sd)
)

test_results[, interval_width_95 := pred_upper_95 - pred_lower_95]
test_results[, interval_width_80 := pred_upper_80 - pred_lower_80]

# Save results
fwrite(test_results, file.path(csv_dir, "bootstrap_prediction_intervals.csv"))

# Summary statistics
interval_summary <- data.table(
  Metric = c("Mean 95% CI Width", "Median 95% CI Width", "Mean 80% CI Width", 
             "Median 80% CI Width", "Mean Prediction SD", "Bootstrap Iterations",
             "Runtime (minutes)"),
  Value = c(round(mean(test_results$interval_width_95), 4),
            round(median(test_results$interval_width_95), 4),
            round(mean(test_results$interval_width_80), 4),
            round(median(test_results$interval_width_80), 4),
            round(mean(test_results$pred_sd), 4),
            N_BOOTSTRAPS,
            round(as.numeric(elapsed), 1))
)
fwrite(interval_summary, file.path(csv_dir, "bootstrap_uncertainty_summary.csv"))

# Visualizations
# 1. Uncertainty by predicted probability
test_results[, pred_bin := cut(pred_mean, breaks = seq(0, 1, 0.1), include.lowest = TRUE)]
uncertainty_by_pred <- test_results[, .(
  mean_width = mean(interval_width_95),
  mean_sd = mean(pred_sd),
  n = .N
), by = pred_bin]

# 2. Calibration of interval coverage
coverage_95 <- mean(test_results$Outcome >= test_results$pred_lower_95 & 
                    test_results$Outcome <= test_results$pred_upper_95)
coverage_80 <- mean(test_results$Outcome >= test_results$pred_lower_80 & 
                    test_results$Outcome <= test_results$pred_upper_80)

# For binary outcomes, coverage is tricky - check if observed rate falls in interval by bins
coverage_by_bin <- test_results[, .(
  observed_rate = mean(Outcome),
  mean_pred = mean(pred_mean),
  mean_lower = mean(pred_lower_95),
  mean_upper = mean(pred_upper_95),
  n = .N
), by = pred_bin][!is.na(pred_bin)]

coverage_by_bin[, covered := observed_rate >= mean_lower & observed_rate <= mean_upper]

################################################################################
# INDIVIDUAL UNCERTAINTY DISTRIBUTIONS (violin plots)
################################################################################

# Select example patients across risk spectrum
set.seed(42)
example_idx <- c(
  sample(which(test_results$pred_mean < 0.05), min(2, sum(test_results$pred_mean < 0.05))),
  sample(which(test_results$pred_mean >= 0.05 & test_results$pred_mean < 0.15), 
         min(2, sum(test_results$pred_mean >= 0.05 & test_results$pred_mean < 0.15))),
  sample(which(test_results$pred_mean >= 0.15 & test_results$pred_mean < 0.30), 
         min(2, sum(test_results$pred_mean >= 0.15 & test_results$pred_mean < 0.30))),
  sample(which(test_results$pred_mean >= 0.30), min(2, sum(test_results$pred_mean >= 0.30)))
)

if (length(example_idx) >= 4) {
  example_preds <- as.data.frame(t(pred_matrix_calibrated[example_idx, ]))
  names(example_preds) <- paste0("Patient ", seq_along(example_idx))
  example_long <- pivot_longer(example_preds, everything(), 
                               names_to = "Patient", values_to = "Predicted_Risk")
  
  example_summary <- data.frame(
    Patient = paste0("Patient ", seq_along(example_idx)),
    Point_Estimate = test_results$pred_mean[example_idx],
    Actual = factor(ifelse(test_results$Outcome[example_idx] == 1, "Positive", "Negative"))
  )
  example_long <- merge(example_long, example_summary, by = "Patient")
  
  p_individual_uncertainty <- ggplot(example_long, aes(x = reorder(Patient, Point_Estimate), 
                                                        y = Predicted_Risk)) +
    geom_violin(fill = "#2E86AB", alpha = 0.4, color = NA) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    geom_point(aes(y = Point_Estimate, color = Actual), size = 3) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text", x = 0.6, y = threshold + 0.03, label = "Decision threshold", 
             hjust = 0, color = "red", size = 3) +
    scale_color_manual(values = c("Positive" = "#E94F37", "Negative" = "#2E86AB")) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    coord_flip() +
    labs(title = "Individual Prediction Uncertainty Distributions",
         subtitle = sprintf("Bootstrap distributions (N=%d) | Violin=uncertainty, dot=point estimate", N_BOOTSTRAPS),
         x = "", y = "Predicted Probability of Positive Culture",
         color = "Actual Outcome") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")
  
  ggsave(file.path(fig_dir, "13_individual_uncertainty_violin.png"), 
         p_individual_uncertainty, width = 10, height = 8, dpi = 300)
}

################################################################################
# CLASSIFICATION INSTABILITY
################################################################################

# Calculate instability: how often does bootstrap classification differ from point estimate?
point_class <- ifelse(test_results$pred_mean >= threshold, 1, 0)
instability <- rowMeans(ifelse(pred_matrix_calibrated >= threshold, 1, 0) != point_class)
test_results[, classification_instability := instability]
test_results[, spans_threshold := (pred_lower_95 < threshold) & (pred_upper_95 > threshold)]

# Update saved results with instability metrics
fwrite(test_results, file.path(csv_dir, "bootstrap_prediction_intervals.csv"))

################################################################################
# DCA x BOOTSTRAP UNCERTAINTY: Threshold Crossing Analysis
################################################################################

# For each threshold, count predictions where CI crosses the threshold
thresholds <- seq(0.01, 0.50, by = 0.01)

crossing_analysis <- lapply(thresholds, function(pt) {
  # 95% CI crosses threshold
  crosses_95 <- test_results$pred_lower_95 < pt & test_results$pred_upper_95 > pt
  
  # 80% CI crosses threshold
  crosses_80 <- test_results$pred_lower_80 < pt & test_results$pred_upper_80 > pt
  
  # Point estimate classification
  above_threshold <- test_results$pred_mean >= pt
  
  # Breakdown by point estimate position
  # "Uncertain positive": point estimate >= threshold, but lower CI < threshold
  uncertain_pos <- above_threshold & test_results$pred_lower_95 < pt
  # "Uncertain negative": point estimate < threshold, but upper CI > threshold  
  uncertain_neg <- !above_threshold & test_results$pred_upper_95 > pt
  
  data.table(
    Threshold = pt,
    N_Cross_95CI = sum(crosses_95),
    Pct_Cross_95CI = 100 * mean(crosses_95),
    N_Cross_80CI = sum(crosses_80),
    Pct_Cross_80CI = 100 * mean(crosses_80),
    N_Uncertain_Pos = sum(uncertain_pos),
    N_Uncertain_Neg = sum(uncertain_neg),
    N_Confident_Pos = sum(above_threshold & !uncertain_pos),
    N_Confident_Neg = sum(!above_threshold & !uncertain_neg)
  )
})

crossing_df <- rbindlist(crossing_analysis)

# Add outcome-stratified analysis
crossing_by_outcome <- lapply(thresholds, function(pt) {
  crosses_95_pos <- test_results[Outcome == 1, pred_lower_95 < pt & pred_upper_95 > pt]
  crosses_95_neg <- test_results[Outcome == 0, pred_lower_95 < pt & pred_upper_95 > pt]
  
  data.table(
    Threshold = pt,
    Pct_Cross_95_TruePos = 100 * mean(crosses_95_pos),
    Pct_Cross_95_TrueNeg = 100 * mean(crosses_95_neg)
  )
})

crossing_outcome_df <- rbindlist(crossing_by_outcome)
crossing_df <- merge(crossing_df, crossing_outcome_df, by = "Threshold")

fwrite(crossing_df, file.path(csv_dir, "dca_uncertainty_crossing.csv"))

# Key thresholds summary
key_thresholds <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
key_summary <- crossing_df[Threshold %in% key_thresholds, .(
  Threshold = paste0(100 * Threshold, "%"),
  `N Uncertain (95% CI)` = N_Cross_95CI,
  `% Uncertain` = round(Pct_Cross_95CI, 1),
  `Confident Pos` = N_Confident_Pos,
  `Uncertain Pos` = N_Uncertain_Pos,
  `Uncertain Neg` = N_Uncertain_Neg,
  `Confident Neg` = N_Confident_Neg
)]

fwrite(key_summary, file.path(csv_dir, "dca_uncertainty_key_thresholds.csv"))

message("Bootstrap uncertainty estimation complete!")
message(sprintf("Results saved to %s", csv_dir))
message(sprintf("At 10%% threshold: %.1f%% of predictions have uncertain classification", 
                crossing_df[Threshold == 0.10]$Pct_Cross_95CI))
message(sprintf("At 20%% threshold: %.1f%% of predictions have uncertain classification", 
                crossing_df[Threshold == 0.20]$Pct_Cross_95CI))

