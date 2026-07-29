library(data.table)
library(xgboost)
library(lubridate)
library(pROC)
library(ggplot2)
library(patchwork)

# ============================================================
# FIGURES 1-3 - reconstructed from the saved model + pipeline
# outputs (the original plotting scripts were not preserved).
#   Figure 1: discrimination + calibration, Oxford temporal hold-out
#             test set (panels A/B) and UCLH external validation
#             cohort (panels C/D, pre-computed by UCLH; see below).
#   Figure 2: SHAP summary plot for the top 20 features.
#   Figure 3: decision curve analysis (A) and resource-neutral
#             reallocation swept across ranking strategies (B).
# Run from the repository root; requires 02_Models/Do_XGBoost_FullCohort.r
# and 02_Models/NBC_9_Apply_Model.r to have already been run, and access
# to the underlying IORD data extract.
# ============================================================

fig_dir <- "Output Files/Figures/Report"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Shared feature engineering ------------------------------
# Mirrors the pipeline in 02_Models/Do_XGBoost_FullCohort.r exactly.

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

engineer_features <- function(dt) {
  dt[, Sex_num := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]
  dt[is.na(SOCIAL_IMDPercentile), SOCIAL_IMDPercentile := median(SOCIAL_IMDPercentile, na.rm = TRUE)]
  dt[, AdmissionType_num := fcase(
    AdmissionMethod == "Emergency", 1L, AdmissionMethod == "Elective", 2L, default = 0L)]
  dt[, SourceType_num := fcase(
    CultureSourceCategory == "Peripheral", 1L,
    CultureSourceCategory %in% c("Central Line", "PICC", "Port/Portacath", "Arterial", "Femoral"), 2L,
    default = 0L)]
  for (col in c("Source_Urinary", "Source_Respiratory", "Source_Skin_Soft_Tissue", "Source_Other",
                "HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before",
                "Carbapenem_Before", "Any_Prior_Other_Pos", "Prior_Urine_Pos",
                "PositiveCulturesInPrior28Days")) {
    dt[is.na(get(col)), (col) := 0L]
  }
  dt[, VS_AnyFever := fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), 1, 0)]
  dt[is.na(LOS_BeforeCulture_days), LOS_BeforeCulture_days := 0]
  for (col in feature_cols) if (!col %in% names(dt)) dt[, (col) := NA_real_]
  dt
}

shapiro_score <- function(age, temp, sbp, wcc, plt, creat) {
  # Modified Shapiro rule (Shapiro*): intravascular-catheter term zeroed
  # throughout, as this variable is not reliably available at triage.
  as.integer(age > 65) +
    fcase(!is.na(temp) & temp >= 39.4, 3L, !is.na(temp) & temp >= 38.3, 1L, default = 0L) +
    fifelse(!is.na(sbp) & sbp < 90, 1L, 0L) +
    fifelse(!is.na(wcc) & wcc > 18, 1L, 0L) +
    fifelse(!is.na(plt) & plt < 150, 1L, 0L) +
    fifelse(!is.na(creat) & creat > 176.8, 1L, 0L)
}

bc_data <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))
bc_data[, CollectionDateTime := ymd_hms(CollectionDateTime)]
bc_data[, Outcome := fifelse(CultivabilityCategory %in%
  c("Gram-positive", "Gram-negative", "Polymicrobial", "Fungi", "Anaerobes", "Other"), 1L, 0L)]
bc_data <- bc_data[CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial",
  "Fungi", "Anaerobes", "Other", "Culture-negative", "Contaminant - single",
  "Contaminant - recurrent", "Polymicrobial - contaminant")]
bc_data <- engineer_features(bc_data)

model       <- xgb.load("Models/XGBoost/FullCohort/xgb_fullcohort.json")
platt_coefs <- readRDS("Models/XGBoost/FullCohort/xgb_fullcohort_platt_calibration.rds")

calibrate <- function(dt) {
  raw   <- predict(model, xgb.DMatrix(as.matrix(dt[, ..feature_cols])))
  logit <- qlogis(pmax(pmin(raw, 0.999), 0.001))
  plogis(platt_coefs[1] + platt_coefs[2] * logit)
}

# ============================================================
# FIGURE 1 - discrimination + calibration
# ============================================================

ox <- fread("Output Files/CSVs/XGBoost/FullCohort/test_predictions.csv")

panel_roc <- function(roc_df, auc_ci, title) {
  ggplot(roc_df, aes(FPR, TPR)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_line(colour = "#2166ac", linewidth = 0.9) +
    annotate("text", x = 0.65, y = 0.1,
             label = sprintf("AUROC %.3f (%.3f-%.3f)", auc_ci[2], auc_ci[1], auc_ci[3]),
             size = 3.4, hjust = 0) +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity") +
    coord_equal() + theme_bw(base_size = 10)
}

panel_calibration <- function(pred, outcome, deciles, title, rug_max) {
  fit  <- loess(outcome ~ pred, span = 0.75)
  grid <- data.frame(pred = seq(min(pred), max(pred), length.out = 200))
  pr   <- predict(fit, newdata = grid, se = TRUE)
  grid$fit <- pmin(pmax(pr$fit, 0), 1)
  grid$lo  <- pmin(pmax(pr$fit - 1.96 * pr$se.fit, 0), 1)
  grid$hi  <- pmin(pmax(pr$fit + 1.96 * pr$se.fit, 0), 1)

  ggplot() +
    geom_ribbon(data = grid, aes(pred, ymin = lo, ymax = hi), fill = "#c6dbef", alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_line(data = grid, aes(pred, fit), colour = "#2166ac", linewidth = 0.9) +
    geom_point(data = deciles, aes(mean_pred, obs_rate), colour = "#2166ac", size = 1.8) +
    geom_rug(data = data.frame(pred = pred[seq(1, length(pred), by = max(1, length(pred) %/% rug_max))]),
             aes(pred), sides = "b", alpha = 0.25, length = unit(0.02, "npc")) +
    labs(title = title, x = "Predicted probability", y = "Observed proportion") +
    coord_cartesian(xlim = c(0, max(grid$pred))) + theme_bw(base_size = 10)
}

ox_roc    <- roc(ox$outcome, ox$pred, quiet = TRUE)
ox_auc    <- ci.auc(ox_roc)
ox_roc_df <- data.frame(FPR = 1 - ox_roc$specificities, TPR = ox_roc$sensitivities)

cal_trunc <- quantile(ox$pred, 0.99)
ox_cal <- ox[pred <= cal_trunc]
ox_cal[, decile := cut(pred, breaks = quantile(pred, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)]
ox_dec <- ox_cal[, .(mean_pred = mean(pred), obs_rate = mean(outcome)), by = decile][order(decile)]

pA <- panel_roc(ox_roc_df, ox_auc, "A. Oxford temporal test - ROC")
pB <- panel_calibration(ox_cal$pred, ox_cal$outcome, ox_dec, "B. Oxford temporal test - calibration", 3000)

# UCLH external validation panels: pre-computed by UCLH applying the
# shared model to their cohort (see manuscript Methods); update this
# path to wherever those result files are placed locally.
ev_dir <- "Output Files/External Validation/CSVs"
figure1 <- pA | pB
if (file.exists(file.path(ev_dir, "ev_model_performance.csv"))) {
  ev_roc   <- fread(file.path(ev_dir, "ev_roc_coords.csv"))
  ev_dec   <- fread(file.path(ev_dir, "ev_calibration_deciles.csv"))
  ev_curve <- fread(file.path(ev_dir, "ev_calibration_curve.csv"))
  ev_perf  <- fread(file.path(ev_dir, "ev_model_performance.csv"))
  ev_auc   <- as.numeric(ev_perf[Metric == "AUROC", Value])
  ev_lo    <- as.numeric(ev_perf[Metric == "AUROC_CI_Lower", Value])
  ev_hi    <- as.numeric(ev_perf[Metric == "AUROC_CI_Upper", Value])

  pC <- panel_roc(ev_roc, c(ev_lo, ev_auc, ev_hi), "C. UCLH external validation - ROC")
  pD <- ggplot() +
    geom_ribbon(data = ev_curve, aes(pred, ymin = lo, ymax = hi), fill = "#c6dbef", alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_line(data = ev_curve, aes(pred, fit), colour = "#2166ac", linewidth = 0.9) +
    geom_point(data = ev_dec, aes(mean_pred, obs_rate), colour = "#2166ac", size = 1.8) +
    labs(title = "D. UCLH external validation - calibration",
         x = "Predicted probability", y = "Observed proportion") +
    coord_cartesian(xlim = c(0, max(ev_curve$pred))) + theme_bw(base_size = 10)

  figure1 <- (pA | pB) / (pC | pD)
}

ggsave(file.path(fig_dir, "Figure1_ROC_Calibration.png"), figure1, width = 9, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "Figure1_ROC_Calibration.pdf"), figure1, width = 9, height = 8)

# ============================================================
# FIGURE 2 - SHAP summary plot
# ============================================================

bc_data[, CollectionDate := as.Date(CollectionDateTime)]
shap_test <- bc_data[CollectionDate >= as.Date("2024-01-01")]
shap_matrix <- as.matrix(shap_test[, ..feature_cols])

contrib <- predict(model, shap_matrix, predcontrib = TRUE)
contrib <- contrib[, feature_cols, drop = FALSE]

top_n <- 20
mean_abs_shap <- sort(colMeans(abs(contrib)), decreasing = TRUE)
top_features  <- names(mean_abs_shap)[seq_len(min(top_n, length(mean_abs_shap)))]

shap_long <- rbindlist(lapply(top_features, function(f) {
  vals <- shap_matrix[, f]
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

# ============================================================
# FIGURE 3 - decision curve analysis (A) and reallocation (B)
# ============================================================

dca <- fread("Output Files/CSVs/XGBoost/FullCohort/dca_comparators.csv")
dca_long <- melt(dca, id.vars = "Threshold",
                  measure.vars = c("NB_XGBoost", "NB_Shapiro", "NB_NEWS", "NB_CRP", "NB_TreatAll", "NB_TreatNone"),
                  variable.name = "Strategy", value.name = "NetBenefit")
dca_long[, Strategy := factor(Strategy,
  levels = c("NB_XGBoost", "NB_Shapiro", "NB_NEWS", "NB_CRP", "NB_TreatAll", "NB_TreatNone"),
  labels = c("XGBoost", "Shapiro*", "NEWS", "CRP", "Treat all", "Treat none"))]

panelA <- ggplot(dca_long, aes(Threshold, NetBenefit, colour = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("XGBoost" = "#2166ac", "Shapiro*" = "#d6604d", "NEWS" = "#4daf4a",
                                  "CRP" = "#984ea3", "Treat all" = "grey50", "Treat none" = "grey80")) +
  scale_linetype_manual(values = c("XGBoost" = 1, "Shapiro*" = 1, "NEWS" = 1, "CRP" = 1,
                                    "Treat all" = 2, "Treat none" = 2)) +
  labs(title = "A. Decision curve analysis", x = "Threshold probability", y = "Net benefit") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.title = element_blank())

bc_2024 <- bc_data[!is.na(CollectionDateTime) & year(CollectionDateTime) == 2024]
bc_2024[, pred_calibrated := calibrate(bc_2024)]
bc_2024[, Shapiro := shapiro_score(AgeAtCulture, VS_Baseline_Temp, VS_Baseline_SBP,
                                    BT_Closest_WCC, BT_Closest_Plt, BT_Closest_Creat)]

nbc <- fread("Output Files/CSVs/NonBC_Admissions/2024/nonbc_admissions_with_predictions.csv",
             colClasses = c(ClusterID = "character"))
nbc <- engineer_features(nbc)
nbc[, pred_calibrated := calibrate(nbc)]
nbc[, Shapiro := shapiro_score(AgeAtCulture, VS_Baseline_Temp, VS_Baseline_SBP,
                                BT_Closest_WCC, BT_Closest_Plt, BT_Closest_Creat)]

n_bc <- nrow(bc_2024)
swap_sizes <- unique(pmin(c(500, 1000, 2000, 3000, 5000, 7500, 10000, 15000, 20000), n_bc))

sweep_strategy <- function(score_bc, score_nbc, label) {
  ord_bc   <- order(score_bc)     # lowest risk sent-culture first (culled)
  ord_nbc  <- order(-score_nbc)   # highest risk untested-admission first (added)
  base_pos <- sum(bc_2024$Outcome)
  rbindlist(lapply(swap_sizes, function(k) {
    culled <- bc_2024[ord_bc[seq_len(k)]]
    added  <- nbc[ord_nbc[seq_len(k)]]
    lost   <- sum(culled$Outcome)
    gained <- sum(added$pred_calibrated)  # XGBoost risk used as outcome proxy throughout
    data.table(Strategy = label, Swaps = k,
               RelIncreasePct = 100 * (gained - lost) / base_pos)
  }))
}

realloc <- rbindlist(list(
  sweep_strategy(bc_2024$pred_calibrated,  nbc$pred_calibrated,  "XGBoost"),
  sweep_strategy(bc_2024$Shapiro,          nbc$Shapiro,          "Shapiro*"),
  sweep_strategy(bc_2024$VS_Baseline_NEWS, nbc$VS_Baseline_NEWS, "NEWS"),
  sweep_strategy(bc_2024$BT_Closest_CRP,   nbc$BT_Closest_CRP,   "CRP")
))

panelB <- ggplot(realloc, aes(Swaps, RelIncreasePct, colour = Strategy)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  geom_vline(xintercept = n_bc, linetype = 3, colour = "grey40") +
  scale_x_log10() +
  scale_colour_manual(values = c("XGBoost" = "#2166ac", "Shapiro*" = "#d6604d",
                                  "NEWS" = "#4daf4a", "CRP" = "#984ea3")) +
  labs(title = "B. Resource-neutral reallocation", x = "Swaps (log scale)",
       y = "Relative increase in positive cultures (%)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.title = element_blank())

figure3 <- panelA + panelB

ggsave(file.path(fig_dir, "Figure3_DCA_Reallocation.png"), figure3, width = 11, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "Figure3_DCA_Reallocation.pdf"), figure3, width = 11, height = 5.5)
fwrite(realloc, "Output Files/CSVs/NonBC_Admissions/2024/reallocation_by_strategy.csv")
