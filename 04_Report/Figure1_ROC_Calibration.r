library(data.table)
library(ggplot2)
library(pROC)
library(patchwork)

# ============================================================
# FIGURE 1 - Discrimination and calibration in temporal and
#            external validation.
# Reconstructed from saved model outputs (the original plotting
# script was not preserved). Panels A/B use the Oxford temporal
# hold-out test set; panels C/D use the UCLH external validation
# cohort (pre-computed by 05_Deployment/external_validation_package).
# Run from the repository root.
# ============================================================

fig_dir <- "Output Files/Figures/Report"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Oxford temporal hold-out test set ----------------------

ox <- fread("Output Files/CSVs/XGBoost/FullCohort/test_predictions.csv")

ox_roc  <- roc(ox$outcome, ox$pred, quiet = TRUE)
ox_auc  <- ci.auc(ox_roc)
ox_roc_df <- data.frame(FPR = 1 - ox_roc$specificities, TPR = ox_roc$sensitivities)

cal_trunc <- quantile(ox$pred, 0.99)
ox_cal <- ox[pred <= cal_trunc]
ox_cal[, decile := cut(pred, breaks = quantile(pred, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)]
ox_dec <- ox_cal[, .(mean_pred = mean(pred), obs_rate = mean(outcome)), by = decile][order(decile)]

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
  fit <- loess(outcome ~ pred, span = 0.75)
  grid <- data.frame(pred = seq(min(pred), max(pred), length.out = 200))
  pr <- predict(fit, newdata = grid, se = TRUE)
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

pA <- panel_roc(ox_roc_df, ox_auc, "A. Oxford temporal test - ROC")
pB <- panel_calibration(ox_cal$pred, ox_cal$outcome, ox_dec, "B. Oxford temporal test - calibration", 3000)

# ---- UCLH external validation cohort -------------------------
# Pre-computed by the shared model package run at UCLH (Apply_XGBoost.r
# / EV_Results_Analysis.r); update this path to the latest results folder.

ev_dir <- "05_Deployment/external_validation_package/UCLH_Results"
ev_roc  <- fread(file.path(ev_dir, "ev_roc_coords.csv"))
ev_dec  <- fread(file.path(ev_dir, "ev_calibration_deciles.csv"))
ev_curve <- fread(file.path(ev_dir, "ev_calibration_curve.csv"))
ev_perf  <- fread(file.path(ev_dir, "ev_model_performance.csv"))
ev_auc   <- as.numeric(ev_perf[Metric == "AUROC", Value])
ev_auc_lo <- as.numeric(ev_perf[Metric == "AUROC_CI_Lower", Value])
ev_auc_hi <- as.numeric(ev_perf[Metric == "AUROC_CI_Upper", Value])

pC <- panel_roc(ev_roc, c(ev_auc_lo, ev_auc, ev_auc_hi), "C. UCLH external validation - ROC")

pD <- ggplot() +
  geom_ribbon(data = ev_curve, aes(pred, ymin = lo, ymax = hi), fill = "#c6dbef", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_line(data = ev_curve, aes(pred, fit), colour = "#2166ac", linewidth = 0.9) +
  geom_point(data = ev_dec, aes(mean_pred, obs_rate), colour = "#2166ac", size = 1.8) +
  labs(title = "D. UCLH external validation - calibration",
       x = "Predicted probability", y = "Observed proportion") +
  coord_cartesian(xlim = c(0, max(ev_curve$pred))) + theme_bw(base_size = 10)

figure1 <- (pA | pB) / (pC | pD)

ggsave(file.path(fig_dir, "Figure1_ROC_Calibration_Composite.png"), figure1, width = 9, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "Figure1_ROC_Calibration_Composite.pdf"), figure1, width = 9, height = 8)
