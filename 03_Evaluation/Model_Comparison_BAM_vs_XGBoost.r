################################################################################
# Formal Statistical Comparison: XGBoost vs BAM (Full Cohort)
#
# AUROC comparison : DeLong test (pROC::roc.test) — accounts for pairing
# AUPRC comparison : paired bootstrap — same resampled indices for both models
# Brier comparison : paired bootstrap on squared errors
################################################################################

library(data.table)
library(pROC)
library(PRROC)
library(mgcv)
library(xgboost)
library(ggplot2)
library(writexl)

set.seed(42)
N_BOOT <- 2000

out_dir <- "Output Files/Publication_Tables_Figures/Metrics"
fig_dir <- "Output Files/Figures/Comparisons"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# 1. LOAD AND PREPARE TEST DATA  (shared pipeline)
################################################################################

bc_data <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))

bc_data <- bc_data[CultivabilityCategory %in% c(
  "Gram-positive", "Gram-negative", "Polymicrobial",
  "Fungi", "Anaerobes", "Other",
  "Culture-negative", "Contaminant - single", "Contaminant - recurrent",
  "Polymicrobial - contaminant"
)]
bc_data[, Outcome := fifelse(CultivabilityCategory %in%
                               c("Gram-positive", "Gram-negative", "Polymicrobial",
                                 "Fungi", "Anaerobes", "Other"), 1L, 0L)]

# ── Shared feature engineering ──────────────────────────────────────────────

median_imd <- median(bc_data$SOCIAL_IMDPercentile, na.rm = TRUE)
bc_data[is.na(SOCIAL_IMDPercentile), SOCIAL_IMDPercentile := median_imd]
bc_data[is.na(PositiveCulturesInPrior28Days), PositiveCulturesInPrior28Days := 0]
bc_data[is.na(LOS_BeforeCulture_days), LOS_BeforeCulture_days := 0]

for (col in c("Source_Urinary", "Source_Respiratory", "Source_Skin_Soft_Tissue", "Source_Other",
              "HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before", "Carbapenem_Before",
              "Any_Prior_Other_Pos", "Prior_Urine_Pos")) {
  bc_data[is.na(get(col)), (col) := 0]
}

# BAM-style factors
bc_data[, Sex_f := factor(fifelse(is.na(DEMOGRAPHICS_LinkedSex), "Unknown", DEMOGRAPHICS_LinkedSex),
                          levels = c("F", "M", "Unknown"))]
bc_data[, AdmissionType_f := factor(fcase(
  is.na(AdmissionMethod), "Unknown",
  AdmissionMethod == "Emergency", "Emergency",
  AdmissionMethod == "Elective", "Elective",
  default = "Other"
), levels = c("Other", "Emergency", "Elective", "Unknown"))]
bc_data[, SourceType_f := factor(fcase(
  is.na(CultureSourceCategory), "Unknown",
  CultureSourceCategory == "Peripheral", "Peripheral",
  CultureSourceCategory %in% c("Central Line", "PICC", "Port/Portacath", "Arterial", "Femoral"), "Central",
  default = "Other"
), levels = c("Other", "Peripheral", "Central", "Unknown"))]

for (col in c("Source_Urinary", "Source_Respiratory", "Source_Skin_Soft_Tissue", "Source_Other")) {
  bc_data[, paste0(col, "_f") := factor(fifelse(get(col) == 1, "Yes", "No"))]
}
for (col in c("HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before", "Carbapenem_Before",
              "Any_Prior_Other_Pos", "Prior_Urine_Pos")) {
  bc_data[, paste0(col, "_f") := factor(fifelse(get(col) == 1, "Yes", "No"))]
}
bc_data[, VS_AnyFever_f   := factor(fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), "Yes", "No"))]
bc_data[, VS_AnyHypotension_f := factor(fifelse(VS_AnyHypotension == 1 & !is.na(VS_AnyHypotension), "Yes", "No"))]

# XGBoost-style numerics
bc_data[, Sex_num := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]
bc_data[, AdmissionType_num := fcase(AdmissionMethod == "Emergency", 1,
                                      AdmissionMethod == "Elective",  2,
                                      default = 0)]
bc_data[, SourceType_num := fcase(
  CultureSourceCategory == "Peripheral", 1,
  CultureSourceCategory %in% c("Central Line", "PICC", "Port/Portacath", "Arterial", "Femoral"), 2,
  default = 0)]
bc_data[, VS_AnyFever     := fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), 1, 0)]
bc_data[, VS_AnyHypotension := fifelse(VS_AnyHypotension == 1 & !is.na(VS_AnyHypotension), 1, 0)]

# Snapshot for XGBoost BEFORE continuous imputation: flags/LOS/prior-cultures are zero-filled
# (exactly as the headline model does) but labs/vitals stay native NA (XGBoost handles missing
# internally). BAM below uses the median-imputed copy since it cannot handle NA. DeLong only needs
# paired predictions over the same test rows, not identical input features.
bc_data[, .row_id := .I]
xgb_native <- copy(bc_data)

# Median imputation for continuous vars (BAM pipeline)
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

# ── Train/test split ─────────────────────────────────────────────────────────
bc_data[, CollectionDate := as.Date(CollectionDateTime)]
split_date  <- as.Date("2024-01-01")
test_data   <- bc_data[CollectionDate >= split_date]
test_label  <- test_data$Outcome

################################################################################
# 2. GET CALIBRATED PREDICTIONS FROM BOTH MODELS
################################################################################

model_dir_xgb <- "Models/XGBoost/FullCohort"
model_dir_bam <- "Models/BAM/FullCohort"

# ── XGBoost ──────────────────────────────────────────────────────────────────
xgb_model    <- xgb.load(file.path(model_dir_xgb, "xgb_fullcohort.json"))
xgb_platt    <- readRDS(file.path(model_dir_xgb, "xgb_fullcohort_platt_calibration.rds"))

feature_cols <- c(
  "AgeAtCulture", "Sex_num",
  "AdmissionType_num", "LOS_BeforeCulture_days",
  "SourceType_num", "PositiveCulturesInPrior28Days",
  "Source_Respiratory",
  "HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before", "Carbapenem_Before",
  "Any_Prior_Other_Pos", "Prior_Urine_Pos",
  "VS_Baseline_Temp", "VS_Baseline_SBP",
  "VS_Baseline_NEWS", "VS_AnyFever", "VS_Baseline_HR",
  "BT_Closest_CRP", "BT_Closest_WCC", "BT_Closest_Neut", "BT_Closest_Lymph",
  "BT_Closest_Alb", "BT_Closest_Plt", "BT_Closest_Creat", "BT_Closest_Bili",
  "BT_Closest_Na", "BT_Closest_Eosin", "BT_Closest_Mono", "BT_Closest_ALP",
  "BG_Closest_Lactate", "BG_Closest_pH"
)

test_matrix    <- as.matrix(xgb_native[.row_id %in% test_data$.row_id][, ..feature_cols])
xgb_raw        <- predict(xgb_model, xgb.DMatrix(data = test_matrix))
xgb_logit      <- qlogis(pmax(pmin(xgb_raw, 0.999), 0.001))
pred_xgb       <- plogis(xgb_platt[1] + xgb_platt[2] * xgb_logit)

# ── BAM ───────────────────────────────────────────────────────────────────────
bam_model      <- readRDS(file.path(model_dir_bam, "bam_fullcohort_model.rds"))
bam_platt      <- readRDS(file.path(model_dir_bam, "bam_fullcohort_platt_calibration.rds"))

test_df        <- as.data.frame(test_data)
bam_raw        <- predict(bam_model, newdata = test_df, type = "response")
bam_logit      <- qlogis(pmax(pmin(bam_raw, 0.999), 0.001))
pred_bam       <- predict(bam_platt, newdata = data.frame(train_pred_logit = bam_logit), type = "response")

################################################################################
# 3. ROC OBJECTS & DELONG TEST FOR AUROC
################################################################################

roc_xgb <- roc(test_label, pred_xgb, quiet = TRUE)
roc_bam <- roc(test_label, pred_bam, quiet = TRUE)

# DeLong — correct for paired predictions (same test set, correlated)
delong <- roc.test(roc_xgb, roc_bam, method = "delong", paired = TRUE)

auroc_xgb    <- as.numeric(auc(roc_xgb))
auroc_bam    <- as.numeric(auc(roc_bam))
ci_xgb       <- ci.auc(roc_xgb)
ci_bam       <- ci.auc(roc_bam)

################################################################################
# 4. PAIRED BOOTSTRAP FOR AUPRC & BRIER
################################################################################

auprc_point <- function(pred, label) {
  valid <- !is.na(pred)
  if (sum(label[valid] == 1) < 1 || sum(label[valid] == 0) < 1) return(NA_real_)
  tryCatch(
    pr.curve(scores.class0 = pred[valid & label == 1],
             scores.class1 = pred[valid & label == 0])$auc.integral,
    error = function(e) NA_real_
  )
}

auprc_xgb_obs <- auprc_point(pred_xgb, test_label)
auprc_bam_obs <- auprc_point(pred_bam, test_label)

brier_xgb_obs <- mean((pred_xgb - test_label)^2)
brier_bam_obs <- mean((pred_bam - test_label)^2)

n_test <- length(test_label)

boot_results <- vapply(seq_len(N_BOOT), function(b) {
  idx <- sample(n_test, replace = TRUE)
  bx  <- pred_xgb[idx];  bb  <- pred_bam[idx];  bl  <- test_label[idx]

  a_xgb <- auprc_point(bx, bl)
  a_bam <- auprc_point(bb, bl)

  brier_xgb <- mean((bx - bl)^2)
  brier_bam <- mean((bb - bl)^2)

  c(auprc_xgb  = a_xgb,
    auprc_bam  = a_bam,
    auprc_diff = a_xgb - a_bam,     # positive = XGBoost better
    brier_xgb  = brier_xgb,
    brier_bam  = brier_bam,
    brier_diff = brier_xgb - brier_bam)  # negative = XGBoost better (lower Brier)
}, numeric(6))

boot_dt <- as.data.table(t(boot_results))

ci_auprc_xgb  <- quantile(boot_dt$auprc_xgb,  c(0.025, 0.975), na.rm = TRUE)
ci_auprc_bam  <- quantile(boot_dt$auprc_bam,  c(0.025, 0.975), na.rm = TRUE)
ci_auprc_diff <- quantile(boot_dt$auprc_diff, c(0.025, 0.975), na.rm = TRUE)
p_auprc       <- mean(boot_dt$auprc_diff <= 0, na.rm = TRUE)  # 1-sided: XGBoost >= BAM
p_auprc_2s    <- 2 * min(p_auprc, 1 - p_auprc)               # 2-sided

ci_brier_xgb  <- quantile(boot_dt$brier_xgb,  c(0.025, 0.975), na.rm = TRUE)
ci_brier_bam  <- quantile(boot_dt$brier_bam,  c(0.025, 0.975), na.rm = TRUE)
ci_brier_diff <- quantile(boot_dt$brier_diff, c(0.025, 0.975), na.rm = TRUE)
p_brier       <- mean(boot_dt$brier_diff >= 0, na.rm = TRUE)  # 1-sided: XGBoost >= BAM
p_brier_2s    <- 2 * min(p_brier, 1 - p_brier)

################################################################################
# 5. RESULTS TABLE
################################################################################

fmt_ci <- function(lo, hi) sprintf("(%.3f–%.3f)", lo, hi)

comparison_table <- data.frame(
  Metric      = c("AUROC", "AUPRC", "Brier Score"),
  XGBoost     = c(sprintf("%.3f %s", auroc_xgb, fmt_ci(ci_xgb[1], ci_xgb[3])),
                  sprintf("%.3f %s", auprc_xgb_obs, fmt_ci(ci_auprc_xgb[1], ci_auprc_xgb[2])),
                  sprintf("%.4f %s", brier_xgb_obs, fmt_ci(ci_brier_xgb[1], ci_brier_xgb[2]))),
  BAM         = c(sprintf("%.3f %s", auroc_bam, fmt_ci(ci_bam[1], ci_bam[3])),
                  sprintf("%.3f %s", auprc_bam_obs, fmt_ci(ci_auprc_bam[1], ci_auprc_bam[2])),
                  sprintf("%.4f %s", brier_bam_obs, fmt_ci(ci_brier_bam[1], ci_brier_bam[2]))),
  Difference_XGB_minus_BAM = c(
    sprintf("%.3f %s", auroc_xgb - auroc_bam,
            fmt_ci(delong$conf.int[1], delong$conf.int[2])),
    sprintf("%.3f %s", auprc_xgb_obs - auprc_bam_obs,
            fmt_ci(ci_auprc_diff[1], ci_auprc_diff[2])),
    sprintf("%.4f %s", brier_xgb_obs - brier_bam_obs,
            fmt_ci(ci_brier_diff[1], ci_brier_diff[2]))
  ),
  P_value     = c(sprintf("%.4f", delong$p.value),
                  sprintf("%.4f", p_auprc_2s),
                  sprintf("%.4f", p_brier_2s)),
  Method      = c("DeLong (paired)", "Paired bootstrap", "Paired bootstrap"),
  stringsAsFactors = FALSE
)

fwrite(comparison_table, file.path(out_dir, "model_comparison_xgboost_vs_bam.csv"))

# Machine-readable version for downstream scripts
comparison_raw <- data.frame(
  Metric             = c("AUROC", "AUPRC", "Brier"),
  XGBoost_Point      = c(auroc_xgb, auprc_xgb_obs, brier_xgb_obs),
  XGBoost_CI_Lo      = c(ci_xgb[1], ci_auprc_xgb[1], ci_brier_xgb[1]),
  XGBoost_CI_Hi      = c(ci_xgb[3], ci_auprc_xgb[2], ci_brier_xgb[2]),
  BAM_Point          = c(auroc_bam, auprc_bam_obs, brier_bam_obs),
  BAM_CI_Lo          = c(ci_bam[1], ci_auprc_bam[1], ci_brier_bam[1]),
  BAM_CI_Hi          = c(ci_bam[3], ci_auprc_bam[2], ci_brier_bam[2]),
  Diff_Point         = c(auroc_xgb - auroc_bam, auprc_xgb_obs - auprc_bam_obs, brier_xgb_obs - brier_bam_obs),
  Diff_CI_Lo         = c(delong$conf.int[1], ci_auprc_diff[1], ci_brier_diff[1]),
  Diff_CI_Hi         = c(delong$conf.int[2], ci_auprc_diff[2], ci_brier_diff[2]),
  P_value            = c(delong$p.value, p_auprc_2s, p_brier_2s),
  stringsAsFactors   = FALSE
)
fwrite(comparison_raw, file.path(out_dir, "model_comparison_raw.csv"))

write_xlsx(
  list(Formatted = comparison_table, Raw = comparison_raw),
  file.path(out_dir, "model_comparison_xgboost_vs_bam.xlsx")
)

################################################################################
# 6. FIGURE — Forest-style comparison plot
################################################################################

# Panel A: point estimates + 95% CIs by metric
plot_df <- data.table(
  Metric  = rep(c("AUROC", "AUPRC", "Brier Score"), each = 2),
  Model   = rep(c("XGBoost", "BAM"), 3),
  Point   = c(auroc_xgb, auroc_bam, auprc_xgb_obs, auprc_bam_obs, brier_xgb_obs, brier_bam_obs),
  Lo      = c(ci_xgb[1], ci_bam[1], ci_auprc_xgb[1], ci_auprc_bam[1], ci_brier_xgb[1], ci_brier_bam[1]),
  Hi      = c(ci_xgb[3], ci_bam[3], ci_auprc_xgb[2], ci_auprc_bam[2], ci_brier_xgb[2], ci_brier_bam[2])
)
plot_df[, Model := factor(Model, levels = c("XGBoost", "BAM"))]
plot_df[, Metric := factor(Metric, levels = c("AUROC", "AUPRC", "Brier Score"))]

# For Brier, lower is better — flag it
plot_df[, Direction := ifelse(Metric == "Brier Score", "lower_better", "higher_better")]

# P-value label positions (top of each facet)
pval_df <- data.table(
  Metric  = factor(c("AUROC", "AUPRC", "Brier Score"),
                   levels = c("AUROC", "AUPRC", "Brier Score")),
  label   = c(sprintf("DeLong p = %s", formatC(delong$p.value, format = "f", digits = 4)),
              sprintf("Bootstrap p = %s", formatC(p_auprc_2s,  format = "f", digits = 4)),
              sprintf("Bootstrap p = %s", formatC(p_brier_2s,  format = "f", digits = 4))),
  y       = c(max(ci_xgb[3], ci_bam[3]) * 1.01,
              max(ci_auprc_xgb[2], ci_auprc_bam[2]) * 1.01,
              max(ci_brier_xgb[2], ci_brier_bam[2]) * 1.01)
)

p_compare <- ggplot(plot_df, aes(x = Model, y = Point, color = Model)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lo, ymax = Hi), width = 0.2, linewidth = 1) +
  geom_text(data = pval_df, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = 3.2, fontface = "italic", hjust = 0.5) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("XGBoost" = "#2E86AB", "BAM" = "#E94F37")) +
  labs(
    title    = "XGBoost vs BAM: Test-Set Performance Comparison",
    subtitle = sprintf("Test set  |  AUROC: DeLong paired test  |  AUPRC & Brier: paired bootstrap (B = %d)", N_BOOT),
    x = "", y = "Metric Value", color = "Model"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    legend.position = "bottom",
    strip.text    = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "xgboost_vs_bam_comparison.png"),
       p_compare, width = 11, height = 5.5, dpi = 300)

# Panel B: bootstrap distribution of AUPRC difference
p_boot_diff <- ggplot(boot_dt[!is.na(auprc_diff)], aes(x = auprc_diff)) +
  geom_histogram(bins = 80, fill = "#2E86AB", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 1) +
  geom_vline(xintercept = auprc_xgb_obs - auprc_bam_obs,
             color = "#E94F37", linewidth = 1.2) +
  geom_vline(xintercept = ci_auprc_diff, linetype = "dotted",
             color = "#E94F37", linewidth = 0.9) +
  annotate("text", x = auprc_xgb_obs - auprc_bam_obs, y = Inf,
           label = sprintf(" Observed\n Δ = %.3f", auprc_xgb_obs - auprc_bam_obs),
           vjust = 1.5, hjust = 0, color = "#E94F37", size = 3.5) +
  labs(
    title    = "Bootstrap Distribution: AUPRC Difference (XGBoost − BAM)",
    subtitle = sprintf("95%% CI: (%.3f, %.3f)  |  2-sided p = %.4f",
                       ci_auprc_diff[1], ci_auprc_diff[2], p_auprc_2s),
    x = "AUPRC(XGBoost) − AUPRC(BAM)", y = "Bootstrap iterations"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "xgboost_vs_bam_auprc_boot_diff.png"),
       p_boot_diff, width = 9, height = 6, dpi = 300)

# Panel C: bootstrap distribution of AUROC difference (from boot samples, not DeLong)
boot_auroc <- vapply(seq_len(N_BOOT), function(b) {
  idx <- sample(n_test, replace = TRUE)
  rx  <- tryCatch(as.numeric(auc(roc(test_label[idx], pred_xgb[idx], quiet = TRUE))), error = function(e) NA_real_)
  rb  <- tryCatch(as.numeric(auc(roc(test_label[idx], pred_bam[idx],  quiet = TRUE))), error = function(e) NA_real_)
  rx - rb
}, numeric(1))

ci_auroc_diff_boot <- quantile(boot_auroc, c(0.025, 0.975), na.rm = TRUE)
p_auroc_boot       <- 2 * min(mean(boot_auroc <= 0, na.rm = TRUE),
                               mean(boot_auroc >= 0, na.rm = TRUE))

p_boot_auroc <- ggplot(data.frame(diff = boot_auroc[!is.na(boot_auroc)]), aes(x = diff)) +
  geom_histogram(bins = 80, fill = "#2E86AB", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 1) +
  geom_vline(xintercept = auroc_xgb - auroc_bam,
             color = "#E94F37", linewidth = 1.2) +
  geom_vline(xintercept = ci_auroc_diff_boot, linetype = "dotted",
             color = "#E94F37", linewidth = 0.9) +
  annotate("text", x = auroc_xgb - auroc_bam, y = Inf,
           label = sprintf(" Observed\n Δ = %.3f", auroc_xgb - auroc_bam),
           vjust = 1.5, hjust = 0, color = "#E94F37", size = 3.5) +
  labs(
    title    = "Bootstrap Distribution: AUROC Difference (XGBoost − BAM)",
    subtitle = sprintf("95%% CI: (%.3f, %.3f)  |  2-sided p = %.4f  (cf. DeLong p = %.4f)",
                       ci_auroc_diff_boot[1], ci_auroc_diff_boot[2],
                       p_auroc_boot, delong$p.value),
    x = "AUROC(XGBoost) − AUROC(BAM)", y = "Bootstrap iterations"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "xgboost_vs_bam_auroc_boot_diff.png"),
       p_boot_auroc, width = 9, height = 6, dpi = 300)
