library(data.table)
library(xgboost)
library(lubridate)
library(ggplot2)
library(patchwork)

# ============================================================
# FIGURE 3 - Comparison with existing risk-stratification tools
#            on the observed 2024 blood culture cohort.
# Reconstructed from saved model + pipeline outputs (the original
# plotting script was not preserved).
#   Panel A: decision curve analysis (net benefit vs threshold)
#            for XGBoost, modified Shapiro rule, NEWS and CRP.
#   Panel B: resource-neutral reallocation swept over swap volumes,
#            ranking sent 2024 cultures vs untested 2024 admissions
#            by each strategy's score. XGBoost-predicted risk is
#            used as the outcome proxy for untested admissions
#            (ground truth unavailable), as described in the paper.
# Run from the repository root; requires 02_Models/Do_XGBoost_FullCohort.r
# and 02_Models/NBC_9_Apply_Model.r to have already been run.
# ============================================================

fig_dir <- "Output Files/Figures/Report"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Panel A: decision curve analysis ------------------------
# NB_XGBoost/NB_CRP/NB_NEWS/NB_Shapiro/NB_TreatAll/NB_TreatNone
# computed in 02_Models/Do_XGBoost_FullCohort.r

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

# ---- Panel B: reallocation swept across ranking strategies ----

bc_data <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))
bc_data[, CollectionDateTime := ymd_hms(CollectionDateTime)]
bc_data[, Outcome := fifelse(CultivabilityCategory %in%
  c("Gram-positive", "Gram-negative", "Polymicrobial", "Fungi", "Anaerobes", "Other"), 1L, 0L)]

bc_data[, Sex_num := fifelse(DEMOGRAPHICS_LinkedSex == "M", 1, 0)]
bc_data[, SOCIAL_IMDPercentile := fifelse(is.na(SOCIAL_IMDPercentile),
  median(SOCIAL_IMDPercentile, na.rm = TRUE), SOCIAL_IMDPercentile)]
bc_data[, AdmissionType_num := fcase(
  AdmissionMethod == "Emergency", 1L, AdmissionMethod == "Elective", 2L, default = 0L)]
bc_data[, SourceType_num := fcase(
  CultureSourceCategory == "Peripheral", 1L,
  CultureSourceCategory %in% c("Central Line", "PICC", "Port/Portacath", "Arterial", "Femoral"), 2L,
  default = 0L)]
for (col in c("Source_Urinary", "Source_Respiratory", "Source_Skin_Soft_Tissue", "Source_Other",
              "HasAntibiotics_Before", "BroadSpectrum_Before", "Ceph3rdGenPlus_Before",
              "Carbapenem_Before", "Any_Prior_Other_Pos", "Prior_Urine_Pos",
              "PositiveCulturesInPrior28Days")) {
  bc_data[is.na(get(col)), (col) := 0L]
}
bc_data[, VS_AnyFever := fifelse(VS_AnyFever_38 == 1 & !is.na(VS_AnyFever_38), 1, 0)]
bc_data[is.na(LOS_BeforeCulture_days), LOS_BeforeCulture_days := 0]

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
for (col in feature_cols) if (!col %in% names(bc_data)) bc_data[, (col) := NA_real_]

bc_2024 <- bc_data[!is.na(CollectionDateTime) & year(CollectionDateTime) == 2024]
rm(bc_data)

model       <- xgb.load("Models/XGBoost/FullCohort/xgb_fullcohort.json")
platt_coefs <- readRDS("Models/XGBoost/FullCohort/xgb_fullcohort_platt_calibration.rds")
pred_raw    <- predict(model, xgb.DMatrix(as.matrix(bc_2024[, ..feature_cols])))
pred_logit  <- qlogis(pmax(pmin(pred_raw, 0.999), 0.001))
bc_2024[, pred_calibrated := plogis(platt_coefs[1] + platt_coefs[2] * pred_logit)]

# Modified Shapiro rule (Shapiro*): intravascular-catheter term zeroed
# throughout, as this variable is not reliably available at triage.
shapiro_score <- function(age, temp, sbp, wcc, plt, creat) {
  as.integer(age > 65) +
    fcase(!is.na(temp) & temp >= 39.4, 3L, !is.na(temp) & temp >= 38.3, 1L, default = 0L) +
    fifelse(!is.na(sbp) & sbp < 90, 1L, 0L) +
    fifelse(!is.na(wcc) & wcc > 18, 1L, 0L) +
    fifelse(!is.na(plt) & plt < 150, 1L, 0L) +
    fifelse(!is.na(creat) & creat > 176.8, 1L, 0L)
}
bc_2024[, Shapiro := shapiro_score(AgeAtCulture, VS_Baseline_Temp, VS_Baseline_SBP,
                                    BT_Closest_WCC, BT_Closest_Plt, BT_Closest_Creat)]

nbc <- fread("Output Files/CSVs/NonBC_Admissions/2024/nonbc_admissions_with_predictions.csv",
             colClasses = c(ClusterID = "character"))
nbc[, Shapiro := shapiro_score(AgeAtCulture, VS_Baseline_Temp, VS_Baseline_SBP,
                                BT_Closest_WCC, BT_Closest_Plt, BT_Closest_Creat)]

n_bc  <- nrow(bc_2024)
swap_sizes <- unique(pmin(c(500, 1000, 2000, 3000, 5000, 7500, 10000, 15000, 20000), n_bc))

sweep_strategy <- function(score_bc, score_nbc, label) {
  ord_bc  <- order(score_bc)          # lowest risk sent-culture first (culled)
  ord_nbc <- order(-score_nbc)        # highest risk untested-admission first (added)
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
