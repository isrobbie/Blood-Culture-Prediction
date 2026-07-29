################################################################################
# Data Completeness by Outcome - All Model Variables
# Computes % completeness for each variable / domain in the final model,
# split by whole cohort, blood culture positives, and negatives.
################################################################################

library(data.table)
library(writexl)

bc_data <- fread("Output Files/CSVs/blood_cultures_with_other_micro.csv",
                 colClasses = c(ClusterID = "character"))

bc_data <- bc_data[CultivabilityCategory %in% c(
  "Gram-positive", "Gram-negative", "Polymicrobial",
  "Culture-negative", "Contaminant - single", "Contaminant - recurrent",
  "Polymicrobial - contaminant"
)]

bc_data[, Outcome := fifelse(
  CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial"), 1L, 0L
)]

pos <- bc_data[Outcome == 1]
neg <- bc_data[Outcome == 0]

n_all <- nrow(bc_data)
n_pos <- nrow(pos)
n_neg <- nrow(neg)

# Helper: non-NA count + pct for a column
completeness <- function(dt, col) {
  v <- dt[[col]]
  n <- sum(!is.na(v) & v != "")
  list(n = n, pct = round(100 * n / nrow(dt), 1))
}

# Infection source: complete = at least one source flag = 1
source_completeness <- function(dt) {
  n <- sum(dt$Source_Urinary == 1 | dt$Source_Respiratory == 1 |
           dt$Source_Skin_Soft_Tissue == 1 | dt$Source_Other == 1, na.rm = TRUE)
  list(n = n, pct = round(100 * n / nrow(dt), 1))
}

# Build variable-level table
# Each entry: Domain, Label, raw column name
var_defs <- list(
  # Demographics
  list(domain = "Demographics",       label = "Age at Culture",              col = "AgeAtCulture"),
  list(domain = "Demographics",       label = "Sex",                         col = "DEMOGRAPHICS_LinkedSex"),
  list(domain = "Demographics",       label = "IMD Percentile",              col = "SOCIAL_IMDPercentile"),
  # Admission & Presentation
  list(domain = "Admission",          label = "Admission Method",            col = "AdmissionMethod"),
  list(domain = "Admission",          label = "LOS Before Culture (days)",   col = "LOS_BeforeCulture_days"),
  # Culture Source
  list(domain = "Culture Source",     label = "Culture Source Category",     col = "CultureSourceCategory"),
  # Prior Blood Cultures
  list(domain = "Prior Blood Cultures", label = "Positive Cultures (Prior 28d)", col = "PositiveCulturesInPrior28Days"),
  list(domain = "Prior Blood Cultures", label = "Gram+ Cultures (Prior 28d)",    col = "GramPositiveCulturesInPrior28Days"),
  list(domain = "Prior Blood Cultures", label = "Gram- Cultures (Prior 28d)",    col = "GramNegativeCulturesInPrior28Days"),
  # Infection Source (complete = at least one source documented, not all-zero)
  list(domain = "Infection Source",   label = "Any Source Documented",       col = "__source_any__"),
  list(domain = "Infection Source",   label = "Source: Urinary",             col = "Source_Urinary"),
  list(domain = "Infection Source",   label = "Source: Respiratory",         col = "Source_Respiratory"),
  list(domain = "Infection Source",   label = "Source: Skin/Soft Tissue",    col = "Source_Skin_Soft_Tissue"),
  list(domain = "Infection Source",   label = "Source: Other",               col = "Source_Other"),
  # Antibiotics (binary flags; NA treated as 0 in model)
  list(domain = "Prior Antibiotics",  label = "Any Antibiotics (Before)",    col = "HasAntibiotics_Before"),
  list(domain = "Prior Antibiotics",  label = "Broad Spectrum (Before)",     col = "BroadSpectrum_Before"),
  list(domain = "Prior Antibiotics",  label = "Ceph 3rd Gen+ (Before)",      col = "Ceph3rdGenPlus_Before"),
  list(domain = "Prior Antibiotics",  label = "Carbapenem (Before)",         col = "Carbapenem_Before"),
  # Prior Organisms
  list(domain = "Prior Organisms",    label = "Any Prior Other Organism",    col = "Any_Prior_Other_Pos"),
  list(domain = "Prior Organisms",    label = "Prior Urine Positive",        col = "Prior_Urine_Pos"),
  # Vital Signs
  list(domain = "Vital Signs",        label = "Heart Rate",                  col = "VS_Baseline_HR"),
  list(domain = "Vital Signs",        label = "Temperature",                 col = "VS_Baseline_Temp"),
  list(domain = "Vital Signs",        label = "Systolic BP",                 col = "VS_Baseline_SBP"),
  list(domain = "Vital Signs",        label = "NEWS Score",                  col = "VS_Baseline_NEWS"),
  list(domain = "Vital Signs",        label = "Any Fever (≥38°C)",           col = "VS_AnyFever_38"),
  list(domain = "Vital Signs",        label = "Any Hypotension",             col = "VS_AnyHypotension"),
  # FBC
  list(domain = "FBC",                label = "WCC",                         col = "BT_Closest_WCC"),
  list(domain = "FBC",                label = "Neutrophils",                 col = "BT_Closest_Neut"),
  list(domain = "FBC",                label = "Lymphocytes",                 col = "BT_Closest_Lymph"),
  list(domain = "FBC",                label = "Haemoglobin",                 col = "BT_Closest_Hb"),
  list(domain = "FBC",                label = "Platelets",                   col = "BT_Closest_Plt"),
  list(domain = "FBC",                label = "Eosinophils",                 col = "BT_Closest_Eosin"),
  list(domain = "FBC",                label = "Monocytes",                   col = "BT_Closest_Mono"),
  list(domain = "FBC",                label = "Haematocrit",                 col = "BT_Closest_HCT"),
  # CRP
  list(domain = "CRP",                label = "CRP",                         col = "BT_Closest_CRP"),
  # U&E
  list(domain = "U&E",                label = "Sodium",                      col = "BT_Closest_Na"),
  list(domain = "U&E",                label = "Creatinine",                  col = "BT_Closest_Creat"),
  list(domain = "U&E",                label = "Urea",                        col = "BT_Closest_Urea"),
  # LFT
  list(domain = "LFT",                label = "Bilirubin",                   col = "BT_Closest_Bili"),
  list(domain = "LFT",                label = "ALP",                         col = "BT_Closest_ALP"),
  list(domain = "LFT",                label = "Albumin",                     col = "BT_Closest_Alb"),
  # Blood Gas
  list(domain = "Blood Gas",          label = "Lactate",                     col = "BG_Closest_Lactate"),
  list(domain = "Blood Gas",          label = "pH",                          col = "BG_Closest_pH"),
  list(domain = "Blood Gas",          label = "pCO2",                        col = "BG_Closest_pCO2")
)

# Build per-variable row
rows <- lapply(var_defs, function(v) {
  col <- v$col
  if (col == "__source_any__") {
    ca <- source_completeness(bc_data)
    cp <- source_completeness(pos)
    cn <- source_completeness(neg)
  } else if (!col %in% names(bc_data)) {
    return(data.table(
      Domain = v$domain, Variable = v$label,
      N_All = NA_integer_, Pct_All = NA_real_,
      N_Positive = NA_integer_, Pct_Positive = NA_real_,
      N_Negative = NA_integer_, Pct_Negative = NA_real_
    ))
  } else {
    ca <- completeness(bc_data, col)
    cp <- completeness(pos, col)
    cn <- completeness(neg, col)
  }
  data.table(
    Domain = v$domain, Variable = v$label,
    N_All = ca$n, Pct_All = ca$pct,
    N_Positive = cp$n, Pct_Positive = cp$pct,
    N_Negative = cn$n, Pct_Negative = cn$pct
  )
})

var_raw <- rbindlist(rows)

################################################################################
# Domain-level summary (all key variables complete within domain)
################################################################################

bc_data[, VS_DomainComplete      := !is.na(VS_Baseline_HR) & !is.na(VS_Baseline_Temp) &
                                    !is.na(VS_Baseline_SBP) & !is.na(VS_Baseline_NEWS)]
bc_data[, FBC_DomainComplete     := !is.na(BT_Closest_WCC) & !is.na(BT_Closest_Neut) &
                                    !is.na(BT_Closest_Lymph) & !is.na(BT_Closest_Hb) & !is.na(BT_Closest_Plt)]
bc_data[, CRP_DomainComplete     := !is.na(BT_Closest_CRP)]
bc_data[, UE_DomainComplete      := !is.na(BT_Closest_Na) & !is.na(BT_Closest_Creat) & !is.na(BT_Closest_Urea)]
bc_data[, LFT_DomainComplete     := !is.na(BT_Closest_Bili) & !is.na(BT_Closest_ALP) & !is.na(BT_Closest_Alb)]
bc_data[, BG_DomainComplete      := !is.na(BG_Closest_Lactate) & !is.na(BG_Closest_pH)]

pos <- bc_data[Outcome == 1]
neg <- bc_data[Outcome == 0]

domain_complete_defs <- list(
  list(label = "Vital Signs (all: HR, Temp, SBP, NEWS)",   col = "VS_DomainComplete"),
  list(label = "FBC (all: WCC, Neut, Lymph, Hb, Plt)",    col = "FBC_DomainComplete"),
  list(label = "CRP",                                       col = "CRP_DomainComplete"),
  list(label = "U&E (all: Na, Creat, Urea)",               col = "UE_DomainComplete"),
  list(label = "LFT (all: Bili, ALP, Alb)",                col = "LFT_DomainComplete"),
  list(label = "Blood Gas (all: Lactate, pH)",             col = "BG_DomainComplete")
)

domain_rows <- lapply(domain_complete_defs, function(d) {
  ca <- list(n = sum(bc_data[[d$col]], na.rm = TRUE),
             pct = round(100 * mean(bc_data[[d$col]], na.rm = TRUE), 1))
  cp <- list(n = sum(pos[[d$col]], na.rm = TRUE),
             pct = round(100 * mean(pos[[d$col]], na.rm = TRUE), 1))
  cn <- list(n = sum(neg[[d$col]], na.rm = TRUE),
             pct = round(100 * mean(neg[[d$col]], na.rm = TRUE), 1))
  data.table(
    Domain = d$label,
    N_All = ca$n, Pct_All = ca$pct,
    N_Positive = cp$n, Pct_Positive = cp$pct,
    N_Negative = cn$n, Pct_Negative = cn$pct
  )
})

domain_raw <- rbindlist(domain_rows)

################################################################################
# Cohort header row
################################################################################

cohort_info <- data.table(
  Metric = c("Total cultures (all)", "Culture positives", "Culture negatives"),
  N      = c(n_all, n_pos, n_neg)
)

################################################################################
# Format helpers
################################################################################

fmt_np <- function(n, pct) sprintf("%s (%s%%)", format(n, big.mark = ","), pct)

var_fmt <- copy(var_raw)
var_fmt[, All       := fmt_np(N_All,      Pct_All)]
var_fmt[, Positives := fmt_np(N_Positive, Pct_Positive)]
var_fmt[, Negatives := fmt_np(N_Negative, Pct_Negative)]
var_fmt <- var_fmt[, .(Domain, Variable, All, Positives, Negatives)]

domain_fmt <- copy(domain_raw)
domain_fmt[, All       := fmt_np(N_All,      Pct_All)]
domain_fmt[, Positives := fmt_np(N_Positive, Pct_Positive)]
domain_fmt[, Negatives := fmt_np(N_Negative, Pct_Negative)]
domain_fmt <- domain_fmt[, .(Domain, All, Positives, Negatives)]

################################################################################
# Save
################################################################################

out_dir <- "Output Files/Publication_Tables_Figures/Per_Model_Summaries"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_xlsx(
  list(
    Cohort_Summary        = as.data.frame(cohort_info),
    Variable_Completeness = as.data.frame(var_fmt),
    Domain_Completeness   = as.data.frame(domain_fmt),
    Variable_Raw          = as.data.frame(var_raw),
    Domain_Raw            = as.data.frame(domain_raw)
  ),
  file.path(out_dir, "Model_Variable_Completeness_ByOutcome.xlsx")
)
