library(data.table)
library(lubridate)

# ============================================================
# ADD PRIOR BLOOD CULTURE HISTORY TO NON-BC ADMISSIONS
# ============================================================
# Equivalent to BC_1 sections 6C/6D. Counts gram-positive and
# gram-negative positive blood cultures in the 28 days prior to
# admission (excluding the 24h immediately before admission).
# ============================================================

nbc_df <- fread("Output Files/CSVs/NonBC_Admissions/nbc_with_demographics.csv",
                colClasses = c(ClusterID = "character"))
nbc_df[, AdmissionDate := as.POSIXct(AdmissionDate, tz = "UTC")]

cohort_ids <- unique(nbc_df$ClusterID)

# ============================================================
# 1. LOAD AND CLASSIFY BLOOD CULTURES
# ============================================================

bug_lookup <- fread("00_Data/lookup_tables/bug_group_table.csv")
setkey(bug_lookup, BugCode)

micro_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_Microbiology.csv",
                  select = c("ClusterID", "AccessionNumber", "CollectionDateTime", "BatTestName", "BugCode"))
bc_micro <- micro_df[BatTestName == "BLOOD CULTURE"]
bc_micro[, ClusterID := as.character(ClusterID)]
bc_micro[, CollectionDateTime := ymd_hms(CollectionDateTime)]
bc_micro <- bc_micro[!is.na(CollectionDateTime)]
bc_micro[, BugGroup := bug_lookup[bc_micro, BugGroup, on = "BugCode"]]

bc_classify <- bc_micro[, .(
  is_gram_pos = any(BugGroup %in% c("S AUREUS", "ENTEROCOCCUS", "S PNEUMONIAE",
                                     "BETA HAEMOLYTIC STREP", "OTHER PATHOGENIC STREPS"), na.rm = TRUE),
  is_gram_neg = any(BugGroup %in% c("E COLI", "KLEBSIELLA", "P AERUGINOSA",
                                     "ENTEROBACTER", "OTHER ENTEROBACTERALES"), na.rm = TRUE)
), by = .(ClusterID, AccessionNumber, CollectionDateTime)]

# One row per unique culture (AccessionNumber) for total count
bc_unique_times <- unique(bc_micro[, .(ClusterID, AccessionNumber, CollectionDateTime)])

# ============================================================
# 2. MATCH TO ADMISSIONS: PRIOR 28 DAYS (EXCLUDING <1 DAY)
# ============================================================

nbc_df[, `:=`(adm_minus_28d = AdmissionDate - days(28),
              adm_minus_1d  = AdmissionDate - days(1))]

# Total prior BCs (all cultures, not just positive)
all_bc <- bc_unique_times[nbc_df,
  on = .(ClusterID, CollectionDateTime >= adm_minus_28d, CollectionDateTime < adm_minus_1d),
  .(AdmissionID = i.AdmissionID),
  allow.cartesian = TRUE, nomatch = 0L]
total_bc_agg <- all_bc[, .(CulturesInPrior28Days = .N), by = AdmissionID]

# Gram-classified positive cultures
prior_bc <- bc_classify[nbc_df,
  on = .(ClusterID, CollectionDateTime >= adm_minus_28d, CollectionDateTime < adm_minus_1d),
  .(AdmissionID = i.AdmissionID, is_gram_pos = x.is_gram_pos, is_gram_neg = x.is_gram_neg),
  allow.cartesian = TRUE, nomatch = 0L]

prior_bc_agg <- prior_bc[, .(
  PositiveCulturesInPrior28Days     = sum(is_gram_pos | is_gram_neg, na.rm = TRUE),
  GramPositiveCulturesInPrior28Days = sum(is_gram_pos, na.rm = TRUE),
  GramNegativeCulturesInPrior28Days = sum(is_gram_neg, na.rm = TRUE)
), by = AdmissionID]

nbc_df <- merge(nbc_df, total_bc_agg, by = "AdmissionID", all.x = TRUE)
nbc_df <- merge(nbc_df, prior_bc_agg, by = "AdmissionID", all.x = TRUE)
for (col in c("CulturesInPrior28Days", "PositiveCulturesInPrior28Days",
              "GramPositiveCulturesInPrior28Days", "GramNegativeCulturesInPrior28Days")) {
  nbc_df[is.na(get(col)), (col) := 0L]
}
nbc_df[, c("adm_minus_28d", "adm_minus_1d") := NULL]

# ============================================================
# 3. SAVE OUTPUT
# ============================================================

fwrite(nbc_df, "Output Files/CSVs/NonBC_Admissions/nbc_with_prior_bc.csv")
