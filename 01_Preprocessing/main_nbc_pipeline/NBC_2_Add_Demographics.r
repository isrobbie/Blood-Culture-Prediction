library(data.table)
library(lubridate)

# ============================================================
# ADD DEMOGRAPHICS TO NON-BC ADMISSIONS
# ============================================================
# Equivalent to BC_2. Adds age, sex, mortality, LOS, IMD,
# and ethnicity. IMD and ethnicity use fallbacks from any
# historical inpatient or ED visit.
# ============================================================

nbc_df <- fread("Output Files/CSVs/NonBC_Admissions/nbc_cohort.csv",
                colClasses = c(ClusterID = "character"))
nbc_df[, AdmissionDate := as.POSIXct(AdmissionDate, tz = "UTC")]
nbc_df[, DischargeDate := as.POSIXct(DischargeDate, tz = "UTC")]

cohort_ids <- unique(nbc_df$ClusterID)

# ============================================================
# 1. JOIN DEMOGRAPHICS
# ============================================================

demographics <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_CohortDemographics.csv",
                      colClasses = c(ClusterID = "character"))
demographics[, `:=`(
  LinkedBirthMonth = ymd_hms(fifelse(LinkedBirthMonth == "NULL", NA_character_, LinkedBirthMonth)),
  LinkedDeathdate  = ymd_hms(fifelse(LinkedDeathdate  == "NULL", NA_character_, LinkedDeathdate))
)]

nbc_df <- merge(nbc_df, demographics[, .(ClusterID, LinkedBirthMonth, LinkedSex, LinkedDeathdate)],
                by = "ClusterID", all.x = TRUE)

# ============================================================
# 2. CALCULATE AGE AND MORTALITY METRICS
# ============================================================

nbc_df[, `:=`(
  DEMOGRAPHICS_LinkedSex = LinkedSex,
  AgeAtCulture           = floor(as.numeric(difftime(AdmissionDate, LinkedBirthMonth, units = "days")) / 365.25),
  DaysToDeath            = as.numeric(difftime(LinkedDeathdate, AdmissionDate, units = "days")),
  LengthOfStay_days      = as.numeric(difftime(DischargeDate, AdmissionDate, units = "days"))
)]
nbc_df[, `:=`(
  Mort_7   = as.integer(!is.na(DaysToDeath) & DaysToDeath >= 0 & DaysToDeath <= 7),
  Mort_28  = as.integer(!is.na(DaysToDeath) & DaysToDeath >= 0 & DaysToDeath <= 28),
  Mort_365 = as.integer(!is.na(DaysToDeath) & DaysToDeath >= 0 & DaysToDeath <= 365)
)]
nbc_df[, c("LinkedBirthMonth", "LinkedSex", "LinkedDeathdate") := NULL]

# ============================================================
# 3. LOAD ED AND INPATIENT DATA FOR SOCIAL DEMOGRAPHICS
# ============================================================

inpat_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_InpatEpisodes.csv",
                  select = c("ClusterID", "IMDScore", "IMDPercentile", "EthnicGroupCode"),
                  colClasses = c(ClusterID = "character"))
ed_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_EDAttendances.csv",
               select = c("ClusterID", "IMDScore", "IMDPercentile", "EthnicGroupCode"),
               colClasses = c(ClusterID = "character"))
inpat_df <- inpat_df[ClusterID %in% cohort_ids]
ed_df    <- ed_df[ClusterID %in% cohort_ids]

# ============================================================
# 4. FILL MISSING IMD FROM ANY HISTORICAL VISIT
# ============================================================

ip_imd <- inpat_df[!is.na(IMDPercentile),
                    .(FB_IMDScore = first(as.numeric(IMDScore)),
                      FB_IMD      = first(as.integer(IMDPercentile))), by = ClusterID]
ed_imd <- ed_df[!is.na(IMDPercentile),
                 .(FB_IMDScore = first(as.numeric(IMDScore)),
                   FB_IMD      = first(as.integer(IMDPercentile))), by = ClusterID]
all_imd <- rbindlist(list(ip_imd, ed_imd))[!duplicated(ClusterID)]

nbc_df[, SOCIAL_IMDPercentile := suppressWarnings(as.integer(fifelse(
  IMDPercentile == "NULL" | is.na(IMDPercentile), NA_character_, as.character(IMDPercentile))))]
nbc_df[, SOCIAL_IMDScore := suppressWarnings(as.numeric(fifelse(
  IMDScore == "NULL" | is.na(IMDScore), NA_character_, as.character(IMDScore))))]

nbc_df <- merge(nbc_df, all_imd, by = "ClusterID", all.x = TRUE)
nbc_df[is.na(SOCIAL_IMDPercentile) & !is.na(FB_IMD),      SOCIAL_IMDPercentile := FB_IMD]
nbc_df[is.na(SOCIAL_IMDScore)      & !is.na(FB_IMDScore),  SOCIAL_IMDScore      := FB_IMDScore]
nbc_df[, c("FB_IMD", "FB_IMDScore") := NULL]

# ============================================================
# 5. FILL MISSING ETHNICITY FROM ANY HISTORICAL VISIT
# ============================================================

ip_eth <- inpat_df[!is.na(EthnicGroupCode) &
                    !(trimws(EthnicGroupCode) %in% c("Z", "99", "NULL", "")),
                    .(FB_Ethnic = first(trimws(EthnicGroupCode))), by = ClusterID]
ed_eth <- ed_df[!is.na(EthnicGroupCode) &
                 !(trimws(EthnicGroupCode) %in% c("Z", "99", "NULL", "")),
                 .(FB_Ethnic = first(trimws(EthnicGroupCode))), by = ClusterID]
all_eth <- rbindlist(list(ip_eth, ed_eth))[!duplicated(ClusterID)]

nbc_df[, SOCIAL_EthnicGroupCode := fifelse(
  is.na(EthnicGroupCode) | trimws(EthnicGroupCode) %in% c("Z", "99", "NULL", ""),
  NA_character_,
  trimws(EthnicGroupCode))]

nbc_df <- merge(nbc_df, all_eth, by = "ClusterID", all.x = TRUE)
nbc_df[(is.na(SOCIAL_EthnicGroupCode) |
        SOCIAL_EthnicGroupCode %in% c("Z", "99", "NULL")) &
       !is.na(FB_Ethnic),
       SOCIAL_EthnicGroupCode := FB_Ethnic]
nbc_df[, FB_Ethnic := NULL]

nbc_df[, SOCIAL_Ethnicity := fcase(
  SOCIAL_EthnicGroupCode %in% c("A", "B", "C"), "White",
  SOCIAL_EthnicGroupCode %in% c("D", "E", "F", "G"), "Mixed",
  SOCIAL_EthnicGroupCode %in% c("H", "J", "K", "L"), "Asian",
  SOCIAL_EthnicGroupCode %in% c("M", "N", "P"), "Black",
  SOCIAL_EthnicGroupCode %in% c("R", "S"), "Other",
  default = "Unknown"
)]

# ============================================================
# 6. SAVE OUTPUT
# ============================================================

fwrite(nbc_df, "Output Files/CSVs/NonBC_Admissions/nbc_with_demographics.csv")
