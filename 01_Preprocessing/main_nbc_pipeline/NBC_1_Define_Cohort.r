library(data.table)
library(lubridate)

# ============================================================
# DEFINE NON-BC EMERGENCY ADMISSION COHORT
# ============================================================
# Equivalent to BC_1 initial extract + BC_3 admission matching.
# Identifies emergency admissions in the study period that did
# NOT have a blood culture within +-24h of admission.
# Anchor: AdmissionDate (vs CollectionDateTime in BC pipeline).
# ============================================================

study_start <- as.POSIXct("2024-01-01", tz = "UTC")
study_end   <- as.POSIXct("2024-12-31 23:59:59", tz = "UTC")

dir.create("Output Files/CSVs/NonBC_Admissions", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. LOAD BLOOD CULTURE TIMES (FOR EXCLUSION)
# ============================================================

micro_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_Microbiology.csv",
                  select = c("ClusterID", "CollectionDateTime", "BatTestName"))
bc_micro <- micro_df[BatTestName == "BLOOD CULTURE"]
bc_micro[, CollectionDateTime := ymd_hms(CollectionDateTime)]
bc_micro <- bc_micro[!is.na(CollectionDateTime)]
bc_times <- unique(bc_micro[, .(ClusterID = as.character(ClusterID), BC_DateTime = CollectionDateTime)])

# ============================================================
# 2. LOAD AND FILTER EMERGENCY ADMISSIONS
# ============================================================

inpat_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_InpatEpisodes.csv",
                  colClasses = c(ClusterID = "character"))
inpat_df[, `:=`(
  AdmissionDate    = suppressWarnings(ymd_hms(AdmissionDate)),
  DischargeDate    = suppressWarnings(ymd_hms(DischargeDate)),
  EpisodeStartDate = suppressWarnings(ymd_hms(EpisodeStartDate))
)]

emergency_codes <- c("21", "22", "23", "24", "25", "2A", "2C", "2D", "28")

# Admitting episode specialty per spell (first episode by EpisodeStartDate)
admitting_ep <- inpat_df[!is.na(EpisodeStartDate)][
  order(ClusterID, SpellID, EpisodeStartDate),
  .(ConsultantMainSpecialtyCode = first(ConsultantMainSpecialtyCode),
    TreatmentFunctionCode       = first(TreatmentFunctionCode),
    LSOA2011Code                = first(LSOA2011Code),
    LSOA2001Code                = first(LSOA2001Code)),
  by = .(ClusterID, SpellID)]

spells <- unique(inpat_df[, .(ClusterID, SpellID, AdmissionDate, DischargeDate, AdmissionMethodCode,
                               IMDScore, IMDPercentile, EthnicGroupCode)],
                  by = c("ClusterID", "SpellID"))
spells <- merge(spells, admitting_ep, by = c("ClusterID", "SpellID"), all.x = TRUE)

spells <- spells[AdmissionMethodCode %in% emergency_codes &
                  !is.na(AdmissionDate) &
                  AdmissionDate >= study_start &
                  AdmissionDate <= study_end]

# ============================================================
# 3. EXCLUDE ADMISSIONS WITH A BLOOD CULTURE WITHIN +-48H
# ============================================================

spells[, `:=`(adm_minus_48h = AdmissionDate - hours(48),
              adm_plus_48h  = AdmissionDate + hours(48))]

has_bc <- bc_times[spells,
                   on = .(ClusterID, BC_DateTime >= adm_minus_48h, BC_DateTime <= adm_plus_48h),
                   nomatch = 0L,
                   .(ClusterID, SpellID, BC_DateTime)][, .(has_bc = 1L), by = .(ClusterID, SpellID)]

spells <- has_bc[spells, on = .(ClusterID, SpellID)]
nbc_df <- spells[is.na(has_bc)]
nbc_df[, c("has_bc", "adm_minus_48h", "adm_plus_48h") := NULL]

nbc_df[, AdmissionID := paste0(ClusterID, "_", SpellID, "_", format(AdmissionDate, "%Y%m%d_%H%M%S"))]

# ============================================================
# 4. SAVE OUTPUT
# ============================================================

fwrite(nbc_df, "Output Files/CSVs/NonBC_Admissions/nbc_cohort.csv")
