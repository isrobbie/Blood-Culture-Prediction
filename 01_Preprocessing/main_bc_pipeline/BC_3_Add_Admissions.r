library(data.table)
library(dplyr)
library(lubridate)
library(writexl)
library(stringr)

# Load Excel output formatting helpers
source("01_Preprocessing/Excel_Output_Helper.r")

# ============================================================
# ADD ADMISSION AND ICU INFORMATION TO BLOOD CULTURES
# ============================================================
# This script adds ED attendance, inpatient admission, and ICU
# information to primary blood cultures 
# ============================================================

# Read blood cultures with demographics (from previous step)
bc_df <- fread("Output Files/CSVs/blood_cultures_with_demographics.csv")
initial_bc_count <- nrow(bc_df)

# Extract cohort ClusterIDs for filtering
cohort_cluster_ids <- unique(bc_df$ClusterID)

# Read other data (using fread for speed), filtered to cohort
ed_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_EDAttendances.csv")
ed_df <- ed_df[ClusterID %in% cohort_cluster_ids]

inpat_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_InpatEpisodes.csv")
inpat_df <- inpat_df[ClusterID %in% cohort_cluster_ids]

ward_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_InpatWardStays.csv")
ward_df <- ward_df[ClusterID %in% cohort_cluster_ids]

outpat_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_OutpatAppointments.csv")
outpat_df <- outpat_df[ClusterID %in% cohort_cluster_ids]

# Fix duplicate column names in source data (ward_df has ClusterID twice at positions 1 and 7)
if(anyDuplicated(names(ward_df)) > 0) {
  # Keep only non-duplicate columns
  ward_df <- ward_df[, !duplicated(names(ward_df)), with = FALSE]
}
if(anyDuplicated(names(ed_df)) > 0) {
  ed_df <- ed_df[, !duplicated(names(ed_df)), with = FALSE]
}
if(anyDuplicated(names(inpat_df)) > 0) {
  inpat_df <- inpat_df[, !duplicated(names(inpat_df)), with = FALSE]
}
if(anyDuplicated(names(outpat_df)) > 0) {
  outpat_df <- outpat_df[, !duplicated(names(outpat_df)), with = FALSE]
}

# ============================================================
# 1. PREPARE DATA
# ============================================================

# Define ICU ward pattern (used multiple times throughout script)
ICU_PATTERN <- "ICU|CTVCC|CTCC|Crit Care|CritCare|HDU|CCU|CIU|SCBU"

# Define wards to EXCLUDE from ICU categorization (even if they match pattern)
EXCLUDED_WARDS <- c("NOC HDU", "Neuro HDU", "H-WD CCU")

# Helper function to identify ICU wards
is_icu_ward <- function(ward_name) {
  case_when(
    is.na(ward_name) ~ 0L,
    ward_name %in% EXCLUDED_WARDS ~ 0L,  # Explicitly exclude these wards
    str_detect(ward_name, ICU_PATTERN) ~ 1L,
    TRUE ~ 0L
  )
}

# Helper function to map admission method codes to readable names
map_admission_method <- function(admission_code) {
  case_when(
    is.na(admission_code) ~ "Not admitted/Unknown",
    admission_code %in% c("11", "12", "13") ~ "Elective",
    admission_code %in% c("21", "22", "23", "24", "25", "2A", "2B", "2D", "28") ~ "Emergency",
    admission_code %in% c("31", "32", "82", "83", "2C") ~ "Maternity",
    TRUE ~ "Other"
  )
}

# Convert ClusterID to character for all datasets
ed_df <- ed_df %>% mutate(ClusterID = as.character(ClusterID))
inpat_df <- inpat_df %>% mutate(ClusterID = as.character(ClusterID))
bc_df <- bc_df %>% mutate(ClusterID = as.character(ClusterID))
ward_df <- ward_df %>% mutate(ClusterID = as.character(ClusterID))
outpat_df <- outpat_df %>% mutate(ClusterID = as.character(ClusterID))

# Convert CollectionDateTime - ensure it's POSIXct
bc_df <- bc_df %>% 
  mutate(CollectionDateTime = if(is.character(CollectionDateTime)) {
    ymd_hms(CollectionDateTime)
  } else {
    as.POSIXct(CollectionDateTime)
  })

# Convert dates in ed_df - filter out invalid dates first  
ed_df <- ed_df %>%
  filter(!is.na(ArrivalDateTime) & !is.na(DepartureDateTime)) %>%
  mutate(across(c(ArrivalDateTime, DepartureDateTime),
                ~ if(is.character(.)) {
                  ymd_hms(if_else(. == "NULL", NA_character_, .))
                } else {
                  as.POSIXct(.)
                }))

# Convert dates in inpat_df
inpat_df <- inpat_df %>%
  mutate(across(c(AdmissionDate, DischargeDate),
                ~ if(is.character(.)) {
                  ymd_hms(if_else(. == "NULL", NA_character_, .))
                } else {
                  as.POSIXct(.)
                }))

# Convert dates in ward_df
ward_df <- ward_df %>%
  mutate(across(c(WardStartDate, WardEndDate),
                ~ if(is.character(.)) {
                  ymd_hms(if_else(. == "NULL", NA_character_, .))
                } else {
                  as.POSIXct(.)
                }))

# Convert dates in outpat_df
outpat_df <- outpat_df %>%
  filter(!is.na(AppointmentDate)) %>%
  mutate(AppointmentDate = if(is.character(AppointmentDate)) {
    ymd_hms(if_else(AppointmentDate == "NULL", NA_character_, AppointmentDate))
  } else {
    as.POSIXct(AppointmentDate)
  })

# Ensure ward_df stays as data.table for non-equi joins
setDT(ward_df)

# ============================================================
# 2. LOAD SPECIALTY MAPPING
# ============================================================

# Load specialty mapping from CSV
specialty_lookup <- fread("00_Data/lookup_tables/specialty_mapping.csv")
setkey(specialty_lookup, SpecialtyCode)

# Function to map specialty code to name
get_specialty_name <- function(spec_codes) {
  specialty_lookup[J(spec_codes), SpecialtyName, on = "SpecialtyCode"]
}

# ============================================================
# 3. JOIN ED ATTENDANCES (IF CULTURE DURING ED VISIT)
# ============================================================

bc_ed_during <- bc_df %>%
  inner_join(
    ed_df %>% select(ClusterID, ArrivalDateTime, DepartureDateTime, AttendanceID),
    by = "ClusterID", 
    relationship = "many-to-many"
  ) %>%
  filter(CollectionDateTime >= ArrivalDateTime & CollectionDateTime <= DepartureDateTime) %>%
  # If multiple ED attendances match, keep the most recent one per specimen
  group_by(BloodCultureID) %>%
  arrange(desc(ArrivalDateTime)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(CultureLocation = "ED",
         TimeInED_minutes = as.numeric(difftime(DepartureDateTime, ArrivalDateTime, units = "mins"))) %>%
  select(BloodCultureID, ClusterID, AccessionNumber, CollectionDateTime, SpecimenFull, AttendanceID, ArrivalDateTime, DepartureDateTime, 
         CultureLocation, TimeInED_minutes)

# ============================================================
# 3B. MATCH OUTPATIENT APPOINTMENTS (FOR UNKNOWN SOURCES)
# ============================================================

# Match blood cultures to outpatient appointments on same DATE
# This helps identify cultures taken during outpatient visits
bc_outpat_match <- bc_df %>%
  mutate(CollectionDate = as.Date(CollectionDateTime)) %>%
  left_join(
    outpat_df %>% 
      mutate(AppointmentDateOnly = as.Date(AppointmentDate)) %>%
      select(ClusterID, AppointmentDateOnly, 
             OutpatAppointment_SpecialtyCode = ConsultantMainSpecialtyCode,
             OutpatAppointment_EpisodeID = EpisodeID),
    by = c("ClusterID", "CollectionDate" = "AppointmentDateOnly"),
    relationship = "many-to-many",
    suffix = c("", "_outpat")
  ) %>%
  # If multiple appointments on same day, keep first one (or keep NA if no match)
  group_by(BloodCultureID) %>%
  arrange(BloodCultureID) %>%
  slice(1) %>%
  ungroup() %>%
  select(BloodCultureID, OutpatAppointment_EpisodeID, OutpatAppointment_SpecialtyCode)

# ============================================================
# 4. ADD ICU AND INPATIENT ADMISSION INFORMATION
# ============================================================

# ICU identification via ward names - create granular ICU timing

# ICU at time of culture: patient must be IN ward when culture taken
# NOTE: Using data.table for complex non-equi join (ward time overlap)
bc_ward_match_at <- as.data.table(bc_df[, c("BloodCultureID", "ClusterID", "AccessionNumber", "CollectionDateTime", "SpecimenFull")])
bc_ward_match_at[, bc_id := .I]

# Join to get wards active at time of culture (data.table non-equi join)
bc_wards_at_culture <- ward_df[bc_ward_match_at,
                                on = .(ClusterID,
                                       WardStartDate <= CollectionDateTime,
                                       WardEndDate > CollectionDateTime),
                                nomatch = 0L,
                                allow.cartesian = TRUE]

# Clean up duplicate columns from non-equi join and convert to tibble
# data.table joins create i.column for right table columns
bc_wards_at_culture <- bc_wards_at_culture %>%
  as_tibble() %>%
  select(-any_of(c("i.WardStartDate", "i.WardEndDate", "i.CollectionDateTime"))) %>%
  mutate(ICU_Ward = is_icu_ward(WardName))

# Create ward names at culture (concatenate all wards, mark ICU with ***)
wards_at_culture_names <- bc_wards_at_culture %>%
  mutate(WardDisplay = if_else(ICU_Ward == 1, paste0("***", WardName, "***"), WardName)) %>%
  group_by(bc_id) %>%
  summarise(Ward_at_Culture = paste(WardDisplay, collapse = "; "), .groups = "drop") %>%
  left_join(bc_ward_match_at %>% select(bc_id, BloodCultureID), 
            by = "bc_id") %>%
  select(-bc_id)

# Aggregate by the bc_id which links back to BloodCultureID
icu_at_culture <- bc_wards_at_culture %>%
  group_by(bc_id) %>%
  summarise(ICU_at_Culture = max(ICU_Ward, na.rm = TRUE), .groups = "drop") %>%
  left_join(bc_ward_match_at %>% select(bc_id, BloodCultureID), 
            by = "bc_id") %>%
  select(-bc_id)

# ICU after culture (0-24h post-culture): any ward that starts within 24h after culture
# NOTE: Using data.table for complex non-equi join (time window)
bc_ward_match_after <- as.data.table(bc_df[, c("BloodCultureID", "ClusterID", "AccessionNumber", "CollectionDateTime", "SpecimenFull")])
bc_ward_match_after[, `:=`(bc_id = .I,
                            culture_plus_24h = CollectionDateTime + hours(24))]

bc_wards_after_culture <- ward_df[bc_ward_match_after,
                                   on = .(ClusterID,
                                          WardStartDate >= CollectionDateTime,
                                          WardStartDate <= culture_plus_24h),
                                   nomatch = 0L,
                                   allow.cartesian = TRUE]

# Clean up duplicate columns from non-equi join and convert to tibble
bc_wards_after_culture <- bc_wards_after_culture %>%
  as_tibble() %>%
  select(-any_of(c("i.WardStartDate", "i.CollectionDateTime", "culture_plus_24h"))) %>%
  mutate(ICU_Ward = is_icu_ward(WardName))

icu_after_culture <- bc_wards_after_culture %>%
  group_by(bc_id) %>%
  summarise(ICU_after_Culture = max(ICU_Ward, na.rm = TRUE), .groups = "drop") %>%
  left_join(bc_ward_match_after %>% 
              select(bc_id, BloodCultureID), 
            by = "bc_id") %>%
  select(-bc_id)

# Match admissions using data.table non-equi joins
# NOTE: Using data.table for complex time window join (24h before admission through discharge)
bc_for_join <- as.data.table(bc_df[, c("BloodCultureID", "ClusterID", "AccessionNumber", "CollectionDateTime", "SpecimenFull")])
bc_for_join[, join_time := CollectionDateTime]

# Set up time windows: 24h before admission through discharge (for spell-level matching)
inpat_for_join <- as.data.table(inpat_df[, c("ClusterID", "SpellID", "AdmissionDate", "DischargeDate", 
                                               "AdmissionMethodCode")])
inpat_for_join[, `:=`(join_time_start = AdmissionDate - hours(24),
                       join_time_end = DischargeDate)]

# Non-equi join on ClusterID and time window (data.table)
bc_matched_admissions <- bc_for_join[inpat_for_join,
                                      on = .(ClusterID,
                                             join_time >= join_time_start,
                                             join_time <= join_time_end),
                                      nomatch = NA,
                                      allow.cartesian = TRUE]

# Clean up duplicate columns from non-equi join
bc_matched_admissions <- bc_matched_admissions %>%
  as_tibble() %>%
  select(-any_of(c("join_time", "i.join_time", "join_time_start", "join_time_end")))

# ============================================================
# SPECIALTY AT TIME OF CULTURE (avoid data leakage)
# ============================================================
# Match specialty from the EPISODE active at time of culture, not a later episode
# This prevents leakage where specialty could be influenced by culture result

episode_for_specialty <- as.data.table(inpat_df[, c("ClusterID", "EpisodeID", "EpisodeStartDate", 
                                                      "EpisodeEndDate", "ConsultantMainSpecialtyCode")])
# Convert dates (same approach as inpat_df earlier)
episode_for_specialty <- episode_for_specialty %>%
  mutate(across(c(EpisodeStartDate, EpisodeEndDate),
                ~ if(is.character(.)) {
                  ymd_hms(if_else(. == "NULL", NA_character_, .))
                } else {
                  as.POSIXct(.)
                }))
setDT(episode_for_specialty)

bc_for_specialty <- as.data.table(bc_df[, c("BloodCultureID", "ClusterID", "CollectionDateTime")])

# Join to get episode active at time of culture
bc_specialty_match <- episode_for_specialty[bc_for_specialty,
                                             on = .(ClusterID,
                                                    EpisodeStartDate <= CollectionDateTime,
                                                    EpisodeEndDate >= CollectionDateTime),
                                             nomatch = NA,
                                             allow.cartesian = TRUE]

# Keep one specialty per blood culture (if multiple, use earliest episode start)
bc_specialty_match <- bc_specialty_match %>%
  as_tibble() %>%
  filter(!is.na(ConsultantMainSpecialtyCode)) %>%
  group_by(BloodCultureID) %>%
  arrange(EpisodeStartDate) %>%
  slice(1) %>%
  ungroup() %>%
  select(BloodCultureID, ConsultantMainSpecialtyCode_AtCulture = ConsultantMainSpecialtyCode)

# Add specialty at time of culture to matched admissions
bc_matched_admissions <- bc_matched_admissions %>%
  left_join(bc_specialty_match, by = "BloodCultureID")

# ICU any time during entire admission - match wards from admission to discharge (data.table)
bc_admission_wards <- ward_df[as.data.table(bc_matched_admissions),
                               on = .(ClusterID,
                                      WardStartDate >= AdmissionDate,
                                      WardStartDate <= DischargeDate),
                               nomatch = 0L,
                               allow.cartesian = TRUE]

# Clean up duplicate columns from non-equi join
bc_admission_wards <- bc_admission_wards %>%
  as_tibble() %>%
  select(-any_of(c("i.WardStartDate", "i.AdmissionDate", "i.DischargeDate"))) %>%
  mutate(ICU_Ward = is_icu_ward(WardName))

# Create list of all wards during admission (mark ICU with ***)
all_wards_during_admission <- bc_admission_wards %>%
  mutate(WardDisplay = if_else(ICU_Ward == 1, paste0("***", WardName, "***"), WardName)) %>%
  group_by(BloodCultureID) %>%
  summarise(
    All_Wards_During_Admission = paste(unique(WardDisplay), collapse = "; "),
    .groups = "drop"
  )

icu_any <- bc_admission_wards %>%
  group_by(BloodCultureID) %>%
  summarise(
    Has_ICU_Ward = max(ICU_Ward, na.rm = TRUE), 
    .groups = "drop"
  )

# Add ICU information and ward names from ward lookups
bc_matched_admissions <- bc_matched_admissions %>%
  left_join(icu_any, by = "BloodCultureID") %>%
  left_join(icu_at_culture, by = "BloodCultureID") %>%
  left_join(icu_after_culture, by = "BloodCultureID") %>%
  left_join(all_wards_during_admission, by = "BloodCultureID")

# For multiple admission matches, choose the admission closest to the culture time
bc_matched_admissions <- bc_matched_admissions %>%
  mutate(TimeDiff = abs(as.numeric(difftime(CollectionDateTime, AdmissionDate, units = "hours")))) %>%
  arrange(BloodCultureID, TimeDiff) %>%
  group_by(BloodCultureID) %>%
  slice(1) %>%
  ungroup() %>%
  select(-TimeDiff) %>%
  select(BloodCultureID, 
         SpellID, AdmissionDate, DischargeDate, 
         AdmissionMethodCode, ConsultantMainSpecialtyCode_AtCulture,
         Has_ICU_Ward, ICU_at_Culture, ICU_after_Culture,
         All_Wards_During_Admission)

# Join admission info back to FULL bc_df (preserves all blood cultures)
bc_df <- bc_df %>%
  left_join(bc_matched_admissions, by = "BloodCultureID")

# Add ED location information (preserve where culture was taken)
ed_locations <- bc_ed_during %>%
  select(BloodCultureID, 
         CultureLocation, AttendanceID, 
         ArrivalDateTime, DepartureDateTime, TimeInED_minutes)

bc_df <- bc_df %>%
  left_join(ed_locations, by = "BloodCultureID") %>%
  left_join(wards_at_culture_names, by = "BloodCultureID") %>%
  left_join(bc_outpat_match, by = "BloodCultureID")

# Culture Location: Where was the blood culture physically taken?
# And create all derived variables
bc_df <- bc_df %>%
  mutate(
    # Culture Location
    CultureLocation = case_when(
      !is.na(CultureLocation) ~ CultureLocation,  # Already marked as ED
      !is.na(SpellID) & !is.na(AdmissionDate) & CollectionDateTime >= AdmissionDate ~ "Inpatient Ward",
      !is.na(OutpatAppointment_SpecialtyCode) ~ "Outpatient",  # Matched to outpatient appointment
      TRUE ~ "Unknown"  # Could not match to any source
    ),
    
    # Time calculations
    HoursFromAdmissionToCulture = as.numeric(difftime(CollectionDateTime, AdmissionDate, units = "hours")),
    HoursFromArrivalToCulture = as.numeric(difftime(CollectionDateTime, ArrivalDateTime, units = "hours")),
    LengthOfStay_days = as.numeric(difftime(DischargeDate, AdmissionDate, units = "days")),
    
    # Length of stay BEFORE culture (days from admission to culture)
    LOS_BeforeCulture_days = as.numeric(difftime(CollectionDateTime, AdmissionDate, units = "days")),
    
    # ICU flags
    ICU_at_Culture = if_else(is.na(ICU_at_Culture) | ICU_at_Culture == 0, 0L, 1L),
    ICU_after_Culture = if_else(is.na(ICU_after_Culture) | ICU_after_Culture == 0, 0L, 1L),
    ICU_any_time = if_else(is.na(Has_ICU_Ward) | Has_ICU_Ward == 0, 0L, 1L),
    
    # Admission Timing: When was culture taken relative to admission? (Detailed breakdown)
    AdmissionTiming = case_when(
      is.na(AdmissionDate) ~ "Not admitted",
      HoursFromAdmissionToCulture < 0 ~ "Admitted within 24h after culture",  # Culture before admission
      HoursFromAdmissionToCulture >= 0 & HoursFromAdmissionToCulture <= 48 ~ "Already admitted 0-48h",
      HoursFromAdmissionToCulture > 48 ~ "Already admitted >48h",
      TRUE ~ "Not admitted"
    ),
    
    # CultureBeforeAdmission: Binary flag for cultures taken before formal admission
    # (e.g., ED patient about to be admitted) - more parsimonious than AdmissionTiming dummies
    # TRUE = culture taken before admission (LOS_BeforeCulture_days < 0)
    CultureBeforeAdmission = fifelse(
      !is.na(AdmissionDate) & HoursFromAdmissionToCulture < 0, 1L, 0L
    ),
    
    # Community-onset vs Hospital-onset (Qingze definition for comparison)
    # Community-onset: Culture taken between 24h BEFORE admission and 48h AFTER admission
    # Hospital-onset: Culture taken more than 48h after admission
    CommunityOnset = case_when(
      is.na(AdmissionDate) ~ NA_character_,  # Exclude episodes without admission data
      HoursFromAdmissionToCulture >= -24 & HoursFromAdmissionToCulture <= 48 ~ "Community-onset",
      HoursFromAdmissionToCulture > 48 ~ "Hospital-onset",
      TRUE ~ NA_character_  # Cultures more than 24h before admission (excluded)
    ),
    
    # Admission Method: How were they admitted?
    AdmissionMethod = map_admission_method(AdmissionMethodCode),
    
    # Consultant Specialty: Use inpatient specialty AT TIME OF CULTURE if available, otherwise use outpatient
    # Note: ConsultantMainSpecialtyCode_AtCulture is from the episode active when culture was taken (avoids leakage)
    ConsultantMainSpecialtyCode_Combined = coalesce(ConsultantMainSpecialtyCode_AtCulture, OutpatAppointment_SpecialtyCode),
    
    # Consultant Specialty: Full name (with specific merges applied)
    ConsultantSpecialtyName_Raw = get_specialty_name(ConsultantMainSpecialtyCode_Combined),
    
    # Consolidate specialties: merge paediatrics, obs/gyn, and small specialties
    ConsultantSpecialtyName = case_when(
      # Merge all paediatric specialties
      ConsultantSpecialtyName_Raw %in% c("Paediatrics", "Paediatric Surgery", 
                                         "Paediatric Cardiology", "Paediatric Neurology",
                                         "Paediatric Dentistry") ~ "Paediatrics (all)",
      
      # Merge obstetrics, gynaecology, and midwifery
      ConsultantSpecialtyName_Raw %in% c("Obstetrics", "Gynaecology", "Midwifery") ~ "Obstetrics & Gynaecology",
      
      # Consolidate small specialties (<100 cases) into "Other"
      ConsultantSpecialtyName_Raw %in% c("Ophthalmology", "Radiology", 
                                         "Clinical Immunology and Allergy", "Nursing",
                                         "Allied Health Professional", "Dermatology",
                                         "Clinical Neurophysiology", "Clinical Genetics",
                                         "Blood Transfusion", "Medical Ophthalmology",
                                         "Histopathology", "Paediatric Dentistry",
                                         "Old Age Psychiatry", "Unknown", "Orthodontics") ~ "Other",
      
      # Keep all other specialties as is
      TRUE ~ ConsultantSpecialtyName_Raw
    ),
    
    # Create composite ICU status for easier interpretation
    ICU_Status = case_when(
      ICU_at_Culture == 1 ~ "ICU at culture",
      ICU_after_Culture == 1 ~ "ICU in 24h after culture",
      ICU_any_time == 1 ~ "ICU other time",
      TRUE ~ "No ICU"
    )
  ) %>%
  select(-Has_ICU_Ward, -ConsultantSpecialtyName_Raw, -ConsultantMainSpecialtyCode_Combined)

# ============================================================
# 4. SAVE OUTPUT
# ============================================================

fwrite(bc_df, "Output Files/CSVs/blood_cultures_with_admissions.csv")

# ============================================================
# 5. SUMMARY STATISTICS
# ============================================================

# CHECKPOINT: Verify no blood cultures were lost
final_bc_count <- nrow(bc_df)

if (initial_bc_count != final_bc_count) {
  warning("Blood culture count mismatch! Lost: ", initial_bc_count - final_bc_count, " blood cultures")
}

# Create admission summary table
summary_table <- bc_df %>%
  group_by(CultureCategory) %>%
  summarise(
    N_Cultures = n(),
    
    # Culture Location
    ED_N = sum(CultureLocation == "ED", na.rm = TRUE),
    ED_Pct = round(100 * mean(CultureLocation == "ED", na.rm = TRUE), 1),
    Inpatient_N = sum(CultureLocation == "Inpatient Ward", na.rm = TRUE),
    Inpatient_Pct = round(100 * mean(CultureLocation == "Inpatient Ward", na.rm = TRUE), 1),
    Outpatient_N = sum(CultureLocation == "Outpatient", na.rm = TRUE),
    Outpatient_Pct = round(100 * mean(CultureLocation == "Outpatient", na.rm = TRUE), 1),
    Unknown_N = sum(CultureLocation == "Unknown", na.rm = TRUE),
    Unknown_Pct = round(100 * mean(CultureLocation == "Unknown", na.rm = TRUE), 1),
    
    # Admission Timing
    NotAdmitted_N = sum(AdmissionTiming == "Not admitted", na.rm = TRUE),
    NotAdmitted_Pct = round(100 * mean(AdmissionTiming == "Not admitted", na.rm = TRUE), 1),
    Admitted_0_48h_N = sum(AdmissionTiming == "Already admitted 0-48h", na.rm = TRUE),
    Admitted_0_48h_Pct = round(100 * mean(AdmissionTiming == "Already admitted 0-48h", na.rm = TRUE), 1),
    Admitted_gt48h_N = sum(AdmissionTiming == "Already admitted >48h", na.rm = TRUE),
    Admitted_gt48h_Pct = round(100 * mean(AdmissionTiming == "Already admitted >48h", na.rm = TRUE), 1),
    
    # Community vs Hospital Onset
    CommunityOnset_N = sum(CommunityOnset == "Community-onset", na.rm = TRUE),
    CommunityOnset_Pct = round(100 * mean(CommunityOnset == "Community-onset", na.rm = TRUE), 1),
    HospitalOnset_N = sum(CommunityOnset == "Hospital-onset", na.rm = TRUE),
    HospitalOnset_Pct = round(100 * mean(CommunityOnset == "Hospital-onset", na.rm = TRUE), 1),
    
    # Admission Method
    Emergency_N = sum(AdmissionMethod == "Emergency", na.rm = TRUE),
    Emergency_Pct = round(100 * mean(AdmissionMethod == "Emergency", na.rm = TRUE), 1),
    Elective_N = sum(AdmissionMethod == "Elective", na.rm = TRUE),
    Elective_Pct = round(100 * mean(AdmissionMethod == "Elective", na.rm = TRUE), 1),
    
    # ICU
    ICU_at_Culture_N = sum(ICU_at_Culture == 1, na.rm = TRUE),
    ICU_at_Culture_Pct = round(100 * mean(ICU_at_Culture == 1, na.rm = TRUE), 1),
    ICU_after_N = sum(ICU_after_Culture == 1, na.rm = TRUE),
    ICU_after_Pct = round(100 * mean(ICU_after_Culture == 1, na.rm = TRUE), 1),
    ICU_any_N = sum(ICU_any_time == 1, na.rm = TRUE),
    ICU_any_Pct = round(100 * mean(ICU_any_time == 1, na.rm = TRUE), 1),
    
    # Length of Stay
    LOS_Median = round(median(LengthOfStay_days, na.rm = TRUE), 1),
    LOS_Mean = round(mean(LengthOfStay_days, na.rm = TRUE), 1),
    
    .groups = "drop"
  ) %>%
  arrange(desc(CultureCategory))

# Add total row
total_row <- bc_df %>%
  summarise(
    CultureCategory = "Total",
    N_Cultures = n(),
    ED_N = sum(CultureLocation == "ED", na.rm = TRUE),
    ED_Pct = round(100 * mean(CultureLocation == "ED", na.rm = TRUE), 1),
    Inpatient_N = sum(CultureLocation == "Inpatient Ward", na.rm = TRUE),
    Inpatient_Pct = round(100 * mean(CultureLocation == "Inpatient Ward", na.rm = TRUE), 1),
    Outpatient_N = sum(CultureLocation == "Outpatient", na.rm = TRUE),
    Outpatient_Pct = round(100 * mean(CultureLocation == "Outpatient", na.rm = TRUE), 1),
    Unknown_N = sum(CultureLocation == "Unknown", na.rm = TRUE),
    Unknown_Pct = round(100 * mean(CultureLocation == "Unknown", na.rm = TRUE), 1),
    NotAdmitted_N = sum(AdmissionTiming == "Not admitted", na.rm = TRUE),
    NotAdmitted_Pct = round(100 * mean(AdmissionTiming == "Not admitted", na.rm = TRUE), 1),
    Admitted_0_48h_N = sum(AdmissionTiming == "Already admitted 0-48h", na.rm = TRUE),
    Admitted_0_48h_Pct = round(100 * mean(AdmissionTiming == "Already admitted 0-48h", na.rm = TRUE), 1),
    Admitted_gt48h_N = sum(AdmissionTiming == "Already admitted >48h", na.rm = TRUE),
    Admitted_gt48h_Pct = round(100 * mean(AdmissionTiming == "Already admitted >48h", na.rm = TRUE), 1),
    CommunityOnset_N = sum(CommunityOnset == "Community-onset", na.rm = TRUE),
    CommunityOnset_Pct = round(100 * mean(CommunityOnset == "Community-onset", na.rm = TRUE), 1),
    HospitalOnset_N = sum(CommunityOnset == "Hospital-onset", na.rm = TRUE),
    HospitalOnset_Pct = round(100 * mean(CommunityOnset == "Hospital-onset", na.rm = TRUE), 1),
    Emergency_N = sum(AdmissionMethod == "Emergency", na.rm = TRUE),
    Emergency_Pct = round(100 * mean(AdmissionMethod == "Emergency", na.rm = TRUE), 1),
    Elective_N = sum(AdmissionMethod == "Elective", na.rm = TRUE),
    Elective_Pct = round(100 * mean(AdmissionMethod == "Elective", na.rm = TRUE), 1),
    ICU_at_Culture_N = sum(ICU_at_Culture == 1, na.rm = TRUE),
    ICU_at_Culture_Pct = round(100 * mean(ICU_at_Culture == 1, na.rm = TRUE), 1),
    ICU_after_N = sum(ICU_after_Culture == 1, na.rm = TRUE),
    ICU_after_Pct = round(100 * mean(ICU_after_Culture == 1, na.rm = TRUE), 1),
    ICU_any_N = sum(ICU_any_time == 1, na.rm = TRUE),
    ICU_any_Pct = round(100 * mean(ICU_any_time == 1, na.rm = TRUE), 1),
    LOS_Median = round(median(LengthOfStay_days, na.rm = TRUE), 1),
    LOS_Mean = round(mean(LengthOfStay_days, na.rm = TRUE), 1)
  )

summary_table <- bind_rows(total_row, summary_table)

# Format for Excel output
summary_table_formatted <- summary_table %>% format_excel_output()

# ============================================================
# 6. WARD STATISTICS
# ============================================================

# Create ward summary - wards at time of culture
ward_stats_at_culture <- bc_wards_at_culture %>%
  filter(!is.na(WardName)) %>%
  group_by(WardName, ICU_Ward) %>%
  summarise(
    N_Patients = n_distinct(ClusterID),
    N_Cultures = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Ward_Type = if_else(ICU_Ward == 1, "ICU", "Non-ICU"),
    WardName_Display = if_else(ICU_Ward == 1, paste0("*** ", WardName, " ***"), WardName)
  ) %>%
  arrange(desc(ICU_Ward), desc(N_Cultures)) %>%
  select(WardName_Display, Ward_Type, N_Patients, N_Cultures)

# Overall ward statistics
overall_ward_stats <- tibble(
  Category = c("Total Wards", "ICU Wards", "Non-ICU Wards", 
               "Cultures with Ward Data", "Cultures in ICU", "Cultures in Non-ICU"),
  Count = c(
    n_distinct(ward_stats_at_culture$WardName_Display),
    sum(ward_stats_at_culture$Ward_Type == "ICU"),
    sum(ward_stats_at_culture$Ward_Type == "Non-ICU"),
    sum(ward_stats_at_culture$N_Cultures),
    sum(ward_stats_at_culture$N_Cultures[ward_stats_at_culture$Ward_Type == "ICU"]),
    sum(ward_stats_at_culture$N_Cultures[ward_stats_at_culture$Ward_Type == "Non-ICU"])
  )
)

# ============================================================
# 7. SPECIALTY BREAKDOWN TABLE
# ============================================================

# Create detailed specialty breakdown
specialty_breakdown <- bc_df %>%
  filter(!is.na(ConsultantSpecialtyName) & ConsultantSpecialtyName != "") %>%
  group_by(ConsultantSpecialtyName) %>%
  summarise(
    Total_N = n(),
    Positive_N = sum(CultureCategory == "Positive", na.rm = TRUE),
    Negative_N = sum(CultureCategory == "Negative (and SC)", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Positive_Pct = 100 * Positive_N / Total_N,
    Negative_Pct = 100 * Negative_N / Total_N
  ) %>%
  arrange(desc(Total_N)) %>%
  format_excel_output()

# Create clearer A&E comparison table
# Filter to ED cultures only for correct percentages
bc_ed <- bc_df %>% filter(CultureLocation == "ED")

# Placeholder: use bc_ed for comparison (actual comparison to all ED attendances would require additional data)
ed_comparison <- bc_ed

ae_comparison_table <- tibble(
  Metric = c("Total N", "Discharged from A&E - N", "Discharged from A&E - Pct", 
             "Admitted to hospital - N", "Admitted to hospital - Pct", "", 
             "Emergency admissions - N", "Emergency admissions - Pct", 
             "ICU any time - N", "ICU any time - Pct", 
             "Median LOS (days)", "Mean LOS (days)"),
  
  Blood_Cultures_in_AE = c(
    format_count(nrow(bc_ed)),
    format_count(sum(bc_ed$AdmissionTiming == "Not admitted", na.rm = TRUE)),
    format_percent(100 * mean(bc_ed$AdmissionTiming == "Not admitted", na.rm = TRUE)),
    format_count(sum(bc_ed$AdmissionTiming != "Not admitted", na.rm = TRUE)),
    format_percent(100 * mean(bc_ed$AdmissionTiming != "Not admitted", na.rm = TRUE)),
    "",
    format_count(sum(bc_ed$AdmissionMethod == "Emergency", na.rm = TRUE)),
    format_percent(100 * mean(bc_ed$AdmissionMethod == "Emergency", na.rm = TRUE)),
    format_count(sum(bc_ed$ICU_any_time == 1, na.rm = TRUE)),
    format_percent(100 * mean(bc_ed$ICU_any_time == 1, na.rm = TRUE)),
    format_number(median(bc_ed$LengthOfStay_days, na.rm = TRUE)),
    format_number(mean(bc_ed$LengthOfStay_days, na.rm = TRUE))
  ),
  
  All_AE_Attendances = c(
    format_count(nrow(ed_comparison)),
    format_count(sum(ed_comparison$AdmissionTiming == "Not admitted", na.rm = TRUE)),
    format_percent(100 * mean(ed_comparison$AdmissionTiming == "Not admitted", na.rm = TRUE)),
    format_count(sum(ed_comparison$AdmissionTiming != "Not admitted", na.rm = TRUE)),
    format_percent(100 * mean(ed_comparison$AdmissionTiming != "Not admitted", na.rm = TRUE)),
    "",
    format_count(sum(ed_comparison$AdmissionMethod == "Emergency", na.rm = TRUE)),
    format_percent(100 * mean(ed_comparison$AdmissionMethod == "Emergency", na.rm = TRUE)),
    "N/A", "N/A", "N/A", "N/A"
  )
)

write_xlsx(list("Summary_Formatted" = summary_table_formatted, 
                "Summary_Raw_Data" = summary_table %>% format_excel_output(),
                "AE_Comparison" = ae_comparison_table,
                "Specialty_Breakdown" = specialty_breakdown,
                "Ward_Summary" = overall_ward_stats %>% format_excel_output(),
                "Ward_Details" = ward_stats_at_culture %>% format_excel_output(),
                "Formatting_Note" = create_formatting_note()),
           "Output Files/Excels/Pipeline/BC_3_Add_Admissions_excel.xlsx")

