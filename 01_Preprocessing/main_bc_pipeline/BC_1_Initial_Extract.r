library(data.table) 
library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(openxlsx)

# Load Excel output formatting helpers
source("01_Preprocessing/Excel_Output_Helper.r")

# ============================================================
# 1. READ AND PREPARE MICROBIOLOGY DATA
# ============================================================

# Use data.table's fread for fast loading of large CSV
micro_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_Microbiology.csv",
                  select = c("ClusterID", "AccessionNumber", "CollectionDateTime", "BatTestName", "BugCode", "BugName", "SpecimenFull")) %>%
  as_tibble()

# Filter and deduplicate 
cultures <- micro_df %>%
  filter(BatTestName == "BLOOD CULTURE", ignore.case = TRUE) %>%
  select(ClusterID, AccessionNumber, CollectionDateTime, SpecimenFull, BugCode, BugName) %>%
  mutate(CollectionDateTime = ymd_hms(CollectionDateTime)) %>%
  filter(!is.na(CollectionDateTime)) %>%  # Remove failed parses
  distinct(ClusterID, AccessionNumber, BugCode, .keep_all = TRUE) %>%
  arrange(ClusterID, CollectionDateTime, AccessionNumber, BugName)

# ============================================================
# 2. LOAD BUG GROUP MAPPING
# ============================================================

# Use data.table for fast keyed join
bug_lookup_dt <- fread("00_Data/lookup_tables/bug_group_table.csv")
setkey(bug_lookup_dt, BugCode)

# Map BugGroup using data.table (fast), then convert back
cultures_dt <- as.data.table(cultures)
cultures_dt[, BugGroup := bug_lookup_dt[cultures_dt, BugGroup, on = "BugCode"]]
cultures <- as_tibble(cultures_dt)

# ============================================================
# 2B. DEFINE BUG GROUP CONSTANTS (USED THROUGHOUT SCRIPT)
# ============================================================

# Contaminant groups (with and without NULL for different uses)
contaminant_only_groups <- c("VIRIDANS AND OTHER STREPS", "CoNS (CONTAMINANT)", "OTHER CONTAMINANT")
contaminant_groups_with_null <- c(contaminant_only_groups, "NULL")

# Pathogen groups
pathogen_groups <- c("S AUREUS", "E COLI", "KLEBSIELLA", "P AERUGINOSA", 
                     "ENTEROCOCCUS", "S PNEUMONIAE", "ENTEROBACTER", 
                     "BETA HAEMOLYTIC STREP", "OTHER ENTEROBACTERALES", 
                     "OTHER PATHOGENIC STREPS", "ANAEROBES", "CANDIDA",
                     "OTHER GRAM-POSITIVE", "OTHER GRAM-NEGATIVE")

# Gram-positive and Gram-negative subgroups
grampos_groups <- c("S AUREUS", "ENTEROCOCCUS", "S PNEUMONIAE", 
                    "BETA HAEMOLYTIC STREP", "OTHER PATHOGENIC STREPS", "OTHER GRAM-POSITIVE")
gramneg_groups <- c("E COLI", "KLEBSIELLA", "P AERUGINOSA", "ENTEROBACTER", "OTHER ENTEROBACTERALES", 
                    "OTHER GRAM-NEGATIVE")

# Non-gram pathogen subgroups (always pathogenic, not gram-stain classifiable)
fungi_groups   <- c("CANDIDA")
anaerobe_groups <- c("ANAEROBES")

# ============================================================
# 3. DEFINE EPISODES FOR CONTAMINANT DETECTION ONLY
# ============================================================
# Episodes are ONLY used to determine single vs recurrent contaminants
# Episode = cultures within 14 days of each other

cultures_dt <- as.data.table(cultures)
setorder(cultures_dt, ClusterID, CollectionDateTime)

# Calculate time since last culture
cultures_dt[, TimeSinceLastCulture := 
              round(as.numeric(difftime(CollectionDateTime,
                                       shift(CollectionDateTime, n = 1, type = "lag"),
                                       units = "days")), 2), 
            by = ClusterID]

# Define new episode: >14 days since last culture
cultures_dt[, NewEpisode := fifelse(is.na(TimeSinceLastCulture) | TimeSinceLastCulture > 14, 1L, 0L)]

# Assign episode numbers within each patient
cultures_dt[NewEpisode == 1, EpisodeNo := seq_len(.N), by = ClusterID]
cultures_dt <- setnafill(cultures_dt, "locf", cols = c("EpisodeNo"))

# Create unique episode ID across all patients (only used for contaminant detection)
cultures_dt[, EpisodeID := .GRP, by = .(ClusterID, EpisodeNo)]

# Convert back to tibble
cultures <- as_tibble(cultures_dt)

# ============================================================
# 4. CREATE WIDE FORMAT (ONE ROW PER BLOOD CULTURE)
# ============================================================

# Get base dataframe with one row per culture event
base_cultures <- cultures %>%
  distinct(ClusterID, AccessionNumber, CollectionDateTime, SpecimenFull) %>%
  arrange(ClusterID, CollectionDateTime, AccessionNumber)

# Pivot bugs wide (group by ClusterID + AccessionNumber)
positive_cultures <- cultures %>%
  filter(!is.na(BugName) & BugName != "" & BugName != "NULL") %>%
  group_by(ClusterID, AccessionNumber) %>%
  mutate(bug_number = row_number()) %>%
  ungroup() %>%
  select(ClusterID, AccessionNumber, CollectionDateTime, SpecimenFull, BugName, BugCode, BugGroup, bug_number) %>%
  pivot_wider(
    id_cols = c(ClusterID, AccessionNumber, CollectionDateTime, SpecimenFull),
    names_from = bug_number,
    values_from = c(BugName, BugCode, BugGroup),
    names_sep = "_"
  )

# Get EpisodeID for each culture (only needed for contaminant detection)
episode_info <- cultures %>%
  group_by(ClusterID, AccessionNumber) %>%
  slice(1) %>%
  ungroup() %>%
  select(ClusterID, AccessionNumber, CollectionDateTime, SpecimenFull, EpisodeID)

# Join everything together
primary_blood_cultures <- base_cultures %>%
  left_join(episode_info, by = c("ClusterID", "AccessionNumber", "CollectionDateTime", "SpecimenFull")) %>%
  left_join(positive_cultures, by = c("ClusterID", "AccessionNumber", "CollectionDateTime", "SpecimenFull")) %>%
  mutate(
    CultureResult = if_else(!is.na(BugName_1) & BugName_1 != "NULL" & BugName_1 != "", 
                           "Positive", "Negative")
  )

# ============================================================
# CREATE UNIQUE COMPOSITE KEY FOR EACH BLOOD CULTURE
# ============================================================
# AccessionNumber alone is NOT unique - it can be reused across different patients
# Even ClusterID + AccessionNumber is not always unique (6 cases in dataset)
# Create stable composite key: ClusterID + AccessionNumber + CollectionDateTime
# This ensures each blood culture has a unique identifier throughout pipeline
primary_blood_cultures <- primary_blood_cultures %>%
  mutate(BloodCultureID = paste0(ClusterID, "_", AccessionNumber, "_", 
                                  format(CollectionDateTime, "%Y%m%d_%H%M%S"))) %>%
  arrange(ClusterID, CollectionDateTime, AccessionNumber)

# ============================================================
# 5. CALCULATE CULTURE-LEVEL SINGLE CONTAMINANT FLAG (OPTIMIZED)
# ============================================================
# For each culture that only grows contaminants, check if those same
# contaminant species appear in any other culture in the EPISODE (14-day window)
# OPTIMIZED: Single pivot, data.table operations, BloodCultureID joins

# Get all BugName and BugGroup columns
bugname_cols <- grep("^BugName_[0-9]+$", names(primary_blood_cultures), value = TRUE)
buggroup_cols <- grep("^BugGroup_[0-9]+$", names(primary_blood_cultures), value = TRUE)

# Single pivot to get all organisms in long format
culture_organisms <- primary_blood_cultures %>%
  select(BloodCultureID, EpisodeID, all_of(bugname_cols), all_of(buggroup_cols)) %>%
  pivot_longer(
    cols = c(all_of(bugname_cols), all_of(buggroup_cols)),
    names_to = c(".value", "bug_number"),
    names_pattern = "(BugName|BugGroup)_(\\d+)"
  ) %>%
  filter(!is.na(BugName) & BugName != "" & BugName != "NULL")

# Convert to data.table for efficient operations
setDT(culture_organisms)

# For each culture, determine if it's contaminant-only
culture_type <- culture_organisms[, .(
  only_contaminants = all(BugGroup %in% contaminant_only_groups)
), by = BloodCultureID]

# For contaminant-only cultures, check if contaminants recur within episode
# Count how many cultures each contaminant species appears in per episode
contaminant_recurrence <- culture_organisms[BugGroup %in% contaminant_only_groups, 
  .(n_cultures = uniqueN(BloodCultureID)), 
  by = .(EpisodeID, BugName)]

# Join recurrence counts back to cultures, then check if ALL contaminants are single
contaminant_status <- culture_organisms[BugGroup %in% contaminant_only_groups][
  contaminant_recurrence, on = .(EpisodeID, BugName)][,
  .(all_contaminants_single = all(n_cultures == 1)),
  by = BloodCultureID]

# Determine SingleContaminant flag: only contaminants AND all are single (non-recurrent)
single_contaminant_flag <- culture_type[contaminant_status, on = "BloodCultureID"][
  only_contaminants == TRUE & all_contaminants_single == TRUE,
  .(BloodCultureID, SingleContaminant = 1L)]

# Join back to primary blood cultures
primary_blood_cultures <- primary_blood_cultures %>%
  left_join(as_tibble(single_contaminant_flag), by = "BloodCultureID") %>%
  mutate(SingleContaminant = if_else(is.na(SingleContaminant), 0L, SingleContaminant))

# ============================================================
# 5B. ADD COLLECTION TO RECEIPT TIME
# ============================================================
# Calculate time from blood culture collection to receipt in lab (hours)
# OPTIMIZED: Single pass through microbiology data, data.table join

# Load ReceiveDateTime from microbiology data
receipt_times <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_Microbiology.csv",
                      select = c("AccessionNumber", "CollectionDateTime", "ReceiveDateTime"))

# Parse datetime columns and calculate time difference
receipt_times[, `:=`(
  CollectionDateTime = ymd_hms(CollectionDateTime),
  ReceiveDateTime = ymd_hms(ReceiveDateTime)
)]

# Remove invalid datetimes and calculate hours
receipt_times <- receipt_times[!is.na(CollectionDateTime) & !is.na(ReceiveDateTime)]
receipt_times[, CollectionToReceipt_Hours := as.numeric(difftime(ReceiveDateTime, CollectionDateTime, units = "hours"))]

# Remove impossible values (negative times or >48 hours - likely data errors)
receipt_times <- receipt_times[CollectionToReceipt_Hours >= 0 & CollectionToReceipt_Hours <= 48]

# Get one value per AccessionNumber (in case of duplicates, take first)
receipt_times_unique <- receipt_times[, .SD[1], by = AccessionNumber][
  , .(AccessionNumber, CollectionToReceipt_Hours)]

# Join to primary blood cultures
primary_blood_cultures <- primary_blood_cultures %>%
  left_join(as_tibble(receipt_times_unique), by = "AccessionNumber")

# ============================================================
# 6. ADD CULTIVABILITY CATEGORIZATION
# ============================================================
# Categorize cultures by cultivability using culture-level contaminant flags

# Categorize cultures based on individual culture contents using vectorized data.table operations
pbc_dt <- as.data.table(primary_blood_cultures)

# Vectorized operations: check each BugGroup column against our defined groups
for (col in buggroup_cols) {
  pbc_dt[, paste0(col, "_pathogen") := get(col) %in% pathogen_groups]
  pbc_dt[, paste0(col, "_contaminant") := get(col) %in% contaminant_only_groups]
  pbc_dt[, paste0(col, "_grampos") := get(col) %in% grampos_groups]
  pbc_dt[, paste0(col, "_gramneg") := get(col) %in% gramneg_groups]
  pbc_dt[, paste0(col, "_fungi") := get(col) %in% fungi_groups]
  pbc_dt[, paste0(col, "_anaerobe") := get(col) %in% anaerobe_groups]
  pbc_dt[, paste0(col, "_valid") := !is.na(get(col)) & get(col) != "" & get(col) != "NULL"]
}

# Aggregate across all BugGroup columns for each culture
pathogen_cols <- paste0(buggroup_cols, "_pathogen")
contaminant_cols <- paste0(buggroup_cols, "_contaminant")
grampos_cols <- paste0(buggroup_cols, "_grampos")
gramneg_cols <- paste0(buggroup_cols, "_gramneg")
fungi_cols <- paste0(buggroup_cols, "_fungi")
anaerobe_cols <- paste0(buggroup_cols, "_anaerobe")
valid_cols <- paste0(buggroup_cols, "_valid")

pbc_dt[, has_pathogen := rowSums(.SD, na.rm = TRUE) > 0, .SDcols = pathogen_cols]
pbc_dt[, has_contaminant := rowSums(.SD, na.rm = TRUE) > 0, .SDcols = contaminant_cols]
pbc_dt[, has_grampos := rowSums(.SD, na.rm = TRUE) > 0, .SDcols = grampos_cols]
pbc_dt[, has_gramneg := rowSums(.SD, na.rm = TRUE) > 0, .SDcols = gramneg_cols]
pbc_dt[, has_fungi := rowSums(.SD, na.rm = TRUE) > 0, .SDcols = fungi_cols]
pbc_dt[, has_anaerobe := rowSums(.SD, na.rm = TRUE) > 0, .SDcols = anaerobe_cols]
# Count only gram-classified pathogens (fungi/anaerobes alongside gram bugs are ignored)
pbc_dt[, n_gram_pathogens := rowSums(.SD, na.rm = TRUE), .SDcols = c(grampos_cols, gramneg_cols)]
pbc_dt[, n_bugs := rowSums(.SD, na.rm = TRUE), .SDcols = valid_cols]

# Categorize using vectorized fifelse (much faster than case_when on large data)
# Use culture-level logic: check if this specific culture's contaminants are single or recurrent
# Categories:
#  - Gram-positive/Gram-negative: any gram-classified pathogen(s) present; non-gram organisms alongside are ignored
#  - Polymicrobial: 2+ gram-classified pathogens (both gram-pos and gram-neg, or multiple of either)
#  - Fungi: Candida with no gram-classified bugs (always pathogenic)
#  - Anaerobes: anaerobic organisms with no gram-classified bugs (always pathogenic)
#  - Polymicrobial - contaminant: n_bugs >= 2, no gram pathogens, contaminants only, no pathogens
#  - Contaminant - single/recurrent: single contaminant, episode-level recurrence determines subtype
#  - Other: remaining non-gram pathogen edge cases (e.g. mixed fungi+anaerobes without gram bugs)
pbc_dt[, CultivabilityCategory := 
  fifelse(n_bugs == 0, "Culture-negative",
  # Gram-classified pathogens (non-gram organisms alongside are ignored)
  fifelse(has_grampos & !has_gramneg & n_gram_pathogens == 1, "Gram-positive",
  fifelse(has_gramneg & !has_grampos & n_gram_pathogens == 1, "Gram-negative",
  # Polymicrobial: requires 2+ gram-classified pathogens
  fifelse((n_gram_pathogens > 1 | (has_grampos & has_gramneg)) & (has_grampos | has_gramneg), "Polymicrobial",
  # Non-gram pathogens: fungi and anaerobes (no gram-classified bugs)
  fifelse(has_fungi & !has_grampos & !has_gramneg, "Fungi",
  fifelse(has_anaerobe & !has_grampos & !has_gramneg, "Anaerobes",
  fifelse(!has_grampos & !has_gramneg & n_bugs >= 2 & has_contaminant & !has_pathogen, "Polymicrobial - contaminant",
  # Single contaminants (episode-level recurrence determines single vs recurrent)
  fifelse(!has_pathogen & has_contaminant & SingleContaminant == 1, "Contaminant - single",
  fifelse(!has_pathogen & has_contaminant & SingleContaminant == 0, "Contaminant - recurrent",
  # Remaining edge cases (e.g. mixed fungi+anaerobes without gram bugs)
  "Other")))))))))]

# Create simplified binary categorization:
# Positive = gram pathogens, fungi, anaerobes, recurrent contaminants, and other non-gram pathogens
# Negative (and SC) = culture-negative, single contaminants, and polymicrobial contaminant
pbc_dt[, CultureCategory := 
  fifelse(CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial",
                                        "Fungi", "Anaerobes", "Other", "Contaminant - recurrent"),
          "Positive",
          "Negative (and SC)")]

# Clean up temporary columns
temp_cols <- c(paste0(buggroup_cols, "_pathogen"),
               paste0(buggroup_cols, "_contaminant"),
               paste0(buggroup_cols, "_grampos"),
               paste0(buggroup_cols, "_gramneg"),
               paste0(buggroup_cols, "_fungi"),
               paste0(buggroup_cols, "_anaerobe"),
               paste0(buggroup_cols, "_valid"),
               "has_pathogen", "has_contaminant", "has_grampos", "has_gramneg",
               "has_fungi", "has_anaerobe", "n_bugs")
pbc_dt[, (temp_cols) := NULL]

# ============================================================
# 6B. CATEGORIZE CULTURE SOURCE
# ============================================================
# Categorize SpecimenFull into standardized categories while keeping original data
# Categories: Peripheral, Central Line, PICC, Port/Portacath, Arterial, Femoral, Unspecified/Other

# Convert to uppercase for case-insensitive matching
pbc_dt[, SpecimenFull_upper := toupper(SpecimenFull)]

# Categorize using vectorized pattern matching (data.table is very fast for this)
pbc_dt[, CultureSourceCategory := fcase(
  # 1. Peripheral - most common, check first
  grepl("PERIPHERAL", SpecimenFull_upper) & 
    !grepl("STEM CELL|PBSC|HPC", SpecimenFull_upper), "Peripheral",
  
  # 2. Port/Portacath - check before general "LINE" patterns
  grepl("PORT|PORTACATH", SpecimenFull_upper) & 
    !grepl("CENTRAL LINE.*PORT", SpecimenFull_upper), "Port/Portacath",
  
  # 3. PICC - check for PIC LINE or PICC LINE (but not PORTACATH)
  grepl("PIC LINE|PICC LINE|PIC$|PICC$", SpecimenFull_upper) & 
    !grepl("PORTACATH", SpecimenFull_upper), "PICC",
  
  # 4. Arterial
  grepl("ARTERIAL", SpecimenFull_upper), "Arterial",
  
  # 5. Femoral
  grepl("FEMORAL", SpecimenFull_upper), "Femoral",
  
  # 6. Central Line - check after PICC and Port
  grepl("CENTRAL LINE|CENTRAL  LINE|CENRAL LINE|CVC|HICKMAN|CLI ", SpecimenFull_upper), "Central Line",
  
  # 7. Unspecified/Other - catch all remaining
  default = "Unspecified/Other"
)]

# Clean up temporary column
pbc_dt[, SpecimenFull_upper := NULL]

# Reorder columns to put BloodCultureID first for easy identification
# (BloodCultureID was created earlier in step 4 for efficient contaminant joins)
setcolorder(pbc_dt, c("BloodCultureID", setdiff(names(pbc_dt), "BloodCultureID")))

# Convert back to tibble
primary_blood_cultures <- as_tibble(pbc_dt)

# ============================================================
# 6C. COUNT BLOOD CULTURES IN PRIOR 28 DAYS (EXCLUDING <24 HOURS)
# ============================================================
# For each blood culture, count how many blood cultures occurred in the 
# 28 days prior for the same patient (ClusterID), EXCLUDING cultures within 24 hours
# OPTIMIZED: Using rolling join instead of nested loops (O(n log n) vs O(n²))

pbc_dt <- as.data.table(primary_blood_cultures)
setorder(pbc_dt, ClusterID, CollectionDateTime)

# Create rolling window: for each culture, find all cultures in prior 28 days (but >24 hours ago)
# Using data.table's binary search with rolling join
pbc_dt[, datetime_minus_28 := CollectionDateTime - days(28)]
pbc_dt[, datetime_minus_1 := CollectionDateTime - days(1)]
pbc_dt[, row_id := .I]

# Self-join: for each row, count how many rows fall in the 28-day window (excluding <24 hours)
# This is O(n log n) instead of O(n²)
count_dt <- pbc_dt[pbc_dt, 
                   on = .(ClusterID, CollectionDateTime >= datetime_minus_28, CollectionDateTime < datetime_minus_1),
                   .N,
                   by = .EACHI]
setnames(count_dt, "N", "CulturesInPrior28Days")

# Merge counts back
pbc_dt[, CulturesInPrior28Days := count_dt$CulturesInPrior28Days]

# Clean up temporary columns
pbc_dt[, c("datetime_minus_28", "datetime_minus_1", "row_id") := NULL]

# Convert back to tibble
primary_blood_cultures <- as_tibble(pbc_dt)

# ============================================================
# 6D. COUNT POSITIVE CULTURES & CULTIVABILITY CATEGORIES IN PRIOR 28 DAYS (EXCLUDING <24 HOURS)
# ============================================================
# For each blood culture, count positive cultures and specific cultivability
# categories in the 28 days prior for the same patient (ClusterID), EXCLUDING cultures within 24 hours
# OPTIMIZED: Using rolling joins instead of nested loops (O(n log n) vs O(n²))

pbc_dt <- as.data.table(primary_blood_cultures)
setorder(pbc_dt, ClusterID, CollectionDateTime)

# Pre-compute indicator columns for each category (vectorized operations)
# NOTE: is_positive_outcome includes all pathogenic organisms (gram-classified, fungi, anaerobes)
pbc_dt[, `:=`(
  is_positive_outcome = (CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial", "Fungi", "Anaerobes", "Other")),
  is_culture_negative = (CultivabilityCategory == "Culture-negative"),
  is_gram_positive = (CultivabilityCategory == "Gram-positive"),
  is_gram_negative = (CultivabilityCategory == "Gram-negative"),
  is_polymicrobial = (CultivabilityCategory == "Polymicrobial"),
  is_fungi = (CultivabilityCategory == "Fungi"),
  is_anaerobe = (CultivabilityCategory == "Anaerobes"),
  is_other = (CultivabilityCategory == "Other"),
  is_contam_recurrent = (CultivabilityCategory == "Contaminant - recurrent"),
  is_contam_single = (CultivabilityCategory == "Contaminant - single"),
  is_polymicrobial_contaminant = (CultivabilityCategory == "Polymicrobial - contaminant")
)]

# Create rolling window (excluding <24 hours)
pbc_dt[, datetime_minus_28 := CollectionDateTime - days(28)]
pbc_dt[, datetime_minus_1 := CollectionDateTime - days(1)]

# Use rolling joins to count each category
# This performs binary search instead of full scan for each row
for (category_info in list(
  list(col = "is_positive_outcome", name = "PositiveCulturesInPrior28Days"),
  list(col = "is_culture_negative", name = "CultureNegativeCulturesInPrior28Days"),
  list(col = "is_gram_positive", name = "GramPositiveCulturesInPrior28Days"),
  list(col = "is_gram_negative", name = "GramNegativeCulturesInPrior28Days"),
  list(col = "is_polymicrobial", name = "PolymicrobialCulturesInPrior28Days"),
  list(col = "is_fungi", name = "FungiCulturesInPrior28Days"),
  list(col = "is_anaerobe", name = "AnaerobesCulturesInPrior28Days"),
  list(col = "is_other", name = "OtherCulturesInPrior28Days"),
  list(col = "is_contam_recurrent", name = "RecurrentContaminantCulturesInPrior28Days"),
  list(col = "is_contam_single", name = "SingleContaminantCulturesInPrior28Days"),
  list(col = "is_polymicrobial_contaminant", name = "PolymicrobialContaminantCulturesInPrior28Days")
)) {
  
  # Self-join with rolling window to sum the indicator column (excluding <24 hours)
  # Note: Do NOT use "i." prefix - we want LEFT table (prior cultures), not RIGHT table (current culture)
  count_dt <- pbc_dt[pbc_dt, 
                     on = .(ClusterID, CollectionDateTime >= datetime_minus_28, CollectionDateTime < datetime_minus_1),
                     .(count = sum(get(category_info$col), na.rm = TRUE)),
                     by = .EACHI]
  
  # Assign back to main table (use 0 for NA values - no prior matches)
  pbc_dt[, (category_info$name) := ifelse(is.na(count_dt$count), 0, count_dt$count)]
}

# ============================================================
# 6E. COUNT ORGANISM-SPECIFIC CULTURES IN PRIOR 28 DAYS
# ============================================================
# For each blood culture, count specific organisms in prior 28 days (excluding <24 hours)
# Gram-Positive: S aureus, Enterococcus, Streptococcus (pathogenic), Other GP
# Gram-Negative: Enterobacterales, Pseudomonas, Other GN (incl. Haemophilus, Neisseria)

# Get BugGroup column names
buggroup_cols <- grep("^BugGroup_[0-9]+$", names(pbc_dt), value = TRUE)

# Define organism categories
grampos_s_aureus <- c("S AUREUS")
grampos_enterococcus <- c("ENTEROCOCCUS")
grampos_streptococcus <- c("S PNEUMONIAE", "BETA HAEMOLYTIC STREP", "OTHER PATHOGENIC STREPS")
grampos_other <- c("OTHER GRAM-POSITIVE")

gramneg_enterobacterales <- c("E COLI", "KLEBSIELLA", "ENTEROBACTER", "OTHER ENTEROBACTERALES")
gramneg_pseudomonas <- c("P AERUGINOSA")
gramneg_other <- c("OTHER GRAM-NEGATIVE")

# Helper: Check if ANY BugGroup column contains organisms from a list
has_organism <- function(dt, organism_list, buggroup_cols) {
  result <- rep(FALSE, nrow(dt))
  for (col in buggroup_cols) {
    if (col %in% names(dt)) {
      result <- result | (dt[[col]] %in% organism_list)
    }
  }
  return(as.integer(result))
}

# Create indicator columns for each organism category
pbc_dt[, `:=`(
  is_prior_s_aureus = has_organism(pbc_dt, grampos_s_aureus, buggroup_cols),
  is_prior_enterococcus = has_organism(pbc_dt, grampos_enterococcus, buggroup_cols),
  is_prior_streptococcus = has_organism(pbc_dt, grampos_streptococcus, buggroup_cols),
  is_prior_other_grampos = has_organism(pbc_dt, grampos_other, buggroup_cols),
  is_prior_enterobacterales = has_organism(pbc_dt, gramneg_enterobacterales, buggroup_cols),
  is_prior_pseudomonas = has_organism(pbc_dt, gramneg_pseudomonas, buggroup_cols),
  is_prior_other_gramneg = has_organism(pbc_dt, gramneg_other, buggroup_cols)
)]

# Count organism-specific cultures in prior 28 days (excluding <24 hours)
for (organism_info in list(
  list(col = "is_prior_s_aureus", name = "Prior_S_aureus_Count"),
  list(col = "is_prior_enterococcus", name = "Prior_Enterococcus_Count"),
  list(col = "is_prior_streptococcus", name = "Prior_Streptococcus_Count"),
  list(col = "is_prior_other_grampos", name = "Prior_Other_GramPos_Count"),
  list(col = "is_prior_enterobacterales", name = "Prior_Enterobacterales_Count"),
  list(col = "is_prior_pseudomonas", name = "Prior_Pseudomonas_Count"),
  list(col = "is_prior_other_gramneg", name = "Prior_Other_GramNeg_Count")
)) {
  
  count_dt <- pbc_dt[pbc_dt, 
                     on = .(ClusterID, CollectionDateTime >= datetime_minus_28, CollectionDateTime < datetime_minus_1),
                     .(count = sum(get(organism_info$col), na.rm = TRUE)),
                     by = .EACHI]
  
  pbc_dt[, (organism_info$name) := ifelse(is.na(count_dt$count), 0, count_dt$count)]
}

# Clean up temporary columns
pbc_dt[, c("datetime_minus_28", "datetime_minus_1",
           "is_positive_outcome", "is_culture_negative", "is_gram_positive", "is_gram_negative",
           "is_polymicrobial", "is_fungi", "is_anaerobe", "is_other",
           "is_contam_recurrent", "is_contam_single", "is_polymicrobial_contaminant",
           "is_prior_s_aureus", "is_prior_enterococcus", "is_prior_streptococcus", "is_prior_other_grampos",
           "is_prior_enterobacterales", "is_prior_pseudomonas", "is_prior_other_gramneg") := NULL]

# Convert back to tibble
primary_blood_cultures <- as_tibble(pbc_dt)

# ============================================================
# 6b. TEMPORAL FEATURES
# ============================================================
# Extract temporal patterns from CollectionDateTime

# Ensure CollectionDateTime is in POSIXct format
primary_blood_cultures$CollectionDateTime <- ymd_hms(primary_blood_cultures$CollectionDateTime, quiet = TRUE)

# Extract hour of day (0-23) and day of week (1=Sunday, 7=Saturday)
primary_blood_cultures$Collection_Hour <- hour(primary_blood_cultures$CollectionDateTime)
primary_blood_cultures$Collection_DayOfWeek <- wday(primary_blood_cultures$CollectionDateTime)

# Out-of-hours indicator: 6pm to 8am on weekdays OR anytime on weekends (Saturday=7, Sunday=1)
# Treats missing as in-hours (conservative)
primary_blood_cultures$OutOfHours <- ifelse(
  !is.na(primary_blood_cultures$Collection_Hour) & 
  !is.na(primary_blood_cultures$Collection_DayOfWeek) & 
  (primary_blood_cultures$Collection_DayOfWeek %in% c(1, 7) |  # Weekend (Sunday=1, Saturday=7)
   primary_blood_cultures$Collection_Hour < 8 |                 # Before 8am on weekdays
   primary_blood_cultures$Collection_Hour >= 18),               # After 6pm on weekdays
  1, 0)
primary_blood_cultures$OutOfHours[is.na(primary_blood_cultures$Collection_Hour)] <- 0

# ============================================================
# 7. SAVE OUTPUT
# ============================================================

fwrite(primary_blood_cultures, "Output Files/CSVs/primary_blood_cultures.csv")

# ============================================================
# 8. GENERATE SUMMARY STATISTICS
# ============================================================

# Helper function to calculate N and percentage
calc_stat <- function(condition) {
  n <- sum(condition, na.rm = TRUE)
  pct <- round(100 * mean(condition, na.rm = TRUE), 1)
  c(n, pct)
}

# --- SECTION 1: Overall Culture Summary ---
total_cultures <- nrow(primary_blood_cultures)
positive_cultures <- sum(primary_blood_cultures$CultureResult == "Positive", na.rm = TRUE)
negative_cultures <- sum(primary_blood_cultures$CultureResult == "Negative", na.rm = TRUE)
single_contaminants <- sum(primary_blood_cultures$SingleContaminant == 1, na.rm = TRUE)
positive_non_contam <- sum(primary_blood_cultures$CultureResult == "Positive" & 
                           primary_blood_cultures$SingleContaminant == 0, na.rm = TRUE)
# Model outcomes (gram-classified + fungi + anaerobes)
model_positive <- sum(primary_blood_cultures$CultivabilityCategory %in% c("Gram-positive", "Gram-negative", "Polymicrobial", "Fungi", "Anaerobes", "Other"), na.rm = TRUE)

overall_summary <- tibble(
  Category = c("All Cultures", 
               "Positive Cultures", 
               "Positive (excl. single contaminants)",
               "Positive (model outcome: Gram-pos + Gram-neg + Polymicrobial + Fungi + Anaerobes)",
               "Negative Cultures", 
               "Single Contaminant Cultures"),
  N = c(total_cultures, 
        positive_cultures,
        positive_non_contam,
        model_positive,
        negative_cultures, 
        single_contaminants),
  Percentage = c(100.0,
                 100 * positive_cultures / total_cultures,
                 100 * positive_non_contam / total_cultures,
                 100 * model_positive / total_cultures,
                 100 * negative_cultures / total_cultures,
                 100 * single_contaminants / total_cultures)
)

# --- SECTION 2: Cultivability Category Breakdown ---
cultivability_summary <- primary_blood_cultures %>%
  group_by(CultivabilityCategory) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(Percentage = round(100 * N / sum(N), 1)) %>%
  arrange(desc(N)) %>%
  rename(Category = CultivabilityCategory)

# Add aggregate rows for gram positive and gram negative (model outcomes)
grampos_count <- sum(primary_blood_cultures$CultivabilityCategory == "Gram-positive", na.rm = TRUE)
gramneg_count <- sum(primary_blood_cultures$CultivabilityCategory == "Gram-negative", na.rm = TRUE)

cultivability_summary <- bind_rows(
  cultivability_summary,
  tibble(Category = "", N = NA_integer_, Percentage = NA_real_),
  tibble(
    Category = "--- Model Outcomes ---",
    N = NA_integer_,
    Percentage = NA_real_
  ),
  tibble(
    Category = "Gram-positive (model outcome)",
    N = grampos_count,
    Percentage = round(100 * grampos_count / total_cultures, 1)
  ),
  tibble(
    Category = "Gram-negative (model outcome)",
    N = gramneg_count,
    Percentage = round(100 * gramneg_count / total_cultures, 1)
  )
)

# --- SECTION 3: Organism Type Breakdown (at Culture Level) ---
# This counts individual organism instances across all positive cultures
# If a culture has 2 organisms, it contributes 2 to the count

# Create long format of all organisms from positive cultures (using buggroup_cols from earlier)
organism_instances <- primary_blood_cultures %>%
  filter(CultureResult == "Positive") %>%
  select(ClusterID, AccessionNumber, CollectionDateTime, SpecimenFull, all_of(buggroup_cols)) %>%
  pivot_longer(cols = all_of(buggroup_cols), 
               names_to = "bug_position", 
               values_to = "BugGroup") %>%
  filter(!is.na(BugGroup) & BugGroup != "" & BugGroup != "NULL")

# Count organism instances
organism_summary <- organism_instances %>%
  group_by(BugGroup) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(Percentage = round(100 * N / sum(N), 1)) %>%
  arrange(desc(N)) %>%
  rename(Organism_Group = BugGroup)

# Add total row
organism_summary <- bind_rows(
  organism_summary,
  tibble(
    Organism_Group = "TOTAL (Organism Instances)",
    N = sum(organism_summary$N),
    Percentage = 100.0
  )
)

# --- SECTION 3b: Multiple Organisms Per Culture Analysis ---
# Check how many cultures have multiple organisms
primary_blood_cultures_analysis <- primary_blood_cultures %>%
  filter(CultureResult == "Positive")

# Count non-NA BugGroup columns per culture
# OPTIMIZED: Extract columns once and reuse for rowSums
bug_group_matrix <- primary_blood_cultures_analysis %>%
  select(all_of(buggroup_cols))

organism_count_per_culture <- primary_blood_cultures_analysis %>%
  mutate(
    n_organisms = rowSums(!is.na(bug_group_matrix) & 
                         bug_group_matrix != "" & 
                         bug_group_matrix != "NULL")
  )

multi_organism_summary <- organism_count_per_culture %>%
  group_by(n_organisms) %>%
  summarise(N_Cultures = n(), .groups = "drop") %>%
  mutate(Percentage = round(100 * N_Cultures / sum(N_Cultures), 1)) %>%
  arrange(n_organisms) %>%
  rename(Organisms_Per_Culture = n_organisms)

# Add summary statistics
multi_organism_summary <- bind_rows(
  multi_organism_summary,
  tibble(
    Organisms_Per_Culture = NA_integer_,
    N_Cultures = NA_integer_,
    Percentage = NA_real_
  ),
  tibble(
    Organisms_Per_Culture = NA_integer_,
    N_Cultures = sum(multi_organism_summary$N_Cultures),
    Percentage = 100.0
  )
) %>%
  mutate(Description = case_when(
    is.na(Organisms_Per_Culture) & is.na(N_Cultures) ~ "",
    is.na(Organisms_Per_Culture) ~ "TOTAL",
    TRUE ~ as.character(Organisms_Per_Culture)
  )) %>%
  select(Description, N_Cultures, Percentage)

# --- SECTION 4: Prior Culture History Analysis ---
# Distribution of prior culture counts in 28 days
prior_culture_distribution <- primary_blood_cultures %>%
  group_by(CulturesInPrior28Days) %>%
  summarise(N_Cultures = n(), .groups = "drop") %>%
  mutate(Percentage = round(100 * N_Cultures / sum(N_Cultures), 1)) %>%
  arrange(CulturesInPrior28Days) %>%
  rename(N_Prior_Cultures_in_28_Days = CulturesInPrior28Days)

# Add summary statistics
prior_culture_distribution <- bind_rows(
  prior_culture_distribution,
  tibble(
    N_Prior_Cultures_in_28_Days = NA_integer_,
    N_Cultures = NA_integer_,
    Percentage = NA_real_
  ),
  tibble(
    N_Prior_Cultures_in_28_Days = NA_integer_,
    N_Cultures = nrow(primary_blood_cultures),
    Percentage = 100.0
  )
) %>%
  mutate(Description = case_when(
    is.na(N_Prior_Cultures_in_28_Days) & is.na(N_Cultures) ~ "",
    is.na(N_Prior_Cultures_in_28_Days) ~ "TOTAL",
    TRUE ~ as.character(N_Prior_Cultures_in_28_Days)
  )) %>%
  select(Description, N_Cultures, Percentage)

# --- SECTION 4B: Prior Positive Culture Analysis ---
# Breakdown by CURRENT culture status: what % had prior positive/negative/no cultures?
# Mutually exclusive categories: positive, only negative, or no cultures (percentages add to 100%)
prior_positive_breakdown <- primary_blood_cultures %>%
  group_by(CultureCategory) %>%
  summarise(
    Total = n(),
    # Any prior positive culture (gram-pos, gram-neg, polymicrobial, other, or recurrent contaminant)
    N_with_Prior_Positive = sum(PositiveCulturesInPrior28Days > 0),
    Pct_with_Prior_Positive = round(100 * sum(PositiveCulturesInPrior28Days > 0) / n(), 1),
    # Only negative prior cultures (no positives, but had cultures)
    N_with_Only_Negative = sum(PositiveCulturesInPrior28Days == 0 & CultureNegativeCulturesInPrior28Days > 0),
    Pct_with_Only_Negative = round(100 * sum(PositiveCulturesInPrior28Days == 0 & CultureNegativeCulturesInPrior28Days > 0) / n(), 1),
    # No prior cultures done
    N_with_No_Cultures_Done = sum(CulturesInPrior28Days == 0),
    Pct_with_No_Cultures_Done = round(100 * sum(CulturesInPrior28Days == 0) / n(), 1),
    # Breakdown of positive prior cultures (not mutually exclusive, so shown separately)
    N_with_Prior_GramPos = sum(GramPositiveCulturesInPrior28Days > 0),
    Pct_with_Prior_GramPos = round(100 * sum(GramPositiveCulturesInPrior28Days > 0) / n(), 1),
    N_with_Prior_GramNeg = sum(GramNegativeCulturesInPrior28Days > 0),
    Pct_with_Prior_GramNeg = round(100 * sum(GramNegativeCulturesInPrior28Days > 0) / n(), 1),
    N_with_Prior_Polymicrobial = sum(PolymicrobialCulturesInPrior28Days > 0),
    Pct_with_Prior_Polymicrobial = round(100 * sum(PolymicrobialCulturesInPrior28Days > 0) / n(), 1),
    N_with_Prior_Other = sum(OtherCulturesInPrior28Days > 0),
    Pct_with_Prior_Other = round(100 * sum(OtherCulturesInPrior28Days > 0) / n(), 1),
    N_with_Prior_RecurContam = sum(RecurrentContaminantCulturesInPrior28Days > 0),
    Pct_with_Prior_RecurContam = round(100 * sum(RecurrentContaminantCulturesInPrior28Days > 0) / n(), 1),
    .groups = "drop"
  ) %>%
  arrange(CultureCategory)

# Add total row
prior_positive_breakdown <- bind_rows(
  prior_positive_breakdown,
  tibble(
    CultureCategory = "TOTAL",
    Total = sum(prior_positive_breakdown$Total),
    N_with_Prior_Positive = sum(prior_positive_breakdown$N_with_Prior_Positive),
    Pct_with_Prior_Positive = round(100 * sum(prior_positive_breakdown$N_with_Prior_Positive) / sum(prior_positive_breakdown$Total), 1),
    N_with_Only_Negative = sum(prior_positive_breakdown$N_with_Only_Negative),
    Pct_with_Only_Negative = round(100 * sum(prior_positive_breakdown$N_with_Only_Negative) / sum(prior_positive_breakdown$Total), 1),
    N_with_No_Cultures_Done = sum(prior_positive_breakdown$N_with_No_Cultures_Done),
    Pct_with_No_Cultures_Done = round(100 * sum(prior_positive_breakdown$N_with_No_Cultures_Done) / sum(prior_positive_breakdown$Total), 1),
    N_with_Prior_GramPos = sum(prior_positive_breakdown$N_with_Prior_GramPos),
    Pct_with_Prior_GramPos = round(100 * sum(prior_positive_breakdown$N_with_Prior_GramPos) / sum(prior_positive_breakdown$Total), 1),
    N_with_Prior_GramNeg = sum(prior_positive_breakdown$N_with_Prior_GramNeg),
    Pct_with_Prior_GramNeg = round(100 * sum(prior_positive_breakdown$N_with_Prior_GramNeg) / sum(prior_positive_breakdown$Total), 1),
    N_with_Prior_Polymicrobial = sum(prior_positive_breakdown$N_with_Prior_Polymicrobial),
    Pct_with_Prior_Polymicrobial = round(100 * sum(prior_positive_breakdown$N_with_Prior_Polymicrobial) / sum(prior_positive_breakdown$Total), 1),
    N_with_Prior_Other = sum(prior_positive_breakdown$N_with_Prior_Other),
    Pct_with_Prior_Other = round(100 * sum(prior_positive_breakdown$N_with_Prior_Other) / sum(prior_positive_breakdown$Total), 1),
    N_with_Prior_RecurContam = sum(prior_positive_breakdown$N_with_Prior_RecurContam),
    Pct_with_Prior_RecurContam = round(100 * sum(prior_positive_breakdown$N_with_Prior_RecurContam) / sum(prior_positive_breakdown$Total), 1)
  )
)

# Breakdown by CURRENT cultivability category: what % had prior cultures in each category?
# Mutually exclusive categories: positive, only negative, or no cultures (percentages add to 100%)
prior_positive_categories <- primary_blood_cultures %>%
  group_by(CultivabilityCategory) %>%
  summarise(
    Total = n(),
    # Any prior positive culture (gram-pos, gram-neg, polymicrobial, other, or recurrent contaminant)
    N_with_Prior_Positive = sum(PositiveCulturesInPrior28Days > 0),
    Pct_with_Prior_Positive = round(100 * sum(PositiveCulturesInPrior28Days > 0) / n(), 1),
    # Only negative prior cultures (no positives, but had cultures)
    N_with_Only_Negative = sum(PositiveCulturesInPrior28Days == 0 & CultureNegativeCulturesInPrior28Days > 0),
    Pct_with_Only_Negative = round(100 * sum(PositiveCulturesInPrior28Days == 0 & CultureNegativeCulturesInPrior28Days > 0) / n(), 1),
    # No prior cultures done
    N_with_No_Cultures_Done = sum(CulturesInPrior28Days == 0),
    Pct_with_No_Cultures_Done = round(100 * sum(CulturesInPrior28Days == 0) / n(), 1),
    # Breakdown of positive prior cultures (not mutually exclusive, so shown separately)
    N_with_Prior_GramPos = sum(GramPositiveCulturesInPrior28Days > 0),
    Pct_with_Prior_GramPos = round(100 * sum(GramPositiveCulturesInPrior28Days > 0) / n(), 1),
    N_with_Prior_GramNeg = sum(GramNegativeCulturesInPrior28Days > 0),
    Pct_with_Prior_GramNeg = round(100 * sum(GramNegativeCulturesInPrior28Days > 0) / n(), 1),
    N_with_Prior_Polymicrobial = sum(PolymicrobialCulturesInPrior28Days > 0),
    Pct_with_Prior_Polymicrobial = round(100 * sum(PolymicrobialCulturesInPrior28Days > 0) / n(), 1),
    N_with_Prior_Other = sum(OtherCulturesInPrior28Days > 0),
    Pct_with_Prior_Other = round(100 * sum(OtherCulturesInPrior28Days > 0) / n(), 1),
    N_with_Prior_RecurContam = sum(RecurrentContaminantCulturesInPrior28Days > 0),
    Pct_with_Prior_RecurContam = round(100 * sum(RecurrentContaminantCulturesInPrior28Days > 0) / n(), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(Total))

# Add total row
prior_positive_categories <- bind_rows(
  prior_positive_categories,
  tibble(
    CultivabilityCategory = "TOTAL",
    Total = sum(prior_positive_categories$Total),
    N_with_Prior_Positive = sum(prior_positive_categories$N_with_Prior_Positive),
    Pct_with_Prior_Positive = round(100 * sum(prior_positive_categories$N_with_Prior_Positive) / sum(prior_positive_categories$Total), 1),
    N_with_Only_Negative = sum(prior_positive_categories$N_with_Only_Negative),
    Pct_with_Only_Negative = round(100 * sum(prior_positive_categories$N_with_Only_Negative) / sum(prior_positive_categories$Total), 1),
    N_with_No_Cultures_Done = sum(prior_positive_categories$N_with_No_Cultures_Done),
    Pct_with_No_Cultures_Done = round(100 * sum(prior_positive_categories$N_with_No_Cultures_Done) / sum(prior_positive_categories$Total), 1),
    N_with_Prior_GramPos = sum(prior_positive_categories$N_with_Prior_GramPos),
    Pct_with_Prior_GramPos = round(100 * sum(prior_positive_categories$N_with_Prior_GramPos) / sum(prior_positive_categories$Total), 1),
    N_with_Prior_GramNeg = sum(prior_positive_categories$N_with_Prior_GramNeg),
    Pct_with_Prior_GramNeg = round(100 * sum(prior_positive_categories$N_with_Prior_GramNeg) / sum(prior_positive_categories$Total), 1),
    N_with_Prior_Polymicrobial = sum(prior_positive_categories$N_with_Prior_Polymicrobial),
    Pct_with_Prior_Polymicrobial = round(100 * sum(prior_positive_categories$N_with_Prior_Polymicrobial) / sum(prior_positive_categories$Total), 1),
    N_with_Prior_Other = sum(prior_positive_categories$N_with_Prior_Other),
    Pct_with_Prior_Other = round(100 * sum(prior_positive_categories$N_with_Prior_Other) / sum(prior_positive_categories$Total), 1),
    N_with_Prior_RecurContam = sum(prior_positive_categories$N_with_Prior_RecurContam),
    Pct_with_Prior_RecurContam = round(100 * sum(prior_positive_categories$N_with_Prior_RecurContam) / sum(prior_positive_categories$Total), 1)
  )
)

# --- SECTION 5: Culture Source Category Breakdown ---
culture_source_summary <- primary_blood_cultures %>%  group_by(CultureSourceCategory) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(Percentage = round(100 * N / sum(N), 1)) %>%
  arrange(desc(N)) %>%
  rename(Category = CultureSourceCategory)

# Add total row
culture_source_summary <- bind_rows(
  culture_source_summary,
  tibble(
    Category = "TOTAL",
    N = sum(culture_source_summary$N),
    Percentage = 100.0
  )
)

# --- SAVE TO EXCEL ---

# Format all summaries before writing to Excel
overall_summary_fmt <- overall_summary %>% format_excel_output()
cultivability_summary_fmt <- cultivability_summary %>% format_excel_output()
culture_source_summary_fmt <- culture_source_summary %>% format_excel_output()
organism_summary_fmt <- organism_summary %>% format_excel_output()
multi_organism_summary_fmt <- multi_organism_summary %>% format_excel_output()
prior_culture_distribution_fmt <- prior_culture_distribution %>% format_excel_output()
prior_positive_breakdown_fmt <- prior_positive_breakdown %>% format_excel_output()
prior_positive_categories_fmt <- prior_positive_categories %>% format_excel_output()

wb <- createWorkbook()
addWorksheet(wb, "Culture_Level_Summary")
addWorksheet(wb, "Cultivability_Breakdown")
addWorksheet(wb, "Culture_Source_Breakdown")
addWorksheet(wb, "Organism_Type_Breakdown")
addWorksheet(wb, "Multi_Organism_Analysis")
addWorksheet(wb, "Prior_Culture_Distribution")
addWorksheet(wb, "Prior_Positive_Breakdown")
addWorksheet(wb, "Prior_Positive_Categories")

writeData(wb, "Culture_Level_Summary", overall_summary_fmt)
writeData(wb, "Cultivability_Breakdown", cultivability_summary_fmt)
writeData(wb, "Culture_Source_Breakdown", culture_source_summary_fmt)
writeData(wb, "Organism_Type_Breakdown", organism_summary_fmt)
writeData(wb, "Multi_Organism_Analysis", multi_organism_summary_fmt)
writeData(wb, "Prior_Culture_Distribution", prior_culture_distribution_fmt)
writeData(wb, "Prior_Positive_Breakdown", prior_positive_breakdown_fmt)
writeData(wb, "Prior_Positive_Categories", prior_positive_categories_fmt)

# Add formatting note
addWorksheet(wb, "Formatting_Note")
writeData(wb, "Formatting_Note", create_formatting_note())

saveWorkbook(wb, "Output Files/Excels/Pipeline/BC_1_Initial_Extract_excel.xlsx", overwrite = TRUE)