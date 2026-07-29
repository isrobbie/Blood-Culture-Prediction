library(data.table)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(purrr)
library(writexl)

# Load Excel output formatting helpers
source("01_Preprocessing/Excel_Output_Helper.r")

# ============================================================
# ADD DEMOGRAPHICS TO BLOOD CULTURE DATA
# ============================================================
# This script adds patient demographics (age, sex, ethnicity,
# death dates) to the primary blood cultures dataset
# ============================================================

# Read primary blood cultures 
bc_df <- fread("Output Files/CSVs/primary_blood_cultures.csv")
initial_bc_count <- nrow(bc_df)

# Extract cohort ClusterIDs for filtering
cohort_cluster_ids <- unique(bc_df$ClusterID)

# Read and process demographics
demographics <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_CohortDemographics.csv")
demographics <- demographics %>%
  mutate(
    # Ensure ClusterID is character for consistent joining
    ClusterID = as.character(ClusterID),
    
    # Replace "NULL" strings with actual NA
    LinkedBirthMonth = na_if(LinkedBirthMonth, "NULL"),
    LinkedDeathdate = na_if(LinkedDeathdate, "NULL"),
    
    # Convert dates
    LinkedBirthMonth = ymd_hms(LinkedBirthMonth),
    LinkedDeathdate = ymd_hms(LinkedDeathdate),
    
    # Convert sex to factor
    LinkedSex = factor(LinkedSex, levels = c("M", "F"))
  )

# Read ED and Inpatient data for social demographics, filtered to cohort
ed_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_EDAttendances.csv")
ed_df <- ed_df[ClusterID %in% cohort_cluster_ids]

inpat_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_InpatEpisodes.csv")
inpat_df <- inpat_df[ClusterID %in% cohort_cluster_ids]

# ============================================================
# 1. PREPARE DATA
# ============================================================

# Ensure ClusterID is character and convert CollectionDateTime
bc_df[, `:=`(ClusterID = as.character(ClusterID),
             CollectionDateTime = ymd_hms(CollectionDateTime))]

# Prepare ED and Inpatient data
ed_df[, ClusterID := as.character(ClusterID)]
inpat_df[, ClusterID := as.character(ClusterID)]

# Convert dates
ed_df[, c("ArrivalDateTime", "DepartureDateTime") := 
       lapply(.SD, function(x) {
         if (is.character(x)) ymd_hms(fifelse(x == "NULL", NA_character_, x))
         else x
       }), 
      .SDcols = c("ArrivalDateTime", "DepartureDateTime")]

inpat_df[, c("AdmissionDate", "DischargeDate") := 
          lapply(.SD, function(x) {
            if (is.character(x)) ymd_hms(fifelse(x == "NULL", NA_character_, x))
            else x
          }), 
         .SDcols = c("AdmissionDate", "DischargeDate")]

# ============================================================
# 2. JOIN DEMOGRAPHICS
# ============================================================

# Use data.table merge for this large join (faster)
bc_df <- merge(bc_df, demographics, by = "ClusterID", all.x = TRUE)

# Add prefixes to demographic columns
setnames(bc_df, 
         c("LinkedSex", "LinkedBirthMonth", "LinkedDeathdate", "SpineCheckDate", "LinkageWarningFlag"),
         c("DEMOGRAPHICS_LinkedSex", "DEMOGRAPHICS_LinkedBirthMonth", "DEMOGRAPHICS_LinkedDeathdate", 
           "DEMOGRAPHICS_SpineCheckDate", "DEMOGRAPHICS_LinkageWarningFlag"))

# ============================================================
# 3. CALCULATE AGE AND MORTALITY METRICS
# ============================================================

bc_df <- bc_df %>%
  mutate(
    AgeAtCulture = floor(as.numeric(difftime(CollectionDateTime, DEMOGRAPHICS_LinkedBirthMonth, units = "days")) / 365.25),
    DaysToDeath = as.numeric(difftime(DEMOGRAPHICS_LinkedDeathdate, CollectionDateTime, units = "days")),
    Mort_7 = as.integer(!is.na(DaysToDeath) & DaysToDeath >= 0 & DaysToDeath <= 7),
    Mort_28 = as.integer(!is.na(DaysToDeath) & DaysToDeath >= 0 & DaysToDeath <= 28),
    Mort_365 = as.integer(!is.na(DaysToDeath) & DaysToDeath >= 0 & DaysToDeath <= 365)
  )

# ============================================================
# 4. ADD SOCIAL DEMOGRAPHICS (ETHNICITY, IMD)
# ============================================================

# Match to ED attendances during culture
bc_ed_demographics <- bc_df %>%
  select(ClusterID, CollectionDateTime) %>%
  inner_join(
    ed_df %>% select(ClusterID, ArrivalDateTime, DepartureDateTime, 
                     EthnicGroupCode, IMDScore, IMDPercentile),
    by = "ClusterID", 
    relationship = "many-to-many"
  ) %>%
  filter(CollectionDateTime >= ArrivalDateTime & CollectionDateTime <= DepartureDateTime) %>%
  group_by(ClusterID, CollectionDateTime) %>%
  arrange(desc(ArrivalDateTime)) %>%
  slice(1) %>%
  ungroup() %>%
  select(ClusterID, CollectionDateTime,
         ED_EthnicGroupCode = EthnicGroupCode, ED_IMDScore = IMDScore, ED_IMDPercentile = IMDPercentile)

# Match to Inpatient episodes
bc_inpat_demographics <- bc_df %>%
  select(ClusterID, CollectionDateTime) %>%
  left_join(
    inpat_df %>% select(ClusterID, AdmissionDate, DischargeDate, 
                        EthnicGroupCode, IMDScore, IMDPercentile),
    by = "ClusterID",
    relationship = "many-to-many"
  ) %>%
  filter(!is.na(AdmissionDate) & 
         CollectionDateTime >= (AdmissionDate - hours(48)) & 
         CollectionDateTime <= DischargeDate) %>%
  group_by(ClusterID, CollectionDateTime) %>%
  arrange(desc(AdmissionDate)) %>%
  slice(1) %>%
  ungroup() %>%
  select(ClusterID, CollectionDateTime,
         INPATIENT_EthnicGroupCode = EthnicGroupCode, 
         INPATIENT_IMDScore = IMDScore, 
         INPATIENT_IMDPercentile = IMDPercentile)

# Merge social demographics
bc_df <- bc_df %>%
  left_join(bc_ed_demographics, by = c("ClusterID", "CollectionDateTime")) %>%
  left_join(bc_inpat_demographics, by = c("ClusterID", "CollectionDateTime")) %>%
  mutate(
    # Combine ED and Inpatient demographics, prioritizing Inpatient data
    SOCIAL_EthnicGroupCode = str_trim(coalesce(INPATIENT_EthnicGroupCode, ED_EthnicGroupCode)),
    SOCIAL_IMDScore = as.numeric(coalesce(INPATIENT_IMDScore, ED_IMDScore)),
    SOCIAL_IMDPercentile = as.integer(coalesce(INPATIENT_IMDPercentile, ED_IMDPercentile)),
    # Add join indicator flags
    HAS_ED_DEMOGRAPHICS = !is.na(ED_EthnicGroupCode) | !is.na(ED_IMDScore),
    HAS_INPATIENT_DEMOGRAPHICS = !is.na(INPATIENT_EthnicGroupCode) | !is.na(INPATIENT_IMDScore)
  ) %>%
  select(-c(INPATIENT_EthnicGroupCode, ED_EthnicGroupCode, 
            INPATIENT_IMDScore, ED_IMDScore,
            INPATIENT_IMDPercentile, ED_IMDPercentile))

# ============================================================
# FILL MISSING IMD FROM ANY HISTORICAL VISIT
# ============================================================
# For patients missing IMD from matched ED/Inpatient visit,
# use any available IMD value from any of their visits

# Get best available IMD per patient from all ED visits
ed_imd_any <- ed_df[!is.na(IMDPercentile), .(
  Fallback_IMDScore = first(IMDScore),
  Fallback_IMDPercentile = first(IMDPercentile)
), by = .(ClusterID = as.character(ClusterID))]

# Get best available IMD per patient from all Inpatient visits  
ip_imd_any <- inpat_df[!is.na(IMDPercentile), .(
  Fallback_IMDScore = first(IMDScore),
  Fallback_IMDPercentile = first(IMDPercentile)
), by = .(ClusterID = as.character(ClusterID))]

# Combine (prefer inpatient)
all_imd <- rbindlist(list(ip_imd_any, ed_imd_any))
all_imd <- all_imd[!duplicated(ClusterID)]

# Convert bc_df to data.table for efficient join
setDT(bc_df)

# Join fallback IMD
bc_df <- merge(bc_df, all_imd, by = "ClusterID", all.x = TRUE)

# Fill missing IMD with fallback values
bc_df[is.na(SOCIAL_IMDPercentile) & !is.na(Fallback_IMDPercentile), 
      `:=`(SOCIAL_IMDScore = Fallback_IMDScore,
           SOCIAL_IMDPercentile = Fallback_IMDPercentile)]

# Clean up fallback columns
bc_df[, c("Fallback_IMDScore", "Fallback_IMDPercentile") := NULL]

# Convert back to data.table
setDT(bc_df)

# ============================================================
# FILL MISSING ETHNICITY FROM ANY HISTORICAL VISIT
# ============================================================
# For patients with Unknown/missing ethnicity from matched visit,
# use any available non-Unknown ethnicity from any of their visits

# Get best available ethnicity per patient from all ED visits
# Exclude Unknown codes (Z, 99, NULL, NA) - trim whitespace first
ed_ethnicity_any <- ed_df[
  !is.na(EthnicGroupCode) & 
  !(trimws(EthnicGroupCode) %in% c("Z", "99", "NULL", "")),
  .(Fallback_EthnicGroupCode = first(trimws(EthnicGroupCode))),
  by = .(ClusterID = as.character(ClusterID))
]

# Get best available ethnicity per patient from all Inpatient visits
ip_ethnicity_any <- inpat_df[
  !is.na(EthnicGroupCode) & 
  !(trimws(EthnicGroupCode) %in% c("Z", "99", "NULL", "")),
  .(Fallback_EthnicGroupCode = first(trimws(EthnicGroupCode))),
  by = .(ClusterID = as.character(ClusterID))
]

# Combine (prefer inpatient)
all_ethnicity <- rbindlist(list(ip_ethnicity_any, ed_ethnicity_any))
all_ethnicity <- all_ethnicity[!duplicated(ClusterID)]

# Join fallback ethnicity
bc_df <- merge(bc_df, all_ethnicity, by = "ClusterID", all.x = TRUE)

# Fill missing/Unknown ethnicity with fallback values
bc_df[
  (is.na(SOCIAL_EthnicGroupCode) | 
   SOCIAL_EthnicGroupCode %in% c("Z", "99", "NULL") |
   SOCIAL_EthnicGroupCode == "") & 
  !is.na(Fallback_EthnicGroupCode),
  SOCIAL_EthnicGroupCode := Fallback_EthnicGroupCode
]

# Clean up fallback column
bc_df[, Fallback_EthnicGroupCode := NULL]

# Group ethnicity into granular categories
bc_df[, SOCIAL_Ethnicity := fcase(
  SOCIAL_EthnicGroupCode %in% c("A", "B", "C"), "White",
  SOCIAL_EthnicGroupCode %in% c("D", "E", "F", "G"), "Mixed",
  SOCIAL_EthnicGroupCode %in% c("H", "J", "K", "L"), "Asian",
  SOCIAL_EthnicGroupCode %in% c("M", "N", "P"), "Black",
  SOCIAL_EthnicGroupCode %in% c("R", "S"), "Other",
  SOCIAL_EthnicGroupCode %in% c("Z", "99", "NULL") | is.na(SOCIAL_EthnicGroupCode), "Unknown",
  default = "Unknown"
)]

# Convert to factor
bc_df[, SOCIAL_Ethnicity := factor(SOCIAL_Ethnicity, 
                                    levels = c("White", "Mixed", "Asian", "Black", "Other", "Unknown"))]

# ============================================================
# 5. SAVE OUTPUT
# ============================================================

fwrite(bc_df, "Output Files/CSVs/blood_cultures_with_demographics.csv")

# ============================================================
# 6. CREATE EXCEL SUMMARY STATISTICS
# ============================================================

# Define culture groups to analyze (using CultureCategory from BC_1 as gospel truth)
culture_groups <- list(
  "All Cultures" = bc_df,
  "Negative (and single contam)" = bc_df %>% filter(CultureCategory == "Negative (and SC)"),
  "Positive Cultures" = bc_df %>% filter(CultureCategory == "Positive")
)

# Helper function to calculate statistics with counts and percentages
# Returns named list with separate N and Pct
calc_stat <- function(condition) {
  list(
    N = sum(condition, na.rm = TRUE),
    Pct = 100 * mean(condition, na.rm = TRUE)
  )
}

# Function to create demographic summary for a group
create_demographic_summary <- function(data, group_name) {
  n_total <- nrow(data)
  
  # Sex distribution
  n_male <- sum(data$DEMOGRAPHICS_LinkedSex == "M", na.rm = TRUE)
  n_female <- sum(data$DEMOGRAPHICS_LinkedSex == "F", na.rm = TRUE)
  n_sex_missing <- sum(is.na(data$DEMOGRAPHICS_LinkedSex))
  
  # Mortality statistics
  mort_7_stats <- calc_stat(data$Mort_7 == 1)
  mort_28_stats <- calc_stat(data$Mort_28 == 1)
  mort_365_stats <- calc_stat(data$Mort_365 == 1)
  
  tibble(
    Group = group_name,
    N_Cultures = n_total,
    N_Cultures_Pct = 100 * n_total / nrow(bc_df),
    
    # Age statistics
    Median_Age = median(data$AgeAtCulture, na.rm = TRUE),
    Q1_Age = quantile(data$AgeAtCulture, 0.25, na.rm = TRUE),
    Q3_Age = quantile(data$AgeAtCulture, 0.75, na.rm = TRUE),
    Mean_Age = mean(data$AgeAtCulture, na.rm = TRUE),
    SD_Age = sd(data$AgeAtCulture, na.rm = TRUE),
    Missing_Age_N = sum(is.na(data$AgeAtCulture)),
    
    # Sex distribution
    Male_N = n_male,
    Male_Pct = 100 * n_male / n_total,
    Female_N = n_female,
    Female_Pct = 100 * n_female / n_total,
    Sex_Missing_N = n_sex_missing,
    Sex_Missing_Pct = 100 * n_sex_missing / n_total,
    
    # Mortality
    Mort_7d_N = mort_7_stats$N,
    Mort_7d_Pct = mort_7_stats$Pct,
    Mort_28d_N = mort_28_stats$N,
    Mort_28d_Pct = mort_28_stats$Pct,
    Mort_365d_N = mort_365_stats$N,
    Mort_365d_Pct = mort_365_stats$Pct
  )
}

# Create summary table
demographics_summary <- bind_rows(
  imap(culture_groups, ~create_demographic_summary(.x, .y))
)

# Format for display using helper functions
demographics_summary_formatted <- demographics_summary %>%
  format_excel_output() %>%
  mutate(
    `Age: Median (Q1-Q3)` = paste0(Median_Age, " (", Q1_Age, "-", Q3_Age, ")"),
    `Age: Mean ± SD` = paste0(Mean_Age, " ± ", SD_Age)
  ) %>%
  select(Group, N_Cultures, N_Cultures_Pct, `Age: Median (Q1-Q3)`, `Age: Mean ± SD`, Missing_Age_N,
         Male_N, Male_Pct, Female_N, Female_Pct, Sex_Missing_N, Sex_Missing_Pct,
         Mort_7d_N, Mort_7d_Pct, Mort_28d_N, Mort_28d_Pct, Mort_365d_N, Mort_365d_Pct)

# Create vertical demographics table (variables as rows, groups as columns)
demographics_vertical <- demographics_summary %>%
  mutate(
    Age_Median_Q1_Q3 = paste0(format_number(Median_Age), " (", format_number(Q1_Age), "-", format_number(Q3_Age), ")"),
    Age_Mean_SD = paste0(format_number(Mean_Age), " ± ", format_number(SD_Age)),
    Male = paste0(format_count(Male_N), " (", format_percent(Male_Pct), ")"),
    Female = paste0(format_count(Female_N), " (", format_percent(Female_Pct), ")"),
    Sex_Missing = paste0(format_count(Sex_Missing_N), " (", format_percent(Sex_Missing_Pct), ")"),
    Mort_7d = paste0(format_count(Mort_7d_N), " (", format_percent(Mort_7d_Pct), ")"),
    Mort_28d = paste0(format_count(Mort_28d_N), " (", format_percent(Mort_28d_Pct), ")"),
    Mort_365d = paste0(format_count(Mort_365d_N), " (", format_percent(Mort_365d_Pct), ")"),
    N_Total = paste0(format_count(N_Cultures), " (", format_percent(N_Cultures_Pct), ")")
  ) %>%
  select(Group, N_Total, Age_Median_Q1_Q3, Age_Mean_SD, Male, Female, Sex_Missing, 
         Mort_7d, Mort_28d, Mort_365d) %>%
  tidyr::pivot_longer(cols = -Group, names_to = "Variable", values_to = "Value") %>%
  tidyr::pivot_wider(names_from = Group, values_from = Value) %>%
  mutate(Variable = factor(Variable, 
                           levels = c("N_Total", "Age_Median_Q1_Q3", "Age_Mean_SD", 
                                    "Male", "Female", "Sex_Missing",
                                    "Mort_7d", "Mort_28d", "Mort_365d"),
                           labels = c("N (%)", "Age: Median (Q1-Q3)", "Age: Mean ± SD",
                                    "Male, N (%)", "Female, N (%)", "Sex Missing, N (%)",
                                    "7-day Mortality, N (%)", "28-day Mortality, N (%)", "365-day Mortality, N (%)"))) %>%
  arrange(Variable)

# Create age distribution table
age_bins <- c(0, 18, 40, 60, 80, Inf)
age_labels <- c("0-17", "18-39", "40-59", "60-79", "80+")

create_age_distribution <- function(data, group_name) {
  data %>%
    mutate(AgeGroup = cut(AgeAtCulture, breaks = age_bins, labels = age_labels, right = FALSE)) %>%
    group_by(AgeGroup) %>%
    summarise(N = n(), .groups = "drop") %>%
    mutate(
      Group = group_name,
      Pct = 100 * N / sum(N)
    ) %>%
    select(Group, AgeGroup, N, Pct)
}

age_distribution <- bind_rows(
  imap(culture_groups, ~create_age_distribution(.x, .y))
) %>%
  pivot_wider(names_from = Group, values_from = c(N, Pct), names_sep = ": ") %>%
  format_excel_output()

# Create ethnicity distribution table
create_ethnicity_distribution <- function(data, group_name) {
  data %>%
    group_by(SOCIAL_Ethnicity) %>%
    summarise(N = n(), .groups = "drop") %>%
    mutate(
      Group = group_name,
      Pct = 100 * N / sum(N)
    ) %>%
    select(Group, SOCIAL_Ethnicity, N, Pct)
}

ethnicity_distribution <- bind_rows(
  imap(culture_groups, ~create_ethnicity_distribution(.x, .y))
) %>%
  pivot_wider(names_from = Group, values_from = c(N, Pct), names_sep = ": ") %>%
  format_excel_output()

# Save to Excel
write_xlsx(
  list(
    "Vertical_Summary" = demographics_vertical,
    "Summary_Formatted" = demographics_summary_formatted,
    "Summary_Raw" = demographics_summary %>% format_excel_output(),
    "Age_Distribution" = age_distribution,
    "Ethnicity_Distribution" = ethnicity_distribution,
    "Formatting_Note" = create_formatting_note()
  ),
  "Output Files/Excels/Pipeline/BC_2_Add_Demographics_excel.xlsx"
)

