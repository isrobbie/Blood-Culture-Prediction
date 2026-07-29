library(data.table)
library(dplyr)
library(lubridate)
library(stringr)
library(openxlsx)
library(writexl)

# Load Excel output formatting helpers
source("01_Preprocessing/Excel_Output_Helper.r")

# ============================================================
# ADD OTHER MICROBIOLOGY RESULTS TO BLOOD CULTURES
# ============================================================
# Adds microbiology results from other sites (Urine, Resp, etc.)
# - Lookback period: 28 days prior to blood culture
# - Only includes POSITIVE cultures with valid bugs
# - Includes BACTERIA and FUNGI; EXCLUDES VIRUSES
# ============================================================

# Load blood culture data from BC_7
bc_df <- fread("Output Files/CSVs/blood_cultures_with_blood_tests.csv",
               colClasses = c(ClusterID = "character"))

# Extract keys for matching (include BloodCultureID for unique identification)
bc_keys <- bc_df[, .(BloodCultureID, ClusterID, AccessionNumber, CollectionDateTime, CultureCategory)]
bc_keys[, CollectionDateTime := as.POSIXct(CollectionDateTime, tz = "UTC")]
cohort_ids <- unique(bc_df$ClusterID)

# Load Microbiology Data
micro_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_Microbiology.csv",
                  select = c("ClusterID", "AccessionNumber", "CollectionDateTime", 
                             "BatTestName", "BugCode", "BugName", "SpecimenFull"))

# ============================================================
# LOAD BUG GROUP MAPPING & DEFINE CONSTANTS
# ============================================================

# Load Bug Group Mapping
bug_lookup_dt <- fread("00_Data/lookup_tables/bug_group_table.csv")
setkey(bug_lookup_dt, BugCode)

# Define Constants 
contaminant_only_groups <- c("VIRIDANS AND OTHER STREPS", "CoNS (CONTAMINANT)", "OTHER CONTAMINANT")
grampos_groups <- c("S AUREUS", "ENTEROCOCCUS", "S PNEUMONIAE", 
                    "BETA HAEMOLYTIC STREP", "OTHER PATHOGENIC STREPS", "OTHER GRAM-POSITIVE")
gramneg_groups <- c("E COLI", "KLEBSIELLA", "P AERUGINOSA", "ENTEROBACTER", "OTHER ENTEROBACTERALES",
                    "OTHER GRAM-NEGATIVE")
pathogen_groups <- c("S AUREUS", "E COLI", "KLEBSIELLA", "P AERUGINOSA", 
                     "ENTEROCOCCUS", "S PNEUMONIAE", "ENTEROBACTER", 
                     "BETA HAEMOLYTIC STREP", "OTHER ENTEROBACTERALES", 
                     "OTHER PATHOGENIC STREPS", "ANAEROBES", "CANDIDA",
                     "OTHER GRAM-POSITIVE", "OTHER GRAM-NEGATIVE")
# ============================================================
# PREPARE MICROBIOLOGY DATA
# ============================================================

# Filter for relevant cultures:
# 1. Not Blood Culture
# 2. Positive (has BugName)
# 3. Valid dates
# 4. Exclude viral PCR/testing (exclude viruses)
# 5. Exclude screening tests (colonization, not infection)
# 6. Exclude screening specimens (rectal swabs, nasal swabs)
viral_tests <- c("PCR VIROLOGY", "FLU & RSV A/B PCR", "RESPIRATORY PCR", 
                 "NOROVIRUS PCR", "MONKEYPOX PCR", "REF LAB NOROVIRUS",
                 "LOWER RESPIRATORY PCR", "REF RESPIRATORY PCR", 
                 "Ref Lab Resp PCR", "FAECES VIRAL NOT DONE")

screening_tests <- c("MRSA SCREEN", "ESBL SCREEN", "CPE SCREEN", 
                     "CANDIDA AURIS SCREEN", "GROUP B STREP SCREEN",
                     "VRE SCREEN", "CRE SCREEN")

other_micro <- micro_df[
  BatTestName != "BLOOD CULTURE" &
  !is.na(BugName) & BugName != "" & BugName != "NULL" &
  !is.na(CollectionDateTime) &
  !(BatTestName %in% viral_tests) &
  !(BatTestName %in% screening_tests) &
  !grepl("virus|VIRUS", BugName, ignore.case = TRUE) &
  !grepl("BLOOD.*FOR.*CULTURE|BLOOD.*CULTURE", SpecimenFull, ignore.case = TRUE) &
  !grepl("^NOSE SWAB$|^NASAL SWAB$|RECTAL SWAB|RECTAL/STOOL SWAB", SpecimenFull, ignore.case = TRUE)
]

# Pre-filter to cohort ClusterIDs before further processing
other_micro <- other_micro[ClusterID %in% cohort_ids]

# Ensure dates are POSIXct (fread already parses correctly, no ymd_hms needed)
other_micro[, CollectionDateTime := as.POSIXct(CollectionDateTime, tz = "UTC")]
other_micro <- other_micro[!is.na(CollectionDateTime)]

# Join with Bug Group
setkey(other_micro, BugCode)
other_micro[, BugGroup := bug_lookup_dt[other_micro, BugGroup, on = "BugCode"]]

# FILTER: BACTERIA and FUNGI 
other_micro <- other_micro[
  !is.na(BugGroup) &
  BugGroup != "" &
  BugGroup != "NULL" &
  BugGroup %in% c(pathogen_groups, contaminant_only_groups)
]

# Deduplicate: One row per organism per culture
other_micro <- unique(other_micro, by = c("ClusterID", "AccessionNumber", "BugName"))

# ============================================================
# CLASSIFY ORGANISMS (GRAM POS/NEG)
# ============================================================

other_micro[, `:=`(
  IsGramPos = BugGroup %in% grampos_groups,
  IsGramNeg = BugGroup %in% gramneg_groups
)]

# Filter to pathogenic bacteria only (exclude fungi and contaminants)
pathogenic_pre <- other_micro[
  BugGroup %in% setdiff(pathogen_groups, "CANDIDA") &
  !(BugGroup %in% contaminant_only_groups)
]

# ============================================================
# CATEGORIZE SOURCES
# ============================================================

pathogenic_pre[, SpecimenUpper := toupper(SpecimenFull)]

pathogenic_pre[, SourceCategory := fcase(
  # URINE
  grepl("URINE|CSU|MSU|CATHETER URINE", SpecimenUpper), "Urine",
  
  # RESPIRATORY (includes cough, nasopharyngeal, mouth swabs; excludes PLEURAL which is sterile site)
  grepl("SPUTUM|BRONCH|LAVAGE|TRACHEAL|RESPIRATORY|THROAT|NASAL|NOSE|SWAB.*VIRAL|COUGH|NASOPHARYNGEAL|NASO PHARYNGEAL|MOUTH SWAB", SpecimenUpper), "Respiratory",
  
  # EYE (check before Skin_SoftTissue to catch conjunctival)
  grepl("EYE|CORNEA|CORNEAL|CONJUNCTIVAL|EYELID", SpecimenUpper), "Eye",
  
  # GENITAL (check before Skin_SoftTissue to catch vulval, penile, cervical)
  grepl("VAGINAL|VAGINA|CERVIX|VULVAL|PENILE|PENIS", SpecimenUpper), "Genital",
  
  # SKIN / SOFT TISSUE / WOUND (includes SWAB with site descriptors, excludes generic SWAB)
  # Matches: WOUND, PUS, TISSUE, BIOPSY, ASPIRATE, or site-specific swabs (both "SWAB;SITE" and "SITE SWAB" formats)
  (grepl("WOUND|SKIN|ULCER|PUS|ABSCESS|TISSUE|BIOPSY|ASPIRATE", SpecimenUpper) | 
   grepl("(^SWAB|SWAB$)", SpecimenUpper) & grepl("FOOT|GROIN|BACK|LEG|ABDOMEN|ARM|HAND|CHEST|STERNUM|THIGH|SHOULDER|ELBOW|WRIST|ANKLE|TOE|FINGER|BUTTOCK|GLUTEAL|PERINEUM|PEG|SITE|AXILLA|UMBILICUS|BREAST", SpecimenUpper)) & 
    !grepl("THROAT|NASAL|NOSE|RECTAL|VAGINAL|URETHRAL|EYE|EAR|COUGH|NASOPHARYNGEAL|MOUTH|VULVAL|PENILE|PENIS|CERVIX|CONJUNCTIVAL|RESPIRATORY|BRONCH|SPUTUM|CORNEA", SpecimenUpper) &
    !(SpecimenUpper == "SWAB") &
    !grepl("^SWAB;OTHER|SWAB OTHER|SWAB:UNSPECIFIED|SWAB;UNSPECIFIED|SWAB;UNKNOWN|SWAB; Other", SpecimenUpper), "Skin_SoftTissue",
  
  # CSF (CEREBROSPINAL FLUID)
  grepl("CSF|CEREBROSPINAL", SpecimenUpper), "CSF",
  
  # FLUID (STERILE SITES - including PLEURAL, excluding CSF)
  grepl("ASCITIC|PERITONEAL|JOINT|SYNOVIAL|SYNOVIUM|BILE|CAPD|KIDNEY PERFUSION|PANCREATIC PERFUSION|SONICATION|BURSA|PERICARDIAL|PLEURAL", SpecimenUpper) |
  (grepl("FLUID", SpecimenUpper) & grepl("ASCITIC|PERITONEAL|JOINT|SYNOVIAL|BILE|CAPD|KIDNEY|PANCREATIC|SONICATION|BURSA|PERICARDIAL|DRAIN|GALL|ABDOMEN|ABDOMINAL|BREAST|SEROMA|BILIARY|PLEURAL", SpecimenUpper)) &
    !grepl("BRONCH|CSF|CEREBROSPINAL", SpecimenUpper), "SterileFluid",
  
  # STOOL / GI
  grepl("STOOL|FAECES|RECTAL|GASTRIC", SpecimenUpper), "Gastrointestinal",
  
  # DEVICE / LINE TIPS
  grepl("TIP|LINE|CATHETER", SpecimenUpper) & !grepl("URINE", SpecimenUpper), "Device_LineTip",
  
  # EAR
  grepl("EAR", SpecimenUpper), "Ear",
  
  # BONE
  grepl("BONE|ACETABULAR|FEMORAL|TIBIAL", SpecimenUpper), "Bone",
  
  # GENERIC/UNSPECIFIED SWAB (catch before default Other) - Reclassified to Skin_SoftTissue
  (SpecimenUpper == "SWAB") | 
  grepl("^SWAB;OTHER|SWAB OTHER|SWAB:UNSPECIFIED|SWAB;UNSPECIFIED|SWAB;UNKNOWN|SWAB; Other|SWAB;OTHER", SpecimenUpper), "Skin_SoftTissue",
  
  # GENERIC FLUID (without descriptors - catch before default Other) - Reclassified to SterileFluid
  (SpecimenUpper == "FLUID"), "SterileFluid",
  
  # DEFAULT
  default = "Other"
)]

# Clean up
pathogenic_pre[, SpecimenUpper := NULL]

# ============================================================
# MATCH TO BLOOD CULTURES (28-DAY LOOKBACK)
# ============================================================

# Define lookback window
bc_keys[, `:=`(
  window_start = CollectionDateTime - days(28),
  window_end = CollectionDateTime
)]

# Perform Non-Equi Join
matched_micro <- pathogenic_pre[bc_keys, 
                             on = .(ClusterID, 
                                    CollectionDateTime >= window_start, 
                                    CollectionDateTime <= window_end),
                             nomatch = 0L,
                             allow.cartesian = TRUE,
                             .(ClusterID, 
                               BloodCultureID = i.BloodCultureID,
                               BC_AccessionNumber = i.AccessionNumber, 
                               BC_Date = i.CollectionDateTime,
                               Micro_AccessionNumber = x.AccessionNumber,
                               Micro_Date = x.CollectionDateTime,
                               SourceCategory,
                               SpecimenFull,
                             BugName,
                             BugGroup,
                             IsGramPos,
                             IsGramNeg)]

# ============================================================
# AGGREGATE PER BLOOD CULTURE
# ============================================================

# Count number of positive cultures per source category
source_counts <- matched_micro[, .(
  # Source Categories
  Prior_Urine_Pos = as.integer(any(SourceCategory == "Urine")),
  Prior_Resp_Pos = as.integer(any(SourceCategory == "Respiratory")),
  Prior_SkinSoftTissue_Pos = as.integer(any(SourceCategory == "Skin_SoftTissue")),
  Prior_SterileFluid_Pos = as.integer(any(SourceCategory == "SterileFluid")),
  Prior_CSF_Pos = as.integer(any(SourceCategory == "CSF")),
  Prior_GI_Pos = as.integer(any(SourceCategory == "Gastrointestinal")),
  Prior_Device_Pos = as.integer(any(SourceCategory == "Device_LineTip")),
  Prior_Genital_Pos = as.integer(any(SourceCategory == "Genital")),
  Prior_Ear_Pos = as.integer(any(SourceCategory == "Ear")),
  Prior_Eye_Pos = as.integer(any(SourceCategory == "Eye")),
  Prior_Bone_Pos = as.integer(any(SourceCategory == "Bone")),
  Prior_Other_Pos = as.integer(any(SourceCategory == "Other")),
  
  # Organism Classifications
  Prior_GramPos_Pos = as.integer(any(IsGramPos)),
  Prior_GramNeg_Pos = as.integer(any(IsGramNeg)),
  
  # Total count of distinct positive cultures (any source)
  Prior_Total_Other_Cultures_Pos = uniqueN(Micro_AccessionNumber)
), by = BloodCultureID]

# 2. Create a flag for "Any Prior Positive Non-Blood Culture"
source_counts[, Any_Prior_Other_Pos := as.integer(Prior_Total_Other_Cultures_Pos > 0)]

# ============================================================
# JOIN BACK TO MAIN DATASET
# ============================================================

setkey(bc_df, BloodCultureID)
setkey(source_counts, BloodCultureID)

bc_df <- source_counts[bc_df, on = .(BloodCultureID)]

# Fill NAs with 0
cols_to_fill <- c("Prior_Urine_Pos", "Prior_Resp_Pos", "Prior_SkinSoftTissue_Pos", 
                  "Prior_SterileFluid_Pos", "Prior_CSF_Pos", "Prior_GI_Pos", "Prior_Device_Pos", 
                  "Prior_Genital_Pos", "Prior_Ear_Pos", "Prior_Eye_Pos", "Prior_Bone_Pos",
                  "Prior_Other_Pos", "Prior_Total_Other_Cultures_Pos", "Any_Prior_Other_Pos",
                  "Prior_GramPos_Pos", "Prior_GramNeg_Pos")

for (col in cols_to_fill) {
  set(bc_df, which(is.na(bc_df[[col]])), col, 0L)
}

# ============================================================
# ============================================================
# ADD COVID-19 TEST RESULT
# ============================================================
# Feature: CoV2_Result_14d
#   1  = most recent COVID PCR within 14 days prior = DETECTED
#   0  = most recent COVID PCR within 14 days prior = NOT DETECTED
#   NA = no test in window (pre-pandemic or not tested)
# ============================================================

covid <- fread(
  "Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_COVIDTests/IORD_BloodCultPatients_57_20251030_COVIDTests.csv",
  select = c("ClusterID", "CollectionDateTime", "Result"),
  colClasses = c(ClusterID = "character")
)
covid[, TestDate := as.Date(as.POSIXct(CollectionDateTime, tz = "UTC"))]
covid <- covid[Result %in% c("DET", "NDET")]
covid[, CoV2_Result_14d := fifelse(Result == "DET", 1L, 0L)]
# One row per (patient, day) - DET beats NDET; required for clean rolling join
covid <- covid[order(ClusterID, TestDate, -CoV2_Result_14d),
               .(CoV2_Result_14d = CoV2_Result_14d[1L]),
               by = .(ClusterID, TestDate)]
setkey(covid, ClusterID, TestDate)

bc_df[, CultureDate := as.Date(as.POSIXct(CollectionDateTime, tz = "UTC"))]
covid_joined <- covid[bc_df[, .(BloodCultureID, ClusterID, CultureDate)],
                      on = .(ClusterID, TestDate = CultureDate),
                      roll = 14, rollends = c(TRUE, FALSE)]
bc_df <- covid_joined[, .(BloodCultureID, CoV2_Result_14d)][bc_df, on = "BloodCultureID"]
bc_df[, CultureDate := NULL]

n_covered <- bc_df[!is.na(CoV2_Result_14d), .N]
message(sprintf("CoV2 coverage: %s / %s cultures (%.1f%%)",
                format(n_covered, big.mark = ","),
                format(nrow(bc_df), big.mark = ","),
                100 * n_covered / nrow(bc_df)))

# ============================================================
# SAVE OUTPUT
# ============================================================

fwrite(bc_df, "Output Files/CSVs/blood_cultures_with_other_micro.csv")

# ============================================================
# CREATE EXCEL SUMMARY
# ============================================================

# Summary Table
summary_stats <- bc_df[, .(
  Category = "All Cultures",
  N = .N,
  
  Any_Prior_Other_Pos_N = sum(Any_Prior_Other_Pos),
  Any_Prior_Other_Pos_Pct = round(100 * mean(Any_Prior_Other_Pos), 1),
  
  # Sources
  Urine_Pos_N = sum(Prior_Urine_Pos),
  Urine_Pos_Pct = round(100 * mean(Prior_Urine_Pos), 1),
  
  Resp_Pos_N = sum(Prior_Resp_Pos),
  Resp_Pos_Pct = round(100 * mean(Prior_Resp_Pos), 1),
  
  Skin_Pos_N = sum(Prior_SkinSoftTissue_Pos),
  Skin_Pos_Pct = round(100 * mean(Prior_SkinSoftTissue_Pos), 1),
  
  Fluid_Pos_N = sum(Prior_SterileFluid_Pos),
  Fluid_Pos_Pct = round(100 * mean(Prior_SterileFluid_Pos), 1),
  
  CSF_Pos_N = sum(Prior_CSF_Pos),
  CSF_Pos_Pct = round(100 * mean(Prior_CSF_Pos), 1),
  
  GI_Pos_N = sum(Prior_GI_Pos),
  GI_Pos_Pct = round(100 * mean(Prior_GI_Pos), 1),
  
  Device_Pos_N = sum(Prior_Device_Pos),
  Device_Pos_Pct = round(100 * mean(Prior_Device_Pos), 1),
  
  Genital_Pos_N = sum(Prior_Genital_Pos),
  Genital_Pos_Pct = round(100 * mean(Prior_Genital_Pos), 1),
  
  Ear_Pos_N = sum(Prior_Ear_Pos),
  Ear_Pos_Pct = round(100 * mean(Prior_Ear_Pos), 1),
  
  Eye_Pos_N = sum(Prior_Eye_Pos),
  Eye_Pos_Pct = round(100 * mean(Prior_Eye_Pos), 1),
  
  Bone_Pos_N = sum(Prior_Bone_Pos),
  Bone_Pos_Pct = round(100 * mean(Prior_Bone_Pos), 1),
  
  Other_Pos_N = sum(Prior_Other_Pos),
  Other_Pos_Pct = round(100 * mean(Prior_Other_Pos), 1),
  
  # Organism Types
  GramPos_Pos_N = sum(Prior_GramPos_Pos),
  GramPos_Pos_Pct = round(100 * mean(Prior_GramPos_Pos), 1),
  
  GramNeg_Pos_N = sum(Prior_GramNeg_Pos),
  GramNeg_Pos_Pct = round(100 * mean(Prior_GramNeg_Pos), 1)
)]

# By Culture Result
summary_by_result <- bc_df[!is.na(CultureCategory), .(
  Category = first(CultureCategory),
  N = .N,
  
  Any_Prior_Other_Pos_N = sum(Any_Prior_Other_Pos),
  Any_Prior_Other_Pos_Pct = round(100 * mean(Any_Prior_Other_Pos), 1),
  
  # Sources
  Urine_Pos_N = sum(Prior_Urine_Pos),
  Urine_Pos_Pct = round(100 * mean(Prior_Urine_Pos), 1),
  
  Resp_Pos_N = sum(Prior_Resp_Pos),
  Resp_Pos_Pct = round(100 * mean(Prior_Resp_Pos), 1),
  
  Skin_Pos_N = sum(Prior_SkinSoftTissue_Pos),
  Skin_Pos_Pct = round(100 * mean(Prior_SkinSoftTissue_Pos), 1),
  
  Fluid_Pos_N = sum(Prior_SterileFluid_Pos),
  Fluid_Pos_Pct = round(100 * mean(Prior_SterileFluid_Pos), 1),
  
  CSF_Pos_N = sum(Prior_CSF_Pos),
  CSF_Pos_Pct = round(100 * mean(Prior_CSF_Pos), 1),
  
  GI_Pos_N = sum(Prior_GI_Pos),
  GI_Pos_Pct = round(100 * mean(Prior_GI_Pos), 1),
  
  Device_Pos_N = sum(Prior_Device_Pos),
  Device_Pos_Pct = round(100 * mean(Prior_Device_Pos), 1),
  
  Genital_Pos_N = sum(Prior_Genital_Pos),
  Genital_Pos_Pct = round(100 * mean(Prior_Genital_Pos), 1),
  
  Ear_Pos_N = sum(Prior_Ear_Pos),
  Ear_Pos_Pct = round(100 * mean(Prior_Ear_Pos), 1),
  
  Eye_Pos_N = sum(Prior_Eye_Pos),
  Eye_Pos_Pct = round(100 * mean(Prior_Eye_Pos), 1),
  
  Bone_Pos_N = sum(Prior_Bone_Pos),
  Bone_Pos_Pct = round(100 * mean(Prior_Bone_Pos), 1),
  
  Other_Pos_N = sum(Prior_Other_Pos),
  Other_Pos_Pct = round(100 * mean(Prior_Other_Pos), 1),
  
  # Organism Types
  GramPos_Pos_N = sum(Prior_GramPos_Pos),
  GramPos_Pos_Pct = round(100 * mean(Prior_GramPos_Pos), 1),
  
  GramNeg_Pos_N = sum(Prior_GramNeg_Pos),
  GramNeg_Pos_Pct = round(100 * mean(Prior_GramNeg_Pos), 1)
), by = CultureCategory][order(-Category)]

final_summary <- rbind(summary_stats, summary_by_result, fill = TRUE)
final_summary_formatted <- as.data.frame(final_summary) %>% format_excel_output()

# ============================================================
# COMMON ORGANISMS BY SOURCE ANALYSIS
# ============================================================

# Calculate total counts per source for percentage calculation
source_totals <- matched_micro[, .(Total_N = .N), by = SourceCategory]

# Count organisms per source
organisms_by_source <- matched_micro[, .(N = .N), by = .(SourceCategory, BugName)]

# Join with totals to calculate percentage
organisms_by_source <- organisms_by_source[source_totals, on = "SourceCategory"]
organisms_by_source[, Percentage := round(100 * N / Total_N, 1)]

# Sort by Source and N (descending)
setorder(organisms_by_source, SourceCategory, -N)

# Filter to top 10 per source
top_organisms_by_source <- organisms_by_source[, head(.SD, 10), by = SourceCategory]

# Format for Excel
top_organisms_formatted <- top_organisms_by_source[, .(
  Source = SourceCategory,
  Organism = BugName,
  Count = N,
  Percentage = Percentage
)] %>% format_excel_output()

# ============================================================
# BREAKDOWN OF "OTHER" CATEGORY
# ============================================================

# Filter to "Other" category only (using pathogenic bacteria only)
other_category_micro <- matched_micro[SourceCategory == "Other"]

# Breakdown by SpecimenFull (actual sample types)
other_by_specimen <- other_category_micro[, .(
  N = .N,
  Unique_Cultures = uniqueN(Micro_AccessionNumber),
  Unique_BCs = uniqueN(BC_AccessionNumber)
), by = SpecimenFull]
setorder(other_by_specimen, -N)

# Breakdown by Organism
other_by_organism <- other_category_micro[, .(
  N = .N,
  Unique_Cultures = uniqueN(Micro_AccessionNumber),
  Unique_BCs = uniqueN(BC_AccessionNumber)
), by = .(BugName, BugGroup)]
setorder(other_by_organism, -N)

# Breakdown by SpecimenFull and Organism (top combinations)
other_by_specimen_organism <- other_category_micro[, .(
  N = .N,
  Unique_Cultures = uniqueN(Micro_AccessionNumber)
), by = .(SpecimenFull, BugName, BugGroup)]
setorder(other_by_specimen_organism, -N)

# Format for Excel
other_by_specimen_formatted <- other_by_specimen[, .(
  Sample_Type = SpecimenFull,
  Count = N,
  Unique_Cultures = Unique_Cultures,
  Unique_BCs = Unique_BCs
)] %>% format_excel_output()

other_by_organism_formatted <- other_by_organism[, .(
  Organism = BugName,
  Bug_Group = BugGroup,
  Count = N,
  Unique_Cultures = Unique_Cultures,
  Unique_BCs = Unique_BCs
)] %>% format_excel_output()

other_by_specimen_organism_formatted <- other_by_specimen_organism[, .(
  Sample_Type = SpecimenFull,
  Organism = BugName,
  Bug_Group = BugGroup,
  Count = N,
  Unique_Cultures = Unique_Cultures
)] %>% format_excel_output()

# Info Table
info_table <- data.frame(
  Item = c("Analysis", "Total blood cultures", "Lookback period", "Included Micro Results", "Filter", "Date"),
  Value = c("Prior Non-Blood Microbiology Results",
            format_count(nrow(bc_df)),
            "28 days",
            "Positive cultures from non-blood sites (Urine, Resp, etc.)",
            "Bacteria and Fungi (Defined pathogen/contaminant groups including OTHER GRAM-POSITIVE and OTHER GRAM-NEGATIVE; excludes CANDIDA, viruses, contaminants, and NULL)",
            as.character(Sys.Date()))
)

# Create Workbook
wb <- createWorkbook()

addWorksheet(wb, "Summary")
writeData(wb, "Summary", 
          data.frame(Title = "Prior Positive Non-Blood Cultures (28 Days)"), 
          startRow = 1, colNames = FALSE)
writeData(wb, "Summary", final_summary_formatted, startRow = 3)

addWorksheet(wb, "Top_Organisms_by_Source")
writeData(wb, "Top_Organisms_by_Source", 
          data.frame(Title = "Top 10 Organisms per Source Category"), 
          startRow = 1, colNames = FALSE)
writeData(wb, "Top_Organisms_by_Source", top_organisms_formatted, startRow = 3)

addWorksheet(wb, "Other_Category_by_Sample")
writeData(wb, "Other_Category_by_Sample", 
          data.frame(Title = "Other Category: Breakdown by Sample Type"), 
          startRow = 1, colNames = FALSE)
writeData(wb, "Other_Category_by_Sample", other_by_specimen_formatted, startRow = 3)

addWorksheet(wb, "Other_Category_by_Organism")
writeData(wb, "Other_Category_by_Organism", 
          data.frame(Title = "Other Category: Breakdown by Organism"), 
          startRow = 1, colNames = FALSE)
writeData(wb, "Other_Category_by_Organism", other_by_organism_formatted, startRow = 3)

addWorksheet(wb, "Other_Category_Combinations")
writeData(wb, "Other_Category_Combinations", 
          data.frame(Title = "Other Category: Sample Type x Organism Combinations"), 
          startRow = 1, colNames = FALSE)
writeData(wb, "Other_Category_Combinations", other_by_specimen_organism_formatted, startRow = 3)

addWorksheet(wb, "Info")
writeData(wb, "Info", info_table)

addWorksheet(wb, "Formatting_Note")
writeData(wb, "Formatting_Note", create_formatting_note())

saveWorkbook(wb, "Output Files/Excels/Pipeline/BC_8_Add_Other_Micro_excel.xlsx", overwrite = TRUE)

