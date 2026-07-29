library(data.table)
library(dplyr)
library(lubridate)
library(openxlsx)
library(writexl)

# Load Excel output formatting helpers
source("01_Preprocessing/Excel_Output_Helper.r")

# ============================================================
# ADD BLOOD TEST RESULTS TO BLOOD CULTURES
# ============================================================
# This script adds blood test information to blood cultures:
# - Matches blood tests within 24-hour time windows (±24 hours)
# - Finds closest values for key parameters
# - Calculates min/max/mean values within window
# ============================================================

# Read blood culture data with blood gas
bc_df <- fread("Output Files/CSVs/blood_cultures_with_blood_gas.csv",
               colClasses = c(ClusterID = "character"))
cohort_cluster_ids <- unique(bc_df$ClusterID)

# ============================================================
# 1. LOAD AND PREPARE BLOOD TEST DATA (OPTIMIZED - FILTER EARLY)
# ============================================================

# Define relevant tests upfront for efficient filtering
# Expanded to include tests from Lancet Digital Health paper
relevant_tests <- c(
  # Existing tests
  "CRP", "WHITE CELLS", "NEUTROPHILS", "LYMPHOCYTES",
  "CREATININE", "UREA", "BILIRUBIN", "ALT", "ALBUMIN",
  "HAEMOGLOBIN", "PLATELETS", "INR", "APTT",
  "SODIUM", "POTASSIUM",
  # ALL WCC components
  "EOSINOPHILS",           # 4.3M rows
  "MONOCYTES",             # Full blood count component
  "BASOPHILS",             # Full blood count component
  # Additional tests from paper (exact names from IORD data)
  "MEAN CELL VOL.",        # 4.4M rows (note the period!)
  "HAEMATOCRIT",           # 4.4M rows
  "ALK.PHOSPHATASE",       # 3.6M rows (note the period!)
  "GLUCOSE",               # 374K rows
  "Fibrinogen",            # 33K rows (capitalized!)
  "D-Dimer"                # 122K rows (capitalized!)
  # Note: FERRITIN and TOTAL PROTEIN not available in IORD data
)

# Load only relevant tests 
bloodtests <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_BloodTests.csv",
                    select = c("ClusterID", "CollectionDateTime", "TestName", "Value"),
                    colClasses = c(ClusterID = "character"))

# Filter to cohort AND relevant tests
bloodtests <- bloodtests[ClusterID %in% cohort_cluster_ids & TestName %in% relevant_tests]

if (nrow(bloodtests) == 0) stop("No blood test data found!")

# Convert datetime
bloodtests[, TestDateTime := as.POSIXct(CollectionDateTime, tz = "UTC")]
bloodtests <- bloodtests[!is.na(TestDateTime)]
bloodtests[, CollectionDateTime := NULL]

# Clean non-numeric values
# Strategy: 
# - Below detection (<X): Replace with X/2 (half the detection limit)
# - Above range (>X): Set to NA (actual value unknown, could be much higher)
# - Flags (*, ]): Set to NA

# Handle below detection limits: <X -> X/2
bloodtests[Value %like% "^<", Value := as.character(as.numeric(gsub("^<", "", Value)) / 2)]

# Handle above range: >X -> NA
bloodtests[Value %like% "^>", Value := NA_character_]

# Handle flags: *, ] -> NA
bloodtests[Value %in% c("*", "]"), Value := NA_character_]

# Convert Value to numeric
bloodtests[, Value := as.numeric(Value)]

# ============================================================
# 2-4. FIND CLOSEST BLOOD TEST PER TYPE (OPTIMIZED ROLLING JOIN)
# ============================================================
# Match blood tests to blood cultures within ±24 hour windows
# Uses rolling join approach like BC_6 for speed (O(n log n) vs O(n*m))

# Prepare blood culture keys
bc_for_join <- bc_df[, .(BloodCultureID, ClusterID, CollectionDateTime)]
bc_for_join[, CollectionDateTime := as.POSIXct(CollectionDateTime, tz = "UTC")]
bc_for_join[, JoinTime := as.numeric(CollectionDateTime)]

# Prepare blood tests with numeric time
bloodtests[, TestNumeric := as.numeric(TestDateTime)]

# Get unique test times per ClusterID per TestName
test_times <- unique(bloodtests[, .(ClusterID, TestName, JoinTime = TestNumeric, ActualTestTime = TestNumeric)])

# Create expanded bc_for_join with all test types (cross join)
test_types <- unique(bloodtests$TestName)
bc_for_join_expanded <- bc_for_join[, .(TestName = test_types), by = .(BloodCultureID, ClusterID, JoinTime)]

# Rolling join to find nearest test per ClusterID per TestName (ONE operation, not 15!)
setkey(test_times, ClusterID, TestName, JoinTime)
setkey(bc_for_join_expanded, ClusterID, TestName, JoinTime)

matched_times <- test_times[bc_for_join_expanded, roll = "nearest", mult = "first"]

# Calculate time difference and filter to -24h to +1h
# TIME WINDOW: 24 hours before, 1 hour after culture
# This avoids data leakage while allowing concurrent workup
matched_times[, HoursFromCulture := (ActualTestTime - JoinTime) / 3600]
matched_times <- matched_times[!is.na(HoursFromCulture) & HoursFromCulture >= -24 & HoursFromCulture <= 1]

# Join back to get actual test values
setkey(bloodtests, ClusterID, TestName, TestNumeric)
setkey(matched_times, ClusterID, TestName, ActualTestTime)
bc_bloodtests <- bloodtests[matched_times, on = .(ClusterID, TestName, TestNumeric = ActualTestTime), nomatch = 0L]

# Clean up columns
bc_bloodtests[, `:=`(TestDateTime = as.POSIXct(TestNumeric, origin = "1970-01-01", tz = "UTC"))]

# Abbreviate test names (exact names from IORD data)
test_abbrev <- c("CRP" = "CRP", "WHITE CELLS" = "WCC", "NEUTROPHILS" = "Neut",
                 "LYMPHOCYTES" = "Lymph", "CREATININE" = "Creat", "UREA" = "Urea",
                 "BILIRUBIN" = "Bili", "ALT" = "ALT", "ALBUMIN" = "Alb",
                 "HAEMOGLOBIN" = "Hb", "PLATELETS" = "Plt", "INR" = "INR",
                 "APTT" = "APTT", "SODIUM" = "Na", "POTASSIUM" = "K",
                 # ALL WCC components
                 "EOSINOPHILS" = "Eosin",
                 "MONOCYTES" = "Mono",
                 "BASOPHILS" = "Baso",
                 # Additional tests (exact IORD names - note periods and capitalization!)
                 "MEAN CELL VOL." = "MCV",        # Has period!
                 "HAEMATOCRIT" = "HCT", 
                 "ALK.PHOSPHATASE" = "ALP",       # Has period!
                 "Fibrinogen" = "Fibrinogen",     # Capitalized!
                 "D-Dimer" = "DDimer",            # Capitalized!
                 "GLUCOSE" = "Glucose")

bc_bloodtests[, TestName_Short := test_abbrev[TestName]]

# Reshape to wide format (rolling join already selected closest measurement)
baseline_summary <- dcast(bc_bloodtests,
                         BloodCultureID ~ TestName_Short,
                         value.var = c("Value", "HoursFromCulture"),
                         fun.aggregate = function(x) x[1])

# Rename columns to match expected format: Value_X -> BT_Closest_X, HoursFromCulture_X -> BT_Closest_X_HoursFromCulture
old_names <- names(baseline_summary)
new_names <- old_names

# Replace Value_ prefix
new_names <- gsub("^Value_", "BT_Closest_", new_names)

# Replace HoursFromCulture_ prefix  
new_names <- gsub("^HoursFromCulture_", "BT_Closest_", new_names)
  
# Add _HoursFromCulture suffix to the hours columns (those that came from HoursFromCulture_)
hours_indices <- grep("^HoursFromCulture_", old_names)
if (length(hours_indices) > 0) {
  new_names[hours_indices] <- paste0(new_names[hours_indices], "_HoursFromCulture")
}

setnames(baseline_summary, old_names, new_names)

# Ensure all expected columns exist (handles test types with no data in cohort)
for (nm in names(test_abbrev)) {
  abbr <- test_abbrev[nm]
  target_val <- paste0("BT_Closest_", abbr)
  target_hrs <- paste0("BT_Closest_", abbr, "_HoursFromCulture")
  if (!target_val %in% names(baseline_summary)) baseline_summary[, (target_val) := NA_real_]
  if (!target_hrs %in% names(baseline_summary)) baseline_summary[, (target_hrs) := NA_real_]
}

# ============================================================
# 5. CALCULATE MEASUREMENT COUNTS
# ============================================================

# Count how many tests have data per culture (from baseline_summary)
# Exclude the _HoursFromCulture columns, only count actual test values
bt_cols <- grep("^BT_Closest_[A-Za-z]+$", names(baseline_summary), value = TRUE)
bt_cols <- bt_cols[!grepl("_HoursFromCulture$", bt_cols)]

measurement_counts <- baseline_summary[, .(
  BloodCultureID,
  BT_NumMeasurements = rowSums(!is.na(.SD))
), .SDcols = bt_cols]

# ============================================================
# 6. JOIN BLOOD TEST DATA BACK TO BLOOD CULTURES
# ============================================================

bc_df[, CollectionDateTime := as.POSIXct(CollectionDateTime, tz = "UTC")]

# Deduplicate baseline_summary and measurement_counts first
baseline_summary <- unique(baseline_summary, by = "BloodCultureID")
measurement_counts <- unique(measurement_counts, by = "BloodCultureID")

setkey(bc_df, BloodCultureID)
setkey(baseline_summary, BloodCultureID)
setkey(measurement_counts, BloodCultureID)

# Left joins to keep all bc_df rows
bc_df <- merge(bc_df, baseline_summary, all.x = TRUE, by = "BloodCultureID")
bc_df <- merge(bc_df, measurement_counts, all.x = TRUE, by = "BloodCultureID")

# Final deduplication check
bc_df <- unique(bc_df, by = "BloodCultureID")

# Set BT_NumMeasurements to 0 where NA (cultures with no matching blood tests)
if ("BT_NumMeasurements" %in% names(bc_df)) {
  bc_df[is.na(BT_NumMeasurements), BT_NumMeasurements := 0L]
} else {
  bc_df[, BT_NumMeasurements := 0L]
}

# ============================================================
# 7. CATEGORIZE BASED ON DATA AVAILABILITY
# ============================================================

bc_df[, BT_DataAvailable := BT_NumMeasurements > 0]

# ============================================================
# 8. SAVE OUTPUT
# ============================================================

# Remove any list columns before saving (can cause fwrite errors)
list_cols <- names(bc_df)[sapply(bc_df, is.list)]
if (length(list_cols) > 0) {
  bc_df[, (list_cols) := NULL]
}

fwrite(bc_df, "Output Files/CSVs/blood_cultures_with_blood_tests.csv")

# ============================================================
# 9. CREATE EXCEL SUMMARY 
# ============================================================

# Helper function to calculate summary for all tests
calc_test_summary <- function(dt, category_name) {
  # Get all BT_Closest_ columns
  bt_cols <- grep("^BT_Closest_[A-Za-z]+$", names(dt), value = TRUE)
  bt_cols <- bt_cols[!grepl("_HoursFromCulture$", bt_cols)]
  
  result <- list(
    Category = category_name,
    N_Cultures = nrow(dt)
  )
  
  # Calculate N and Pct for each test
  for (col in bt_cols) {
    test_name <- gsub("^BT_Closest_", "", col)
    result[[paste0(test_name, "_N")]] <- sum(!is.na(dt[[col]]), na.rm = TRUE)
    result[[paste0(test_name, "_Pct")]] <- round(100 * mean(!is.na(dt[[col]]), na.rm = TRUE), 1)
  }
  
  return(as.data.table(result))
}

# Total summary (using helper function)
total_summary <- calc_test_summary(bc_df, "Total")

# By Positive/Negative (using helper function)
by_pos_neg_list <- lapply(unique(bc_df[!is.na(CultureCategory), CultureCategory]), function(cat) {
  calc_test_summary(bc_df[CultureCategory == cat], cat)
})
by_pos_neg <- rbindlist(by_pos_neg_list)
# Sort: Positive first, then Negative
by_pos_neg <- by_pos_neg[order(-Category)]

# Combine summaries (fill=TRUE to handle column differences)
combined_summary <- rbindlist(list(total_summary, by_pos_neg), fill = TRUE)

# Format output (convert to data.frame for Excel formatting compatibility)
combined_summary_formatted <- as.data.frame(combined_summary) %>% format_excel_output()

# Create info table
info_table <- data.frame(
  Item = c("Analysis", "Total cultures", "Cultures with blood tests", "Coverage", "Date", "Note"),
  Value = c("Blood Test Parameters",
            format_count(nrow(bc_df)), 
            format_count(sum(bc_df$BT_DataAvailable, na.rm = TRUE)),
            paste0(round(100 * mean(bc_df$BT_DataAvailable, na.rm = TRUE), 1), "%"),
            as.character(Sys.Date()),
            "Baseline values: closest measurement within ±24 hours per test")
)

# Create workbook with multiple sheets
wb <- createWorkbook()

# Sheet 1: Blood Test Summary
addWorksheet(wb, "Blood_Test_Summary")
writeData(wb, "Blood_Test_Summary", 
          data.frame(Title = "Blood Test Parameters by Positive/Negative Cultures"), 
          startRow = 1, colNames = FALSE)
writeData(wb, "Blood_Test_Summary", combined_summary_formatted, startRow = 3)

# Sheet 2: Info
addWorksheet(wb, "Info")
writeData(wb, "Info", info_table)

# Sheet 3: Formatting Note
addWorksheet(wb, "Formatting_Note")
writeData(wb, "Formatting_Note", create_formatting_note())

saveWorkbook(wb, "Output Files/Excels/Pipeline/BC_7_Add_Blood_Tests_excel.xlsx", overwrite = TRUE)

