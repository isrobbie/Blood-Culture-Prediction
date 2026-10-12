library(data.table)
library(dplyr)
library(lubridate)
library(openxlsx)
library(writexl)

source("Active Scripts/Prepare/Excel_Output_Helper.r")

# ============================================================
# TRIAL: ADD BLOOD TEST RESULTS TO BLOOD CULTURES
# ============================================================
# Identical to main BC_7 except: tests must be >=2 hours BEFORE culture
# Window: -48h to -2h (no concurrent or post-culture tests)
# Uses shifted-LOCF rolling join to find the most recent valid test
# ============================================================

LOOKBACK_HOURS <- 48
trial_prefix   <- "trial_"

bc_df <- fread(paste0("Output Files/CSVs/", trial_prefix, "blood_cultures_with_blood_gas.csv"),
               colClasses = c(ClusterID = "character"))
cohort_cluster_ids <- unique(bc_df$ClusterID)

# ============================================================
# 1. LOAD AND PREPARE BLOOD TEST DATA (OPTIMIZED - FILTER EARLY)
# ============================================================

relevant_tests <- c(
  "CRP", "WHITE CELLS", "NEUTROPHILS", "LYMPHOCYTES",
  "CREATININE", "UREA", "BILIRUBIN", "ALT", "ALBUMIN",
  "HAEMOGLOBIN", "PLATELETS", "INR", "APTT",
  "SODIUM", "POTASSIUM",
  "EOSINOPHILS",
  "MONOCYTES",
  "BASOPHILS",
  "MEAN CELL VOL.",
  "HAEMATOCRIT",
  "ALK.PHOSPHATASE",
  "GLUCOSE",
  "Fibrinogen",
  "D-Dimer"
)

bloodtests <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_BloodTests.csv",
                    select = c("ClusterID", "CollectionDateTime", "TestName", "Value"),
                    colClasses = c(ClusterID = "character"))

bloodtests <- bloodtests[ClusterID %in% cohort_cluster_ids & TestName %in% relevant_tests]

if (nrow(bloodtests) == 0) stop("No blood test data found!")

bloodtests[, TestDateTime := as.POSIXct(CollectionDateTime, tz = "UTC")]
bloodtests <- bloodtests[!is.na(TestDateTime)]
bloodtests[, CollectionDateTime := NULL]

bloodtests[Value %like% "^<", Value := as.character(as.numeric(gsub("^<", "", Value)) / 2)]
bloodtests[Value %like% "^>", Value := NA_character_]
bloodtests[Value %in% c("*", "]"), Value := NA_character_]
bloodtests[, Value := as.numeric(Value)]

# ============================================================
# 2-4. FIND CLOSEST BLOOD TEST PER TYPE (SHIFTED LOCF)
# ============================================================
# TRIAL: shift culture time back by 2 hours, then roll=TRUE finds the
# latest test <= (CultureTime - 2h). Then filter to within 48h.

bc_for_join <- bc_df[, .(BloodCultureID, ClusterID, CollectionDateTime)]
bc_for_join[, CollectionDateTime := as.POSIXct(CollectionDateTime, tz = "UTC")]
bc_for_join[, CultureTimeNum := as.numeric(CollectionDateTime)]
bc_for_join[, JoinTime := CultureTimeNum - 2*3600]

bloodtests[, TestNumeric := as.numeric(TestDateTime)]

test_times <- unique(bloodtests[, .(ClusterID, TestName, JoinTime = TestNumeric, ActualTestTime = TestNumeric)])

test_types <- unique(bloodtests$TestName)
bc_for_join_expanded <- bc_for_join[, .(TestName = test_types), 
                                     by = .(BloodCultureID, ClusterID, CultureTimeNum, JoinTime)]

setkey(test_times, ClusterID, TestName, JoinTime)
setkey(bc_for_join_expanded, ClusterID, TestName, JoinTime)

matched_times <- test_times[bc_for_join_expanded, roll = TRUE, mult = "first"]

matched_times[, HoursFromCulture := (ActualTestTime - CultureTimeNum) / 3600]
matched_times <- matched_times[!is.na(HoursFromCulture) & HoursFromCulture >= -LOOKBACK_HOURS & HoursFromCulture <= -2]

setkey(bloodtests, ClusterID, TestName, TestNumeric)
setkey(matched_times, ClusterID, TestName, ActualTestTime)
bc_bloodtests <- bloodtests[matched_times, on = .(ClusterID, TestName, TestNumeric = ActualTestTime), nomatch = 0L]

bc_bloodtests[, `:=`(TestDateTime = as.POSIXct(TestNumeric, origin = "1970-01-01", tz = "UTC"))]

test_abbrev <- c("CRP" = "CRP", "WHITE CELLS" = "WCC", "NEUTROPHILS" = "Neut",
                 "LYMPHOCYTES" = "Lymph", "CREATININE" = "Creat", "UREA" = "Urea",
                 "BILIRUBIN" = "Bili", "ALT" = "ALT", "ALBUMIN" = "Alb",
                 "HAEMOGLOBIN" = "Hb", "PLATELETS" = "Plt", "INR" = "INR",
                 "APTT" = "APTT", "SODIUM" = "Na", "POTASSIUM" = "K",
                 "EOSINOPHILS" = "Eosin",
                 "MONOCYTES" = "Mono",
                 "BASOPHILS" = "Baso",
                 "MEAN CELL VOL." = "MCV",
                 "HAEMATOCRIT" = "HCT", 
                 "ALK.PHOSPHATASE" = "ALP",
                 "Fibrinogen" = "Fibrinogen",
                 "D-Dimer" = "DDimer",
                 "GLUCOSE" = "Glucose")

bc_bloodtests[, TestName_Short := test_abbrev[TestName]]

baseline_summary <- dcast(bc_bloodtests,
                         BloodCultureID ~ TestName_Short,
                         value.var = c("Value", "HoursFromCulture"),
                         fun.aggregate = function(x) x[1])

old_names <- names(baseline_summary)
new_names <- old_names

new_names <- gsub("^Value_", "BT_Closest_", new_names)
new_names <- gsub("^HoursFromCulture_", "BT_Closest_", new_names)

hours_indices <- grep("^HoursFromCulture_", old_names)
if (length(hours_indices) > 0) {
  new_names[hours_indices] <- paste0(new_names[hours_indices], "_HoursFromCulture")
}

setnames(baseline_summary, old_names, new_names)

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

baseline_summary <- unique(baseline_summary, by = "BloodCultureID")
measurement_counts <- unique(measurement_counts, by = "BloodCultureID")

setkey(bc_df, BloodCultureID)
setkey(baseline_summary, BloodCultureID)
setkey(measurement_counts, BloodCultureID)

bc_df <- merge(bc_df, baseline_summary, all.x = TRUE, by = "BloodCultureID")
bc_df <- merge(bc_df, measurement_counts, all.x = TRUE, by = "BloodCultureID")

bc_df <- unique(bc_df, by = "BloodCultureID")

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

list_cols <- names(bc_df)[sapply(bc_df, is.list)]
if (length(list_cols) > 0) {
  bc_df[, (list_cols) := NULL]
}

fwrite(bc_df, paste0("Output Files/CSVs/", trial_prefix, "blood_cultures_with_blood_tests.csv"))

# ============================================================
# 9. CREATE EXCEL SUMMARY 
# ============================================================

calc_test_summary <- function(dt, category_name) {
  bt_cols <- grep("^BT_Closest_[A-Za-z]+$", names(dt), value = TRUE)
  bt_cols <- bt_cols[!grepl("_HoursFromCulture$", bt_cols)]
  
  result <- list(
    Category = category_name,
    N_Cultures = nrow(dt)
  )
  
  for (col in bt_cols) {
    test_name <- gsub("^BT_Closest_", "", col)
    result[[paste0(test_name, "_N")]] <- sum(!is.na(dt[[col]]), na.rm = TRUE)
    result[[paste0(test_name, "_Pct")]] <- round(100 * mean(!is.na(dt[[col]]), na.rm = TRUE), 1)
  }
  
  return(as.data.table(result))
}

total_summary <- calc_test_summary(bc_df, "Total")

by_pos_neg_list <- lapply(unique(bc_df[!is.na(CultureCategory), CultureCategory]), function(cat) {
  calc_test_summary(bc_df[CultureCategory == cat], cat)
})
by_pos_neg <- rbindlist(by_pos_neg_list)
by_pos_neg <- by_pos_neg[order(-Category)]

combined_summary <- rbindlist(list(total_summary, by_pos_neg), fill = TRUE)

combined_summary_formatted <- as.data.frame(combined_summary) %>% format_excel_output()

info_table <- data.frame(
  Item = c("Analysis", "Total cultures", "Cultures with blood tests", "Coverage", "Date", "Note"),
  Value = c("TRIAL: Blood Test Parameters (-48h to -2h before culture)",
            format_count(nrow(bc_df)), 
            format_count(sum(bc_df$BT_DataAvailable, na.rm = TRUE)),
            paste0(round(100 * mean(bc_df$BT_DataAvailable, na.rm = TRUE), 1), "%"),
            as.character(Sys.Date()),
            "Baseline values: closest measurement within -48h to -2h before culture")
)

wb <- createWorkbook()

addWorksheet(wb, "Blood_Test_Summary")
writeData(wb, "Blood_Test_Summary", 
          data.frame(Title = "TRIAL: Blood Test Parameters by Positive/Negative Cultures (-48h to -2h before culture)"), 
          startRow = 1, colNames = FALSE)
writeData(wb, "Blood_Test_Summary", combined_summary_formatted, startRow = 3)

addWorksheet(wb, "Info")
writeData(wb, "Info", info_table)

addWorksheet(wb, "Formatting_Note")
writeData(wb, "Formatting_Note", create_formatting_note())

saveWorkbook(wb, "Output Files/Excels/Trial_BC_7_Add_Blood_Tests_excel.xlsx", overwrite = TRUE)
