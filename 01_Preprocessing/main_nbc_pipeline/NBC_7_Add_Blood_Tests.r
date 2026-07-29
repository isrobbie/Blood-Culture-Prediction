library(data.table)
library(lubridate)

# ============================================================
# ADD BLOOD TEST RESULTS TO NON-BC ADMISSIONS
# ============================================================
# Equivalent to BC_7. Uses a rolling nearest join per test type
# to find the closest value within +-24h of AdmissionDate.
# ============================================================

nbc_df <- fread("Output Files/CSVs/NonBC_Admissions/nbc_with_blood_gas.csv",
                colClasses = c(ClusterID = "character"))
nbc_df[, AdmissionDate := as.POSIXct(AdmissionDate, tz = "UTC")]

cohort_ids <- unique(nbc_df$ClusterID)

# ============================================================
# 1. LOAD AND PREPARE BLOOD TEST DATA
# ============================================================

relevant_tests <- c("CRP", "WHITE CELLS", "NEUTROPHILS", "LYMPHOCYTES",
                    "ALBUMIN", "PLATELETS", "CREATININE", "BILIRUBIN",
                    "UREA", "HAEMOGLOBIN", "SODIUM",
                    "EOSINOPHILS", "MONOCYTES", "ALK.PHOSPHATASE", "HAEMATOCRIT",
                    "ALT", "APTT", "BASOPHILS", "D-Dimer", "Fibrinogen",
                    "GLUCOSE", "INR", "POTASSIUM", "MEAN CELL VOL.")

bloodtests <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_BloodTests.csv",
                    select = c("ClusterID", "CollectionDateTime", "TestName", "Value"),
                    colClasses = c(ClusterID = "character"))
bloodtests <- bloodtests[ClusterID %in% cohort_ids & TestName %in% relevant_tests]
if (nrow(bloodtests) == 0) stop("No blood test data found!")
bloodtests[, TestDateTime := suppressWarnings(ymd_hms(CollectionDateTime, tz = "UTC"))]
bloodtests <- bloodtests[!is.na(TestDateTime)]
bloodtests[, CollectionDateTime := NULL]

bloodtests[Value %like% "^<", Value := as.character(as.numeric(gsub("^<", "", Value)) / 2)]
bloodtests[Value %like% "^>", Value := NA_character_]
bloodtests[Value %in% c("*", "]"), Value := NA_character_]
bloodtests[, Value := as.numeric(Value)]
bloodtests[, TestNumeric := as.numeric(TestDateTime)]

# ============================================================
# 2. ROLLING NEAREST JOIN PER TEST TYPE (+-24H)
# ============================================================

nbc_bt <- nbc_df[, .(AdmissionID, ClusterID, JoinTime = as.numeric(AdmissionDate))]
test_types  <- unique(bloodtests$TestName)
nbc_bt_exp  <- nbc_bt[, .(TestName = test_types), by = .(AdmissionID, ClusterID, JoinTime)]

test_times <- unique(bloodtests[, .(ClusterID, TestName, JoinTime = TestNumeric, ActualTestTime = TestNumeric)])
setkey(test_times, ClusterID, TestName, JoinTime)
setkey(nbc_bt_exp, ClusterID, TestName, JoinTime)

bt_matched <- test_times[nbc_bt_exp, roll = "nearest", mult = "first"]
bt_matched[, HoursFromAdmission := (ActualTestTime - JoinTime) / 3600]
bt_matched <- bt_matched[!is.na(HoursFromAdmission) & HoursFromAdmission >= -24 & HoursFromAdmission <= 1]

setkey(bloodtests, ClusterID, TestName, TestNumeric)
setkey(bt_matched, ClusterID, TestName, ActualTestTime)
bt_values <- bloodtests[bt_matched, on = .(ClusterID, TestName, TestNumeric = ActualTestTime), nomatch = 0L]

# ============================================================
# 3. RESHAPE AND RENAME TO BT_CLOSEST_* COLUMNS
# ============================================================

test_abbrev <- c("CRP" = "CRP", "WHITE CELLS" = "WCC", "NEUTROPHILS" = "Neut",
                 "LYMPHOCYTES" = "Lymph", "ALBUMIN" = "Alb", "PLATELETS" = "Plt",
                 "CREATININE" = "Creat", "BILIRUBIN" = "Bili", "UREA" = "Urea",
                 "HAEMOGLOBIN" = "Hb", "SODIUM" = "Na",
                 "EOSINOPHILS" = "Eosin", "MONOCYTES" = "Mono",
                 "ALK.PHOSPHATASE" = "ALP", "HAEMATOCRIT" = "HCT",
                 "ALT" = "ALT", "APTT" = "APTT", "BASOPHILS" = "Baso",
                 "D-Dimer" = "DDimer", "Fibrinogen" = "Fibrinogen",
                 "GLUCOSE" = "Glucose", "INR" = "INR",
                 "POTASSIUM" = "K", "MEAN CELL VOL." = "MCV")
bt_values[, TestName_Short := test_abbrev[TestName]]

bt_wide <- dcast(bt_values, AdmissionID ~ TestName_Short,
                 value.var = c("Value", "HoursFromAdmission"),
                 fun.aggregate = function(x) x[1])

# Rename: Value_X -> BT_Closest_X, HoursFromAdmission_X -> BT_Closest_X_HoursFromAdmission
old_names <- names(bt_wide)
new_names <- old_names
new_names <- gsub("^Value_", "BT_Closest_", new_names)
new_names <- gsub("^HoursFromAdmission_", "BT_Closest_", new_names)
hours_indices <- grep("^HoursFromAdmission_", old_names)
if (length(hours_indices) > 0) {
  new_names[hours_indices] <- paste0(new_names[hours_indices], "_HoursFromAdmission")
}
setnames(bt_wide, old_names, new_names)

# Ensure all expected columns exist (handles test types with no data in cohort)
for (nm in names(test_abbrev)) {
  abbr <- test_abbrev[nm]
  target_val <- paste0("BT_Closest_", abbr)
  target_hrs <- paste0("BT_Closest_", abbr, "_HoursFromAdmission")
  if (!target_val %in% names(bt_wide)) bt_wide[, (target_val) := NA_real_]
  if (!target_hrs %in% names(bt_wide)) bt_wide[, (target_hrs) := NA_real_]
}

bt_out_cols <- c("AdmissionID",
                 paste0("BT_Closest_", unname(test_abbrev)),
                 paste0("BT_Closest_", unname(test_abbrev), "_HoursFromAdmission"))
bt_out <- bt_wide[, ..bt_out_cols]

# Count distinct test types with data per admission (matches BC_7 approach)
bt_value_cols <- grep("^BT_Closest_[A-Za-z]+$", names(bt_out), value = TRUE)
bt_out[, BT_NumMeasurements := rowSums(!is.na(.SD)), .SDcols = bt_value_cols]
bt_out[, BT_DataAvailable := BT_NumMeasurements > 0]

nbc_df <- merge(nbc_df, bt_out, by = "AdmissionID", all.x = TRUE)
nbc_df[is.na(BT_NumMeasurements), BT_NumMeasurements := 0L]
nbc_df[is.na(BT_DataAvailable),   BT_DataAvailable   := FALSE]

rm(bloodtests, bt_matched, bt_values, bt_wide)

# ============================================================
# 4. SAVE OUTPUT
# ============================================================

fwrite(nbc_df, "Output Files/CSVs/NonBC_Admissions/nbc_with_blood_tests.csv")
