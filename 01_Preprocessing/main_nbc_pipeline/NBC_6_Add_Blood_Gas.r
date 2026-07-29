library(data.table)
library(lubridate)

# ============================================================
# ADD BLOOD GAS RESULTS TO NON-BC ADMISSIONS
# ============================================================
# Equivalent to BC_6. Uses a rolling nearest join to find the
# closest blood gas within +-24h of AdmissionDate.
# ============================================================

nbc_df <- fread("Output Files/CSVs/NonBC_Admissions/nbc_with_vitals.csv",
                colClasses = c(ClusterID = "character"))
nbc_df[, AdmissionDate := as.POSIXct(AdmissionDate, tz = "UTC")]

cohort_ids <- unique(nbc_df$ClusterID)

# ============================================================
# 1. LOAD AND PREPARE BLOOD GAS DATA
# ============================================================

bloodgases <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_BloodGases.csv",
                    select = c("ClusterID", "PerformedDateTime", "EventName", "EventResult"),
                    colClasses = c(ClusterID = "character"))
bloodgases <- bloodgases[ClusterID %in% cohort_ids]
bloodgases[, PerformedDateTime := suppressWarnings(ymd_hms(PerformedDateTime))]
bloodgases <- bloodgases[!is.na(PerformedDateTime)]
bloodgases <- bloodgases[EventName %in% c("pH (BG)", "cLAC (BG)", "pCO2 POC",
                                          "BE Std (BG)", "Glucose (BG)", "Bicarb (BG)", "Hb (BG)")]
bloodgases[, `:=`(PerformedNumeric = as.numeric(PerformedDateTime),
                   EventResult = as.numeric(EventResult))]

bloodgases[, EventName_Short := fcase(
  EventName == "pH (BG)",      "pH",
  EventName == "cLAC (BG)",    "Lactate",
  EventName == "pCO2 POC",     "pCO2",
  EventName == "BE Std (BG)",  "BaseExcess",
  EventName == "Glucose (BG)", "Glucose",
  EventName == "Bicarb (BG)",  "Bicarb",
  EventName == "Hb (BG)",      "Hb",
  default = "Other")]

# ============================================================
# 2. ROLLING NEAREST JOIN TO ADMISSION (+-24H)
# ============================================================

nbc_bg <- nbc_df[, .(AdmissionID, ClusterID, JoinTime = as.numeric(AdmissionDate))]

gas_times <- unique(bloodgases[, .(ClusterID, JoinTime = PerformedNumeric, ActualGasTime = PerformedNumeric)])
setkey(gas_times, ClusterID, JoinTime)
setkey(nbc_bg, ClusterID, JoinTime)

bg_matched <- gas_times[nbc_bg, roll = "nearest", mult = "first"]
bg_matched[, HoursFromAdmission := (ActualGasTime - JoinTime) / 3600]
bg_matched <- bg_matched[!is.na(HoursFromAdmission) & HoursFromAdmission >= -24 & HoursFromAdmission <= 1]

setkey(bloodgases, ClusterID, PerformedNumeric)
setkey(bg_matched, ClusterID, ActualGasTime)
bg_values <- bloodgases[bg_matched, on = .(ClusterID, PerformedNumeric = ActualGasTime), nomatch = 0L]

bg_wide <- dcast(bg_values, AdmissionID + HoursFromAdmission ~ EventName_Short, value.var = "EventResult", fun.aggregate = function(x) x[1])

bg_cols_map <- c("Lactate"    = "BG_Closest_Lactate",
                  "pH"         = "BG_Closest_pH",
                  "pCO2"       = "BG_Closest_pCO2",
                  "BaseExcess" = "BG_Closest_BaseExcess",
                  "Glucose"    = "BG_Closest_Glucose",
                  "Bicarb"     = "BG_Closest_Bicarb",
                  "Hb"         = "BG_Closest_Hb")
for (old_nm in names(bg_cols_map)) {
  new_nm <- bg_cols_map[old_nm]
  if (old_nm %in% names(bg_wide)) setnames(bg_wide, old_nm, new_nm)
  else bg_wide[, (new_nm) := NA_real_]
}
bg_out <- bg_wide[, c("AdmissionID", "HoursFromAdmission", unname(bg_cols_map)), with = FALSE]
setnames(bg_out, "HoursFromAdmission", "BG_HoursFromAdmission")

bg_out <- bg_out[!(is.na(BG_Closest_pH) & is.na(BG_Closest_Lactate) & is.na(BG_Closest_BaseExcess) &
                   is.na(BG_Closest_Glucose) & is.na(BG_Closest_Bicarb) & is.na(BG_Closest_pCO2) &
                   is.na(BG_Closest_Hb))]

nbc_df <- merge(nbc_df, bg_out, by = "AdmissionID", all.x = TRUE)
rm(bloodgases, nbc_bg, bg_matched, bg_values, bg_wide, bg_out)

# ============================================================
# 3. ACID-BASE CLASSIFICATION (matching BC_6)
# ============================================================

nbc_df[, BG_AcidBase_Status := fcase(
  is.na(BG_Closest_pH) | is.na(BG_Closest_Bicarb) | is.na(BG_Closest_pCO2), NA_character_,
  BG_Closest_pH >= 7.35 & BG_Closest_pH <= 7.45,                                "Normal",
  BG_Closest_pH < 7.35 & BG_Closest_Bicarb < 22 & BG_Closest_pCO2 <= 6.0,      "Metabolic Acidosis",
  BG_Closest_pH < 7.35 & BG_Closest_Bicarb < 22 & BG_Closest_pCO2 > 6.0,       "Mixed Acidosis",
  BG_Closest_pH < 7.35 & BG_Closest_pCO2 > 6.0  & BG_Closest_Bicarb >= 22,     "Respiratory Acidosis",
  BG_Closest_pH < 7.35 & BG_Closest_Bicarb < 22,                                "Metabolic Acidosis",
  BG_Closest_pH < 7.35 & BG_Closest_pCO2 > 6.0,                                 "Respiratory Acidosis",
  BG_Closest_pH < 7.35,                                                           "Acidosis (unclassified)",
  BG_Closest_pH > 7.45 & BG_Closest_Bicarb > 26 & BG_Closest_pCO2 >= 4.7,      "Metabolic Alkalosis",
  BG_Closest_pH > 7.45 & BG_Closest_pCO2 < 4.7  & BG_Closest_Bicarb <= 26,     "Respiratory Alkalosis",
  BG_Closest_pH > 7.45 & BG_Closest_Bicarb > 26 & BG_Closest_pCO2 < 4.7,       "Mixed Alkalosis",
  BG_Closest_pH > 7.45 & BG_Closest_Bicarb > 26,                                "Metabolic Alkalosis",
  BG_Closest_pH > 7.45 & BG_Closest_pCO2 < 4.7,                                 "Respiratory Alkalosis",
  BG_Closest_pH > 7.45,                                                           "Alkalosis (unclassified)",
  default = "Unknown"
)]

# ============================================================
# 4. AVAILABILITY FLAGS
# ============================================================

nbc_df[, `:=`(
  BG_DataAvailable        = !is.na(BG_Closest_pH) | !is.na(BG_Closest_Lactate),
  BG_pH_Available         = !is.na(BG_Closest_pH),
  BG_Lactate_Available    = !is.na(BG_Closest_Lactate),
  BG_BaseExcess_Available = !is.na(BG_Closest_BaseExcess),
  BG_Glucose_Available    = !is.na(BG_Closest_Glucose),
  BG_Bicarb_Available     = !is.na(BG_Closest_Bicarb),
  BG_pCO2_Available       = !is.na(BG_Closest_pCO2),
  BG_Hb_Available         = !is.na(BG_Closest_Hb),
  BG_AcidBase_Available   = !is.na(BG_AcidBase_Status)
)]

# ============================================================
# 5. SAVE OUTPUT
# ============================================================

fwrite(nbc_df, "Output Files/CSVs/NonBC_Admissions/nbc_with_blood_gas.csv")

