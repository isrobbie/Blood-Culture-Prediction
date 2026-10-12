library(data.table)
library(dplyr)
library(lubridate)
library(openxlsx)

source("Active Scripts/Prepare/Excel_Output_Helper.r")

# ============================================================
# TRIAL: ADD BLOOD GAS RESULTS TO BLOOD CULTURES
# ============================================================
# Identical to main BC_6 except: gases must be >=5 minutes BEFORE culture
# Window: -48h to -5min (no concurrent or post-culture gases)
# Uses shifted-LOCF rolling join to find the most recent valid gas
# ============================================================

LOOKBACK_HOURS <- 48
trial_prefix   <- "trial_"

bc_df <- fread(paste0("Output Files/CSVs/", trial_prefix, "blood_cultures_with_vital_signs.csv"),
               colClasses = c(ClusterID = "character"))
initial_bc_count <- nrow(bc_df)

cohort_cluster_ids <- unique(bc_df$ClusterID)

# ============================================================
# 1. LOAD AND PREPARE BLOOD GAS DATA
# ============================================================

bloodgases <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_BloodGases.csv",
                    select = c("ClusterID", "PerformedDateTime", "EventName", "EventResult"),
                    colClasses = c(ClusterID = "character"))
bloodgases <- bloodgases[ClusterID %in% cohort_cluster_ids]
bloodgases[, PerformedDateTime := suppressWarnings(ymd_hms(PerformedDateTime))]
bloodgases <- bloodgases[!is.na(PerformedDateTime)]

bloodgases <- bloodgases[EventName %in% c("pH (BG)", "cLAC (BG)", "BE Std (BG)",
                                           "Glucose (BG)", "Bicarb (BG)", "pCO2 POC", "Hb (BG)")]

bloodgases[, `:=`(
  PerformedNumeric = as.numeric(PerformedDateTime),
  EventResult      = as.numeric(EventResult)
)]

bloodgases[, EventName_Short := fcase(
  EventName == "pH (BG)",      "pH",
  EventName == "cLAC (BG)",    "Lactate",
  EventName == "BE Std (BG)",  "BaseExcess",
  EventName == "Glucose (BG)", "Glucose",
  EventName == "Bicarb (BG)",  "Bicarb",
  EventName == "pCO2 POC",     "pCO2",
  EventName == "Hb (BG)",      "Hb",
  default = "Other"
)]

# ============================================================
# 2. MATCH BLOOD GASES TO BLOOD CULTURES (SHIFTED LOCF)
# ============================================================
# TRIAL: shift culture time back by 5 minutes, then roll=TRUE finds the
# latest gas <= (CultureTime - 5min). Then filter to within 48h.

bc_for_join <- bc_df[, .(BloodCultureID, ClusterID,
                          CultureTimeNum = as.numeric(as.POSIXct(CollectionDateTime, tz = "UTC")))]
bc_for_join[, JoinTime := CultureTimeNum - (5/60)*3600]

# Anchor rolling join on model-relevant types only (pH, Lactate, pCO2).
model_gas_types <- c("pH (BG)", "cLAC (BG)", "pCO2 POC")
gas_times <- unique(bloodgases[EventName %in% model_gas_types,
                               .(ClusterID, JoinTime = PerformedNumeric, ActualGasTime = PerformedNumeric)])
setkey(gas_times, ClusterID, JoinTime)
setkey(bc_for_join, ClusterID, JoinTime)

matched_times <- gas_times[bc_for_join, roll = TRUE, mult = "first"]
matched_times[, HoursFromCulture := (ActualGasTime - CultureTimeNum) / 3600]
matched_times <- matched_times[!is.na(HoursFromCulture) & HoursFromCulture >= -LOOKBACK_HOURS & HoursFromCulture <= -(5/60)]

setkey(bloodgases, ClusterID, PerformedNumeric)
setkey(matched_times, ClusterID, ActualGasTime)
bc_bloodgas <- bloodgases[matched_times, on = .(ClusterID, PerformedNumeric = ActualGasTime), nomatch = 0L]

# ============================================================
# 3. RESHAPE TO WIDE FORMAT
# ============================================================

bc_bloodgas_wide <- dcast(bc_bloodgas,
                           BloodCultureID + HoursFromCulture ~ EventName_Short,
                           value.var = "EventResult",
                           fun.aggregate = function(x) x[1])

bg_cols_map <- c(pH         = "BG_Closest_pH",
                 Lactate    = "BG_Closest_Lactate",
                 BaseExcess = "BG_Closest_BaseExcess",
                 Glucose    = "BG_Closest_Glucose",
                 Bicarb     = "BG_Closest_Bicarb",
                 pCO2       = "BG_Closest_pCO2",
                 Hb         = "BG_Closest_Hb")
for (old_nm in names(bg_cols_map)) {
  new_nm <- bg_cols_map[old_nm]
  if (old_nm %in% names(bc_bloodgas_wide)) setnames(bc_bloodgas_wide, old_nm, new_nm)
  else bc_bloodgas_wide[, (new_nm) := NA_real_]
}

bg_out <- bc_bloodgas_wide[, c("BloodCultureID", "HoursFromCulture", unname(bg_cols_map)), with = FALSE]
setnames(bg_out, "HoursFromCulture", "BG_HoursFromCulture")

bg_out <- bg_out[!(is.na(BG_Closest_pH) & is.na(BG_Closest_Lactate) & is.na(BG_Closest_BaseExcess) &
                   is.na(BG_Closest_Glucose) & is.na(BG_Closest_Bicarb) & is.na(BG_Closest_pCO2) &
                   is.na(BG_Closest_Hb))]

# ============================================================
# 4. JOIN BLOOD GAS DATA BACK TO BLOOD CULTURES
# ============================================================

bc_df[, CollectionDateTime := as.POSIXct(CollectionDateTime, tz = "UTC")]
setkey(bc_df, BloodCultureID)
setkey(bg_out, BloodCultureID)
bc_df <- bg_out[bc_df]

# ============================================================
# 5. ACID-BASE CLASSIFICATION
# ============================================================

bc_df[, BG_AcidBase_Status := fcase(
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
# 6. AVAILABILITY FLAGS
# ============================================================

bc_df[, `:=`(
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

rm(bloodgases, bc_for_join, gas_times, matched_times, bc_bloodgas, bc_bloodgas_wide, bg_out)

# ============================================================
# 7. SAVE OUTPUT
# ============================================================

fwrite(bc_df, paste0("Output Files/CSVs/", trial_prefix, "blood_cultures_with_blood_gas.csv"))

# ============================================================
# 8. CREATE EXCEL SUMMARY
# ============================================================
total_summary <- bc_df %>%
  summarise(
    Category = "Total",
    N_Cultures = n(),
    BG_Coverage_Pct = round(100 * mean(BG_DataAvailable, na.rm = TRUE), 1),
    pH_N = sum(BG_pH_Available, na.rm = TRUE),
    pH_Pct = round(100 * mean(BG_pH_Available, na.rm = TRUE), 1),
    pH_Median = round(median(BG_Closest_pH, na.rm = TRUE), 3),
    pH_IQR = paste0(round(quantile(BG_Closest_pH, 0.25, na.rm = TRUE), 3), "-",
                    round(quantile(BG_Closest_pH, 0.75, na.rm = TRUE), 3)),
    Lactate_N = sum(BG_Lactate_Available, na.rm = TRUE),
    Lactate_Pct = round(100 * mean(BG_Lactate_Available, na.rm = TRUE), 1),
    Lactate_Median = round(median(BG_Closest_Lactate, na.rm = TRUE), 2),
    Lactate_IQR = paste0(round(quantile(BG_Closest_Lactate, 0.25, na.rm = TRUE), 2), "-",
                         round(quantile(BG_Closest_Lactate, 0.75, na.rm = TRUE), 2)),
    BaseExcess_N = sum(BG_BaseExcess_Available, na.rm = TRUE),
    BaseExcess_Pct = round(100 * mean(BG_BaseExcess_Available, na.rm = TRUE), 1),
    BaseExcess_Median = round(median(BG_Closest_BaseExcess, na.rm = TRUE), 2),
    BaseExcess_IQR = paste0(round(quantile(BG_Closest_BaseExcess, 0.25, na.rm = TRUE), 2), "-",
                            round(quantile(BG_Closest_BaseExcess, 0.75, na.rm = TRUE), 2)),
    Glucose_N = sum(BG_Glucose_Available, na.rm = TRUE),
    Glucose_Pct = round(100 * mean(BG_Glucose_Available, na.rm = TRUE), 1),
    Glucose_Median = round(median(BG_Closest_Glucose, na.rm = TRUE), 1),
    Glucose_IQR = paste0(round(quantile(BG_Closest_Glucose, 0.25, na.rm = TRUE), 1), "-",
                         round(quantile(BG_Closest_Glucose, 0.75, na.rm = TRUE), 1)),
    Bicarb_Median = round(median(BG_Closest_Bicarb, na.rm = TRUE), 1),
    pCO2_Median = round(median(BG_Closest_pCO2, na.rm = TRUE), 1),
    Hb_Median = round(median(BG_Closest_Hb, na.rm = TRUE), 1)
  )

by_pos_neg <- bc_df %>%
  filter(!is.na(CultureCategory)) %>%
  group_by(CultureCategory) %>%
  summarise(
    Category = first(CultureCategory),
    N_Cultures = n(),
    BG_Coverage_Pct = round(100 * mean(BG_DataAvailable, na.rm = TRUE), 1),
    pH_N = sum(BG_pH_Available, na.rm = TRUE),
    pH_Pct = round(100 * mean(BG_pH_Available, na.rm = TRUE), 1),
    pH_Median = round(median(BG_Closest_pH, na.rm = TRUE), 3),
    pH_IQR = paste0(round(quantile(BG_Closest_pH, 0.25, na.rm = TRUE), 3), "-",
                    round(quantile(BG_Closest_pH, 0.75, na.rm = TRUE), 3)),
    Lactate_N = sum(BG_Lactate_Available, na.rm = TRUE),
    Lactate_Pct = round(100 * mean(BG_Lactate_Available, na.rm = TRUE), 1),
    Lactate_Median = round(median(BG_Closest_Lactate, na.rm = TRUE), 2),
    Lactate_IQR = paste0(round(quantile(BG_Closest_Lactate, 0.25, na.rm = TRUE), 2), "-",
                         round(quantile(BG_Closest_Lactate, 0.75, na.rm = TRUE), 2)),
    BaseExcess_N = sum(BG_BaseExcess_Available, na.rm = TRUE),
    BaseExcess_Pct = round(100 * mean(BG_BaseExcess_Available, na.rm = TRUE), 1),
    BaseExcess_Median = round(median(BG_Closest_BaseExcess, na.rm = TRUE), 2),
    BaseExcess_IQR = paste0(round(quantile(BG_Closest_BaseExcess, 0.25, na.rm = TRUE), 2), "-",
                            round(quantile(BG_Closest_BaseExcess, 0.75, na.rm = TRUE), 2)),
    Glucose_N = sum(BG_Glucose_Available, na.rm = TRUE),
    Glucose_Pct = round(100 * mean(BG_Glucose_Available, na.rm = TRUE), 1),
    Glucose_Median = round(median(BG_Closest_Glucose, na.rm = TRUE), 1),
    Glucose_IQR = paste0(round(quantile(BG_Closest_Glucose, 0.25, na.rm = TRUE), 1), "-",
                         round(quantile(BG_Closest_Glucose, 0.75, na.rm = TRUE), 1)),
    Bicarb_Median = round(median(BG_Closest_Bicarb, na.rm = TRUE), 1),
    pCO2_Median = round(median(BG_Closest_pCO2, na.rm = TRUE), 1),
    Hb_Median = round(median(BG_Closest_Hb, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  select(-CultureCategory) %>%
  arrange(desc(Category))

combined_summary <- bind_rows(total_summary, by_pos_neg)

acidbase_by_pos_neg <- bc_df %>%
  filter(BG_AcidBase_Available & !is.na(CultureCategory)) %>%
  group_by(CultureCategory, BG_AcidBase_Status) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(Percentage = round(100 * Count / sum(Count), 1)) %>%
  ungroup()

acidbase_by_pos_neg <- acidbase_by_pos_neg %>%
  rename(Category = CultureCategory)

acidbase_total <- bc_df %>%
  filter(BG_AcidBase_Available & !is.na(CultureCategory)) %>%
  group_by(BG_AcidBase_Status) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(
    Category = "Total",
    Percentage = round(100 * Count / sum(Count), 1)
  ) %>%
  select(Category, BG_AcidBase_Status, Count, Percentage)

acidbase_by_pos_neg <- bind_rows(acidbase_total, acidbase_by_pos_neg)

status_order <- c("Normal", 
                  "Metabolic Acidosis", "Respiratory Acidosis", "Mixed Acidosis", "Acidosis (unclassified)",
                  "Metabolic Alkalosis", "Respiratory Alkalosis", "Mixed Alkalosis", "Alkalosis (unclassified)",
                  "Unknown")

acidbase_by_pos_neg <- acidbase_by_pos_neg %>%
  mutate(
    BG_AcidBase_Status = factor(BG_AcidBase_Status, levels = status_order),
    Category_Order = case_when(
      Category == "Total" ~ 1,
      Category == "Positive" ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(Category_Order, BG_AcidBase_Status) %>%
  select(-Category_Order) %>%
  mutate(BG_AcidBase_Status = as.character(BG_AcidBase_Status))

combined_summary_formatted <- combined_summary %>% format_excel_output()
acidbase_formatted <- acidbase_by_pos_neg %>% format_excel_output()

wb <- createWorkbook()
addWorksheet(wb, "Blood Gas Summary")
writeData(wb, "Blood Gas Summary", data.frame(Title = "TRIAL: Blood Gas Parameters by Positive/Negative Cultures (-48h to -5min before culture)"), startRow = 1, colNames = FALSE)
writeData(wb, "Blood Gas Summary", combined_summary_formatted, startRow = 3)

addWorksheet(wb, "Acid-Base Status")
writeData(wb, "Acid-Base Status", data.frame(Title = "TRIAL: Acid-Base Classification by Positive/Negative Cultures"), startRow = 1, colNames = FALSE)
writeData(wb, "Acid-Base Status", acidbase_formatted, startRow = 3)

addWorksheet(wb, "Formatting_Note")
writeData(wb, "Formatting_Note", create_formatting_note())

saveWorkbook(wb, "Output Files/Excels/Trial_BC_6_Add_Blood_Gas_excel.xlsx", overwrite = TRUE)
