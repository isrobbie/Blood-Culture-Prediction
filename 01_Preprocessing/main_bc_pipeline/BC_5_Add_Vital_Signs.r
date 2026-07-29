library(data.table)
library(dplyr)
library(lubridate)
library(stringr)
library(openxlsx)

source("01_Preprocessing/Excel_Output_Helper.r")

# ============================================================
# ADD VITAL SIGNS AND NEWS SCORES TO BLOOD CULTURES
# ============================================================
# Uses rolling join approach (like BC_6/BC_7) for speed

# Read blood culture data
bc_df <- fread("Output Files/CSVs/blood_cultures_with_infection_source.csv",
               colClasses = c(ClusterID = "character"))
initial_bc_count <- nrow(bc_df)
cohort_cluster_ids <- unique(bc_df$ClusterID)

# ============================================================
# 1. LOAD AND PREPARE VITAL SIGNS DATA
# ============================================================

vitalsigns1 <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_VitalSigns1.csv",
                     select = c("ClusterID", "PerformedDateTime", "EventName", "EventResult"),
                     colClasses = c(ClusterID = "character"))
vitalsigns2 <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_VitalSigns2.csv",
                     select = c("ClusterID", "PerformedDateTime", "EventName", "EventResult"),
                     colClasses = c(ClusterID = "character"))
oxygen_delivery <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_OxygenDelivery.csv",
                        select = c("ClusterID", "PerformedDateTime", "EventName", "EventResult"),
                        colClasses = c(ClusterID = "character"))

vitalsigns <- rbindlist(list(vitalsigns1, vitalsigns2, oxygen_delivery))
rm(vitalsigns1, vitalsigns2, oxygen_delivery)

vitalsigns <- vitalsigns[ClusterID %in% cohort_cluster_ids]
vitalsigns[, PerformedDateTime := as.POSIXct(PerformedDateTime, tz = "UTC")]
vitalsigns <- vitalsigns[!is.na(PerformedDateTime)]

# Filter to relevant vital signs
vitalsigns <- vitalsigns[
  EventName %in% c("Heart Rate", "Respiratory Rate", "Oxygen Saturation", 
                   "Oxygen L/min delivered", "Systolic Blood Pressure", 
                   "Diastolic Blood Pressure", "AVPU Scale", "Temperature Tympanic")
]

# Abbreviate event names early
vitalsigns[, EventName_Short := fcase(
  EventName == "Heart Rate", "HR",
  EventName == "Respiratory Rate", "RR",
  EventName == "Temperature Tympanic", "Temp",
  EventName == "Oxygen Saturation", "O2Sat",
  EventName == "Oxygen L/min delivered", "O2Flow",
  EventName == "Systolic Blood Pressure", "SBP",
  EventName == "Diastolic Blood Pressure", "DBP",
  EventName == "AVPU Scale", "AVPU",
  default = "Other"
)]

# ============================================================
# 2. MATCH VITAL SIGNS USING ROLLING JOIN (FAST)
# ============================================================

# Prepare blood culture data
bc_for_join <- bc_df[, .(ClusterID, BloodCultureID, AccessionNumber, CollectionDateTime, SpecimenFull, EpisodeID)]
bc_for_join[, CollectionDateTime := as.POSIXct(CollectionDateTime, tz = "UTC")]
bc_for_join[, JoinTime := as.numeric(CollectionDateTime)]

# Prepare vital signs with numeric time
vitalsigns[, PerformedNumeric := as.numeric(PerformedDateTime)]

# Expand bc_for_join with all vital types
vital_types <- unique(vitalsigns$EventName_Short)
bc_for_join_expanded <- bc_for_join[, .(EventName_Short = vital_types), 
                                     by = .(ClusterID, BloodCultureID, AccessionNumber, CollectionDateTime, 
                                            SpecimenFull, EpisodeID, JoinTime)]

bc_for_join_expanded[, `:=`(
  window_start = JoinTime - 24*3600,
  window_end   = JoinTime + 1*3600
)]

# Prepare vitals for join - save actual time before non-equi join overwrites it
vitalsigns[, ActualVitalTime := PerformedNumeric]
vitalsigns_join <- vitalsigns[, .(ClusterID, EventName_Short, PerformedNumeric, ActualVitalTime, EventResult)]

matched_vitals <- vitalsigns_join[bc_for_join_expanded,
                                   on = .(ClusterID, EventName_Short, 
                                          PerformedNumeric >= window_start,
                                          PerformedNumeric <= window_end),
                                   nomatch = 0L,
                                   allow.cartesian = TRUE]

matched_vitals[, HoursFromCulture := (ActualVitalTime - JoinTime) / 3600]

# For aggregation: keep ALL vitals in window
bc_vitals_all <- matched_vitals[, .(
  ClusterID, BloodCultureID, AccessionNumber, CollectionDateTime, SpecimenFull, EpisodeID, EventName_Short,
  ActualVitalTime, EventResult, HoursFromCulture
)]
bc_vitals_all[, `:=`(
  PerformedDateTime = as.POSIXct(ActualVitalTime, origin = "1970-01-01", tz = "UTC"),
  DaysFromCulture = HoursFromCulture / 24
)]

# For baseline: find closest vital per culture per EventName_Short
bc_vitals <- matched_vitals[, .SD[which.min(abs(HoursFromCulture))], 
                             by = .(ClusterID, BloodCultureID, AccessionNumber, CollectionDateTime, 
                                    SpecimenFull, EpisodeID, EventName_Short)]
bc_vitals[, `:=`(
  PerformedDateTime = as.POSIXct(ActualVitalTime, origin = "1970-01-01", tz = "UTC"),
  DaysFromCulture = HoursFromCulture / 24
)]

# ============================================================
# 3. RESHAPE TO WIDE FORMAT
# ============================================================

# For aggregation: use ALL vitals in window
bc_vitals_wide_all <- dcast(bc_vitals_all, 
                            ClusterID + BloodCultureID + AccessionNumber + CollectionDateTime + SpecimenFull + EpisodeID + 
                            PerformedDateTime + HoursFromCulture + DaysFromCulture ~ EventName_Short, 
                            value.var = "EventResult",
                            fun.aggregate = function(x) x[1])

process_vitals_wide <- function(dt) {
  # Convert numeric columns
  potential_numeric_cols <- c("HR", "RR", "Temp", "O2Sat", "O2Flow", "SBP", "DBP")
  numeric_cols <- intersect(potential_numeric_cols, names(dt))
  if (length(numeric_cols) > 0) {
    dt[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]
  }
  
  # Process AVPU
  if ("AVPU" %in% names(dt)) {
    dt[, AVPU := fifelse(str_detect(AVPU, "^A|^Alert"), 0L, 1L)]
  }
  
  # Process O2Flow
  if ("O2Flow" %in% names(dt)) {
    dt[, O2orAir := fifelse(!is.na(O2Flow) & O2Flow > 0, 1L, 0L)]
  } else {
    dt[, O2orAir := NA_integer_]
  }
  return(dt)
}

bc_vitals_wide_all <- process_vitals_wide(bc_vitals_wide_all)

# ============================================================
# 4. COMPUTE NEWS SCORES
# ============================================================

compute_news_scores <- function(dt) {
  if ("RR" %in% names(dt)) {
    dt[, RR_score := fcase(is.na(RR), NA_integer_, RR <= 8, 3L, RR <= 11, 1L, RR <= 20, 0L, RR <= 24, 2L, default = 3L)]
  } else { dt[, RR_score := NA_integer_] }
  
  if ("O2Sat" %in% names(dt)) {
    dt[, O2Sat_score := fcase(is.na(O2Sat), NA_integer_, O2Sat <= 91, 3L, O2Sat <= 93, 2L, O2Sat <= 95, 1L, default = 0L)]
  } else { dt[, O2Sat_score := NA_integer_] }
  
  if ("O2orAir" %in% names(dt)) {
    dt[, O2orAir_score := fcase(is.na(O2orAir), NA_integer_, O2orAir == 1, 2L, default = 0L)]
  } else { dt[, O2orAir_score := NA_integer_] }
  
  if ("SBP" %in% names(dt)) {
    dt[, SBP_score := fcase(is.na(SBP), NA_integer_, SBP <= 90, 3L, SBP <= 100, 2L, SBP <= 110, 1L, SBP <= 219, 0L, default = 3L)]
  } else { dt[, SBP_score := NA_integer_] }
  
  if ("HR" %in% names(dt)) {
    dt[, HR_score := fcase(is.na(HR), NA_integer_, HR <= 40, 3L, HR <= 50, 1L, HR <= 90, 0L, HR <= 110, 1L, HR <= 130, 2L, default = 3L)]
  } else { dt[, HR_score := NA_integer_] }
  
  if ("AVPU" %in% names(dt)) {
    dt[, AVPU_score := fcase(is.na(AVPU), 0L, AVPU == 0, 0L, default = 3L)]
  } else { dt[, AVPU_score := 0L] }
  
  if ("Temp" %in% names(dt)) {
    dt[, Temp_score := fcase(is.na(Temp), NA_integer_, Temp <= 35, 3L, Temp <= 36, 1L, Temp <= 38, 0L, Temp <= 39, 1L, default = 2L)]
  } else { dt[, Temp_score := NA_integer_] }
  
  score_cols <- c("RR_score", "O2Sat_score", "O2orAir_score", "SBP_score", "HR_score", "AVPU_score", "Temp_score")
  dt[, NEWS := rowSums(.SD, na.rm = TRUE), .SDcols = score_cols]
  dt[, n_scores_available := rowSums(!is.na(.SD)), .SDcols = score_cols]
  dt[n_scores_available < 3, NEWS := NA_integer_]
  return(dt)
}

bc_vitals_wide_all <- compute_news_scores(bc_vitals_wide_all)

# ============================================================
# 5. FIND BASELINE VITALS (CLOSEST PER VITAL TYPE - INDEPENDENT)
# ============================================================

# bc_vitals already has the closest measurement per vital type per culture
# Now extract each vital's value and time independently

# Process bc_vitals for baseline extraction
bc_vitals_processed <- copy(bc_vitals)

# Process AVPU before numeric coercion (text "Alert" etc. must be matched first)
bc_vitals_processed[EventName_Short == "AVPU", EventResult := fifelse(
  grepl("^A|^Alert", EventResult, ignore.case = TRUE), "0", "1")]

# Process O2Flow before numeric coercion
bc_vitals_processed[EventName_Short == "O2Flow", EventResult := fifelse(
  !is.na(suppressWarnings(as.numeric(EventResult))) & as.numeric(EventResult) > 0, "1", "0")]

bc_vitals_processed[, EventResult := as.numeric(EventResult)]

# Reshape to get each vital's baseline value (one column per vital type)
baseline_values <- dcast(bc_vitals_processed, 
                         ClusterID + BloodCultureID ~ EventName_Short,
                         value.var = "EventResult",
                         fun.aggregate = function(x) x[1])

# Also get hours from culture for each vital type
baseline_hours <- dcast(bc_vitals_processed,
                        ClusterID + BloodCultureID ~ EventName_Short,
                        value.var = "HoursFromCulture",
                        fun.aggregate = function(x) x[1])

# Create baseline summary with per-vital times
baseline_vitals_summary <- baseline_values[, .(
  ClusterID, BloodCultureID,
  VS_Baseline_HR = HR, VS_Baseline_RR = RR, VS_Baseline_Temp = Temp,
  VS_Baseline_O2Sat = O2Sat, VS_Baseline_SBP = SBP, VS_Baseline_DBP = DBP,
  VS_Baseline_AVPU = AVPU
)]

# Add O2orAir from O2Flow (0 = room air when no O2Flow data, matching EV approach)
if ("O2Flow" %in% names(baseline_values)) {
  baseline_vitals_summary[, VS_Baseline_O2orAir := fifelse(!is.na(baseline_values$O2Flow) & baseline_values$O2Flow > 0, 1L, 0L)]
} else {
  baseline_vitals_summary[, VS_Baseline_O2orAir := 0L]
}

# Add hours from culture for HR (representative - could add per-vital if needed)
if ("HR" %in% names(baseline_hours)) {
  baseline_vitals_summary[, VS_Baseline_HoursFromCulture := baseline_hours$HR]
} else {
  baseline_vitals_summary[, VS_Baseline_HoursFromCulture := NA_real_]
}

# Compute NEWS from baseline values (using the baseline vitals which may be from different times)
baseline_vitals_summary[, `:=`(
  RR_score = fcase(is.na(VS_Baseline_RR), NA_integer_, VS_Baseline_RR <= 8, 3L, VS_Baseline_RR <= 11, 1L, VS_Baseline_RR <= 20, 0L, VS_Baseline_RR <= 24, 2L, default = 3L),
  O2Sat_score = fcase(is.na(VS_Baseline_O2Sat), NA_integer_, VS_Baseline_O2Sat <= 91, 3L, VS_Baseline_O2Sat <= 93, 2L, VS_Baseline_O2Sat <= 95, 1L, default = 0L),
  O2orAir_score = fcase(is.na(VS_Baseline_O2orAir), NA_integer_, VS_Baseline_O2orAir == 1, 2L, default = 0L),
  SBP_score = fcase(is.na(VS_Baseline_SBP), NA_integer_, VS_Baseline_SBP <= 90, 3L, VS_Baseline_SBP <= 100, 2L, VS_Baseline_SBP <= 110, 1L, VS_Baseline_SBP <= 219, 0L, default = 3L),
  HR_score = fcase(is.na(VS_Baseline_HR), NA_integer_, VS_Baseline_HR <= 40, 3L, VS_Baseline_HR <= 50, 1L, VS_Baseline_HR <= 90, 0L, VS_Baseline_HR <= 110, 1L, VS_Baseline_HR <= 130, 2L, default = 3L),
  AVPU_score = fcase(is.na(VS_Baseline_AVPU), 0L, VS_Baseline_AVPU == 0, 0L, default = 3L),
  Temp_score = fcase(is.na(VS_Baseline_Temp), NA_integer_, VS_Baseline_Temp <= 35, 3L, VS_Baseline_Temp <= 36, 1L, VS_Baseline_Temp <= 38, 0L, VS_Baseline_Temp <= 39, 1L, default = 2L)
)]

score_cols <- c("RR_score", "O2Sat_score", "O2orAir_score", "SBP_score", "HR_score", "AVPU_score", "Temp_score")
baseline_vitals_summary[, VS_Baseline_NEWS := rowSums(.SD, na.rm = TRUE), .SDcols = score_cols]
baseline_vitals_summary[, n_scores := rowSums(!is.na(.SD)), .SDcols = score_cols]
baseline_vitals_summary[n_scores < 3, VS_Baseline_NEWS := NA_integer_]
baseline_vitals_summary[, c(score_cols, "n_scores") := NULL]

# ============================================================
# 6. AGGREGATE VITALS (USE ALL VITALS IN WINDOW)
# ============================================================

vitals_aggregated <- bc_vitals_wide_all[, .(
  VS_NumMeasurements = .N,
  VS_HR_Max = max(HR, na.rm = TRUE), VS_HR_Min = min(HR, na.rm = TRUE), VS_HR_Mean = mean(HR, na.rm = TRUE),
  VS_RR_Max = max(RR, na.rm = TRUE), VS_RR_Min = min(RR, na.rm = TRUE), VS_RR_Mean = mean(RR, na.rm = TRUE),
  VS_Temp_Max = max(Temp, na.rm = TRUE), VS_Temp_Min = min(Temp, na.rm = TRUE), VS_Temp_Mean = mean(Temp, na.rm = TRUE),
  VS_O2Sat_Max = max(O2Sat, na.rm = TRUE), VS_O2Sat_Min = min(O2Sat, na.rm = TRUE), VS_O2Sat_Mean = mean(O2Sat, na.rm = TRUE),
  VS_SBP_Max = max(SBP, na.rm = TRUE), VS_SBP_Min = min(SBP, na.rm = TRUE), VS_SBP_Mean = mean(SBP, na.rm = TRUE),
  VS_DBP_Max = max(DBP, na.rm = TRUE), VS_DBP_Min = min(DBP, na.rm = TRUE), VS_DBP_Mean = mean(DBP, na.rm = TRUE),
  VS_NEWS_Max = max(NEWS, na.rm = TRUE), VS_NEWS_Min = min(NEWS, na.rm = TRUE), VS_NEWS_Mean = mean(NEWS, na.rm = TRUE),
  VS_AnyO2Use = as.integer(any(O2orAir == 1, na.rm = TRUE)),
  VS_AnyFever_38 = as.integer(any(Temp >= 38.0, na.rm = TRUE)),
  VS_AnyFever_38_3 = as.integer(any(Temp >= 38.3, na.rm = TRUE)),
  VS_AnyHypothermia = as.integer(any(Temp < 36.0, na.rm = TRUE)),
  VS_AnyHypotension = as.integer(any(SBP < 90, na.rm = TRUE)),
  VS_AnyTachycardia = as.integer(any(HR > 100, na.rm = TRUE)),
  VS_AnyTachypnea = as.integer(any(RR > 20, na.rm = TRUE)),
  VS_AnyHypoxia = as.integer(any(O2Sat < 95, na.rm = TRUE))
), by = .(ClusterID, BloodCultureID)]

# Replace Inf/-Inf with NA
numeric_cols_agg <- names(vitals_aggregated)[sapply(vitals_aggregated, is.numeric)]
for (col in numeric_cols_agg) {
  vitals_aggregated[is.infinite(get(col)) | is.nan(get(col)), (col) := NA]
}

# ============================================================
# 7. JOIN BACK TO BLOOD CULTURES
# ============================================================

bc_df <- as_tibble(bc_df) %>%
  mutate(CollectionDateTime = as.POSIXct(CollectionDateTime, tz = "UTC"))

bc_df <- bc_df %>%
  left_join(baseline_vitals_summary %>% as_tibble(), 
            by = c("ClusterID", "BloodCultureID")) %>%
  left_join(vitals_aggregated %>% as_tibble(), 
            by = c("ClusterID", "BloodCultureID"))

bc_df <- bc_df %>%
  mutate(
    VS_DataAvailable = !is.na(VS_NumMeasurements) & VS_NumMeasurements > 0,
    VS_NEWSAvailable = !is.na(VS_Baseline_NEWS),
    VS_SIRS_Count = as.integer(!is.na(VS_Baseline_Temp) & (VS_Baseline_Temp < 36 | VS_Baseline_Temp > 38)) +
                    as.integer(!is.na(VS_Baseline_HR) & VS_Baseline_HR > 90) +
                    as.integer(!is.na(VS_Baseline_RR) & VS_Baseline_RR > 20),
    VS_qSOFA_Count = as.integer(!is.na(VS_Baseline_RR) & VS_Baseline_RR >= 22) +
                     as.integer(!is.na(VS_Baseline_SBP) & VS_Baseline_SBP <= 100) +
                     as.integer(!is.na(VS_Baseline_AVPU) & VS_Baseline_AVPU > 0)
  )

fwrite(bc_df, "Output Files/CSVs/blood_cultures_with_vital_signs.csv")

# ============================================================
# 8. SUMMARY STATISTICS
# ============================================================

final_bc_count <- nrow(bc_df)
if (initial_bc_count != final_bc_count) {
  warning("Blood culture count mismatch! Lost: ", initial_bc_count - final_bc_count)
}

create_vitals_summary <- function(data, group_name) {
  calc_stat <- function(cond) list(N = sum(cond, na.rm = TRUE), Pct = round(100 * mean(cond, na.rm = TRUE), 1))
  tibble(
    Group = group_name, N_Cultures = nrow(data),
    With_VitalSigns_N = calc_stat(data$VS_DataAvailable == TRUE)$N,
    With_VitalSigns_Pct = calc_stat(data$VS_DataAvailable == TRUE)$Pct,
    With_NEWS_N = calc_stat(data$VS_NEWSAvailable == TRUE)$N,
    With_NEWS_Pct = calc_stat(data$VS_NEWSAvailable == TRUE)$Pct,
    Baseline_HR_Median = round(median(data$VS_Baseline_HR, na.rm = TRUE), 1),
    Baseline_RR_Median = round(median(data$VS_Baseline_RR, na.rm = TRUE), 1),
    Baseline_Temp_Median = round(median(data$VS_Baseline_Temp, na.rm = TRUE), 2),
    Baseline_SBP_Median = round(median(data$VS_Baseline_SBP, na.rm = TRUE), 1),
    Baseline_NEWS_Median = round(median(data$VS_Baseline_NEWS, na.rm = TRUE), 1),
    Fever_38_Pct = calc_stat(data$VS_AnyFever_38 == 1)$Pct,
    Hypothermia_Pct = calc_stat(data$VS_AnyHypothermia == 1)$Pct,
    Hypotension_Pct = calc_stat(data$VS_AnyHypotension == 1)$Pct,
    Tachycardia_Pct = calc_stat(data$VS_AnyTachycardia == 1)$Pct,
    SIRS_2plus_Pct = calc_stat(data$VS_SIRS_Count >= 2)$Pct,
    qSOFA_2plus_Pct = calc_stat(data$VS_qSOFA_Count >= 2)$Pct
  )
}

summary_table <- bind_rows(
  create_vitals_summary(bc_df, "All Cultures"),
  create_vitals_summary(bc_df %>% filter(CultureCategory == "Negative (and SC)"), "Negative"),
  create_vitals_summary(bc_df %>% filter(CultureCategory == "Positive"), "Positive")
)

wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_table %>% format_excel_output(), startRow = 1)
addWorksheet(wb, "Formatting_Note")
writeData(wb, "Formatting_Note", create_formatting_note())
saveWorkbook(wb, "Output Files/Excels/Pipeline/BC_5_Add_Vital_Signs_excel.xlsx", overwrite = TRUE)
