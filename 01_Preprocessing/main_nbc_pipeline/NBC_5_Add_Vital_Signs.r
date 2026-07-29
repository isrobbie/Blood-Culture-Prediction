library(data.table)
library(lubridate)
library(stringr)

# ============================================================
# ADD VITAL SIGNS TO NON-BC ADMISSIONS
# ============================================================
# Equivalent to BC_5. Matches vitals in a +-24h window around
# AdmissionDate. Computes baseline (closest) values, aggregated
# flags (fever, hypotension, etc.), and NEWS score.
# ============================================================

nbc_df <- fread("Output Files/CSVs/NonBC_Admissions/nbc_with_infection_source.csv",
                colClasses = c(ClusterID = "character"))
nbc_df[, AdmissionDate := as.POSIXct(AdmissionDate, tz = "UTC")]

cohort_ids <- unique(nbc_df$ClusterID)

# ============================================================
# 1. LOAD AND PREPARE VITAL SIGNS DATA
# ============================================================

vs1 <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_VitalSigns1.csv",
             select = c("ClusterID", "PerformedDateTime", "EventName", "EventResult"),
             colClasses = c(ClusterID = "character"))
vs2 <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_VitalSigns2.csv",
             select = c("ClusterID", "PerformedDateTime", "EventName", "EventResult"),
             colClasses = c(ClusterID = "character"))
o2d <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_OxygenDelivery.csv",
             select = c("ClusterID", "PerformedDateTime", "EventName", "EventResult"),
             colClasses = c(ClusterID = "character"))

vitals <- rbindlist(list(vs1, vs2, o2d))
rm(vs1, vs2, o2d)
vitals <- vitals[ClusterID %in% cohort_ids]
vitals[, PerformedDateTime := suppressWarnings(ymd_hms(PerformedDateTime))]
vitals <- vitals[!is.na(PerformedDateTime)]

vitals <- vitals[EventName %in% c("Heart Rate", "Respiratory Rate", "Oxygen Saturation",
                                   "Oxygen L/min delivered", "Systolic Blood Pressure",
                                   "Diastolic Blood Pressure", "AVPU Scale", "Temperature Tympanic")]

vitals[, EventName_Short := fcase(
  EventName == "Heart Rate",             "HR",
  EventName == "Respiratory Rate",       "RR",
  EventName == "Temperature Tympanic",   "Temp",
  EventName == "Oxygen Saturation",      "O2Sat",
  EventName == "Oxygen L/min delivered", "O2Flow",
  EventName == "Systolic Blood Pressure","SBP",
  EventName == "Diastolic Blood Pressure","DBP",
  EventName == "AVPU Scale",             "AVPU",
  default = "Other"
)]

# ============================================================
# 2. MATCH VITALS TO ADMISSIONS (+-24H WINDOW)
# ============================================================

nbc_vs <- nbc_df[, .(AdmissionID, ClusterID, JoinTime = as.numeric(AdmissionDate))]
vitals[, PerformedNumeric := as.numeric(PerformedDateTime)]
vitals[, ActualVitalTime  := PerformedNumeric]

vital_types <- unique(vitals$EventName_Short)
nbc_vs_exp  <- nbc_vs[, .(EventName_Short = vital_types), by = .(AdmissionID, ClusterID, JoinTime)]
nbc_vs_exp[, `:=`(window_start = JoinTime - 24*3600, window_end = JoinTime + 1*3600)]

vs_join <- vitals[, .(ClusterID, EventName_Short, PerformedNumeric, ActualVitalTime, EventResult)]

matched_vs <- vs_join[nbc_vs_exp,
  on = .(ClusterID, EventName_Short, PerformedNumeric >= window_start, PerformedNumeric <= window_end),
  nomatch = 0L, allow.cartesian = TRUE]
matched_vs[, HoursFromAdmission := (ActualVitalTime - JoinTime) / 3600]

# ============================================================
# 3. AGGREGATED FLAGS (MIN/MAX/MEAN AND CLINICAL FLAGS)
# ============================================================

closest_vs <- matched_vs[, .SD[which.min(abs(HoursFromAdmission))],
                          by = .(AdmissionID, EventName_Short)]

all_vs_wide <- dcast(matched_vs, AdmissionID + ActualVitalTime ~ EventName_Short,
                     value.var = "EventResult", fun.aggregate = function(x) x[1])
for (col in c("HR","RR","Temp","O2Sat","O2Flow","SBP","DBP")) {
  if (col %in% names(all_vs_wide)) all_vs_wide[, (col) := as.numeric(get(col))]
}
if ("AVPU"   %in% names(all_vs_wide)) all_vs_wide[, AVPU   := fifelse(str_detect(AVPU, "^A|^Alert"), 0L, 1L)]
if ("O2Flow" %in% names(all_vs_wide)) all_vs_wide[, O2orAir := fifelse(!is.na(O2Flow) & O2Flow > 0, 1L, 0L)]

# Compute per-row NEWS for each timepoint (matches BC_5 approach)
if ("RR"      %in% names(all_vs_wide)) all_vs_wide[, RR_score    := fcase(is.na(RR),    NA_integer_, RR    <= 8,  3L, RR    <= 11, 1L, RR    <= 20, 0L, RR    <= 24, 2L, default = 3L)]
if ("O2Sat"   %in% names(all_vs_wide)) all_vs_wide[, O2Sat_score := fcase(is.na(O2Sat), NA_integer_, O2Sat <= 91, 3L, O2Sat <= 93, 2L, O2Sat <= 95, 1L, default = 0L)]
if ("O2orAir" %in% names(all_vs_wide)) all_vs_wide[, O2orAir_score := fcase(is.na(O2orAir), NA_integer_, O2orAir == 1, 2L, default = 0L)]
if ("SBP"     %in% names(all_vs_wide)) all_vs_wide[, SBP_score   := fcase(is.na(SBP),  NA_integer_, SBP   <= 90, 3L, SBP   <= 100, 2L, SBP   <= 110, 1L, SBP <= 219, 0L, default = 3L)]
if ("HR"      %in% names(all_vs_wide)) all_vs_wide[, HR_score    := fcase(is.na(HR),   NA_integer_, HR    <= 40, 3L, HR    <= 50,  1L, HR    <= 90,  0L, HR  <= 110, 1L, HR <= 130, 2L, default = 3L)]
if ("AVPU"    %in% names(all_vs_wide)) all_vs_wide[, AVPU_score  := fcase(is.na(AVPU), 0L, AVPU == 0, 0L, default = 3L)]
if ("Temp"    %in% names(all_vs_wide)) all_vs_wide[, Temp_score  := fcase(is.na(Temp), NA_integer_, Temp  <= 35, 3L, Temp  <= 36,  1L, Temp  <= 38,  0L, Temp <= 39, 1L, default = 2L)]
news_score_cols <- intersect(c("RR_score","O2Sat_score","O2orAir_score","SBP_score","HR_score","AVPU_score","Temp_score"), names(all_vs_wide))
all_vs_wide[, NEWS_row := rowSums(.SD, na.rm = TRUE), .SDcols = news_score_cols]
all_vs_wide[, n_scores_row := rowSums(!is.na(.SD)), .SDcols = news_score_cols]
all_vs_wide[n_scores_row < 3, NEWS_row := NA_integer_]
all_vs_wide[, c(news_score_cols, "n_scores_row") := NULL]

vs_agg <- all_vs_wide[, .(
  VS_NumMeasurements = .N,
  VS_HR_Max  = max(HR,    na.rm = TRUE), VS_HR_Min  = min(HR,    na.rm = TRUE), VS_HR_Mean  = mean(HR,    na.rm = TRUE),
  VS_RR_Max  = max(RR,    na.rm = TRUE), VS_RR_Min  = min(RR,    na.rm = TRUE), VS_RR_Mean  = mean(RR,    na.rm = TRUE),
  VS_Temp_Max = max(Temp, na.rm = TRUE), VS_Temp_Min = min(Temp, na.rm = TRUE), VS_Temp_Mean = mean(Temp, na.rm = TRUE),
  VS_O2Sat_Max = max(O2Sat, na.rm = TRUE), VS_O2Sat_Min = min(O2Sat, na.rm = TRUE), VS_O2Sat_Mean = mean(O2Sat, na.rm = TRUE),
  VS_SBP_Max = max(SBP,  na.rm = TRUE), VS_SBP_Min = min(SBP,  na.rm = TRUE), VS_SBP_Mean = mean(SBP,  na.rm = TRUE),
  VS_DBP_Max = max(DBP,  na.rm = TRUE), VS_DBP_Min = min(DBP,  na.rm = TRUE), VS_DBP_Mean = mean(DBP,  na.rm = TRUE),
  VS_NEWS_Max  = max(NEWS_row,  na.rm = TRUE),
  VS_NEWS_Min  = min(NEWS_row,  na.rm = TRUE),
  VS_NEWS_Mean = mean(NEWS_row, na.rm = TRUE),
  VS_AnyO2Use       = as.integer(any(O2orAir == 1, na.rm = TRUE)),
  VS_AnyFever_38    = as.integer(any(Temp >= 38.0, na.rm = TRUE)),
  VS_AnyFever_38_3  = as.integer(any(Temp >= 38.3, na.rm = TRUE)),
  VS_AnyHypothermia = as.integer(any(Temp < 36.0,  na.rm = TRUE)),
  VS_AnyHypotension = as.integer(any(SBP < 90,     na.rm = TRUE)),
  VS_AnyTachycardia = as.integer(any(HR > 100,      na.rm = TRUE)),
  VS_AnyTachypnea   = as.integer(any(RR > 20,       na.rm = TRUE)),
  VS_AnyHypoxia     = as.integer(any(O2Sat < 95,    na.rm = TRUE))
), by = AdmissionID]

agg_num_cols <- names(vs_agg)[sapply(vs_agg, is.numeric)]
for (col in agg_num_cols) vs_agg[is.infinite(get(col)) | is.nan(get(col)), (col) := NA]

# ============================================================
# 4. BASELINE VITALS AND NEWS SCORE
# ============================================================

# Process text-valued vitals before numeric coercion (AVPU text must be matched first)
closest_vs[EventName_Short == "AVPU",   EventResult := fifelse(grepl("^A|^Alert", EventResult, ignore.case = TRUE), "0", "1")]
closest_vs[EventName_Short == "O2Flow", EventResult := fifelse(!is.na(suppressWarnings(as.numeric(EventResult))) & as.numeric(EventResult) > 0, "1", "0")]
closest_vs[, EventResult := as.numeric(EventResult)]

baseline_vs     <- dcast(closest_vs, AdmissionID ~ EventName_Short, value.var = "EventResult", fun.aggregate = function(x) x[1])
baseline_vs_hrs <- dcast(closest_vs, AdmissionID ~ EventName_Short, value.var = "HoursFromAdmission", fun.aggregate = function(x) x[1])

for (col in c("HR","RR","Temp","O2Sat","SBP","DBP","AVPU","O2Flow")) {
  if (!col %in% names(baseline_vs)) baseline_vs[, (col) := NA_real_]
}
baseline_vs[, O2orAir := fifelse(!is.na(O2Flow) & O2Flow > 0, 1L, 0L)]

baseline_vs[, `:=`(
  RR_score    = fcase(is.na(RR),    NA_integer_, RR    <= 8,  3L, RR    <= 11, 1L, RR    <= 20, 0L, RR    <= 24, 2L, default = 3L),
  O2Sat_score = fcase(is.na(O2Sat), NA_integer_, O2Sat <= 91, 3L, O2Sat <= 93, 2L, O2Sat <= 95, 1L, default = 0L),
  O2orAir_score = fcase(is.na(O2orAir), NA_integer_, O2orAir == 1, 2L, default = 0L),
  SBP_score   = fcase(is.na(SBP),  NA_integer_, SBP   <= 90, 3L, SBP   <= 100, 2L, SBP   <= 110, 1L, SBP <= 219, 0L, default = 3L),
  HR_score    = fcase(is.na(HR),   NA_integer_, HR    <= 40, 3L, HR    <= 50,  1L, HR    <= 90,  0L, HR  <= 110, 1L, HR <= 130, 2L, default = 3L),
  AVPU_score  = fcase(is.na(AVPU), 0L, AVPU == 0, 0L, default = 3L),
  Temp_score  = fcase(is.na(Temp), NA_integer_, Temp  <= 35, 3L, Temp  <= 36,  1L, Temp  <= 38,  0L, Temp <= 39, 1L, default = 2L)
)]
score_cols <- c("RR_score","O2Sat_score","O2orAir_score","SBP_score","HR_score","AVPU_score","Temp_score")
baseline_vs[, VS_Baseline_NEWS := rowSums(.SD, na.rm = TRUE), .SDcols = score_cols]
baseline_vs[, n_scores := rowSums(!is.na(.SD)), .SDcols = score_cols]
baseline_vs[n_scores < 3, VS_Baseline_NEWS := NA_integer_]

setnames(baseline_vs,
  old = c("HR","RR","Temp","O2Sat","SBP","DBP","AVPU","O2orAir"),
  new = c("VS_Baseline_HR","VS_Baseline_RR","VS_Baseline_Temp","VS_Baseline_O2Sat",
          "VS_Baseline_SBP","VS_Baseline_DBP","VS_Baseline_AVPU","VS_Baseline_O2orAir"),
  skip_absent = TRUE)

if ("HR" %in% names(baseline_vs_hrs)) {
  baseline_vs[, VS_Baseline_HoursFromAdmission := baseline_vs_hrs$HR]
} else {
  baseline_vs[, VS_Baseline_HoursFromAdmission := NA_real_]
}

baseline_vs_out <- baseline_vs[, .(AdmissionID,
  VS_Baseline_HR, VS_Baseline_RR, VS_Baseline_Temp, VS_Baseline_O2Sat,
  VS_Baseline_SBP, VS_Baseline_DBP, VS_Baseline_AVPU, VS_Baseline_O2orAir,
  VS_Baseline_HoursFromAdmission, VS_Baseline_NEWS)]

nbc_df <- merge(nbc_df, baseline_vs_out, by = "AdmissionID", all.x = TRUE)
nbc_df <- merge(nbc_df, vs_agg,          by = "AdmissionID", all.x = TRUE)

nbc_df[, `:=`(
  VS_DataAvailable = !is.na(VS_NumMeasurements) & VS_NumMeasurements > 0,
  VS_NEWSAvailable = !is.na(VS_Baseline_NEWS),
  VS_SIRS_Count    = as.integer(!is.na(VS_Baseline_Temp) & (VS_Baseline_Temp < 36 | VS_Baseline_Temp > 38)) +
                     as.integer(!is.na(VS_Baseline_HR)   &  VS_Baseline_HR   > 90) +
                     as.integer(!is.na(VS_Baseline_RR)   &  VS_Baseline_RR   > 20),
  VS_qSOFA_Count   = as.integer(!is.na(VS_Baseline_RR)   &  VS_Baseline_RR   >= 22) +
                     as.integer(!is.na(VS_Baseline_SBP)  &  VS_Baseline_SBP  <= 100) +
                     as.integer(!is.na(VS_Baseline_AVPU) &  VS_Baseline_AVPU > 0)
)]

rm(vitals, matched_vs, closest_vs, all_vs_wide, baseline_vs, baseline_vs_hrs, vs_agg)

# ============================================================
# 5. SAVE OUTPUT
# ============================================================

fwrite(nbc_df, "Output Files/CSVs/NonBC_Admissions/nbc_with_vitals.csv")
