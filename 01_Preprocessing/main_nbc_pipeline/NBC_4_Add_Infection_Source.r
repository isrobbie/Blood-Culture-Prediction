library(data.table)
library(lubridate)
library(stringr)

# ============================================================
# ADD INFECTION SOURCE FROM ANTIBIOTICS TO NON-BC ADMISSIONS
# ============================================================
# Equivalent to BC_4. Matches EPR and ITU antibiotic prescriptions
# in the 14 days before admission (+1h). Uses BERT-classified
# indications for infection source, and classifies antibiotic
# classes. Window anchor: AdmissionDate.
# ============================================================

nbc_df <- fread("Output Files/CSVs/NonBC_Admissions/nbc_with_prior_bc.csv",
                colClasses = c(ClusterID = "character"))
nbc_df[, AdmissionDate := as.POSIXct(AdmissionDate, tz = "UTC")]

cohort_ids <- unique(nbc_df$ClusterID)

# ============================================================
# 1. LOAD AND PREPARE ANTIBIOTIC DATA
# ============================================================

# Using updated add_fields file which includes ScheduledStart/StopDateTime for TTA duration (mirrors BC_4)
epr_abx <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_EPRPrescriptions_add_fields/IORD_BloodCultPatients_57_20251030_EPRPrescriptionsAntibiotics_add_fields.csv",
                 select = c("ClusterID", "PrescriptionDateTime", "PrescriptionType", "Drug", "Route",
                            "Indication", "ScheduledStartDateTime", "ScheduledStopDateTime",
                            "FirstAdministrationDateTime", "LastAdministrationDateTime",
                            "NumberOfDosesAdministered"),
                 colClasses = c(ClusterID = "character"))
epr_abx <- epr_abx[ClusterID %in% cohort_ids]
epr_abx[, `:=`(
  PrescriptionDateTime        = suppressWarnings(ymd_hms(PrescriptionDateTime)),
  ScheduledStartDateTime      = suppressWarnings(ymd_hms(ScheduledStartDateTime)),
  ScheduledStopDateTime       = suppressWarnings(ymd_hms(ScheduledStopDateTime)),
  FirstAdministrationDateTime = suppressWarnings(ymd_hms(FirstAdministrationDateTime)),
  LastAdministrationDateTime  = suppressWarnings(ymd_hms(LastAdministrationDateTime)),
  NumberOfDosesAdministered   = suppressWarnings(as.numeric(NumberOfDosesAdministered)),
  Indication = tolower(trimws(Indication)),
  IsTTO = (PrescriptionType == "TTO"),
  Source = "EPR"
)]
# For TTAs without administration dates: use scheduled dates so they have a real duration
# (mirrors BC_4: patients take TTAs at home so admin dates absent; long courses valid, cap at 365 days)
epr_abx[IsTTO == TRUE & is.na(FirstAdministrationDateTime) & !is.na(ScheduledStartDateTime), `:=`(
  FirstAdministrationDateTime = ScheduledStartDateTime,
  LastAdministrationDateTime  = fifelse(
    !is.na(ScheduledStopDateTime) & as.numeric(difftime(ScheduledStopDateTime, ScheduledStartDateTime, units="days")) <= 365,
    ScheduledStopDateTime,
    ScheduledStartDateTime + days(365)
  )
)]
# Remaining prescriptions without administration dates: point-in-time at PrescriptionDateTime
epr_abx[is.na(FirstAdministrationDateTime), `:=`(
  FirstAdministrationDateTime = PrescriptionDateTime,
  LastAdministrationDateTime  = PrescriptionDateTime
)]

itu_abx <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_ITUPrescriptions.csv",
                 colClasses = c(ClusterID = "character"))
itu_abx <- itu_abx[ClusterID %in% cohort_ids]
itu_abx[, `:=`(
  FirstAdministrationDateTime = suppressWarnings(ymd_hms(FirstAdministrationDateTime)),
  LastAdministrationDateTime  = suppressWarnings(ymd_hms(LastAdministrationDateTime)),
  NumberOfDosesAdministered   = suppressWarnings(as.numeric(NumberOfDosesAdministered)),
  Indication = NA_character_,
  IsTTO = FALSE,
  Source = "ITU"
)]
itu_abx <- itu_abx[!is.na(FirstAdministrationDateTime),
                    .(ClusterID, Drug, Route, Indication, Source, IsTTO,
                      NumberOfDosesAdministered,
                      FirstAdministrationDateTime, LastAdministrationDateTime)]

# ============================================================
# 2. JOIN BERT INFECTION SOURCE CLASSIFICATIONS
# ============================================================

bert_class <- fread("Input Files/bert_classified_indications_update.csv")
setDT(bert_class)
bert_class[, indication_lower := tolower(Indication)]

epr_abx[!is.na(Indication), Indication := gsub("\\r\\n|\\n|\\r|\\s+", " ", Indication)]
epr_abx[!is.na(Indication), Indication := trimws(Indication)]
epr_abx[Indication %in% c("", "null", "na", "n/a", "none"), Indication := NA_character_]

epr_abx <- bert_class[epr_abx, on = c(indication_lower = "Indication")]

bert_cols <- c("urinary", "respiratory", "abdominal", "neurological",
               "skin_soft_tissue", "ent", "orthopaedic", "other_specific", "no_specific_source")

all_abx <- rbind(
  epr_abx[, c(list(ClusterID = ClusterID, Drug = Drug, Route = Route, Source = Source,
                   IsTTO = IsTTO, NumberOfDosesAdministered = NumberOfDosesAdministered,
                   PrescriptionDateTime = PrescriptionDateTime,
                   FirstAdministrationDateTime = FirstAdministrationDateTime,
                   LastAdministrationDateTime  = LastAdministrationDateTime),
               lapply(bert_cols, function(col) if (col %in% names(epr_abx)) get(col) else NA_integer_)) |>
             setNames(c("ClusterID","Drug","Route","Source","IsTTO","NumberOfDosesAdministered",
                        "PrescriptionDateTime",
                        "FirstAdministrationDateTime","LastAdministrationDateTime", bert_cols))],
  itu_abx[, c(list(ClusterID = ClusterID, Drug = Drug, Route = Route, Source = Source,
                   IsTTO = IsTTO, NumberOfDosesAdministered = NumberOfDosesAdministered,
                   PrescriptionDateTime = as.POSIXct(NA_real_, origin = "1970-01-01"),
                   FirstAdministrationDateTime = FirstAdministrationDateTime,
                   LastAdministrationDateTime  = LastAdministrationDateTime),
               setNames(lapply(bert_cols, function(col) NA_integer_), bert_cols))],
  fill = TRUE)

# ============================================================
# 3. CLASSIFY ANTIBIOTIC DRUGS
# ============================================================

# Normalize drug names
all_abx[, Drug_lower := tolower(trimws(Drug))]
all_abx[, Drug_lower := gsub("\\s*\\d+(\\.\\d+)?\\s*(mg|g|ml|mcg|units?|iu|tablets?|capsules?|doses?)\\s*/?\\s*\\d*(\\.\\d+)?\\s*(mg|g|ml|mcg|units?|iu|tablets?|capsules?|doses?)?", "", Drug_lower, ignore.case = TRUE)]
all_abx[, Drug_lower := gsub("\\([^)]*\\)", "", Drug_lower)]
all_abx[, Drug_lower := gsub("[+\\-\\/]", " ", Drug_lower)]
all_abx[, Drug_lower := gsub("\\s*\\([^)]*$", "", Drug_lower)]
all_abx[, Drug_lower := gsub("\\s*(injection|suspension|topical|oral|iv|im|sc|infusion|tablet|tablets|capsule|capsules|cream|ointment|gel|solution|powder|for injection|concentrate|pcc)\\s*$", "", Drug_lower, ignore.case = TRUE)]
all_abx[, Drug_lower := trimws(gsub("\\s+", " ", Drug_lower))]
all_abx[Drug_lower %like% "^augmentin",                   Drug_lower := "co-amoxiclav"]
all_abx[Drug_lower %like% "^tazocin" | Drug_lower %like% "^zocin", Drug_lower := "piperacillin tazobactam"]
all_abx[Drug_lower %like% "^zavicefta",                   Drug_lower := "ceftazidime avibactam"]
all_abx[Drug_lower %like% "^zyvox",                       Drug_lower := "linezolid"]
all_abx[Drug_lower %like% "^vancocin",                    Drug_lower := "vancomycin"]
all_abx[Drug_lower %like% "^flagyl",                      Drug_lower := "metronidazole"]
all_abx[Drug_lower %like% "^septrim" | Drug_lower %like% "^bactrim", Drug_lower := "co-trimoxazole"]
all_abx[Drug_lower %like% "amoxicillin.*clavulan",        Drug_lower := "co-amoxiclav"]
all_abx[Drug_lower %like% "trimethoprim.*sulfa",          Drug_lower := "co-trimoxazole"]
all_abx[Drug_lower %like% "piperacillin.*tazobactam",     Drug_lower := "piperacillin tazobactam"]

# Drug class patterns
broad_spectrum_pattern      <- paste(c("piperacillin","tazocin","co-amoxiclav","ceftriaxone","ceftazidime","cefepime","cefotaxime","meropenem","imipenem","ertapenem","ciprofloxacin","levofloxacin","moxifloxacin","aztreonam","temocillin","tigecycline","linezolid","daptomycin","colistin","cefiderocol","avibactam","vaborbactam","relebactam"), collapse = "|")
beta_lactamase_combo_pattern <- paste(c("co-amoxiclav","piperacillin.*tazobactam","tazocin","tazobactam"), collapse = "|")
narrow_penicillin_pattern   <- paste(c("amoxicillin","flucloxacillin","benzylpenicillin","phenoxymethylpenicillin","penicillin","pivmecillinam","procaine"), collapse = "|")
ceph_1st_2nd_pattern        <- paste(c("cefazolin","cefuroxime","cefoxitin","cefradine","cefixime","cefalexin"), collapse = "|")
ceph_3rd_plus_pattern       <- paste(c("ceftriaxone","ceftazidime","cefotaxime","cefepime","ceftaroline","ceftolozane","cefiderocol"), collapse = "|")
carbapenem_pattern          <- paste(c("meropenem","imipenem","ertapenem"), collapse = "|")
fluoroquinolone_pattern     <- paste(c("ciprofloxacin","levofloxacin","moxifloxacin","ofloxacin","delafloxacin"), collapse = "|")
aminoglycoside_pattern      <- paste(c("gentamicin","amikacin","tobramycin","neomycin"), collapse = "|")
glycopeptide_pattern        <- paste(c("vancomycin","teicoplanin","dalbavancin"), collapse = "|")
macrolide_pattern           <- paste(c("erythromycin","clarithromycin","azithromycin"), collapse = "|")
tetracycline_pattern        <- paste(c("doxycycline","tetracycline","tigecycline","minocycline","lymecycline","oxytetracycline"), collapse = "|")
metronidazole_pattern       <- "metronidazole"
nitrofurantoin_pattern      <- "nitrofurantoin"
co_trimoxazole_pattern      <- "co-trimoxazole"

all_antibiotic_patterns <- paste(c(
  "co-amoxiclav","piperacillin","tazocin","amoxicillin","flucloxacillin","benzylpenicillin",
  "phenoxymethylpenicillin","penicillin","pivmecillinam",
  "cefazolin","cefuroxime","cefoxitin","cefradine","cefixime","cefalexin",
  "ceftriaxone","ceftazidime","cefotaxime","cefepime","ceftaroline","ceftolozane","cefiderocol",
  "meropenem","imipenem","ertapenem","ciprofloxacin","levofloxacin","moxifloxacin","ofloxacin",
  "gentamicin","amikacin","tobramycin","vancomycin","teicoplanin","dalbavancin",
  "erythromycin","clarithromycin","azithromycin","doxycycline","tetracycline","tigecycline",
  "minocycline","lymecycline","metronidazole","nitrofurantoin","co-trimoxazole",
  "clindamycin","rifampicin","rifampin","fosfomycin","aztreonam","colistin","fidaxomicin",
  "daptomycin","linezolid","temocillin","trimethoprim",
  "isoniazid","ethambutol","pyrazinamide","dapsone","pristinamycin",
  "polymyxin","quinupristin","dalfopristin"
), collapse = "|")

antiviral_antifungal_pattern <- paste(c(
  # Antivirals
  "aciclovir","valaciclovir","valganciclovir","ganciclovir","oseltamivir","remdesivir",
  "foscarnet","zanamivir","cidofovir","famciclovir","maribavir",
  # Antifungals
  "fluconazole","voriconazole","posaconazole","itraconazole","isavuconazole",
  "micafungin","caspofungin","anidulafungin","amphotericin","terbinafine",
  "flucytosine","griseofulvin","rezafungin",
  # Topical antibiotics, rarely used systemics, and urinary antiseptics
  "fusidic","fusidate","mupirocin","chloramphenicol","methenamine",
  # Antiparasitics and antiprotozoals
  "pentamidine","ivermectin","mebendazole","albendazole","praziquantel",
  "nitazoxanide","stibogluconate","tinidazole",
  # Antimalarials
  "artemether","lumefantrine","artesunate","atovaquone","proguanil",
  "chloroquine","mefloquine","primaquine","pyrimethamine","mepacrine",
  # Antiretrovirals
  "abacavir","lamivudine","tenofovir","emtricitabine","efavirenz",
  # Non-antimicrobials
  "dalteparin","enoxaparin","heparin"
), collapse = "|")

all_abx <- all_abx[
  (Source == "ITU" & grepl(all_antibiotic_patterns, Drug_lower, ignore.case = TRUE)) |
  (Source == "EPR" & !grepl(antiviral_antifungal_pattern, Drug_lower, ignore.case = TRUE))
]

all_abx[, `:=`(
  IsBroadSpectrum      = as.integer(grepl(broad_spectrum_pattern,      Drug_lower)),
  IsBetaLactamaseCombo = as.integer(grepl(beta_lactamase_combo_pattern, Drug_lower)),
  IsNarrowPenicillin   = as.integer(grepl(narrow_penicillin_pattern,   Drug_lower)),
  IsCeph1st2ndGen      = as.integer(grepl(ceph_1st_2nd_pattern,        Drug_lower)),
  IsCeph3rdGenPlus     = as.integer(grepl(ceph_3rd_plus_pattern,       Drug_lower)),
  IsCarbapenem         = as.integer(grepl(carbapenem_pattern,           Drug_lower)),
  IsFluoroquinolone    = as.integer(grepl(fluoroquinolone_pattern,     Drug_lower)),
  IsAminoglycoside     = as.integer(grepl(aminoglycoside_pattern,      Drug_lower)),
  IsGlycopeptide       = as.integer(grepl(glycopeptide_pattern,        Drug_lower)),
  IsMacrolide          = as.integer(grepl(macrolide_pattern,           Drug_lower)),
  IsTetracycline       = as.integer(grepl(tetracycline_pattern,        Drug_lower)),
  IsMetronidazole      = as.integer(grepl(metronidazole_pattern,       Drug_lower)),
  IsNitrofurantoin     = as.integer(grepl(nitrofurantoin_pattern,      Drug_lower)),
  IsCoTrimoxazole      = as.integer(grepl(co_trimoxazole_pattern,      Drug_lower))
)]

all_abx[, AntibioticCategory := fcase(
  IsBetaLactamaseCombo == 1, "Beta-lactamase combo",
  IsCarbapenem         == 1, "Carbapenem",
  IsGlycopeptide       == 1, "Glycopeptide",
  IsCeph3rdGenPlus     == 1, "3rd gen+ cephalosporin",
  IsCeph1st2ndGen      == 1, "1st/2nd gen cephalosporin",
  IsFluoroquinolone    == 1, "Fluoroquinolone",
  IsAminoglycoside     == 1, "Aminoglycoside",
  IsNarrowPenicillin   == 1, "Narrow penicillin",
  IsMacrolide          == 1, "Macrolide",
  IsTetracycline       == 1, "Tetracycline",
  IsMetronidazole      == 1, "Metronidazole",
  IsNitrofurantoin     == 1, "Nitrofurantoin",
  IsCoTrimoxazole      == 1, "Co-trimoxazole",
  default = "Other"
)]

all_abx[, Route_lower := tolower(Route)]
all_abx[, RouteCategory := fcase(
  grepl("iv|intravenous|infusion", Route_lower, ignore.case = TRUE), "IV",
  grepl("oral|po|enteral|nasogastric|gastrostomy|jejunostomy|peg|gastrojejunostomy", Route_lower, ignore.case = TRUE), "Oral",
  default = "Other"
)]

# ============================================================
# 4. MATCH ANTIBIOTICS TO ADMISSIONS (14 DAYS BEFORE TO +1H)
# ============================================================

nbc_windows <- nbc_df[, .(AdmissionID, ClusterID, AdmissionDate,
                           window_start = AdmissionDate - days(14),
                           window_end   = AdmissionDate + hours(1))]

min_admission_date <- min(nbc_df$AdmissionDate, na.rm = TRUE)
max_admission_date <- max(nbc_df$AdmissionDate, na.rm = TRUE)
all_abx <- all_abx[FirstAdministrationDateTime >= (min_admission_date - days(14)) &
                   FirstAdministrationDateTime <= (max_admission_date + hours(1))]

all_abx[, `:=`(OrigFirst = FirstAdministrationDateTime, OrigLast = LastAdministrationDateTime)]
setkey(all_abx, ClusterID, FirstAdministrationDateTime, LastAdministrationDateTime)
setkey(nbc_windows, ClusterID, window_start, window_end)

abx_matched <- all_abx[nbc_windows,
  on = .(ClusterID, FirstAdministrationDateTime <= window_end, LastAdministrationDateTime >= window_start),
  nomatch = 0L, allow.cartesian = TRUE]
abx_matched[, `:=`(FirstAdministrationDateTime = OrigFirst, LastAdministrationDateTime = OrigLast)]

# Leakage guard: TTAs use scheduled-start as their effective start, but a discharge order
# written at/after admission (with a back-dated scheduled start) was not known at admission
# time, so it must not count as "before" (mirrors BC_4).
abx_matched[, TTO_OrderedAfterAdmission := IsTTO == TRUE & !is.na(PrescriptionDateTime) &
                                           PrescriptionDateTime >= AdmissionDate]
abx_matched[, `:=`(
  IsAdminBefore = LastAdministrationDateTime  > (AdmissionDate - days(14)) &
                  FirstAdministrationDateTime  < AdmissionDate &
                  !TTO_OrderedAfterAdmission,
  IsAdmin_24h   = LastAdministrationDateTime  >= (AdmissionDate - hours(24)) &
                  FirstAdministrationDateTime  < AdmissionDate &
                  !TTO_OrderedAfterAdmission,
  IsAdmin_1wk   = LastAdministrationDateTime  >= (AdmissionDate - days(7)) &
                  FirstAdministrationDateTime  < AdmissionDate &
                  !TTO_OrderedAfterAdmission,
  IsAdmin_2wk   = LastAdministrationDateTime  > (AdmissionDate - days(14)) &
                  FirstAdministrationDateTime  < AdmissionDate &
                  !TTO_OrderedAfterAdmission
)]

# ============================================================
# 5. AGGREGATE PER ADMISSION
# ============================================================

abx_agg <- abx_matched[IsAdminBefore == TRUE, .(
  Source_Urinary          = as.integer(any(urinary == 1,          na.rm = TRUE)),
  Source_Respiratory      = as.integer(any(respiratory == 1,      na.rm = TRUE)),
  Source_Abdominal        = as.integer(any(abdominal == 1,        na.rm = TRUE)),
  Source_CNS              = as.integer(any(neurological == 1,     na.rm = TRUE)),
  Source_Skin_Soft_Tissue = as.integer(any(skin_soft_tissue == 1, na.rm = TRUE)),
  Source_ENT              = as.integer(any(ent == 1,              na.rm = TRUE)),
  Source_Orthopaedic      = as.integer(any(orthopaedic == 1,      na.rm = TRUE)),
  Source_Other            = as.integer(any(other_specific == 1,   na.rm = TRUE)),
  Source_Unspecific       = as.integer(any(no_specific_source == 1, na.rm = TRUE)),
  HasAntibiotics_Before      = 1L,
  BroadSpectrum_Before       = as.integer(any(IsBroadSpectrum      == 1)),
  BetaLactamaseCombo_Before  = as.integer(any(IsBetaLactamaseCombo == 1)),
  NarrowPenicillin_Before    = as.integer(any(IsNarrowPenicillin   == 1)),
  Ceph1st2ndGen_Before       = as.integer(any(IsCeph1st2ndGen      == 1)),
  Ceph3rdGenPlus_Before      = as.integer(any(IsCeph3rdGenPlus     == 1)),
  Carbapenem_Before          = as.integer(any(IsCarbapenem         == 1)),
  Fluoroquinolone_Before     = as.integer(any(IsFluoroquinolone    == 1)),
  Aminoglycoside_Before      = as.integer(any(IsAminoglycoside     == 1)),
  Glycopeptide_Before        = as.integer(any(IsGlycopeptide       == 1)),
  Macrolide_Before           = as.integer(any(IsMacrolide          == 1)),
  Tetracycline_Before        = as.integer(any(IsTetracycline       == 1)),
  Metronidazole_Before       = as.integer(any(IsMetronidazole      == 1)),
  Nitrofurantoin_Before      = as.integer(any(IsNitrofurantoin     == 1)),
  CoTrimoxazole_Before       = as.integer(any(IsCoTrimoxazole      == 1)),
  Other_Before               = as.integer(any(AntibioticCategory == "Other")),
  HasAntibiotics_24h_Before         = as.integer(any(IsAdmin_24h)),
  NumDistinctAntibiotics_24h_Before = uniqueN(Drug_lower[IsAdmin_24h]),
  HasAntibiotics_1wk_Before         = as.integer(any(IsAdmin_1wk)),
  NumDistinctAntibiotics_1wk_Before = uniqueN(Drug_lower[IsAdmin_1wk]),
  HasAntibiotics_2wk_Before         = as.integer(any(IsAdmin_2wk)),
  NumDosesAdministered_WeekBefore   = sum(NumberOfDosesAdministered[IsAdminBefore], na.rm = TRUE),
  DaysContinuousAbx_PreCulture = {
    if (all(!IsAdminBefore)) {
      0L
    } else {
      adm_dt <- first(AdmissionDate)
      abx_before <- data.table(
        Start = FirstAdministrationDateTime[IsAdminBefore],
        End   = LastAdministrationDateTime[IsAdminBefore]
      )
      abx_before <- abx_before[order(Start)]
      continuous_start <- adm_dt
      for (i in nrow(abx_before):1) {
        if (abx_before$End[i] >= continuous_start - hours(48)) {
          continuous_start <- min(continuous_start, abx_before$Start[i])
        }
      }
      as.integer(ceiling(difftime(adm_dt, continuous_start, units = "days")))
    }
  },
  Abx_Route_Category = {
    if (any(RouteCategory == "IV")) {
      "IV"
    } else {
      routes <- RouteCategory[!is.na(RouteCategory)]
      if (length(routes) == 0) NA_character_
      else names(table(routes))[which.max(table(routes))]
    }
  },
  HasTTO_Before = as.integer(any(IsTTO & IsAdminBefore, na.rm = TRUE))
), by = AdmissionID]

specific_src <- c("Source_Urinary","Source_Respiratory","Source_Abdominal","Source_CNS",
                  "Source_Skin_Soft_Tissue","Source_ENT","Source_Orthopaedic","Source_Other")
abx_agg[, n_specific_src := rowSums(as.matrix(.SD), na.rm = TRUE), .SDcols = specific_src]
abx_agg[, MultipleSourcesFlag := as.integer(n_specific_src >= 2)]
abx_agg[n_specific_src > 0, Source_Unspecific := 0L]
abx_agg[, n_specific_src := NULL]

nbc_df <- merge(nbc_df, abx_agg, by = "AdmissionID", all.x = TRUE)

zero_cols <- c("Source_Urinary","Source_Respiratory","Source_Abdominal","Source_CNS",
               "Source_Skin_Soft_Tissue","Source_ENT","Source_Orthopaedic","Source_Other",
               "Source_Unspecific","MultipleSourcesFlag",
               "HasAntibiotics_Before","BroadSpectrum_Before","BetaLactamaseCombo_Before",
               "NarrowPenicillin_Before","Ceph1st2ndGen_Before","Ceph3rdGenPlus_Before",
               "Carbapenem_Before","Fluoroquinolone_Before","Aminoglycoside_Before",
               "Glycopeptide_Before","Macrolide_Before","Tetracycline_Before",
               "Metronidazole_Before","Nitrofurantoin_Before","CoTrimoxazole_Before",
               "Other_Before",               "HasAntibiotics_24h_Before","NumDistinctAntibiotics_24h_Before",
               "HasAntibiotics_1wk_Before","NumDistinctAntibiotics_1wk_Before",
               "HasAntibiotics_2wk_Before","NumDosesAdministered_WeekBefore",
               "DaysContinuousAbx_PreCulture","HasTTO_Before")
for (col in zero_cols) nbc_df[is.na(get(col)), (col) := 0L]

# ============================================================
# 6. SAVE OUTPUT
# ============================================================

fwrite(nbc_df, "Output Files/CSVs/NonBC_Admissions/nbc_with_infection_source.csv")
