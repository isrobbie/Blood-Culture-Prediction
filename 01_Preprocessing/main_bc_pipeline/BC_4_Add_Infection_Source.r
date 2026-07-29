
library(data.table)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(openxlsx)

# Load Excel output formatting helpers
source("01_Preprocessing/Excel_Output_Helper.r")

# ============================================================
# STEP 1: EXTRACT AND CLASSIFY ANTIBIOTIC INDICATIONS
# ============================================================

# Load blood culture cohort ClusterIDs first for filtering
bc_episodes <- fread("Output Files/CSVs/blood_cultures_with_admissions.csv", 
                     select = "ClusterID",
                     colClasses = c(ClusterID = "character"))
cohort_cluster_ids <- unique(bc_episodes$ClusterID)

# Use data.table for big CSV loading, filtered to cohort
# Using updated add_fields file which includes ScheduledStart/StopDateTime for TTA duration
epr_abx <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_EPRPrescriptions_add_fields/IORD_BloodCultPatients_57_20251030_EPRPrescriptionsAntibiotics_add_fields.csv",
                 select = c("ClusterID", "PrescriptionDateTime", "PrescriptionType", "Drug", "Route", "Indication",
                            "ScheduledStartDateTime", "ScheduledStopDateTime",
                            "FirstAdministrationDateTime", "LastAdministrationDateTime", "NumberOfDosesAdministered"),
                 colClasses = c(ClusterID = "character"))
epr_abx <- epr_abx[ClusterID %in% cohort_cluster_ids]
epr_abx[, `:=`(
  PrescriptionDateTime         = as.POSIXct(PrescriptionDateTime,         tz = "UTC", format = "%Y-%m-%d %H:%M:%S"),
  ScheduledStartDateTime       = suppressWarnings(as.POSIXct(ScheduledStartDateTime,       tz = "UTC", format = "%Y-%m-%d %H:%M:%S")),
  ScheduledStopDateTime        = suppressWarnings(as.POSIXct(ScheduledStopDateTime,        tz = "UTC", format = "%Y-%m-%d %H:%M:%S")),
  FirstAdministrationDateTime  = suppressWarnings(as.POSIXct(FirstAdministrationDateTime,  tz = "UTC", format = "%Y-%m-%d %H:%M:%S")),
  LastAdministrationDateTime   = suppressWarnings(as.POSIXct(LastAdministrationDateTime,   tz = "UTC", format = "%Y-%m-%d %H:%M:%S")),
  Indication                   = tolower(trimws(Indication)),
  NumberOfDosesAdministered    = suppressWarnings(as.numeric(NumberOfDosesAdministered)),
  Source = "EPR"
)]

# Set TTO flag based on PrescriptionType (more robust for future data)
epr_abx[, IsTTO := (PrescriptionType == "TTO")]

# For TTAs without administration dates: use scheduled dates so they have a real duration
# (patients take TTAs at home so admin dates are typically absent, but scheduled dates reflect intended course)
# Long courses (>28 days) are clinically valid: PCP/CMV prophylaxis, TB, osteomyelitis, asplenism penicillin etc.
# Only cap at 365 days to exclude the small number of obvious data-entry placeholders (999 days, 30-year entries).
epr_abx[IsTTO == TRUE & is.na(FirstAdministrationDateTime) & !is.na(ScheduledStartDateTime), `:=`(
  FirstAdministrationDateTime = ScheduledStartDateTime,
  LastAdministrationDateTime  = fifelse(
    !is.na(ScheduledStopDateTime) & as.numeric(difftime(ScheduledStopDateTime, ScheduledStartDateTime, units="days")) <= 365,
    ScheduledStopDateTime,
    ScheduledStartDateTime + days(365)
  )
)]

# For all remaining prescriptions without administration dates, use PrescriptionDateTime
# Treat these as point-in-time prescriptions (no duration)
epr_abx[is.na(FirstAdministrationDateTime), `:=`(
  FirstAdministrationDateTime = PrescriptionDateTime,
  LastAdministrationDateTime  = PrescriptionDateTime
)]

itu_abx <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_ITUPrescriptions.csv",
                 colClasses = c(ClusterID = "character"))
itu_abx <- itu_abx[ClusterID %in% cohort_cluster_ids]

# ITU doesn't have Indication or PrescriptionDateTime - set to NA
itu_abx[, `:=`(
  PrescriptionDateTime = NA_POSIXct_,
  PrescriptionType = NA_character_,
  FirstAdministrationDateTime = as.POSIXct(FirstAdministrationDateTime, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"),
  LastAdministrationDateTime = as.POSIXct(LastAdministrationDateTime, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"),
  Indication = NA_character_,
  NumberOfDosesAdministered = as.numeric(NumberOfDosesAdministered),
  Source = "ITU",
  IsTTO = FALSE
)]
itu_abx <- itu_abx[!is.na(FirstAdministrationDateTime), 
                   .(ClusterID, PrescriptionDateTime, PrescriptionType, Drug, Route, Indication, Source,
                     FirstAdministrationDateTime, LastAdministrationDateTime, NumberOfDosesAdministered, IsTTO)]

# Note: ITU file contains all medications, not just antibiotics

# Use data.table for fast string cleaning (only for non-NA indications)
epr_abx[!is.na(Indication), Indication := gsub("\\r\\n|\\n|\\r|\\s+", " ", Indication)]
epr_abx[!is.na(Indication), Indication := trimws(Indication)]
epr_abx[Indication %in% c("", "null", "na", "n/a", "none"), Indication := NA_character_]

# Use data.table for big CSVs
bert_class <- fread("Input Files/bert_classified_indications_update.csv")
setDT(bert_class)
bert_class[, indication_lower := tolower(Indication)]

# Get unique indications from our antibiotic data
unique_indications <- unique(epr_abx$Indication)

# Keep only relevant classifications (inner join on indications we actually have)
classified <- bert_class[indication_lower %in% unique_indications]

# ============================================================
# STEP 2: LINK ANTIBIOTICS TO INDIVIDUAL BLOOD CULTURES (±14 DAYS)
# ============================================================

# Reload full bc_episodes - each row is already a unique blood culture
bc_episodes <- fread("Output Files/CSVs/blood_cultures_with_admissions.csv", colClasses = c(ClusterID = "character"))
setDT(bc_episodes)
bc_episodes[, CollectionDateTime := as.POSIXct(CollectionDateTime, tz = "UTC")]

# Join antibiotic indications with classifications to bring in source columns
epr_abx <- classified[epr_abx, on = c(indication_lower = "Indication")]

# Combine EPR and ITU data (ITU will be filtered to antibiotics later)
all_abx <- rbind(epr_abx, itu_abx, fill = TRUE)

# Define source categories (these columns now exist in epr_abx from the join above)
source_categories <- c("urinary", "respiratory", "abdominal", "neurological", 
                      "skin_soft_tissue", "ent", "orthopaedic", "other_specific")

# ============================================================
# STEP 2A: CLASSIFY ANTIBIOTICS - NEW CLASSIFICATION SYSTEM
# ============================================================

# Beta-lactamase inhibitor combinations 
beta_lactamase_combos <- c("co-amoxiclav", "piperacillin.*tazobactam", "tazocin", "tazobactam")

# Narrow penicillins 
narrow_penicillins <- c("amoxicillin", "flucloxacillin", "benzylpenicillin", 
                       "phenoxymethylpenicillin", "penicillin", "pivmecillinam", "procaine")

# 1st and 2nd generation cephalosporins 
ceph_1st_2nd_gen <- c("cefazolin", "cefuroxime", "cefoxitin", "cefradine", "cefixime", "cefalexin")

# 3rd generation and beyond cephalosporins 
ceph_3rd_gen_plus <- c("ceftriaxone", "ceftazidime", "cefotaxime", "cefepime", 
                      "ceftaroline", "ceftolozane", "cefiderocol")

# Carbapenems 
carbapenems <- c("meropenem", "imipenem", "ertapenem")

# Fluoroquinolones 
fluoroquinolones <- c("ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin", "delafloxacin")

# Aminoglycosides 
aminoglycosides <- c("gentamicin", "amikacin", "tobramycin", "neomycin")

# Glycopeptides 
glycopeptides <- c("vancomycin", "teicoplanin", "dalbavancin")

# Macrolides 
macrolides <- c("erythromycin", "clarithromycin", "azithromycin")

# Tetracyclines 
tetracyclines <- c("doxycycline", "tetracycline", "tigecycline", "minocycline", 
                  "lymecycline", "oxytetracycline")

# Metronidazole 
metronidazole <- c("metronidazole")

# Nitrofurantoin
nitrofurantoin <- c("nitrofurantoin")

# Co-trimoxazole 
co_trimoxazole <- c("co-trimoxazole")

# Broad-spectrum antibiotics (kept separately as clinically useful)
broad_spectrum <- c("piperacillin", "tazocin", "co-amoxiclav",
                   "ceftriaxone", "ceftazidime", "cefepime", "cefotaxime",
                   "meropenem", "imipenem", "ertapenem",
                   "ciprofloxacin", "levofloxacin", "moxifloxacin",
                   "aztreonam", "temocillin", "tigecycline", "linezolid", "daptomycin",
                   "colistin", "cefiderocol", "avibactam", "vaborbactam", "relebactam")

# Pre-compile regex patterns for efficiency
beta_lactamase_combo_pattern <- paste(beta_lactamase_combos, collapse = "|")
narrow_penicillin_pattern <- paste(narrow_penicillins, collapse = "|")
ceph_1st_2nd_pattern <- paste(ceph_1st_2nd_gen, collapse = "|")
ceph_3rd_plus_pattern <- paste(ceph_3rd_gen_plus, collapse = "|")
carbapenem_pattern <- paste(carbapenems, collapse = "|")
fluoroquinolone_pattern <- paste(fluoroquinolones, collapse = "|")
aminoglycoside_pattern <- paste(aminoglycosides, collapse = "|")
glycopeptide_pattern <- paste(glycopeptides, collapse = "|")
macrolide_pattern <- paste(macrolides, collapse = "|")
tetracycline_pattern <- paste(tetracyclines, collapse = "|")
metronidazole_pattern <- paste(metronidazole, collapse = "|")
nitrofurantoin_pattern <- paste(nitrofurantoin, collapse = "|")
co_trimoxazole_pattern <- paste(co_trimoxazole, collapse = "|")
broad_spectrum_pattern <- paste(broad_spectrum, collapse = "|")

# Filter ITU prescriptions to only include antibiotics (ITU file contains all medications)
all_antibiotic_patterns <- paste(c(
  beta_lactamase_combos, "piperacillin", narrow_penicillins, ceph_1st_2nd_gen, ceph_3rd_gen_plus,
  carbapenems, fluoroquinolones, aminoglycosides, glycopeptides, macrolides,
  tetracyclines, metronidazole, nitrofurantoin, co_trimoxazole,
  # Other antibiotics (will be classified as "Other" category)
  "clindamycin", "rifampicin", "rifampin", "isoniazid", "ethambutol", "pyrazinamide",
  "fosfomycin", "aztreonam", "colistin", "fidaxomicin",
  "dapsone", "daptomycin", "linezolid", "temocillin",
  "cefiderocol", "avibactam", "vaborbactam", "relebactam",
  "trimethoprim", "polymyxin", "quinupristin", "dalfopristin"
), collapse = "|")

# Pattern to exclude non-antibiotics from EPR
# (antivirals, antifungals, antiparasitics, antimalarials, antiretrovirals, topical antibiotics, and non-antimicrobials)
antiviral_antifungal_pattern <- paste(c(
  # Antivirals
  "aciclovir", "valaciclovir", "valganciclovir", "ganciclovir", "oseltamivir", "remdesivir",
  "foscarnet", "zanamivir", "cidofovir", "famciclovir", "maribavir",
  # Antifungals
  "fluconazole", "voriconazole", "posaconazole", "itraconazole", "isavuconazole",
  "micafungin", "caspofungin", "anidulafungin", "amphotericin", "terbinafine",
  "flucytosine", "griseofulvin", "rezafungin",
  # Topical antibiotics, rarely used systemics, and urinary antiseptics
  "fusidic", "fusidate", "mupirocin", "chloramphenicol", "methenamine",
  # Antiparasitics and antiprotozoals
  "pentamidine", "ivermectin", "mebendazole", "albendazole", "praziquantel",
  "nitazoxanide", "stibogluconate", "tinidazole",
  # Antimalarials
  "artemether", "lumefantrine", "artesunate", "atovaquone", "proguanil",
  "chloroquine", "mefloquine", "primaquine", "pyrimethamine", "mepacrine",
  # Antiretrovirals
  "abacavir", "lamivudine", "tenofovir", "emtricitabine", "efavirenz",
  # Non-antimicrobials
  "dalteparin", "enoxaparin", "heparin"
), collapse = "|")

# Normalize drug names before filtering (catches brand-named ITU prescriptions)
all_abx[, Drug_lower := tolower(trimws(Drug))]
all_abx[, Drug_lower := gsub("\\s*\\d+(\\.\\d+)?\\s*(mg|g|ml|mcg|units?|iu|tablets?|capsules?|doses?)\\s*/?\\s*\\d*(\\.\\d+)?\\s*(mg|g|ml|mcg|units?|iu|tablets?|capsules?|doses?)?", "", Drug_lower, ignore.case = TRUE)]
all_abx[, Drug_lower := gsub("\\([^)]*\\)", "", Drug_lower)]
all_abx[, Drug_lower := gsub("[+\\-\\/]", " ", Drug_lower)]
all_abx[, Drug_lower := gsub("\\s*\\([^)]*$", "", Drug_lower)]
all_abx[, Drug_lower := gsub("\\s*(injection|suspension|topical|oral|iv|im|sc|infusion|tablet|tablets|capsule|capsules|cream|ointment|gel|solution|powder|for injection|concentrate|pcc)\\s*$", "", Drug_lower, ignore.case = TRUE)]
all_abx[, Drug_lower := trimws(gsub("\\s+", " ", Drug_lower))]

all_abx[Drug_lower %like% "^augmentin", Drug_lower := "co-amoxiclav"]
all_abx[Drug_lower %like% "^tazocin" | Drug_lower %like% "^zocin", Drug_lower := "piperacillin tazobactam"]
all_abx[Drug_lower %like% "^zavicefta", Drug_lower := "ceftazidime avibactam"]
all_abx[Drug_lower %like% "^zyvox", Drug_lower := "linezolid"]
all_abx[Drug_lower %like% "^vancocin", Drug_lower := "vancomycin"]
all_abx[Drug_lower %like% "^flagyl", Drug_lower := "metronidazole"]
all_abx[Drug_lower %like% "^septrim" | Drug_lower %like% "^bactrim", Drug_lower := "co-trimoxazole"]
all_abx[Drug_lower %like% "amoxicillin.*clavulan", Drug_lower := "co-amoxiclav"]
all_abx[Drug_lower %like% "trimethoprim.*sulfa", Drug_lower := "co-trimoxazole"]
all_abx[Drug_lower %like% "piperacillin.*tazobactam", Drug_lower := "piperacillin tazobactam"]

# Filter medications after normalization (catches brand-named ITU prescriptions)
# - ITU: keep only antibiotics (ITU contains all meds)
# - EPR: exclude antivirals/antifungals but keep unmatched antibiotics in "Other" category
all_abx <- all_abx[
  (Source == "ITU" & grepl(all_antibiotic_patterns, Drug_lower, ignore.case = TRUE)) |
  (Source == "EPR" & !grepl(antiviral_antifungal_pattern, Drug_lower, ignore.case = TRUE))
]

# Classify antibiotics (do once before join)
all_abx[, `:=`(
  IsBetaLactamaseCombo = as.integer(grepl(beta_lactamase_combo_pattern, Drug_lower)),
  IsNarrowPenicillin = as.integer(grepl(narrow_penicillin_pattern, Drug_lower)),
  IsCeph1st2ndGen = as.integer(grepl(ceph_1st_2nd_pattern, Drug_lower)),
  IsCeph3rdGenPlus = as.integer(grepl(ceph_3rd_plus_pattern, Drug_lower)),
  IsCarbapenem = as.integer(grepl(carbapenem_pattern, Drug_lower)),
  IsFluoroquinolone = as.integer(grepl(fluoroquinolone_pattern, Drug_lower)),
  IsAminoglycoside = as.integer(grepl(aminoglycoside_pattern, Drug_lower)),
  IsGlycopeptide = as.integer(grepl(glycopeptide_pattern, Drug_lower)),
  IsMacrolide = as.integer(grepl(macrolide_pattern, Drug_lower)),
  IsTetracycline = as.integer(grepl(tetracycline_pattern, Drug_lower)),
  IsMetronidazole = as.integer(grepl(metronidazole_pattern, Drug_lower)),
  IsNitrofurantoin = as.integer(grepl(nitrofurantoin_pattern, Drug_lower)),
  IsCoTrimoxazole = as.integer(grepl(co_trimoxazole_pattern, Drug_lower)),
  IsBroadSpectrum = as.integer(grepl(broad_spectrum_pattern, Drug_lower))
)]


# Create time windows for each blood culture (-14 days to +1 hour)
bc_windows <- bc_episodes[, .(
  BloodCultureID, ClusterID, SpecimenFull, EpisodeID, CollectionDateTime,
  window_start = CollectionDateTime - days(14),
  window_end = CollectionDateTime + hours(1)
)]

# OPTIMIZATION: Pre-filter antibiotics to possible date range before expensive join
# Get overall date range with buffer
min_culture_date <- min(bc_episodes$CollectionDateTime, na.rm = TRUE)
max_culture_date <- max(bc_episodes$CollectionDateTime, na.rm = TRUE)
all_abx <- all_abx[FirstAdministrationDateTime >= (min_culture_date - days(14)) & 
                   FirstAdministrationDateTime <= (max_culture_date + hours(1))]

# Save original administration dates before join (they get overwritten)
all_abx[, `:=`(
  OrigFirstAdmin = FirstAdministrationDateTime,
  OrigLastAdmin = LastAdministrationDateTime
)]

# Set keys for faster joining
setkey(all_abx, ClusterID, FirstAdministrationDateTime, LastAdministrationDateTime)
setkey(bc_windows, ClusterID, window_start, window_end)

abx_bc <- all_abx[bc_windows,
                  on = .(ClusterID, 
                        FirstAdministrationDateTime <= window_end,
                        LastAdministrationDateTime >= window_start),
                  nomatch = 0L, allow.cartesian = TRUE]

# Restore original administration dates
abx_bc[, `:=`(
  FirstAdministrationDateTime = OrigFirstAdmin,
  LastAdministrationDateTime = OrigLastAdmin
)]

setDT(abx_bc)

# Create mutually exclusive antibiotic category for each prescription
# Priority order: Beta-lactamase combos > Carbapenems > Glycopeptides > 3rd gen+ ceph > 
# 1st/2nd gen ceph > Fluoroquinolones > Aminoglycosides > Narrow penicillins > 
# Macrolides > Tetracyclines > Metronidazole > Nitrofurantoin > Co-trimoxazole > Other
abx_bc[, AntibioticCategory := fcase(
  IsBetaLactamaseCombo == 1, "Beta-lactamase combo",
  IsCarbapenem == 1, "Carbapenem",
  IsGlycopeptide == 1, "Glycopeptide",
  IsCeph3rdGenPlus == 1, "3rd gen+ cephalosporin",
  IsCeph1st2ndGen == 1, "1st/2nd gen cephalosporin",
  IsFluoroquinolone == 1, "Fluoroquinolone",
  IsAminoglycoside == 1, "Aminoglycoside",
  IsNarrowPenicillin == 1, "Narrow penicillin",
  IsMacrolide == 1, "Macrolide",
  IsTetracycline == 1, "Tetracycline",
  IsMetronidazole == 1, "Metronidazole",
  IsNitrofurantoin == 1, "Nitrofurantoin",
  IsCoTrimoxazole == 1, "Co-trimoxazole",
  default = "Other"
)]

# ============================================================
# STEP 2B: SPLIT INTO BEFORE AND AFTER BLOOD CULTURE
# ============================================================

# Determine timing relative to blood culture based on administration period
# An antibiotic can be in both Before and After if it spans the BC time
# Leakage guard: TTAs use scheduled-start as their effective start, but a discharge order
# written at/after the culture (with a back-dated scheduled start) was not known at collection
# time, so it must not count as "before".
abx_bc[, TTO_OrderedAfterCulture := IsTTO == TRUE & !is.na(PrescriptionDateTime) &
                                    PrescriptionDateTime >= CollectionDateTime]
abx_bc[, `:=`(
  IsAdminBefore = LastAdministrationDateTime > (CollectionDateTime - days(14)) & 
                  FirstAdministrationDateTime < CollectionDateTime &
                  !TTO_OrderedAfterCulture,
  IsAdminAfter = FirstAdministrationDateTime < (CollectionDateTime + hours(1)) & 
                 LastAdministrationDateTime >= CollectionDateTime
)]

# Pre-compute route category for each antibiotic (faster than doing it in aggregation)
# Categories: IV, Oral (including enteral routes), Other (IM, SC, topical, etc.)
abx_bc[, Route_lower := tolower(Route)]
abx_bc[, IsIV := grepl("iv|intravenous|infusion", Route_lower, ignore.case = TRUE)]
abx_bc[, IsOral := grepl("oral|po|enteral|nasogastric|gastrostomy|jejunostomy|peg|gastrojejunostomy", Route_lower, ignore.case = TRUE)]
abx_bc[, RouteCategory := fcase(
  IsIV, "IV",
  IsOral, "Oral",
  default = "Other"
)]

# Add time-binned flags for different periods BEFORE blood culture
# Each flag indicates if the antibiotic administration overlapped with that time period
abx_bc[, `:=`(
  IsAdmin_24h_Before = FirstAdministrationDateTime < CollectionDateTime & 
                       LastAdministrationDateTime >= (CollectionDateTime - hours(24)) &
                       !TTO_OrderedAfterCulture,
  IsAdmin_1wk_Before = FirstAdministrationDateTime < CollectionDateTime & 
                       LastAdministrationDateTime >= (CollectionDateTime - days(7)) &
                       !TTO_OrderedAfterCulture,
  IsAdmin_2wk_Before = FirstAdministrationDateTime < CollectionDateTime & 
                       LastAdministrationDateTime >= (CollectionDateTime - days(14)) &
                       !TTO_OrderedAfterCulture
)]

# Calculate ALL metrics in a single pass through the data for efficiency
# Aggregate from abx_bc (only blood cultures with antibiotics)
bc_abx_metrics <- abx_bc[, .(
  # BEFORE metrics (14 days)
  HasAntibiotics_Before = as.integer(any(IsAdminBefore)),
  BroadSpectrum_Before = as.integer(any(IsBroadSpectrum == 1 & IsAdminBefore)),
  
  # BEFORE (14 days)
  BetaLactamaseCombo_Before = as.integer(any(IsBetaLactamaseCombo == 1 & IsAdminBefore)),
  NarrowPenicillin_Before = as.integer(any(IsNarrowPenicillin == 1 & IsAdminBefore)),
  Ceph1st2ndGen_Before = as.integer(any(IsCeph1st2ndGen == 1 & IsAdminBefore)),
  Ceph3rdGenPlus_Before = as.integer(any(IsCeph3rdGenPlus == 1 & IsAdminBefore)),
  Carbapenem_Before = as.integer(any(IsCarbapenem == 1 & IsAdminBefore)),
  Fluoroquinolone_Before = as.integer(any(IsFluoroquinolone == 1 & IsAdminBefore)),
  Aminoglycoside_Before = as.integer(any(IsAminoglycoside == 1 & IsAdminBefore)),
  Glycopeptide_Before = as.integer(any(IsGlycopeptide == 1 & IsAdminBefore)),
  Macrolide_Before = as.integer(any(IsMacrolide == 1 & IsAdminBefore)),
  Tetracycline_Before = as.integer(any(IsTetracycline == 1 & IsAdminBefore)),
  Metronidazole_Before = as.integer(any(IsMetronidazole == 1 & IsAdminBefore)),
  Nitrofurantoin_Before = as.integer(any(IsNitrofurantoin == 1 & IsAdminBefore)),
  CoTrimoxazole_Before = as.integer(any(IsCoTrimoxazole == 1 & IsAdminBefore)),
  Other_Before = as.integer(any(AntibioticCategory == "Other" & IsAdminBefore)),
  
  # AFTER metrics (14 days)
  HasAntibiotics_After = as.integer(any(IsAdminAfter)),
  BroadSpectrum_After = as.integer(any(IsBroadSpectrum == 1 & IsAdminAfter)),
  
  # AFTER (14 days)
  BetaLactamaseCombo_After = as.integer(any(IsBetaLactamaseCombo == 1 & IsAdminAfter)),
  NarrowPenicillin_After = as.integer(any(IsNarrowPenicillin == 1 & IsAdminAfter)),
  Ceph1st2ndGen_After = as.integer(any(IsCeph1st2ndGen == 1 & IsAdminAfter)),
  Ceph3rdGenPlus_After = as.integer(any(IsCeph3rdGenPlus == 1 & IsAdminAfter)),
  Carbapenem_After = as.integer(any(IsCarbapenem == 1 & IsAdminAfter)),
  Fluoroquinolone_After = as.integer(any(IsFluoroquinolone == 1 & IsAdminAfter)),
  Aminoglycoside_After = as.integer(any(IsAminoglycoside == 1 & IsAdminAfter)),
  Glycopeptide_After = as.integer(any(IsGlycopeptide == 1 & IsAdminAfter)),
  Macrolide_After = as.integer(any(IsMacrolide == 1 & IsAdminAfter)),
  Tetracycline_After = as.integer(any(IsTetracycline == 1 & IsAdminAfter)),
  Metronidazole_After = as.integer(any(IsMetronidazole == 1 & IsAdminAfter)),
  Nitrofurantoin_After = as.integer(any(IsNitrofurantoin == 1 & IsAdminAfter)),
  CoTrimoxazole_After = as.integer(any(IsCoTrimoxazole == 1 & IsAdminAfter)),
  Other_After = as.integer(any(AntibioticCategory == "Other" & IsAdminAfter)),
  
  # Time-binned metrics (24h, 1 week)
  HasAntibiotics_24h_Before = as.integer(any(IsAdmin_24h_Before)),
  NumDistinctAntibiotics_24h_Before = uniqueN(Drug_lower[IsAdmin_24h_Before]),
  HasAntibiotics_1wk_Before = as.integer(any(IsAdmin_1wk_Before)),
  NumDistinctAntibiotics_1wk_Before = uniqueN(Drug_lower[IsAdmin_1wk_Before]),
  HasAntibiotics_2wk_Before = as.integer(any(IsAdmin_2wk_Before)),
  
  # Days of continuous antibiotic pre-culture (WITH GAP CHECKING)
  # Checks for gaps >48h in antibiotic coverage
  DaysContinuousAbx_PreCulture = {
    if(all(!IsAdminBefore)) {
      0L
    } else {
      culture_dt <- first(CollectionDateTime)
      # Get antibiotics that were active before culture
      abx_before <- data.table(
        Start = FirstAdministrationDateTime[IsAdminBefore],
        End = LastAdministrationDateTime[IsAdminBefore]
      )
      # Sort by start time
      abx_before <- abx_before[order(Start)]
      
      if(nrow(abx_before) == 0) {
        0L
      } else {
        # Find the most recent continuous period ending at or before culture
        # Work backwards from culture time
        continuous_start <- culture_dt
        
        for(i in nrow(abx_before):1) {
          # Does this antibiotic reach back in time from our continuous period?
          # Allow 48h gap (antibiotics might not be perfectly recorded)
          if(abx_before$End[i] >= continuous_start - hours(48)) {
            # Extend continuous period backwards
            continuous_start <- min(continuous_start, abx_before$Start[i])
          }
        }
        
        # Calculate days from continuous start to culture
        as.integer(ceiling(difftime(culture_dt, continuous_start, units = "days")))
      }
    }
  },
  

  # Sum doses for antibiotics given in the 14 days before culture
  NumDosesAdministered_WeekBefore = sum(NumberOfDosesAdministered[IsAdminBefore], na.rm = TRUE),
  
  # Categorize route based on antibiotics given before culture (pre-computed for speed)
  # If IV present, mark as IV; otherwise use most common route category
  Abx_Route_Category = {
    if(any(IsAdminBefore)) {
      if(any(IsIV[IsAdminBefore], na.rm = TRUE)) {
        "IV"
      } else {
        # Find most common route among antibiotics given before
        routes <- RouteCategory[IsAdminBefore]
        routes <- routes[!is.na(routes)]
        if(length(routes) == 0) {
          NA_character_
        } else {
          # Return most frequent route (excluding IV which we already checked)
          route_counts <- table(routes)
          names(route_counts)[which.max(route_counts)]
        }
      }
    } else {
      NA_character_
    }
  },
  
  # TTO metrics - track prescriptions to take home
  HasTTO_Before = as.integer(any(IsTTO & IsAdminBefore)),
  HasTTO_After = as.integer(any(IsTTO & IsAdminAfter))
), by = BloodCultureID]


bc_sources_indication_final <- abx_bc[IsAdminBefore == TRUE, .(
  urinary          = as.integer(any(urinary == 1,          na.rm = TRUE)),
  respiratory      = as.integer(any(respiratory == 1,      na.rm = TRUE)),
  abdominal        = as.integer(any(abdominal == 1,        na.rm = TRUE)),
  neurological     = as.integer(any(neurological == 1,     na.rm = TRUE)),
  skin_soft_tissue = as.integer(any(skin_soft_tissue == 1, na.rm = TRUE)),
  ent              = as.integer(any(ent == 1,              na.rm = TRUE)),
  orthopaedic      = as.integer(any(orthopaedic == 1,      na.rm = TRUE)),
  other_specific   = as.integer(any(other_specific == 1,   na.rm = TRUE)),
  no_specific_source = as.integer(any(no_specific_source == 1, na.rm = TRUE))
), by = BloodCultureID]

# Calculate sum of specific sources for each row
bc_sources_indication_final[, has_specific := rowSums(.SD), .SDcols = source_categories]

# If any specific source is present, set no_specific_source to 0
bc_sources_indication_final[has_specific > 0, no_specific_source := 0]

# Calculate MultipleSourcesFlag (2+ specific sources)
bc_sources_indication_final[, MultipleSourcesFlag := as.integer(has_specific >= 2)]

# Remove temporary column
bc_sources_indication_final[, has_specific := NULL]

# Rename columns to be more readable
setnames(bc_sources_indication_final,
         old = c("urinary", "respiratory", "abdominal", "neurological", 
                 "skin_soft_tissue", "ent", "orthopaedic", "other_specific", "no_specific_source"),
         new = c("Source_Urinary", "Source_Respiratory", "Source_Abdominal", "Source_CNS",
                 "Source_Skin_Soft_Tissue", "Source_ENT", "Source_Orthopaedic", "Source_Other", "Source_Unspecific"))

# ============================================================
# STEP 3: COMBINE AND SAVE
# ============================================================

# Join source information and antibiotic metrics to all blood cultures using BloodCultureID
# Using data.table joins for performance
setDT(bc_with_sources <- bc_episodes)
bc_with_sources <- bc_sources_indication_final[bc_with_sources, on = "BloodCultureID"]
bc_with_sources <- bc_abx_metrics[bc_with_sources, on = "BloodCultureID"]

# Define source column names (matching the renamed columns from bc_sources_indication_final)
source_cols <- c("Source_Urinary", "Source_Respiratory", "Source_Abdominal", "Source_CNS",
                 "Source_Skin_Soft_Tissue", "Source_ENT", "Source_Orthopaedic", "Source_Other", "Source_Unspecific")

# Define antibiotic column names - NEW CLASSIFICATION SYSTEM
abx_cols_before <- c("HasAntibiotics_Before", "BroadSpectrum_Before",
                    "BetaLactamaseCombo_Before", "NarrowPenicillin_Before",
                    "Ceph1st2ndGen_Before", "Ceph3rdGenPlus_Before",
                    "Carbapenem_Before", "Fluoroquinolone_Before",
                    "Aminoglycoside_Before", "Glycopeptide_Before",
                    "Macrolide_Before", "Tetracycline_Before",
                    "Metronidazole_Before", "Nitrofurantoin_Before", "CoTrimoxazole_Before", "Other_Before")

abx_cols_after <- c("HasAntibiotics_After", "BroadSpectrum_After",
                   "BetaLactamaseCombo_After", "NarrowPenicillin_After",
                   "Ceph1st2ndGen_After", "Ceph3rdGenPlus_After",
                   "Carbapenem_After", "Fluoroquinolone_After",
                   "Aminoglycoside_After", "Glycopeptide_After",
                   "Macrolide_After", "Tetracycline_After",
                   "Metronidazole_After", "Nitrofurantoin_After", "CoTrimoxazole_After", "Other_After")

# Define time-binned antibiotic column names
abx_cols_timebin <- c("HasAntibiotics_24h_Before", "NumDistinctAntibiotics_24h_Before",
                      "HasAntibiotics_1wk_Before", "NumDistinctAntibiotics_1wk_Before",
                      "HasAntibiotics_2wk_Before")

# Define TTO column names
tta_cols <- c("HasTTO_Before", "HasTTO_After")

# Define additional numeric columns that should be 0 when no antibiotics
abx_cols_numeric <- c("DaysContinuousAbx_PreCulture", "NumDosesAdministered_WeekBefore")

# Replace NAs with 0 using data.table syntax (much faster)
# Note: Abx_Route_Category is intentionally left as NA for blood cultures without antibiotics
na_cols <- c(source_cols, "MultipleSourcesFlag", abx_cols_before, abx_cols_after, abx_cols_timebin, tta_cols, abx_cols_numeric)
for (col in na_cols) {
  bc_with_sources[is.na(get(col)), (col) := 0]
}

fwrite(bc_with_sources, "Output Files/CSVs/blood_cultures_with_infection_source.csv")

# ============================================================
# STEP 4: SUMMARY EXCEL - BY BLOOD CULTURE
# ============================================================

# Calculate if blood culture has any source identified (excluding Unspecific)
# Using data.table for efficiency
source_cols_specific <- setdiff(source_cols, "Source_Unspecific")
bc_with_sources[, HasAnySource := rowSums(.SD, na.rm = TRUE) > 0, .SDcols = source_cols_specific]
bc_with_sources[, HasOnlyUnspecific := Source_Unspecific == 1 & !HasAnySource]
bc_with_sources[, HasNoSource := !HasAnySource & Source_Unspecific == 0]

# Overall summary
summary <- data.frame(
  Metric = c("Total blood cultures", "With specific infection source", "Only unspecific source", "No source identified"),
  Total = c(
    nrow(bc_with_sources),
    sum(bc_with_sources$HasAnySource, na.rm = TRUE),
    sum(bc_with_sources$HasOnlyUnspecific, na.rm = TRUE),
    sum(bc_with_sources$HasNoSource, na.rm = TRUE)
  ),
  Positive = c(
    sum(bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$HasAnySource & bc_with_sources$CultureCategory == "Positive", na.rm = TRUE),
    sum(bc_with_sources$HasOnlyUnspecific & bc_with_sources$CultureCategory == "Positive", na.rm = TRUE),
    sum(bc_with_sources$HasNoSource & bc_with_sources$CultureCategory == "Positive", na.rm = TRUE)
  ),
  `Negative (and single contam)` = c(
    sum(bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$HasAnySource & bc_with_sources$CultureCategory == "Negative (and SC)", na.rm = TRUE),
    sum(bc_with_sources$HasOnlyUnspecific & bc_with_sources$CultureCategory == "Negative (and SC)", na.rm = TRUE),
    sum(bc_with_sources$HasNoSource & bc_with_sources$CultureCategory == "Negative (and SC)", na.rm = TRUE)
  ),
  check.names = FALSE
) %>%
  mutate(
    Total_Pct = 100 * Total / first(Total),
    Positive_Pct = 100 * Positive / first(Positive),
    Negative_Pct = 100 * .data[["Negative (and single contam)"]] / first(.data[["Negative (and single contam)"]])
  )

# Breakdown by individual infection sources (each source gets its own row)
source_breakdown <- bc_with_sources %>%
  select(all_of(source_cols), MultipleSourcesFlag, CultureCategory) %>%
  pivot_longer(cols = c(all_of(source_cols), MultipleSourcesFlag),
               names_to = "Source", values_to = "HasSource") %>%
  filter(HasSource == 1) %>%
  count(Source, CultureCategory) %>%
  pivot_wider(names_from = CultureCategory, values_from = n, values_fill = 0) %>%
  mutate(Total = Positive + `Negative (and SC)`) %>%
  select(Source, Total, Positive, `Negative (and SC)`) %>%
  arrange(desc(Total)) %>%
  mutate(
    Total_Pct = 100 * Total / nrow(bc_with_sources),
    Positive_Pct = 100 * Positive / sum(bc_with_sources$CultureCategory == "Positive"),
    Negative_SC_Pct = 100 * `Negative (and SC)` / sum(bc_with_sources$CultureCategory == "Negative (and SC)")
  )


# Antibiotic usage summary 
abx_before_summary <- data.frame(
  Metric = c(
    "Any antibiotic before BC",
    "Broad-spectrum antibiotic before BC",
    "Beta-lactamase combo before BC",
    "Narrow penicillin before BC",
    "1st/2nd gen cephalosporin before BC",
    "3rd gen+ cephalosporin before BC",
    "Carbapenem before BC",
    "Fluoroquinolone before BC",
    "Aminoglycoside before BC",
    "Glycopeptide before BC",
    "Macrolide before BC",
    "Tetracycline before BC",
    "Metronidazole before BC",
    "Nitrofurantoin before BC",
    "Co-trimoxazole before BC",
    "Other antibiotic before BC",
    "Days of continuous antibiotic pre-culture",
    "Mean days of continuous antibiotic pre-culture",
    "Total doses administered (week before)",
    "Mean doses administered (week before)",
    "Route: IV",
    "Route: Oral",
    "TTO prescription before BC"
  ),
  Total = c(
    sum(bc_with_sources$HasAntibiotics_Before),
    sum(bc_with_sources$BroadSpectrum_Before),
    sum(bc_with_sources$BetaLactamaseCombo_Before),
    sum(bc_with_sources$NarrowPenicillin_Before),
    sum(bc_with_sources$Ceph1st2ndGen_Before),
    sum(bc_with_sources$Ceph3rdGenPlus_Before),
    sum(bc_with_sources$Carbapenem_Before),
    sum(bc_with_sources$Fluoroquinolone_Before),
    sum(bc_with_sources$Aminoglycoside_Before),
    sum(bc_with_sources$Glycopeptide_Before),
    sum(bc_with_sources$Macrolide_Before),
    sum(bc_with_sources$Tetracycline_Before),
    sum(bc_with_sources$Metronidazole_Before),
    sum(bc_with_sources$Nitrofurantoin_Before),
    sum(bc_with_sources$CoTrimoxazole_Before),
    sum(bc_with_sources$Other_Before),
    sum(bc_with_sources$DaysContinuousAbx_PreCulture > 0, na.rm = TRUE),
    round(mean(bc_with_sources$DaysContinuousAbx_PreCulture, na.rm = TRUE), 1),
    sum(bc_with_sources$NumDosesAdministered_WeekBefore, na.rm = TRUE),
    round(mean(bc_with_sources$NumDosesAdministered_WeekBefore, na.rm = TRUE), 1),
    sum(bc_with_sources$Abx_Route_Category == "IV", na.rm = TRUE),
    sum(bc_with_sources$Abx_Route_Category == "Oral", na.rm = TRUE),
    sum(bc_with_sources$HasTTO_Before)
  ),
  Positive = c(
    sum(bc_with_sources$HasAntibiotics_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$BroadSpectrum_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$BetaLactamaseCombo_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$NarrowPenicillin_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Ceph1st2ndGen_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Ceph3rdGenPlus_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Carbapenem_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Fluoroquinolone_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Aminoglycoside_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Glycopeptide_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Macrolide_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Tetracycline_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Metronidazole_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Nitrofurantoin_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$CoTrimoxazole_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Other_Before & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$DaysContinuousAbx_PreCulture > 0 & bc_with_sources$CultureCategory == "Positive", na.rm = TRUE),
    round(mean(bc_with_sources$DaysContinuousAbx_PreCulture[bc_with_sources$CultureCategory == "Positive"], na.rm = TRUE), 1),
    sum(bc_with_sources$NumDosesAdministered_WeekBefore[bc_with_sources$CultureCategory == "Positive"], na.rm = TRUE),
    round(mean(bc_with_sources$NumDosesAdministered_WeekBefore[bc_with_sources$CultureCategory == "Positive"], na.rm = TRUE), 1),
    sum(bc_with_sources$Abx_Route_Category == "IV" & bc_with_sources$CultureCategory == "Positive", na.rm = TRUE),
    sum(bc_with_sources$Abx_Route_Category == "Oral" & bc_with_sources$CultureCategory == "Positive", na.rm = TRUE),
    sum(bc_with_sources$HasTTO_Before & bc_with_sources$CultureCategory == "Positive")
  ),
  `Negative (and single contam)` = c(
    sum(bc_with_sources$HasAntibiotics_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$BroadSpectrum_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$BetaLactamaseCombo_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$NarrowPenicillin_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Ceph1st2ndGen_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Ceph3rdGenPlus_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Carbapenem_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Fluoroquinolone_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Aminoglycoside_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Glycopeptide_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Macrolide_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Tetracycline_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Metronidazole_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Nitrofurantoin_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$CoTrimoxazole_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Other_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$DaysContinuousAbx_PreCulture > 0 & bc_with_sources$CultureCategory == "Negative (and SC)", na.rm = TRUE),
    round(mean(bc_with_sources$DaysContinuousAbx_PreCulture[bc_with_sources$CultureCategory == "Negative (and SC)"], na.rm = TRUE), 1),
    sum(bc_with_sources$NumDosesAdministered_WeekBefore[bc_with_sources$CultureCategory == "Negative (and SC)"], na.rm = TRUE),
    round(mean(bc_with_sources$NumDosesAdministered_WeekBefore[bc_with_sources$CultureCategory == "Negative (and SC)"], na.rm = TRUE), 1),
    sum(bc_with_sources$Abx_Route_Category == "IV" & bc_with_sources$CultureCategory == "Negative (and SC)", na.rm = TRUE),
    sum(bc_with_sources$Abx_Route_Category == "Oral" & bc_with_sources$CultureCategory == "Negative (and SC)", na.rm = TRUE),
    sum(bc_with_sources$HasTTO_Before & bc_with_sources$CultureCategory == "Negative (and SC)")
  ),
  check.names = FALSE
) %>%
  mutate(
    Total_Pct = ifelse(row_number() %in% c(18, 20, 22, 23), NA, 100 * Total / nrow(bc_with_sources)),
    Positive_Pct = ifelse(row_number() %in% c(18, 20, 22, 23), NA, 100 * Positive / sum(bc_with_sources$CultureCategory == "Positive")),
    Negative_Pct = ifelse(row_number() %in% c(18, 20, 22, 23), NA, 100 * .data[["Negative (and single contam)"]] / sum(bc_with_sources$CultureCategory == "Negative (and SC)"))
  )

# Antibiotic usage summary - AFTER blood culture (NEW CLASSIFICATION)
abx_after_summary <- data.frame(
  Metric = c(
    "Any antibiotic after BC",
    "Broad-spectrum antibiotic after BC",
    "Beta-lactamase combo after BC",
    "Narrow penicillin after BC",
    "1st/2nd gen cephalosporin after BC",
    "3rd gen+ cephalosporin after BC",
    "Carbapenem after BC",
    "Fluoroquinolone after BC",
    "Aminoglycoside after BC",
    "Glycopeptide after BC",
    "Macrolide after BC",
    "Tetracycline after BC",
    "Metronidazole after BC",
    "Nitrofurantoin after BC",
    "Co-trimoxazole after BC",
    "Other antibiotic after BC",
    "TTO prescription after BC"
  ),
  Total = c(
    sum(bc_with_sources$HasAntibiotics_After),
    sum(bc_with_sources$BroadSpectrum_After),
    sum(bc_with_sources$BetaLactamaseCombo_After),
    sum(bc_with_sources$NarrowPenicillin_After),
    sum(bc_with_sources$Ceph1st2ndGen_After),
    sum(bc_with_sources$Ceph3rdGenPlus_After),
    sum(bc_with_sources$Carbapenem_After),
    sum(bc_with_sources$Fluoroquinolone_After),
    sum(bc_with_sources$Aminoglycoside_After),
    sum(bc_with_sources$Glycopeptide_After),
    sum(bc_with_sources$Macrolide_After),
    sum(bc_with_sources$Tetracycline_After),
    sum(bc_with_sources$Metronidazole_After),
    sum(bc_with_sources$Nitrofurantoin_After),
    sum(bc_with_sources$CoTrimoxazole_After),
    sum(bc_with_sources$Other_After),
    sum(bc_with_sources$HasTTO_After)
  ),
  Positive = c(
    sum(bc_with_sources$HasAntibiotics_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$BroadSpectrum_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$BetaLactamaseCombo_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$NarrowPenicillin_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Ceph1st2ndGen_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Ceph3rdGenPlus_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Carbapenem_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Fluoroquinolone_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Aminoglycoside_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Glycopeptide_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Macrolide_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Tetracycline_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Metronidazole_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Nitrofurantoin_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$CoTrimoxazole_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$Other_After & bc_with_sources$CultureCategory == "Positive"),
    sum(bc_with_sources$HasTTO_After & bc_with_sources$CultureCategory == "Positive")
  ),
  `Negative (and single contam)` = c(
    sum(bc_with_sources$HasAntibiotics_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$BroadSpectrum_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$BetaLactamaseCombo_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$NarrowPenicillin_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Ceph1st2ndGen_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Ceph3rdGenPlus_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Carbapenem_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Fluoroquinolone_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Aminoglycoside_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Glycopeptide_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Macrolide_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Tetracycline_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Metronidazole_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Nitrofurantoin_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$CoTrimoxazole_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$Other_After & bc_with_sources$CultureCategory == "Negative (and SC)"),
    sum(bc_with_sources$HasTTO_After & bc_with_sources$CultureCategory == "Negative (and SC)")
  ),
  check.names = FALSE
) %>%
  mutate(
    Total_Pct = 100 * Total / nrow(bc_with_sources),
    Positive_Pct = 100 * Positive / sum(bc_with_sources$CultureCategory == "Positive"),
    Negative_Pct = 100 * .data[["Negative (and single contam)"]] / sum(bc_with_sources$CultureCategory == "Negative (and SC)")
  )

# Time-binned antibiotic summary - BEFORE blood culture
abx_timebin_summary <- data.frame(
  Metric = c(
    "Any antibiotic 24h before BC",
    "Mean distinct antibiotics 24h before BC",
    "Any antibiotic 1 week before BC",
    "Mean distinct antibiotics 1 week before BC"
  ),
  Total = c(
    sum(bc_with_sources$HasAntibiotics_24h_Before),
    mean(bc_with_sources$NumDistinctAntibiotics_24h_Before),
    sum(bc_with_sources$HasAntibiotics_1wk_Before),
    mean(bc_with_sources$NumDistinctAntibiotics_1wk_Before)
  ),
  Positive = c(
    sum(bc_with_sources$HasAntibiotics_24h_Before & bc_with_sources$CultureCategory == "Positive"),
    mean(bc_with_sources$NumDistinctAntibiotics_24h_Before[bc_with_sources$CultureCategory == "Positive"]),
    sum(bc_with_sources$HasAntibiotics_1wk_Before & bc_with_sources$CultureCategory == "Positive"),
    mean(bc_with_sources$NumDistinctAntibiotics_1wk_Before[bc_with_sources$CultureCategory == "Positive"])
  ),
  `Negative (and single contam)` = c(
    sum(bc_with_sources$HasAntibiotics_24h_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    round(mean(bc_with_sources$NumDistinctAntibiotics_24h_Before[bc_with_sources$CultureCategory == "Negative (and SC)"]), 2),
    sum(bc_with_sources$HasAntibiotics_1wk_Before & bc_with_sources$CultureCategory == "Negative (and SC)"),
    round(mean(bc_with_sources$NumDistinctAntibiotics_1wk_Before[bc_with_sources$CultureCategory == "Negative (and SC)"]), 2)
  ),
  check.names = FALSE
) %>%
  mutate(
    Total_Pct = ifelse(row_number() %% 2 == 0, NA, 100 * Total / nrow(bc_with_sources)),
    Positive_Pct = ifelse(row_number() %% 2 == 0, NA, 100 * Positive / sum(bc_with_sources$CultureCategory == "Positive")),
    Negative_Pct = ifelse(row_number() %% 2 == 0, NA, 100 * .data[["Negative (and single contam)"]] / sum(bc_with_sources$CultureCategory == "Negative (and SC)"))
  )

# Format all summary tables before writing to Excel
summary_formatted <- summary %>% format_excel_output()
source_breakdown_formatted <- source_breakdown %>% format_excel_output()
abx_before_summary_formatted <- abx_before_summary %>% format_excel_output()
abx_after_summary_formatted <- abx_after_summary %>% format_excel_output()
abx_timebin_summary_formatted <- abx_timebin_summary %>% format_excel_output()

wb <- createWorkbook()
addWorksheet(wb, "Summary")

# Overall summary
writeData(wb, "Summary", summary_formatted, startRow = 1)

# Infection source breakdown
source_start <- nrow(summary_formatted) + 3
writeData(wb, "Summary", data.frame(Section = "Infection Source Breakdown"), 
          startRow = source_start, colNames = FALSE)
writeData(wb, "Summary", source_breakdown_formatted, startRow = source_start + 1)

# Antibiotic usage BEFORE blood culture
abx_before_start <- source_start + 1 + nrow(source_breakdown_formatted) + 2
writeData(wb, "Summary", data.frame(Section = "Antibiotic Usage - BEFORE Blood Culture (14 days)"), 
          startRow = abx_before_start, colNames = FALSE)
writeData(wb, "Summary", abx_before_summary_formatted, startRow = abx_before_start + 1)

# Antibiotic usage AFTER blood culture
abx_after_start <- abx_before_start + 1 + nrow(abx_before_summary_formatted) + 2
writeData(wb, "Summary", data.frame(Section = "Antibiotic Usage - AFTER Blood Culture (14 days)"), 
          startRow = abx_after_start, colNames = FALSE)
writeData(wb, "Summary", abx_after_summary_formatted, startRow = abx_after_start + 1)

# Time-binned antibiotic usage
abx_timebin_start <- abx_after_start + 1 + nrow(abx_after_summary_formatted) + 2
writeData(wb, "Summary", data.frame(Section = "Antibiotic Usage - Time-binned BEFORE Blood Culture"), 
          startRow = abx_timebin_start, colNames = FALSE)
writeData(wb, "Summary", abx_timebin_summary_formatted, startRow = abx_timebin_start + 1)

# Add formatting note
addWorksheet(wb, "Formatting_Note")
writeData(wb, "Formatting_Note", create_formatting_note())

saveWorkbook(wb, "Output Files/Excels/Pipeline/BC_4_Add_Infection_Source_excel.xlsx", overwrite = TRUE)
