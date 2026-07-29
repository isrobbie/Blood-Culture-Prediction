library(data.table)
library(lubridate)

# ============================================================
# ADD OTHER MICROBIOLOGY RESULTS TO NON-BC ADMISSIONS
# ============================================================
# Equivalent to BC_8. Looks back 28 days before admission for
# positive non-blood-culture microbiology results.
# Uses full exclusion criteria (viral tests, screening tests,
# blood culture specimens, nasal/rectal swabs) matching BC_8.
# All flags including Any_Prior_Other_Pos and Prior_Urine_Pos
# are derived from the same fully-filtered dataset.
# ============================================================

nbc_df <- fread("Output Files/CSVs/NonBC_Admissions/nbc_with_blood_tests.csv",
                colClasses = c(ClusterID = "character"))
nbc_df[, AdmissionDate := as.POSIXct(AdmissionDate, tz = "UTC")]

cohort_ids <- unique(nbc_df$ClusterID)

# ============================================================
# 1. LOAD MICROBIOLOGY AND BUG LOOKUP
# ============================================================

bug_lookup <- fread("00_Data/lookup_tables/bug_group_table.csv")
setkey(bug_lookup, BugCode)

micro_df <- fread("Input Files/IORD CSVs/IORD_BloodCultPatients_57_20251030_Microbiology.csv",
                  select = c("ClusterID", "AccessionNumber", "CollectionDateTime", "BatTestName",
                             "BugCode", "BugName", "SpecimenFull"))

pathogen_groups <- c("S AUREUS", "E COLI", "KLEBSIELLA", "P AERUGINOSA",
                     "ENTEROCOCCUS", "S PNEUMONIAE", "ENTEROBACTER",
                     "BETA HAEMOLYTIC STREP", "OTHER ENTEROBACTERALES",
                     "OTHER PATHOGENIC STREPS", "ANAEROBES", "CANDIDA",
                     "OTHER GRAM-POSITIVE", "OTHER GRAM-NEGATIVE")
contaminant_groups <- c("VIRIDANS AND OTHER STREPS", "CoNS (CONTAMINANT)", "OTHER CONTAMINANT")

# ============================================================
# 2. EXTENDED PER-SOURCE AND GRAM CLASSIFICATION
# ============================================================
# Uses full exclusion criteria matching BC_8: excludes viral tests,
# screening tests, blood culture specimens, and nasal/rectal swabs.
# Prior_Urine_Pos and Any_Prior_Other_Pos are derived from this
# same filtered dataset for consistency.

nbc_df[, adm_minus_28d := AdmissionDate - days(28)]

viral_tests_bc8     <- c("PCR VIROLOGY","FLU & RSV A/B PCR","RESPIRATORY PCR",
                         "NOROVIRUS PCR","MONKEYPOX PCR","REF LAB NOROVIRUS",
                         "LOWER RESPIRATORY PCR","REF RESPIRATORY PCR",
                         "Ref Lab Resp PCR","FAECES VIRAL NOT DONE")
screening_tests_bc8 <- c("MRSA SCREEN","ESBL SCREEN","CPE SCREEN",
                         "CANDIDA AURIS SCREEN","GROUP B STREP SCREEN",
                         "VRE SCREEN","CRE SCREEN")
grampos_groups_bc8  <- c("S AUREUS","ENTEROCOCCUS","S PNEUMONIAE","BETA HAEMOLYTIC STREP","OTHER PATHOGENIC STREPS","OTHER GRAM-POSITIVE")
gramneg_groups_bc8  <- c("E COLI","KLEBSIELLA","P AERUGINOSA","ENTEROBACTER","OTHER ENTEROBACTERALES","OTHER GRAM-NEGATIVE")
all_groups_bc8      <- c(pathogen_groups, contaminant_groups)

other_micro_full <- micro_df[
  BatTestName != "BLOOD CULTURE" &
  !is.na(BugName) & BugName != "" & BugName != "NULL" &
  !is.na(CollectionDateTime) &
  !(BatTestName %in% viral_tests_bc8) &
  !(BatTestName %in% screening_tests_bc8) &
  !grepl("virus|VIRUS", BugName, ignore.case = TRUE) &
  !grepl("BLOOD.*FOR.*CULTURE|BLOOD.*CULTURE", SpecimenFull, ignore.case = TRUE) &
  !grepl("^NOSE SWAB$|^NASAL SWAB$|RECTAL SWAB|RECTAL/STOOL SWAB", SpecimenFull, ignore.case = TRUE)]

other_micro_full[, ClusterID := as.character(ClusterID)]
other_micro_full <- other_micro_full[ClusterID %in% cohort_ids]
other_micro_full[, CollectionDateTime := as.POSIXct(CollectionDateTime, tz = "UTC")]
other_micro_full <- other_micro_full[!is.na(CollectionDateTime)]
other_micro_full[, BugGroup_ext := bug_lookup[other_micro_full, BugGroup, on = "BugCode"]]
other_micro_full <- other_micro_full[!is.na(BugGroup_ext) & BugGroup_ext != "" &
                                      BugGroup_ext != "NULL" & BugGroup_ext %in% all_groups_bc8]
other_micro_full <- unique(other_micro_full, by = c("ClusterID", "AccessionNumber", "BugName"))

other_micro_full[, `:=`(
  IsGramPos_ext = BugGroup_ext %in% grampos_groups_bc8,
  IsGramNeg_ext = BugGroup_ext %in% gramneg_groups_bc8
)]

pathogenic_full <- other_micro_full[
  BugGroup_ext %in% setdiff(pathogen_groups, "CANDIDA") &
  !(BugGroup_ext %in% contaminant_groups)]

pathogenic_full[, SpecimenUpper := toupper(SpecimenFull)]
pathogenic_full[, SourceCategory := fcase(
  grepl("URINE|CSU|MSU|CATHETER URINE", SpecimenUpper), "Urine",
  grepl("SPUTUM|BRONCH|LAVAGE|TRACHEAL|RESPIRATORY|THROAT|NASAL|NOSE|SWAB.*VIRAL|COUGH|NASOPHARYNGEAL|NASO PHARYNGEAL|MOUTH SWAB", SpecimenUpper), "Respiratory",
  grepl("EYE|CORNEA|CORNEAL|CONJUNCTIVAL|EYELID", SpecimenUpper), "Eye",
  grepl("VAGINAL|VAGINA|CERVIX|VULVAL|PENILE|PENIS", SpecimenUpper), "Genital",
  (grepl("WOUND|SKIN|ULCER|PUS|ABSCESS|TISSUE|BIOPSY|ASPIRATE", SpecimenUpper) |
   (grepl("(^SWAB|SWAB$)", SpecimenUpper) &
    grepl("FOOT|GROIN|BACK|LEG|ABDOMEN|ARM|HAND|CHEST|STERNUM|THIGH|SHOULDER|ELBOW|WRIST|ANKLE|TOE|FINGER|BUTTOCK|GLUTEAL|PERINEUM|PEG|SITE|AXILLA|UMBILICUS|BREAST", SpecimenUpper))) &
   !grepl("THROAT|NASAL|NOSE|RECTAL|VAGINAL|URETHRAL|EYE|EAR|COUGH|NASOPHARYNGEAL|MOUTH|VULVAL|PENILE|PENIS|CERVIX|CONJUNCTIVAL|RESPIRATORY|BRONCH|SPUTUM|CORNEA", SpecimenUpper) &
   !(SpecimenUpper == "SWAB") &
   !grepl("^SWAB;OTHER|SWAB OTHER|SWAB:UNSPECIFIED|SWAB;UNSPECIFIED|SWAB;UNKNOWN|SWAB; Other", SpecimenUpper), "Skin_SoftTissue",
  grepl("CSF|CEREBROSPINAL", SpecimenUpper), "CSF",
  grepl("ASCITIC|PERITONEAL|JOINT|SYNOVIAL|SYNOVIUM|BILE|CAPD|KIDNEY PERFUSION|PANCREATIC PERFUSION|SONICATION|BURSA|PERICARDIAL|PLEURAL", SpecimenUpper) |
  (grepl("FLUID", SpecimenUpper) &
   grepl("ASCITIC|PERITONEAL|JOINT|SYNOVIAL|BILE|CAPD|KIDNEY|PANCREATIC|SONICATION|BURSA|PERICARDIAL|DRAIN|GALL|ABDOMEN|ABDOMINAL|BREAST|SEROMA|BILIARY|PLEURAL", SpecimenUpper) &
   !grepl("BRONCH|CSF|CEREBROSPINAL", SpecimenUpper)), "SterileFluid",
  grepl("STOOL|FAECES|RECTAL|GASTRIC", SpecimenUpper), "Gastrointestinal",
  grepl("TIP|LINE|CATHETER", SpecimenUpper) & !grepl("URINE", SpecimenUpper), "Device_LineTip",
  grepl("EAR", SpecimenUpper), "Ear",
  grepl("BONE|ACETABULAR|FEMORAL|TIBIAL", SpecimenUpper), "Bone",
  (SpecimenUpper == "SWAB") |
  grepl("^SWAB;OTHER|SWAB OTHER|SWAB:UNSPECIFIED|SWAB;UNSPECIFIED|SWAB;UNKNOWN|SWAB; Other|SWAB;OTHER", SpecimenUpper), "Skin_SoftTissue",
  (SpecimenUpper == "FLUID"), "SterileFluid",
  default = "Other"
)]
pathogenic_full[, SpecimenUpper := NULL]

matched_full <- pathogenic_full[nbc_df,
  on = .(ClusterID, CollectionDateTime >= adm_minus_28d, CollectionDateTime <= AdmissionDate),
  .(AdmissionID = i.AdmissionID, SourceCategory = x.SourceCategory,
    AccessionNumber = x.AccessionNumber,
    IsGramPos_ext = x.IsGramPos_ext, IsGramNeg_ext = x.IsGramNeg_ext),
  nomatch = 0L, allow.cartesian = TRUE]

other_agg_full <- matched_full[, .(
  Prior_Urine_Pos          = as.integer(any(SourceCategory == "Urine")),
  Prior_Resp_Pos           = as.integer(any(SourceCategory == "Respiratory")),
  Prior_SkinSoftTissue_Pos = as.integer(any(SourceCategory == "Skin_SoftTissue")),
  Prior_SterileFluid_Pos   = as.integer(any(SourceCategory == "SterileFluid")),
  Prior_CSF_Pos            = as.integer(any(SourceCategory == "CSF")),
  Prior_GI_Pos             = as.integer(any(SourceCategory == "Gastrointestinal")),
  Prior_Device_Pos         = as.integer(any(SourceCategory == "Device_LineTip")),
  Prior_Genital_Pos        = as.integer(any(SourceCategory == "Genital")),
  Prior_Ear_Pos            = as.integer(any(SourceCategory == "Ear")),
  Prior_Eye_Pos            = as.integer(any(SourceCategory == "Eye")),
  Prior_Bone_Pos           = as.integer(any(SourceCategory == "Bone")),
  Prior_Other_Pos          = as.integer(any(SourceCategory == "Other")),
  Prior_Total_Other_Cultures_Pos = uniqueN(AccessionNumber),
  Prior_GramPos_Pos        = as.integer(any(IsGramPos_ext)),
  Prior_GramNeg_Pos        = as.integer(any(IsGramNeg_ext))
), by = AdmissionID]

other_agg_full[, Any_Prior_Other_Pos := as.integer(Prior_Total_Other_Cultures_Pos > 0)]

nbc_df <- merge(nbc_df, other_agg_full, by = "AdmissionID", all.x = TRUE)
for (col in c("Any_Prior_Other_Pos","Prior_Urine_Pos","Prior_Resp_Pos",
              "Prior_SkinSoftTissue_Pos","Prior_SterileFluid_Pos",
              "Prior_CSF_Pos","Prior_GI_Pos","Prior_Device_Pos","Prior_Genital_Pos",
              "Prior_Ear_Pos","Prior_Eye_Pos","Prior_Bone_Pos","Prior_Other_Pos",
              "Prior_Total_Other_Cultures_Pos","Prior_GramPos_Pos","Prior_GramNeg_Pos")) {
  nbc_df[is.na(get(col)), (col) := 0L]
}
nbc_df[, adm_minus_28d := NULL]

rm(micro_df, other_micro_full, pathogenic_full, matched_full)

# ============================================================
# 4. SAVE OUTPUT
# ============================================================

fwrite(nbc_df, "Output Files/CSVs/NonBC_Admissions/nbc_with_other_micro.csv")
