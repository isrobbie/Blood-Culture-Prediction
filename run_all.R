# Runs the full pipeline in order from the repo root: Rscript run_all.R
# Requires the IORD extract in "Input Files/" (see README) and the packages listed there.

stages <- c(
  sprintf("01_Preprocessing/main_bc_pipeline/BC_%d_%s.r", 1:8,
          c("Initial_Extract", "Add_Demographics", "Add_Admissions", "Add_Infection_Source",
            "Add_Vital_Signs", "Add_Blood_Gas", "Add_Blood_Tests", "Add_OtherInfections")),
  sprintf("01_Preprocessing/front_door_pipeline/Trial_BC_%d_%s.r", 5:8,
          c("Add_Vital_Signs", "Add_Blood_Gas", "Add_Blood_Tests", "Add_OtherInfections")),
  sprintf("01_Preprocessing/main_nbc_pipeline/NBC_%d_%s.r", 1:8,
          c("Define_Cohort", "Add_Demographics", "Add_Prior_BC", "Add_Infection_Source",
            "Add_Vital_Signs", "Add_Blood_Gas", "Add_Blood_Tests", "Add_Other_Micro")),
  "02_Models/Do_XGBoost_FullCohort.r",
  "02_Models/Do_BAM_FullCohort.r",
  "02_Models/Do_XGBoost_GramPredictFullCohort.r",
  "02_Models/NBC_9_Apply_Model.r",
  list.files("03_Evaluation", pattern = "\\.r$", full.names = TRUE),
  "04_Report/Generate_Figures.r"
)

for (s in stages) {
  message("Running ", s)
  source(s, local = new.env())
}
