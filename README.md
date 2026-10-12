# Blood Culture Prediction

Code to reproduce the analysis in:

> Forrest-Hammond R, Gupta R, McVean G, Noursadeghi M, O'Grady J, Samuels THA, Eyre DW.
> **Machine learning models to improve targeting of blood culture testing.**
> medRxiv 2026. https://doi.org/10.1101/2026.07.17.26358320

## Abstract

Bloodstream infections are a major cause of mortality, yet blood cultures have low positivity (<10%) and turnaround times of 24–48 hours. Many are taken from patients at low risk, while some infections are diagnosed late or missed entirely. Here we develop and externally validate machine learning models to improve targeting of blood culture testing using routinely collected clinical and laboratory data. In 294,064 cultures from a multi-site NHS trust (positivity 5.6%), XGBoost discriminated pathogenic positivity in a temporal hold-out test set (AUROC 0.853, 95% CI 0.846–0.860; 0.876 in emergency department patients), with good calibration (slope 1.046), and performance was maintained in external validation of 37,326 encounters from University College London Hospitals (0.847, 0.839–0.856). In a simulated resource-neutral reallocation, replacing the 10,000 lowest-risk cultures with the highest-risk untested admissions yielded 627 additional positives (28.3% relative increase). This model could improve blood culture yield without increasing testing volume, supporting resource-neutral diagnostic stewardship.

## Repo structure

```
.
├── 00_Data                    # Lookup/reference tables (no patient data)
├── 01_Preprocessing
│   ├── main_bc_pipeline       # Cohort: blood cultures taken (BC_1 .. BC_8)
│   ├── front_door_pipeline    # Same cohort with front-door time windows (Trial_BC_5 .. 8)
│   └── main_nbc_pipeline      # Cohort: admissions without a culture (NBC_1 .. NBC_8)
├── 02_Models                  # XGBoost, BAM and GramPredict model training
├── 03_Evaluation              # Bootstrap CIs, subgroup, SHAP, front-door and deduplication analyses
├── 04_Report                  # Generate_Figures.r produces Figures 1-3
└── run_all.R                  # Runs every stage in order
```

## Data

Patient-level data (IORD, an NHS research database subject to governance approval) are not included and cannot be made public. Data are available to researchers whose proposal meets IORD's ethical and governance requirements, under a data access agreement (see manuscript Data sharing statement).

The pipeline expects the IORD extract in `Input Files/IORD CSVs/` (file names as referenced in `BC_1`–`BC_8` and `NBC_1`–`NBC_8`) and the BERT-classified antibiotic indications in `Input Files/bert_classified_indications_update.csv`. Variable definitions are in Supplementary Table 2b.

The trained XGBoost model and Platt calibration are available from the authors on request.

## Setup

R 4.3+ with:
`data.table`, `dplyr`, `tidyr`, `purrr`, `stringr`, `lubridate`, `xgboost`, `mgcv`, `pROC`, `PRROC`, `ggplot2`, `patchwork`, `writexl`, `openxlsx`, `SHAPforxgboost`.

## Reproducing the results

From the repo root, with the IORD extract in place:

```
Rscript run_all.R
```

This runs preprocessing, model training, evaluation and figures in order. Outputs are written to `Output Files/` and models to `Models/`.

| Output | Script |
|---|---|
| Table 2, Figure 1 | `02_Models/Do_XGBoost_FullCohort.r`, `04_Report/Generate_Figures.r` |
| BAM comparison | `02_Models/Do_BAM_FullCohort.r`, `03_Evaluation/Model_Comparison_BAM_vs_XGBoost.r` |
| Figure 2, Supplementary Table 8 (SHAP) | `04_Report/Generate_Figures.r`, `03_Evaluation/Do_XGBoost_SHAP_Consistency.r` |
| Front-door analysis | `01_Preprocessing/front_door_pipeline`, `03_Evaluation/Do_Evaluate_MainModel_on_StrictCohort.r` |
| Figure 3, Table 3 (reallocation) | `02_Models/NBC_9_Apply_Model.r`, `03_Evaluation/NBC_10_Reallocation_Analysis.r` |
| Supplementary Table 9 (subgroups, completeness, confusion matrices) | `03_Evaluation/Do_XGBoost_Supp_Subgroup_Tables.r` |

## Cite

```bibtex
@article{forresthammond2026blood,
  title={Machine learning models to improve targeting of blood culture testing},
  author={Forrest-Hammond, Robbie and Gupta, Rishi and McVean, Gil and Noursadeghi, Mahdad and O'Grady, Justin and Samuels, Thomas HA and Eyre, David W},
  journal={medRxiv},
  year={2026},
  doi={10.1101/2026.07.17.26358320}
}
```
