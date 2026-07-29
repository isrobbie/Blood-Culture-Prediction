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
│   └── main_nbc_pipeline      # Cohort: admissions without a culture (NBC_1 .. NBC_8)
├── 02_Models                  # XGBoost, BAM and GramPredict model training
├── 03_Evaluation              # Bootstrap CIs, subgroup/strict-cohort checks, model comparison
├── 04_Report                  # Figure-generating scripts (ROC/calibration, SHAP, DCA/reallocation)
└── 05_Deployment
    └── external_validation_package  # Self-contained package to apply the trained model to a new site
```

Each pipeline stage reads the previous stage's output and writes to `Output Files/CSVs/...`; run scripts in numeric order from the repo root. Paths in scripts are relative to the repo root.

## Data

Patient-level data (IORD, an NHS research database subject to governance approval) are not included and cannot be made public. Scripts here operate on de-identified extracts with the column names referenced throughout the code; see the manuscript Methods and `05_Deployment/external_validation_package/Data_Dictionary/` for feature definitions. Data are available on reasonable request via IORD, subject to approval (see manuscript for details).

Free-text antibiotic indications are classified into infection-source categories using a fine-tuned Bio+ClinicalBERT model (`05_Deployment/external_validation_package/Use_Code/Classify_Indications_BERT.py`) rather than being distributed as text, to avoid releasing clinical free text.

## Setup

R 4.3+ with the packages below (also listed in `install_packages.R`):

```r
source("install_packages.R")
```

Model deployment/inference additionally needs Python with `torch`, `transformers`, `datasets`, `tqdm`, `pandas` (see header of `Classify_Indications_BERT.py`).

## Applying the trained model elsewhere

`05_Deployment/external_validation_package/` is self-contained: it holds the trained XGBoost model, Platt calibration coefficients, feature data dictionary and the code used to apply the model and analyse results (as done for the UCLH external validation). Start from `Use_Code/Apply_XGBoost.r`.

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

## License

MIT (see `LICENSE`). This covers the code only; the underlying patient data are governed separately by IORD.
