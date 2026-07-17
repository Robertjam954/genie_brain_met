# Architecture and Data Flow

This page describes the structure of the pipeline, how data moves through it, and where
each component lives.

---

## High-level flow

```
raw private cBioPortal / GENIE / GDC exports
  |
  v
[1] Pull / extract (R)
    pull_genie_data.R
    -> data/extracted_variables_of_interest.xlsx
  |
  v
[2] Harmonize to canonical schema (Python, per cohort)
    harmonize_genie.py
    harmonize_tcga.py
    harmonize_breast_msk_2018.py
    -> src/exploratory data analysis/extracted_variables_<cohort>_data.csv
    -> src/exploratory data analysis/extracted_variables_<cohort>_top_genes.txt
    -> src/exploratory data analysis/extracted_variables_<cohort>_gene_prev_brain_met.csv
    -> src/exploratory data analysis/extracted_variables_<cohort>_dictionary.csv
  |
  v
[3] Enrich (Python)
    enrich_harmonized.py
    -> rewrites the cohort CSVs in place (alias columns, competing-event endpoints,
       GENIE regimen flags)
  |
  v
[4] Stage analytic frames (manual copy)
    data/processed/extracted_variables_<cohort>_data.csv   <-- _lib.load_cohort() reads here
    data/processed/extracted_variables_<cohort>_top_genes.txt
  |
  v
[5] Aim 1 - association analysis (Python)
    gene_prevtable_oncoprint_forest.py --cohort <cohort>
    aim1_top10_pq_table.py --cohort <cohort>
    univariate and multivariate regression df and tables.py --cohort <cohort>
    -> src/modeling/<cohort>/aim1/ (CSV tables)
    -> manuscript components/<cohort>/aim1/ (figures/PDFs)
  |
  v
[6] Aim 2 - time to brain met (Python)
    proportional_hazard_afttest.py --cohort <cohort>
    -> src/modeling/<cohort>/aim2/
    -> manuscript components/<cohort>/aim2/
  |
  v
[7] Aim 3 - overall survival (Python)
    aim3_os.py --cohort <cohort>
    -> src/modeling/<cohort>/aim3/
    -> manuscript components/<cohort>/aim3/
  |
  v
[8] ML survival benchmark (Python + R)
    xgb_aft_preprocessing_...evaluate.py
    xgb_aft_feature_...importance.py
    xgb_aft_shap_feature_importance.py
    shap analysis and plot generation.R
    -> reports/figures/ml_benchmark/
```

---

## Important location nuance

The harmonize scripts **write** their outputs to `src/exploratory data analysis/`, but
`_lib.load_cohort()` **reads** from `data/processed/`. These are different directories.
`data/processed/` is not present in the repo. The analytic CSVs must be copied there
before running any Aim analysis.

---

## Repository layout

```
src/
  data collection and processing/   Raw pull (R) + harmonization ETL (Python) + schema spec
  exploratory data analysis/        Aim 1 association analyses; also holds the harmonized
                                    extracted_variables_<cohort>_* CSV caches
  modeling/                         _lib.py shared backbone + Aim 2/3 survival + XGBoost-AFT
  data reports/                     Generated CSV result tables (genie/aim1 currently committed)
notebooks/                          Descriptive / OS / time-to-brain-met plotting notebooks
reports/figures/                    Rendered PNG figures (aim1, aim2, aim3, ml_benchmark)
manuscript components/              Publication-figure PDFs/PNGs (genie/aim1 currently committed)
docs/                               results_dashboard.html + manuscript-planning PDFs + plans
references/                         Source-paper PDFs and study materials
data/                               Raw / master workbooks (git-ignored)
archive/                            Archived AI-workflow tutorial content
wiki/                               This wiki
```

---

## Key components

### `src/modeling/_lib.py` - Shared backbone

Central configuration and utility module shared by all Aim scripts.

- Defines `CohortSpec` for each of the three cohorts (`genie`, `tcga`, `msk18`)
- `load_cohort(cohort)` - reads `data/processed/extracted_variables_<cohort>_data.csv`
- `load_top_genes(cohort)` - reads the top-gene list for a cohort
- `prep_covariates()` - one-hot encoding of clinical covariates
- `drop_low_variance()` - drops columns with near-zero variance (small-cohort safety)
- `ensure_dirs()` - creates output directories as needed
- `PATHWAY_COLS` - the ten Sanchez-Vega oncogenic pathway column names

### `src/data collection and processing/harmonization_spec.md` - Schema contract

Documents the canonical analytic frame that all three cohorts must conform to. Defines:
- One row per `SAMPLE_ID`
- Exact column names and types
- Recoding rules for race, ethnicity, grade, stage, receptor subtype
- Organ-specific metastasis flags
- OS / PFS / time-to-brain-met endpoints
- Brain-met cohort derivation logic
- Mutation aggregation from MAF
- gnomeR gene-binary matrix columns
- Ten Sanchez-Vega pathway columns

### Cohort keys

Scripts that support multiple cohorts accept `--cohort {genie,tcga,msk18,all}`.

| Key | Cohort | Data file |
|-----|--------|-----------|
| `genie` | GENIE BPC BRCA | `extracted_variables_genie_data.csv` |
| `tcga` | TCGA BRCA | `extracted_variables_tcga_data.csv` |
| `msk18` | Breast MSK 2018 | `extracted_variables_breast_msk_2018_data.csv` |

---

## Technologies

| Language | Key libraries |
|----------|--------------|
| Python | pandas, numpy, lifelines (KM, Cox PH, AFT, Fine-Gray), statsmodels, scipy, scikit-learn, xgboost (AFT), matplotlib, seaborn, shap |
| R | dplyr, readr, tidyr, purrr, stringr, data.table, openxlsx, broom, lubridate, xgboost, SHAPforxgboost |
