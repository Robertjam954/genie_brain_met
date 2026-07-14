# Architecture

This project is a reproducible clinical-genomics research pipeline that studies
breast-cancer brain metastasis across three harmonized cancer cohorts: GENIE BPC
BRCA (`genie`), TCGA (`tcga`), and Breast MSK 2018 (`msk18`). It is organized as a
linear data-engineering and statistical-modeling pipeline rather than a deployed
service: raw public cBioPortal / GENIE datasets are pulled, harmonized to a single
canonical schema, then run through descriptive, regression, survival, and machine-learning
analyses whose outputs feed a manuscript.

## High-level flow

```
raw data (cBioPortal / GENIE)
  -> pull / extract        (R)
  -> harmonize to schema   (Python + R)
  -> enrich (pathways)     (R)
  -> analytic CSVs         data/processed/extracted_variables_<cohort>_data.csv
  -> Aim 1 / 2 / 3 analyses (Python)
  -> tables + figures      src/data reports/, manuscript components/
```

## Main components

### 1. Data collection and processing (`src/data collection and processing/`)
- `pull_genie_data.R` - pulls GENIE / TCGA / MSK IMPACT study folders from a local
  data root, joins HUGO gene annotations, and emits extracted variable workbooks.
- `harmonize_genie.py`, `harmonize_tcga.py`, `harmonize_breast_msk_2018.py` - per-cohort
  ETL that builds `extracted_variables_<cohort>_data.csv`. Each aggregates the mutation
  MAF to a per-sample summary (non-silent filter, mutation counts, alt-count max), recodes
  clinical variables, defines the brain-met cohort, derives top-mutated genes, and adds
  time-to-event columns.
- `harmonization_spec.md` - the canonical contract: one row per `SAMPLE_ID`, exact column
  names, source-column candidates, and recoding rules (race/ethnicity/grade/stage/receptor
  subtype, organ-specific metastasis flags, OS / PFS / time-to-brain-met endpoints). This
  spec is what makes the three cohorts conformable.
- `add_pathways_genie_bpc.R` / `enrich_harmonized.py` - attach Sanchez-Vega oncogenic
  pathway indicators (`pathway_RTK/RAS`, `pathway_PI3K`, `pathway_p53`, etc.).
- `missing_analysis_table.py`, `create_deidentified_dataset.R`, date-handling templates -
  supporting QC and de-identification.

### 2. Exploratory data analysis (`src/exploratory data analysis/`)
- `gene_prevtable_oncoprint_forest.py` - Aim 1 core: 2x2 gene/pathway prevalence between
  brain-met and non-brain-met groups, Fisher exact tests, prevalence ratios with
  log-normal CIs, forest plots, and oncoprints.
- `univariate and multivariate regression df and tables.py` - logistic regression of
  `any_brain_met` on top genes and pathways (statsmodels), with multiple-testing correction.
- `aim1_top10_pq_table.py`, R correlation/regression scripts, and the per-cohort
  `*_retrieved_data/` and `msk_retrieved_data/` extracted CSV caches.

### 3. Modeling (`src/modeling/`)
- `_lib.py` - shared backbone. Defines the three `CohortSpec`s, loads analytic CSVs from
  `data/processed/`, reads top-gene lists, builds output dirs, and `prep_covariates()` /
  `drop_low_variance()` for one-hot encoding and small-cohort safety.
- `proportional_hazard_afttest.py` - Aim 2: time to brain metastasis (no-CNS-at-dx cohort).
  KM + log-rank, Cox PH with scaled Schoenfeld diagnostics, and AFT distribution selection
  (Weibull / lognormal / log-logistic by AIC).
- `aim3_os.py` - Aim 3: overall survival in the brain-mets-ever cohort, same KM/Cox/AFT
  machinery with a minimum-events guard (`MIN_EVENTS = 25`).
- XGBoost AFT scripts (`xgb_aft_preprocessing_..._evaluate.py`,
  `xgb_aft_feature_..._importance.py`, `xgb_aft_shap_feature_importance.py`) - survival
  ML using XGBoost's Accelerated Failure Time objective, concordance-index evaluation, and
  SHAP-based feature importance. `shap analysis and plot generation.R` reproduces SHAP
  plots via `SHAPforxgboost`.
- `survival_and_xgb_analysis.ipynb` + `README_survival_AFT_pipeline.md` document the
  combined KM / Cox / Fine-Gray / XGBoost-AFT run and its CLI flags.

### 4. Notebooks and reporting
- `notebooks/` - descriptive, OS, and time-to-brain-met plotting notebooks, retrieval /
  metadata recoding, and a Cox-vs-AFT diagnostic notebook, plus a dataset-metadata template.
- `src/data reports/genie/aim1/` - generated CSV result tables (top-10 gene, univariate /
  multivariate logistic, pathway prevalence, gene-vs-pathway comparison).
- `manuscript components/genie/aim1/` - publication PDFs/figures mirroring those tables.
- `docs/`, `references/` - manuscript structure, revision plans, and source-paper PDFs.

## Data flow and contracts

The hinge of the architecture is the canonical analytic frame
`data/processed/extracted_variables_<cohort>_data.csv` plus
`extracted_variables_<cohort>_top_genes.txt`. ETL writes them; `_lib.load_cohort()` and
`load_top_genes()` read them. Because all three cohorts share canonical columns
(`any_brain_met`, `brain_met_at_dx`, `OS_months`, `os_status_bin`, `top5_any_mutated`,
`receptor_primary_cat`, `grade_ord`, the `pathway_*` set), the same Aim 1/2/3 scripts run
unchanged across cohorts by passing a cohort key.

## Key technologies

- Python: pandas, numpy, lifelines (KM, Cox PH, Schoenfeld, AFT fitters, Fine-Gray),
  statsmodels, scipy, scikit-learn, xgboost (AFT objective), matplotlib, seaborn.
- R: tidyverse (dplyr, readr, tidyr, purrr, stringr), data.table, openxlsx, broom,
  lubridate, xgboost, SHAPforxgboost.
- Data sources: GENIE BPC, TCGA, and MSK cBioPortal study exports (MAF mutation files plus
  clinical tables).

Auxiliary top-level folders (`self-documenting-ai-agent/`, `claude-md-memory-workflow/`,
`context-engineering-workflow.md`) are AI-engineering workflow tutorials kept alongside the
research code and are not part of the analytic pipeline.
