# Product

> **Status caveat.** This is an active research work-in-progress. The specific results shown
> in `docs/results_dashboard.html` and the methods narratives have not been independently
> re-verified and may contain errors; treat all reported numbers as provisional. This
> document describes intended capabilities and methods, not confirmed findings.

## What this project is

A research analytics pipeline investigating the genomic and clinical drivers of
**breast-cancer brain metastasis**. It harmonizes three public clinical-genomic cancer
cohorts - GENIE BPC BRCA, TCGA, and Breast MSK 2018 - into a single comparable analytic
schema and runs a consistent set of statistical and machine-learning analyses across all
of them. The end product is a set of reproducible tables and figures that feed a scientific
manuscript on which gene mutations, oncogenic pathways, and clinical features are associated
with developing brain metastases and with survival once they occur.

## What it does

The work is organized around three analytic aims, each runnable per cohort:

### Aim 1 - Association with brain metastasis
- Compares prevalence of the top mutated genes and Sanchez-Vega oncogenic pathways between
  patients who develop CNS/brain metastasis and those who do not.
- Computes 2x2 contingency tables, Fisher exact tests, prevalence ratios with log-normal
  95% CIs, and multiple-testing-corrected p/q values.
- Fits univariate and multivariate logistic regressions of `any_brain_met` on genes and
  pathways.
- Produces forest plots, oncoprints, and prevalence comparison tables.

### Aim 2 - Time to brain metastasis
- In the cohort without CNS involvement at diagnosis, models time from diagnosis to brain
  metastasis.
- Kaplan-Meier curves stratified by mutation status with log-rank tests, Cox proportional
  hazards models with scaled Schoenfeld diagnostics, and Accelerated Failure Time (AFT)
  models with distribution selection (Weibull / lognormal / log-logistic by AIC).

### Aim 3 - Overall survival
- In the brain-mets-ever cohort, models overall survival from diagnosis using the same
  KM / log-rank / Cox PH / AFT machinery, with a minimum-event-count guard before fitting.

### Machine-learning survival modeling
- XGBoost with the Accelerated Failure Time objective for survival prediction, evaluated by
  concordance index, with feature-importance and SHAP-based explainability outputs (both
  Python and R).

## Key features

- **Cross-cohort harmonization.** A documented canonical schema (`harmonization_spec.md`)
  defines one row per sample with standardized clinical recodings (race, ethnicity, grade,
  stage, receptor subtype), organ-specific metastasis flags, and OS / PFS / time-to-brain-met
  endpoints, so the same analysis runs identically on GENIE, TCGA, and MSK 2018.
- **Mutation and pathway feature engineering.** Per-sample mutation summaries derived from
  MAF files (non-silent filtering, mutation counts, top-gene indicators) plus Sanchez-Vega
  oncogenic pathway flags.
- **Classical and ML survival analysis side by side.** Cox PH (with assumption diagnostics),
  AFT, and XGBoost-AFT, allowing comparison of interpretable and predictive approaches.
- **Manuscript-oriented outputs.** Analyses emit CSV result tables (Aim scripts write to
  `src/modeling/<cohort>/aimN/`; a GENIE Aim 1 subset is committed under `src/data reports/`)
  and publication figures/PDFs (`manuscript components/`, `reports/figures/`). Output
  coverage in the repo is currently partial - primarily GENIE Aim 1 tables plus rendered
  figures - so not every aim/cohort listed above has committed result files.
- **Centralized configuration.** Cohort definitions, gene lists, and covariate handling are
  centralized in shared code (`src/modeling/_lib.py`) and driven by cohort keys, so the same
  analysis runs across cohorts. Note the repo does not yet pin a Python/R environment
  (no `requirements.txt` / `environment.yml` / `renv.lock`), which is a reproducibility gap.

## Intended users and use cases

- **Cancer genomics / translational researchers** examining biomarkers of brain metastasis
  in breast cancer.
- **Biostatisticians and computational biologists** who need a harmonized, multi-cohort
  framework for survival and association analyses.
- **Manuscript authors** producing the tables and figures for a brain-metastasis publication.

This is a research codebase, not a clinical or production decision-support tool. It operates
on de-identified public datasets and is intended for analysis and publication, not for
guiding individual patient care.
