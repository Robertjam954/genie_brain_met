# Product

> **Status caveat.** This is an active research work-in-progress. The specific results shown
> in `docs/results_dashboard.html` and the methods narratives have not been independently
> re-verified and may contain errors; treat all reported numbers as provisional. This
> document describes intended capabilities and methods, not confirmed findings.

## What this project is

A research analytics pipeline investigating the genomic and clinical drivers of
**breast-cancer brain metastasis** in the **GENIE BPC BRCA** cohort. Genes and clinical risk
factors shape both the likelihood and the timing of brain metastasis, and some carry far more
weight than others; this project quantifies that variable impact by pairing classical
association and survival statistics with nonlinear time-to-event models and modern machine
learning (XGBoost with an Accelerated Failure Time objective). The end product is a set of
reproducible tables and figures that feed a scientific manuscript on which gene mutations,
oncogenic pathways, and clinical features are associated with developing brain metastases and
with survival once they occur.

## What it does

The work is organized around three analytic aims:

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

- **Documented analytic schema.** A written spec (`harmonization_spec.md`) defines one row
  per sample with standardized clinical recodings (race, ethnicity, grade, stage, receptor
  subtype), organ-specific metastasis flags, and OS / PFS / time-to-brain-met endpoints,
  giving every downstream analysis a single, well-defined GENIE BPC frame to read.
- **Mutation and pathway feature engineering.** Per-sample mutation summaries derived from
  MAF files (non-silent filtering, mutation counts, top-gene indicators) plus Sanchez-Vega
  oncogenic pathway flags.
- **Classical and ML survival analysis side by side.** Cox PH (with assumption diagnostics),
  AFT, and XGBoost-AFT, allowing comparison of interpretable and predictive approaches.
- **Manuscript-oriented outputs.** Analyses emit CSV result tables (Aim scripts write to
  `src/modeling/genie/aimN/`; a GENIE Aim 1 subset is committed under `src/data reports/`)
  and publication figures/PDFs (`manuscript components/`, `reports/figures/`). Output
  coverage in the repo is currently partial - primarily GENIE Aim 1 tables plus rendered
  figures - so not every aim listed above has committed result files.
- **Centralized configuration.** The cohort definition, gene lists, and covariate handling
  are centralized in shared code (`src/modeling/_lib.py`), so every aim reads the same frame
  the same way. Note the repo does not yet pin a Python/R environment
  (no `requirements.txt` / `environment.yml` / `renv.lock`), which is a reproducibility gap.

## Intended users and use cases

- **Cancer genomics / translational researchers** examining biomarkers of brain metastasis
  in breast cancer.
- **Biostatisticians and computational biologists** who need a worked framework for survival
  and association analyses on genomic cohort data.
- **Manuscript authors** producing the tables and figures for a brain-metastasis publication.

This is a research codebase, not a clinical or production decision-support tool. It operates
on de-identified public datasets and is intended for analysis and publication, not for
guiding individual patient care.

## Related reference workflow

The repo also carries the CBBio candidate driver-gene workflow (Leila Mirsadeghi;
`notebooks/Workflow_CBBio.pdf`, run notes in `docs/CBBio_workflow.md`) as reference material:
four driver callers (MutSigCV, OncodriveCLUST, OncodriveFM, NetBox) feed one `-log10(p)`
feature each into an SVM/ANN/RF ensemble with an algebraic combiner to prioritize candidate
driver genes. Its method citations and tool sources are in `references/REFERENCES.md` (§2-§3).

## References

See [`references/REFERENCES.md`](references/REFERENCES.md) for the full bibliography: source
papers, the AFT-in-XGBoost basis for the ML modeling, the CBBio driver-gene method and its
cited papers, and the data sources.
