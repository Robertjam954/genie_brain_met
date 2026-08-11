# GENIE BPC - Breast-Cancer Brain Metastasis

Genes and clinical risk factors shape both the likelihood and the timing of breast-cancer
brain metastasis, and some carry far more weight than others. This project quantifies that
**variable impact** in the **GENIE BPC BRCA** cohort: it pairs classical association and
survival statistics with nonlinear time-to-event models and modern machine learning
(**XGBoost with an Accelerated Failure Time objective**) to rank which gene mutations,
oncogenic pathways, and clinical features most drive brain metastasis and survival after it.

> **Status / verification note.** This repository is an active research work-in-progress.
> The narrative results in `docs/results_dashboard.html` and the methods write-ups have
> **not been independently re-verified** and may contain errors; treat every specific
> number (cohort sizes, prevalences, p/q-values, hazard ratios, concordance indices) as
> provisional pending author review. This README describes the **pipeline structure**
> (what each script reads, writes, and does) rather than asserting any scientific finding.

## Analytic aims

| Aim | Question | Methods (as implemented) |
|-----|----------|--------------------------|
| **1 - Association** | Which genes / oncogenic pathways are enriched in brain-met patients? | 2x2 contingency tables, Fisher exact, prevalence ratios with log-normal CIs, multiple-testing correction, logistic regression, forest plots, oncoprints |
| **2 - Time to brain met** | Time from diagnosis to brain metastasis, in the no-CNS-at-diagnosis cohort | Kaplan-Meier + log-rank, Cox PH + scaled Schoenfeld diagnostics, AFT distribution selection by AIC |
| **3 - Overall survival** | Overall survival in the brain-mets-ever cohort | Same KM / Cox / AFT machinery, with a minimum-events guard |
| **ML benchmark** | Predictive survival modeling | XGBoost with the AFT objective, concordance-index evaluation, SHAP-based feature importance (Python + R) |

## Repository layout

```
src/
  data collection and processing/   Raw pull (R) + data-prep ETL (Python) + schema spec
  exploratory data analysis/        Aim 1 association analyses; also holds the prepared
                                    extracted_variables_genie_* CSV caches
  modeling/                         _lib.py shared backbone + Aim 2/3 survival + XGBoost-AFT
  data reports/                     Generated CSV result tables (currently genie/aim1 only)
notebooks/                          Descriptive / OS / time-to-brain-met plotting notebooks
reports/figures/                    Rendered PNG figures (aim1, aim2, aim3, ml_benchmark)
manuscript components/             Publication-figure PDFs/PNGs (currently genie/aim1 only)
docs/                               results_dashboard.html + manuscript-planning PDFs +
                                    plans + CBBio_workflow.md (driver-gene workflow notes)
references/                         Source-paper PDFs, study materials, and REFERENCES.md
data/                              Raw / master workbooks (git-ignored; see Data below)
```

See [`PRODUCT.md`](PRODUCT.md) and [`ARCHITECTURE.md`](ARCHITECTURE.md) for the product and
architecture overview, [`CONTRIBUTING.md`](CONTRIBUTING.md) for developer conventions,
[`references/REFERENCES.md`](references/REFERENCES.md) for the full bibliography (source
papers, the CBBio driver-gene method, and the papers/tools it cites), and
`docs/results_dashboard.html` for the (unverified) figure dashboard. A complete
file-by-file inventory is in [Repository file inventory](#repository-file-inventory) below.

## Prerequisites and setup

The Python survival + ML stack is pinned in [`requirements.txt`](requirements.txt) (pip) and
[`environment.yml`](environment.yml) (conda), a coherent Python 3.9-3.11 constellation. There
is no R `renv.lock` yet; the R packages below are still installed manually.

### Python

```zsh
# pip
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# or conda / mamba (env name matches the notebook kernel `genie-brainmet`)
conda env create -f environment.yml && conda activate genie-brainmet
```

`requirements.txt` covers everything the scripts and notebooks import: pandas, numpy,
scipy, scikit-learn, statsmodels, lifelines, scikit-survival, xgboost, lightgbm, optuna,
shap, matplotlib, seaborn, and openpyxl. Python 3.9+ is required (`_lib.py` uses
`from __future__ import annotations` and modern type hints); XGBoost is pinned >= 1.6 for
the AFT objective. If you are extending the author's existing `tcga-analysis` conda env
rather than building a fresh venv, reconcile the pins against that env first
(see `docs/plans/finish-project-run-notebooks-plan.md`).

### R

```r
install.packages(c("dplyr", "readr", "tidyr", "purrr", "stringr", "data.table",
                   "openxlsx", "broom", "lubridate", "xgboost", "SHAPforxgboost"))
```

### Data

Raw and derived data are **git-ignored** (see `.gitignore`: `data/`, `*.csv`, `*.xlsx`,
`**/extracted_variables_*`, the `*_retrieved_data/` folders). This project operates on
de-identified public datasets for research; it is **not** a clinical tool. Never commit raw
or identifiable patient data.

The pipeline reads raw cohort exports from **local, machine-specific absolute paths**
(the R pull scripts point at `/Users/robertjames/loc/data private/...`; the notebooks
hard-code a `PROJ` root under the author's iCloud `Documents`). You will need to adjust
those paths for your own machine.

> **Prerequisite gap.** `_lib.load_cohort()` expects the prepared analytic frame at
> `data/processed/extracted_variables_genie_data.csv`, but this directory is **not
> present in the repo**, and the harmonize scripts actually *write* their outputs into
> `src/exploratory data analysis/`. In practice the analytic CSVs must be placed under
> `data/processed/` for the modeling scripts to find them. Do not assume the ETL has been
> run end-to-end here; verify the inputs exist before running any analysis stage.

## Pipeline run order

The pipeline runs stage by stage. Paths below are relative to the repo root; several scripts
use absolute, machine-specific paths internally and will need editing.

1. **Raw pull (R)** - read raw cBioPortal / GENIE / GDC exports and gene annotations, emit a
   consolidated workbook.
   ```zsh
   Rscript "src/data collection and processing/pull_genie_data.R"
   ```
   Input: private cBioPortal study folders + `data/hugo_symbols.xlsx`.
   Output: `data/extracted_variables_of_interest.xlsx`.
   (`extract_variables_of_interest.R` appears to be a near-duplicate of this script.)

2. **Build the analytic frame (Python)** - `harmonize_genie.py` reads intermediate master
   CSVs from `data/processed/` (e.g. `genie_bpc_v1_sample_master_full.csv`,
   `genie_bpc_v1_mutations.csv`) - which are produced upstream (partly in R / `gnomeR`) and
   are **not** generated by any script visible in this repo - aggregates the mutation MAF to
   per-sample features, recodes the clinical variables, defines the brain-met cohort, and
   derives the top-mutated genes and time-to-event endpoints. It writes the analytic frame
   plus a top-genes list and a recoding dictionary.
   ```zsh
   python "src/data collection and processing/harmonize_genie.py"
   ```
   Outputs (written to `src/exploratory data analysis/`):
   `extracted_variables_genie_data.csv`, `_top_genes.txt`,
   `_gene_prev_brain_met.csv`, `_dictionary.csv`.

3. **Enrich (Python)** - `enrich_harmonized.py` rewrites the three
   `extracted_variables_genie_data.csv` frames in place, adding long-form alias columns,
   competing-event time-to-event columns for Aim 2, and (GENIE only) regimen-derived
   treatment flags read from `data/processed/genie_bpc_v1_regimens.csv`.
   ```zsh
   python "src/data collection and processing/enrich_harmonized.py"
   ```
   Note: `add_pathways_genie_bpc.R` is the upstream R step that builds the gene-binary
   matrix and the 10 Sanchez-Vega `pathway_*` columns (via gnomeR, per
   `harmonization_spec.md` sections 13-14) and writes
   `data/processed/genie_bpc_v1_sample_master_full.csv` - the clinical + gene_binary +
   pathways master that `harmonize_genie.py` reads. Run it before Stage 2 if that master
   does not already carry pathway columns.

4. **Make the analytic frame discoverable to modeling.** Ensure the prepared
   `extracted_variables_genie_data.csv` and `_top_genes.txt` files are present under
   `data/processed/`, because `_lib.load_cohort()` / `load_top_genes()` read from there.

5. **Aim 1 - association (Python).** Scripts take a `--cohort genie` argument.
   ```zsh
   python "src/exploratory data analysis/gene_prevtable_oncoprint_forest.py" --cohort genie
   python "src/exploratory data analysis/aim1_top10_pq_table.py" --cohort genie
   python "src/exploratory data analysis/univariate and multivariate regression df and tables.py" --cohort genie
   ```
   CSV tables are written under `src/modeling/genie/aim1/`; rendered tables/figures under
   `manuscript components/genie/aim1/`.

6. **Aim 2 / Aim 3 - survival (Python).**
   ```zsh
   python "src/modeling/proportional_hazard_afttest.py" --cohort genie   # Aim 2, no-CNS-at-dx cohort
   python "src/modeling/aim3_os.py" --cohort genie                       # Aim 3, brain-mets-ever cohort
   ```
   Both use `_lib.prep_covariates` / `drop_low_variance` and a minimum-events guard, and
   write to `src/modeling/genie/aimN/` (CSVs) and `manuscript components/genie/aimN/`
   (figures).

7. **ML survival benchmark (Python + R).** The XGBoost-AFT trainer
   (`xgb_aft_preprocessing_feature_constuction_train_validate_evaluate.py`) and its SHAP
   explainer (`xgb_aft_shap_feature_importance.py`) are `_lib`-based and take
   `--cohort`/`--aim`; `shap analysis and plot generation.R` reproduces SHAP plots via
   `SHAPforxgboost`. `stratifiedKM_CoxFG_feature_prep_AFT.py` is an older standalone driver
   still on the legacy schema. See `src/modeling/README_survival_AFT_pipeline.md` for the
   full breakdown of what is ported vs legacy.

There is no automated test suite. Validate changes by re-running the relevant stage and
sanity-checking the printed diagnostics and output tables/figures.

## Repository file inventory

Every tracked file needed to understand and run the project, grouped by directory. Raw and
derived data are git-ignored (see [Data](#data)) and are not listed here.

### Prep docs (repo root)

- `README.md` - this file: pipeline structure, setup, run order, inventory.
- `PRODUCT.md` - product/aims overview.
- `ARCHITECTURE.md` - pipeline structure, inputs/outputs, data contracts.
- `CONTRIBUTING.md` - developer conventions and setup.
- `STATUS.md` - live completion checklist.
- `requirements.txt`, `environment.yml` - pinned Python (pip) and conda environments.
- `finish-project-run-notebooks-plan.md` - plan for running the remaining notebooks.
- `.gitignore` - excludes `data/`, CSV/XLSX, and extracted-variable caches.

### `src/data collection and processing/` - ETL

- `pull_genie_data.R` - pull raw GENIE BPC study folders + HUGO annotations -> `data/extracted_variables_of_interest.xlsx`.
- `extract_variables_of_interest.R` - near-duplicate of `pull_genie_data.R`.
- `harmonize_genie.py` - build the GENIE analytic frame (mutation aggregation, clinical recoding, cohort/endpoints).
- `harmonize_tcga.py`, `harmonize_breast_msk_2018.py` - harmonizers for the TCGA and MSK-2018 cohorts.
- `enrich_harmonized.py` - add alias, competing-event, and regimen-derived columns in place.
- `add_pathways_genie_bpc.R` - gene-binary + Sanchez-Vega pathway builder (currently a 0-byte stub).
- `harmonization_spec.md` - the canonical analytic-schema contract (one row per `SAMPLE_ID`).
- `missing_analysis_table.py` - missingness audit CSV + rendered table (`--cohort`).
- `create_deidentified_dataset.R` - de-identification (references offline Windows paths).
- `date_handling_template.R` - date-handling helper template.
- `metadata template.pdf` - dataset-metadata template.

### `src/exploratory data analysis/` - Aim 1 association

- `gene_prevtable_oncoprint_forest.py` - gene/pathway prevalence, Fisher exact, prevalence ratios, forest plots, oncoprints.
- `univariate and multivariate regression df and tables.py` - logistic regression of `any_brain_met` (statsmodels).
- `aim1_top10_pq_table.py` - top-10 gene p/q table.
- `calculatecorrelation_regressions.R` - correlation / regression (R).
- `genomic exploratory analysis with tables and plots.rtf`, `retrieve tcga genie data.rtf` - exploratory notes.

### `src/modeling/` - survival + ML

- `_lib.py` - shared backbone: `CohortSpec`, `load_cohort()`, `load_top_genes()`, `prep_covariates()`, `drop_low_variance()`.
- `proportional_hazard_afttest.py` - Aim 2 (time to brain met): KM + log-rank, Cox PH + Schoenfeld, AFT by AIC.
- `aim3_os.py` - Aim 3 (overall survival), same machinery with a minimum-events guard.
- `xgb_aft_preprocessing_feature_constuction_train_validate_evaluate.py` - XGBoost-AFT trainer.
- `xgb_aft_shap_feature_importance.py` - SHAP feature importance for the XGBoost-AFT model.
- `stratifiedKM_CoxFG_feature_prep_AFT.py` - older standalone KM/Cox/Fine-Gray/AFT driver (legacy schema).
- `shap analysis and plot generation.R` - SHAP plots via `SHAPforxgboost`.
- `survival_and_xgb_analysis.ipynb` - combined KM/Cox/Fine-Gray/XGBoost-AFT notebook.
- `README_survival_AFT_pipeline.md` - modeling pipeline notes (ported vs legacy).

### `notebooks/`

- `descriptiveplots.ipynb`, `osplots.ipynb`, `timebrainmettplots.ipynb`, `coxvsaft_error_os.ipynb` - plotting/diagnostic notebooks.
- `retrieval_metadata_recode.ipynb` - retrieval + metadata recoding.
- `data merge and preprocessing markdowm.md`, `project notebook template.md` - notebook docs/templates.
- `Workflow_CBBio.pdf` - candidate driver-gene workflow diagram (see below).

### CBBio driver-gene workflow

- `src/Workflow_CBBio.pdf`, `notebooks/Workflow_CBBio.pdf` - the CBBio workflow diagram.
- `docs/CBBio_workflow.md` - run notes reproduced from the CBBio v0.1 package (tools, commands, ensemble).
- `references/REFERENCES.md` §2-§3 - the CBBio method, its cited papers, and the caller-tool sources.

### `docs/`, `references/`, `reports/`, `manuscript components/`

- `docs/results_dashboard.html`, `docs/index.html` - the (unverified) figure dashboard.
- `docs/plans/finish-project-run-notebooks-plan.md` - notebook-run plan.
- `docs/*.pdf` - manuscript revision plan and intro/methods structure edits.
- `references/*.pdf` - source-paper PDFs and study materials (indexed in `references/REFERENCES.md` §1).
- `references/REFERENCES.md` - full bibliography.
- `reports/figures/{aim1,aim2,aim3,ml_benchmark}/*.png` - rendered figures.
- `manuscript components/genie/aim1/*.pdf` - publication PDFs for the GENIE Aim 1 tables.

### `wiki/`

- `Home.md`, `Getting-Started.md`, `Analysis-Aims.md`, `Architecture-and-Data-Flow.md`,
  `Data-Schema-and-Harmonization.md`, `Pipeline-Run-Guide.md`, `Known-Issues-and-Gaps.md`,
  `Contributing.md` - project wiki.

### `archive/` - auxiliary content (not part of the analytic pipeline)

The `archive/self-documenting-ai-agent/`, `archive/claude-md-memory-workflow/`,
`archive/context-engineering-workflow.md`, and `archive/plan-template.md` files are
AI-engineering workflow scaffolding kept alongside the research code (see
`references/REFERENCES.md` §5). They are unrelated to the brain-metastasis analysis.

## References

Full bibliography: [`references/REFERENCES.md`](references/REFERENCES.md). Key items:

- **AFT survival modeling** - Barnwal et al., "Survival regression with accelerated failure
  time model in XGBoost" (basis for the `xgb_aft_*` scripts).
- **Brain-metastasis genomics** - meta-analysis of breast-cancer brain-metastasis genomics
  (motivates Aim 1).
- **CBBio candidate driver-gene method** (Leila Mirsadeghi) - four driver callers
  (MutSigCV, OncodriveCLUST, OncodriveFM, NetBox) -> per-gene `-log10(p)` features ->
  SVM/ANN/RF ensemble with an algebraic combiner; workflow in
  [`docs/CBBio_workflow.md`](docs/CBBio_workflow.md), citations in `references/REFERENCES.md`
  §2, tool sources in §3.
- **Data** - AACR Project GENIE BPC BRCA (via cBioPortal); Sanchez-Vega oncogenic pathways.
