# Contributing

This is a research analytics codebase combining Python and R. The conventions below are
inferred from the existing scripts; follow them to keep results reproducible.

## Developer setup

### Python
```zsh
python3 -m venv .venv
source .venv/bin/activate
pip install pandas numpy lifelines statsmodels scipy scikit-learn xgboost matplotlib seaborn
```
Python 3.9+ is recommended (the code uses `from __future__ import annotations` and modern
type hints). XGBoost must be >= 1.6 for the AFT objective.

### R
Install the packages used across the R scripts:
```r
install.packages(c("dplyr", "readr", "tidyr", "purrr", "stringr", "data.table",
                   "openxlsx", "broom", "lubridate", "xgboost", "SHAPforxgboost"))
```

### Data
Raw GENIE / TCGA / MSK cBioPortal exports are not in the repo. ETL scripts read from a local
data root and write canonical analytic frames to `data/processed/`:
`extracted_variables_<cohort>_data.csv` and `extracted_variables_<cohort>_top_genes.txt`.
Generate these before running any Aim analysis. Never commit raw or identifiable patient data;
use the de-identification script and keep private data outside the tracked tree.

## Running and testing

The pipeline is run stage by stage:
1. Pull / extract (R): `Rscript "src/data collection and processing/pull_genie_data.R"`
2. Harmonize (Python): run the per-cohort `harmonize_*.py` scripts.
3. Enrich pathways (R): `add_pathways_genie_bpc.R` / `enrich_harmonized.py`.
4. Analyze (Python): run the Aim 1/2/3 scripts; most take a `--cohort {genie,tcga,msk18}`
   argument and write to `src/data reports/` and `manuscript components/`.

The survival + XGBoost pipeline has its own CLI documented in
`src/modeling/README_survival_AFT_pipeline.md` (e.g. `--xlsx`, `--time-col`, `--event-col`,
`--strata`, `--topn-genes`). There is no automated test suite; validate changes by re-running
the relevant stage and confirming output tables/figures are sane (row counts, event counts,
p-values, concordance indices printed to stdout).

## Coding guidelines

- **Single source of truth for cohorts.** Add or modify cohorts via `CohortSpec` in
  `src/modeling/_lib.py`; use `load_cohort()` / `load_top_genes()` rather than re-reading CSVs.
- **Respect the canonical schema.** Keep column names consistent with `harmonization_spec.md`.
  If a source lacks a variable, create the column as NA rather than dropping it. Never invent
  values; record mapping decisions in a side-car file.
- **Headless plotting.** Set `matplotlib.use("Agg")` before importing `pyplot` (as existing
  scripts do) so figures render without a display.
- **Reproducible stats.** Apply multiple-testing correction where the existing code does,
  use the small-cohort guards (`drop_low_variance`, `MIN_EVENTS`), and keep Cox PH diagnostics
  (Schoenfeld) alongside fits.
- **Module docstrings.** Start each script with a docstring describing inputs, outputs, and
  the analytic step, matching the existing style.
- **Style.** Snake_case for Python, tidyverse style for R, absolute/explicit output paths,
  and informative `print()` progress lines.
- **Text:** use single hyphens (-), not em dashes, in code comments and docs.

Commit only code and generated result tables/figures intended for the manuscript - never
raw patient-level data.
