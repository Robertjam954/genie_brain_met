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
Raw GENIE BPC cBioPortal exports are not in the repo (and are git-ignored). The R
pull scripts and the notebooks reference local, machine-specific absolute paths - you must
adjust these. The harmonize scripts read intermediate master CSVs from `data/processed/`
(produced upstream, partly in R / `gnomeR`, and not scripted in this repo) and **write** the
canonical analytic frames (`extracted_variables_genie_data.csv`,
`extracted_variables_genie_top_genes.txt`, etc.) into `src/exploratory data analysis/`.
However, `_lib.load_cohort()` / `load_top_genes()` **read** those frames from
`data/processed/`. Because those locations differ - and `data/processed/` is not present in
the repo - you must ensure the prepared frames are placed under `data/processed/` before
running any Aim analysis. Never commit raw or identifiable patient data; use the
de-identification script and keep private data outside the tracked tree.

### Environment
There is no pinned environment spec in the repo (no `requirements.txt`, `environment.yml`, or
`renv.lock`) - see README for the manual install list. Adding a pinned spec is a welcome
contribution.

## Running and testing

The pipeline is run stage by stage:
1. Pull / extract (R): `Rscript "src/data collection and processing/pull_genie_data.R"`
   (emits `data/extracted_variables_of_interest.xlsx`).
2. Build the analytic frame (Python): run `harmonize_genie.py`. Pathway and gene-binary
   columns are attached during data prep / in the upstream masters - **not** by
   `add_pathways_genie_bpc.R`, which is an empty 0-byte stub.
3. Enrich (Python): `enrich_harmonized.py` (adds alias, competing-event, and GENIE regimen
   columns; rewrites the cohort CSVs in place).
4. Stage the analytic frames under `data/processed/` so `_lib` can find them.
5. Analyze (Python): run the Aim 1/2/3 scripts; they take a `--cohort genie`
   argument and write CSVs to `src/modeling/genie/aimN/` and figures to
   `manuscript components/genie/aimN/`.

The survival + XGBoost documentation in `src/modeling/README_survival_AFT_pipeline.md`
describes a standalone `survival_and_xgb_analysis.py` CLI (`--xlsx`, `--time-col`,
`--event-col`, `--strata`, `--topn-genes`) that does not match the current `_lib`-based
scripts or the canonical schema; treat it as an older/separate variant and reconcile before
relying on it. There is no automated test suite; validate changes by re-running the relevant
stage and confirming output tables/figures are sane (row counts, event counts, p-values,
concordance indices printed to stdout).

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

## References

The full bibliography lives in [`references/REFERENCES.md`](references/REFERENCES.md) (source
papers, the CBBio driver-gene method and its cited papers, and driver-caller tool sources).
When adding a method or tool, cite it there and link from the relevant prep doc rather than
inlining full references.
