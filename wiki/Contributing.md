# Contributing

Guidelines for contributing to this research codebase. These conventions are inferred
from existing scripts and are intended to keep results reproducible.

---

## Developer setup

See [Getting Started](Getting-Started.md) for full setup instructions.

Quick summary:

```zsh
# Python
python3 -m venv .venv && source .venv/bin/activate
pip install pandas numpy lifelines statsmodels scipy scikit-learn xgboost matplotlib seaborn

# R
Rscript -e 'install.packages(c("dplyr","readr","tidyr","purrr","stringr","data.table","openxlsx","broom","lubridate","xgboost","SHAPforxgboost"))'
```

---

## Coding conventions

### Single source of truth for cohorts

- Add or modify cohorts via `CohortSpec` in `src/modeling/_lib.py`
- Use `load_cohort()` / `load_top_genes()` to read data, not bare `pd.read_csv()`
- Pass cohort key (`genie`, `tcga`, `msk18`) to scripts via `--cohort`

### Respect the canonical schema

- Keep column names consistent with `src/data collection and processing/harmonization_spec.md`
- If a source lacks a variable, create the column as NA rather than dropping it
- Never invent values; record mapping decisions in a side-car `harmonization_decisions.csv`

### Python style

- Python 3.9+ features (`from __future__ import annotations`, modern type hints)
- Snake_case for all identifiers
- Start each script with a module docstring describing inputs, outputs, and the analytic step
- Set `matplotlib.use("Agg")` before importing `pyplot` (scripts must render headless)
- Use informative `print()` progress lines
- Use single hyphens (-) in comments and docs, not em dashes

### R style

- Tidyverse style (dplyr pipelines, readr, purrr)
- Absolute / explicit output paths

### Statistical conventions

- Apply multiple-testing correction (Benjamini-Hochberg) wherever the existing code does
- Use `drop_low_variance` and the `MIN_EVENTS = 25` guard for small-cohort safety
- Include Schoenfeld diagnostics alongside Cox PH fits

---

## Running the pipeline

The pipeline is run stage by stage; see [Pipeline Run Guide](Pipeline-Run-Guide.md).

1. Pull / extract (R)
2. Harmonize (Python, per cohort)
3. Enrich (Python)
4. Stage analytic frames under `data/processed/`
5. Run Aim 1/2/3 and ML analyses (Python)

There is no automated test suite. Validate changes by re-running the relevant stage and
confirming output tables/figures are sane (row counts, event counts, p-values, C-indices).

---

## What to commit

**Commit only:**
- Source code (`.py`, `.R`, `.ipynb`, `.md`)
- Generated result tables and figures intended for the manuscript
  (`src/modeling/<cohort>/aimN/`, `manuscript components/`, `reports/figures/`)

**Never commit:**
- Raw or identifiable patient data
- `data/` directory contents
- `extracted_variables_*` CSV files
- `*_retrieved_data/` folders
- Private / machine-specific data files

The `.gitignore` enforces these rules. Use `create_deidentified_dataset.R` and keep
private data outside the tracked tree.

---

## Adding a new cohort

1. Add a new `CohortSpec` entry in `src/modeling/_lib.py`
2. Write a `harmonize_<cohort>.py` script following the existing pattern and the schema
   in `harmonization_spec.md`
3. Run the harmonize script and copy outputs to `data/processed/`
4. Verify the cohort loads: `_lib.load_cohort("<key>")`
5. Run Aim 1 with `--cohort <key>` and check output tables/figures

---

## Adding a new analysis

1. Import `_lib` for cohort loading, covariate prep, and output-dir management
2. Accept `--cohort {genie,tcga,msk18,all}` as a CLI argument (use `argparse`)
3. Write CSV result tables to `src/modeling/<cohort>/aimN/` via `_lib.ensure_dirs()`
4. Write figures to `manuscript components/<cohort>/aimN/`
5. Add a module docstring describing inputs, outputs, and the analytic step

---

## Reproducibility checklist before submitting changes

- [ ] Script runs to completion without errors on a clean checkout (after data placement)
- [ ] Output CSVs have expected row counts and no unexpected empty columns
- [ ] Figures render correctly (plausible KM curves, forest plots, concordance values)
- [ ] No raw patient data was committed (`git status` shows only code/figure files)
- [ ] Module docstring updated if inputs or outputs changed
