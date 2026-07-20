# Getting Started

This page covers environment setup, data placement, and the prerequisites needed before
running any part of the pipeline.

---

## Prerequisites

### Python (3.9+)

Python 3.9 or later is required. The shared library (`_lib.py`) uses
`from __future__ import annotations` and modern type hints. XGBoost must be >= 1.6 for
the AFT survival objective.

```zsh
python3 -m venv .venv
source .venv/bin/activate
pip install pandas numpy scipy scikit-learn statsmodels lifelines \
            xgboost matplotlib seaborn openpyxl shap
```

Some scripts and notebooks additionally use:

```zsh
pip install lightgbm optuna scikit-survival
```

### R

Install the packages used across all R scripts:

```r
install.packages(c(
  "dplyr", "readr", "tidyr", "purrr", "stringr",
  "data.table", "openxlsx", "broom", "lubridate",
  "xgboost", "SHAPforxgboost"
))
```

> **Note.** There is no pinned environment spec in the repo (`requirements.txt`,
> `environment.yml`, or `renv.lock`). Adding one is a tracked gap - see
> [Known Issues and Gaps](Known-Issues-and-Gaps.md).

---

## Data setup

### What is git-ignored

All raw and derived data are git-ignored:

- `data/` directory (all contents)
- `*.csv`, `*.xlsx`
- `**/extracted_variables_*`
- `*_retrieved_data/` folders

**Never commit raw or identifiable patient data.**

### Where the raw data comes from

The pipeline is built on three de-identified public cancer-genomics datasets. You must
download them separately and place them at the paths the scripts expect.

| Cohort | Source | Required files |
|--------|--------|---------------|
| GENIE BPC BRCA | AACR Project GENIE BPC | MAF + clinical tables; referenced from `data/processed/genie_bpc_v1_*.csv` |
| TCGA BRCA | GDC / cBioPortal | cBioPortal study export |
| Breast MSK 2018 | cBioPortal | cBioPortal study export |

> The R pull scripts reference **machine-specific absolute paths** under
> `/Users/robertjames/loc/data private/...`. Adjust these for your machine before running.

### Required directory structure before running analyses

```
data/
  processed/
    extracted_variables_genie_data.csv
    extracted_variables_genie_top_genes.txt
    extracted_variables_tcga_data.csv
    extracted_variables_tcga_top_genes.txt
    extracted_variables_breast_msk_2018_data.csv
    extracted_variables_breast_msk_2018_top_genes.txt
    genie_bpc_v1_sample_master_full.csv
    genie_bpc_v1_mutations.csv
    genie_bpc_v1_regimens.csv
  hugo_symbols.xlsx
```

The `extracted_variables_*` CSVs are produced by the harmonize scripts (step 2 of the
pipeline) and written to `src/exploratory data analysis/`. You must copy or symlink
them into `data/processed/` so that `_lib.load_cohort()` can find them.
See [Pipeline Run Guide](Pipeline-Run-Guide.md) for details.

---

## Notebooks

The notebooks in `notebooks/` hard-code a `PROJ` root pointing to the author's iCloud
`Documents` folder. Update the `PROJ` variable at the top of each notebook to point to
your local clone before running.

---

## Verifying the setup

After installing dependencies and placing data, test that the shared library loads:

```python
import sys
sys.path.insert(0, "src/modeling")
import _lib
df = _lib.load_cohort("genie")
print(df.shape, df.columns.tolist()[:5])
```

If `data/processed/extracted_variables_genie_data.csv` exists, this should print the
frame shape and the first five column names without error.
