# Survival and ML Analysis Pipeline

This pipeline performs survival analysis and an AFT-based machine learning model on the merged GENIE dataset.

## What it does

- Kaplan–Meier curves overall and stratified by SUBTYPE and MANTIS_BIN
- Univariate Cox PH for gene indicators (`G__*`) and baseline covariates
- Multivariate Cox PH with PH diagnostics (Schoenfeld tests)
- Fine–Gray competing-risk model (if a competing event column exists)
- XGBoost Accelerated Failure Time (AFT) model with feature importance

## Inputs

- `datasets_analysis_dictionary/merged_genie.xlsx` (Excel; contains time-to-event and event columns)
- Optional: `output/fine_gray_python/per_gene_cs_cox_sksurv_bootstrap.csv` (provides top-N genes)

## Outputs

Written under `output/survival/`:
- `km/` — PNG plots of Kaplan–Meier curves
- `cox/` — `cox_univariate.csv`, `cox_multivariate_summary.csv`, and stratified results
- `finegray/` — `finegray_summary.csv` (if applicable)
- `xgb/` — `xgb_aft_model.json`, `xgb_feature_importance.csv`, `xgb_aft_metrics.json`
- `logs/` — Cox PH diagnostics text files
- `run_summary.json` — overview of detected columns and steps performed

## Install dependencies (macOS, zsh)

```zsh
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run

```zsh
# Basic run (auto-detect time/event columns)
python survival_and_xgb_analysis.py \
  --xlsx datasets_analysis_dictionary/merged_genie.xlsx

# Explicit columns (example)
python survival_and_xgb_analysis.py \
  --xlsx datasets_analysis_dictionary/merged_genie.xlsx \
  --time-col DFS_MONTHS \
  --event-col DFS_STATUS \
  --competing-col OS_STATUS \
  --strata SUBTYPE,MANTIS_BIN \
  --topn-genes 10
```

## Column expectations

- Time-to-event: tries `DFS_MONTHS`, `DFS_TIME`, `TIME_TO_DFS`, `TIME_TO_EVENT`, `MONTHS_TO_DFS`, `PFS_MONTHS`, `OS_MONTHS` (case-insensitive)
- Event indicator (1=event, 0=censored): tries `DFS_STATUS`, `DFS_EVENT`, `EVENT`, `PFS_EVENT`, `OS_EVENT`, `OS_STATUS`
- Competing event for Fine–Gray (optional): tries `DEATH_EVENT`, `OS_STATUS`, `DECEASED`, `COMPETING_EVENT`
- Gene indicators: columns that start with `G__`; if missing, script will derive them from a HUGO-like column using the top genes CSV

## Notes

- Fine–Gray uses `lifelines.FineAndGrayFitter`; if lifelines is not installed with this class, the step will be skipped.
- XGBoost AFT requires xgboost>=1.6; labels are encoded as (lower, upper) bounds with +inf for right-censored observations.
- For Cox PH, we include baseline covariates and select up to 10 genes with p<0.1 from the univariate screen.

## SHAP explainability (R)

An R helper script is included to reproduce SHAP plots similar to your example.

- Script: `shap_feature_explainer.R`
- It installs SHAPforxgboost (CRAN first; falls back to GitHub if needed), trains a small demo model if no CSV is provided, and saves plots.

Run (demo data):

```zsh
Rscript shap_feature_explainer.R --outdir output/shap
```

Run (your CSV):

```zsh
Rscript shap_feature_explainer.R \
  --csv path/to/your/features.csv \
  --outcome diffcwv \
  --outdir output/shap
```

Outputs:
- `shap_ranked_features.csv`
- `shap_summary.png`, `shap_summary_light.png`
- `shap_dependence_<topfeat>.png` (+ colored variant when available)
- `shap_interaction_<feat1>_<feat2>.png` (if computed)
- `shap_force_plot.png`, `shap_force_plot_bygroup.png`

## SHAP-like summary for XGB AFT (Python)

For the survival AFT model trained by `survival_and_xgb_analysis.py`, use `xgb_aft_shap_summary.py` to compute contribution matrices and a top-N bar plot.

Run:

```zsh
python xgb_aft_shap_summary.py \
  --xlsx datasets_analysis_dictionary/merged_genie.xlsx \
  --model-json output/survival/xgb/xgb_aft_model.json \
  --per-gene-csv output/fine_gray_python/per_gene_cs_cox_sksurv_bootstrap.csv \
  --topn-genes 10 \
  --outdir output/survival/xgb
```

Outputs:
- `output/survival/xgb/shap_contribs.csv` (N x (F+1), last is bias)
- `output/survival/xgb/shap_feature_ranking.csv`
- `output/survival/xgb/shap_topN_bar.png`

## Troubleshooting installs

If dependency installation failed earlier, ensure you install from the directory that contains `requirements.txt`:

```zsh
cd /Users/robertjames/Downloads/Transfer
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
