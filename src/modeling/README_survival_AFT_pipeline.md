# Survival and ML modeling scripts (`src/modeling/`)

This directory holds two related-but-distinct families of survival code. Know which one
you are running.

> **History note.** An earlier version of this README described a single standalone
> `survival_and_xgb_analysis.py` CLI driven by a `merged_genie.xlsx` workbook with
> `DFS_MONTHS`/`DFS_STATUS`/`SUBTYPE`/`MANTIS_BIN`/`G__*` columns. That file name and that
> schema are **not** what the repo currently contains. This page has been reconciled to the
> scripts that actually exist; the older column contract survives only inside the
> standalone benchmark scripts noted below and does not match the canonical analytic frame.

---

## A. Canonical Aim 2 / Aim 3 survival (the current pipeline)

These are the scripts the [Pipeline Run Guide](../../wiki/Pipeline-Run-Guide.md) drives.
They share the `_lib.py` backbone, take `--cohort {genie,tcga,msk18}`, read the canonical
`data/processed/extracted_variables_<cohort>_data.csv` frame via `_lib.load_cohort()`, and
use the canonical schema (`OS_months`, `os_status_bin`, `tt_brain_met_mos`,
`brain_met_event`, `top5_any_mutated`, `receptor_primary_cat`, bare-HUGO gene columns).

| Script | Aim | Cohort filter | Endpoint |
|---|---|---|---|
| `proportional_hazard_afttest.py` | Aim 2 - time to brain met | `brain_met_at_dx == 0` | `tt_brain_met_mos` / `brain_met_event` |
| `aim3_os.py` | Aim 3 - overall survival | `any_brain_met == 1` | `OS_months` / `os_status_bin` |

Each does: KM stratified by `top5_any_mutated` + log-rank; Cox PH with `_lib.prep_covariates`
/ `drop_low_variance` covariates and scaled Schoenfeld PH diagnostics; AFT distribution
selection (Weibull / lognormal / log-logistic by AIC) then a full AFT fit. Both enforce
`MIN_EVENTS = 25` before fitting Cox / AFT.

```zsh
python proportional_hazard_afttest.py --cohort genie
python aim3_os.py --cohort genie          # repeat with --cohort tcga / msk18
```

Outputs: CSV tables under `src/modeling/<cohort>/aim{2,3}/`; figures under
`manuscript components/<cohort>/aim{2,3}/`.

---

## B. Standalone XGBoost-AFT ML benchmark

These implement the machine-learning survival benchmark (Stage 8 of the run guide). They
are older, standalone scripts that do **not** yet use `_lib` or the canonical schema - they
still reference the legacy `merged_genie.xlsx` / `DFS_*` / `SUBTYPE` / `MANTIS_BIN` / `G__*`
contract and, in one case, a placeholder input path. Treat them as a reference
implementation to port onto the canonical frame; see "Known gaps" below.

| Script | Role |
|---|---|
| `stratifiedKM_CoxFG_feature_prep_AFT.py` | The combined "survival + ML" driver (the script the old README called `survival_and_xgb_analysis.py`): KM, uni/multivariate Cox PH, Fine-Gray competing risks, and the XGBoost-AFT fit. Reads `merged_genie.xlsx`; writes under `output/survival/`. |
| `xgb_aft_preprocessing_feature_constuction_train_validate_evaluate.py` | Standalone feature construction + XGBoost-AFT train/validate/evaluate (concordance index). Input path is a `path/to/cleaned_data.csv` placeholder that must be set. |
| `xgb_aft_shap_feature_importance.py` | SHAP-like feature contributions for a trained XGB-AFT model (`pred_contribs`); writes `shap_contribs.csv`, `shap_feature_ranking.csv`, `shap_topN_bar.png`. |
| `xgb_aft_feature_.processing_feature_processing_feature_importance.py` | Near-duplicate of the SHAP script above (same docstring/outputs); consolidate the two. |
| `shap analysis and plot generation.R` | R SHAP plots via `SHAPforxgboost` (summary / dependence / force plots). |

XGBoost-AFT label encoding: intervals `(lower, upper)` with `+inf` upper bound for
right-censored rows; XGBoost must be `>= 1.6` for the `survival:aft` objective.

```zsh
python stratifiedKM_CoxFG_feature_prep_AFT.py            # combined survival + XGB-AFT
python xgb_aft_preprocessing_feature_constuction_train_validate_evaluate.py
python xgb_aft_shap_feature_importance.py
Rscript "shap analysis and plot generation.R"
```

Outputs: `reports/figures/ml_benchmark/` (and, for the standalone driver, `output/survival/`).

---

## Install

Dependencies are pinned at the repo root, not here:

```zsh
python3 -m venv .venv && source .venv/bin/activate
pip install -r ../../requirements.txt     # or: conda env create -f ../../environment.yml
```

`lifelines` provides the Cox/KM/AFT fitters and (where available) `FineAndGrayFitter`; if
that class is absent in your lifelines version the Fine-Gray step is skipped with a message.

---

## Known gaps (to reconcile before relying on family B)

1. **Schema mismatch.** The family-B scripts expect `merged_genie.xlsx` and
   `DFS_*`/`SUBTYPE`/`MANTIS_BIN`/`G__*` columns. The canonical frame is
   `extracted_variables_<cohort>_data.csv` with `OS_months`/`os_status_bin`/
   `tt_brain_met_mos`/`brain_met_event`, `receptor_primary_cat`, bare-HUGO gene columns,
   and `pathway_*` indicators. Port family B onto `_lib.load_cohort()` + the canonical
   columns so it runs on current data.
2. **Placeholder input.** `xgb_aft_preprocessing_feature_constuction_train_validate_evaluate.py`
   opens `path/to/cleaned_data.csv` - wire it to `_lib.load_cohort()` (or a real path).
3. **Duplicate SHAP scripts.** `xgb_aft_shap_feature_importance.py` and
   `xgb_aft_feature_.processing_feature_processing_feature_importance.py` are near-identical;
   keep one.
