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

## B. XGBoost-AFT ML benchmark

The machine-learning survival benchmark (Stage 8 of the run guide). The XGB-AFT trainer
and its SHAP explainer have been **ported onto the canonical `_lib` schema** (`--cohort` /
`--aim`, `load_cohort()`, canonical endpoints and covariates), so they share the feature
definition and run on the same `extracted_variables_<cohort>_data.csv` frame as Aim 2/3.
The combined KM/Cox/Fine-Gray driver is still on the legacy contract (see "Known gaps").

| Script | Role |
|---|---|
| `xgb_aft_preprocessing_feature_constuction_train_validate_evaluate.py` | XGBoost-AFT train/validate/evaluate on the canonical frame. `--cohort {genie,tcga,msk18,all}` `--aim {aim3,aim2}`. Features = `prep_covariates` clinical block + top-10 bare-HUGO genes + `pathway_*` + genomic burden. Exposes `assemble_xy` / `ENDPOINTS` for reuse. |
| `xgb_aft_shap_feature_importance.py` | SHAP-like contributions (`pred_contribs`) for the trained model; imports `assemble_xy` from the trainer so the explanation matrix matches training exactly. |
| `stratifiedKM_CoxFG_feature_prep_AFT.py` | Legacy combined "survival + ML" driver (the script the old README called `survival_and_xgb_analysis.py`): KM, uni/multivariate Cox PH, Fine-Gray, XGB-AFT. Still reads `merged_genie.xlsx` and writes under `output/survival/`. |
| `shap analysis and plot generation.R` | R SHAP plots via `SHAPforxgboost` (summary / dependence / force plots). |

XGBoost-AFT label encoding: intervals `(lower, upper)` with `+inf` upper bound for
right-censored rows; XGBoost must be `>= 1.6` for the `survival:aft` objective. The trainer
enforces the same `MIN_EVENTS = 25` guard as Aim 2/3.

```zsh
python xgb_aft_preprocessing_feature_constuction_train_validate_evaluate.py --cohort genie --aim aim3
python xgb_aft_shap_feature_importance.py --cohort genie --aim aim3
Rscript "shap analysis and plot generation.R"
```

Outputs: model JSON / metrics / feature-importance CSV under
`src/modeling/<cohort>/ml_benchmark/`; figures under `reports/figures/ml_benchmark/`.

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

## Known gaps

Resolved:
- **XGB-AFT trainer + SHAP ported** to `_lib.load_cohort()` and the canonical schema
  (`--cohort`/`--aim`); the `path/to/cleaned_data.csv` placeholder is gone, and the two
  scripts now share one feature definition.
- **Duplicate SHAP script removed** (`xgb_aft_feature_.processing_feature_processing_feature_importance.py`).

Remaining:
1. **Legacy driver not yet ported.** `stratifiedKM_CoxFG_feature_prep_AFT.py` still expects
   `merged_genie.xlsx` and `DFS_*`/`SUBTYPE`/`MANTIS_BIN`/`G__*` columns. Its KM/Cox/Fine-Gray
   pieces overlap the canonical Aim 2/3 scripts (family A); either port it onto
   `_lib.load_cohort()` or retire it in favor of family A + the ported XGB-AFT benchmark.
2. **Not yet executed on real data.** The ported scripts byte-compile but have not been run
   against a cohort frame here; validate C-index / feature ranking on first real run.
