"""XGBoost-AFT survival benchmark on the canonical analytic frame.

Ported onto the `_lib` cohort schema (was a standalone script that opened a
`path/to/cleaned_data.csv` placeholder with a legacy `G__*`/`SUBTYPE`/`MANTIS_BIN`
column contract). It now reads `data/processed/extracted_variables_<cohort>_data.csv`
via `_lib.load_cohort()` and uses the same canonical endpoints and covariates as the
Aim 2 / Aim 3 survival scripts.

Endpoints (choose with --aim):
  aim3 (default): overall survival in the brain-mets-ever cohort
                  (any_brain_met == 1; OS_months / os_status_bin)
  aim2          : time to brain met in the no-CNS-at-dx cohort
                  (brain_met_at_dx == 0; tt_brain_met_mos / brain_met_event)

Features: the `_lib.prep_covariates` clinical block + top-N bare-HUGO gene indicators
+ Sanchez-Vega pathway indicators + genomic burden (mutation_count_all_sites_sum,
t_alt_count_max), reduced with `_lib.drop_low_variance`.

The `assemble_xy` and `EndpointSpec` helpers are imported by
`xgb_aft_shap_feature_importance.py` so training and explanation share one feature
definition. Run from `src/modeling/` (like the other modeling scripts).

Outputs (per cohort/aim):
  src/modeling/<cohort>/ml_benchmark/<label>_xgb_aft_model.json
  src/modeling/<cohort>/ml_benchmark/<label>_xgb_features.txt
  src/modeling/<cohort>/ml_benchmark/<label>_xgb_metrics.json
  src/modeling/<cohort>/ml_benchmark/<label>_xgb_feature_importance.csv
  reports/figures/ml_benchmark/<label>_xgb_feature_importance.png
"""
from __future__ import annotations

import argparse
import json
import warnings
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split

from _lib import (
    COHORTS, PATHWAY_COLS, PROJECT_ROOT, drop_low_variance, load_cohort,
    load_top_genes, prep_covariates,
)

MIN_EVENTS = 25
CONTINUOUS_GENOMIC = ["mutation_count_all_sites_sum", "t_alt_count_max"]
FIG_DIR = PROJECT_ROOT / "reports" / "figures" / "ml_benchmark"


@dataclass
class EndpointSpec:
    label: str
    cohort_filter: str      # column that must equal filter_value
    filter_value: int
    time_col: str
    event_col: str
    include_stage_iv: bool


ENDPOINTS: dict[str, EndpointSpec] = {
    "aim3": EndpointSpec("aim3_os", "any_brain_met", 1, "OS_months",
                         "os_status_bin", include_stage_iv=False),
    "aim2": EndpointSpec("aim2", "brain_met_at_dx", 0, "tt_brain_met_mos",
                         "brain_met_event", include_stage_iv=True),
}


def model_dir(cohort: str) -> Path:
    d = PROJECT_ROOT / "src" / "modeling" / cohort / "ml_benchmark"
    d.mkdir(parents=True, exist_ok=True)
    return d


def assemble_xy(cohort: str, aim: str) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """Build the (X features, time, event) triple for a cohort/endpoint.

    Row-aligned and NaN-preserving (XGBoost treats NaN as missing). Shared by the
    trainer and the SHAP script so both use an identical feature definition.
    """
    spec = ENDPOINTS[aim]
    df = load_cohort(cohort)
    base = df[df[spec.cohort_filter] == spec.filter_value].copy()

    # Clinical covariate block (one-hot encoded, numeric-coerced)
    cov = prep_covariates(base, include_stage_iv=spec.include_stage_iv)

    # Top-N bare-HUGO gene indicators present in the frame
    genes = load_top_genes(cohort, n=10)
    gene_cols = [g for g in genes if g in base.columns]
    gene_block = base[gene_cols].apply(pd.to_numeric, errors="coerce")

    # Sanchez-Vega pathway indicators present in the frame
    pw_cols = [c for c in PATHWAY_COLS if c in base.columns]
    pw_block = base[pw_cols].apply(pd.to_numeric, errors="coerce")

    # Continuous genomic burden
    cont_cols = [c for c in CONTINUOUS_GENOMIC if c in base.columns]
    cont_block = base[cont_cols].apply(pd.to_numeric, errors="coerce")

    X = pd.concat(
        [b.reset_index(drop=True) for b in (cov, gene_block, pw_block, cont_block)],
        axis=1,
    )
    X = drop_low_variance(X)

    time = pd.to_numeric(base[spec.time_col], errors="coerce").reset_index(drop=True)
    event = pd.to_numeric(base[spec.event_col], errors="coerce").reset_index(drop=True)

    # Keep rows with a valid, positive time and a known event
    ok = time.notna() & event.notna() & (time > 0)
    return X.loc[ok].reset_index(drop=True), time[ok].reset_index(drop=True), \
        event[ok].astype(int).reset_index(drop=True)


def _aft_bounds(time: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    y_lower = time.astype(float).copy()
    y_upper = time.astype(float).copy()
    y_upper[event == 0] = np.inf   # right-censored -> upper bound +inf
    return y_lower, y_upper


def run(cohort: str, aim: str) -> None:
    spec = ENDPOINTS[aim]
    X, time, event = assemble_xy(cohort, aim)
    n_events = int(event.sum())
    print(f"[{cohort}/{spec.label}] n={len(X)}  events={n_events}  features={X.shape[1]}")

    out_dir = model_dir(cohort)
    if n_events < MIN_EVENTS:
        msg = f"Skipped XGB-AFT - {n_events} events (< {MIN_EVENTS})."
        (out_dir / f"{spec.label}_xgb_skipped.txt").write_text(msg + "\n")
        print(f"   SKIP: {msg}")
        return

    feat_names = list(X.columns)
    y_lower, y_upper = _aft_bounds(time.to_numpy(), event.to_numpy())

    X_tr, X_va, lo_tr, lo_va, up_tr, up_va, ev_tr, ev_va = train_test_split(
        X.to_numpy(dtype=float), y_lower, y_upper, event.to_numpy(),
        test_size=0.2, random_state=42,
        stratify=event if event.nunique() > 1 else None,
    )
    print(f"   train={len(X_tr)}  val={len(X_va)}")

    dtrain = xgb.DMatrix(X_tr, feature_names=feat_names)
    dtrain.set_float_info("label_lower_bound", lo_tr)
    dtrain.set_float_info("label_upper_bound", up_tr)
    dval = xgb.DMatrix(X_va, feature_names=feat_names)
    dval.set_float_info("label_lower_bound", lo_va)
    dval.set_float_info("label_upper_bound", up_va)

    params = {
        "objective": "survival:aft",
        "eval_metric": "aft-nloglik",
        "aft_loss_distribution": "normal",
        "aft_loss_distribution_scale": 1.0,
        "tree_method": "hist",
        "learning_rate": 0.05,
        "max_depth": 6,
        "seed": 42,
    }

    # Pick boosting rounds by CV, with a fixed fallback if CV cannot run.
    try:
        cv = xgb.cv(params, dtrain, num_boost_round=1000, nfold=5,
                    early_stopping_rounds=10, metrics="aft-nloglik", seed=42)
        best_rounds = len(cv)
    except Exception as e:
        best_rounds = 200
        print(f"   WARN: xgb.cv failed ({e}); using {best_rounds} rounds")
    print(f"   best_rounds={best_rounds}")

    bst = xgb.train(params, dtrain, num_boost_round=best_rounds,
                    evals=[(dval, "validation")], early_stopping_rounds=10,
                    verbose_eval=False)

    pred_va = bst.predict(dval)
    # Higher predicted log-time = longer survival = lower risk; negate for concordance.
    c_index = concordance_index(lo_va, -pred_va, ev_va)
    predicted_time = np.exp(pred_va)
    rmse = float(np.sqrt(np.mean(
        (lo_va[ev_va == 1] - predicted_time[ev_va == 1]) ** 2))) \
        if (ev_va == 1).any() else float("nan")
    print(f"   val C-index={c_index:.3f}  RMSE(uncensored)={rmse:.3f}")

    # ----- Artifacts -----
    bst.save_model(str(out_dir / f"{spec.label}_xgb_aft_model.json"))
    (out_dir / f"{spec.label}_xgb_features.txt").write_text(
        "\n".join(feat_names) + "\n")
    json.dump(
        {"cohort": cohort, "aim": spec.label, "n": int(len(X)), "events": n_events,
         "n_features": X.shape[1], "best_rounds": int(best_rounds),
         "val_c_index": float(c_index), "val_rmse_uncensored": rmse},
        open(out_dir / f"{spec.label}_xgb_metrics.json", "w"), indent=2)

    gain = bst.get_score(importance_type="gain")
    imp = (pd.DataFrame({"feature": list(gain.keys()),
                         "gain": list(gain.values())})
           .sort_values("gain", ascending=False))
    imp.to_csv(out_dir / f"{spec.label}_xgb_feature_importance.csv", index=False)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    top = imp.head(20)
    plt.figure(figsize=(10, max(4, 0.4 * len(top))))
    plt.barh(top["feature"][::-1], top["gain"][::-1], color="#1f77b4")
    plt.xlabel("Gain"); plt.title(f"XGB-AFT feature importance - {cohort}/{spec.label}")
    plt.tight_layout()
    plt.savefig(FIG_DIR / f"{cohort}_{spec.label}_xgb_feature_importance.png", dpi=150)
    plt.close()
    print(f"   wrote artifacts to {out_dir} and figure to {FIG_DIR}")


def main() -> None:
    p = argparse.ArgumentParser(description="XGBoost-AFT survival benchmark")
    p.add_argument("--cohort", choices=list(COHORTS.keys()) + ["all"], default="all")
    p.add_argument("--aim", choices=list(ENDPOINTS.keys()), default="aim3")
    args = p.parse_args()
    cohorts = list(COHORTS.keys()) if args.cohort == "all" else [args.cohort]
    for c in cohorts:
        print(f"\n=== XGB-AFT benchmark cohort={c} aim={args.aim} ===")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run(c, args.aim)


if __name__ == "__main__":
    main()
