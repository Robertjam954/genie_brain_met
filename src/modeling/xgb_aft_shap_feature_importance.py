"""SHAP-like feature contributions for the trained XGBoost-AFT survival model.

Ported onto the canonical `_lib` schema. Reuses the exact feature definition from
`xgb_aft_preprocessing_feature_constuction_train_validate_evaluate.py` (via
`assemble_xy`) so the explanation matrix matches training, loads the saved model,
and computes XGBoost `pred_contribs`.

Run from `src/modeling/` after the trainer, e.g.:
    python xgb_aft_shap_feature_importance.py --cohort genie --aim aim3

Inputs  (per cohort/aim, produced by the trainer):
    src/modeling/<cohort>/ml_benchmark/<label>_xgb_aft_model.json
    src/modeling/<cohort>/ml_benchmark/<label>_xgb_features.txt
Outputs (same dir + reports/figures/ml_benchmark/):
    <label>_shap_contribs.csv         (N x (F+1); last column is the bias term)
    <label>_shap_feature_ranking.csv  (feature, mean_abs_contrib)
    <cohort>_<label>_shap_bar.png     (top-N bar chart)
"""
from __future__ import annotations

import argparse
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb

from _lib import COHORTS
from xgb_aft_preprocessing_feature_constuction_train_validate_evaluate import (
    ENDPOINTS, FIG_DIR, assemble_xy, model_dir,
)


def run(cohort: str, aim: str, topn_plot: int = 20) -> None:
    spec = ENDPOINTS[aim]
    out_dir = model_dir(cohort)
    model_path = out_dir / f"{spec.label}_xgb_aft_model.json"
    feat_path = out_dir / f"{spec.label}_xgb_features.txt"
    if not model_path.exists() or not feat_path.exists():
        print(f"   SKIP {cohort}/{spec.label}: model or feature list missing "
              f"(run the trainer first)")
        return

    feat_names = [ln.strip() for ln in feat_path.read_text().splitlines() if ln.strip()]

    # Rebuild features with the identical definition, then align to the trained order.
    X, _time, _event = assemble_xy(cohort, aim)
    X = X.reindex(columns=feat_names)  # add any missing trained cols as NaN, drop extras

    bst = xgb.Booster()
    bst.load_model(str(model_path))

    dmat = xgb.DMatrix(X.to_numpy(dtype=float), feature_names=feat_names)
    contrib = bst.predict(dmat, pred_contribs=True)  # (n, F+1); last col = bias

    pd.DataFrame(contrib, columns=feat_names + ["bias"]).to_csv(
        out_dir / f"{spec.label}_shap_contribs.csv", index=False)

    abs_mean = np.abs(contrib[:, :-1]).mean(axis=0)
    ranking = (pd.DataFrame({"feature": feat_names, "mean_abs_contrib": abs_mean})
               .sort_values("mean_abs_contrib", ascending=False))
    ranking.to_csv(out_dir / f"{spec.label}_shap_feature_ranking.csv", index=False)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    top = ranking.head(topn_plot)
    plt.figure(figsize=(10, max(4, 0.4 * len(top))))
    plt.barh(top["feature"][::-1], top["mean_abs_contrib"][::-1], color="#1f77b4")
    plt.xlabel("Mean |contribution|")
    plt.title(f"XGB-AFT SHAP-like importance - {cohort}/{spec.label}")
    plt.tight_layout()
    plt.savefig(FIG_DIR / f"{cohort}_{spec.label}_shap_bar.png", dpi=150)
    plt.close()
    print(f"   [{cohort}/{spec.label}] wrote SHAP outputs to {out_dir}")


def main() -> None:
    p = argparse.ArgumentParser(description="XGB-AFT SHAP-like feature importance")
    p.add_argument("--cohort", choices=list(COHORTS.keys()) + ["all"], default="all")
    p.add_argument("--aim", choices=list(ENDPOINTS.keys()), default="aim3")
    p.add_argument("--topn-plot", type=int, default=20)
    args = p.parse_args()
    cohorts = list(COHORTS.keys()) if args.cohort == "all" else [args.cohort]
    for c in cohorts:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run(c, args.aim, topn_plot=args.topn_plot)


if __name__ == "__main__":
    main()
