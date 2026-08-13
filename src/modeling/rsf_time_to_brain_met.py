"""Steps 6-7: predictive modeling of time to brain metastasis, and interpretation.

Step 6 - Random Survival Forest
  - restrict to the incident cohort (brain_met_at_dx == 0)
  - one-hot encode the categorical covariates
  - restrict gene features to the genes selected by the step-4 risk models
    (`select_risk_model_genes.py`), falling back to all G_top10_* columns
  - fit a RandomSurvivalForest (n_estimators=200, min_samples_leaf=5,
    max_features='sqrt') and tune it with 5-fold cross-validation
  - evaluate on a held-out test split: Harrell + IPCW concordance index,
    time-dependent (cumulative/dynamic) AUC, and the integrated Brier score

Step 7 - Interpretation
  - permutation variable importance on the test split
  - partial dependence of predicted risk on the top predictors

Death is treated as censoring here: a Random Survival Forest models a single
right-censored endpoint, so the competing-risk interpretation stays with the
Fine-Gray model from step 4.

Usage:
    python3 "src/modeling/rsf_time_to_brain_met.py" --cohort genie
    python3 "src/modeling/rsf_time_to_brain_met.py" --data path/to/frame.csv \
        --genes-file path/to/selected_genes.txt --outdir OUT --figdir FIG
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.inspection import permutation_importance
from sklearn.model_selection import GridSearchCV, KFold, train_test_split
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import (
    concordance_index_censored,
    concordance_index_ipcw,
    cumulative_dynamic_auc,
    integrated_brier_score,
)
from sksurv.util import Surv

from _lib import COHORTS, MODEL_OUT_DIR, FIG_OUT_DIR, load_cohort

TIME_COL = "tt_brain_met_mos"
EVENT_COL = "brain_met_event"
CATEGORICAL = [
    "receptor_primary_cat", "stage_dx_cat", "grade_ord", "smoking_status",
    "insurance", "race_clean", "ethnicity_clean", "sex", "SEQ_ASSAY_ID",
]
NUMERIC = ["age_dx_num", "mutation_count", "tmb"]
GENE_PREFIX = "G_top10_"

MIN_EVENTS = 25
RSF_BASE = dict(n_estimators=200, min_samples_leaf=5, max_features="sqrt",
                n_jobs=-1, random_state=42)
PARAM_GRID = {
    "min_samples_leaf": [3, 5, 10],
    "min_samples_split": [6, 10],
    "max_features": ["sqrt", 0.5],
}


# ------------------------------------------------------------------ data -----

def load_frame(cohort: str, data_path: str | None) -> pd.DataFrame:
    if data_path:
        p = Path(data_path)
        if p.suffix.lower() in {".xlsx", ".xls"}:
            return pd.read_excel(p)
        return pd.read_csv(p, low_memory=False)
    return load_cohort(cohort)


def read_gene_list(path: Path | None, df: pd.DataFrame) -> list[str]:
    """Genes selected in step 5; fall back to every gene-indicator column."""
    if path is not None and path.exists():
        genes = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
        present = [g for g in genes if g in df.columns]
        missing = sorted(set(genes) - set(present))
        if missing:
            print(f"   note: selected genes absent from the frame: {', '.join(missing)}")
        if present:
            print(f"   using {len(present)} risk-model-selected genes")
            return present
        print("   note: none of the selected genes are present - falling back to all genes")
    else:
        print("   note: no selected-genes file - falling back to all gene columns")
    return [c for c in df.columns if c.startswith(GENE_PREFIX)]


def build_design(df: pd.DataFrame, genes: list[str]) -> tuple[pd.DataFrame, np.ndarray]:
    """One-hot encode categoricals, assemble the feature matrix, and build the
    structured survival array. Rows with a missing time/event are dropped."""
    base = df[(df.get("brain_met_at_dx", 0).fillna(0) == 0)].copy()
    time = pd.to_numeric(base[TIME_COL], errors="coerce")
    event = pd.to_numeric(base[EVENT_COL], errors="coerce")
    ok = time.notna() & (time > 0) & event.notna()
    base, time, event = base[ok], time[ok], event[ok]

    num_cols = [c for c in NUMERIC if c in base.columns]
    cat_cols = [c for c in CATEGORICAL if c in base.columns]
    gene_cols = [g for g in genes if g in base.columns]

    X = base[num_cols].apply(pd.to_numeric, errors="coerce")
    for c in num_cols:
        X[c] = X[c].fillna(X[c].median())
    if cat_cols:
        cats = base[cat_cols].astype("string").fillna("Unknown")
        X = pd.concat([X, pd.get_dummies(cats, columns=cat_cols, drop_first=True)], axis=1)
    if gene_cols:
        g = base[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        X = pd.concat([X, g], axis=1)

    X = X.astype(float)
    # Drop constant and near-empty one-hot columns to keep the forest stable.
    keep = [c for c in X.columns if X[c].nunique(dropna=False) > 1]
    dropped = sorted(set(X.columns) - set(keep))
    if dropped:
        print(f"   dropped {len(dropped)} constant feature(s)")
    X = X[keep]

    y = Surv.from_arrays(event=event.astype(bool).to_numpy(), time=time.to_numpy())
    return X, y


def eval_times(y_train: np.ndarray, y_test: np.ndarray,
               n: int = 12) -> np.ndarray:
    """Grid of evaluation times valid for both IPCW censoring and Brier scores:
    strictly inside the follow-up support of train and test."""
    t_train = y_train["time"]
    t_test = y_test["time"]
    lo = max(np.percentile(t_test[y_test["event"]], 10), t_train.min() + 1e-6)
    hi = min(np.percentile(t_test[y_test["event"]], 90), t_train.max(), t_test.max())
    hi = min(hi, np.max(t_test) * 0.999)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.array([])
    return np.linspace(lo, hi, n)


def survival_matrix(model: RandomSurvivalForest, X: pd.DataFrame,
                    times: np.ndarray) -> np.ndarray:
    """Predicted survival probabilities at `times` (n_samples x n_times)."""
    fns = model.predict_survival_function(X)
    return np.row_stack([[fn(t) for t in times] for fn in fns])


# ------------------------------------------------------------- interpretation -

def plot_importance(imp: pd.DataFrame, path: Path, cohort: str, top: int = 20) -> None:
    d = imp.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, max(4, 0.32 * len(d))))
    ax.barh(d["feature"], d["importance_mean"], xerr=d["importance_std"],
            color="#4477aa")
    ax.set_xlabel("Permutation importance (drop in C-index)")
    ax.set_title(f"Variable importance - time to brain metastasis ({cohort})")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def partial_dependence(model: RandomSurvivalForest, X: pd.DataFrame, feature: str,
                       grid_points: int = 12) -> pd.DataFrame:
    """Marginal mean predicted risk score as `feature` is swept over its grid,
    holding the observed joint distribution of the other features fixed."""
    col = X[feature]
    uniq = np.unique(col.dropna().to_numpy())
    if len(uniq) <= grid_points:
        grid = uniq
    else:
        grid = np.quantile(col.dropna(), np.linspace(0.05, 0.95, grid_points))
        grid = np.unique(grid)
    rows = []
    for v in grid:
        Xg = X.copy()
        Xg[feature] = v
        rows.append({"feature": feature, "value": float(v),
                     "mean_risk": float(np.mean(model.predict(Xg)))})
    return pd.DataFrame(rows)


def plot_partial_dependence(pdps: list[pd.DataFrame], path: Path, cohort: str) -> None:
    n = len(pdps)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow), squeeze=False)
    for ax, pdp in zip(axes.ravel(), pdps):
        name = pdp["feature"].iloc[0]
        if len(pdp) <= 2:
            ax.bar([str(v) for v in pdp["value"]], pdp["mean_risk"], color="#4477aa",
                   width=0.5)
        else:
            ax.plot(pdp["value"], pdp["mean_risk"], marker="o", color="#4477aa")
        ax.set_title(name, fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("mean predicted risk")
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle(f"Partial dependence of RSF risk on top predictors - {cohort}")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------- main -----

def run(cohort: str, data_path: str | None, genes_file: str | None,
        outdir: Path, figdir: Path, test_size: float, seed: int,
        n_estimators: int, tune: bool) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    print(f"[step 6] cohort={cohort}")
    df = load_frame(cohort, data_path)
    print(f"   frame: {len(df)} rows x {df.shape[1]} cols")

    gpath = Path(genes_file) if genes_file else (
        MODEL_OUT_DIR / cohort / "risk_models" / "selected_genes.txt")
    genes = read_gene_list(gpath, df)

    X, y = build_design(df, genes)
    n_events = int(y["event"].sum())
    print(f"   design: {X.shape[0]} rows x {X.shape[1]} features, {n_events} events")
    if n_events < MIN_EVENTS:
        msg = f"SKIP RSF: {n_events} brain-met events (< {MIN_EVENTS})"
        print(f"   {msg}")
        (outdir / "rsf_skipped.txt").write_text(msg + "\n")
        return

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=seed, stratify=y["event"])
    print(f"   split: train={len(X_train)} (events={int(y_train['event'].sum())})  "
          f"test={len(X_test)} (events={int(y_test['event'].sum())})")

    base = dict(RSF_BASE, n_estimators=n_estimators, random_state=seed)
    if tune:
        print("   5-fold CV hyperparameter search ...")
        gs = GridSearchCV(
            RandomSurvivalForest(**base), PARAM_GRID,
            cv=KFold(n_splits=5, shuffle=True, random_state=seed),
            n_jobs=1, refit=True, verbose=0,
        )
        gs.fit(X_train, y_train)
        model = gs.best_estimator_
        best_params = gs.best_params_
        cv_score = float(gs.best_score_)
        cv_table = pd.DataFrame(gs.cv_results_)[
            ["params", "mean_test_score", "std_test_score", "rank_test_score"]
        ].sort_values("rank_test_score")
        cv_table.to_csv(outdir / "rsf_cv_results.csv", index=False)
        print(f"   best params: {best_params}  (CV C-index {cv_score:.4f})")
    else:
        model = RandomSurvivalForest(**base).fit(X_train, y_train)
        best_params = {k: v for k, v in base.items() if k in PARAM_GRID}
        cv_score = float("nan")

    # ---- evaluation ----
    risk_test = model.predict(X_test)
    c_harrell = concordance_index_censored(
        y_test["event"], y_test["time"], risk_test)[0]
    metrics: dict[str, object] = {
        "cohort": cohort,
        "n_train": int(len(X_train)),
        "n_test": int(len(X_test)),
        "n_features": int(X.shape[1]),
        "n_events_total": n_events,
        "best_params": best_params,
        "cv_c_index_mean": cv_score,
        "test_c_index_harrell": float(c_harrell),
    }

    times = eval_times(y_train, y_test)
    if times.size:
        try:
            metrics["test_c_index_ipcw"] = float(
                concordance_index_ipcw(y_train, y_test, risk_test, tau=times[-1])[0])
        except Exception as exc:  # pragma: no cover - depends on censoring pattern
            print(f"   IPCW C-index unavailable: {exc}")
        try:
            auc_t, auc_mean = cumulative_dynamic_auc(y_train, y_test, risk_test, times)
            pd.DataFrame({"time": times, "auc": auc_t}).to_csv(
                outdir / "rsf_time_dependent_auc.csv", index=False)
            metrics["test_mean_time_dependent_auc"] = float(auc_mean)
        except Exception as exc:  # pragma: no cover
            print(f"   time-dependent AUC unavailable: {exc}")
            auc_t = None
        try:
            surv = survival_matrix(model, X_test, times)
            metrics["test_integrated_brier_score"] = float(
                integrated_brier_score(y_train, y_test, surv, times))
            pd.DataFrame({"time": times,
                          "mean_predicted_survival": surv.mean(axis=0)}).to_csv(
                outdir / "rsf_predicted_survival.csv", index=False)
        except Exception as exc:  # pragma: no cover
            print(f"   Brier score unavailable: {exc}")
        if auc_t is not None:
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.plot(times, auc_t, marker="o", color="#cc6677")
            ax.axhline(0.5, ls="--", c="grey", lw=1)
            ax.set_xlabel("Months from diagnosis")
            ax.set_ylabel("Time-dependent AUC")
            ax.set_title(f"RSF time-dependent AUC - {cohort} "
                         f"(mean {metrics.get('test_mean_time_dependent_auc', float('nan')):.3f})")
            fig.tight_layout()
            fig.savefig(figdir / "rsf_time_dependent_auc.png", dpi=150)
            plt.close(fig)
    else:
        print("   note: no valid evaluation time grid - skipping AUC / Brier")

    print(f"   test C-index (Harrell) = {c_harrell:.4f}")
    for k in ("test_c_index_ipcw", "test_mean_time_dependent_auc",
              "test_integrated_brier_score"):
        if k in metrics:
            print(f"   {k} = {metrics[k]:.4f}")

    # ---- step 7: interpretation ----
    print("[step 7] permutation importance + partial dependence")
    pi = permutation_importance(model, X_test, y_test, n_repeats=10,
                               random_state=seed, n_jobs=1)
    imp = pd.DataFrame({
        "feature": X_test.columns,
        "importance_mean": pi.importances_mean,
        "importance_std": pi.importances_std,
    }).sort_values("importance_mean", ascending=False)
    imp.to_csv(outdir / "rsf_variable_importance.csv", index=False)
    plot_importance(imp, figdir / "rsf_variable_importance.png", cohort)
    print("   top predictors: " + ", ".join(imp["feature"].head(5)))

    top_feats = [f for f in imp["feature"].head(6)]
    pdps = [partial_dependence(model, X_test, f) for f in top_feats]
    pd.concat(pdps, ignore_index=True).to_csv(
        outdir / "rsf_partial_dependence.csv", index=False)
    plot_partial_dependence(pdps, figdir / "rsf_partial_dependence.png", cohort)

    (outdir / "rsf_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(f"[done]   outputs in {outdir}")
    print(f"         figures in {figdir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default="genie", choices=sorted(COHORTS))
    ap.add_argument("--data", default=None,
                    help="analytic frame (.csv/.xlsx); default = the cohort's processed CSV")
    ap.add_argument("--genes-file", default=None,
                    help="selected_genes.txt from select_risk_model_genes.py")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--figdir", default=None)
    ap.add_argument("--test-size", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-estimators", type=int, default=200)
    ap.add_argument("--no-tune", action="store_true",
                    help="skip the 5-fold CV search and fit the base configuration")
    args = ap.parse_args()

    outdir = Path(args.outdir) if args.outdir else MODEL_OUT_DIR / args.cohort / "rsf"
    figdir = Path(args.figdir) if args.figdir else FIG_OUT_DIR / args.cohort / "rsf"
    run(args.cohort, args.data, args.genes_file, outdir, figdir,
        args.test_size, args.seed, args.n_estimators, not args.no_tune)


if __name__ == "__main__":
    main()
