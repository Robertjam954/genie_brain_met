"""Aim 3 (OS only): overall survival within the brain-mets-ever cohort.

Per cohort: KM stratified by top5_any_mutated + log-rank; Cox PH with covariates
+ scaled Schoenfeld diagnostics; AFT distribution selection (Weibull / lognormal
/ log-logistic) via AIC; full AFT fit at chosen distribution. PFS endpoints
dropped per user spec - OS only across all 3 cohorts.
"""
from __future__ import annotations

import argparse
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import (
    CoxPHFitter, KaplanMeierFitter, LogLogisticAFTFitter, LogNormalAFTFitter,
    WeibullAFTFitter,
)
from lifelines.statistics import proportional_hazard_test, logrank_test

from _lib import COHORTS, ensure_dirs, drop_low_variance, load_cohort, prep_covariates


MIN_EVENTS = 25


def aim3(cohort: str) -> None:
    df = load_cohort(cohort)
    model_dir, fig_dir = ensure_dirs(cohort, "aim3")

    base = df[df["any_brain_met"] == 1].copy()
    print(f"[{cohort}] brain-mets-ever cohort n={len(base)}  OS events={int(base['os_status_bin'].sum())}")

    # ----- KM + log-rank -----
    km_df = base[["OS_months", "os_status_bin", "top5_any_mutated"]].dropna()
    km_df = km_df[km_df["OS_months"] > 0]
    fig, ax = plt.subplots(figsize=(7, 5))
    for grp, sub in km_df.groupby("top5_any_mutated"):
        kmf = KaplanMeierFitter()
        kmf.fit(sub["OS_months"], sub["os_status_bin"],
                label=f"top5_any_mutated={int(grp)}  (n={len(sub)}, events={int(sub.os_status_bin.sum())})")
        kmf.plot_survival_function(ax=ax, ci_show=True)
    g0 = km_df[km_df["top5_any_mutated"] == 0]
    g1 = km_df[km_df["top5_any_mutated"] == 1]
    if len(g0) and len(g1):
        lr = logrank_test(g0["OS_months"], g1["OS_months"],
                          g0["os_status_bin"], g1["os_status_bin"])
        ax.set_title(f"Aim 3 KM overall survival - {cohort}\nlog-rank p={lr.p_value:.4g}")
    else:
        ax.set_title(f"Aim 3 KM overall survival - {cohort}")
    ax.set_xlabel("Months from dx"); ax.set_ylabel("P(alive)")
    fig.tight_layout()
    out_fig = fig_dir / "aim3_os_km.png"
    fig.savefig(out_fig, dpi=150); plt.close(fig)
    print(f"   wrote {out_fig}")

    events = int(base["os_status_bin"].dropna().sum())
    if events < MIN_EVENTS:
        print(f"   SKIP Cox/AFT: only {events} OS events (< {MIN_EVENTS})")
        (model_dir / "aim3_skipped.txt").write_text(
            f"Skipped Cox/AFT - {events} OS events (< {MIN_EVENTS}).\n"
        )
        return

    # ----- Cox PH (no stage_iv_bin per README Sec 9.2) -----
    cov = prep_covariates(base, include_stage_iv=False)
    cph_df = pd.concat([
        base[["OS_months", "os_status_bin"]].reset_index(drop=True),
        cov.reset_index(drop=True),
    ], axis=1).dropna()
    cph_df = cph_df[cph_df["OS_months"] > 0]
    cph_cov = cph_df.drop(columns=["OS_months", "os_status_bin"])
    cph_cov = drop_low_variance(cph_cov)
    fit_df = pd.concat([
        cph_df[["OS_months", "os_status_bin"]].reset_index(drop=True),
        cph_cov.reset_index(drop=True),
    ], axis=1)
    print(f"   Cox: n={len(fit_df)} covariates={cph_cov.shape[1]}")

    cph = CoxPHFitter(penalizer=0.01)
    cph.fit(fit_df, duration_col="OS_months", event_col="os_status_bin",
            show_progress=False)
    cox_summary = cph.summary.reset_index().rename(columns={"covariate": "term"})
    cox_summary.to_csv(model_dir / "aim3_os_cox.csv", index=False)
    print(f"   wrote aim3_os_cox.csv  C-index={cph.concordance_index_:.3f}")

    try:
        zph = proportional_hazard_test(cph, fit_df, time_transform="rank")
        zph.summary.reset_index().rename(columns={"index": "term"}).to_csv(
            model_dir / "aim3_os_cox_zph.csv", index=False)
        print(f"   wrote aim3_os_cox_zph.csv")
    except Exception as e:
        print(f"   WARN: zph failed: {e}")

    # ----- AFT distribution selection -----
    aft_results = []
    fitters = {
        "weibull": WeibullAFTFitter(penalizer=0.01),
        "lognormal": LogNormalAFTFitter(penalizer=0.01),
        "loglogistic": LogLogisticAFTFitter(penalizer=0.01),
    }
    for name, f in fitters.items():
        try:
            f.fit(fit_df, duration_col="OS_months",
                  event_col="os_status_bin", show_progress=False)
            aft_results.append({"distribution": name, "AIC": f.AIC_,
                                "log_likelihood": f.log_likelihood_,
                                "concordance": f.concordance_index_})
        except Exception as e:
            print(f"   WARN: AFT {name} failed: {e}")
            aft_results.append({"distribution": name, "AIC": np.nan,
                                "log_likelihood": np.nan, "concordance": np.nan})
    aft_cmp = pd.DataFrame(aft_results).sort_values("AIC")
    aft_cmp.to_csv(model_dir / "aim3_os_aft_dist_comparison.csv", index=False)
    print(f"   AFT dist comparison:\n{aft_cmp.to_string(index=False)}")

    if aft_cmp["AIC"].notna().any():
        best = aft_cmp.iloc[0]["distribution"]
        best_fitter = fitters[best]
        try:
            aft_summary = best_fitter.summary.reset_index().rename(
                columns={"covariate": "term"})
            aft_summary["distribution"] = best
            aft_summary.to_csv(model_dir / "aim3_os_aft.csv", index=False)
            print(f"   wrote aim3_os_aft.csv  best={best}")
        except Exception as e:
            print(f"   WARN: AFT summary failed: {e}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", choices=list(COHORTS.keys()) + ["all"], default="all")
    args = p.parse_args()
    cohorts = list(COHORTS.keys()) if args.cohort == "all" else [args.cohort]
    for c in cohorts:
        print(f"\n=== Aim 3 OS cohort={c} ===")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            aim3(c)


if __name__ == "__main__":
    main()
