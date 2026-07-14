"""Aim 1: top-5 gene + Sanchez-Vega pathway comparison between CNS-met and
no-CNS-met groups, plus forest plot of prevalence ratios and an oncoprint of
the top-5 split by any_brain_met.

Per-cohort top-5 is read from extracted_variables_<cohort>_top_genes.txt
(this was derived in the brain-mets-ever sub-cohort during ETL).
"""
from __future__ import annotations

import argparse
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from _lib import (
    COHORTS, PATHWAY_COLS, ensure_dirs, load_cohort, load_top_genes,
)


def two_by_two(df: pd.DataFrame, gene_col: str) -> dict:
    """2x2 table for a binary gene column vs any_brain_met."""
    sub = df[[gene_col, "any_brain_met"]].dropna()
    sub = sub[sub[gene_col].isin([0, 1]) & sub["any_brain_met"].isin([0, 1])]
    a = int(((sub["any_brain_met"] == 1) & (sub[gene_col] == 1)).sum())  # met+ gene+
    b = int(((sub["any_brain_met"] == 1) & (sub[gene_col] == 0)).sum())  # met+ gene-
    c = int(((sub["any_brain_met"] == 0) & (sub[gene_col] == 1)).sum())  # met- gene+
    d = int(((sub["any_brain_met"] == 0) & (sub[gene_col] == 0)).sum())  # met- gene-

    p1 = a / (a + b) if (a + b) > 0 else np.nan
    p0 = c / (c + d) if (c + d) > 0 else np.nan
    # Prevalence ratio with log-normal 95% CI
    if a > 0 and c > 0:
        log_pr = np.log(p1 / p0)
        se = np.sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + d))
        pr = float(np.exp(log_pr))
        pr_lo = float(np.exp(log_pr - 1.96 * se))
        pr_hi = float(np.exp(log_pr + 1.96 * se))
    else:
        pr, pr_lo, pr_hi = np.nan, np.nan, np.nan

    try:
        _, pval = fisher_exact([[a, b], [c, d]])
    except Exception:
        pval = np.nan
    return {
        "feature": gene_col,
        "n_brain_met": a + b,
        "n_no_brain_met": c + d,
        "prev_brain_met": p1,
        "prev_no_brain_met": p0,
        "prev_ratio": pr,
        "pr_lo": pr_lo,
        "pr_hi": pr_hi,
        "fisher_p": pval,
        "count_brain_met_mut": a,
        "count_no_brain_met_mut": c,
    }


def aim1(cohort: str) -> None:
    df = load_cohort(cohort)
    top5 = load_top_genes(cohort, n=5)
    print(f"[{cohort}] top-5 from disk: {top5}")
    print(f"[{cohort}] n={len(df)}  brain-met={int(df['any_brain_met'].sum())}")

    model_dir, fig_dir = ensure_dirs(cohort, "aim1")

    # Build feature list: top-5 gene cols (bare HUGO) + pathway cols present
    features = []
    for g in top5:
        if g not in df.columns:
            print(f"   WARN: gene col {g!r} not found in cohort {cohort}, skipping")
            continue
        features.append(("gene", g))
    for p in PATHWAY_COLS:
        if p in df.columns:
            features.append(("pathway", p))
        else:
            print(f"   WARN: pathway col {p!r} not found in cohort {cohort}, skipping")

    # Two-prop / Fisher per feature
    rows = []
    for kind, col in features:
        r = two_by_two(df, col)
        r["kind"] = kind
        rows.append(r)
    res = pd.DataFrame(rows)

    # Multiple-testing correction WITHIN each kind (top-5 = 5 tests; pathways = 10)
    res["p_holm"] = np.nan
    res["p_bh"] = np.nan
    for k in res["kind"].unique():
        mask = res["kind"] == k
        ps = res.loc[mask, "fisher_p"].fillna(1.0).values
        if len(ps) > 0:
            _, holm, _, _ = multipletests(ps, method="holm")
            _, bh, _, _ = multipletests(ps, method="fdr_bh")
            res.loc[mask, "p_holm"] = holm
            res.loc[mask, "p_bh"] = bh

    out_csv = model_dir / "aim1_gene_pathway_comparison.csv"
    res = res[[
        "kind", "feature", "n_brain_met", "n_no_brain_met",
        "count_brain_met_mut", "count_no_brain_met_mut",
        "prev_brain_met", "prev_no_brain_met",
        "prev_ratio", "pr_lo", "pr_hi",
        "fisher_p", "p_holm", "p_bh",
    ]]
    res.to_csv(out_csv, index=False)
    print(f"   wrote {out_csv}")

    # Forest plot of prevalence ratios
    plot_df = res.dropna(subset=["prev_ratio"]).copy()
    if len(plot_df) > 0:
        plot_df = plot_df.sort_values(["kind", "prev_ratio"], ascending=[True, True])
        fig, ax = plt.subplots(figsize=(7, max(3, 0.35 * len(plot_df))))
        y = np.arange(len(plot_df))
        ax.errorbar(plot_df["prev_ratio"], y,
                    xerr=[plot_df["prev_ratio"] - plot_df["pr_lo"],
                          plot_df["pr_hi"] - plot_df["prev_ratio"]],
                    fmt="o", color="black", ecolor="gray", capsize=2)
        ax.axvline(1, color="red", linestyle="--", linewidth=0.8)
        labels = [f"{r['feature']}  (p={r['fisher_p']:.3g})" for _, r in plot_df.iterrows()]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xscale("log")
        ax.set_xlabel("Prevalence ratio (brain-met / no-brain-met), 95% CI")
        ax.set_title(f"Aim 1 prevalence ratios - cohort: {cohort}")
        fig.tight_layout()
        out_fig = fig_dir / "aim1_forest.png"
        fig.savefig(out_fig, dpi=150)
        plt.close(fig)
        print(f"   wrote {out_fig}")

    # Oncoprint of top-5 genes split by any_brain_met
    gene_cols = [g for g in top5 if g in df.columns]
    if gene_cols:
        sub = df[gene_cols + ["any_brain_met"]].dropna(subset=["any_brain_met"]).copy()
        sub["any_brain_met"] = sub["any_brain_met"].astype(int)
        sub = sub.sort_values(["any_brain_met"] + gene_cols, ascending=[True] + [False] * len(gene_cols))
        mat = sub[gene_cols].fillna(0).astype(int).values.T  # genes x samples
        fig, (ax_top, ax) = plt.subplots(
            2, 1, figsize=(12, max(2.5, 0.4 * len(gene_cols) + 1.5)),
            gridspec_kw={"height_ratios": [0.5, 4]}, sharex=True,
        )
        ax_top.imshow(sub["any_brain_met"].values.reshape(1, -1),
                      aspect="auto", cmap="Reds", vmin=0, vmax=1)
        ax_top.set_yticks([0]); ax_top.set_yticklabels(["any_brain_met"], fontsize=8)
        ax_top.set_xticks([])
        ax.imshow(mat, aspect="auto", cmap="Greys", vmin=0, vmax=1)
        ax.set_yticks(range(len(gene_cols)))
        ax.set_yticklabels(gene_cols, fontsize=9)
        ax.set_xlabel(f"samples (n={len(sub)}, sorted by any_brain_met then top-5 status)")
        ax.set_xticks([])
        fig.suptitle(f"Aim 1 oncoprint top-5 - cohort: {cohort}", fontsize=10)
        fig.tight_layout()
        out_fig = fig_dir / "aim1_oncoprint.png"
        fig.savefig(out_fig, dpi=150)
        plt.close(fig)
        print(f"   wrote {out_fig}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", choices=list(COHORTS.keys()) + ["all"], default="all")
    args = p.parse_args()
    cohorts = list(COHORTS.keys()) if args.cohort == "all" else [args.cohort]
    for c in cohorts:
        print(f"\n=== Aim 1 cohort={c} ===")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            aim1(c)


if __name__ == "__main__":
    main()
