"""Aim 1 supplementary: top-10 gene + Sanchez-Vega pathway tables with
p and q values, brain-met vs no-brain-met. Renders both CSV and a polished
PDF/PNG table in the spirit of TCGAsurvival MIA correlation tables.

Reads:
  data/processed/extracted_variables_<cohort>_data.csv
  data/processed/extracted_variables_<cohort>_top_genes.txt   (10 lines)

Writes (per cohort):
  src/modeling/<cohort>/aim1/aim1_top10_gene_pq.csv
  src/modeling/<cohort>/aim1/aim1_pathway_pq.csv
  manuscript components/<cohort>/aim1/aim1_top10_gene_pq.png
  manuscript components/<cohort>/aim1/aim1_pathway_pq.png
  manuscript components/<cohort>/aim1/aim1_top10_gene_pq.pdf
  manuscript components/<cohort>/aim1/aim1_pathway_pq.pdf
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

from _lib import COHORTS, PATHWAY_COLS, ensure_dirs, load_cohort, load_top_genes


def two_by_two(df: pd.DataFrame, col: str) -> dict:
    sub = df[[col, "any_brain_met"]].dropna()
    sub = sub[sub[col].isin([0, 1]) & sub["any_brain_met"].isin([0, 1])]
    a = int(((sub["any_brain_met"] == 1) & (sub[col] == 1)).sum())
    b = int(((sub["any_brain_met"] == 1) & (sub[col] == 0)).sum())
    c = int(((sub["any_brain_met"] == 0) & (sub[col] == 1)).sum())
    d = int(((sub["any_brain_met"] == 0) & (sub[col] == 0)).sum())

    n1 = a + b
    n0 = c + d
    p1 = a / n1 if n1 else np.nan
    p0 = c / n0 if n0 else np.nan
    if a > 0 and c > 0:
        log_pr = np.log(p1 / p0)
        se = np.sqrt(1 / a - 1 / n1 + 1 / c - 1 / n0)
        pr = float(np.exp(log_pr))
        pr_lo = float(np.exp(log_pr - 1.96 * se))
        pr_hi = float(np.exp(log_pr + 1.96 * se))
    else:
        pr, pr_lo, pr_hi = np.nan, np.nan, np.nan
    odds, pval = fisher_exact([[a, b], [c, d]])
    return {
        "feature": col,
        "n_brain_met": n1,
        "n_no_brain_met": n0,
        "mut_brain_met": a,
        "mut_no_brain_met": c,
        "prev_brain_met": p1,
        "prev_no_brain_met": p0,
        "prev_ratio": pr,
        "pr_lo": pr_lo,
        "pr_hi": pr_hi,
        "odds_ratio": float(odds) if np.isfinite(odds) else np.nan,
        "fisher_p": float(pval),
    }


def add_corrections(res: pd.DataFrame) -> pd.DataFrame:
    ps = res["fisher_p"].fillna(1.0).values
    _, p_holm, _, _ = multipletests(ps, method="holm")
    _, p_bh, _, _ = multipletests(ps, method="fdr_bh")
    res = res.copy()
    res["p_holm"] = p_holm
    res["q_bh"] = p_bh
    return res


def fmt_p(p: float) -> str:
    if pd.isna(p):
        return "-"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def fmt_pct(p: float) -> str:
    if pd.isna(p):
        return "-"
    return f"{100 * p:.1f}%"


def render_table_png(res: pd.DataFrame, title: str, out_png, out_pdf,
                     feature_label: str = "Gene") -> None:
    cols = [
        feature_label,
        "Mut in brain-met (n/N)",
        "% brain-met",
        "Mut in no-brain-met (n/N)",
        "% no-brain-met",
        "Prev. ratio (95% CI)",
        "Fisher p",
        "Holm p",
        "BH q",
    ]
    rows = []
    for _, r in res.iterrows():
        nm1 = f"{r['mut_brain_met']}/{r['n_brain_met']}"
        nm0 = f"{r['mut_no_brain_met']}/{r['n_no_brain_met']}"
        pr_str = (
            f"{r['prev_ratio']:.2f} ({r['pr_lo']:.2f}-{r['pr_hi']:.2f})"
            if pd.notna(r["prev_ratio"])
            else "-"
        )
        rows.append([
            str(r["feature"]).replace("pathway_", ""),
            nm1,
            fmt_pct(r["prev_brain_met"]),
            nm0,
            fmt_pct(r["prev_no_brain_met"]),
            pr_str,
            fmt_p(r["fisher_p"]),
            fmt_p(r["p_holm"]),
            fmt_p(r["q_bh"]),
        ])

    nrows = len(rows)
    fig_w = 13
    fig_h = 1.2 + 0.42 * nrows
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    table = ax.table(
        cellText=rows, colLabels=cols, loc="center", cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.35)
    # Header styling
    for col_idx in range(len(cols)):
        cell = table[(0, col_idx)]
        cell.set_facecolor("#1f3a5f")
        cell.set_text_props(color="white", weight="bold")
    # Highlight significant BH q (<= 0.05) and Holm p (<= 0.05)
    bh_col = cols.index("BH q")
    holm_col = cols.index("Holm p")
    fp_col = cols.index("Fisher p")
    for i, r in enumerate(res.itertuples(index=False), start=1):
        if pd.notna(r.q_bh) and r.q_bh <= 0.05:
            table[(i, bh_col)].set_facecolor("#d4f4d2")
        if pd.notna(r.p_holm) and r.p_holm <= 0.05:
            table[(i, holm_col)].set_facecolor("#d4f4d2")
        if pd.notna(r.fisher_p) and r.fisher_p <= 0.05:
            table[(i, fp_col)].set_facecolor("#d4f4d2")
        # Alternate row shading
        if i % 2 == 0:
            for j in range(len(cols)):
                if table[(i, j)].get_facecolor() == (1.0, 1.0, 1.0, 1.0):
                    table[(i, j)].set_facecolor("#f4f6fa")

    ax.set_title(title, fontsize=12, weight="bold", loc="left", pad=14)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def aim1_top10(cohort: str) -> None:
    df = load_cohort(cohort)
    top10 = load_top_genes(cohort, n=10)
    n_total = len(df)
    n_bm = int(df["any_brain_met"].sum())
    print(f"[{cohort}] top-10 from disk: {top10}")
    print(f"[{cohort}] n={n_total}  brain-met={n_bm}  no-brain-met={n_total - n_bm}")

    model_dir, fig_dir = ensure_dirs(cohort, "aim1")

    # --- Genes ---
    gene_rows = []
    for g in top10:
        if g not in df.columns:
            print(f"   WARN: gene col {g!r} missing, skipping")
            continue
        gene_rows.append(two_by_two(df, g))
    gene_res = add_corrections(pd.DataFrame(gene_rows))
    gene_res = gene_res.sort_values("prev_brain_met", ascending=False).reset_index(drop=True)

    gene_csv = model_dir / "aim1_top10_gene_pq.csv"
    gene_res.to_csv(gene_csv, index=False)
    print(f"   wrote {gene_csv}")

    gene_png = fig_dir / "aim1_top10_gene_pq.png"
    gene_pdf = fig_dir / "aim1_top10_gene_pq.pdf"
    render_table_png(
        gene_res,
        title=(f"Aim 1 - Top 10 genes by brain-met prevalence (cohort={cohort}; "
               f"brain-met n={n_bm}, no-brain-met n={n_total - n_bm}); "
               f"Fisher exact, Holm and BH-FDR across 10 tests"),
        out_png=gene_png,
        out_pdf=gene_pdf,
        feature_label="Gene",
    )
    print(f"   wrote {gene_png}")
    print(f"   wrote {gene_pdf}")

    # --- Pathways ---
    pw_rows = []
    for p in PATHWAY_COLS:
        if p not in df.columns:
            print(f"   WARN: pathway col {p!r} missing, skipping")
            continue
        pw_rows.append(two_by_two(df, p))
    pw_res = add_corrections(pd.DataFrame(pw_rows))
    pw_res = pw_res.sort_values("prev_brain_met", ascending=False).reset_index(drop=True)

    pw_csv = model_dir / "aim1_pathway_pq.csv"
    pw_res.to_csv(pw_csv, index=False)
    print(f"   wrote {pw_csv}")

    pw_png = fig_dir / "aim1_pathway_pq.png"
    pw_pdf = fig_dir / "aim1_pathway_pq.pdf"
    render_table_png(
        pw_res,
        title=(f"Aim 1 - Sanchez-Vega pathways by brain-met prevalence (cohort={cohort}; "
               f"brain-met n={n_bm}, no-brain-met n={n_total - n_bm}); "
               f"Fisher exact, Holm and BH-FDR across 10 tests"),
        out_png=pw_png,
        out_pdf=pw_pdf,
        feature_label="Pathway",
    )
    print(f"   wrote {pw_png}")
    print(f"   wrote {pw_pdf}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", choices=list(COHORTS.keys()) + ["all"], default="genie")
    args = p.parse_args()
    cohorts = list(COHORTS.keys()) if args.cohort == "all" else [args.cohort]
    for c in cohorts:
        print(f"\n=== Aim 1 top-10 + pathway p/q tables cohort={c} ===")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            aim1_top10(c)


if __name__ == "__main__":
    main()
