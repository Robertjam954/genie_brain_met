"""Missingness audit for the finalized analytic frame.

Reports per-variable NA counts and percentages across:
  - the full frame
  - any_brain_met == 1 sub-cohort
  - any_brain_met == 0 sub-cohort

Covered variable groups (Table 1 + modeling covariates + outcome + time/event):
  - identifiers / cohort flags
  - demographics
  - tumor / receptor / stage / grade
  - assay
  - mutation burden
  - outcome / time-to-event columns (Aim 2 + Aim 3)
  - top-10 gene indicators + pathway indicators

Writes:
  src/modeling/<cohort>/aim1/aim1_missingness_audit.csv
  manuscript components/<cohort>/aim1/aim1_missingness_audit.{png,pdf}
"""
from __future__ import annotations

import argparse
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _lib import COHORTS, PATHWAY_COLS, ensure_dirs, load_cohort, load_top_genes


VARIABLE_GROUPS: list[tuple[str, list[str]]] = [
    ("Identifier / cohort", [
        "record_id", "SAMPLE_ID", "ca_seq",
        "any_brain_met", "brain_met_at_dx",
    ]),
    ("Demographics", [
        "age_dx_num", "age_cat",
        "race_clean", "ethnicity_clean",
    ]),
    ("Tumor / receptor / stage / grade", [
        "receptor_primary_cat", "her2_status_bin",
        "stage_diag_group", "stage_iv_bin",
        "grade_ord", "sample_type_bin",
    ]),
    ("Assay", [
        "SEQ_ASSAY_ID",
    ]),
    ("Mutation burden", [
        "mutation_count_all_sites_sum",
        "t_alt_count_max",
        "mutation_count_q",
        "t_alt_count_q",
    ]),
    ("Outcome / time-to-event (Aim 2)", [
        "tt_brain_met_mos", "brain_met_event",
    ]),
    ("Outcome / time-to-event (Aim 3)", [
        "OS_months", "os_status_bin",
        "PFS_imaging_months", "pfs_i_event_bin",
        "PFS_medonc_months", "pfs_m_event_bin",
    ]),
]


def _is_missing(s: pd.Series) -> pd.Series:
    """A value is missing if NaN, or an empty / sentinel string."""
    if s.dtype == object:
        sl = s.astype(str).str.strip().str.lower()
        sentinel = {"", "nan", "na", "n/a", "none", "null", "unknown", "not collected",
                    "not applicable", "not assessed", "not reported"}
        return s.isna() | sl.isin(sentinel)
    return s.isna()


def audit(df: pd.DataFrame, var_groups: list[tuple[str, list[str]]]) -> pd.DataFrame:
    n_total = len(df)
    bm_mask = df["any_brain_met"].astype("Int64") == 1
    nbm_mask = df["any_brain_met"].astype("Int64") == 0
    n_bm = int(bm_mask.sum())
    n_nbm = int(nbm_mask.sum())

    rows = []
    for group, cols in var_groups:
        for c in cols:
            if c not in df.columns:
                rows.append({
                    "group": group, "variable": c,
                    "present": False,
                    "n_total": n_total,
                    "n_missing_total": np.nan, "pct_missing_total": np.nan,
                    "n_missing_bm": np.nan, "pct_missing_bm": np.nan,
                    "n_missing_nbm": np.nan, "pct_missing_nbm": np.nan,
                    "dtype": "",
                    "n_unique_observed": np.nan,
                })
                continue
            s = df[c]
            miss = _is_missing(s)
            n_miss_t = int(miss.sum())
            n_miss_bm = int((miss & bm_mask).sum())
            n_miss_nbm = int((miss & nbm_mask).sum())
            rows.append({
                "group": group, "variable": c,
                "present": True,
                "n_total": n_total,
                "n_missing_total": n_miss_t,
                "pct_missing_total": 100 * n_miss_t / n_total if n_total else np.nan,
                "n_missing_bm": n_miss_bm,
                "pct_missing_bm": 100 * n_miss_bm / n_bm if n_bm else np.nan,
                "n_missing_nbm": n_miss_nbm,
                "pct_missing_nbm": 100 * n_miss_nbm / n_nbm if n_nbm else np.nan,
                "dtype": str(s.dtype),
                "n_unique_observed": int(s[~miss].nunique()),
            })
    out = pd.DataFrame(rows)
    out.attrs.update({"n_total": n_total, "n_bm": n_bm, "n_nbm": n_nbm})
    return out


def _fmt_pct(n: float, pct: float) -> str:
    if pd.isna(n) or pd.isna(pct):
        return "-"
    return f"{int(n):>4d} ({pct:5.1f}%)"


def render_table(audit_df: pd.DataFrame, title: str, out_png, out_pdf) -> None:
    cols = [
        "Group", "Variable", "Present", "Dtype", "Unique obs",
        "Missing total (n, %)",
        "Missing brain-met (n, %)",
        "Missing no-brain-met (n, %)",
    ]
    rows = []
    for _, r in audit_df.iterrows():
        if not r["present"]:
            rows.append([
                r["group"], r["variable"], "MISSING COL", "-", "-", "-", "-", "-",
            ])
            continue
        rows.append([
            r["group"], r["variable"], "yes", r["dtype"], int(r["n_unique_observed"]),
            _fmt_pct(r["n_missing_total"], r["pct_missing_total"]),
            _fmt_pct(r["n_missing_bm"], r["pct_missing_bm"]),
            _fmt_pct(r["n_missing_nbm"], r["pct_missing_nbm"]),
        ])

    nrows = len(rows)
    fig_w = 16
    fig_h = 1.4 + 0.36 * nrows
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    table = ax.table(cellText=rows, colLabels=cols, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 1.35)
    for j in range(len(cols)):
        cell = table[(0, j)]
        cell.set_facecolor("#1f3a5f")
        cell.set_text_props(color="white", weight="bold")

    # Color cells: red for absent column; tint by total missingness severity
    for i, r in enumerate(audit_df.itertuples(index=False), start=1):
        if not r.present:
            for j in range(len(cols)):
                table[(i, j)].set_facecolor("#fad6d6")
            continue
        pct = r.pct_missing_total
        if pd.notna(pct):
            if pct >= 50:
                color = "#fad6d6"
            elif pct >= 20:
                color = "#ffe7b3"
            elif pct >= 5:
                color = "#fff7d6"
            elif pct > 0:
                color = "#eef7ee"
            else:
                color = "#d4f4d2"
            for j in range(5, 8):
                table[(i, j)].set_facecolor(color)
        if i % 2 == 0:
            for j in range(5):
                if table[(i, j)].get_facecolor() == (1.0, 1.0, 1.0, 1.0):
                    table[(i, j)].set_facecolor("#f4f6fa")

    # Group separator visual (bold first cell when group changes)
    last_group = None
    for i, r in enumerate(audit_df.itertuples(index=False), start=1):
        if r.group != last_group:
            table[(i, 0)].set_text_props(weight="bold")
            last_group = r.group
        else:
            table[(i, 0)].get_text().set_text("")

    ax.set_title(title, fontsize=12, weight="bold", loc="left", pad=14)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def run(cohort: str) -> None:
    df = load_cohort(cohort)
    n_total = len(df)
    n_bm = int(df["any_brain_met"].sum())
    n_nbm = n_total - n_bm
    print(f"[{cohort}] n_total={n_total}  brain-met={n_bm}  no-brain-met={n_nbm}")

    # Add genomic feature group from the on-disk top-10 and the canonical pathways
    top10 = load_top_genes(cohort, n=10)
    groups = list(VARIABLE_GROUPS) + [
        ("Top-10 genes (mutation indicator)", top10),
        ("Sanchez-Vega pathways", PATHWAY_COLS),
    ]

    audit_df = audit(df, groups)

    model_dir, fig_dir = ensure_dirs(cohort, "aim1")
    out_csv = model_dir / "aim1_missingness_audit.csv"
    audit_df.to_csv(out_csv, index=False)
    print(f"   wrote {out_csv}")

    # Console summary - just the rows with > 0% missing in the full frame
    print("\n   Variables with >0% missing (full frame):")
    summary = audit_df[audit_df["present"] & (audit_df["n_missing_total"] > 0)]\
        [["group", "variable", "n_missing_total", "pct_missing_total",
          "pct_missing_bm", "pct_missing_nbm"]]
    if len(summary) == 0:
        print("      (none)")
    else:
        for _, r in summary.iterrows():
            print(
                f"      [{r['group']}] {r['variable']:<32s}  "
                f"total={int(r['n_missing_total']):>4d} ({r['pct_missing_total']:5.1f}%)   "
                f"bm={r['pct_missing_bm']:5.1f}%   no-bm={r['pct_missing_nbm']:5.1f}%"
            )

    # Console list of absent variables
    absent = audit_df[~audit_df["present"]]
    if len(absent) > 0:
        print("\n   Expected variables MISSING from the frame:")
        for _, r in absent.iterrows():
            print(f"      [{r['group']}] {r['variable']}")

    render_table(
        audit_df,
        title=(f"Aim 1 - Missingness audit "
               f"(cohort={cohort}; n={n_total}, brain-met={n_bm}, no-brain-met={n_nbm})"),
        out_png=fig_dir / "aim1_missingness_audit.png",
        out_pdf=fig_dir / "aim1_missingness_audit.pdf",
    )
    print(f"   wrote {fig_dir / 'aim1_missingness_audit.png'}")
    print(f"   wrote {fig_dir / 'aim1_missingness_audit.pdf'}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", choices=list(COHORTS.keys()) + ["all"], default="genie")
    args = p.parse_args()
    cohorts = list(COHORTS.keys()) if args.cohort == "all" else [args.cohort]
    for c in cohorts:
        print(f"\n=== Missingness audit cohort={c} ===")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run(c)


if __name__ == "__main__":
    main()
