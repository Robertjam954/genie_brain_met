"""Aim 1 supplementary: univariate logistic regression of `any_brain_met` on
each top-10 gene and each Sanchez-Vega pathway, then a multivariate logistic
regression with the univariate-significant features.

Outcome:  any_brain_met (1 = brain mets ever, 0 = no brain mets)
Exposure: bare-HUGO gene indicators (top 10) + 10 pathway_<X> indicators

Univariate model:    logit(P(any_brain_met)) ~ feature
Multivariate model:  logit(P(any_brain_met)) ~ all features with univariate p < 0.05

Writes:
  src/modeling/<cohort>/aim1/aim1_univariate_logistic.csv
  src/modeling/<cohort>/aim1/aim1_multivariate_logistic.csv
  manuscript components/<cohort>/aim1/aim1_univariate_logistic.{png,pdf}
  manuscript components/<cohort>/aim1/aim1_multivariate_logistic.{png,pdf}
"""
from __future__ import annotations

import argparse
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from _lib import COHORTS, PATHWAY_COLS, ensure_dirs, load_cohort, load_top_genes


UNIVAR_ALPHA = 0.05


def _fit_logit(y: pd.Series, X: pd.DataFrame) -> sm.discrete.discrete_model.BinaryResultsWrapper:
    Xc = sm.add_constant(X.astype(float), has_constant="add")
    model = sm.Logit(y.astype(int).values, Xc.values)
    return model.fit(disp=False, maxiter=200)


def univariate(df: pd.DataFrame, feats: list[tuple[str, str]]) -> pd.DataFrame:
    rows = []
    y = df["any_brain_met"].astype(int)
    for kind, col in feats:
        s = df[col]
        sub_mask = s.isin([0, 1]) & y.isin([0, 1])
        sub_y = y[sub_mask]
        sub_x = s[sub_mask].astype(int)
        n_pos = int(sub_x.sum())
        n_neg = int(len(sub_x) - n_pos)
        if n_pos < 5 or n_neg < 5 or sub_y.nunique() < 2:
            rows.append({
                "kind": kind, "feature": col,
                "n": int(len(sub_x)), "n_mutated": n_pos,
                "or": np.nan, "or_lo": np.nan, "or_hi": np.nan,
                "coef": np.nan, "se": np.nan, "p": np.nan,
                "note": "skipped: low cell count",
            })
            continue
        try:
            res = _fit_logit(sub_y, sub_x.to_frame(name=col))
            coef = res.params[1]
            se = res.bse[1]
            p = res.pvalues[1]
            or_ = float(np.exp(coef))
            or_lo = float(np.exp(coef - 1.96 * se))
            or_hi = float(np.exp(coef + 1.96 * se))
            rows.append({
                "kind": kind, "feature": col,
                "n": int(len(sub_x)), "n_mutated": n_pos,
                "or": or_, "or_lo": or_lo, "or_hi": or_hi,
                "coef": float(coef), "se": float(se), "p": float(p),
                "note": "",
            })
        except Exception as e:
            rows.append({
                "kind": kind, "feature": col,
                "n": int(len(sub_x)), "n_mutated": n_pos,
                "or": np.nan, "or_lo": np.nan, "or_hi": np.nan,
                "coef": np.nan, "se": np.nan, "p": np.nan,
                "note": f"fit_error: {type(e).__name__}",
            })
    res_df = pd.DataFrame(rows)
    # BH within each kind (so 10 genes and 10 pathways are corrected separately)
    res_df["q_bh"] = np.nan
    for k in res_df["kind"].unique():
        m = (res_df["kind"] == k) & res_df["p"].notna()
        if m.sum() == 0:
            continue
        ps = res_df.loc[m, "p"].values
        _, q, _, _ = multipletests(ps, method="fdr_bh")
        res_df.loc[m, "q_bh"] = q
    return res_df


def multivariate(df: pd.DataFrame, feats: list[str]) -> pd.DataFrame:
    y = df["any_brain_met"].astype(int)
    X = df[feats].copy()
    for c in feats:
        X[c] = pd.to_numeric(X[c], errors="coerce")
    mask = y.isin([0, 1]) & X.notna().all(axis=1) & X.isin([0, 1]).all(axis=1)
    y_use = y[mask]
    X_use = X[mask].astype(int)
    if len(y_use) == 0 or y_use.nunique() < 2:
        raise RuntimeError("multivariate: no usable rows")
    res = _fit_logit(y_use, X_use)
    # First row of params is the intercept
    coefs = res.params[1:]
    ses = res.bse[1:]
    ps = res.pvalues[1:]
    rows = []
    for i, name in enumerate(feats):
        coef = float(coefs[i])
        se = float(ses[i])
        rows.append({
            "feature": name,
            "or": float(np.exp(coef)),
            "or_lo": float(np.exp(coef - 1.96 * se)),
            "or_hi": float(np.exp(coef + 1.96 * se)),
            "coef": coef,
            "se": se,
            "p": float(ps[i]),
        })
    out = pd.DataFrame(rows)
    out.attrs["n"] = int(len(y_use))
    out.attrs["n_brain_met"] = int(y_use.sum())
    out.attrs["pseudo_r2"] = float(res.prsquared)
    out.attrs["llf"] = float(res.llf)
    out.attrs["llnull"] = float(res.llnull)
    return out


def fmt_p(p: float) -> str:
    if pd.isna(p):
        return "-"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def render_table(res: pd.DataFrame, title: str, out_png, out_pdf,
                 cols_spec: list[tuple[str, str]],
                 sig_cols: list[str] | None = None,
                 sort_by: str | None = None) -> None:
    if sort_by and sort_by in res.columns:
        res = res.sort_values(sort_by, na_position="last").reset_index(drop=True)
    rows = []
    for _, r in res.iterrows():
        out_row = []
        for label, attr in cols_spec:
            if attr == "feature":
                out_row.append(str(r["feature"]).replace("pathway_", ""))
            elif attr == "or_ci":
                if pd.notna(r.get("or")) and pd.notna(r.get("or_lo")) and pd.notna(r.get("or_hi")):
                    out_row.append(f"{r['or']:.2f} ({r['or_lo']:.2f}-{r['or_hi']:.2f})")
                else:
                    out_row.append("-")
            elif attr in {"p", "q_bh"}:
                out_row.append(fmt_p(r.get(attr, np.nan)))
            elif attr == "kind":
                out_row.append(str(r.get("kind", "")))
            elif attr == "n_pair":
                out_row.append(f"{int(r['n_mutated'])}/{int(r['n'])}")
            else:
                v = r.get(attr, "")
                out_row.append(str(v) if pd.notna(v) else "-")
        rows.append(out_row)

    nrows = len(rows)
    ncols = len(cols_spec)
    fig_w = 1.6 * ncols
    fig_h = 1.2 + 0.42 * nrows
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    table = ax.table(
        cellText=rows,
        colLabels=[c[0] for c in cols_spec],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.4)
    for j in range(ncols):
        cell = table[(0, j)]
        cell.set_facecolor("#1f3a5f")
        cell.set_text_props(color="white", weight="bold")
    sig_cols = sig_cols or []
    sig_col_idx = [i for i, (label, attr) in enumerate(cols_spec) if attr in sig_cols]
    for i, r in enumerate(res.itertuples(index=False), start=1):
        for j_attr in sig_col_idx:
            attr = cols_spec[j_attr][1]
            v = getattr(r, attr, np.nan)
            if pd.notna(v) and v <= 0.05:
                table[(i, j_attr)].set_facecolor("#d4f4d2")
        if i % 2 == 0:
            for j in range(ncols):
                if table[(i, j)].get_facecolor() == (1.0, 1.0, 1.0, 1.0):
                    table[(i, j)].set_facecolor("#f4f6fa")
    ax.set_title(title, fontsize=12, weight="bold", loc="left", pad=14)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def run(cohort: str) -> None:
    df = load_cohort(cohort)
    top10 = load_top_genes(cohort, n=10)
    n_total = len(df)
    n_bm = int(df["any_brain_met"].sum())
    print(f"[{cohort}] n={n_total}  brain-met={n_bm}  no-brain-met={n_total - n_bm}")

    feats = []
    for g in top10:
        if g in df.columns:
            feats.append(("gene", g))
        else:
            print(f"   WARN: gene col {g!r} missing")
    for p in PATHWAY_COLS:
        if p in df.columns:
            feats.append(("pathway", p))
        else:
            print(f"   WARN: pathway col {p!r} missing")

    model_dir, fig_dir = ensure_dirs(cohort, "aim1")

    # --- Univariate ---
    uni = univariate(df, feats)
    uni_csv = model_dir / "aim1_univariate_logistic.csv"
    uni.to_csv(uni_csv, index=False)
    print(f"   wrote {uni_csv}")

    cols_spec_uni = [
        ("Kind", "kind"),
        ("Feature", "feature"),
        ("Mut/N", "n_pair"),
        ("OR (95% CI)", "or_ci"),
        ("p", "p"),
        ("BH q", "q_bh"),
    ]
    render_table(
        uni,
        title=(f"Aim 1 - Univariate logistic regression of any_brain_met "
               f"(cohort={cohort}; brain-met={n_bm}/{n_total})"),
        out_png=fig_dir / "aim1_univariate_logistic.png",
        out_pdf=fig_dir / "aim1_univariate_logistic.pdf",
        cols_spec=cols_spec_uni,
        sig_cols=["p", "q_bh"],
        sort_by=None,
    )
    print(f"   wrote {fig_dir / 'aim1_univariate_logistic.png'}")

    # --- Multivariate with univariate-significant features ---
    sig_mask = uni["p"].notna() & (uni["p"] < UNIVAR_ALPHA)
    sig_feats = uni.loc[sig_mask, "feature"].tolist()
    print(f"   univariate-significant (p<{UNIVAR_ALPHA}) features: {sig_feats}")

    if len(sig_feats) < 1:
        print("   no significant features; skipping multivariate")
        return
    if len(sig_feats) == 1:
        print("   only 1 significant feature; multivariate is identical to univariate")

    multi = multivariate(df, sig_feats)
    multi_csv = model_dir / "aim1_multivariate_logistic.csv"
    multi.to_csv(multi_csv, index=False)
    print(
        f"   wrote {multi_csv}  (n={multi.attrs['n']}, brain-met={multi.attrs['n_brain_met']}, "
        f"McFadden pseudo-R^2={multi.attrs['pseudo_r2']:.3f})"
    )

    cols_spec_mv = [
        ("Feature", "feature"),
        ("Adjusted OR (95% CI)", "or_ci"),
        ("p", "p"),
    ]
    render_table(
        multi,
        title=(f"Aim 1 - Multivariate logistic regression of any_brain_met "
               f"(cohort={cohort}; n={multi.attrs['n']}, brain-met={multi.attrs['n_brain_met']}; "
               f"features with univariate p<{UNIVAR_ALPHA})"),
        out_png=fig_dir / "aim1_multivariate_logistic.png",
        out_pdf=fig_dir / "aim1_multivariate_logistic.pdf",
        cols_spec=cols_spec_mv,
        sig_cols=["p"],
        sort_by="p",
    )
    print(f"   wrote {fig_dir / 'aim1_multivariate_logistic.png'}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", choices=list(COHORTS.keys()) + ["all"], default="genie")
    args = p.parse_args()
    cohorts = list(COHORTS.keys()) if args.cohort == "all" else [args.cohort]
    for c in cohorts:
        print(f"\n=== Aim 1 logistic regressions cohort={c} ===")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run(c)


if __name__ == "__main__":
    main()
