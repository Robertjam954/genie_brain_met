"""Step 5: feature selection from the step-4 risk models.

Selection rule (as specified in the analysis order of operations): keep a gene if

    q < 0.05 in the Fine-Gray subdistribution model
      OR
    (hazard_ratio > 1 AND p < 0.05) in the Cox PH overall-survival model

Reads the per-gene tables written by `finegray_cox_risk_models.R`:
    <risk_dir>/finegray_gene_subhazards.csv   gene, subhazard_ratio, p_value, q_value
    <risk_dir>/cox_os_gene_hazards.csv        gene, hazard_ratio, p_value, q_value

Writes:
    <risk_dir>/selected_genes.csv   one row per gene with the reason it was kept
    <risk_dir>/selected_genes.txt   bare gene list, consumed by rsf_time_to_brain_met.py

Usage:
    python3 "src/modeling/select_risk_model_genes.py" --cohort genie
    python3 "src/modeling/select_risk_model_genes.py" --risk-dir path/to/risk_models
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from _lib import MODEL_OUT_DIR

FG_FILE = "finegray_gene_subhazards.csv"
COX_FILE = "cox_os_gene_hazards.csv"
OUT_CSV = "selected_genes.csv"
OUT_TXT = "selected_genes.txt"

FG_Q_MAX = 0.05
COX_P_MAX = 0.05
COX_HR_MIN = 1.0


def _read_optional(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        print(f"   note: {path.name} not found in {path.parent} - skipping that criterion")
        return None
    df = pd.read_csv(path)
    if "gene" not in df.columns:
        raise ValueError(f"{path} has no 'gene' column")
    return df


def select_genes(
    fg: pd.DataFrame | None,
    cox: pd.DataFrame | None,
    fg_q_max: float = FG_Q_MAX,
    cox_p_max: float = COX_P_MAX,
    cox_hr_min: float = COX_HR_MIN,
) -> pd.DataFrame:
    """Apply the step-5 rule and return one row per selected gene."""
    rows: dict[str, dict] = {}

    if fg is not None and "q_value" in fg.columns:
        hit = fg[pd.to_numeric(fg["q_value"], errors="coerce") < fg_q_max]
        for _, r in hit.iterrows():
            rows.setdefault(r["gene"], {"gene": r["gene"]}).update({
                "finegray_subhazard_ratio": r.get("subhazard_ratio"),
                "finegray_p_value": r.get("p_value"),
                "finegray_q_value": r.get("q_value"),
                "selected_by_finegray": True,
            })

    if cox is not None and {"hazard_ratio", "p_value"} <= set(cox.columns):
        hr = pd.to_numeric(cox["hazard_ratio"], errors="coerce")
        p = pd.to_numeric(cox["p_value"], errors="coerce")
        hit = cox[(hr > cox_hr_min) & (p < cox_p_max)]
        for _, r in hit.iterrows():
            rows.setdefault(r["gene"], {"gene": r["gene"]}).update({
                "cox_hazard_ratio": r.get("hazard_ratio"),
                "cox_p_value": r.get("p_value"),
                "cox_q_value": r.get("q_value"),
                "selected_by_cox": True,
            })

    if not rows:
        return pd.DataFrame(columns=["gene", "selected_by_finegray", "selected_by_cox"])

    out = pd.DataFrame(list(rows.values()))
    for col in ("selected_by_finegray", "selected_by_cox"):
        if col not in out.columns:
            out[col] = False
        out[col] = out[col].fillna(False).astype(bool)
    lead = ["gene", "selected_by_finegray", "selected_by_cox"]
    out = out[lead + [c for c in out.columns if c not in lead]]
    sort_key = out["finegray_q_value"] if "finegray_q_value" in out.columns else out["gene"]
    return out.assign(_k=sort_key).sort_values("_k", na_position="last").drop(columns="_k")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default="genie",
                    help="cohort key; risk dir defaults to src/modeling/<cohort>/risk_models")
    ap.add_argument("--risk-dir", default=None,
                    help="directory holding the step-4 per-gene tables")
    ap.add_argument("--fg-q-max", type=float, default=FG_Q_MAX)
    ap.add_argument("--cox-p-max", type=float, default=COX_P_MAX)
    ap.add_argument("--cox-hr-min", type=float, default=COX_HR_MIN)
    args = ap.parse_args()

    risk_dir = Path(args.risk_dir) if args.risk_dir else MODEL_OUT_DIR / args.cohort / "risk_models"
    if not risk_dir.exists():
        raise SystemExit(f"risk-model directory not found: {risk_dir}\n"
                         f"run finegray_cox_risk_models.R first")
    print(f"[step 5] risk-model tables in {risk_dir}")

    fg = _read_optional(risk_dir / FG_FILE)
    cox = _read_optional(risk_dir / COX_FILE)
    if fg is None and cox is None:
        raise SystemExit("neither per-gene table is present - nothing to select from")

    sel = select_genes(fg, cox, args.fg_q_max, args.cox_p_max, args.cox_hr_min)
    sel.to_csv(risk_dir / OUT_CSV, index=False)
    (risk_dir / OUT_TXT).write_text("\n".join(sel["gene"].astype(str)) + ("\n" if len(sel) else ""))

    n_fg = int(sel.get("selected_by_finegray", pd.Series(dtype=bool)).sum())
    n_cox = int(sel.get("selected_by_cox", pd.Series(dtype=bool)).sum())
    print(f"   selected {len(sel)} genes  (Fine-Gray q<{args.fg_q_max}: {n_fg}; "
          f"Cox HR>{args.cox_hr_min} & p<{args.cox_p_max}: {n_cox})")
    for g in sel["gene"]:
        print(f"     {g}")
    print(f"   wrote {risk_dir / OUT_CSV}")
    print(f"   wrote {risk_dir / OUT_TXT}")


if __name__ == "__main__":
    main()
