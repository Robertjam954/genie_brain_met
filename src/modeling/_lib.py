"""Shared utilities for the brain-mets analytic pipeline.

Reads the extracted_variables_<cohort>_data.csv frames produced by ETL and
returns harmonized analytic frames. All three cohorts share the canonical
column names (any_brain_met, brain_met_at_dx, OS_months, os_status_bin,
top5_any_mutated, receptor_primary_cat, ...).
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "processed"
MODEL_OUT_DIR = PROJECT_ROOT / "src" / "modeling"
FIG_OUT_DIR = PROJECT_ROOT / "manuscript components"

PATHWAY_COLS = [
    "pathway_RTK/RAS", "pathway_Nrf2", "pathway_PI3K", "pathway_TGFB",
    "pathway_p53", "pathway_Wnt", "pathway_Myc", "pathway_Cell cycle",
    "pathway_Hippo", "pathway_Notch",
]


@dataclass
class CohortSpec:
    name: str
    data_file: str
    top_genes_file: str


COHORTS: dict[str, CohortSpec] = {
    "genie": CohortSpec(
        name="genie",
        data_file="extracted_variables_genie_data.csv",
        top_genes_file="extracted_variables_genie_top_genes.txt",
    ),
    "tcga": CohortSpec(
        name="tcga",
        data_file="extracted_variables_tcga_data.csv",
        top_genes_file="extracted_variables_tcga_top_genes.txt",
    ),
    "msk18": CohortSpec(
        name="msk18",
        data_file="extracted_variables_breast_msk_2018_data.csv",
        top_genes_file="extracted_variables_breast_msk_2018_top_genes.txt",
    ),
}


def load_cohort(cohort: str) -> pd.DataFrame:
    spec = COHORTS[cohort]
    path = DATA_DIR / spec.data_file
    df = pd.read_csv(path, low_memory=False)
    return df


def load_top_genes(cohort: str, n: int = 5) -> list[str]:
    spec = COHORTS[cohort]
    path = DATA_DIR / spec.top_genes_file
    genes: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genes.append(line)
            if len(genes) >= n:
                break
    return genes


def ensure_dirs(cohort: str, aim: str) -> tuple[Path, Path]:
    model_dir = MODEL_OUT_DIR / cohort / aim
    fig_dir = FIG_OUT_DIR / cohort / aim
    model_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    return model_dir, fig_dir


def safe_int(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").astype("Int64")


def prep_covariates(df: pd.DataFrame, include_stage_iv: bool = True) -> pd.DataFrame:
    """One-hot encode categoricals and coerce numerics for Cox/AFT covariates.

    Returns a covariate-only frame; the caller supplies time/event columns.
    """
    cov_cols = [
        "top5_any_mutated", "receptor_primary_cat", "grade_ord", "age_dx_num",
        "race_clean", "ethnicity_clean",
    ]
    if include_stage_iv and "stage_iv_bin" in df.columns:
        cov_cols.append("stage_iv_bin")
    if "SEQ_ASSAY_ID" in df.columns:
        cov_cols.append("SEQ_ASSAY_ID")

    cov = df[cov_cols].copy()
    cov["top5_any_mutated"] = pd.to_numeric(cov["top5_any_mutated"], errors="coerce")
    if "stage_iv_bin" in cov.columns:
        cov["stage_iv_bin"] = pd.to_numeric(cov["stage_iv_bin"], errors="coerce")
    cov["age_dx_num"] = pd.to_numeric(cov["age_dx_num"], errors="coerce")

    # Categoricals -> one-hot, drop first level
    cat_cols = [c for c in ["receptor_primary_cat", "grade_ord", "race_clean",
                            "ethnicity_clean", "SEQ_ASSAY_ID"] if c in cov.columns]
    cov = pd.get_dummies(cov, columns=cat_cols, drop_first=True, dummy_na=False)
    # Cast bool dummies to int
    bool_cols = cov.select_dtypes(include=["bool"]).columns
    for c in bool_cols:
        cov[c] = cov[c].astype(int)
    return cov


def drop_low_variance(cov: pd.DataFrame, min_unique: int = 2) -> pd.DataFrame:
    """Drop constant columns and columns with one-hot level n < 5 (to avoid
    Cox/AFT singularities in small cohorts)."""
    keep = []
    for c in cov.columns:
        s = cov[c]
        nu = s.dropna().nunique()
        if nu < min_unique:
            continue
        if set(s.dropna().unique()).issubset({0, 1}) and s.sum() < 5:
            continue
        keep.append(c)
    return cov[keep]
