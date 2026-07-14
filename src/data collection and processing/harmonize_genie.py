"""Build extracted_variables_genie_data.csv: analytic dataset for the
GENIE BPC brain-met analysis.

Inputs:
  data/processed/genie_bpc_v1_sample_master_full.csv  (clinical + gene_binary + pathways)
  data/processed/genie_bpc_v1_mutations.csv           (full MAF, for proper aggregation)

Pipeline:
  1. Aggregate MAF to per-sample mutation summary (count, t_alt_count_max, genes)
  2. Recode clinical variables per analytic plan
  3. Define brain-met cohort: dist_mets_brain_cns == 1 OR DMETS_DX_BRAIN == 'Yes'
  4. Find top 10 most-mutated genes in brain-met group; create top-5 binary
  5. Add time-to-event columns for OS, PFS_i, PFS_m, time-to-brain-met
  6. Write extracted_variables_genie_data.csv + top genes list + data dictionary

Cohort key: dist_mets_brain_cns == 1 (cancer-level flag, captures mets at dx OR
during study; matches user's spec).
"""

from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path(
    "/Users/robertjames/Documents/Documents - Robert’s iMac/Research Projects/"
    "MSKCC Research Fellowship/Projects/genie_tcga_impact_brain_mets"
)
PROC = PROJ / "data/processed"
OUT_DIR = PROJ / "src" / "exploratory data analysis"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT_DIR / "extracted_variables_genie_data.csv"
OUT_TOP = OUT_DIR / "extracted_variables_genie_top_genes.txt"
OUT_DICT = OUT_DIR / "extracted_variables_genie_dictionary.csv"


def main() -> None:
    print("=== load inputs ===")
    df = pd.read_csv(PROC / "genie_bpc_v1_sample_master_full.csv",
                     low_memory=False)
    print(f"  sample_master_full: {df.shape}")
    maf = pd.read_csv(PROC / "genie_bpc_v1_mutations.csv", low_memory=False)
    print(f"  mutations MAF:      {maf.shape}")

    # -------- 1. MAF -> per-sample summary ---------
    print("\n=== aggregate MAF to per-sample ===")
    # Filter to non-silent, somatic-style (Mutation_Status often missing in
    # GENIE BPC; keep all rows that pass Variant_Classification filter)
    DROP_VC = {"Silent", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
               "Intron", "RNA", "IGR"}
    maf_f = maf[~maf["Variant_Classification"].isin(DROP_VC)].copy()
    maf_f = maf_f[maf_f["hugoGeneSymbol"].notna()
                  & maf_f["sampleId"].notna()]
    print(f"  after non-silent filter: {len(maf_f)} (from {len(maf)})")

    mut_summary = (maf_f.assign(t_alt_count=pd.to_numeric(maf_f.get("t_alt_count"),
                                                          errors="coerce"))
                        .groupby("sampleId")
                        .agg(mutation_count_all_sites_sum=("hugoGeneSymbol", "size"),
                             t_alt_count_max=("t_alt_count", "max"),
                             genes=("hugoGeneSymbol",
                                    lambda s: ";".join(sorted(s.dropna().unique()))))
                        .reset_index()
                        .rename(columns={"sampleId": "SAMPLE_ID"}))
    print(f"  mut_summary: {mut_summary.shape}")

    df = df.merge(mut_summary, on="SAMPLE_ID", how="left")
    df["mutation_count_all_sites_sum"] = df["mutation_count_all_sites_sum"].fillna(0)
    df["genes"] = df["genes"].fillna("")

    # -------- 2. Variable recoding ---------
    print("\n=== recode clinical variables ===")
    dict_rows = []

    def add_dict(var, original_levels, mapped_levels, note=""):
        for o, m in zip(original_levels, mapped_levels):
            dict_rows.append({"variable": var, "original": str(o),
                              "mapped": str(m), "note": note})

    # 2a. sample_type_bin: 1 Primary, 2 Metastasis
    if "SAMPLE_TYPE_DETAILED" in df.columns:
        s = df["SAMPLE_TYPE_DETAILED"].astype(str).str.lower()
        df["sample_type_bin"] = np.where(s.str.contains("primary"), 1,
                                  np.where(s.str.contains("metast"), 2,
                                           np.float64("nan")))
        add_dict("sample_type_bin",
                 ["Primary tumor", "Metastasis"], [1, 2],
                 "from SAMPLE_TYPE_DETAILED")

    # 2b. age_cat: <50, 50-70, >70 (from age_dx)
    if "age_dx" in df.columns:
        df["age_dx_num"] = pd.to_numeric(df["age_dx"], errors="coerce")
        df["age_cat"] = pd.cut(df["age_dx_num"],
                               bins=[-np.inf, 50, 70, np.inf],
                               labels=["<50", "50-70", ">70"], right=False)
        add_dict("age_cat",
                 ["<50", "50-70 (50 ≤ age ≤ 70)", ">70"],
                 ["<50", "50-70", ">70"], "from age_dx")

    # 2c. grade ordered: Low/Intermediate/High (from ca_grade)
    if "ca_grade" in df.columns:
        g = df["ca_grade"].astype(str)
        df["grade_ord"] = np.where(g.str.contains("Low", case=False, na=False)
                                   | g.str.match(r"^\s*I\s*$"), "Low",
                            np.where(g.str.contains("Intermediate", case=False, na=False)
                                     | g.str.contains(r"^\s*II\s*$", na=False), "Intermediate",
                              np.where(g.str.contains("High", case=False, na=False)
                                       | g.str.contains(r"^\s*III\s*$", na=False), "High",
                                None)))
        df["grade_ord"] = pd.Categorical(df["grade_ord"],
                                         categories=["Low", "Intermediate", "High"],
                                         ordered=True)
        add_dict("grade_ord", ["I/Low", "II/Intermediate", "III/High"],
                 ["Low", "Intermediate", "High"], "from ca_grade")

    # 2d. stage groups (from stage_dx)
    if "stage_dx" in df.columns:
        s = df["stage_dx"].astype(str).str.upper().str.strip()
        df["stage_diag_group"] = np.where(s.isin(["I", "IA", "IB", "STAGE I"]), "Stage I",
                                   np.where(s.isin(["II", "IIA", "IIB", "STAGE II"]), "Stage II",
                                     np.where(s.isin(["III", "IIIA", "IIIB", "IIIC", "STAGE III"]), "Stage III",
                                       np.where(s.isin(["IV", "STAGE IV"]), "Stage IV",
                                         None))))
        df["stage_diag_group"] = pd.Categorical(df["stage_diag_group"],
                                                categories=["Stage I", "Stage II", "Stage III", "Stage IV"],
                                                ordered=True)
        add_dict("stage_diag_group", ["I/IA/IB", "II/IIA/IIB", "III/IIIA-C", "IV"],
                 ["Stage I", "Stage II", "Stage III", "Stage IV"], "from stage_dx")

    # 2e. stage_dx_iv binary (already present, just normalize to 0/1)
    if "stage_dx_iv" in df.columns:
        df["stage_iv_bin"] = (df["stage_dx_iv"].astype(str)
                               .str.strip().str.lower()
                               .isin(["stage iv", "iv", "yes", "true", "1"])
                               .astype(int))
        add_dict("stage_iv_bin", ["Stage IV", "Stage I-III"], [1, 0],
                 "from stage_dx_iv")

    # 2f. receptor_primary_cat: bca_subtype already 4-class HR/HER2
    if "bca_subtype" in df.columns:
        sub = df["bca_subtype"].astype(str)
        df["receptor_primary_cat"] = np.where(
            sub.isin(["HR+, HER2-", "HR+/HER2-"]), "HR+/HER2-",
            np.where(sub.isin(["HR+, HER2+", "HR+/HER2+"]), "HR+/HER2+",
              np.where(sub.isin(["HR-, HER2+", "HR-/HER2+"]), "HR-/HER2+",
                np.where(sub.isin(["TNBC", "Triple Negative", "HR-, HER2-", "HR-/HER2-"]),
                         "Triple Negative", None))))
        df["receptor_primary_cat"] = pd.Categorical(df["receptor_primary_cat"],
                                                    categories=["HR+/HER2-", "HR+/HER2+",
                                                                "HR-/HER2+", "Triple Negative"])
        add_dict("receptor_primary_cat",
                 ["HR+, HER2-", "HR+, HER2+", "HR-, HER2+", "TNBC"],
                 ["HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "Triple Negative"],
                 "from bca_subtype")

    # 2g. race: collapse PRIMARY_RACE
    if "PRIMARY_RACE" in df.columns:
        r = df["PRIMARY_RACE"].astype(str)
        df["race_clean"] = np.where(r.str.contains("White|Middle East", case=False, na=False),
                                    "White",
                              np.where(r.str.contains("Black|African", case=False, na=False),
                                       "Black",
                                np.where(r.str.contains("Asian|Indian|Chinese", case=False, na=False),
                                         "Asian",
                                  np.where(r.str.contains("Native|Alaska", case=False, na=False),
                                           "Native American", None))))
        add_dict("race_clean",
                 ["White / Middle Eastern", "Black or African American", "Asian", "Native American / Alaska"],
                 ["White", "Black", "Asian", "Native American"], "from PRIMARY_RACE")

    # 2h. ethnicity: Hispanic / Non-Hispanic
    if "ETHNICITY" in df.columns:
        e = df["ETHNICITY"].astype(str)
        df["ethnicity_clean"] = np.where(
            e.str.contains("Non-Spanish|Non-Hispanic", case=False, na=False), "Non-Hispanic",
            np.where(e.str.contains("Hispanic|Latino|Cuban|Mexican|Puerto", case=False, na=False),
                     "Hispanic", None))
        add_dict("ethnicity_clean",
                 ["Non-Spanish; Non-Hispanic", "Hispanic / Latino variants"],
                 ["Non-Hispanic", "Hispanic"], "from ETHNICITY")

    # 2i. HER2 binary (from ca_bca_her_summ)
    if "ca_bca_her_summ" in df.columns:
        h = df["ca_bca_her_summ"].astype(str)
        df["her2_status_bin"] = np.where(
            h.str.contains("positive|amplif|equivocal_pos", case=False, na=False), 1,
            np.where(h.str.contains("negative|not amplif", case=False, na=False), 0, np.float64("nan")))
        add_dict("her2_status_bin", ["Positive/elevated/amplified", "Negative"],
                 [1, 0], "from ca_bca_her_summ")

    # 2j. OS status factor (from os_dx_status)
    if "os_dx_status" in df.columns:
        s = df["os_dx_status"].astype(str)
        df["os_status_f"] = np.where(s.str.startswith("0"), "Alive",
                              np.where(s.str.startswith("1"), "Deceased", None))
        df["os_status_bin"] = pd.to_numeric(
            df["os_dx_status"].astype(str).str[0], errors="coerce")

    # 2k. PFS event (imaging)
    if "pfs_i_adv_status" in df.columns:
        df["pfs_i_event_bin"] = pd.to_numeric(
            df["pfs_i_adv_status"].astype(str).str[0], errors="coerce")
    if "pfs_m_adv_status" in df.columns:
        df["pfs_m_event_bin"] = pd.to_numeric(
            df["pfs_m_adv_status"].astype(str).str[0], errors="coerce")

    # 2l. mutation count quartiles
    df["mutation_count_q"] = pd.qcut(df["mutation_count_all_sites_sum"],
                                     q=4, labels=["Q1", "Q2", "Q3", "Q4"],
                                     duplicates="drop")

    # 2m. t_alt_count quartiles
    if "t_alt_count_max" in df.columns:
        df["t_alt_count_q"] = pd.qcut(df["t_alt_count_max"], q=4,
                                      labels=["Q1", "Q2", "Q3", "Q4"],
                                      duplicates="drop")

    # -------- 3. brain-met cohort indicator ---------
    print("\n=== build brain-met cohort indicator ===")
    # Any brain met (at dx or during study). dist_mets_brain_cns is the
    # cancer-level cumulative flag from cancer_level_dataset_index; it includes
    # at-dx mets. DMETS_DX_BRAIN is from clinical_patient (at-dx only).
    brain_overall = df["dist_mets_brain_cns"].fillna(0).astype(int) == 1
    brain_at_dx = (df.get("DMETS_DX_BRAIN", pd.Series([""]*len(df)))
                     .astype(str).str.lower() == "yes")
    df["any_brain_met"] = (brain_overall | brain_at_dx).astype(int)
    df["brain_met_at_dx"] = brain_at_dx.astype(int)

    print(f"  any_brain_met:    {df['any_brain_met'].sum()} samples ({df.loc[df['any_brain_met']==1,'record_id'].nunique()} patients)")
    print(f"  brain_met_at_dx:  {df['brain_met_at_dx'].sum()} samples")
    print(f"  no brain met:     {(df['any_brain_met']==0).sum()} samples")

    # met_loc factor: Brain / Other / None
    other_organ_cols = [c for c in df.columns if c.startswith("dist_mets_")
                        and c not in {"dist_mets_brain_cns"}
                        and not c.startswith("dx_to_dist_mets_")]
    has_other_met = df[other_organ_cols].fillna(0).astype(int).max(axis=1) == 1
    df["met_loc"] = np.where(df["any_brain_met"] == 1, "Brain",
                       np.where(has_other_met, "Other", "None"))
    df["met_loc"] = pd.Categorical(df["met_loc"],
                                   categories=["Brain", "Other", "None"])

    # -------- 4. Top 10 / Top 5 genes in brain-met cohort ---------
    print("\n=== top 10 mutated genes in brain-met cohort ===")
    # Identify gene-binary mutation cols: in the gnomeR output, mutation cols
    # are HUGO symbols from the MAF (no suffix). Restrict to columns whose
    # names appear in the MAF gene list AND that are integer 0/1 typed.
    maf_gene_set = set(maf["hugoGeneSymbol"].dropna().astype(str).unique())
    mut_gene_cols = [c for c in df.columns
                     if c in maf_gene_set
                     and pd.api.types.is_integer_dtype(df[c])
                     and not c.endswith((".Amp", ".Del", ".fus"))]
    print(f"  mutation gene cols (from MAF + int dtype): {len(mut_gene_cols)}")

    brain_df = df[df["any_brain_met"] == 1]
    if len(brain_df) > 0:
        gene_prev = (brain_df[mut_gene_cols].sum()
                              .sort_values(ascending=False)
                              .to_frame("n_brain_met_samples_mutated"))
        gene_prev["pct_brain_met"] = gene_prev["n_brain_met_samples_mutated"] / len(brain_df)
        gene_prev["n_total_samples_mutated"] = df[gene_prev.index.tolist()].sum().values
        gene_prev["pct_total"] = gene_prev["n_total_samples_mutated"] / len(df)
        gene_prev = gene_prev.reset_index().rename(columns={"index": "gene"})
        top10 = gene_prev.head(10)
        top10_genes = top10["gene"].tolist()
        top5_genes = top10_genes[:5]
        print("  top 10:")
        for _, row in top10.iterrows():
            print(f"    {row['gene']:8s} brain-met: {int(row['n_brain_met_samples_mutated']):>4} "
                  f"({row['pct_brain_met']*100:5.1f}%)   "
                  f"total: {int(row['n_total_samples_mutated']):>4} ({row['pct_total']*100:5.1f}%)")
        # Save full gene prevalence and top lists
        gene_prev.to_csv(OUT_DIR / "extracted_variables_genie_gene_prev_brain_met.csv", index=False)
        with open(OUT_TOP, "w") as f:
            f.write("# Top 10 genes by mutation prevalence in any_brain_met == 1 cohort\n")
            for g in top10_genes:
                f.write(g + "\n")

        # Binary indicators
        for g in top10_genes:
            df[f"G_top10_{g}"] = df[g].astype(int)
        df["top5_any_mutated"] = (df[top5_genes].sum(axis=1) > 0).astype(int)
        df["top10_any_mutated"] = (df[top10_genes].sum(axis=1) > 0).astype(int)
        df["top5_n_mutated"] = df[top5_genes].sum(axis=1).astype(int)
        df["top10_n_mutated"] = df[top10_genes].sum(axis=1).astype(int)
        print(f"  top5_any_mutated:  {df['top5_any_mutated'].sum()} samples")
        print(f"  top10_any_mutated: {df['top10_any_mutated'].sum()} samples")

    # -------- 5. Time-to-event columns ready for modeling ---------
    print("\n=== time-to-event columns ===")
    # Brain-met time-to-event: in patients without brain met, use overall OS time
    # as censoring time.
    df["tt_brain_met_mos"] = np.where(
        df["any_brain_met"] == 1,
        pd.to_numeric(df.get("dx_to_dist_mets_brain_cns_mos"), errors="coerce"),
        pd.to_numeric(df.get("tt_os_dx_mos"), errors="coerce")
    )
    df["brain_met_event"] = df["any_brain_met"]
    print(f"  tt_brain_met_mos: {df['tt_brain_met_mos'].notna().sum()} non-NA")
    print(f"  brain_met_event:  {df['brain_met_event'].sum()} events")

    # OS / PFS already present as tt_os_dx_mos / os_dx_status etc; just rename
    # standard analytic names for clarity
    rename_for_analysis = {
        "tt_os_dx_mos": "OS_months",
        "tt_pfs_i_adv_mos": "PFS_imaging_months",
        "tt_pfs_m_adv_mos": "PFS_medonc_months",
        "dx_to_dist_mets_brain_cns_mos": "time_to_brain_met_mos",
    }
    for old, new in rename_for_analysis.items():
        if old in df.columns and new not in df.columns:
            df[new] = df[old]

    # -------- 6. Write output ---------
    print("\n=== write ===")
    df.to_csv(OUT_CSV, index=False)
    print(f"  wrote: {OUT_CSV}  ({df.shape[0]} x {df.shape[1]})  "
          f"{OUT_CSV.stat().st_size/1e6:.1f} MB")

    pd.DataFrame(dict_rows).to_csv(OUT_DICT, index=False)
    print(f"  wrote: {OUT_DICT}  ({len(dict_rows)} dictionary rows)")

    print("\n=== summary table for analytic plan ===")
    if "top5_any_mutated" in df.columns:
        ct = pd.crosstab(df["any_brain_met"], df["top5_any_mutated"],
                         margins=True, margins_name="total")
        print("any_brain_met x top5_any_mutated (samples):")
        print(ct)
    print("\nmet_loc distribution:")
    print(df["met_loc"].value_counts(dropna=False))
    print("\nreceptor_primary_cat x any_brain_met (samples):")
    print(pd.crosstab(df.get("receptor_primary_cat"), df["any_brain_met"]))


if __name__ == "__main__":
    main()
