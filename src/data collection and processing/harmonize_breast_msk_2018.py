"""
Pull breast_msk_2018 (Razavi 2018 / MSK 2017-2018 BRCA cohort, N=1918) from
cBioPortal via pyBioPortal, and harmonize to the canonical schema in
harmonization_spec.md (column names match extracted_variables_genie_data.csv).

Outputs (data/processed/):
  extracted_variables_breast_msk_2018_data.csv
  extracted_variables_breast_msk_2018_top_genes.txt
  extracted_variables_breast_msk_2018_gene_prev_brain_met.csv
  extracted_variables_breast_msk_2018_dictionary.csv

Raw cache (data/raw/cbioportal/breast_msk_2018/) - re-runs skip re-download.
"""

from __future__ import annotations
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

# pyBioPortal + requests are only needed on cache miss. Import lazily so a
# fully-cached re-run works without the network deps installed.
sys.path.insert(0, "/Users/robertjames/Documents/Documents - Robert’s iMac/"
                  "pyBioPortal-master")
samples = clinical_data = mutations = cna_mod = None
def _ensure_pybioportal():
    global samples, clinical_data, mutations, cna_mod
    if samples is None:
        from pybioportal import (samples as _s, clinical_data as _c,
                                  mutations as _m,
                                  discrete_copy_number_alterations as _cna)
        samples, clinical_data, mutations, cna_mod = _s, _c, _m, _cna

STUDY = "breast_msk_2018"
ROOT = Path("/Users/robertjames/Documents/Documents - Robert’s iMac/"
            "Research Projects/MSKCC Research Fellowship/Projects/"
            "genie_tcga_impact_brain_mets")
RAW = ROOT / f"data/raw/cbioportal/{STUDY}"
RAW.mkdir(parents=True, exist_ok=True)
PROC = ROOT / "data/processed"
OUT_DIR = ROOT / "src" / "exploratory data analysis"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_MAIN = OUT_DIR / f"extracted_variables_{STUDY}_data.csv"
OUT_TOPGENES = OUT_DIR / f"extracted_variables_{STUDY}_top_genes.txt"
OUT_GENEPREV = OUT_DIR / f"extracted_variables_{STUDY}_gene_prev_brain_met.csv"
OUT_DICT = OUT_DIR / f"extracted_variables_{STUDY}_dictionary.csv"

DAYS_PER_MONTH = 30.4375
DAYS_PER_YEAR = 365.25

SPEC_17_SITES = ["abdomen","adrenal","bone","bone_marrow","brain_cns","breast",
                 "head_and_neck","liver","lymph_nodes","other","pelvis",
                 "pericardial_and_malignant_pericardial_effusion",
                 "peritoneum_and_malignant_peritoneal_effusion",
                 "pleura_and_malignant_pleural_effusion","pulmonary","skin","thorax"]
ORGAN_DX = ["BRAIN","BONE","LIVER","LUNG","ADRENAL","LYMPH","PLEURA","SUBC_TISSUE","OTHER"]
# SAMPLE_SITE → organ regex
SITE_RGX = {
    "BRAIN":       r"brain|cns|cerebr|menin|epidural",
    "BONE":        r"bone",
    "LIVER":       r"liver|hepat",
    "LUNG":        r"lung|pulmon",
    "ADRENAL":     r"adrenal",
    "LYMPH":       r"lymph",
    "PLEURA":      r"pleur",
    "SUBC_TISSUE": r"skin|subc|soft tissue|connective|chest wall",
}

dict_rows: list[dict] = []
def note(variable, original, mapped, msg):
    dict_rows.append({"variable":variable, "original":original,
                      "mapped":mapped, "note":msg})

# ============================================================
# PULL (cached)
# ============================================================
def cached(name, fn):
    p = RAW / name
    if p.exists():
        print(f"   [cache] {p.name}")
        return pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p, low_memory=False)
    print(f"   [pull]  {name}")
    _ensure_pybioportal()
    df = fn()
    if p.suffix == ".parquet":
        df.to_parquet(p)
    else:
        df.to_csv(p, index=False)
    return df

print(">>> Pulling samples + clinical…")
samp_df = cached("samples.csv",
                 lambda: pd.DataFrame(samples.get_all_samples_in_study(study_id=STUDY)))
print(f"   N samples: {len(samp_df)}")

def pull_clin(level):
    import requests
    r = requests.get(
        f"https://www.cbioportal.org/api/studies/{STUDY}/clinical-data",
        params={"clinicalDataType": level, "projection": "SUMMARY", "pageSize": 1000000},
        timeout=120)
    r.raise_for_status()
    return pd.DataFrame(r.json())

clin_sample = cached("clinical_data_sample.parquet", lambda: pull_clin("SAMPLE"))
clin_patient = cached("clinical_data_patient.parquet", lambda: pull_clin("PATIENT"))
print(f"   clinical SAMPLE rows: {len(clin_sample)};  PATIENT rows: {len(clin_patient)}")
print(f"   unique SAMPLE attrs: {clin_sample['clinicalAttributeId'].nunique()}; "
      f"PATIENT attrs: {clin_patient['clinicalAttributeId'].nunique()}")

print(">>> Pulling mutations…")
mut_df = cached("mutations.parquet",
                lambda: pd.DataFrame(mutations.fetch_muts_in_mol_prof(
                    molecular_profile_id=f"{STUDY}_mutations",
                    sample_list_id=f"{STUDY}_sequenced",
                    projection="DETAILED")))
print(f"   mutation rows: {len(mut_df)}")

print(">>> Pulling discrete CNA…")
cna_df = cached("cna.parquet",
                lambda: pd.DataFrame(cna_mod.fetch_discrete_copy_numbers_in_molecular_profile(
                    molecular_profile_id=f"{STUDY}_cna",
                    sample_list_id=f"{STUDY}_cna",
                    discrete_cn_evt_type="HOMDEL_AND_AMP",
                    projection="DETAILED")))
print(f"   CNA rows: {len(cna_df)}")

# ============================================================
# RESHAPE clinical: long → wide
# ============================================================
print(">>> Reshaping clinical to wide…")
samp_wide = (clin_sample.pivot_table(index=["sampleId","patientId"],
                                      columns="clinicalAttributeId",
                                      values="value", aggfunc="first")
                        .reset_index())
pat_wide  = (clin_patient.pivot_table(index="patientId",
                                       columns="clinicalAttributeId",
                                       values="value", aggfunc="first")
                         .reset_index())
print(f"   samp_wide shape: {samp_wide.shape};  pat_wide shape: {pat_wide.shape}")

clin = samp_wide.merge(pat_wide, on="patientId", how="left", suffixes=("","__pat"))
print(f"   joined clinical shape: {clin.shape}")

# ============================================================
# Determine brain-met status (patient-level, broadcast to samples)
# ============================================================
print(">>> Deriving brain-met flag from SAMPLE_SITE…")
samp_wide["_brain_site"] = samp_wide.get("SAMPLE_SITE", pd.Series([np.nan]*len(samp_wide))).astype(str).str.contains(
    SITE_RGX["BRAIN"], case=False, na=False, regex=True)
patient_brain = samp_wide.groupby("patientId")["_brain_site"].max().rename("any_brain_met_pat").reset_index()
clin = clin.merge(patient_brain, on="patientId", how="left")
clin["any_brain_met"] = clin["any_brain_met_pat"].fillna(False).astype(int)
note("any_brain_met","SAMPLE_SITE matches /brain|cns|cerebr|menin|epidural/i",
     "1 if patient had ANY brain-site sample",
     "Brain biopsy → patient had brain met (sampling time)")
n_brain = int(patient_brain["any_brain_met_pat"].fillna(False).sum())
print(f"   patients with brain-met: {n_brain}   samples in brain-met patients: {int(clin['any_brain_met'].sum())}")

# ============================================================
# Build mutation features (MAF aggregation §12) + germline/silent filter
# ============================================================
print(">>> Building mutation features…")
hg_col = ("gene_hugoGeneSymbol" if "gene_hugoGeneSymbol" in mut_df.columns
          else "hugoGeneSymbol" if "hugoGeneSymbol" in mut_df.columns
          else "Hugo_Symbol")
sid_col = "sampleId"
vc_col = "mutationType" if "mutationType" in mut_df.columns else "Variant_Classification"
ms_col = "mutationStatus" if "mutationStatus" in mut_df.columns else "Mutation_Status"
ta_col = "tumorAltCount" if "tumorAltCount" in mut_df.columns else ("t_alt_count" if "t_alt_count" in mut_df.columns else None)

DROP_VC = {"Silent","silent","3'UTR","5'UTR","3'Flank","5'Flank","Intron","intron","RNA","IGR"}
xf = mut_df.dropna(subset=[hg_col, sid_col]).copy()
xf = xf[~xf[vc_col].isin(DROP_VC)]
print(f"   MAF after filter: {len(xf)} rows; {xf[sid_col].nunique()} unique samples")

agg_kwargs = {
    "mutation_count_all_sites_sum": (hg_col, "size"),
    "genes": (hg_col, lambda s: ";".join(sorted(s.dropna().unique()))),
}
if ta_col is not None:
    agg_kwargs["t_alt_count_max"] = (ta_col, lambda s: pd.to_numeric(s, errors="coerce").max())
mut_feat = xf.groupby(sid_col).agg(**agg_kwargs).reset_index().rename(columns={sid_col:"SAMPLE_ID"})

# Reindex to ALL samples (samples with no muts → 0)
all_samples = clin[["sampleId"]].rename(columns={"sampleId":"SAMPLE_ID"}).drop_duplicates()
mut_feat = all_samples.merge(mut_feat, on="SAMPLE_ID", how="left")
mut_feat["mutation_count_all_sites_sum"] = mut_feat["mutation_count_all_sites_sum"].fillna(0).astype(int)
if "t_alt_count_max" not in mut_feat.columns:
    mut_feat["t_alt_count_max"] = np.nan
mut_feat["t_alt_count_max"] = pd.to_numeric(mut_feat["t_alt_count_max"], errors="coerce")
mut_feat["genes"] = mut_feat["genes"].fillna("")

def safe_qcut(s, q=4):
    s = pd.to_numeric(s, errors="coerce")
    try:
        return pd.qcut(s, q=q, labels=[f"Q{i+1}" for i in range(q)], duplicates="drop")
    except ValueError:
        return pd.Series([np.nan]*len(s), index=s.index)
mut_feat["mutation_count_q"] = safe_qcut(mut_feat["mutation_count_all_sites_sum"])
mut_feat["t_alt_count_q"]    = safe_qcut(mut_feat["t_alt_count_max"])

# ============================================================
# Build gene-binary matrix (mutations + CNA) per gnomeR §13
# ============================================================
print(">>> Building gene binary (mutations + CNA Amp/Del)…")
# Mutation binary
xf_for_bin = xf[xf[ms_col].astype(str).str.lower().eq("somatic") | xf[ms_col].astype(str).str.lower().eq("somatic")]
# More permissive: keep all if mutationStatus is missing or =="SOMATIC"/"Somatic"
xf_for_bin = xf[xf[ms_col].astype(str).str.upper().isin(["SOMATIC","UNKNOWN","NA"]) | xf[ms_col].isna()]
mut_bin = (xf_for_bin.assign(val=1)
                     .pivot_table(index=sid_col, columns=hg_col,
                                  values="val", aggfunc="max", fill_value=0)
                     .reset_index()
                     .rename(columns={sid_col:"SAMPLE_ID"}))
mut_bin.columns.name = None
mut_bin = all_samples.merge(mut_bin, on="SAMPLE_ID", how="left").fillna(0)
gene_mut_cols = [c for c in mut_bin.columns if c != "SAMPLE_ID"]
mut_bin[gene_mut_cols] = mut_bin[gene_mut_cols].astype("int8")
print(f"   mutation binary: {len(gene_mut_cols)} gene cols")

# CNA: long → Amp / Del per gnomeR (high-level only per spec §13)
cna_hg = ("gene_hugoGeneSymbol" if "gene_hugoGeneSymbol" in cna_df.columns
          else "hugoGeneSymbol" if "hugoGeneSymbol" in cna_df.columns
          else "gene")
cna_val = "alteration"
cna_sid = "sampleId"
cna_df["_alt"] = pd.to_numeric(cna_df[cna_val], errors="coerce")
# high-level: 2 (Amp), -2 (Del); spec §13 says high_level_cna_only=TRUE
amp = cna_df[cna_df["_alt"].eq(2)]
del_ = cna_df[cna_df["_alt"].eq(-2)]
amp_w = (amp.assign(col=lambda d: d[cna_hg]+".Amp", val=1)
            .pivot_table(index=cna_sid, columns="col", values="val", aggfunc="max", fill_value=0))
del_w = (del_.assign(col=lambda d: d[cna_hg]+".Del", val=1)
            .pivot_table(index=cna_sid, columns="col", values="val", aggfunc="max", fill_value=0))
cna_bin = amp_w.join(del_w, how="outer").fillna(0).reset_index().rename(columns={cna_sid:"SAMPLE_ID"})
cna_bin.columns.name = None
cna_bin = all_samples.merge(cna_bin, on="SAMPLE_ID", how="left").fillna(0)
gene_cna_cols = [c for c in cna_bin.columns if c != "SAMPLE_ID"]
cna_bin[gene_cna_cols] = cna_bin[gene_cna_cols].astype("int8")
print(f"   CNA binary: {len(gene_cna_cols)} gene/event cols")

# ============================================================
# Apply Sanchez-Vega pathways (10 default) - light Python re-impl
# ============================================================
PATHWAYS = {
    "RTK/RAS":  ["ABL1","ALK","BRAF","CBL","ERBB2","ERBB3","ERBB4","EGFR","FGFR1","FGFR2","FGFR3","FGFR4",
                  "FLT3","HRAS","KIT","KRAS","MAP2K1","MAP2K2","MAPK1","MET","NF1","NRAS","NTRK1","NTRK2","NTRK3",
                  "PDGFRA","PDGFRB","PTPN11","RAF1","RET","ROS1","SOS1"],
    "Nrf2":     ["CUL3","KEAP1","NFE2L2"],
    "PI3K":     ["AKT1","AKT2","AKT3","INPP4B","MTOR","PIK3CA","PIK3CB","PIK3R1","PIK3R2","PPP2R1A","PTEN",
                  "RICTOR","RPTOR","STK11","TSC1","TSC2"],
    "TGFB":     ["ACVR1B","ACVR2A","SMAD2","SMAD3","SMAD4","TGFBR1","TGFBR2"],
    "p53":      ["ATM","CHEK2","MDM2","MDM4","RPS6KA3","TP53","TP53BP1"],
    "Wnt":      ["AMER1","APC","ARID1A","ARID2","AXIN1","AXIN2","CTNNB1","DKK1","DKK2","DKK3","DKK4",
                  "RNF43","SOX17","TCF7L1","TCF7L2","TLE1","TLE2","TLE3","TLE4","ZNRF3"],
    "Myc":      ["MAX","MGA","MLX","MLXIP","MLXIPL","MNT","MXD1","MXD3","MXD4","MXI1","MYC","MYCL","MYCN"],
    "Cell cycle":["CCND1","CCND2","CCND3","CCNE1","CDK4","CDK6","CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C",
                   "E2F1","E2F3","RB1"],
    "Hippo":    ["CRB1","CRB2","CRB3","DCHS1","DCHS2","FAT1","FAT2","FAT3","FAT4","HMCN1","LATS1","LATS2",
                  "LLGL1","LLGL2","MOB1A","MOB1B","NF2","PARD3","PARD6A","PARD6B","PARD6G","SAV1","SCRIB",
                  "STK3","STK4","TAOK1","TAOK2","TAOK3","TEAD4","TJP1","TJP2","TJP3","WWC1","YAP1"],
    "Notch":    ["APH1A","ARRDC1","CIR1","CNTN6","CREBBP","CTBP1","CTBP2","DLL1","DLL3","DLL4","DTX1","DTX2",
                  "DTX3","DTX3L","DTX4","EP300","FBXW7","FHL1","HDAC2","HES1","HES2","HES3","HES4","HES5","HEYL",
                  "JAG1","JAG2","KAT2B","KDM5A","LFNG","MAML1","MAML2","MAML3","MFNG","NCOR1","NCOR2","NOTCH1",
                  "NOTCH2","NOTCH3","NOTCH4","NRARP","NUMB","NUMBL","PSEN1","PSEN2","PSENEN","RBPJ","RBX1","RFNG",
                  "SAP30","SNW1","SPEN","TLE1","TLE2","TLE3","TLE4"],
}

# Build pathway 0/1 per sample
def pathway_hit(sample_row, pathway_genes):
    for g in pathway_genes:
        if g in sample_row.index and sample_row[g] == 1: return 1
        amp = f"{g}.Amp"; deln = f"{g}.Del"
        if amp in sample_row.index and sample_row[amp] == 1: return 1
        if deln in sample_row.index and sample_row[deln] == 1: return 1
    return 0

# Merge mut + cna binaries on SAMPLE_ID before pathway calc
binmat = mut_bin.merge(cna_bin, on="SAMPLE_ID", how="left").fillna(0)
print(f"   binary matrix: {binmat.shape}")
# Compute pathway cols
for pname, pgenes in PATHWAYS.items():
    mut_hits = pd.Series([0]*len(binmat), index=binmat.index, dtype="int8")
    for g in pgenes:
        if g in binmat.columns:
            mut_hits = mut_hits | binmat[g].astype("int8")
        amp = f"{g}.Amp"; deln = f"{g}.Del"
        if amp in binmat.columns: mut_hits = mut_hits | binmat[amp].astype("int8")
        if deln in binmat.columns: mut_hits = mut_hits | binmat[deln].astype("int8")
    binmat[f"pathway_{pname}"] = mut_hits.astype("int8")
print(f"   added 10 pathway cols")

# Drop all-zero gene cols (gnomeR default)
all_zero_drop = [c for c in binmat.columns if c != "SAMPLE_ID" and not c.startswith("pathway_")
                 and binmat[c].sum() == 0]
binmat = binmat.drop(columns=all_zero_drop)
print(f"   dropped {len(all_zero_drop)} all-zero gene columns")
gene_cols_final = [c for c in binmat.columns if c != "SAMPLE_ID"]

# ============================================================
# Build canonical OUT
# ============================================================
print(">>> Building canonical OUT frame…")
out = pd.DataFrame()
out["SAMPLE_ID"]   = clin["sampleId"]
out["record_id"]   = clin["patientId"]
out["ca_seq"]      = np.nan
out["cpt_number"]  = np.nan

def col(name):
    return clin[name] if name in clin.columns else pd.Series([np.nan]*len(clin))

def yn1(v):
    if pd.isna(v): return pd.NA
    s = str(v).strip().lower()
    if s.startswith("y") or s in {"1","true","positive","metastasis","metastasis_yes"}: return 1
    if s.startswith("n") or s in {"0","false","negative"}: return 0
    return pd.NA

def first_char_bin(v):
    if pd.isna(v): return pd.NA
    s = str(v).strip()
    if not s: return pd.NA
    try: return int(s[0])
    except: return pd.NA

# Demographics
out["SEX"]                  = col("SEX")
out["BIRTH_YEAR"]           = np.nan
out["CENTER"]               = "MSK"
out["PRIMARY_RACE"]         = np.nan
out["ETHNICITY"]            = np.nan
out["naaccr_ethnicity_code"]= np.nan
out["age_dx"]               = pd.to_numeric(col("INVASIVE_CARCINOMA_DX_AGE"), errors="coerce")
out["age_dx_num"]           = out["age_dx"]
out["age_last_fu_yrs"]      = pd.to_numeric(col("LAST_CONTACT_DAYS_TO"), errors="coerce") / DAYS_PER_YEAR
out["age_death_yrs"]        = np.nan
out["race_clean"]           = np.nan
out["ethnicity_clean"]      = np.nan

def age_cat(v):
    v = pd.to_numeric(v, errors="coerce")
    if pd.isna(v): return np.nan
    if v < 50: return "<50"
    if v <= 70: return "50-70"
    return ">70"
out["age_cat"] = out["age_dx_num"].map(age_cat)

# Menopausal (bonus)
out["menopausal_status"] = col("MENOPAUSAL_STATUS_AT_DIAGNOSIS")
note("menopausal_status","MENOPAUSAL_STATUS_AT_DIAGNOSIS","direct","Bonus var not in canonical spec")

# Breast subtype
recep = col("OVERALL_RECEPTOR_STATUS_PATIENT").astype(str)
def map_overall(v):
    if pd.isna(v): return np.nan
    s = str(v).strip().replace(" ","").replace("/","").replace(",","").upper()
    if "HR+HER2-" in s or s == "HR+HER2-": return "HR+, HER2-"
    if "HR+HER2+" in s: return "HR+, HER2+"
    if "HR-HER2+" in s: return "HR-, HER2+"
    if "TNBC" in s or "TRIPLENEGATIVE" in s or "HR-HER2-" in s: return "TNBC"
    return np.nan
out["bca_subtype"] = recep.map(map_overall)
def recep_cat(v):
    if pd.isna(v): return np.nan
    s = str(v).replace(" ","").replace(",","/").replace("/","").upper().replace("HER2","HER2")
    if "HR+/HER2-" in str(v) or v == "HR+, HER2-": return "HR+/HER2-"
    if "HR+/HER2+" in str(v) or v == "HR+, HER2+": return "HR+/HER2+"
    if "HR-/HER2+" in str(v) or v == "HR-, HER2+": return "HR-/HER2+"
    if "TNBC" in str(v).upper() or "TRIPLE" in str(v).upper(): return "Triple Negative"
    return np.nan
out["receptor_primary_cat"] = recep.map(recep_cat)

out["ca_bca_er"]            = col("ER_STATUS_PRIMARY")
out["ca_bca_pr"]            = col("PR_STATUS_PRIMARY")
out["ca_bca_her2ihc_val"]   = col("HER2_IHC_VALUE_PRIMARY")
out["ca_bca_her2ihc_intp"]  = col("HER2_IHC_PRIMARY")
out["ca_bca_her_summ"]      = col("OVERALL_HER2_STATUS_PATIENT")
out["ca_bca_oncotypedx"]    = np.nan
out["ca_bca_mgene"]         = np.nan
out["ca_bca_mgeneresult"]   = np.nan
out["ca_bca_herish"]        = col("HER2_FISH_STATUS")

def her2_bin(v):
    if pd.isna(v): return pd.NA
    s = str(v).lower()
    if re.search(r"positive|amplif|equivocal_pos", s): return 1
    if re.search(r"negative|not amplif|not_amplif", s): return 0
    return pd.NA
out["her2_status_bin"] = out["ca_bca_her_summ"].map(her2_bin).astype("Int64")

# Histology / grade
out["naaccr_histology_cd"]   = np.nan
out["ca_histology"]          = col("TUMOR_SAMPLE_HISTOLOGY")
out["ca_hist_brca"]          = col("TUMOR_SAMPLE_HISTOLOGY")
out["ca_hist_adeno_squamous"]= np.nan
out["ca_grade"]              = col("OVERALL_TUMOR_GRADE")
out["primary_nuclear_grade"] = col("PRIMARY_NUCLEAR_GRADE")
note("primary_nuclear_grade","PRIMARY_NUCLEAR_GRADE","direct","Bonus var")

def grade_ord(v):
    if pd.isna(v): return np.nan
    s = str(v).strip()
    if re.search(r"\blow\b|^\s*I\s*$|^G?1$|grade\s*1|well", s, re.I): return "Low"
    if re.search(r"intermediate|moderately|^\s*II\s*$|^G?2$|grade\s*2", s, re.I): return "Intermediate"
    if re.search(r"\bhigh\b|poorly|^\s*III\s*$|^G?3$|^G?4$|grade\s*3|grade\s*4", s, re.I): return "High"
    return np.nan
out["grade_ord"] = out["ca_grade"].map(grade_ord)

# Stage
out["best_ajcc_stage_cd"] = col("STAGE_AT_DIAGNOSIS")
out["stage_dx"]           = col("STAGE_AT_DIAGNOSIS")
def stage_group(v):
    if pd.isna(v): return np.nan
    s = str(v).upper().strip()
    s = re.sub(r"STAGE\s+", "", s)
    if s in {"0","IS","TIS"}: return "Stage I"
    if s in {"I","IA","IB","IC"}: return "Stage I"
    if s in {"II","IIA","IIB","IIC"}: return "Stage II"
    if s in {"III","IIIA","IIIB","IIIC"}: return "Stage III"
    if s in {"IV","IVA","IVB","IVC"}: return "Stage IV"
    return np.nan
out["stage_diag_group"] = out["stage_dx"].map(stage_group)
out["stage_iv_bin"]     = out["stage_diag_group"].eq("Stage IV").astype("Int64")
out["stage_dx_iv"]      = out["stage_iv_bin"].map({1:"Stage IV", 0:"Stage I-III"})
out["ca_path_t_stage"]  = col("T_STAGE")
out["ca_path_n_stage"]  = col("N_STAGE")
out["naaccr_path_m_cd"] = col("M_STAGE")
out["ca_path_group_stage"] = out["stage_diag_group"]
for c in ["ca_clin_t_stage","ca_clin_n_stage","naaccr_clin_stage_cd","naaccr_path_stage_cd",
          "naaccr_seer_sum_stage","ca_tx_pre_path_stage","naaccr_tnm_path_desc",
          "naaccr_path_t_cd","naaccr_path_n_cd",
          "ca_path_tis_det","ca_path_t1_det","ca_path_t2_det","ca_path_t3_det","ca_path_t4_det",
          "ca_path_n0_det","ca_path_n1_det","ca_path_n2_det","ca_path_n3_det"]:
    out[c] = np.nan

# DMETS at dx: we don't have organ-by-organ dx-time fields; treat brain at sampling as proxy for "any time".
# Derive DMETS_DX_* from SAMPLE_SITE (= site of biopsy, which is a met if SAMPLE_TYPE=Metastasis)
sample_site = col("SAMPLE_SITE").astype(str)
sample_type = col("SAMPLE_TYPE").astype(str)
is_met = sample_type.str.lower().str.contains("metastasis", na=False)

for organ, rgx in SITE_RGX.items():
    hit = (sample_site.str.contains(rgx, case=False, na=False, regex=True) & is_met).astype("Int64")
    out[f"DMETS_DX_{organ}"] = hit
# OTHER: met-site but none of specific organs
any_specific = np.zeros(len(out), dtype=bool)
for o in [x for x in SITE_RGX]:
    any_specific |= (out[f"DMETS_DX_{o}"].fillna(0).astype(int).to_numpy() == 1)
out["DMETS_DX_OTHER"] = ((is_met.to_numpy() & ~any_specific).astype(int))
note("DMETS_DX_*","SAMPLE_SITE & SAMPLE_TYPE","biopsy site of metastatic sample",
     "Razavi cohort: brain biopsy → patient had brain met")

out["ca_dmets_yn"] = col("METASTATIC_DZ_FUP").map(yn1).astype("Int64")
out["CA_DMETS_YN"] = out["ca_dmets_yn"]

for i in range(1, 11):
    out[f"ca_first_dmets{i}"] = np.nan

# Mets during study: only brain_cns populated; use patient-level any_brain_met
# Note: METASTATIC_RECURRENCE_TIME_MONTHS is time-to-first-met (any site)
recur_mos = pd.to_numeric(col("METASTATIC_RECURRENCE_TIME_MONTHS"), errors="coerce")
note("dist_mets_brain_cns","SAMPLE_SITE matches /brain.../ (any of patient's samples)",
     "patient-level broadcast","Razavi cohort biopsy-site only")
for site in SPEC_17_SITES:
    out[f"dist_mets_{site}"]            = np.nan
    out[f"dx_to_dist_mets_{site}_days"] = np.nan
    out[f"dx_to_dist_mets_{site}_mos"]  = np.nan
    out[f"dx_to_dist_mets_{site}_yrs"]  = np.nan
out["dist_mets_brain_cns"]              = clin["any_brain_met"].astype("Int64")
out["dx_to_dist_mets_brain_cns_mos"]    = recur_mos.where(clin["any_brain_met"] == 1, np.nan)
out["dx_to_dist_mets_brain_cns_days"]   = out["dx_to_dist_mets_brain_cns_mos"] * DAYS_PER_MONTH
out["dx_to_dist_mets_brain_cns_yrs"]    = out["dx_to_dist_mets_brain_cns_mos"] / 12.0

out["dmets_post_dx"]    = out["ca_dmets_yn"]
out["dx_to_dmets_days"] = recur_mos * DAYS_PER_MONTH
out["dx_to_dmets_mos"]  = recur_mos
out["dx_to_dmets_yrs"]  = recur_mos / 12.0

# Brain-met cohort
brain_at_dx = ((out["DMETS_DX_BRAIN"].fillna(0).astype(int) == 1) &
                (recur_mos.fillna(99999) <= 1)).astype(int)   # met within 1 month of dx
out["brain_met_at_dx"]  = brain_at_dx
out["any_brain_met"]    = (clin["any_brain_met"].fillna(0).astype(int) == 1).astype(int)
out["brain_met_event"]  = out["any_brain_met"]

def met_loc(row):
    if row["any_brain_met"] == 1: return "Brain"
    for o in ["BONE","LIVER","LUNG","ADRENAL","LYMPH","PLEURA","SUBC_TISSUE","OTHER"]:
        v = row.get(f"DMETS_DX_{o}")
        if pd.notna(v) and v == 1: return "Other"
    if pd.notna(row.get("CA_DMETS_YN")) and row.get("CA_DMETS_YN") == 1: return "Other"
    return "None"
out["met_loc"] = out.apply(met_loc, axis=1)

# Survival
os_status = col("OS_STATUS")
os_months = pd.to_numeric(col("OS_MONTHS"), errors="coerce")
out["OS_STATUS"]   = os_status
out["OS_MONTHS"]   = os_months
out["os_dx_status"] = os_status
out["tt_os_dx_days"] = os_months * DAYS_PER_MONTH
out["tt_os_dx_mos"]  = os_months
out["tt_os_dx_yrs"]  = os_months / 12.0
out["OS_months"]     = os_months
out["os_status_bin"] = os_status.map(first_char_bin).astype("Int64")
out["os_status_f"]   = out["os_status_bin"].map({0:"Alive", 1:"Deceased"})

out["os_adv_status"] = np.nan
out["tt_os_adv_days"] = np.nan
out["tt_os_adv_mos"]  = np.nan
out["tt_os_adv_yrs"]  = np.nan

# PFS not separated in MSK 2018; map nothing → NA
out["PFS_I_ADV_STATUS"]   = np.nan
out["PFS_I_ADV_MONTHS"]   = np.nan
out["pfs_i_adv_status"]   = np.nan
out["tt_pfs_i_adv_days"]  = np.nan
out["tt_pfs_i_adv_mos"]   = np.nan
out["tt_pfs_i_adv_yrs"]   = np.nan
out["PFS_imaging_months"] = np.nan
out["pfs_i_event_bin"]    = pd.NA
out["PFS_M_ADV_STATUS"]   = np.nan
out["PFS_M_ADV_MONTHS"]   = np.nan
out["pfs_m_adv_status"]   = np.nan
out["tt_pfs_m_adv_days"]  = np.nan
out["tt_pfs_m_adv_mos"]   = np.nan
out["tt_pfs_m_adv_yrs"]   = np.nan
out["PFS_medonc_months"]  = np.nan
out["pfs_m_event_bin"]    = pd.NA
out["pfs_cohort"]         = np.nan
out["pfs_i_or_m_adv_status"]  = np.nan
out["tt_pfs_i_or_m_adv_days"] = np.nan
out["tt_pfs_i_or_m_adv_mos"]  = np.nan
out["tt_pfs_i_or_m_adv_yrs"]  = np.nan
out["pfs_i_and_m_adv_status"] = np.nan
out["tt_pfs_i_and_m_adv_days"] = np.nan
out["tt_pfs_i_and_m_adv_mos"]  = np.nan
out["tt_pfs_i_and_m_adv_yrs"]  = np.nan

# DFS (from cBio)
dfs_event = col("DFS_EVENT")
dfs_mos   = pd.to_numeric(col("DFS_MONTHS"), errors="coerce")
out["dfs_status"] = dfs_event
out["dfs_months"] = dfs_mos
out["dss_status"] = np.nan
out["dss_months"] = np.nan

# Time to brain met
out["tt_brain_met_mos"] = np.where(
    out["any_brain_met"] == 1,
    recur_mos,                      # use time to first met as best proxy
    os_months,
)
out["time_to_brain_met_mos"] = out["tt_brain_met_mos"]

# Sample / panel metadata
out["ONCOTREE_CODE"]         = col("ONCOTREE_CODE")
out["SAMPLE_TYPE_DETAILED"]  = col("SAMPLE_TYPE")
out["sample_type"]           = col("SAMPLE_TYPE")
out["SEQ_ASSAY_ID"]          = "MSK-IMPACT"
out["AGE_AT_SEQUENCING"]     = pd.to_numeric(col("INVASIVE_CARCINOMA_DX_AGE"), errors="coerce")
out["PDL1_POSITIVE_ANY"]     = np.nan
out["PDL1_TESTING"]          = np.nan
out["CPT_SEQ_DATE"]          = col("NGS_SAMPLE_COLLECTION_TIME")
out["SAMPLE_SITE"]           = col("SAMPLE_SITE")
note("SAMPLE_SITE","SAMPLE_SITE","direct","Razavi-specific biopsy site")
for nacol in ["cohort","institution","release_version","cpt_n_ca_seq","cpt_order_int",
              "cpt_seq_date","dob_cpt_report_days","dob_cpt_report_mos","dob_cpt_report_yrs",
              "cpt_report_post_death","cpt_report_post_last_alive","dx_cpt_rep_days",
              "dx_cpt_rep_mos","dx_cpt_rep_yrs","dx_path_proc_cpt_days","dx_path_proc_cpt_mos",
              "dx_path_proc_cpt_yrs","path_proc_cpt_rep_days","path_proc_cpt_rep_mos",
              "path_proc_cpt_rep_yrs","path_proc_number","path_rep_number","cpt_oncotree_code",
              "cpt_seq_assay_id","mutations","cna","cohort__master","institution__master",
              "release_version__master","tr_eligible","redcap_ca_index","dob_ca_dx_days",
              "dob_ca_dx_mos","dob_ca_dx_yrs","ca_dx_how","dob_next_ca_days","dob_next_ca_mos",
              "dob_next_ca_yrs","first_index_ca_days","first_index_ca_mos","first_index_ca_yrs",
              "naaccr_first_contact_int","ca_d_site","ca_type","naaccr_behavior_cd",
              "naaccr_laterality_cd"]:
    out[nacol] = np.nan

def sample_type_bin(v):
    if pd.isna(v): return pd.NA
    s = str(v).lower()
    if "primary" in s: return 1
    if "metast" in s: return 2
    return pd.NA
out["sample_type_bin"] = out["SAMPLE_TYPE_DETAILED"].map(sample_type_bin).astype("Int64")

# Non-index cancer
out["non_idx_n_cancers"] = 0
out["non_idx_any_heme"]  = pd.NA
out["non_idx_any_brain"] = pd.NA
out["non_idx_types"]     = ""
out["had_prior_non_breast_cancer"] = col("PRIOR_BREAST_PRIMARY").map(
    lambda v: pd.NA if pd.isna(v) else (1 if str(v).strip().lower().startswith("y") else 0)
).astype("Int64")

# Chemo/radiation - MSK 2018 cohort doesn't have treatment fields; all NA
out["hx_chemo_neoadj_bin"]   = pd.NA
out["hx_chemo_any_bin"]      = pd.NA
out["hx_chemo_adjuvant_bin"] = pd.NA
out["hx_radiation_bin"]      = pd.NA
note("hx_chemo_*","not in breast_msk_2018","NA",
     "Razavi cohort doesn't include treatment fields in cBioPortal export")

# Regimens parity
out["N_REGIMENS_PT"] = np.nan
out["n_regimens_pt"] = np.nan
out["ca_n_regimens"] = np.nan
out["n_cpt_pt"]      = pd.to_numeric(col("SAMPLE_COUNT"), errors="coerce")

# Bonus: prior local recurrence
out["PRIOR_LOCAL_RECURRENCE"] = col("PRIOR_LOCAL_RECURRENCE")

# Genomic features attach
out = out.merge(mut_feat, on="SAMPLE_ID", how="left")
out = out.merge(binmat, on="SAMPLE_ID", how="left")
print(f"   OUT post-genomic attach shape: {out.shape}")

# ============================================================
# Top-gene pipeline (§15) - cohort-specific
# ============================================================
print(">>> Computing top-gene pipeline…")
mut_only = [c for c in gene_cols_final
            if not (c.endswith(".Amp") or c.endswith(".Del") or c.endswith(".fus")
                    or c.startswith("pathway_"))]
brain = out[out["any_brain_met"] == 1]
total_prev = out[mut_only].sum(axis=0).rename("n_total_samples_mutated")
brain_prev = brain[mut_only].sum(axis=0).rename("n_brain_met_samples_mutated")
prev_df = pd.concat([brain_prev, total_prev], axis=1).reset_index().rename(columns={"index":"gene"})
prev_df["pct_brain_met"] = prev_df["n_brain_met_samples_mutated"] / max(len(brain),1)
prev_df["pct_total"]     = prev_df["n_total_samples_mutated"] / max(len(out),1)
prev_df = prev_df.sort_values("n_brain_met_samples_mutated", ascending=False)
prev_df.to_csv(OUT_GENEPREV, index=False)

TOP10 = prev_df["gene"].head(10).tolist()
TOP5  = TOP10[:5]
print(f"   TOP10 (brain-met samples, N={len(brain)}): {TOP10}")
with open(OUT_TOPGENES, "w") as fh:
    fh.write("\n".join(TOP10) + "\n")

for g in TOP10:
    out[f"G_top10_{g}"] = out[g].fillna(0).astype(int)
out["top5_n_mutated"]   = out[TOP5].fillna(0).sum(axis=1).astype(int)
out["top10_n_mutated"]  = out[TOP10].fillna(0).sum(axis=1).astype(int)
out["top5_any_mutated"] = (out["top5_n_mutated"] > 0).astype(int)
out["top10_any_mutated"]= (out["top10_n_mutated"] > 0).astype(int)

# Save
out.to_csv(OUT_MAIN, index=False)
pd.DataFrame(dict_rows).to_csv(OUT_DICT, index=False)

# ============================================================
# Validation
# ============================================================
print("\n========== VALIDATION ==========")
print(f"TOTAL ROWS: {len(out)}")
print(f"TOTAL COLS: {out.shape[1]}")
print(f"\nCohort split:")
print(f"  any_brain_met == 1:   {(out['any_brain_met']==1).sum()} samples "
      f"({out.loc[out['any_brain_met']==1,'record_id'].nunique()} unique patients)")
print(f"  brain_met_at_dx == 1: {(out['brain_met_at_dx']==1).sum()} samples")
print(f"  no brain met:         {(out['any_brain_met']==0).sum()} samples")
print(f"\nReceptor subtype x brain met:")
print(pd.crosstab(out["receptor_primary_cat"], out["any_brain_met"], dropna=False))
print(f"\nTop 10 mutated genes in any_brain_met cohort:")
tp = prev_df.head(10)[["gene","pct_brain_met","pct_total"]].copy()
tp["pct_brain_met"] = (tp["pct_brain_met"]*100).round(1).astype(str) + "%"
tp["pct_total"]     = (tp["pct_total"]*100).round(1).astype(str) + "%"
print(tp.to_string(index=False))
print(f"\nWrote {OUT_MAIN.name}")
print(f"Wrote {OUT_TOPGENES.name}")
print(f"Wrote {OUT_GENEPREV.name}")
print(f"Wrote {OUT_DICT.name}")
