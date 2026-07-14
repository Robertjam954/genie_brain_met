"""
Build extracted_variables_tcga_data.csv per harmonization_spec.md.

Column names mirror extracted_variables_genie_data.csv exactly (canonical
schema), so the two files can be row-stacked for cross-cohort analysis.

Sources:
  M = complete_tcga_data_merged_with_pathways.csv  (clinical + gene binary + pathways)
  C = brca_tcga_pan_can_atlas_2018_cna_additional_clinical.csv (GDC clinical export)
  X = brca_tcga_pan_can_atlas_2018_mutations.csv  (cBio API MAF)

Outputs (data/processed/):
  extracted_variables_tcga_data.csv
  extracted_variables_tcga_top_genes.txt
  extracted_variables_tcga_gene_prev_brain_met.csv
  extracted_variables_tcga_dictionary.csv
"""

from __future__ import annotations
import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/Users/robertjames/Documents/Documents - Robert’s iMac/"
            "Research Projects/MSKCC Research Fellowship/Projects/"
            "genie_tcga_impact_brain_mets")
PROC = ROOT / "data/processed"
OUT_DIR = ROOT / "src" / "exploratory data analysis"
OUT_DIR.mkdir(parents=True, exist_ok=True)
M_PATH = PROC / "complete_tcga_data_merged_with_pathways.csv"
C_PATH = PROC / "brca_tcga_pan_can_atlas_2018_cna_additional_clinical.csv"
X_PATH = PROC / "brca_tcga_pan_can_atlas_2018_mutations.csv"
# Brain mets are ONLY in the GDC follow_up.tsv for TCGA-BRCA (n=6).
F_PATH = ROOT / "data/external/clinical.project-tcga-brca.2026-06-08/follow_up.tsv"
OUT_MAIN = OUT_DIR / "extracted_variables_tcga_data.csv"
OUT_TOPGENES = OUT_DIR / "extracted_variables_tcga_top_genes.txt"
OUT_GENEPREV = OUT_DIR / "extracted_variables_tcga_gene_prev_brain_met.csv"
OUT_DICT = OUT_DIR / "extracted_variables_tcga_dictionary.csv"

DAYS_PER_MONTH = 30.4375
DAYS_PER_YEAR = 365.25

ORGAN_DX = ["BRAIN", "BONE", "LIVER", "LUNG", "ADRENAL",
             "LYMPH", "PLEURA", "SUBC_TISSUE", "OTHER"]
SITE_PATTERNS = {
    "BRAIN":       r"brain|cns|cerebr|menin",
    "BONE":        r"bone",
    "LIVER":       r"liver|hepat",
    "LUNG":        r"lung|pulmon",
    "ADRENAL":     r"adrenal",
    "LYMPH":       r"lymph",
    "PLEURA":      r"pleur",
    "SUBC_TISSUE": r"skin|subc|soft tissue|connective",
}
SPEC_17_SITES = ["abdomen","adrenal","bone","bone_marrow","brain_cns","breast",
                 "head_and_neck","liver","lymph_nodes","other","pelvis",
                 "pericardial_and_malignant_pericardial_effusion",
                 "peritoneum_and_malignant_peritoneal_effusion",
                 "pleura_and_malignant_pleural_effusion","pulmonary","skin","thorax"]

dict_rows: list[dict] = []
def note(variable, original, mapped, msg):
    dict_rows.append({"variable": variable, "original": original,
                      "mapped": mapped, "note": msg})

# ============================================================
# Load
# ============================================================
print(">>> Loading M (merged TCGA)…")
m = pd.read_csv(M_PATH, low_memory=False)
print(f"   M shape: {m.shape}")
print(">>> Loading C (GDC clinical)…")
c = pd.read_csv(C_PATH, low_memory=False, na_values=["'--", "--", "Not Reported"])
print(f"   C shape: {c.shape}")
print(">>> Loading X (mutations MAF)…")
x = pd.read_csv(X_PATH, low_memory=False)
print(f"   X shape: {x.shape}")

# ============================================================
# Dedupe C to per-patient + roll up chemo
# ============================================================
c["patient_barcode"] = c["cases.submitter_id"].astype(str).str.strip()
treat_cols = [col for col in c.columns if col.startswith("treatments.")]
print(f"   {len(treat_cols)} treatments.* cols in C")

def first_non_null(s):
    s = s.dropna()
    return s.iloc[0] if len(s) else np.nan

non_treat = [col for col in c.columns if not col.startswith("treatments.")]
c_pat = c.groupby("patient_barcode")[non_treat].agg(first_non_null).reset_index(drop=True)

ta = c[["patient_barcode"] + treat_cols].copy()
def gstr(name):
    return (ta[name].fillna("").astype(str).str.lower()
            if name in ta.columns else pd.Series([""]*len(ta), index=ta.index))
ttype  = gstr("treatments.treatment_type")
tor    = gstr("treatments.treatment_or_therapy")
intent = gstr("treatments.treatment_intent_type")

chemo_pat = re.compile(r"chemo|pharmaceutical|cytotox", re.I)
is_chemo = (ttype.str.contains(chemo_pat, na=False) |
            (tor.str.contains("yes", na=False) & ttype.str.contains(chemo_pat, na=False)))
is_neo = intent.str.contains("neoadj", na=False) & is_chemo
is_adj = intent.str.contains("adjuv", na=False) & is_chemo & ~is_neo

ta["_chemo"] = is_chemo.astype(int)
ta["_neo"]   = is_neo.astype(int)
ta["_adj"]   = is_adj.astype(int)

tx_roll = ta.groupby("patient_barcode").agg(
    hx_chemo_any_bin     = ("_chemo", "max"),
    hx_chemo_neoadj_C    = ("_neo", "max"),
    hx_chemo_adjuvant_bin= ("_adj", "max"),
).reset_index()
for k in ["hx_chemo_any_bin","hx_chemo_neoadj_C","hx_chemo_adjuvant_bin"]:
    tx_roll[k] = tx_roll[k].fillna(0).astype(int)

c_pat = c_pat.merge(tx_roll, on="patient_barcode", how="left")
print(f"   C deduped: {len(c_pat)} patients; chemo_any>=1: {int(c_pat['hx_chemo_any_bin'].fillna(0).sum())}")

# ============================================================
# Pull brain-mets recurrence from GDC follow_up.tsv
# (the 3 source files don't have it; this is required for any_brain_met)
# ============================================================
if F_PATH.exists():
    print(">>> Loading F (GDC follow_up) for brain-mets recurrence…")
    f = pd.read_csv(F_PATH, sep="\t", low_memory=False,
                    na_values=["'--", "--", "Not Reported"])
    print(f"   F shape: {f.shape}")
    f["patient_barcode"] = f["cases.submitter_id"].astype(str).str.strip()
    site_col = "follow_ups.progression_or_recurrence_anatomic_site"
    days_col = "follow_ups.days_to_progression"
    days_col2 = "follow_ups.days_to_recurrence"
    days_first = "follow_ups.days_to_first_event"

    is_brain = f[site_col].astype(str).str.contains(r"brain|cns|cerebr|menin",
                                                     case=False, na=False, regex=True)
    fb = f[is_brain].copy()
    print(f"   brain recurrences in follow_up: {len(fb)} rows / {fb['patient_barcode'].nunique()} patients")

    days_progression = pd.to_numeric(fb[days_col], errors="coerce")
    days_recurrence  = pd.to_numeric(fb[days_col2], errors="coerce") if days_col2 in fb.columns else pd.Series(dtype=float)
    days_first_event = pd.to_numeric(fb[days_first], errors="coerce") if days_first in fb.columns else pd.Series(dtype=float)
    fb["_min_days"] = pd.concat([days_progression, days_recurrence, days_first_event], axis=1).min(axis=1)

    brain_fu = fb.groupby("patient_barcode").agg(
        dist_mets_brain_cns_FU = ("_min_days", lambda s: 1),
        dx_to_dist_mets_brain_cns_days_FU = ("_min_days", "min"),
    ).reset_index()
else:
    print("   (follow_up.tsv not found - dist_mets_brain_cns will be NA)")
    brain_fu = pd.DataFrame(columns=["patient_barcode","dist_mets_brain_cns_FU",
                                     "dx_to_dist_mets_brain_cns_days_FU"])

# ============================================================
# Join M ↔ C ↔ follow_up brain mets
# ============================================================
m["patient_barcode"] = m["PATIENT_ID"].astype(str).str.strip()
overlap = set(m["patient_barcode"]) & set(c_pat["patient_barcode"])
print(f"   M∩C overlap: {len(overlap)}/{m['patient_barcode'].nunique()}")
mc = m.merge(c_pat, on="patient_barcode", how="left", suffixes=("", "__c"))
mc = mc.merge(brain_fu, on="patient_barcode", how="left")
n_brain_join = int(mc["dist_mets_brain_cns_FU"].fillna(0).sum())
print(f"   joined MC shape: {mc.shape};  brain-met patients linked: {n_brain_join}")

# ============================================================
# Mutation features from MAF (§12)
# ============================================================
print(">>> Building mutation features from MAF…")
hg = "hugoGeneSymbol" if "hugoGeneSymbol" in x.columns else "Hugo_Symbol"
sid = "sampleId" if "sampleId" in x.columns else "Tumor_Sample_Barcode"
vc  = "Variant_Classification" if "Variant_Classification" in x.columns else "mutationType"
ta_col = "t_alt_count" if "t_alt_count" in x.columns else None
gaf = "gnomAD_AF" if "gnomAD_AF" in x.columns else None

DROP_VC = {"Silent","3'UTR","5'UTR","3'Flank","5'Flank","Intron","RNA","IGR"}
xf = x.dropna(subset=[hg, sid]).copy()
xf = xf[~xf[vc].isin(DROP_VC)]
if gaf is not None:
    pre = len(xf)
    af_num = pd.to_numeric(xf[gaf], errors="coerce")
    xf = xf[~(af_num.fillna(0) >= 0.001)]
    note("MAF_filter", gaf, "<0.001", f"germline filter dropped {pre - len(xf)} rows")
print(f"   MAF after filter: {len(xf)} rows, {xf[sid].nunique()} unique sampleIds")

agg_kwargs = {
    "mutation_count_all_sites_sum": (hg, "size"),
    "genes": (hg, lambda s: ";".join(sorted(s.dropna().unique()))),
}
if ta_col is not None:
    agg_kwargs["t_alt_count_max"] = (ta_col, lambda s: pd.to_numeric(s, errors="coerce").max())

mut_feat = xf.groupby(sid).agg(**agg_kwargs).reset_index().rename(columns={sid: "SAMPLE_ID"})
all_samples = mc[["SAMPLE_ID"]].drop_duplicates()
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
# Helpers
# ============================================================
def yn1(v):
    if pd.isna(v): return pd.NA
    s = str(v).strip().lower()
    if s.startswith("y") or s in {"1","true","positive"}: return 1
    if s.startswith("n") or s in {"0","false","negative"}: return 0
    return pd.NA

def first_char_bin(v):
    if pd.isna(v): return pd.NA
    s = str(v).strip()
    if not s: return pd.NA
    try: return int(s[0])
    except: return pd.NA

def get(col_name, default=np.nan):
    return mc[col_name] if col_name in mc.columns else pd.Series([default]*len(mc), index=mc.index)

# ============================================================
# Build canonical OUT frame
# ============================================================
print(">>> Building canonical OUT frame…")
out = pd.DataFrame(index=mc.index)

# ---- IDs ----
out["SAMPLE_ID"]   = mc["SAMPLE_ID"]
out["record_id"]   = mc["PATIENT_ID"]
out["ca_seq"]      = np.nan; note("ca_seq", "", "NA", "TCGA: 1 cancer/patient")
out["cpt_number"]  = np.nan; note("cpt_number", "", "NA", "TCGA not panel-based")

# ---- Demographics ----
out["SEX"]          = get("SEX")
out["BIRTH_YEAR"]   = pd.to_numeric(get("demographic.year_of_birth"), errors="coerce")
out["CENTER"]       = get("TISSUE_SOURCE_SITE")
out["PRIMARY_RACE"] = get("RACE")
out["ETHNICITY"]    = get("ETHNICITY")
out["naaccr_ethnicity_code"] = np.nan
out["age_dx"]       = pd.to_numeric(get("AGE"), errors="coerce")
out["age_dx_num"]   = out["age_dx"]
out["age_last_fu_yrs"] = pd.to_numeric(get("DAYS_LAST_FOLLOWUP"), errors="coerce") / DAYS_PER_YEAR
out["age_death_yrs"] = np.nan

def race_clean(v):
    if pd.isna(v): return np.nan
    s = str(v)
    if re.search(r"white|middle east", s, re.I): return "White"
    if re.search(r"black|african", s, re.I): return "Black"
    if re.search(r"asian|indian|chinese", s, re.I): return "Asian"
    if re.search(r"native|alaska", s, re.I): return "Native American"
    return np.nan
def eth_clean(v):
    if pd.isna(v): return np.nan
    s = str(v)
    if re.search(r"non-spanish|non-hispanic|not hispanic", s, re.I): return "Non-Hispanic"
    if re.search(r"hispanic|latino|cuban|mexican|puerto", s, re.I): return "Hispanic"
    return np.nan
out["race_clean"]      = out["PRIMARY_RACE"].map(race_clean)
out["ethnicity_clean"] = out["ETHNICITY"].map(eth_clean)

def age_cat(v):
    v = pd.to_numeric(v, errors="coerce")
    if pd.isna(v): return np.nan
    if v < 50: return "<50"
    if v <= 70: return "50-70"
    return ">70"
out["age_cat"] = out["age_dx_num"].map(age_cat)

# ---- Breast subtype / receptor ----
def pam50_to_bca(s):
    if pd.isna(s): return np.nan
    sl = str(s).strip().lower()
    if "luma" in sl: return "HR+, HER2-"
    if "lumb" in sl: return "HR+, HER2-"
    if "her2" in sl: return "HR-, HER2+"
    if "basal" in sl: return "TNBC"
    if "normal" in sl: return np.nan
    return np.nan
out["bca_subtype"] = get("SUBTYPE").map(pam50_to_bca)
note("bca_subtype","SUBTYPE (PAM50)","HR+/-, HER2+/-",
     "LumA/LumB→HR+/HER2-, Her2→HR-/HER2+, Basal→TNBC; Normal→NA")

def recep_cat(v):
    if pd.isna(v): return np.nan
    s = str(v).replace(" ","").replace(",","/")
    if s == "HR-/HER2-": return "Triple Negative"
    return {"HR+/HER2-":"HR+/HER2-","HR+/HER2+":"HR+/HER2+","HR-/HER2+":"HR-/HER2+",
            "TNBC":"Triple Negative"}.get(s, np.nan)
out["receptor_primary_cat"] = out["bca_subtype"].map(recep_cat)

for col in ["ca_bca_er","ca_bca_pr","ca_bca_her2ihc_val","ca_bca_her2ihc_intp",
            "ca_bca_her_summ","ca_bca_oncotypedx","ca_bca_mgene","ca_bca_mgeneresult",
            "ca_bca_herish"]:
    out[col] = np.nan
note("ca_bca_*","NA","NA","No IHC fields in TCGA; SUBTYPE proxy via bca_subtype")
out["her2_status_bin"] = pd.NA
note("her2_status_bin","SUBTYPE","NA","PAM50 LumB mixed - can't reliably derive HER2 IHC")

# ---- Histology / grade ----
out["naaccr_histology_cd"]   = np.nan
out["ca_histology"]          = get("ICD_O_3_HISTOLOGY")
out["ca_hist_brca"]          = np.nan
out["ca_hist_adeno_squamous"]= np.nan
out["ca_grade"]              = get("GRADE")

def grade_ord(v):
    if pd.isna(v): return np.nan
    s = str(v).strip()
    if re.search(r"\blow\b|^\s*I\s*$|^G?1$|grade\s*1", s, re.I): return "Low"
    if re.search(r"intermediate|^\s*II\s*$|^G?2$|grade\s*2", s, re.I): return "Intermediate"
    if re.search(r"\bhigh\b|^\s*III\s*$|^G?3$|^G?4$|grade\s*3|grade\s*4", s, re.I): return "High"
    return np.nan
out["grade_ord"] = out["ca_grade"].map(grade_ord)

# ---- Stage ----
out["best_ajcc_stage_cd"] = get("AJCC_PATHOLOGIC_TUMOR_STAGE")
out["stage_dx"]           = get("AJCC_PATHOLOGIC_TUMOR_STAGE")
def stage_group(v):
    if pd.isna(v): return np.nan
    s = str(v).upper().strip()
    s = re.sub(r"STAGE\s+", "", s)
    if s in {"I","IA","IB","IC"}: return "Stage I"
    if s in {"II","IIA","IIB","IIC"}: return "Stage II"
    if s in {"III","IIIA","IIIB","IIIC"}: return "Stage III"
    if s in {"IV","IVA","IVB","IVC"}: return "Stage IV"
    return np.nan
out["stage_diag_group"] = out["stage_dx"].map(stage_group)
out["stage_iv_bin"]     = out["stage_diag_group"].eq("Stage IV").astype("Int64")
out["stage_dx_iv"]      = out["stage_iv_bin"].map({1:"Stage IV", 0:"Stage I-III"})
out["ca_path_t_stage"]  = get("PATH_T_STAGE")
out["ca_path_n_stage"]  = get("PATH_N_STAGE")
out["naaccr_path_m_cd"] = get("PATH_M_STAGE")
out["ca_path_group_stage"] = out["stage_diag_group"]
for col in ["ca_clin_t_stage","ca_clin_n_stage","naaccr_clin_stage_cd","naaccr_path_stage_cd",
            "naaccr_seer_sum_stage","ca_tx_pre_path_stage","naaccr_tnm_path_desc",
            "naaccr_path_t_cd","naaccr_path_n_cd",
            "ca_path_tis_det","ca_path_t1_det","ca_path_t2_det","ca_path_t3_det","ca_path_t4_det",
            "ca_path_n0_det","ca_path_n1_det","ca_path_n2_det","ca_path_n3_det"]:
    out[col] = np.nan

# ---- Mets at dx (DMETS_DX_*) ----
soi = get("diagnoses.sites_of_involvement").fillna("").astype(str)
mad = get("diagnoses.metastasis_at_diagnosis")
has_text = (soi.str.len() > 0).to_numpy()
text_idx = np.where(has_text)[0]

for organ in [o for o in ORGAN_DX if o != "OTHER"]:
    pat = SITE_PATTERNS[organ]
    hit = soi.str.contains(pat, case=False, regex=True, na=False).to_numpy()
    s = pd.Series(pd.NA, index=out.index, dtype="Int64")
    s.iloc[text_idx] = hit[text_idx].astype(int)
    out[f"DMETS_DX_{organ}"] = s

any_specific = np.zeros(len(out), dtype=bool)
for o in [x for x in ORGAN_DX if x != "OTHER"]:
    any_specific |= (out[f"DMETS_DX_{o}"].fillna(0).astype(int).to_numpy() == 1)
s = pd.Series(pd.NA, index=out.index, dtype="Int64")
s.iloc[text_idx] = (~any_specific[text_idx]).astype(int)
out["DMETS_DX_OTHER"] = s

out["ca_dmets_yn"] = mad.map(yn1).astype("Int64")
out["CA_DMETS_YN"] = out["ca_dmets_yn"]

for i in range(1, 11):
    out[f"ca_first_dmets{i}"] = np.nan

# ---- Mets during study (17 sites, days/mos/yrs) ----
# Only brain_cns is populated for TCGA (from follow_up.tsv); other organs NA.
note("dist_mets_*", "follow_up.tsv (brain_cns only)", "brain_cns derived",
     "TCGA other organs NA; brain_cns from GDC follow_up (n=6)")
for site in SPEC_17_SITES:
    out[f"dist_mets_{site}"]                = np.nan
    out[f"dx_to_dist_mets_{site}_days"]     = np.nan
    out[f"dx_to_dist_mets_{site}_mos"]      = np.nan
    out[f"dx_to_dist_mets_{site}_yrs"]      = np.nan

# Populate brain_cns from the follow-up join
fu_flag = pd.to_numeric(mc.get("dist_mets_brain_cns_FU"), errors="coerce")
fu_days = pd.to_numeric(mc.get("dx_to_dist_mets_brain_cns_days_FU"), errors="coerce")
out["dist_mets_brain_cns"] = fu_flag.fillna(0).astype("Int64")
out["dx_to_dist_mets_brain_cns_days"] = fu_days
out["dx_to_dist_mets_brain_cns_mos"]  = fu_days / DAYS_PER_MONTH
out["dx_to_dist_mets_brain_cns_yrs"]  = fu_days / DAYS_PER_YEAR

por = get("diagnoses.progression_or_recurrence")
out["dmets_post_dx"]   = por.map(yn1).astype("Int64")
dtr = pd.to_numeric(get("diagnoses.days_to_recurrence"), errors="coerce")
out["dx_to_dmets_days"] = dtr
out["dx_to_dmets_mos"]  = dtr / DAYS_PER_MONTH
out["dx_to_dmets_yrs"]  = dtr / DAYS_PER_YEAR

# ---- Brain-met cohort (derived: at-dx OR during-study) ----
brain_at_dx = (out["DMETS_DX_BRAIN"].fillna(0).astype(int) == 1).astype(int)
brain_fu_int = (out["dist_mets_brain_cns"].fillna(0).astype(int) == 1).astype(int)
out["brain_met_at_dx"] = brain_at_dx
out["any_brain_met"]   = ((brain_at_dx == 1) | (brain_fu_int == 1)).astype(int)
note("any_brain_met","DMETS_DX_BRAIN | dist_mets_brain_cns","at-dx OR during-study",
     "dist_mets_brain_cns NA in TCGA-only sources")

def met_loc(row):
    if row["any_brain_met"] == 1: return "Brain"
    for o in ["BONE","LIVER","LUNG","ADRENAL","LYMPH","PLEURA","SUBC_TISSUE","OTHER"]:
        v = row.get(f"DMETS_DX_{o}")
        if pd.notna(v) and v == 1: return "Other"
    if pd.notna(row.get("CA_DMETS_YN")) and row.get("CA_DMETS_YN") == 1: return "Other"
    return "None"
out["met_loc"]         = out.apply(met_loc, axis=1)
out["brain_met_event"] = out["any_brain_met"]

# ---- Survival ----
os_status = get("OS_STATUS")
os_months = pd.to_numeric(get("OS_MONTHS"), errors="coerce")
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

pfs_status = get("PFS_STATUS")
pfs_months = pd.to_numeric(get("PFS_MONTHS"), errors="coerce")
out["PFS_I_ADV_STATUS"]   = pfs_status
out["PFS_I_ADV_MONTHS"]   = pfs_months
out["pfs_i_adv_status"]   = pfs_status
out["tt_pfs_i_adv_days"]  = pfs_months * DAYS_PER_MONTH
out["tt_pfs_i_adv_mos"]   = pfs_months
out["tt_pfs_i_adv_yrs"]   = pfs_months / 12.0
out["PFS_imaging_months"] = pfs_months
out["pfs_i_event_bin"]    = pfs_status.map(first_char_bin).astype("Int64")

out["PFS_M_ADV_STATUS"]   = np.nan
out["PFS_M_ADV_MONTHS"]   = np.nan
out["pfs_m_adv_status"]   = np.nan
out["tt_pfs_m_adv_days"]  = np.nan
out["tt_pfs_m_adv_mos"]   = np.nan
out["tt_pfs_m_adv_yrs"]   = np.nan
out["PFS_medonc_months"]  = np.nan
out["pfs_m_event_bin"]    = pd.NA

out["pfs_cohort"] = np.nan
out["pfs_i_or_m_adv_status"]  = np.nan
out["tt_pfs_i_or_m_adv_days"] = np.nan
out["tt_pfs_i_or_m_adv_mos"]  = np.nan
out["tt_pfs_i_or_m_adv_yrs"]  = np.nan
out["pfs_i_and_m_adv_status"] = np.nan
out["tt_pfs_i_and_m_adv_days"] = np.nan
out["tt_pfs_i_and_m_adv_mos"]  = np.nan
out["tt_pfs_i_and_m_adv_yrs"]  = np.nan

# tt_brain_met_mos: if event, use the brain-recurrence time; else censor at OS
brain_cns_mos = out["dx_to_dist_mets_brain_cns_mos"]
out["tt_brain_met_mos"] = np.where(
    out["any_brain_met"] == 1,
    brain_cns_mos.where(brain_cns_mos.notna(), np.nan),
    os_months,
)
out["time_to_brain_met_mos"] = out["tt_brain_met_mos"]

# DFS / DSS (TCGA bonus)
out["dfs_status"] = get("DFS_STATUS")
out["dfs_months"] = pd.to_numeric(get("DFS_MONTHS"), errors="coerce")
out["dss_status"] = get("DSS_STATUS")
out["dss_months"] = pd.to_numeric(get("DSS_MONTHS"), errors="coerce")

# ---- Sample / panel metadata ----
out["ONCOTREE_CODE"]         = get("ONCOTREE_CODE")
out["SAMPLE_TYPE_DETAILED"]  = get("SAMPLE_TYPE")
out["sample_type"]           = get("SAMPLE_TYPE")
out["SEQ_ASSAY_ID"]          = np.nan
out["AGE_AT_SEQUENCING"]     = pd.to_numeric(get("AGE"), errors="coerce")
out["PDL1_POSITIVE_ANY"]     = np.nan
out["PDL1_TESTING"]          = np.nan
out["CPT_SEQ_DATE"]          = np.nan
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

# ---- Non-index cancer (TCGA mostly 1-cancer) ----
out["non_idx_n_cancers"] = 0
out["non_idx_any_heme"]  = pd.NA
out["non_idx_any_brain"] = pd.NA
out["non_idx_types"]     = ""
prior_dx = get("PRIOR_DX")
out["had_prior_non_breast_cancer"] = prior_dx.map(
    lambda v: pd.NA if pd.isna(v) else (0 if str(v).lower().startswith("no") else 1)
).astype("Int64")

# ---- Mutation features attach ----
out = out.merge(mut_feat, on="SAMPLE_ID", how="left")

# ---- Chemo / radiation hx (your addition) ----
out["hx_chemo_neoadj_bin"]   = get("HISTORY_NEOADJUVANT_TRTYN").map(yn1).astype("Int64")
hxc = mc.get("hx_chemo_any_bin")
out["hx_chemo_any_bin"]      = pd.to_numeric(hxc, errors="coerce").astype("Int64") if hxc is not None else pd.NA
hxa = mc.get("hx_chemo_adjuvant_bin")
out["hx_chemo_adjuvant_bin"] = pd.to_numeric(hxa, errors="coerce").astype("Int64") if hxa is not None else pd.NA
out["hx_radiation_bin"]      = get("RADIATION_THERAPY").map(yn1).astype("Int64")
note("hx_chemo_neoadj_bin","HISTORY_NEOADJUVANT_TRTYN","Y/N -> 1/0","TCGA neoadjuvant yes/no")
note("hx_chemo_any_bin","GDC treatments.treatment_type ~ /chemo|pharmaceutical|cytotox/","1 if any chemo row","Rolled from GDC treatments by patient")
note("hx_chemo_adjuvant_bin","GDC treatments.treatment_intent_type ~ /adjuv/","1 if adjuvant chemo","Excludes neoadjuvant")
note("hx_radiation_bin","RADIATION_THERAPY","Y/N -> 1/0","")

# Regimens parity placeholders
out["N_REGIMENS_PT"] = np.nan
out["n_regimens_pt"] = np.nan
out["ca_n_regimens"] = np.nan
out["n_cpt_pt"] = 1

# ---- Gene binary + pathways (copy from M) ----
KNOWN_CLINICAL = {
    "PATIENT_ID","SAMPLE_ID","patient_barcode","SUBTYPE","AGE","SEX","RACE","ETHNICITY",
    "OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS","PFS_STATUS","PFS_MONTHS",
    "DSS_STATUS","DSS_MONTHS","AJCC_PATHOLOGIC_TUMOR_STAGE","PATH_T_STAGE","PATH_N_STAGE",
    "PATH_M_STAGE","ICD_O_3_HISTOLOGY","ICD_O_3_SITE","ONCOTREE_CODE","CANCER_TYPE",
    "CANCER_TYPE_DETAILED","TUMOR_TYPE","GRADE","TISSUE_PROSPECTIVE_COLLECTION_INDICATOR",
    "TISSUE_RETROSPECTIVE_COLLECTION_INDICATOR","TISSUE_SOURCE_SITE_CODE","TUMOR_TISSUE_SITE",
    "ANEUPLOIDY_SCORE","SAMPLE_TYPE","MSI_SCORE_MANTIS","MSI_SENSOR_SCORE","SOMATIC_STATUS",
    "TMB_NONSYNONYMOUS","TISSUE_SOURCE_SITE","TBL_SCORE","OTHER_PATIENT_ID",
    "CANCER_TYPE_ACRONYM","AJCC_STAGING_EDITION","BUFFA_HYPOXIA_SCORE","DAYS_LAST_FOLLOWUP",
    "DAYS_TO_BIRTH","DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS","FORM_COMPLETION_DATE",
    "GENETIC_ANCESTRY_LABEL","HISTORY_NEOADJUVANT_TRTYN","ICD_10","IN_PANCANPATHWAYS_FREEZE",
    "INFORMED_CONSENT_VERIFIED","NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT",
    "PERSON_NEOPLASM_CANCER_STATUS","PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT","PRIOR_DX",
    "RADIATION_THERAPY","RAGNUM_HYPOXIA_SCORE","SAMPLE_COUNT","WEIGHT","WINTER_HYPOXIA_SCORE",
    "hx_chemo_any_bin","hx_chemo_neoadj_C","hx_chemo_adjuvant_bin",
}

gene_cols = []
for col in m.columns:
    if col in KNOWN_CLINICAL: continue
    if col.startswith("pathway_"):
        gene_cols.append(col); continue
    if col.endswith(".Amp") or col.endswith(".Del") or col.endswith(".fus"):
        gene_cols.append(col); continue
    if isinstance(col, str) and col.isupper() and not col.endswith("_ID") and not col.endswith("_STATUS"):
        if re.fullmatch(r"[A-Z0-9][A-Z0-9\-]*", col):
            gene_cols.append(col)
gene_existing = [g for g in gene_cols if g not in out.columns]
print(f"   attaching {len(gene_existing)} gene/pathway columns")
out = out.merge(m[["SAMPLE_ID"] + gene_existing], on="SAMPLE_ID", how="left")

# ---- Top-gene pipeline (cohort-specific) ----
print(">>> Computing top-gene pipeline…")
mut_only = [c for c in gene_existing
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
print(f"   TOP10 (brain-met cohort, N={len(brain)}): {TOP10}")
with open(OUT_TOPGENES, "w") as fh:
    fh.write("\n".join(TOP10) + "\n")

for g in TOP10:
    out[f"G_top10_{g}"] = out[g].fillna(0).astype(int)
out["top5_n_mutated"]   = out[TOP5].fillna(0).sum(axis=1).astype(int)
out["top10_n_mutated"]  = out[TOP10].fillna(0).sum(axis=1).astype(int)
out["top5_any_mutated"] = (out["top5_n_mutated"] > 0).astype(int)
out["top10_any_mutated"]= (out["top10_n_mutated"] > 0).astype(int)

# ---- Save ----
out.to_csv(OUT_MAIN, index=False)
pd.DataFrame(dict_rows).to_csv(OUT_DICT, index=False)

# ============================================================
# Validation (§17)
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
print(f"\nChemo / radiation hx (sample counts):")
for c in ["hx_chemo_neoadj_bin","hx_chemo_any_bin","hx_chemo_adjuvant_bin","hx_radiation_bin"]:
    s = out[c]
    print(f"  {c}: 1={int((s==1).sum())}, 0={int((s==0).sum())}, NA={int(s.isna().sum())}")

print(f"\nWrote {OUT_MAIN.name}")
print(f"Wrote {OUT_TOPGENES.name}")
print(f"Wrote {OUT_GENEPREV.name}")
print(f"Wrote {OUT_DICT.name}")
