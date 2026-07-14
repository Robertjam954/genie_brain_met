"""Post-harmonization enrichment for the 3 cohort frames.

Adds three layers on top of extracted_variables_<cohort>_data.csv:

1. Long-form alias columns - the human-readable GENIE BPC display names
   (e.g. `patient_id`, `year_of_birth`, `ajcc_stage`, `histology`,
   `overall_survival_months`, `distant_mets_brain`) point to the existing
   canonical columns (`record_id`, `BIRTH_YEAR`, `stage_dx`, `ca_histology`,
   `OS_months`, `dist_mets_brain_cns`). A column truly missing in a cohort
   is created and filled with NaN per the spec.

2. Competing-event time-to-event for brain mets (Aim 2 supplementary):
   `tt_compete_first_dmets_mos` = earliest dx_to_dist_mets_<organ>_mos across
   non-brain organs (NaN if no non-brain met)
   `compete_first_dmets_event` = 1 if any non-brain dist_mets_<organ> is 1
   `tt_brain_met_competing_mos` = min(brain time, first non-brain time, OS time)
   `brain_met_competing_status` = 0 censor / 1 brain / 2 competing event

3. Regimen-derived treatment flags (GENIE BPC only - TCGA/MSK18 get NaN):
   chemo_overall, endocrine_overall, aromatase_inhibitor_overall,
   tamoxifen_overall, cdk4_6_inhibitor_overall, akt_inhibitor_overall,
   neoadj_chemo_or_rt, received_adjuvant_ai, received_adjuvant_tam_therapy.
   Computed from genie_bpc_v1_regimens.csv via drug-name regex.

Reads/writes:
  src/exploratory data analysis/extracted_variables_<cohort>_data.csv
  (in-place rewrite). Also appends rows to the dictionary.csv.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/robertjames/Documents/Documents - Robert’s iMac/"
            "Research Projects/MSKCC Research Fellowship/Projects/"
            "genie_tcga_impact_brain_mets")
PROC = ROOT / "data" / "processed"
OUT_DIR = ROOT / "src" / "exploratory data analysis"

ALL_DIST_SITES = [
    "abdomen", "adrenal", "bone", "bone_marrow", "brain_cns", "breast",
    "head_and_neck", "liver", "lymph_nodes", "other", "pelvis",
    "pericardial_and_malignant_pericardial_effusion",
    "peritoneum_and_malignant_peritoneal_effusion",
    "pleura_and_malignant_pleural_effusion", "pulmonary", "skin", "thorax",
]
NON_BRAIN_SITES = [s for s in ALL_DIST_SITES if s != "brain_cns"]


# ---------- (1) Long-form alias map ----------
# Each entry: long-form-name -> ordered list of candidate canonical sources.
# First found candidate is used; if none present, the long-form column is
# created and filled with NaN (per harmonization spec).
ALIAS_MAP: dict[str, list[str]] = {
    # Identifiers
    "study_id": ["STUDY_ID", "cohort", "study"],
    "patient_id": ["record_id", "PATIENT_ID", "patient_id"],
    "sample_id": ["SAMPLE_ID", "sample_id"],
    # Demographics
    "center": ["CENTER", "institution"],
    "sex": ["SEX", "naaccr_sex_code"],
    "year_of_birth": ["BIRTH_YEAR", "birth_year", "year_birth"],
    "primary_race": ["PRIMARY_RACE", "naaccr_race_code_primary"],
    "secondary_race": ["SECONDARY_RACE", "naaccr_race_code_secondary"],
    "tertiary_race": ["TERTIARY_RACE", "naaccr_race_code_tertiary"],
    "ethnicity_category": ["ETHNICITY", "naaccr_ethnicity_code"],
    "region_of_patient_origin": ["region_of_patient_origin"],
    # Oncology history
    "age_at_diagnosis": ["age_dx_num", "age_dx", "AGE_AT_DIAGNOSIS"],
    "age_at_initial_diagnosis_years": ["age_dx_num", "age_dx"],
    "age_at_primary_diagnosis": ["age_dx_num", "age_dx"],
    "was_a_biopsy_performed_of_metastatic_site": ["was_a_biopsy_performed_of_metastatic_site"],
    "oncotree_code": ["ONCOTREE_CODE", "cpt_oncotree_code", "oncotree_code"],
    "cancer_type": ["CANCER_TYPE"],
    "cancer_type_detailed": ["CANCER_TYPE_DETAILED"],
    "tumor_registry_icd_o_3_behavior_code": ["naaccr_behavior_cd"],
    "additional_breast_cancer_diagnoses": ["additional_breast_cancer_diagnoses"],
    "other_oncotree_diagnosis": ["other_oncotree_diagnosis"],
    "other_oncotree_diagnosis_1": ["other_oncotree_diagnosis_1"],
    "other_oncotree_diagnosis_2": ["other_oncotree_diagnosis_2"],
    "other_oncotree_diagnosis_3": ["other_oncotree_diagnosis_3"],
    "number_of_other_invasive_cancer_diagnosis": ["n_cancers_index", "non_idx_n_cancers"],
    "other_invasive_cancer_diagnosis": ["non_idx_types"],
    "number_of_cancers_any_type": ["n_cancers"],
    "number_of_bpc_project_cancers_index_cancers": ["n_cancers_index"],
    # Tumor characteristics
    "ki67_percent": ["ki67_percent", "KI67_PERCENT"],
    "oncotype_dx_score": ["ca_bca_oncotypedx", "CA_BCA_ONCOTYPEDX"],
    "ajcc_stage": ["stage_dx", "STAGE_DX", "best_ajcc_stage_cd"],
    "stage_at_diagnosis": ["stage_dx", "STAGE_DX"],
    "grade": ["ca_grade", "CA_GRADE", "OVERALL_TUMOR_GRADE"],
    "her2_status": ["ca_bca_her_summ", "CA_BCA_HER_SUMM", "HER2_STATUS"],
    "hormone_receptor_status": ["bca_subtype"],
    "er_pr_receptor_change": ["er_pr_receptor_change"],
    "discordance_receptor_status_with_primary": ["discordance_receptor_status_with_primary"],
    "her2_receptor_status_at_metastatic_diagnosis": ["her2_receptor_status_at_metastatic_diagnosis"],
    "her2_receptor_status_at_initial_diagnosis": ["ca_bca_her_summ"],
    "hr_status_at_metastatic_diagnosis": ["hr_status_at_metastatic_diagnosis"],
    "hr_status_at_initial_diagnosis": ["bca_subtype"],
    "site_of_sample_tested": ["SAMPLE_TYPE_DETAILED", "sample_site"],
    "stage_at_initial_breast_cancer_diagnosis": ["stage_dx"],
    "subtype_of_metastatic_diagnosis": ["subtype_of_metastatic_diagnosis"],
    "subtype_of_initial_diagnosis": ["bca_subtype"],
    "histology": ["ca_histology", "CA_HISTOLOGY"],
    "histology_category": ["ca_hist_brca"],
    # Tumor + sequence characteristics
    "sequence_assay_id": ["SEQ_ASSAY_ID", "cpt_seq_assay_id"],
    "age_at_which_sequencing_was_reported": ["AGE_AT_SEQUENCING"],
    "age_at_sample_collection_4": ["AGE_AT_SEQUENCING"],
    "age_at_sequencing": ["AGE_AT_SEQUENCING"],
    "sequenced_sample": ["SAMPLE_ID"],
    "how_long_sample_was_collected_after_starting_first_line_systemic_therapy_30_days_30_days_treatment_naive": [
        "how_long_sample_was_collected_after_starting_first_line_systemic_therapy_30_days_30_days_treatment_naive"
    ],
    "sample_site": ["SAMPLE_TYPE_DETAILED", "SAMPLE_SITE"],
    "number_of_samples_per_patient": ["number_of_samples_per_patient"],
    "sample_type": ["SAMPLE_TYPE_DETAILED", "SAMPLE_TYPE"],
    "sequencing_method": ["SEQ_ASSAY_ID"],
    # Genetic alterations
    "fraction_genome_altered": ["FRACTION_GENOME_ALTERED", "fraction_genome_altered"],
    "mutation_count": ["mutation_count_all_sites_sum"],
    "mismatch_repair_mmr_testing_at_time_of_sample_acquisition": ["mismatch_repair_mmr_testing_at_time_of_sample_acquisition"],
    "msi_h_test_result_at_time_of_sample_acquisition": ["msi_h_test_result_at_time_of_sample_acquisition"],
    "microsatellite_instability_msi_testing_at_time_of_sample_acquisition": ["microsatellite_instability_msi_testing_at_time_of_sample_acquisition"],
    "tumor_mutational_burden": ["TMB_NONSYNONYMOUS", "TMB"],
    "tmb_nonsynonymous": ["TMB_NONSYNONYMOUS"],
    "x1p_19q_codeletion": ["x1p_19q_codeletion"],
    "idh1_2_mutation": ["idh1_2_mutation"],
    "integrated_histomolecular_group": ["integrated_histomolecular_group"],
    "mgmt_methylation_status": ["mgmt_methylation_status"],
    # Survival outcomes
    "patients_age_of_death_in_days": ["age_death_yrs"],
    "age_at_first_distant_metastasis_in_days": ["dx_to_dmets_days"],
    "time_from_initial_diagnosis_to_metastatic_diagnosis_months": ["dx_to_dmets_mos"],
    "cause_of_death": ["cause_of_death"],
    "patients_vital_status": ["os_dx_status", "os_status_f"],
    "vital_status": ["os_status_f", "os_dx_status"],
    "overall_survival_months": ["OS_months", "tt_os_dx_mos"],
    "overall_survival_status": ["os_dx_status", "OS_STATUS"],
    "patients_age_of_last_follow_up_in_days": ["last_anyvisit_int", "last_alive_int"],
    "interval_in_days_from_dob_to_date_of_last_contact": ["dob_lastalive_int", "last_alive_int"],
    "interval_in_days_from_dob_to_dod": ["hybrid_death_int"],
    "sample_class": ["sample_class"],
    "year_of_last_contact": ["year_of_last_contact"],
    "year_of_death": ["year_of_death"],
    "pfs_i_from_diagnosis_status": ["pfs_i_adv_status"],
    "pfs_m_from_diagnosis_status": ["pfs_m_adv_status"],
    # LRR / metastatic time
    "age_at_metastatic_diagnosis_years": ["dx_to_dmets_yrs"],
    "lrr_date_interval_days": ["lrr_date_interval_days"],
    "lrr_radiation_therapy": ["lrr_radiation_therapy"],
    "site_of_lrr": ["site_of_lrr"],
    "distant_metastatic_disease_interval": ["dx_to_dmets_days"],
    # Distant mets at dx (organ-level)
    "distant_mets_adrenal": ["dist_mets_adrenal", "DMETS_DX_ADRENAL"],
    "distant_mets_bone": ["dist_mets_bone", "DMETS_DX_BONE"],
    "distant_mets_brain": ["dist_mets_brain_cns", "DMETS_DX_BRAIN"],
    "distant_mets_liver": ["dist_mets_liver", "DMETS_DX_LIVER"],
    "distant_mets_lung": ["dist_mets_pulmonary", "DMETS_DX_LUNG"],
    "distant_mets_lymph_nodes": ["dist_mets_lymph_nodes", "DMETS_DX_LYMPH"],
    "distant_mets_other": ["dist_mets_other", "DMETS_DX_OTHER"],
    "distant_mets_pleura": ["dist_mets_pleura_and_malignant_pleural_effusion", "DMETS_DX_PLEURA"],
    "distant_mets_subcutaneous_tissue": ["dist_mets_skin", "DMETS_DX_SUBC_TISSUE"],
    "age_at_metastatic_diagnosis": ["dx_to_dmets_mos"],
    "met_site_liver": ["dist_mets_liver"],
    "has_this_patient_experienced_a_loco_regional_recurrence_lrr": ["has_this_patient_experienced_a_loco_regional_recurrence_lrr"],
    "met_site_lung": ["dist_mets_pulmonary"],
    "met_site_lymph_node": ["dist_mets_lymph_nodes"],
    "total_number_of_therapies_received_in_metastatic_disease_treatment": [
        "total_number_of_therapies_received_in_metastatic_disease_treatment", "ca_n_regimens"],
    "interval_between_sequencing_and_metastatis_daignosis": ["interval_between_sequencing_and_metastatis_daignosis"],
    "met_site_multiple_sites": ["met_site_multiple_sites"],
    "met_site_other_site": ["dist_mets_other"],
    "met_site_bone_only": ["met_site_bone_only"],
    "met_site_bone": ["dist_mets_bone"],
    "met_site_brain": ["dist_mets_brain_cns"],
    "met_site_soft_tissue": ["dist_mets_skin"],
    "met_site_visceral": ["met_site_visceral"],
    "sites_of_distant_metastasis_at_the_time_of_cancer_diagnosis_stage_iv_patients": [
        "sites_of_distant_metastasis_at_the_time_of_cancer_diagnosis_stage_iv_patients"],
    "year_of_next_generation_sequencing": ["year_of_next_generation_sequencing", "CPT_SEQ_DATE"],
    # Patient id supp table (genieBPC merge key)
    "patient_id_in_supp_table_1": ["record_id", "PATIENT_ID"],
}


def add_aliases(df: pd.DataFrame, dict_rows: list[dict]) -> pd.DataFrame:
    cols = set(df.columns)
    for long_name, candidates in ALIAS_MAP.items():
        if long_name in cols:
            continue  # already present, don't overwrite
        src = next((c for c in candidates if c in cols), None)
        if src is not None:
            df[long_name] = df[src]
            dict_rows.append({
                "variable": long_name, "original": src, "mapped": "alias",
                "note": "long-form alias for canonical source",
            })
        else:
            df[long_name] = np.nan
            dict_rows.append({
                "variable": long_name, "original": "(absent in cohort)",
                "mapped": "NaN", "note": "long-form var not present in source",
            })
    return df


# ---------- (2) Competing-event TTE ----------
def add_competing_event_tte(df: pd.DataFrame, dict_rows: list[dict]) -> pd.DataFrame:
    """Brain-met TTE with competing events from earliest non-brain dist_mets."""
    # Earliest non-brain dist_mets time
    non_brain_time_cols = [f"dx_to_dist_mets_{s}_mos" for s in NON_BRAIN_SITES
                          if f"dx_to_dist_mets_{s}_mos" in df.columns]
    non_brain_event_cols = [f"dist_mets_{s}" for s in NON_BRAIN_SITES
                          if f"dist_mets_{s}" in df.columns]

    if not non_brain_time_cols or not non_brain_event_cols:
        df["tt_compete_first_dmets_mos"] = np.nan
        df["compete_first_dmets_event"] = np.nan
        df["tt_brain_met_competing_mos"] = pd.to_numeric(df.get("tt_brain_met_mos"), errors="coerce")
        df["brain_met_competing_status"] = pd.to_numeric(df.get("brain_met_event", 0), errors="coerce")
        dict_rows.append({"variable": "tt_brain_met_competing_mos",
                          "original": "tt_brain_met_mos",
                          "mapped": "brain-only TTE (no non-brain sites in cohort)",
                          "note": ""})
        return df

    times = df[non_brain_time_cols].apply(pd.to_numeric, errors="coerce")
    events = df[non_brain_event_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    # mask times where event=0 so they don't count
    site_to_time = {s: f"dx_to_dist_mets_{s}_mos" for s in NON_BRAIN_SITES
                    if f"dx_to_dist_mets_{s}_mos" in df.columns
                    and f"dist_mets_{s}" in df.columns}
    site_to_event = {s: f"dist_mets_{s}" for s in NON_BRAIN_SITES
                    if f"dist_mets_{s}" in df.columns
                    and f"dx_to_dist_mets_{s}_mos" in df.columns}
    masked = pd.DataFrame({s: np.where(events[site_to_event[s]] == 1,
                                       times[site_to_time[s]], np.nan)
                          for s in site_to_time})
    df["tt_compete_first_dmets_mos"] = masked.min(axis=1)
    df["compete_first_dmets_event"] = (events[list(site_to_event.values())].max(axis=1) == 1).astype(int)

    # Composite multi-state status: 0 censor, 1 brain, 2 other met
    brain_time = pd.to_numeric(df.get("dx_to_dist_mets_brain_cns_mos"), errors="coerce")
    brain_event = pd.to_numeric(df.get("dist_mets_brain_cns", 0), errors="coerce").fillna(0)
    other_time = df["tt_compete_first_dmets_mos"]
    other_event = df["compete_first_dmets_event"].astype(float)
    os_time = pd.to_numeric(df.get("tt_os_dx_mos"), errors="coerce")

    # Pick the earliest of brain / other; if neither, censor at OS time
    candidates = pd.DataFrame({
        "brain": brain_time.where(brain_event == 1, np.nan),
        "other": other_time.where(other_event == 1, np.nan),
        "censor": os_time,
    })
    df["tt_brain_met_competing_mos"] = candidates.min(axis=1)
    first = candidates.idxmin(axis=1)
    df["brain_met_competing_status"] = first.map({"brain": 1, "other": 2, "censor": 0}).fillna(0).astype(int)

    dict_rows += [
        {"variable": "tt_compete_first_dmets_mos",
         "original": "min(dx_to_dist_mets_<organ>_mos) across non-brain sites where dist_mets_<organ>==1",
         "mapped": "float months", "note": "earliest non-brain distant met"},
        {"variable": "compete_first_dmets_event",
         "original": "any(dist_mets_<organ>==1) across non-brain sites",
         "mapped": "0/1", "note": "competing event indicator"},
        {"variable": "tt_brain_met_competing_mos",
         "original": "min(brain_time if event, other_time if event, OS_time)",
         "mapped": "float months",
         "note": "multi-state time-to-first-event for Aim 2 competing-risks model"},
        {"variable": "brain_met_competing_status",
         "original": "argmin among brain/other/censor",
         "mapped": "0 censor / 1 brain / 2 other",
         "note": "competing-risks status for Aim 2 (cmprsk / FineGray)"},
    ]
    return df


# ---------- (3) GENIE BPC regimen-derived treatment flags ----------
DRUG_PATTERNS = {
    "chemo_overall": r"capecitabine|cyclophosphamide|paclitaxel|docetaxel|epirubicin|doxorubicin|carboplatin|cisplatin|gemcitabine|5-fluorouracil|fluorouracil|eribulin|vinorelbine|methotrexate|ifosfamide|etoposide|irinotecan|oxaliplatin|sacituzumab|trastuzumab[- ]?deruxtecan|sn-38",
    "endocrine_overall": r"tamoxifen|fulvestrant|toremifene|raloxifene|anastrozole|letrozole|exemestane|goserelin|leuprolide",
    "aromatase_inhibitor_overall": r"anastrozole|letrozole|exemestane",
    "tamoxifen_overall": r"tamoxifen|toremifene|raloxifene",
    "cdk4_6_inhibitor_overall": r"palbociclib|ribociclib|abemaciclib",
    "akt_inhibitor_overall": r"capivasertib|ipatasertib|miransertib|afuresertib",
    "her2_targeted_overall": r"trastuzumab|pertuzumab|lapatinib|neratinib|tucatinib|t-dm1|ado-trastuzumab|emtansine|trastuzumab[- ]?deruxtecan",
    "parp_inhibitor_overall": r"olaparib|talazoparib|niraparib|rucaparib|veliparib",
    "immuno_overall": r"pembrolizumab|atezolizumab|nivolumab|ipilimumab|durvalumab",
}


def add_genie_regimen_flags(df: pd.DataFrame, dict_rows: list[dict]) -> pd.DataFrame:
    """Compute per-patient regimen-derived flags from regimens.csv. Patient-level
    broadcast to sample rows via record_id."""
    if "record_id" not in df.columns:
        for col in DRUG_PATTERNS:
            df[col] = np.nan
        return df

    reg_path = PROC / "genie_bpc_v1_regimens.csv"
    if not reg_path.exists():
        for col in DRUG_PATTERNS:
            df[col] = np.nan
        return df

    reg = pd.read_csv(reg_path, low_memory=False, usecols=["record_id", "regimen_drugs"])
    reg["regimen_drugs_lc"] = reg["regimen_drugs"].astype(str).str.lower()

    pt_flags = pd.DataFrame({"record_id": reg["record_id"].unique()})
    for col, pattern in DRUG_PATTERNS.items():
        mask = reg["regimen_drugs_lc"].str.contains(pattern, regex=True, na=False)
        any_per_pt = reg.assign(_hit=mask.astype(int)).groupby("record_id")["_hit"].max()
        pt_flags[col] = pt_flags["record_id"].map(any_per_pt).fillna(0).astype(int)
        dict_rows.append({
            "variable": col, "original": f"regex on regimen_drugs: {pattern[:60]}...",
            "mapped": "0/1", "note": "any regimen with matching drug (GENIE BPC only)",
        })

    df = df.merge(pt_flags, on="record_id", how="left")
    for col in DRUG_PATTERNS:
        df[col] = df[col].fillna(0).astype(int)
    return df


# ---------- driver ----------
COHORTS = {
    "genie": "extracted_variables_genie_data.csv",
    "tcga": "extracted_variables_tcga_data.csv",
    "msk18": "extracted_variables_breast_msk_2018_data.csv",
}


def run(cohort: str) -> None:
    fname = OUT_DIR / COHORTS[cohort]
    dict_fname = fname.with_name(fname.stem.replace("_data", "_dictionary") + ".csv")
    print(f"\n=== enrich {cohort}: {fname.name} ===")
    df = pd.read_csv(fname, low_memory=False)
    n_before = df.shape[1]
    dict_rows: list[dict] = []

    df = add_aliases(df, dict_rows)
    df = add_competing_event_tte(df, dict_rows)
    if cohort == "genie":
        df = add_genie_regimen_flags(df, dict_rows)
    else:
        # TCGA / MSK18: BPC-specific treatment flags not derivable -> NaN cols
        for col in DRUG_PATTERNS:
            df[col] = np.nan

    n_after = df.shape[1]
    print(f"  added {n_after - n_before} new columns ({n_before} -> {n_after})")

    df.to_csv(fname, index=False)
    print(f"  wrote {fname.name}  shape={df.shape}")

    # Append to dictionary
    if dict_fname.exists():
        existing = pd.read_csv(dict_fname)
        out = pd.concat([existing, pd.DataFrame(dict_rows)], ignore_index=True)
    else:
        out = pd.DataFrame(dict_rows)
    out.to_csv(dict_fname, index=False)
    print(f"  wrote {dict_fname.name}  +{len(dict_rows)} dict rows")


def main() -> None:
    for c in COHORTS:
        run(c)


if __name__ == "__main__":
    main()
