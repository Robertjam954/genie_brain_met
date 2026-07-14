# Brain-met analytic dataset - harmonization spec

This document tells an LLM exactly how to extract and recode a new clinical-genomic
cancer cohort so it matches the canonical schema used in the GENIE BPC BRCA brain-met
analysis. Apply this to a NEW source dataset; output must be one row per SAMPLE_ID with
all the variables below.

If a source column does not exist in the new dataset, leave the derived variable as NA
for all rows but still create the column (so frames are conformable). Never invent
values. Document any mapping ambiguity in a side-car `harmonization_decisions.csv` file
with columns `variable, source_col, decision, sample_value`.

---

## 1. Required identifiers (column names in output)

| Output col | Required | Description |
|---|---|---|
| `SAMPLE_ID` | yes | Sample-level barcode (1 row per SAMPLE_ID) |
| `record_id` | yes | Patient-level id (synonyms: PATIENT_ID, Patient_Id) |
| `ca_seq` | if exists | Cancer sequence number per patient (for multi-cancer patients) |
| `cpt_number` | if exists | Panel test number |

Output grain: ONE ROW PER SAMPLE_ID. If a patient has multiple samples, replicate
patient-level clinical attributes across their sample rows. If multiple cancers per
patient exist, attach the INDEX cancer (`ca_seq == 1` or equivalent).

---

## 2. Clinical demographics + breast subtype

### Demographics

| Output col | Source col candidates | Derivation |
|---|---|---|
| `SEX` | `SEX`, `gender`, `sex` | Direct copy. Canonical values: `Female`, `Male`, `Other`, NA |
| `BIRTH_YEAR` | `BIRTH_YEAR`, `year_birth` | Direct copy (int) |
| `CENTER` | `CENTER`, `institution`, `site` | Direct copy |
| `PRIMARY_RACE` | `PRIMARY_RACE`, `naaccr_race_code_primary`, `race` | Direct copy of raw label |
| `ETHNICITY` | `ETHNICITY`, `naaccr_ethnicity_code` | Direct copy of raw label |
| `race_clean` | `PRIMARY_RACE` | Regex collapse (see Section 8) |
| `ethnicity_clean` | `ETHNICITY` | Regex collapse (see Section 8) |

### Breast subtype + receptor

| Output col | Source col candidates | Derivation |
|---|---|---|
| `bca_subtype` | `bca_subtype`, `BCA_SUBTYPE`, derived from ER/PR/HER2 | Canonical: `HR+, HER2-`, `HR+, HER2+`, `HR-, HER2+`, `TNBC` |
| `receptor_primary_cat` | `bca_subtype` | Recode to: `HR+/HER2-`, `HR+/HER2+`, `HR-/HER2+`, `Triple Negative` (see Section 8) |
| `her2_status_bin` | `ca_bca_her_summ`, `HER2_STATUS`, `OVERALL_HER2_STATUS` | 1 if matches `positive\|amplif\|equivocal_pos`; 0 if `negative\|not amplif`; else NA |

---

## 3. Histology / grade

| Output col | Source col candidates | Derivation |
|---|---|---|
| `ca_histology` | `ca_histology`, `CA_HISTOLOGY`, `histology` | Direct copy |
| `ca_hist_brca` | `ca_hist_brca` | Breast-specific histology code |
| `ca_grade` | `ca_grade`, `CA_GRADE`, `OVERALL_TUMOR_GRADE` | Direct copy |
| `grade_ord` | `ca_grade` | Ordered factor `Low < Intermediate < High` (see Section 8) |

---

## 4. Stage at diagnosis

| Output col | Source col candidates | Derivation |
|---|---|---|
| `stage_dx` | `stage_dx`, `STAGE_DX`, `STAGE_AT_DIAGNOSIS`, `best_ajcc_stage_cd` | Direct copy |
| `stage_dx_iv` | `stage_dx_iv` (if exists) OR derive: `stage_dx == "Stage IV"` | Canonical: `Stage IV` / `Stage I-III` |
| `stage_diag_group` | `stage_dx` | Ordered factor mapping (see Section 8): `Stage I`, `Stage II`, `Stage III`, `Stage IV` |
| `stage_iv_bin` | `stage_dx_iv` | int 0/1: 1 if `lower() in {"stage iv","iv","yes","true","1"}` |

---

## 5. Mets at diagnosis (organ-specific binary flags)

ALL of the following are binary indicators captured at time of dx. Required if source has them.

```
DMETS_DX_BRAIN, DMETS_DX_BONE, DMETS_DX_LIVER, DMETS_DX_LUNG,
DMETS_DX_ADRENAL, DMETS_DX_LYMPH, DMETS_DX_PLEURA, DMETS_DX_SUBC_TISSUE,
DMETS_DX_OTHER, CA_DMETS_YN
```

If using cancer-level features (preferred), also extract:

```
ca_dmets_yn, ca_first_dmets1, ca_first_dmets2, ..., ca_first_dmets10
```

---

## 6. Met DEVELOPMENT during study (cancer-level, organ-specific)

These cover mets developing AFTER dx. For each of the 17 sites below, output two columns:
`dist_mets_<site>` (binary 0/1) and `dx_to_dist_mets_<site>_mos` (float, months from dx).

Canonical site list:

```
abdomen, adrenal, bone, bone_marrow, brain_cns, breast, head_and_neck, liver,
lymph_nodes, other, pelvis, pericardial_and_malignant_pericardial_effusion,
peritoneum_and_malignant_peritoneal_effusion,
pleura_and_malignant_pleural_effusion, pulmonary, skin, thorax
```

Also include:

| Output col | Source col candidates | Derivation |
|---|---|---|
| `dmets_post_dx` | `dmets_post_dx` | 0/1 any post-dx met |
| `dx_to_dmets_mos` | `dx_to_dmets_mos` | time from dx to FIRST distant met (months) |

---

## 7. Brain-met cohort (DERIVED - the analytic cohort split)

| Output col | Derivation |
|---|---|
| `any_brain_met` | `(dist_mets_brain_cns == 1) OR (DMETS_DX_BRAIN.lower() == "yes")`. int 0/1 |
| `brain_met_at_dx` | `DMETS_DX_BRAIN.lower() == "yes"`. int 0/1 |
| `met_loc` | factor (Brain / Other / None): `Brain` if `any_brain_met == 1`; else `Other` if ANY other `dist_mets_<organ>` is 1; else `None`. Use unordered factor with levels `["Brain", "Other", "None"]` |
| `tt_brain_met_mos` | `dx_to_dist_mets_brain_cns_mos` if `any_brain_met == 1`, else `tt_os_dx_mos` (overall follow-up time used as censoring time). float |
| `brain_met_event` | copy of `any_brain_met` (for survival analysis) |

---

## 8. Categorical recodings (exact mappings)

### `race_clean` (from raw race string)
| Raw value (regex match, case insensitive) | -> |
|---|---|
| `White` or `Middle East` | `White` |
| `Black` or `African` | `Black` |
| `Asian`, `Indian`, `Chinese` | `Asian` |
| `Native`, `Alaska` | `Native American` |
| (no match) | NA |

### `ethnicity_clean`
| Raw (regex, case insensitive) | -> |
|---|---|
| `Non-Spanish` or `Non-Hispanic` | `Non-Hispanic` |
| `Hispanic`, `Latino`, `Cuban`, `Mexican`, `Puerto` | `Hispanic` |
| (no match) | NA |

### `age_cat` (from continuous age at diagnosis)
| Source value | -> |
|---|---|
| `< 50` | `<50` |
| `50 <= age <= 70` | `50-70` |
| `> 70` | `>70` |
Ordered factor: `<50 < 50-70 < >70`.

### `grade_ord` (from `ca_grade`)
Regex (case insensitive) on `ca_grade` string:
| Raw match | -> |
|---|---|
| `Low` or `^\s*I\s*$` | `Low` |
| `Intermediate` or `^\s*II\s*$` | `Intermediate` |
| `High` or `^\s*III\s*$` | `High` |
Ordered factor: `Low < Intermediate < High`.

### `stage_diag_group` (from `stage_dx`)
Apply `.upper().strip()` first, then map:
| Raw | -> |
|---|---|
| `I`, `IA`, `IB`, `STAGE I` | `Stage I` |
| `II`, `IIA`, `IIB`, `STAGE II` | `Stage II` |
| `III`, `IIIA`, `IIIB`, `IIIC`, `STAGE III` | `Stage III` |
| `IV`, `STAGE IV` | `Stage IV` |
Ordered factor.

### `receptor_primary_cat` (from `bca_subtype`)
| Source | -> |
|---|---|
| `HR+, HER2-` or `HR+/HER2-` | `HR+/HER2-` |
| `HR+, HER2+` or `HR+/HER2+` | `HR+/HER2+` |
| `HR-, HER2+` or `HR-/HER2+` | `HR-/HER2+` |
| `TNBC`, `Triple Negative`, `HR-, HER2-`, `HR-/HER2-` | `Triple Negative` |
Unordered factor with levels: `["HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "Triple Negative"]`.

### `sample_type_bin`
| Source `SAMPLE_TYPE_DETAILED` (case insensitive contains) | -> |
|---|---|
| `primary` | 1 |
| `metast` | 2 |
| (other) | NA |

### `her2_status_bin`
| Source `ca_bca_her_summ` (regex, case insensitive) | -> |
|---|---|
| `positive` or `amplif` or `equivocal_pos` | 1 |
| `negative` or `not amplif` | 0 |
| (other) | NA |

---

## 9. Survival endpoints

### Overall survival (from diagnosis)

| Output col | Source col candidates | Derivation |
|---|---|---|
| `os_dx_status` | `os_dx_status`, `OS_STATUS` | cBioPortal coding `"0:LIVING"`, `"1:DECEASED"` |
| `tt_os_dx_mos` | `tt_os_dx_mos`, `OS_MONTHS` | float, months from dx |
| `OS_months` | `tt_os_dx_mos` | rename for analysis (duplicate) |
| `os_status_bin` | `os_dx_status` | `int(first_char)` -> 0/1 |
| `os_status_f` | `os_dx_status` | factor: first char `0` -> `Alive`, `1` -> `Deceased` |

### OS from advanced disease

| Output col | Source col candidates |
|---|---|
| `os_adv_status`, `tt_os_adv_mos` | direct copies |

### PFS (imaging-based, advanced disease)

| Output col | Source col candidates | Derivation |
|---|---|---|
| `pfs_i_adv_status` | `pfs_i_adv_status`, `PFS_I_ADV_STATUS` | cBioPortal `"0:CENSORED"` / `"1:PROGRESSION"` |
| `tt_pfs_i_adv_mos` | `tt_pfs_i_adv_mos`, `PFS_I_ADV_MONTHS` | float |
| `PFS_imaging_months` | `tt_pfs_i_adv_mos` | rename for analysis |
| `pfs_i_event_bin` | `pfs_i_adv_status` | `int(first_char)` -> 0/1 |

### PFS (medical oncology, advanced disease)

| Output col | Source col candidates | Derivation |
|---|---|---|
| `pfs_m_adv_status` | `pfs_m_adv_status`, `PFS_M_ADV_STATUS` | cBioPortal `"0:CENSORED"` / `"1:PROGRESSION"` |
| `tt_pfs_m_adv_mos` | `tt_pfs_m_adv_mos`, `PFS_M_ADV_MONTHS` | float |
| `PFS_medonc_months` | `tt_pfs_m_adv_mos` | rename |
| `pfs_m_event_bin` | `pfs_m_adv_status` | `int(first_char)` -> 0/1 |

### PFS combined (OR / AND of imaging + medonc)

`pfs_cohort`, `pfs_i_or_m_adv_status`, `tt_pfs_i_or_m_adv_mos`,
`pfs_i_and_m_adv_status`, `tt_pfs_i_and_m_adv_mos` - direct copies if present.

---

## 10. Sample / panel-test metadata

| Output col | Source col candidates |
|---|---|
| `ONCOTREE_CODE`, `SEQ_ASSAY_ID`, `AGE_AT_SEQUENCING`, `SAMPLE_TYPE_DETAILED`, `PDL1_POSITIVE_ANY`, `PDL1_TESTING`, `CPT_SEQ_DATE` | direct copies |
| `cpt_order_int`, `cpt_seq_date`, `dx_cpt_rep_mos`, `dx_cpt_rep_days`, `cpt_oncotree_code`, `cpt_seq_assay_id`, `mutations` (panel), `cna` (panel) | direct copies |

---

## 11. Non-index cancer summary (per patient)

If `cancer_level_dataset_non_index.csv` (or equivalent) exists, build per-patient
summary and broadcast to all that patient's sample rows:

| Output col | Derivation |
|---|---|
| `non_idx_n_cancers` | `groupby(record_id).size()` from non-index file; fill 0 if patient absent |
| `non_idx_any_heme` | any row where `ca_heme_malig == "Yes"` OR `ca_heme_type != ""` |
| `non_idx_any_brain` | any row where `ca_d_site` or `ca_histology` regex-matches `brain` |
| `non_idx_types` | semicolon-joined sorted unique `ca_type` strings |
| `had_prior_non_breast_cancer` | `non_idx_n_cancers > 0` (bool) |

---

## 12. Mutation features (per-sample summary from MAF)

Input MAF (cBioPortal columns): `sampleId` (or `Tumor_Sample_Barcode`),
`hugoGeneSymbol` (or `Hugo_Symbol`), `Variant_Classification`, `Mutation_Status`,
`t_alt_count`.

### Pre-filter (apply BEFORE aggregating)
Drop rows where `Variant_Classification` is in:
```
{"Silent", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "Intron", "RNA", "IGR"}
```
Drop rows where `hugoGeneSymbol` or `sampleId` is NA.

If `gnomAD_AF` column exists, also drop rows with `gnomAD_AF >= 0.001` (germline filter).

### Per-sample aggregation

| Output col | Derivation |
|---|---|
| `mutation_count_all_sites_sum` | `groupby(sampleId).size()` after filter. Fill 0 for samples not in MAF |
| `t_alt_count_max` | `groupby(sampleId)["t_alt_count"].max()` (numeric coerce) |
| `genes` | `groupby(sampleId).agg(lambda s: ";".join(sorted(s.unique())))` |
| `mutation_count_q` | `pd.qcut(mutation_count_all_sites_sum, q=4, labels=["Q1","Q2","Q3","Q4"])` |
| `t_alt_count_q` | `pd.qcut(t_alt_count_max, q=4, labels=["Q1","Q2","Q3","Q4"])` |

---

## 13. Gene-binary matrix (per-sample x per-gene)

Use gnomeR (R) or equivalent. Settings to match this cohort:

```r
gnomeR::create_gene_binary(
  samples             = all_sample_ids_from_clinical,
  mutation            = maf,                  # cBioPortal long format
  cna                 = cna_long,             # long: sampleId, hugoGeneSymbol, alteration
  fusion              = sv,                   # cBioPortal SV long format
  mut_type            = "somatic_only",
  include_silent      = FALSE,
  snp_only            = FALSE,
  high_level_cna_only = TRUE,                 # only +2 (Amp) and -2 (Del)
  specify_panel       = "no",                 # multi-panel cohort
  recode_aliases      = "no"
)
```

Output columns:
- bare HUGO symbol -> mutation 0/1 (e.g. `TP53`, `PIK3CA`)
- `<GENE>.Amp` -> 0/1 (high-level amplification)
- `<GENE>.Del` -> 0/1 (deep deletion)
- `<GENE>.fus` -> 0/1 (any SV involving the gene)

DROP columns where the column sum across all samples is 0 (gnomeR does this by default).

Samples missing from MAF/CNA/SV are kept with all-zero rows (use the `samples = ` arg).

---

## 14. Pathway annotation

Apply `gnomeR::add_pathways()` with the 10 default Sanchez-Vega pathways. Output 10 binary
columns:
```
pathway_RTK/RAS, pathway_Nrf2, pathway_PI3K, pathway_TGFB, pathway_p53,
pathway_Wnt, pathway_Myc, pathway_Cell cycle, pathway_Hippo, pathway_Notch
```
A pathway is 1 if ANY gene in that pathway has a mutation, .Amp, .Del, or .fus event.

---

## 15. Top-gene pipeline (analytic indicator)

1. Restrict to `any_brain_met == 1` samples (the brain-met cohort).
2. For each mutation gene col (bare HUGO), compute `n_brain_met_samples_mutated = col.sum()`.
3. Sort descending; take top 10 -> `TOP10_GENES`.
4. Take first 5 of TOP10 -> `TOP5_GENES`.
5. For each gene `G` in TOP10, create `G_top10_<G>` = direct copy of gene col `G` (int 0/1).
6. Compute:
   - `top5_any_mutated = (TOP5 cols sum > 0)` int 0/1
   - `top10_any_mutated = (TOP10 cols sum > 0)` int 0/1
   - `top5_n_mutated = TOP5 cols sum` int
   - `top10_n_mutated = TOP10 cols sum` int

This list of 10 genes is COHORT-SPECIFIC. Save the top-10 gene list as a separate file
(`<study>_top_genes.txt`) so it can be reused/audited.

NOTE: Do not hard-code TP53/PIK3CA/etc. from the GENIE cohort; derive top-10 fresh per
new cohort.

---

## 16. Required outputs

Save 4 files in the same directory:

1. `extracted_variables_<study>_data.csv` - the harmonized analytic frame
2. `extracted_variables_<study>_top_genes.txt` - one gene per line, top 10
3. `extracted_variables_<study>_gene_prev_brain_met.csv` - full gene prevalence ranking with cols `gene, n_brain_met_samples_mutated, pct_brain_met, n_total_samples_mutated, pct_total`
4. `extracted_variables_<study>_dictionary.csv` - recoding audit trail, cols `variable, original, mapped, note`

---

## 17. Validation checks (LLM should print these after building)

```
TOTAL ROWS: N_samples
TOTAL COLS: N_cols (expect ~250 clinical + N_gene_binary + 10 pathway + ~10 top-gene)

Cohort split:
  any_brain_met == 1: N samples (M unique patients)
  brain_met_at_dx == 1: N samples
  no brain met:        N samples

Receptor subtype x brain met (sample counts):
  HR+/HER2-:        ...
  HR+/HER2+:        ...
  HR-/HER2+:        ...
  Triple Negative:  ...

Top 10 mutated genes in any_brain_met cohort (rate vs whole cohort):
  GENE  brain_met_pct  total_pct
  ...
```

---

## 18. Column-name compatibility notes

cBioPortal uses inconsistent casing across files. The output spec uses these canonical names.
If your source file uses a synonym, map it before joining:

| Canonical (this spec) | Synonyms |
|---|---|
| `SAMPLE_ID` | `Tumor_Sample_Barcode`, `Sample_Id`, `sampleId`, `cpt_genie_sample_id` |
| `record_id` | `PATIENT_ID`, `Patient_Id`, `patientId` |
| `hugoGeneSymbol` | `Hugo_Symbol`, `hugo_symbol`, `gene_name`, `Gene_Symbol` |
| `Tumor_Sample_Barcode` (in MAF) | `sampleId` (cBioPortal API) |

When reading cBioPortal `.txt` clinical files, skip rows starting with `#` (4 commented
header lines before the actual column header).

---

## 19. Modeling-ready endpoints (final reminder)

These are the time-to-event variables for survival modeling. Pair each `*_months` with its event indicator. Three primary endpoints for the brain-met analysis: **time to brain met**, **OS**, **PFS**.

| Time | Event | Filter | Use for |
|---|---|---|---|
| `tt_brain_met_mos` | `brain_met_event` | none | brain-met survival analysis (primary endpoint) |
| `OS_months` (= `tt_os_dx_mos`) | `os_status_bin` | none | overall survival from dx |
| `PFS_imaging_months` (= `tt_pfs_i_adv_mos`) | `pfs_i_event_bin` | advanced-disease cohort only | PFS imaging |
| `PFS_medonc_months` (= `tt_pfs_m_adv_mos`) | `pfs_m_event_bin` | advanced-disease cohort only | PFS medonc |

Stratification candidates for KM / Cox PH / AFT models:
- `top5_any_mutated`, `top10_any_mutated`
- `receptor_primary_cat`
- `bca_subtype`
- `stage_iv_bin`, `stage_diag_group`
- `grade_ord`
- `age_cat`
- Individual `pathway_<X>` cols
- `met_loc`

Continuous covariates: `age_dx_num`, `mutation_count_all_sites_sum`, `t_alt_count_max`.
