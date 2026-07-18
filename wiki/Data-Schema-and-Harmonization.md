# Data Schema and Harmonization

The canonical analytic frame is defined in
`src/data collection and processing/harmonization_spec.md`. This page summarizes the key
schema elements. For the full specification, read that file directly.

---

## Overview

The pipeline produces one canonical analytic CSV per cohort:

```
extracted_variables_<cohort>_data.csv
```

- One **row per `SAMPLE_ID`**
- All three cohorts share the same column names
- Columns that don't exist in a given source are created as NA (frames must be conformable)
- Recoding decisions are recorded in `extracted_variables_<cohort>_dictionary.csv`

---

## Required output files (per cohort)

| File | Contents |
|------|---------|
| `extracted_variables_<cohort>_data.csv` | Harmonized analytic frame |
| `extracted_variables_<cohort>_top_genes.txt` | Top-10 mutated genes in the brain-met cohort (one per line) |
| `extracted_variables_<cohort>_gene_prev_brain_met.csv` | Gene-prevalence ranking with columns `gene`, `n_brain_met_samples_mutated`, `pct_brain_met`, `n_total_samples_mutated`, `pct_total` |
| `extracted_variables_<cohort>_dictionary.csv` | Recoding audit trail with columns `variable`, `original`, `mapped`, `note` |

---

## Required identifiers

| Column | Description |
|--------|-------------|
| `SAMPLE_ID` | Sample-level barcode; the grain of the frame |
| `record_id` | Patient-level ID |
| `ca_seq` | Cancer sequence number (if available) |
| `cpt_number` | Panel test number (if available) |

---

## Key clinical columns

### Demographics

| Column | Values |
|--------|--------|
| `SEX` | `Female`, `Male`, `Other`, NA |
| `BIRTH_YEAR` | Integer |
| `CENTER` | Institution / site |
| `race_clean` | `White`, `Black`, `Asian`, `Native American`, NA |
| `ethnicity_clean` | `Non-Hispanic`, `Hispanic`, NA |
| `age_cat` | Ordered: `<50 < 50-70 < >70` |

### Breast subtype / receptor

| Column | Values |
|--------|--------|
| `receptor_primary_cat` | `HR+/HER2-`, `HR+/HER2+`, `HR-/HER2+`, `Triple Negative` |
| `bca_subtype` | Raw subtype string |
| `her2_status_bin` | 0/1 |

### Grade and stage

| Column | Values |
|--------|--------|
| `grade_ord` | Ordered: `Low < Intermediate < High` |
| `stage_diag_group` | Ordered: `Stage I < II < III < IV` |
| `stage_iv_bin` | 0/1 |

---

## Brain-met cohort derivation

The primary analytic split is based on these derived columns:

| Column | Derivation |
|--------|-----------|
| `any_brain_met` | 1 if `dist_mets_brain_cns == 1` OR `DMETS_DX_BRAIN == "yes"` |
| `brain_met_at_dx` | 1 if `DMETS_DX_BRAIN == "yes"` |
| `met_loc` | `Brain` / `Other` / `None` (unordered factor) |
| `tt_brain_met_mos` | Months from dx to brain met (or OS time as censoring) |
| `brain_met_event` | Copy of `any_brain_met` |

Cohort restrictions used per aim:

| Aim | Filter |
|-----|--------|
| Aim 1 | Full cohort |
| Aim 2 | `brain_met_at_dx == 0` (no CNS at diagnosis) |
| Aim 3 | `any_brain_met == 1` (brain-mets-ever) |

---

## Survival endpoints

| Time column | Event column | Description |
|------------|-------------|-------------|
| `tt_brain_met_mos` | `brain_met_event` | Time to brain metastasis |
| `OS_months` (= `tt_os_dx_mos`) | `os_status_bin` | Overall survival from diagnosis |
| `PFS_imaging_months` | `pfs_i_event_bin` | PFS imaging (advanced-disease cohort) |
| `PFS_medonc_months` | `pfs_m_event_bin` | PFS medical oncology (advanced-disease cohort) |

---

## Mutation features

Derived from MAF files after filtering:

- Drop silent and non-coding variant classifications: `Silent`, `3'UTR`, `5'UTR`, `3'Flank`, `5'Flank`, `Intron`, `RNA`, `IGR`
- Drop rows with `gnomAD_AF >= 0.001` if the column exists (germline filter)

| Column | Description |
|--------|-------------|
| `mutation_count_all_sites_sum` | Total filtered mutation count per sample |
| `t_alt_count_max` | Maximum alt allele count per sample |
| `mutation_count_q` | Quartile of mutation count (Q1-Q4) |
| `t_alt_count_q` | Quartile of alt count (Q1-Q4) |

---

## Gene-binary matrix

Produced by gnomeR (`create_gene_binary`) with settings:

- `mut_type = "somatic_only"`, `include_silent = FALSE`
- `high_level_cna_only = TRUE` (+2 Amp, -2 Del only)
- Columns: bare HUGO symbol (mutation), `<GENE>.Amp`, `<GENE>.Del`, `<GENE>.fus`
- Samples absent from MAF/CNA/SV get all-zero rows

---

## Pathway columns

Ten Sanchez-Vega oncogenic pathways from `gnomeR::add_pathways()`:

```
pathway_RTK/RAS  pathway_Nrf2  pathway_PI3K  pathway_TGFB  pathway_p53
pathway_Wnt      pathway_Myc   pathway_Cell cycle  pathway_Hippo  pathway_Notch
```

A pathway is 1 if any gene in that pathway has a mutation, Amp, Del, or fusion event.
These are defined in `_lib.PATHWAY_COLS`.

---

## Top-gene pipeline

1. Restrict to `any_brain_met == 1` samples
2. Count per-gene mutation rate in the brain-met cohort
3. Take top 10 genes -> `TOP10_GENES`; first 5 -> `TOP5_GENES`
4. Create `G_top10_<G>` indicator columns for each gene in TOP10
5. Compute `top5_any_mutated`, `top10_any_mutated`, `top5_n_mutated`, `top10_n_mutated`

> The top-10 gene list is **cohort-specific** and derived fresh for each cohort. Do not
> hard-code GENIE-derived gene names for TCGA or MSK.

---

## Column-name synonyms

| Canonical | Synonyms in source files |
|-----------|--------------------------|
| `SAMPLE_ID` | `Tumor_Sample_Barcode`, `Sample_Id`, `sampleId`, `cpt_genie_sample_id` |
| `record_id` | `PATIENT_ID`, `Patient_Id`, `patientId` |
| `hugoGeneSymbol` | `Hugo_Symbol`, `hugo_symbol`, `gene_name`, `Gene_Symbol` |

When reading cBioPortal `.txt` files, skip rows starting with `#` (4 commented header
lines precede the actual column header row).
