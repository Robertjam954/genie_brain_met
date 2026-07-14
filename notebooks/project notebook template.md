---
# Dataset metadata (YAML front-matter) - fill these fields for each dataset
dataset_name: "<short dataset name>"
dataset_id: "<unique id or slug>"        # e.g. studyname_v1
dataset_version: "<version>"             # semantic or numeric version
description: |
  <One-paragraph description of the dataset: what it contains and purpose>
source: "<origin or provider>"            # e.g. internal database, SEER, public API
source_url: "<url or accession id>"
date_acquired: "YYYY-MM-DD"
license: "<license or usage terms>"
access_restrictions: "<open/restricted/controlled>"
sensitive: false                          # true if contains PHI or sensitive data
data_root: "data/<raw|external|processed>" # recommended top-level storage
files:
  - path: "<relative/path/to/file.ext>"
    filename: "<filename.ext>"
    format: "<format>"                   # csv, parquet, xlsx, json, etc.
    size_bytes: <int>
    checksum: "<sha256|md5>:<value>"
    rows: <int>                            # optional, number of rows/records
schema: |
  - name: <column_name>
    type: <data_type>                     # integer, float, string, datetime, categorical
    description: "<short description>"
    nullable: true/false
key_identifiers: ["<id_col1>", "<id_col2>"]
primary_key: "<primary_key_column>"
sample_size:
  observations: <int>
  units: "<rows/patients/samples>"
date_range:
  start: "YYYY-MM-DD"
  end: "YYYY-MM-DD"
update_frequency: "<one-time|daily|weekly|monthly|as-needed>"
preprocessing_notes: |
  - <brief steps applied to raw data: cleaning, encoding, derivations>
missing_value_strategy: "<how missingness handled>"
transformations: |
  - <transformation 1: e.g. log-transform of column X>
provenance: |
  - created_by: "<script or notebook path>"
  - created_on: "YYYY-MM-DD"
  - commit: "<git sha or tag>"                # link data generation to code version
contact:
  name: "<owner name>"
  email: "<owner@example.org>"
  role: "<data steward / analyst>"
manifest_file: "<path/to/manifest.json|csv>"
notes: |
  <freeform notes, include known data quality issues, known biases, citations>
tags: ["<tag1>", "<tag2>"]
---

# Dataset README / Summary (human-readable)

## Quick summary
- Name: **<dataset_name>**
- Version: **<dataset_version>**
- Source: **<source>**
- Acquired: **<date_acquired>**
- Size: **<observations>** rows

## Files included
- `path/to/file.ext` — <one-line description>

## Schema (high level)
- `<column_name>` — <type> — <short description>

## Key fields and identifiers
- Primary key: `<primary_key>`
- ID columns: `<id_col1>, <id_col2>`

## Preprocessing and transformations
- Steps applied: (cleaning, imputation, encoding, feature creation)

## Missingness and quality notes
- Missing-rate per key variable: `<var1>: 3%, var2: 0%`
- Known issues: `<describe>`

## Provenance and reproducibility
- Script/notebook used to create dataset: `<path>`
- Repo commit: `<git-sha>`
- Manifest: `<manifest_file>`

## Access and licensing
- License: `<license>`
- Access restrictions: `<access_restrictions>`

## Contact
- `<owner name> (<owner@example.org>)`

---

# How to use this template
1. Copy this file into the dataset folder or repository root and fill the YAML front-matter.
2. Keep the `manifest_file` and `checksum` up to date when files change.
3. Link the `commit` field to the code that produces the dataset so others can reproduce it.

# Minimal example
---
dataset_name: "brain_mets_clinical"
dataset_id: "brain_mets_clinical_v1"
dataset_version: "1.0"
description: |
  Clinical and genomic data for brain metastasis patients collected from internal registry.
source: "Internal registry"
source_url: ""
date_acquired: "2024-12-15"
license: "internal-research-use"
access_restrictions: "controlled"
sensitive: true
data_root: "data/raw"
files:
  - path: "data/raw/clinical.csv"
    filename: "clinical.csv"
    format: "csv"
    size_bytes: 1234567
    checksum: "sha256:abcdef..."
    rows: 1024
schema: |
  - name: patient_id
    type: string
    description: "de-identified patient identifier"
    nullable: false
  - name: age
    type: integer
    description: "age at diagnosis"
    nullable: false
key_identifiers: ["patient_id"]
primary_key: "patient_id"
sample_size:
  observations: 1024
  units: "patients"
date_range:
  start: "2010-01-01"
  end: "2023-12-31"
update_frequency: "as-needed"
preprocessing_notes: |
  - Drop duplicate patient rows keeping first
  - Encode categorical variables with label encoding
missing_value_strategy: "impute_median_for_numeric, mode_for_categorical"
provenance: |
  - created_by: "src/data collection and processing/ingest_clinical.py"
  - created_on: "2024-12-15"
  - commit: "abc123def"
contact:
  name: "Data Owner"
  email: "owner@example.org"
  role: "data steward"
manifest_file: "manifests/brain_mets_manifest.csv"
notes: |
  - Sensitive clinical variables present; dataset access controlled.
  - See manuscript for cohort definition.
tags: ["clinical", "genomics", "brain-metastasis"]
---
