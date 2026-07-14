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
| Task | Description | Example Tool |
|------|-------------|--------------|
| Replace or drop NA fields | Handle missing values: impute, drop a row, or drop a column | Pandas |
| Remove duplicates | Drop redundant rows | SQL DISTINCT |
| Feature scaling | Normalize numerical data | Scikit-learn |
| Encoding | Categorize nominal and ordinal factors into appropriate unordered and ordered factor levels respectively | Pandas, Scikit-learn |
| Output | Clean, well-structured dataset ready for analysis |  |


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
def fetch_molecular_data(entrez_gene_ids=None, molecular_profile_ids=None, sample_molecular_identifiers=None, projection="SUMMARY"):
    """
    Fetch molecular data. \n
    :param entrez_gene_ids: List of Entrez Gene IDs (e.g., ["672", "675"]). \n
    :type entrez_gene_ids: list of str \n
    :param molecular_profile_ids: List of MolecularProfile IDs (e.g., ["brca_tcga_mrna", "acc_tcga_rna_seq_v2_mrna"]). \n
    :type molecular_profile_ids: list of str \n
    :param sample_molecular_identifiers: List of Molecular Profile ID and Sample ID pairs. \n
    :type sample_molecular_identifiers: list of dict \n
        Each dict should have the following format: \n
            sample_molecular_identifiers = [
                                           {"molecular_profile_id": "brca_tcga_mrna", 
                                            "sample_ids": ["TCGA-AR-A1AR-01","TCGA-BH-A1EO-01"]},
                                           {"molecular_profile_id": "acc_tcga_rna_seq_v2_mrna", 
                                            "sample_ids": ["TCGA-OR-A5J1-01","TCGA-OR-A5J2"]}
                                           ]
    :param projection: Level of detail of the response. \n
        Possible values: \n
            - "DETAILED": Detailed information.
            - "ID": Information with only IDs.
            - "META": Metadata information.
            - "SUMMARY": Summary information (default).
    :type projection: str \n
    :returns: A DataFrame containing molecular data. \n
    :rtype: pandas.DataFrame \n
    """
    endpoint = "/molecular-data/fetch"
    params = {"projection": projection}

    molecular_data_filter = {}

    if entrez_gene_ids:
        molecular_data_filter['entrezGeneIds'] = entrez_gene_ids
    
    if molecular_profile_ids:
        molecular_data_filter['molecularProfileIds'] = molecular_profile_ids

    if sample_molecular_identifiers:
        molecular_data_filter['sampleMolecularIdentifiers'] = []

        for item in sample_molecular_identifiers:
            molec_prof_id = item["molecular_profile_id"]
            sample_ids = item["sample_ids"]

            for sample_id in sample_ids:
                identifier = {
                    "molecularProfileId": molec_prof_id,
                    "sampleId": sample_id
                }
                molecular_data_filter["sampleMolecularIdentifiers"].append(identifier)

    response = requests.post(f"{base_url}{endpoint}", params=params, json=molecular_data_filter)
    return process_response(response, "Failed to fetch molecular data.")
        
def get_all_molecular_data_in_molecular_profile(molecular_profile_id, sample_list_id, entrez_gene_id, projection="SUMMARY"):
    """
    Get all molecular data in a molecular profile for a specific gene. \n
    :param molecular_profile_id: Molecular Profile ID (e.g., "acc_tcga_rna_seq_v2_mrna"). \n
    :type molecular_profile_id: str \n
    :param sample_list_id: Sample List ID (e.g., "acc_tcga_all"). \n
    :type sample_list_id: str \n
    :param entrez_gene_id: Entrez Gene ID (e.g., "1"). \n
    :type entrez_gene_id: str \n
    :param projection: Level of detail of the response. \n
        Possible values: \n
            - "DETAILED": Detailed information.
            - "ID": Information with only IDs.
            - "META": Metadata information.
            - "SUMMARY": Summary information (default).
    :type projection: str \n
    :returns: A DataFrame containing molecular data for the specified gene. \n
    :rtype: pandas.DataFrame \n
    """
    endpoint = f"/molecular-profiles/{molecular_profile_id}/molecular-data"
    params = {
        "entrezGeneId": entrez_gene_id,
        "projection": projection,
        "sampleListId": sample_list_id
    }

    response = requests.get(f"{base_url}{endpoint}", params=params)
    return process_response(response, "Failed to get molecular data in molecular profile.")

def fetch_all_molecular_data_in_molecular_profile(molecular_profile_id, entrez_gene_ids = None, sample_ids = None, 
                                                  sample_list_id = None, projection="SUMMARY"):
    """
    Fetch molecular data in a molecular profile for a list of genes. \n
    :param molecular_profile_id: Molecular Profile ID (e.g., "acc_tcga_rna_seq_v2_mrna"). \n
    :type molecular_profile_id: str \n
    :param entrez_gene_ids: List of Entrez Gene IDs (e.g., ["672","675"]). \n
    :type entrez_gene_ids: list of str \n
    :param sample_ids: List of Sample IDs (e.g., ["TCGA-AR-A1AR-01","TCGA-BH-A1EO-01"]). \n
    :type sample_ids: list of str \n
    :param sample_list_id: Sample List ID (e.g., "brca_tcga_all"). \n
    :type sample_list_id: str \n
    :param projection: Level of detail of the response. \n
        Possible values: \n
            - "DETAILED": Detailed information.
            - "ID": Information with only IDs.
            - "META": Metadata information.
            - "SUMMARY": Summary information (default).
    :type projection: str \n
    :returns: A DataFrame containing molecular data for the specified genes. \n
    :rtype: pandas.DataFrame \n
    """
    endpoint = f"/molecular-profiles/{molecular_profile_id}/molecular-data/fetch"
    params = {
        "projection": projection,
    }

    molecular_data_filter = {}

    if entrez_gene_ids:
        molecular_data_filter['entrezGeneIds'] = entrez_gene_ids
    
    if sample_ids:
        molecular_data_filter['sampleIds'] = sample_ids

    if sample_list_id:
        molecular_data_filter['sampleListId'] = sample_list_id

    response = requests.post(f"{base_url}{endpoint}", params=params, json=molecular_data_filter)
    return process_response(response, "Failed to fetch molecular data in molecular profile.")
