# Status

Live checklist of what is left to complete on this project, in the same checkbox
format the portfolio tracking dashboard reads. Check a box (`[ ]` -> `[x]`) as
you finish.

> Project: **GENIE BPC BrCa Brain-Metastasis Analysis** · Stage: **data / modeling**
> Cohorts: GENIE BPC (primary), TCGA-BRCA + MSK-IMPACT 2018 (replication).

## 1. Scoping
- [x] Aims written with per-aim cohort restriction (see README §1)
- [x] Primary outcome defined (CNS metastasis ever)

## 2. Data collection & processing
- [ ] Run the ETL: harmonize scripts never executed; `data/processed/` empty (all Aims depend on it)
- [ ] Implement `src/data collection and processing/add_pathways_genie_bpc.R` (currently a 0-byte stub)
- [x] Raw GENIE BPC release landed in `data/raw/` (git-ignored)
- [ ] Missingness audit finalized across the three cohorts

## 3. Exploratory data analysis
- [x] Aim 1 gene-prevalence tables + oncoprint/forest (GENIE)
- [ ] Descriptive/EDA parity for TCGA + MSK cohorts

## 4. Modeling
- [ ] Generate missing outputs: Aim 1 for TCGA+MSK, and all Aim 2 & 3 survival / XGBoost-AFT results
- [ ] Validation across replication cohorts

## 5. Reporting & repo hygiene
- [ ] Add requirements.txt / environment.yml with pinned versions
- [ ] Flesh out root README with setup + pipeline run order
- [x] Move/archive AI-workflow tutorial folders cluttering the research repo (-> `archive/`)
- [ ] Publication figures finalized in `reports/figures/`
- [ ] Reproducibility check: fresh clone runs the pipeline end-to-end
