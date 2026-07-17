# Known Issues and Gaps

This page tracks current issues, broken paths, and remaining work. See also `STATUS.md`
in the repo root for the live task checklist.

---

## Critical gaps (pipeline will not run end-to-end without addressing these)

### 1. ETL has not been executed in this checkout

`data/processed/` is absent. All Aim scripts depend on the harmonized analytic frames
being present there. The harmonize scripts write to `src/exploratory data analysis/`,
not `data/processed/`, so the two locations must be reconciled manually.

**Fix:** Run harmonize scripts, then copy outputs to `data/processed/`. See
[Pipeline Run Guide - Stage 4](Pipeline-Run-Guide.md).

### 2. Machine-specific absolute paths in R scripts and notebooks

The R pull scripts (`pull_genie_data.R`, `extract_variables_of_interest.R`) reference
`/Users/robertjames/loc/data private/...`. The notebooks hard-code a `PROJ` root under
the author's iCloud `Documents` folder.

**Fix:** Update these paths to your local data root before running.

### 3. `data/processed/` not present in the repository

The directory is git-ignored and not created by any script in this repo. It must be
created manually and populated with the harmonized frames.

### 4. No pinned environment spec

There is no `requirements.txt`, `environment.yml`, or `renv.lock`. Package versions
may drift and break the pipeline.

**Fix:** After installing dependencies that work, freeze them:
```zsh
pip freeze > requirements.txt
```

---

## Known bugs and inconsistencies

### `add_pathways_genie_bpc.R` is an empty stub

`src/data collection and processing/add_pathways_genie_bpc.R` is a 0-byte file. It does
nothing. Pathway columns (`pathway_RTK/RAS`, etc.) are attached during the harmonization
step or in upstream master files, not by this script.

### `README_survival_AFT_pipeline.md` describes a different pipeline variant

`src/modeling/README_survival_AFT_pipeline.md` documents a standalone
`survival_and_xgb_analysis.py` CLI with `--xlsx`, `--time-col`, `--event-col`,
`--strata`, `--topn-genes` arguments and column names (`DFS_MONTHS`, `SUBTYPE`,
`MANTIS_BIN`, `G__*`) that do not match the current `_lib`-based canonical schema.
Reconcile before relying on it.

### Harmonize scripts write to wrong directory for modeling

The harmonize scripts write outputs to `src/exploratory data analysis/` but
`_lib.load_cohort()` reads from `data/processed/`. This disconnect requires a manual
copy step after harmonization.

### `create_deidentified_dataset.R` references offline Windows paths

This script (`src/data collection and processing/create_deidentified_dataset.R`)
contains hardcoded Windows paths. It must be updated for any other machine or OS.

---

## Completeness of committed results

Only a subset of expected outputs are currently committed:

| Aim | Cohort | Status |
|-----|--------|--------|
| Aim 1 | GENIE | Result tables committed under `src/data reports/genie/aim1/` |
| Aim 1 | TCGA | Not committed |
| Aim 1 | MSK 2018 | Not committed |
| Aim 2 | All cohorts | Not committed |
| Aim 3 | All cohorts | Not committed |
| ML benchmark | All cohorts | Figures in `reports/figures/ml_benchmark/` (partially committed) |

---

## Remaining work (from STATUS.md)

- [ ] Run the ETL: harmonize scripts never executed; `data/processed/` empty
- [ ] Implement `add_pathways_genie_bpc.R` (currently a 0-byte stub)
- [ ] Missingness audit finalized across the three cohorts
- [ ] Descriptive/EDA parity for TCGA + MSK cohorts
- [ ] Generate missing outputs: Aim 1 for TCGA/MSK, all Aim 2 & 3 results
- [ ] Validation across replication cohorts
- [ ] Add `requirements.txt` / `environment.yml` with pinned versions
- [ ] Publication figures finalized in `reports/figures/`
- [ ] Reproducibility check: fresh clone runs the pipeline end-to-end

---

## Not part of the analytic pipeline

The top-level `archive/` directory (formerly `self-documenting-ai-agent/`,
`claude-md-memory-workflow/`) and the `context-engineering-workflow.md`,
`plan-template.md`, and `.claude/` files are AI-engineering workflow scaffolding and
are unrelated to the brain-metastasis analysis.
