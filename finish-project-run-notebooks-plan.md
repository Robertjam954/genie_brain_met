---
title: Finish the brain-metastasis project by executing the analysis notebooks
version: 1.0
date_created: 2026-07-13
last_updated: 2026-07-13
---
# Implementation Plan: Execute the GENIE/TCGA/MSK brain-met analysis notebooks end-to-end

Bring the project to a "finished" state by running the analysis notebooks to
completion against the real processed data, producing the Aim 1/2/3 tables,
figures, and the ML survival benchmark. The notebooks are the deliverable; they
read cached analytic CSVs and emit manuscript components.

## Architecture and design

### Where things live (verified)
- **Notebooks** (the copy to run): `notebooks/*.ipynb` in this repo.
- **Data + outputs root (`PROJ`)**: notebooks hard-code
  `PROJ = /Users/robertjames/Documents/Documents - Robert's iMac/Research Projects/MSKCC Research Fellowship/Projects/genie_tcga_impact_brain_mets`.
  They read `PROJ/data/processed/*` and **write** to `PROJ/data/interim/`,
  `PROJ/manuscript components/figures/`, and `PROJ/src/modeling/genie/aimN/`.
  They `sys.path.insert(PROJ/src/modeling)` and import the **iCloud** `_lib.py`
  (confirmed present), not this repo's copy.
- **Data**: `PROJ/data/processed/` holds the full set — `extracted_variables_<cohort>_data.csv`
  (genie 7 MB, tcga 124 MB, msk 6 MB), `*_top_genes.txt`, `*_gene_prev_brain_met.csv`,
  the `genie_bpc_v1_*` masters, and `complete_*` files. Files are readable (materialized).

### Blockers found
1. **No Python environment has the survival stack.** Every conda env is missing
   `lifelines, statsmodels, shap, seaborn, openpyxl, lightgbm, optuna, sksurv`.
   Full third-party import set across the notebooks:
   `pandas numpy scipy scikit-learn matplotlib seaborn statsmodels lifelines
   xgboost lightgbm optuna shap scikit-survival(sksurv) openpyxl`.
2. **The registered `python3` Jupyter kernel is dead** — it points at
   `~/Documents/llm_summarization/.venv/bin/python`, which no longer exists.
3. **`retrieval_metadata_recode.ipynb` references missing raw data**
   (`/Users/robertjames/loc/data private/genie_bpc_datav1/`, absent). It is a
   documentation-only notebook ("does not write data") and is expected to fail
   partway; it is not required to produce Aim outputs.
4. **iCloud eviction risk**: large files under `~/Documents` can be dataless;
   materialize the big CSVs (esp. tcga 124 MB) with `dd .. of=/dev/null` before pandas reads.

### Environment strategy
Extend the existing `miniforge3/envs/tcga-analysis` env (already has
pandas/numpy/scipy/sklearn/matplotlib/xgboost/ipykernel/nbconvert) by adding the
missing packages, then register it as a dedicated Jupyter kernel
`genie-brainmet`. Execute each notebook with `jupyter nbconvert --to notebook
--execute` pinned to that kernel. Per project convention, pin scientific stack
versions compatible with the existing env; do not upgrade base.

### Execution order (dependency-aware)
1. `descriptiveplots.ipynb`      - Table 1 / descriptive (Aim 0)
2. `timebrainmettplots.ipynb`    - Aim 2, time to brain met
3. `osplots.ipynb`               - Aim 3, overall survival
4. `coxvsaft_error_os.ipynb`     - ML survival benchmark (XGB-AFT + RSF/GBSA/LGBM)
5. `src/modeling/survival_and_xgb_analysis.ipynb` - combined modeling (if in scope)
- `retrieval_metadata_recode.ipynb` - run **last / optional**; documentation-only, likely errors on missing raw data.

Notebooks are "final reflections" of the `src/modeling` scripts reading cached
CSVs, so they are largely independent; order above is for logical review, not hard coupling.

## Tasks

- [ ] **Env**: add `statsmodels seaborn openpyxl lifelines lightgbm optuna shap scikit-survival` to `tcga-analysis`.
- [ ] **Kernel**: `python -m ipykernel install --user --name genie-brainmet`.
- [ ] **Materialize** large processed CSVs under `PROJ/data/processed` (dd to /dev/null).
- [ ] **Smoke test**: import `_lib`, `load_cohort('genie')`, confirm shape, before full runs.
- [ ] **Execute** notebooks 1-4 in order via nbconvert (`--execute`, save executed copy).
- [ ] Optionally execute `survival_and_xgb_analysis.ipynb` and `retrieval_metadata_recode.ipynb`.
- [ ] **Verify** each: no unhandled exception; expected figures/CSVs written under `PROJ`.
- [ ] **Report** per-notebook status (cells run, outputs written, any skipped cells).

## Open questions
1. **Output location**: notebooks write into the iCloud research project
   (`PROJ`), not this repo. Proceed writing there, or redirect outputs into the repo?
2. **Scope**: run all notebooks, or skip the doc-only `retrieval_metadata_recode.ipynb`
   (and the `src/modeling` combined notebook)?
3. **Env**: extend `tcga-analysis` in place, or build a fresh isolated env?
