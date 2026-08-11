# Analysis order of operations

The eight-step protocol for the brain-metastasis analysis, and the scripts that implement
each step. Steps 1-4 run in R (`survival` only, no tidyverse dependency); steps 5-7 run in
Python. Every stage writes CSV tables under `src/modeling/<cohort>/<stage>/` and figures
under `manuscript components/<cohort>/<stage>/`.

This protocol is a complement to, not a replacement for, the Aim 1/2/3 scripts described in
[`../README.md`](../README.md#pipeline-run-order): it adds Table 1, a PRISMA-style cohort
flow chart, a competing-risk (Fine-Gray) brain-metastasis model, risk-model-driven feature
selection, and a Random Survival Forest for time to brain metastasis.

## Steps and implementing files

| Step | What | File |
|------|------|------|
| 1 | Data import - read the analytic frame (`.csv` or `.xlsx`): clinical covariates, mutation indicators, durations, vital status | `src/modeling/_lib.R` (`read_any`) |
| 2 | Preprocessing - categoricals to factors, `survival_time` / `time_to_brain_mets`, binary `event_dead` / `brain_mets_flag` | `src/modeling/_lib.R` (`preprocess`) |
| 3 | Descriptive analysis - Table 1 (overall and by brain-met status) and Figure 1 (PRISMA-style cohort flow chart) | `src/exploratory data analysis/table1_prisma_descriptive.R` |
| 4A | Brain-metastasis risk - cumulative incidence with death as the competing risk, multivariable Fine-Gray subdistribution hazards, per-gene adjusted subhazard ratios with BH q-values | `src/modeling/finegray_cox_risk_models.R` |
| 4B | Overall survival - Kaplan-Meier by receptor subtype with a log-rank test, multivariable Cox PH (with scaled Schoenfeld diagnostics), per-gene adjusted hazard ratios | `src/modeling/finegray_cox_risk_models.R` |
| 5 | Feature selection - genes with `q < 0.05` in Fine-Gray, or `HR > 1` and `p < 0.05` in Cox PH | `src/modeling/select_risk_model_genes.py` |
| 6 | Predictive modeling of time to brain met - one-hot encoding, Random Survival Forest, 5-fold CV tuning, C-index / time-dependent AUC / Brier score | `src/modeling/rsf_time_to_brain_met.py` |
| 7 | Interpretation - permutation variable importance and partial dependence for the top predictors | `src/modeling/rsf_time_to_brain_met.py` |
| 8 | Reproducibility - this document, `requirements.txt`, `environment.yml`, and the R package list in [`../README.md`](../README.md#r) | - |

## Run order

```zsh
# Steps 1-3: import, preprocess, Table 1 + PRISMA flow chart
Rscript "src/exploratory data analysis/table1_prisma_descriptive.R" --cohort genie

# Step 4: Fine-Gray competing-risk + Cox PH risk models
Rscript "src/modeling/finegray_cox_risk_models.R" --cohort genie

# Step 5: gene selection from the step-4 tables
python "src/modeling/select_risk_model_genes.py" --cohort genie

# Steps 6-7: Random Survival Forest + interpretation
python "src/modeling/rsf_time_to_brain_met.py" --cohort genie
```

All four accept `--cohort {genie,tcga,msk18}`. The R scripts additionally accept
`--data PATH` (a `.csv`, `.tsv`, or `.xlsx` analytic frame), `--sheet`, `--outdir`,
`--figdir`, `--min-events`, and `--strata`; the Python scripts accept `--risk-dir`,
`--data`, `--genes-file`, `--outdir`, `--figdir`, `--test-size`, `--seed`,
`--n-estimators`, and `--no-tune`. With no `--data`, the frame is read from
`data/processed/extracted_variables_<cohort>_data.csv`, the same location
`_lib.load_cohort()` uses - see the prerequisite note in
[`../README.md`](../README.md#data).

## Inputs

Canonical columns consumed from the harmonized analytic frame (only the columns present in
a cohort are used, so a cohort missing e.g. `insurance` still runs):

- Durations / status: `OS_months` (or `tt_os_dx_mos`), `os_status_bin`, `tt_brain_met_mos`
  (or `time_to_brain_met_mos`), `brain_met_event`, `any_brain_met`, `brain_met_at_dx`
- Clinical covariates: `age_dx_num`, `sex`, `smoking_status`, `stage_dx_cat`,
  `receptor_primary_cat`, `insurance`, `grade_ord`, `race_clean`, `ethnicity_clean`,
  `mutation_count`, `SEQ_ASSAY_ID`
- Gene alterations: the `G_top10_*` indicator columns

## Outputs

### Step 3 - `src/modeling/<cohort>/descriptive/`

| File | Contents |
|------|----------|
| `table1_overall.csv` | Table 1, whole cohort (median [IQR] for numerics, n (%) per level for factors) |
| `table1_by_brain_met.csv` | Table 1 stratified by brain-met status, with Wilcoxon / chi-square (Fisher when expected counts < 5) p-values |
| `cohort_flow_counts.csv` | the counts behind Figure 1 |
| `preprocessed_summary.csv` | audit of the derived durations and event flags |
| `manuscript components/<cohort>/descriptive/figure1_prisma_flowchart.png` | Figure 1 |

### Step 4 - `src/modeling/<cohort>/risk_models/`

| File | Contents |
|------|----------|
| `cif_brain_met_overall.csv` | Aalen-Johansen cumulative incidence of brain met and of death |
| `cif_brain_met_by_receptor_primary_cat.csv` | the same, per receptor subtype |
| `finegray_multivariable.csv` | clinical-covariate subhazard ratios with 95% CIs |
| `finegray_gene_subhazards.csv` | per-gene sHR, 95% CI, p, BH q (each gene adjusted for the clinical covariates) |
| `km_os_by_receptor_primary_cat.csv`, `logrank_os_by_receptor_primary_cat.txt` | KM estimates and the log-rank test |
| `cox_os_multivariable.csv` | Cox PH hazard ratios with 95% CIs and per-term Schoenfeld p-values |
| `cox_os_gene_hazards.csv` | per-gene HR, 95% CI, p, BH q |
| `risk_model_summary.txt` | cohort sizes, event counts, global PH test, model-skip reasons |
| `manuscript components/<cohort>/risk_models/figure2_cif_brain_met.png` | Figure 2, cumulative incidence by subtype |
| `manuscript components/<cohort>/risk_models/figure3_km_os.png` | Figure 3, KM overall survival by subtype |

### Step 5 - `src/modeling/<cohort>/risk_models/`

`selected_genes.csv` (one row per gene, recording which criterion kept it and the
underlying estimates) and `selected_genes.txt` (the bare list consumed by step 6).

### Steps 6-7 - `src/modeling/<cohort>/rsf/`

| File | Contents |
|------|----------|
| `rsf_metrics.json` | best hyperparameters, CV C-index, test Harrell and IPCW C-index, mean time-dependent AUC, integrated Brier score |
| `rsf_cv_results.csv` | the full 5-fold grid-search table |
| `rsf_time_dependent_auc.csv` | cumulative/dynamic AUC per evaluation time |
| `rsf_predicted_survival.csv` | mean predicted survival at each evaluation time |
| `rsf_variable_importance.csv` | permutation importance (drop in C-index), mean and SD over 10 repeats |
| `rsf_partial_dependence.csv` | marginal mean predicted risk over each top predictor's grid |
| `manuscript components/<cohort>/rsf/*.png` | time-dependent AUC, variable importance, partial dependence |

## Modeling notes

- **Fine-Gray without extra packages.** `survival::finegray()` builds the risk-weighted data
  set and `coxph(..., weights = fgwt)` fits the subdistribution hazards, so `survival` (a
  base-recommended package) is the only R dependency. `cmprsk::crr` is an equivalent route.
- **Competing risk = death.** Step 4A builds a three-level status (0 censored / 1 brain met /
  2 death without brain met) at the first observed event time, so death no longer inflates
  the brain-metastasis incidence the way naive Kaplan-Meier does.
- **Death is censoring in step 6.** A Random Survival Forest models a single right-censored
  endpoint, so the RSF treats death as censoring; the competing-risk interpretation stays
  with the Fine-Gray model. Read the RSF as a discrimination benchmark, not as a
  competing-risk estimator.
- **Per-gene models, not one saturated model.** Each gene enters its own model alongside the
  clinical covariates. With ~10 gene indicators and cohort-scale event counts, a single
  model containing every gene is unstable; per-gene models with BH correction across genes
  are the standard compromise, and they are what the step-5 rule expects.
- **Small-cohort guards.** Fine-Gray and Cox PH are skipped below `--min-events` (default 25)
  events, the RSF below 25 events, and cumulative-incidence strata below 10 patients or 3
  events. Skips are recorded in `risk_model_summary.txt` / `rsf_skipped.txt` rather than
  failing the run. Sparse factor levels (n < 5) are collapsed into `Other` instead of the
  covariate being dropped.
- **Evaluation-time grid.** Time-dependent AUC and the Brier score are evaluated on 12 times
  between the 10th and 90th percentile of the test-set event times, clipped to the training
  follow-up so the IPCW censoring weights stay defined.

## Reproducibility

- Python: `requirements.txt` / `environment.yml`. Steps 6-7 need `scikit-survival` (pinned
  at 0.22.2), `scikit-learn`, `pandas`, `numpy`, and `matplotlib`.
- R: `survival` (ships with R) is sufficient. `readxl` or `openxlsx` is needed only for
  `.xlsx` input.
- Seeds: the RSF split, forest, and permutation importance all take `--seed` (default 42).
- Validation to date: these four scripts have been exercised end-to-end on a synthetic
  GENIE-shaped frame (600 patients, planted `PIK3CA` / `TP53` brain-met signal, which the
  Fine-Gray model, the step-5 selection, and the RSF importance ranking all recover), plus a
  60-patient frame with missing covariate columns to check the skip guards. They have **not**
  been run on the real cohorts here, because the analytic frames are git-ignored and were not
  available.
