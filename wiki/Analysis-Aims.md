# Analysis Aims

The pipeline implements four analytic aims across up to three cohorts. Each aim can be
run independently once the harmonized analytic frames are in place.

> **Status caveat.** Committed result tables and figures are currently limited to GENIE
> Aim 1. TCGA / MSK and Aim 2/3 results are not yet present in the repo.

---

## Cohort restrictions

| Aim | Patient filter |
|-----|---------------|
| Aim 1 | Full cohort |
| Aim 2 | `brain_met_at_dx == 0` (no CNS involvement at diagnosis) |
| Aim 3 | `any_brain_met == 1` (brain metastasis ever) |
| ML benchmark | Same as Aim 2/3 |

---

## Aim 1 - Association with brain metastasis

**Question:** Which genes and oncogenic pathways are enriched in patients who develop
brain metastases compared to those who do not?

### Methods

- 2x2 contingency tables (gene/pathway mutated vs not, brain-met vs not)
- Fisher exact tests
- Prevalence ratios with log-normal 95% confidence intervals
- Multiple-testing correction (Benjamini-Hochberg FDR)
- Univariate logistic regression of `any_brain_met` on each gene/pathway
- Multivariate logistic regression of `any_brain_met` on top genes + clinical covariates

### Scripts

| Script | Output |
|--------|--------|
| `src/exploratory data analysis/gene_prevtable_oncoprint_forest.py --cohort <cohort>` | Forest plots, oncoprints, prevalence tables |
| `src/exploratory data analysis/aim1_top10_pq_table.py --cohort <cohort>` | Top-10 gene p/q table |
| `src/exploratory data analysis/univariate and multivariate regression df and tables.py --cohort <cohort>` | Logistic regression tables |

### Outputs

- CSV tables: `src/modeling/<cohort>/aim1/`
- Figures / PDFs: `manuscript components/<cohort>/aim1/`
- Committed results: `src/data reports/genie/aim1/` (GENIE only)

---

## Aim 2 - Time to brain metastasis

**Question:** In patients without CNS involvement at diagnosis, what is the time from
diagnosis to brain metastasis, and which features predict it?

**Cohort filter:** `brain_met_at_dx == 0`

### Methods

- Kaplan-Meier curves stratified by top-gene mutation status with log-rank tests
- Cox proportional hazards model with scaled Schoenfeld diagnostics (PH assumption check)
- Accelerated Failure Time (AFT) models: Weibull, lognormal, log-logistic; distribution
  selected by AIC

### Script

```zsh
python "src/modeling/proportional_hazard_afttest.py" --cohort <cohort>
```

### Outputs

- `src/modeling/<cohort>/aim2/` (CSV tables)
- `manuscript components/<cohort>/aim2/` (figures)

---

## Aim 3 - Overall survival

**Question:** In the brain-mets-ever cohort, what predicts overall survival from
diagnosis?

**Cohort filter:** `any_brain_met == 1`

### Methods

Same KM / Cox PH / AFT machinery as Aim 2, with a minimum-events guard:
`MIN_EVENTS = 25`. Fitting is skipped if the number of events is below this threshold.

### Script

```zsh
python "src/modeling/aim3_os.py" --cohort <cohort>
```

### Outputs

- `src/modeling/<cohort>/aim3/` (CSV tables)
- `manuscript components/<cohort>/aim3/` (figures)

---

## ML survival benchmark

**Question:** How well can an XGBoost-AFT model predict survival, and which features
drive predictions?

### Methods

- XGBoost with the Accelerated Failure Time (AFT) objective (`xgb:survival:aft`)
- Evaluated by concordance index (C-index)
- Feature importance and SHAP-based explainability (Python + R)

### Scripts

| Script | Purpose |
|--------|---------|
| `src/modeling/xgb_aft_preprocessing_feature_constuction_train_validate_evaluate.py` | Feature construction, training, evaluation |
| `src/modeling/xgb_aft_feature_.processing_feature_processing_feature_importance.py` | Feature importance |
| `src/modeling/xgb_aft_shap_feature_importance.py` | SHAP analysis (Python) |
| `src/modeling/shap analysis and plot generation.R` | SHAP plots (R, via SHAPforxgboost) |
| `src/modeling/stratifiedKM_CoxFG_feature_prep_AFT.py` | Feature prep + stratified KM + Cox Fine-Gray |

### Outputs

- `reports/figures/ml_benchmark/`

> **Note.** `src/modeling/README_survival_AFT_pipeline.md` and
> `src/modeling/survival_and_xgb_analysis.ipynb` describe a standalone CLI
> (`survival_and_xgb_analysis.py`) with column names (`DFS_MONTHS`, `SUBTYPE`, `G__*`)
> that do not match the current `_lib`-based canonical schema. This appears to document
> an older/separate variant. Reconcile before relying on that README.

---

## Covariates used across aims

All Aim scripts use `_lib.prep_covariates()` and `_lib.drop_low_variance()` for
consistent one-hot encoding and small-cohort safety.

**Stratification candidates for KM / Cox / AFT:**

| Variable | Type |
|----------|------|
| `top5_any_mutated`, `top10_any_mutated` | Binary |
| `receptor_primary_cat` | Categorical (4 levels) |
| `stage_iv_bin`, `stage_diag_group` | Binary / ordered |
| `grade_ord` | Ordered (Low / Intermediate / High) |
| `age_cat` | Ordered (<50 / 50-70 / >70) |
| `pathway_*` (10 columns) | Binary |
| `met_loc` | Categorical (Brain / Other / None) |

**Continuous covariates:**

- `age_dx_num`
- `mutation_count_all_sites_sum`
- `t_alt_count_max`
