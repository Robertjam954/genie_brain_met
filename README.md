# GENIE / TCGA / IMPACT — Breast-Cancer Brain Metastasis

A reproducible clinical-genomics research pipeline investigating the genomic and
clinical drivers of **breast-cancer brain metastasis**. Three public clinical-genomic
cohorts — GENIE BPC BRCA, TCGA, and Breast MSK 2018 — are harmonized to a single
canonical schema and run through a consistent set of statistical and machine-learning
analyses feeding a scientific manuscript.

## Analytic aims

| Aim | Question | Methods |
|-----|----------|---------|
| **1 — Association** | Which genes / oncogenic pathways are enriched in brain-met patients? | Fisher exact, prevalence ratios (log-normal 95% CI), BH correction, logistic regression, forest plots, oncoprints |
| **2 — Time to brain met** | Time from diagnosis to brain metastasis (no-CNS-at-dx cohort) | Kaplan-Meier + log-rank, Cox PH + Schoenfeld, AFT distribution selection |
| **3 — Overall survival & PFS** | Survival in the brain-mets-ever cohort | Same KM / Cox / AFT machinery, minimum-events guard |
| **ML benchmark** | Predictive survival modeling | XGBoost-AFT (primary), Random Survival Forest, Gradient-Boosted Survival, LightGBM; Harrell + IPCW C-index; SHAP explainability |

## Repository layout

```
notebooks/        Executed analysis notebooks (Aims 1–3 + ML benchmark)
src/
  data collection and processing/   Harmonization ETL + canonical schema spec
  exploratory data analysis/        Aim 1 association analyses
  modeling/                         Aim 2/3 survival + XGBoost-AFT modeling
docs/
  results_dashboard.html            Interactive browser of all 32 figures
  plans/                            Execution / implementation plans
manuscript components/              Publication figures (PDF/PNG)
references/                         Source papers, manuscript structure
```

See [`PRODUCT.md`](PRODUCT.md) and [`ARCHITECTURE.md`](ARCHITECTURE.md) for the full
product and architecture overview, and [`docs/results_dashboard.html`](docs/results_dashboard.html)
for the interactive figure dashboard.

## Data

Raw and derived data (cBioPortal / GENIE exports, harmonized analytic CSVs) are **not**
tracked in git — they are large and de-identified but out of scope for VCS. The pipeline
reads them from a local `data/processed/` root. This project operates on de-identified
public datasets for research and manuscript preparation; it is **not** a clinical tool.
