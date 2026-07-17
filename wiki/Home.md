# GENIE / TCGA / IMPACT - Breast-Cancer Brain Metastasis: Project Wiki

Welcome to the project wiki for the **GENIE / TCGA / IMPACT Breast-Cancer Brain Metastasis** analysis pipeline.

> **Status notice.** This repository is an active research work-in-progress. All reported
> numbers (cohort sizes, p/q-values, hazard ratios, concordance indices) are provisional
> and have not been independently re-verified. This wiki documents pipeline structure and
> conventions, not confirmed scientific findings.

---

## What is this project?

A clinical-genomics research pipeline investigating the genomic and clinical drivers of
**breast-cancer brain metastasis**. Three public clinical-genomic cohorts are harmonized
to a single canonical schema and analyzed with a consistent statistical and machine-learning
pipeline whose outputs feed a scientific manuscript.

| Cohort | Key | Source |
|--------|-----|--------|
| GENIE BPC BRCA | `genie` | AACR Project GENIE BPC |
| TCGA BRCA | `tcga` | The Cancer Genome Atlas |
| Breast MSK 2018 | `msk18` | cBioPortal / MSK-IMPACT |

---

## Wiki pages

| Page | Description |
|------|-------------|
| [Getting Started](Getting-Started.md) | Environment setup, prerequisites, and data placement |
| [Architecture and Data Flow](Architecture-and-Data-Flow.md) | High-level pipeline architecture and file layout |
| [Pipeline Run Guide](Pipeline-Run-Guide.md) | Step-by-step instructions to run the full pipeline |
| [Data Schema and Harmonization](Data-Schema-and-Harmonization.md) | Canonical schema, column names, and recoding rules |
| [Analysis Aims](Analysis-Aims.md) | Details on Aim 1 (association), Aim 2 (time-to-brain-met), Aim 3 (OS), and ML benchmark |
| [Known Issues and Gaps](Known-Issues-and-Gaps.md) | Current gaps, broken paths, and work remaining |
| [Contributing](Contributing.md) | Coding guidelines, commit rules, and developer conventions |

---

## Quick-start summary

```zsh
# 1 - Pull raw data (R)
Rscript "src/data collection and processing/pull_genie_data.R"

# 2 - Harmonize each cohort (Python)
python "src/data collection and processing/harmonize_genie.py"
python "src/data collection and processing/harmonize_tcga.py"
python "src/data collection and processing/harmonize_breast_msk_2018.py"

# 3 - Enrich harmonized frames (Python)
python "src/data collection and processing/enrich_harmonized.py"

# 4 - Copy analytic CSVs to data/processed/ (manual step)

# 5 - Run Aim 1 association analysis
python "src/exploratory data analysis/gene_prevtable_oncoprint_forest.py" --cohort genie

# 6 - Run survival analyses
python "src/modeling/proportional_hazard_afttest.py" --cohort genie   # Aim 2
python "src/modeling/aim3_os.py" --cohort genie                        # Aim 3
```

See the [Pipeline Run Guide](Pipeline-Run-Guide.md) for full details.

---

## Key documents in the repo

| File | Purpose |
|------|---------|
| `README.md` | Repo overview, layout, setup, run order |
| `ARCHITECTURE.md` | Architecture and component descriptions |
| `PRODUCT.md` | Product/analytic aims overview |
| `CONTRIBUTING.md` | Developer conventions |
| `STATUS.md` | Live task checklist |
| `src/modeling/_lib.py` | Shared pipeline backbone (cohort specs, loaders, covariate prep) |
| `src/data collection and processing/harmonization_spec.md` | Canonical schema contract |
| `docs/results_dashboard.html` | Rendered figure dashboard (provisional) |
