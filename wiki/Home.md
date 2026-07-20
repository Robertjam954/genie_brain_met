# GENIE Brain Metastasis - Project Wiki

Welcome to the wiki for the **GENIE-only breast-cancer brain metastasis** analysis
pipeline.

> **Status notice.** This repository is an active research work-in-progress. Reported
> values are provisional and pending independent verification.

---

## Scope

This repository is now scoped to **GENIE BPC BRCA only**. TCGA and MSK replication
cohorts are out of scope for this repo/wiki.

---

## Wiki pages

| Page | Description |
|------|-------------|
| [Getting Started](Getting-Started.md) | Environment setup, prerequisites, and data placement |
| [Architecture and Data Flow](Architecture-and-Data-Flow.md) | High-level pipeline architecture and file layout |
| [Pipeline Run Guide](Pipeline-Run-Guide.md) | Step-by-step instructions to run the full GENIE pipeline |
| [Data Schema and Harmonization](Data-Schema-and-Harmonization.md) | Canonical schema and GENIE harmonization notes |
| [Analysis Aims](Analysis-Aims.md) | Aim 1, Aim 2, Aim 3, and ML benchmark for GENIE |
| [Known Issues and Gaps](Known-Issues-and-Gaps.md) | Current gaps and remaining work |
| [Contributing](Contributing.md) | Coding guidelines and contributor workflow |

---

## Quick-start summary

```zsh
# 1 - Pull raw data (R)
Rscript "src/data collection and processing/pull_genie_data.R"

# 2 - Harmonize GENIE (Python)
python "src/data collection and processing/harmonize_genie.py"

# 3 - Enrich GENIE analytic frame (Python)
python "src/data collection and processing/enrich_harmonized.py"

# 4 - Ensure analytic files are available in data/processed/

# 5 - Run Aim 1 association analysis
python "src/exploratory data analysis/gene_prevtable_oncoprint_forest.py" --cohort genie

# 6 - Run survival analyses
python "src/modeling/proportional_hazard_afttest.py" --cohort genie
python "src/modeling/aim3_os.py" --cohort genie
```

See the [Pipeline Run Guide](Pipeline-Run-Guide.md) for full details.
