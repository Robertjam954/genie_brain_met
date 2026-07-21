# GENIE BPC - Breast-Cancer Brain Metastasis

Quantifies how genes and clinical risk factors shape the likelihood and timing of breast-cancer brain metastasis using GENIE BPC data. Research analysis feeding a manuscript.

## Quick Facts
- **Type**: Research / data-science analysis (feeds a manuscript, not a deployed service)
- **Language**: Python (see per-notebook / per-script imports; `src/` holds the pipeline)
- **Deploy**: Azure (see project deployment notes); Azure subscription is currently on billing hold - deploys are blocked until reactivated.
- **No formal test suite** - do not invent commands that don't exist.

## Key Directories
- `src/data collection and processing/` - ingestion + cleaning
- `src/exploratory data analysis/` - EDA
- `src/modeling/` - survival / risk modeling
- `src/services/` - shared helpers
- `data/` - datasets (treat as sensitive; do not commit raw patient data)
- `notebooks/` - exploratory notebooks
- `reports/` - figures and generated artifacts
- `manuscript/` - manuscript drafts and figures
- `PRODUCT.md`, `ARCHITECTURE.md`, `STATUS.md` - product/architecture/status context

## Working Rules
- **PHI / patient data**: never print, commit, or send patient-level data externally. Everything under `data/` is sensitive.
- **Missingness before EDA**: run a missingness audit on any finalized dataset before EDA / Table 1 / modeling; the analyst decides imputation. See the `data-analysis-hygiene` skill.
- **Reuse cached outputs**: integrate from existing extraction/feature/metric outputs rather than re-running.
- **Figures**: follow the `publication-figures` skill (vector PDF + 600 DPI PNG, Arial, Okabe-Ito, grayscale-safe, error bars + n + significance).
- **Prose**: single hyphen (-), never em dashes.
- **Workflow**: follow `context-engineering-workflow` (curate context -> plan -> implement) for new work.
- **Models**: LLM/pipeline work uses `claude-sonnet-5`.
