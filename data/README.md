# data/

**Code-only repository.** Under the GENIE BPC data-use agreement, no
patient-level or controlled-access data is committed. Only this README and the
folder structure (`.gitkeep`) are tracked; every real data file is git-ignored.

| Folder | Purpose |
|--------|---------|
| `raw/` | Original release drops (GENIE BPC v1.0, TCGA-BRCA, MSK-IMPACT 2018), immutable. |
| `external/` | Third-party reference files (hotspots, hugo symbols, data guides). |
| `interim/` | Intermediate harmonized/merged frames, regenerable. |
| `processed/` | Final analytic frames the analysis reads from. |

Flow: `raw` (+ `external`) -> `interim` -> `processed`. The analytic frame
(`extracted_variables_genie_data.csv`) is built by the scripts in
`src/data collection and processing/` - see the README pipeline order.
