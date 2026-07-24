# References

Bibliography for the GENIE BPC breast-cancer brain-metastasis project. It collects
(1) the source-paper PDFs kept in this `references/` folder, (2) the candidate
driver-gene ("CBBio") method and the papers it cites, and (3) the external
driver-caller software the CBBio workflow depends on. Citations are provided so the
prep docs (`README.md`, `PRODUCT.md`, `ARCHITECTURE.md`, `CONTRIBUTING.md`) can point
here instead of restating full references inline.

## 1. Source-paper PDFs in this repo

Kept under `references/` (not part of the analytic pipeline; background reading):

- `Genomics of Breast Cancer Brain Metastases- A Meta-Analysis and Therapeutic Implications.pdf`
  - Meta-analysis of the genomic landscape of breast-cancer brain metastases; motivates
    the Aim 1 association analysis.
- `Survival regression with accelerated failure time model in XGBoost.pdf`
  - Barnwal, Cho, Hocking. Survival regression with the Accelerated Failure Time (AFT)
    objective in XGBoost. *Journal of Computational and Graphical Statistics* (2022).
    Basis for the `xgb_aft_*` modeling scripts (`survival:aft` objective).
- `transformer survival analysis and xgb.pdf`
  - Background on transformer-based and gradient-boosted survival modeling.
- `reactome_pathway.pdf`
  - Reactome pathway reference used alongside the Sanchez-Vega oncogenic pathway set.
- `Optimization_An Evaluation of Large Language Models on Text Summarization Tasks Using Prompt .pdf`
  - LLM text-summarization evaluation (adjacent-work reference).
- `🧬 Breast Cancer Brain Metastasis Study Materials.pdf`
  - Collected study materials / notes for the brain-metastasis project.

## 2. Candidate driver-gene method (CBBio)

The candidate driver-gene fusion workflow is specified by the diagram
`notebooks/Workflow_CBBio.pdf` (duplicated at `src/Workflow_CBBio.pdf`). It is the
**CBBio v0.1** method by **Leila Mirsadeghi**: four driver-caller tools rank genes by
p-value, each caller contributes one `-log10(p)` feature per gene, and an ensemble of
non-linear SVM + ANN + Random Forest classifiers combines the features via an algebraic
(average-score) combiner to output a final candidate-driver-gene set. Labels are 0
(passenger) / 1 (driver).

Method package documents (provided by the author): the CBBio `ReadMe.txt`, the
`Command_Line_Software_Tools` run notes, and the thesis `Additional_Files_*.xlsx`
supplementary tables. These are external artifacts and are not committed to this repo
(large binaries / data); see `docs/CBBio_workflow.md` for the run notes reproduced from
them.

### Papers cited by the CBBio method

1. N. Wagle et al., "The Metastatic Breast Cancer (MBC) project: Accelerating
   translational research through direct patient engagement." *American Society of
   Clinical Oncology*, 2017. (Initial input cohort.)
2. E. Cerami et al., "The cBio cancer genomics portal: an open platform for exploring
   multidimensional cancer genomics data." *Cancer Discovery / AACR*, 2012.
3. J. Gao et al., "Integrative analysis of complex cancer genomics and clinical profiles
   using the cBioPortal," *Science Signaling*, vol. 6, no. 269, pp. pl1, 2013.
4. M. S. Lawrence et al., "Mutational heterogeneity in cancer and the search for new
   cancer-associated genes," *Nature*, vol. 499, no. 7457, pp. 214-218, 2013. (MutSigCV.)
5. D. Tamborero, A. Gonzalez-Perez, N. Lopez-Bigas, "OncodriveCLUST: exploiting the
   positional clustering of somatic mutations to identify cancer genes,"
   *Bioinformatics*, vol. 29, no. 18, pp. 2238-2244, 2013.
6. A. Gonzalez-Perez, N. Lopez-Bigas, "Functional impact bias reveals cancer drivers,"
   *Nucleic Acids Research*, vol. 40, no. 21, pp. e169, 2012. (OncodriveFM.)
7. E. Cerami, E. Demir, N. Schultz, B. S. Taylor, C. Sander, "Automated network analysis
   identifies core pathways in glioblastoma," *PLoS One*, vol. 5, no. 2, p. e8918, 2010.
   (NetBox.)
8. L. Rokach, "Ensemble-based classifiers," *Artificial Intelligence Review*, vol. 33,
   no. 1-2, pp. 1-39, 2010. (Algebraic combiners.)
9. C. Cortes, V. Vapnik, "Support-vector networks," *Machine Learning*, vol. 20, no. 3,
   pp. 273-297, 1995.
10. F. Rosenblatt, "The perceptron: a probabilistic model for information storage and
    organization in the brain," *Psychological Review*, vol. 65, no. 6, p. 386, 1958.
11. J. Schmidhuber, "Deep learning in neural networks: An overview," *Neural Networks*,
    vol. 61, pp. 85-117, 2015.
12. T. K. Ho, "Random decision forests," in *Proc. 3rd Int. Conf. Document Analysis and
    Recognition*, 1995, vol. 1, pp. 278-282.
13. R. Polikar, "Ensemble based systems in decision making," *IEEE Circuits and Systems
    Magazine*, vol. 6, no. 3, pp. 21-45, 2006.

## 3. Driver-caller software (CBBio dependencies)

Feature extraction ranks genes by p-value with four external tools. Sources and example
commands are reproduced in `docs/CBBio_workflow.md`.

| Tool | Version | Requirements | Source |
|------|---------|--------------|--------|
| MutSigCV | 1.4 | MATLAB (full toolbox) or MCR; reference coverage/covariate files | https://software.broadinstitute.org/cancer/cga/mutsig |
| OncodriveCLUST | - | Python 3 (numpy, scipy, pandas, statsmodels); Anaconda 3 | https://bitbucket.org/bbglab/oncodriveclust |
| OncodriveFM | - | Python 3 (numpy, scipy, pandas, statsmodels); Anaconda 3 | https://bitbucket.org/bbglab/oncodrivefm |
| NetBox | 1.0 | Java 1.5+; Python 2.5+ | http://cbio.mskcc.org/netbox |

Ensemble ML methods use scikit-learn (https://scikit-learn.org/stable/).

## 4. Data sources

- AACR Project GENIE BPC (Biopharma Collaborative) breast-cancer (BRCA) release,
  via cBioPortal. HUGO gene annotations from `data/hugo_symbols.xlsx`.
- Sanchez-Vega et al. ten oncogenic signaling pathways (used for the `pathway_*` feature
  set; see `harmonization_spec.md` sections 13-14).

## 5. Workflow methodology (background)

- `archive/context-engineering-workflow.md` - the context-engineering prep-docs +
  plan-then-implement workflow these docs follow.
- "The Agentic Operating Model" (author-provided PDF, not committed) - background on the
  agentic operating model that informs the AI-engineering workflow scaffolding kept in
  `archive/`.
