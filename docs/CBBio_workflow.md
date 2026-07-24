# CBBio candidate driver-gene workflow

Run notes for the candidate driver-gene fusion workflow diagrammed in
`notebooks/Workflow_CBBio.pdf` (= `src/Workflow_CBBio.pdf`). This reproduces the
author-provided **CBBio v0.1** package notes (`ReadMe.txt`, `Command_Line_Software_Tools`)
so the workflow is documented in-repo. See `references/REFERENCES.md` §2-§3 for the method
citations and tool sources.

> This is a reference workflow. The four caller tools and the ensemble scripts below are
> external to this repo (large binaries / MATLAB / legacy Python 2) and are **not** committed
> here. The `data/` inputs are git-ignored. Treat this as the specification to follow, not a
> runnable pipeline in this checkout.

## Overview

CBBio is a fusion system for Whole-Exome Sequencing (WES) / MAF (`.maf`) data that identifies
candidate driver genes from gene mutations. The initial input in the source study is a
metastatic breast-cancer cohort (Wagle et al. 2017), available via cBioPortal.

1. **Feature extraction** - four caller tools rank genes by p-value; each caller assigns one
   number (feature) per gene:
   - MutSigCV v1.4, OncodriveCLUST, OncodriveFM, NetBox 1.0.
2. **Ensemble learning** - three classifiers (non-linear SVM, ANN, Random Forest) label genes
   0 (passenger) / 1 (driver) and predict a per-gene score.
3. **Algebraic combiner** - the ensemble averages the three per-gene scores for the final
   decision, producing the candidate-driver-gene set.

## Feature-extraction tools

### MutSigCV v1.4 (MATLAB)

Requires a MATLAB full-toolbox license (or MCR) and the reference input files. Source and
reference files: https://software.broadinstitute.org/cancer/cga/mutsig

```
% in MATLAB:
addpath <MutSigCV dir>
MutSigCV( ...
  '<cohort>/data_mutations_extended.maf', ...
  '<MutSigCV>/Required inputs/exome_full192.coverage.txt', ...
  '<MutSigCV>/Required inputs/gene.covariates.txt', ...
  '<out>/output.txt', ...
  '<MutSigCV>/Required inputs/mutation_type_dictionary_file.txt', ...
  '<MutSigCV>/Required inputs/chr_files_hg19')
```

### OncodriveCLUST (Python 3 / Anaconda 3)

Source: https://bitbucket.org/bbglab/oncodriveclust

```
oncodriveclust -m 3 --cgc <data>/CGC_phenotype.tsv \
  <examples>/nonsyn-mbca.txt <examples>/syn-mbca.txt \
  <data>/gene_transcripts.tsv
```

### OncodriveFM (Python 3 / Anaconda 3)

Source: https://bitbucket.org/bbglab/oncodrivefm

```
pip install OncodriveFM
oncodrivefm -e median -m <data>/ensg_kegg.tsv <data>/MBCA1.tdm
```

### NetBox 1.0 (Java 1.5+, Python 2.5+)

Source: http://cbio.mskcc.org/netbox

```
cd <netbox>/brca_data
python netAnalyze.py <netbox>/brca_data/netbox1.props
```

## Ensemble machine-learning methods

Requires Python 2.5+ (Anaconda 2) and scikit-learn. Replace the default feature file (`.csv`),
positive/negative label files (`.txt`), and test file (`.txt`) with your own, preserving file
names and formats.

```
# SVM
python svm_gene_classification.py -cvr 100
# ANN
python nn_gene_classification.py -cvr 100
# Random Forest
python rf_gene_classification.py -cvr 100
# Ensemble (algebraic combiner)
python all_gene_classification.py -cvr 100
```

Each run writes two files: `*_predictions_train.csv` (labels + scores for train genes) and
`*_predictions_test.csv` (labels + scores for test genes). Genes are prioritized by predicted
importance. Rank and interrogate the final candidate driver genes using TCGA / cBioPortal
analytical tools.

## Relationship to this repo

- `notebooks/Workflow_CBBio.pdf`, `src/Workflow_CBBio.pdf` - the workflow diagram.
- `references/REFERENCES.md` - method citations (§2) and tool sources (§3).
- The four callers, their reference files, a cohort MAF, and the ensemble scripts are
  external prerequisites supplied on the author's machine, not tracked here.
