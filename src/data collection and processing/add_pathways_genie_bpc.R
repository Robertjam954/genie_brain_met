# Build the GENIE BPC BRCA genomic feature block: a per-sample gene-binary matrix
# plus the 10 Sanchez-Vega oncogenic pathway indicators, merged onto the clinical
# sample master.
#
# This is the upstream step that produces `genie_bpc_v1_sample_master_full.csv`
# (clinical + gene_binary + pathways), which `harmonize_genie.py` consumes. It
# implements sections 13 (gene-binary matrix) and 14 (pathway annotation) of
# `harmonization_spec.md` using gnomeR.
#
# Inputs (cBioPortal / GENIE BPC exports; adjust the machine-specific paths below):
#   - MAF                : data_mutations_extended.txt  (long-format mutations)
#   - CNA (wide)         : data_CNA.txt                 (genes x samples, -2..2)
#   - SV / fusions       : data_sv.txt                  (Site1/Site2 hugo symbols)
#   - clinical master    : genie_bpc_v1_sample_master_clinical.csv (one row per SAMPLE_ID)
#
# Output:
#   - <proc_dir>/genie_bpc_v1_sample_master_full.csv
#
# gnomeR settings match harmonization_spec.md section 13:
#   mut_type="somatic_only", include_silent=FALSE, snp_only=FALSE,
#   high_level_cna_only=TRUE (+2 Amp / -2 Del only), specify_panel="no",
#   recode_aliases="no". Samples absent from MAF/CNA/SV are kept as all-zero rows.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(rlang)
  library(tibble)
  library(data.table)
  library(gnomeR)
})

# ---------------- Config (edit these for your machine) ----------------
data_root <- "/Users/robertjames/loc/data private/genie_bpc_datav1"
proc_dir  <- "data/processed"

maf_file      <- file.path(data_root, "data_mutations_extended.txt")
cna_file      <- file.path(data_root, "data_CNA.txt")
sv_file       <- file.path(data_root, "data_sv.txt")
clinical_file <- file.path(proc_dir,  "genie_bpc_v1_sample_master_clinical.csv")
out_file      <- file.path(proc_dir,  "genie_bpc_v1_sample_master_full.csv")

# Canonical pathway columns expected downstream (see src/modeling/_lib.py:PATHWAY_COLS).
PATHWAY_COLS <- c("pathway_RTK/RAS", "pathway_Nrf2", "pathway_PI3K", "pathway_TGFB",
                  "pathway_p53", "pathway_Wnt", "pathway_Myc", "pathway_Cell cycle",
                  "pathway_Hippo", "pathway_Notch")

# ---------------- Read helpers ----------------
read_tsv_q <- function(p, ...) {
  if (!file.exists(p)) {
    cat(sprintf("  WARN: missing %s (skipping)\n", p))
    return(NULL)
  }
  suppressWarnings(suppressMessages(
    readr::read_tsv(p, comment = "#", show_col_types = FALSE,
                    progress = FALSE, na = c("", "NA"), ...)
  ))
}

# Lowercase names so gnomeR gets the snake_case column names it expects.
lc_names <- function(df) {
  names(df) <- tolower(names(df))
  df
}

# ---------------- Clinical master (defines the sample universe) ----------------
cat("== clinical master ==\n")
clinical <- readr::read_csv(clinical_file, show_col_types = FALSE, progress = FALSE)
if (!"SAMPLE_ID" %in% names(clinical)) {
  # tolerate common synonyms (see harmonization_spec.md section 18)
  syn <- intersect(c("Tumor_Sample_Barcode", "Sample_Id", "sampleId",
                     "cpt_genie_sample_id"), names(clinical))
  if (length(syn) == 0) stop("clinical master has no SAMPLE_ID column or known synonym")
  clinical <- dplyr::rename(clinical, SAMPLE_ID = !!rlang::sym(syn[1]))
}
all_samples <- unique(as.character(clinical$SAMPLE_ID))
cat(sprintf("  clinical: %d rows, %d unique SAMPLE_ID\n",
            nrow(clinical), length(all_samples)))

# ---------------- Mutations (MAF, long) ----------------
cat("== mutations ==\n")
maf_raw <- read_tsv_q(maf_file)
mutation <- NULL
if (!is.null(maf_raw)) {
  maf_raw <- lc_names(maf_raw)
  # gnomeR expects sample_id + hugo_symbol; map from cBioPortal names.
  if (!"sample_id" %in% names(maf_raw) && "tumor_sample_barcode" %in% names(maf_raw))
    maf_raw <- dplyr::rename(maf_raw, sample_id = tumor_sample_barcode)
  mutation <- maf_raw
  cat(sprintf("  maf: %d rows\n", nrow(mutation)))
}

# ---------------- CNA (wide -> long) ----------------
cat("== cna (wide -> long) ==\n")
cna_long <- NULL
if (file.exists(cna_file)) {
  cna_wide <- data.table::fread(cna_file, sep = "\t", na.strings = c("", "NA"))
  meta <- intersect(c("Hugo_Symbol", "Entrez_Gene_Id", "Cytoband"), names(cna_wide))
  cna_long <- data.table::melt(
    cna_wide, id.vars = meta,
    variable.name = "sample_id", value.name = "alteration",
    variable.factor = FALSE, value.factor = FALSE
  )
  data.table::setnames(cna_long, "Hugo_Symbol", "hugo_symbol", skip_absent = TRUE)
  cna_long <- cna_long[!is.na(alteration)]
  cna_long[, alteration := as.numeric(alteration)]
  cna_long[, sample_id := as.character(sample_id)]
  # high_level_cna_only in create_gene_binary keeps only +/-2; pre-trim to shrink.
  cna_long <- cna_long[abs(alteration) >= 2, c("sample_id", "hugo_symbol", "alteration")]
  cat(sprintf("  cna long (|val|>=2): %d rows\n", nrow(cna_long)))
} else {
  cat(sprintf("  WARN: missing %s (skipping CNA)\n", cna_file))
}

# ---------------- SV / fusions ----------------
cat("== sv / fusions ==\n")
fusion <- NULL
sv_raw <- read_tsv_q(sv_file)
if (!is.null(sv_raw)) {
  sv_raw <- lc_names(sv_raw)
  # gnomeR fusion input: sample_id, site_1_hugo_symbol, site_2_hugo_symbol
  ren <- c(site1_hugo_symbol = "site_1_hugo_symbol",
           site2_hugo_symbol = "site_2_hugo_symbol")
  for (src in names(ren)) {
    if (src %in% names(sv_raw))
      sv_raw <- dplyr::rename(sv_raw, !!ren[[src]] := !!rlang::sym(src))
  }
  fusion <- sv_raw
  cat(sprintf("  sv: %d rows\n", nrow(fusion)))
}

# ---------------- Gene-binary matrix (spec section 13) ----------------
cat("== gene binary (gnomeR::create_gene_binary) ==\n")
gene_binary <- gnomeR::create_gene_binary(
  samples             = all_samples,
  mutation            = mutation,
  cna                 = cna_long,
  fusion              = fusion,
  mut_type            = "somatic_only",
  include_silent      = FALSE,
  snp_only            = FALSE,
  high_level_cna_only = TRUE,
  specify_panel       = "no",
  recode_aliases      = "no"
)
cat(sprintf("  gene_binary: %d samples x %d gene cols\n",
            nrow(gene_binary), ncol(gene_binary)))

# ---------------- Pathways (spec section 14) ----------------
cat("== add_pathways (Sanchez-Vega) ==\n")
gene_binary <- gnomeR::add_pathways(gene_binary)  # default 10 Sanchez-Vega pathways

added <- intersect(PATHWAY_COLS, names(gene_binary))
missing_pw <- setdiff(PATHWAY_COLS, names(gene_binary))
cat(sprintf("  pathways added: %d/%d\n", length(added), length(PATHWAY_COLS)))
if (length(missing_pw) > 0)
  cat(sprintf("  WARN: expected pathway cols not produced: %s\n",
              paste(missing_pw, collapse = ", ")))

# gnomeR returns sample ids as rownames; lift them into a SAMPLE_ID column.
gene_binary <- tibble::rownames_to_column(as.data.frame(gene_binary), var = "SAMPLE_ID")

# ---------------- Merge onto clinical + write ----------------
cat("== merge + write ==\n")
full <- dplyr::left_join(clinical, gene_binary, by = "SAMPLE_ID")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(full, out_file)

cat(sprintf("\nWrote: %s\n", out_file))
cat(sprintf("  %d rows x %d cols (clinical + gene_binary + %d pathways)\n",
            nrow(full), ncol(full), length(added)))
