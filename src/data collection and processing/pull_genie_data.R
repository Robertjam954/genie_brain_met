suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(data.table)
  library(openxlsx)
})

data_root <- "/Users/robertjames/loc/data private"
proj_data <- "/Users/robertjames/Documents/Documents - Robert’s iMac/Research Projects/MSKCC Research Fellowship/Projects/genie_bpc_brain_mets/data"
out_file  <- file.path(proj_data, "extracted_variables_of_interest.xlsx")
hugo_file <- file.path(proj_data, "hugo_symbols.xlsx")
merged_impact_file <- "/Users/robertjames/Documents/GitHub/IMPACT/merged_impact_patient_level_deid.xlsx"

ds <- list(
  genie  = file.path(data_root, "genie_bpc_datav1"),
  tcga   = file.path(data_root, "lgg_tcga_pan_can_atlas_2018"),
  impact = file.path(data_root, "msk_impact_2017")
)
study_ids <- c(
  genie  = "genie_bpc_brca_v1",
  tcga   = "lgg_tcga_pan_can_atlas_2018",
  impact = "msk_impact_2017"
)

read_tsv_q <- function(p, ...) {
  if (!file.exists(p)) return(NULL)
  suppressWarnings(suppressMessages(
    readr::read_tsv(p, comment = "#", show_col_types = FALSE,
                    progress = FALSE, na = c("","NA"), ...)
  ))
}
read_csv_q <- function(p, ...) {
  if (!file.exists(p)) return(NULL)
  suppressWarnings(suppressMessages(
    readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                    na = c("","NA"), ...)
  ))
}

# ---------------- HUGO annotation ----------------
cat("== reading hugo annotation ==\n")
hugo_full <- openxlsx::read.xlsx(hugo_file, sheet = "hugo")
hugo_lookup <- hugo_full[, c("hgnc_id","symbol","name","locus_type","location",
                             "location_sortable","gene_group","gene_group_id",
                             "entrez_id","ensembl_gene_id")]
hugo_lookup <- as.data.table(hugo_lookup)

attach_hugo <- function(dt, sym_col, prefix = NULL) {
  # left join dt with hugo_lookup on sym_col == symbol
  h <- hugo_lookup
  if (!is.null(prefix)) {
    setnames(h, setdiff(names(h), "symbol"),
             paste0(prefix, setdiff(names(h), "symbol")))
  }
  dt <- merge(dt, h, by.x = sym_col, by.y = "symbol",
              all.x = TRUE, sort = FALSE)
  dt
}

# ---------------- CLINICAL ----------------
cat("== clinical ==\n")

## GENIE
g_pt   <- read_tsv_q(file.path(ds$genie, "data_clinical_patient.txt"))
g_supp <- read_tsv_q(file.path(ds$genie, "data_clinical_supp_survival.txt"))
g_supp_trt <- read_tsv_q(file.path(ds$genie, "data_clinical_supp_survival_treatment.txt"))
g_pl   <- read_csv_q(file.path(ds$genie, "patient_level_dataset.csv"))
g_cl   <- read_csv_q(file.path(ds$genie, "cancer_level_dataset_index.csv"))
g_cl_ni <- read_csv_q(file.path(ds$genie, "cancer_level_dataset_non_index.csv"))

# normalize GENIE CSVs: uppercase column names, rename record_id -> PATIENT_ID
norm_genie_csv <- function(df) {
  if (is.null(df)) return(NULL)
  names(df) <- toupper(names(df))
  if ("RECORD_ID" %in% names(df)) df <- dplyr::rename(df, PATIENT_ID = RECORD_ID)
  df
}
g_pl <- norm_genie_csv(g_pl)
g_cl <- norm_genie_csv(g_cl)
g_cl_ni <- norm_genie_csv(g_cl_ni)

# tag non-index columns with prefix so they don't collide
if (!is.null(g_cl_ni)) {
  non_pid <- setdiff(names(g_cl_ni), "PATIENT_ID")
  names(g_cl_ni)[match(non_pid, names(g_cl_ni))] <- paste0("NONIDX_", non_pid)
}

dedup_by_pid <- function(df) {
  if (is.null(df) || !"PATIENT_ID" %in% names(df)) return(df)
  df[!duplicated(df$PATIENT_ID), , drop = FALSE]
}
g_cl    <- dedup_by_pid(g_cl)
g_cl_ni <- dedup_by_pid(g_cl_ni)
g_pl    <- dedup_by_pid(g_pl)

merge_pid <- function(x, y, suf) {
  if (is.null(y)) return(x)
  z <- dplyr::left_join(x, y, by = "PATIENT_ID", suffix = c("", suf))
  z[, !grepl(paste0(gsub("\\.","\\\\.",suf),"$"), names(z)), drop = FALSE]
}
genie_clin <- g_pt |>
  merge_pid(g_supp,     ".supp") |>
  merge_pid(g_supp_trt, ".supptrt") |>
  merge_pid(g_pl,       ".pl") |>
  merge_pid(g_cl,       ".cl") |>
  merge_pid(g_cl_ni,    ".cl_ni")
genie_clin$study_id <- study_ids[["genie"]]

## TCGA LGG
t_pt  <- read_tsv_q(file.path(ds$tcga, "data_clinical_patient.txt"))
t_hyp <- read_tsv_q(file.path(ds$tcga, "data_clinical_supp_hypoxia.txt"))
tcga_clin <- t_pt |> merge_pid(t_hyp, ".hyp")
tcga_clin$study_id <- study_ids[["tcga"]]

## IMPACT - use merged xlsx as the rich source, plus the cBioPortal patient table
i_pt <- read_tsv_q(file.path(ds$impact, "data_clinical_patient.txt"))
i_pt$study_id <- study_ids[["impact"]]
i_pt$source <- "msk_impact_2017_cbioportal"

i_merged <- openxlsx::read.xlsx(merged_impact_file, sheet = 1)
i_merged$PATIENT_ID <- i_merged$deid_patient_id
i_merged$study_id <- study_ids[["impact"]]
i_merged$source <- "merged_impact_patient_level_deid"

# Combine all clinical
clinical <- dplyr::bind_rows(genie_clin, tcga_clin, i_pt, i_merged)
# put leading cols up front
front <- c("study_id","source","PATIENT_ID")
clinical <- clinical[, c(intersect(front, names(clinical)),
                         setdiff(names(clinical), front))]
cat(sprintf("  clinical: %d rows x %d cols\n", nrow(clinical), ncol(clinical)))

# ---------------- SAMPLE ----------------
cat("== sample ==\n")
g_s <- read_tsv_q(file.path(ds$genie, "data_clinical_sample.txt"))
g_cpt <- read_csv_q(file.path(ds$genie, "cancer_panel_test_level_dataset.csv"))
if (!is.null(g_cpt)) {
  names(g_cpt) <- toupper(names(g_cpt))
  if ("CPT_GENIE_SAMPLE_ID" %in% names(g_cpt))
    g_cpt <- dplyr::rename(g_cpt, SAMPLE_ID = CPT_GENIE_SAMPLE_ID)
  g_cpt <- g_cpt[!duplicated(g_cpt$SAMPLE_ID), ]
}
g_gm <- read_tsv_q(file.path(ds$genie, "data_gene_matrix.txt"))
genie_samp <- g_s
if (!is.null(g_cpt))
  genie_samp <- dplyr::left_join(genie_samp, g_cpt, by = "SAMPLE_ID", suffix = c("",".cpt"))
if (!is.null(g_gm))
  genie_samp <- dplyr::left_join(genie_samp, g_gm,  by = "SAMPLE_ID", suffix = c("",".gm"))
genie_samp <- genie_samp[, !grepl("\\.cpt$|\\.gm$", names(genie_samp))]
genie_samp$study_id <- study_ids[["genie"]]

t_s <- read_tsv_q(file.path(ds$tcga, "data_clinical_sample.txt"))
t_s$study_id <- study_ids[["tcga"]]

i_s <- read_tsv_q(file.path(ds$impact, "data_clinical_sample.txt"))
i_s$study_id <- study_ids[["impact"]]

sample_df <- dplyr::bind_rows(genie_samp, t_s, i_s)
sample_df <- sample_df[, c(intersect(c("study_id","PATIENT_ID","SAMPLE_ID"), names(sample_df)),
                           setdiff(names(sample_df), c("study_id","PATIENT_ID","SAMPLE_ID")))]
cat(sprintf("  sample: %d rows x %d cols\n", nrow(sample_df), ncol(sample_df)))

# ---------------- SEGMENTS ----------------
cat("== segments ==\n")
seg_files <- list(
  genie  = file.path(ds$genie,  "data_cna_hg19.seg.txt"),
  tcga   = file.path(ds$tcga,   "data_cna_hg19.seg"),
  impact = file.path(ds$impact, "data_cna_hg19.seg")
)
read_seg <- function(k) {
  p <- seg_files[[k]]; if (!file.exists(p)) return(NULL)
  s <- fread(p, sep = "\t", colClasses = c(ID = "character", chrom = "character"))
  s[, study_id := study_ids[[k]]]
  setcolorder(s, c("study_id", setdiff(names(s), "study_id")))
  s
}
segments <- rbindlist(lapply(names(seg_files), read_seg), use.names = TRUE, fill = TRUE)
cat(sprintf("  segments: %d rows x %d cols\n", nrow(segments), ncol(segments)))

# ---------------- MUTATIONS ----------------
cat("== mutations ==\n")
mut_paths <- list(
  genie  = file.path(ds$genie,  "data_mutations_extended.txt"),
  tcga   = file.path(ds$tcga,   "data_mutations.txt"),
  impact = file.path(ds$impact, "data_mutations.txt")
)
mut_cols <- c("Tumor_Sample_Barcode","Hugo_Symbol","Entrez_Gene_Id","Transcript_ID","RefSeq",
              "NCBI_Build","Chromosome","Start_Position","End_Position","Strand",
              "Variant_Classification","Variant_Type","Consequence",
              "Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","dbSNP_RS",
              "HGVSc","HGVSp","HGVSp_Short","Protein_position","Codons",
              "Polyphen_Prediction","Polyphen_Score","SIFT_Prediction","SIFT_Score",
              "gnomAD_AF","Hotspot","Mutation_Status","FILTER",
              "t_ref_count","t_alt_count","n_ref_count","n_alt_count",
              "t_depth","n_depth","mutationInCis_Flag")
read_maf <- function(k) {
  p <- mut_paths[[k]]; if (!file.exists(p)) return(NULL)
  # IMPACT has a leading "#sequenced_samples:..." line; skip=#-comments via fill+skip
  hdr_line <- system(sprintf("grep -n -m1 '^Hugo_Symbol' %s", shQuote(p)), intern = TRUE)
  skip_n <- as.integer(sub(":.*","", hdr_line)) - 1L
  m <- fread(p, sep = "\t", skip = skip_n, fill = TRUE,
             quote = "", na.strings = c("","NA"))
  keep <- intersect(mut_cols, names(m))
  m <- m[, ..keep]
  m[, study_id := study_ids[[k]]]
  m
}
mut_list <- lapply(names(mut_paths), read_maf)
mutations <- rbindlist(mut_list, use.names = TRUE, fill = TRUE)
# join hugo
mutations <- attach_hugo(mutations, "Hugo_Symbol")
front <- c("study_id","Tumor_Sample_Barcode","Hugo_Symbol","Entrez_Gene_Id",
           "entrez_id","ensembl_gene_id","location","gene_group")
setcolorder(mutations, c(intersect(front, names(mutations)),
                         setdiff(names(mutations), front)))
cat(sprintf("  mutations: %d rows x %d cols\n", nrow(mutations), ncol(mutations)))

# ---------------- CNA (gene-level, filtered to |val|>=2) ----------------
cat("== cna (gene-level, filtered to |val|>=2) ==\n")
cna_paths <- list(
  genie  = file.path(ds$genie,  "data_CNA.txt"),
  tcga   = file.path(ds$tcga,   "data_cna.txt"),
  impact = file.path(ds$impact, "data_cna.txt")
)
read_cna <- function(k) {
  p <- cna_paths[[k]]; if (!file.exists(p)) return(NULL)
  cna <- fread(p, sep = "\t", na.strings = c("","NA"))
  meta <- intersect(c("Hugo_Symbol","Entrez_Gene_Id","Cytoband"), names(cna))
  id_vars <- meta
  measure_vars <- setdiff(names(cna), id_vars)
  long <- melt(cna, id.vars = id_vars, measure.vars = measure_vars,
               variable.name = "SAMPLE_ID", value.name = "cna_value",
               value.factor = FALSE)
  long <- long[!is.na(cna_value) & abs(as.numeric(cna_value)) >= 2]
  long[, cna_value := as.integer(cna_value)]
  long[, SAMPLE_ID := as.character(SAMPLE_ID)]
  long[, study_id := study_ids[[k]]]
  long
}
cna_list <- lapply(names(cna_paths), read_cna)
cna <- rbindlist(cna_list, use.names = TRUE, fill = TRUE)
cna <- attach_hugo(cna, "Hugo_Symbol")
front <- c("study_id","SAMPLE_ID","Hugo_Symbol","Entrez_Gene_Id","cna_value",
           "entrez_id","ensembl_gene_id","location","gene_group")
setcolorder(cna, c(intersect(front, names(cna)),
                   setdiff(names(cna), front)))
cat(sprintf("  cna: %d rows x %d cols\n", nrow(cna), ncol(cna)))

# ---------------- SV ----------------
cat("== sv ==\n")
sv_paths <- list(
  genie  = file.path(ds$genie,  "data_sv.txt"),
  tcga   = file.path(ds$tcga,   "data_sv.txt"),
  impact = file.path(ds$impact, "data_sv.txt")
)
sv_cols <- c("Sample_Id","SV_Status",
             "Site1_Hugo_Symbol","Site1_Entrez_Gene_Id","Site1_Ensembl_Transcript_Id",
             "Site1_Chromosome","Site1_Position","Site1_Region","Site1_Description",
             "Site2_Hugo_Symbol","Site2_Entrez_Gene_Id","Site2_Ensembl_Transcript_Id",
             "Site2_Chromosome","Site2_Position","Site2_Region","Site2_Description",
             "Site2_Effect_on_Frame","Site2_Effect_On_Frame",
             "Class","Event_Info","Connection_Type","Breakpoint_Type",
             "SV_Length","NCBI_Build","Annotation",
             "Tumor_Read_Count","Tumor_Split_Read_Count","Tumor_Paired_End_Read_Count",
             "Normal_Read_Count","Normal_Variant_Count")
read_sv <- function(k) {
  p <- sv_paths[[k]]; if (!file.exists(p)) return(NULL)
  s <- fread(p, sep = "\t", na.strings = c("","NA"))
  keep <- intersect(sv_cols, names(s))
  s <- s[, ..keep]
  s[, study_id := study_ids[[k]]]
  s
}
sv <- rbindlist(lapply(names(sv_paths), read_sv), use.names = TRUE, fill = TRUE)
# join hugo for Site1 (with prefix)
sv <- merge(sv,
            setNames(hugo_lookup, c("hgnc_id","Site1_Hugo_Symbol","name","locus_type",
                                    "location","location_sortable","gene_group",
                                    "gene_group_id","entrez_id","ensembl_gene_id")),
            by = "Site1_Hugo_Symbol", all.x = TRUE, sort = FALSE)
setnames(sv,
         c("hgnc_id","name","locus_type","location","location_sortable",
           "gene_group","gene_group_id","entrez_id","ensembl_gene_id"),
         paste0("site1_", c("hgnc_id","name","locus_type","location","location_sortable",
                            "gene_group","gene_group_id","entrez_id","ensembl_gene_id")))
sv <- merge(sv,
            setNames(hugo_lookup, c("hgnc_id","Site2_Hugo_Symbol","name","locus_type",
                                    "location","location_sortable","gene_group",
                                    "gene_group_id","entrez_id","ensembl_gene_id")),
            by = "Site2_Hugo_Symbol", all.x = TRUE, sort = FALSE)
setnames(sv,
         c("hgnc_id","name","locus_type","location","location_sortable",
           "gene_group","gene_group_id","entrez_id","ensembl_gene_id"),
         paste0("site2_", c("hgnc_id","name","locus_type","location","location_sortable",
                            "gene_group","gene_group_id","entrez_id","ensembl_gene_id")))
front <- c("study_id","Sample_Id","Site1_Hugo_Symbol","Site2_Hugo_Symbol","Class",
           "Event_Info","SV_Status",
           "site1_entrez_id","site1_ensembl_gene_id","site1_location","site1_gene_group",
           "site2_entrez_id","site2_ensembl_gene_id","site2_location","site2_gene_group")
setcolorder(sv, c(intersect(front, names(sv)),
                  setdiff(names(sv), front)))
cat(sprintf("  sv: %d rows x %d cols\n", nrow(sv), ncol(sv)))

# ---------------- WRITE WORKBOOK ----------------
cat("== writing workbook ==\n")
dir.create(proj_data, recursive = TRUE, showWarnings = FALSE)
wb <- openxlsx::createWorkbook()
add <- function(name, df) {
  openxlsx::addWorksheet(wb, name)
  # openxlsx doesn't like data.table; ensure data.frame
  openxlsx::writeData(wb, name, as.data.frame(df), keepNA = FALSE)
  openxlsx::freezePane(wb, name, firstRow = TRUE)
}
add("clinical", clinical)
add("sample", sample_df)
add("segments", segments)
add("mutations", mutations)
add("cna",       cna)
add("sv",        sv)
add("hugo_annotation", hugo_lookup)
openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
cat("\nWrote:", out_file, "\n")
cat("Sheet row counts:\n")
print(c(clinical = nrow(clinical), sample = nrow(sample_df),
        segments = nrow(segments), mutations = nrow(mutations),
        cna = nrow(cna), sv = nrow(sv),
        hugo_annotation = nrow(hugo_lookup)))
cat("\nClinical rows by study:\n");  print(table(clinical$study_id, useNA = "ifany"))
cat("Sample rows by study:\n");     print(table(sample_df$study_id, useNA = "ifany"))
cat("Mutations rows by study:\n");  print(table(mutations$study_id, useNA = "ifany"))
cat("CNA rows by study:\n");        print(table(cna$study_id, useNA = "ifany"))
cat("SV rows by study:\n");         print(table(sv$study_id, useNA = "ifany"))
