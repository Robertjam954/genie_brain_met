#!/usr/bin/env Rscript
## Create a de-identified version of merged_impact_patient_level.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

outdir <- file.path("C:", "Users", "jamesr4", "OneDrive - Memorial Sloan Kettering Cancer Center", "Documents", "Research", "Projects", "impact", "Output")
infile <- file.path(outdir, 'merged_impact_patient_level.csv')
outfile <- file.path(outdir, 'merged_impact_patient_level_deid.csv')
if (!file.exists(infile)) stop('Input merged file not found: ', infile)

message('Reading: ', infile)
df <- read_csv(infile, show_col_types = FALSE)

# Choose a patient identifier to pseudonymize (prefer merge_key -> SAMPLE_ID -> Study ID)
id_candidates <- c('merge_key','SAMPLE_ID','Study ID','Sample ID','Sample_ID')
id_col <- intersect(id_candidates, names(df))
id_col <- if (length(id_col)>0) id_col[1] else NULL

if (!is.null(id_col)){
  message('Using id column: ', id_col)
  unique_ids <- unique(df[[id_col]])
  # create pseudonyms P000001, P000002 ... by unique id
  map <- tibble::tibble(orig = unique_ids, deid = paste0('P', stringr::str_pad(seq_along(unique_ids), width=6, pad='0')))
  df <- df %>% left_join(map, by = setNames('orig', id_col))
  df <- df %>% rename(deid_patient_id = deid)
} else {
  # fallback: create per-row pseudonym
  message('No explicit id column found; creating per-row pseudonym')
  df <- df %>% mutate(deid_patient_id = paste0('P', stringr::str_pad(row_number(), width=6, pad='0')))
}

# Remove direct identifiers and potentially identifying free-text columns
drop_cols <- c('merge_key','SAMPLE_ID','Study ID','Sample ID','Sample_ID','Sample coverage','Sample ID.1')
drop_cols <- intersect(drop_cols, names(df))
if (length(drop_cols)>0) df <- df %>% select(-all_of(drop_cols))

# Remove exact age fields while keeping age category
age_exact_cols <- c('Invasive_carcinoma_dx_age','Age at Which Sequencing was Reported (Years)','Patient Current Age')
age_exact_cols <- intersect(age_exact_cols, names(df))
if (length(age_exact_cols)>0){
  df <- df %>% select(-all_of(age_exact_cols))
}

# Cap any remaining age-like numeric columns at 90 (for >=90)
age_like <- names(df)[str_detect(names(df), regex('age|Age|AGE', ignore_case=TRUE))]
for (c in age_like){
  if (is.numeric(df[[c]])){
    df[[c]] <- ifelse(!is.na(df[[c]]) & df[[c]] >= 90, 90, df[[c]])
  }
}

# Remove or redact any columns that are likely PHI: exact dates, contact info
phi_patterns <- c('contact','address','phone','email','name','patient name')
phi_cols <- names(df)[str_detect(tolower(names(df)), paste(phi_patterns, collapse='|'))]
if (length(phi_cols)>0){
  message('Dropping potential PHI columns: ', paste(phi_cols, collapse=', '))
  df <- df %>% select(-all_of(phi_cols))
}

# Reorder to put deid id first and drop internal mapping objects
df_deid <- df %>% select(deid_patient_id, everything())

message('Writing deidentified dataset to: ', outfile)
readr::write_csv(df_deid, outfile)
message('Done. De-identified file created: ', outfile)
