# Shared R helpers for the brain-mets analytic pipeline (steps 1-2 of the
# analysis order of operations: data import and preprocessing).
#
# Sourced by:
#   src/exploratory data analysis/table1_prisma_descriptive.R
#   src/modeling/finegray_cox_risk_models.R
#
# Base R only. Excel input additionally needs `readxl` or `openxlsx`.

# ------------------------------------------------------------------ CLI ------

#' Minimal `--key value` command-line parser.
#' `defaults` is a named list; unknown keys are an error.
parse_args <- function(args, defaults) {
  out <- defaults
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("expected --option, got: ", key)
    name <- gsub("-", "_", substring(key, 3))
    if (!name %in% names(defaults)) stop("unknown argument: ", key)
    if (i == length(args)) stop("missing value for ", key)
    out[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

#' Repository root, inferred from the running script's path.
#' `depth` = directories between the script and the repo root.
project_root <- function(depth = 2) {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 1) {
    # Rscript escapes spaces in --file= as `~+~`.
    self <- gsub("~\\+~", " ", sub("^--file=", "", hit))
    root <- do.call(file.path, c(list(dirname(self)), as.list(rep("..", depth))))
    if (dir.exists(root)) return(normalizePath(root))
  }
  normalizePath(getwd())
}

# -------------------------------------------------------- Step 1: import -----

#' Read a .csv / .tsv / .xlsx analytic frame into a data.frame.
read_any <- function(path, sheet = 1) {
  if (!file.exists(path)) stop("input not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (requireNamespace("readxl", quietly = TRUE)) {
      sh <- suppressWarnings(as.integer(sheet))
      if (is.na(sh)) sh <- sheet
      return(as.data.frame(readxl::read_excel(path, sheet = sh)))
    }
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      return(openxlsx::read.xlsx(path, sheet = sheet))
    }
    stop("reading ", ext, " needs the readxl or openxlsx package")
  }
  sep <- if (ext == "tsv") "\t" else ","
  read.csv(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Default analytic-frame path for a cohort key.
default_data_path <- function(root, cohort) {
  file.path(root, "data", "processed",
            sprintf("extracted_variables_%s_data.csv", cohort))
}

# --------------------------------------------------- Step 2: preprocessing ---

num <- function(x) suppressWarnings(as.numeric(as.character(x)))

#' First of `candidates` present in `df`; NULL (or error) when none are.
first_present <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, names(df))
  if (!length(hit)) {
    if (required) {
      stop("none of these columns are present: ", paste(candidates, collapse = ", "))
    }
    return(NULL)
  }
  hit[[1]]
}

CATEGORICAL_COVARIATES <- c(
  "receptor_primary_cat", "race_clean", "ethnicity_clean", "sex",
  "smoking_status", "stage_dx_cat", "insurance", "grade_ord", "SEQ_ASSAY_ID"
)
NUMERIC_COVARIATES <- c("age_dx_num", "mutation_count", "tmb")
MISSING_TOKENS <- c("", "NA", "N/A", "Unknown", "unknown", "UNKNOWN", "Not Available")

#' Coerce categoricals to factors, derive `survival_time` /
#' `time_to_brain_mets`, and create the `event_dead` / `brain_mets_flag` flags.
#' Only columns present in the cohort are touched.
preprocess <- function(df) {
  for (cc in intersect(CATEGORICAL_COVARIATES, names(df))) {
    v <- as.character(df[[cc]])
    v[v %in% MISSING_TOKENS] <- NA
    df[[cc]] <- factor(v)
  }
  for (nc in intersect(NUMERIC_COVARIATES, names(df))) {
    df[[nc]] <- num(df[[nc]])
  }

  # The harmonized frame carries months-scale durations already; dates are only
  # parsed when a duration column is absent.
  df$survival_time <- num(df[[first_present(df, c("OS_months", "tt_os_dx_mos"))]])
  df$time_to_brain_mets <- num(df[[first_present(
    df, c("tt_brain_met_mos", "time_to_brain_met_mos")
  )]])

  df$event_dead <- as.integer(num(df[[first_present(
    df, c("os_status_bin", "death_event", "DEATH_EVENT")
  )]]) == 1)
  df$brain_mets_flag <- as.integer(num(df[[first_present(
    df, c("any_brain_met", "brain_met_event", "dist_mets_brain_cns")
  )]]) == 1)
  df$event_dead[is.na(df$event_dead)] <- 0L
  df$brain_mets_flag[is.na(df$brain_mets_flag)] <- 0L

  df$brain_met_at_dx <- if ("brain_met_at_dx" %in% names(df)) {
    as.integer(num(df$brain_met_at_dx) == 1)
  } else {
    NA_integer_
  }

  # Brain-met event indicator restricted to the incident (post-dx) endpoint.
  bm_event_col <- first_present(df, c("brain_met_event"), required = FALSE)
  df$brain_met_event <- if (is.null(bm_event_col)) {
    df$brain_mets_flag
  } else {
    ev <- as.integer(num(df[[bm_event_col]]) == 1)
    ev[is.na(ev)] <- 0L
    ev
  }
  df
}

#' Gene-alteration indicator columns in the canonical frame.
gene_columns <- function(df, prefix = "G_top10_") {
  grep(paste0("^", prefix), names(df), value = TRUE)
}

#' Clinical covariates available for regression, dropping levels that are too
#' sparse to fit (< `min_level_n` in any retained level).
usable_covariates <- function(df, candidates, min_level_n = 5) {
  keep <- character(0)
  for (cc in intersect(candidates, names(df))) {
    v <- df[[cc]]
    if (is.factor(v)) {
      v <- droplevels(v)
      if (nlevels(v) < 2) next
      tab <- table(v)
      if (any(tab < min_level_n)) {
        # Collapse rare levels into "Other" rather than dropping the covariate.
        rare <- names(tab)[tab < min_level_n]
        v <- as.character(v)
        v[v %in% rare] <- "Other"
        v <- factor(v)
        if (nlevels(v) < 2) next
        df[[cc]] <- v
      } else {
        df[[cc]] <- v
      }
    } else {
      if (length(unique(v[!is.na(v)])) < 2) next
    }
    keep <- c(keep, cc)
  }
  list(covariates = keep, data = df)
}
