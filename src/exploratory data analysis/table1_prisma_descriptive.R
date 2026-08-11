#!/usr/bin/env Rscript
# Steps 1-3 of the analysis order of operations: data import, preprocessing, and
# descriptive analysis.
#
# Inputs
#   --data   analytic frame (.csv or .xlsx). Default:
#            data/processed/extracted_variables_genie_data.csv
#   --sheet  worksheet name/index when --data is an .xlsx (default 1)
#   --cohort cohort key used for output paths (default genie)
#   --outdir output root (default src/modeling/<cohort>/descriptive)
#   --figdir figure root (default manuscript components/<cohort>/descriptive)
#
# Outputs
#   table1_overall.csv            Table 1, whole cohort
#   table1_by_brain_met.csv       Table 1 stratified by brain-met status + p-values
#   figure1_prisma_flowchart.png  PRISMA-style cohort flow diagram
#   cohort_flow_counts.csv        the counts behind the flow diagram
#   preprocessed_summary.csv      audit of the derived factors/flags/TTE columns
#
# Base R + `survival` only (no tidyverse dependency) so the script runs on a
# stock R install. Excel input needs `readxl` or `openxlsx`.
suppressWarnings(suppressPackageStartupMessages(library(survival)))

# Shared import/preprocessing helpers (steps 1-2).
.self_dir <- local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a) == 1) dirname(gsub("~\\+~", " ", sub("^--file=", "", a))) else getwd()
})
source(file.path(.self_dir, "..", "modeling", "_lib.R"), chdir = FALSE)


# ---------------------------------------------- Step 3: descriptive ----------

#' Summarize one variable, overall and by `group` (NULL = overall only).
#' Numerics report median [IQR]; factors report n (%) per level.
summarize_var <- function(df, var, group = NULL) {
  x <- df[[var]]
  rows <- list()
  add <- function(label, level, values, p = NA_real_) {
    rows[[length(rows) + 1]] <<- c(
      list(variable = label, level = level), values, list(p_value = p)
    )
  }
  groups <- if (is.null(group)) list(Overall = seq_len(nrow(df))) else split(
    seq_len(nrow(df)), factor(df[[group]])
  )

  if (is.numeric(x)) {
    vals <- lapply(groups, function(idx) {
      v <- x[idx]
      v <- v[!is.na(v)]
      if (!length(v)) return("-")
      sprintf("%.1f [%.1f-%.1f]", median(v), quantile(v, .25), quantile(v, .75))
    })
    p <- NA_real_
    if (!is.null(group) && length(groups) == 2) {
      ok <- tryCatch(
        wilcox.test(x ~ factor(df[[group]]))$p.value, error = function(e) NA_real_
      )
      p <- ok
    }
    add(paste0(var, ", median [IQR]"), "", vals, p)
  } else {
    fx <- factor(x)
    p <- NA_real_
    if (!is.null(group) && nlevels(fx) >= 2) {
      tab <- table(fx, factor(df[[group]]))
      p <- tryCatch(
        if (any(suppressWarnings(chisq.test(tab)$expected) < 5)) {
          fisher.test(tab, simulate.p.value = TRUE)$p.value
        } else {
          chisq.test(tab)$p.value
        },
        error = function(e) NA_real_
      )
    }
    first <- TRUE
    for (lev in levels(fx)) {
      vals <- lapply(groups, function(idx) {
        v <- fx[idx]
        denom <- sum(!is.na(v))
        n <- sum(v == lev, na.rm = TRUE)
        if (!denom) return("-")
        sprintf("%d (%.1f%%)", n, 100 * n / denom)
      })
      add(paste0(var, ", n (%)"), lev, vals, if (first) p else NA_real_)
      first <- FALSE
    }
    n_miss <- sum(is.na(fx))
    if (n_miss > 0) {
      vals <- lapply(groups, function(idx) as.character(sum(is.na(fx[idx]))))
      add(paste0(var, ", n (%)"), "(missing)", vals, NA_real_)
    }
  }
  do.call(rbind, lapply(rows, function(r) {
    as.data.frame(r, stringsAsFactors = FALSE, check.names = FALSE)
  }))
}

TABLE1_VARS <- c(
  "age_dx_num", "sex", "smoking_status", "stage_dx_cat",
  "receptor_primary_cat", "insurance", "race_clean", "ethnicity_clean",
  "grade_ord", "mutation_count"
)

build_table1 <- function(df, group = NULL) {
  vars <- intersect(TABLE1_VARS, names(df))
  missing <- setdiff(TABLE1_VARS, vars)
  if (length(missing)) {
    message("  note: Table 1 variables absent from this cohort: ",
            paste(missing, collapse = ", "))
  }
  parts <- lapply(vars, function(v) summarize_var(df, v, group))
  out <- do.call(rbind, parts)
  # Header row with group sizes.
  if (is.null(group)) {
    hdr <- data.frame(variable = "N", level = "", Overall = as.character(nrow(df)),
                      p_value = NA_real_, stringsAsFactors = FALSE,
                      check.names = FALSE)
  } else {
    counts <- table(factor(df[[group]]))
    hdr <- data.frame(variable = "N", level = "", stringsAsFactors = FALSE,
                      check.names = FALSE)
    for (g in names(counts)) hdr[[g]] <- as.character(counts[[g]])
    hdr$p_value <- NA_real_
  }
  out <- out[, c("variable", "level",
                 setdiff(names(out), c("variable", "level", "p_value")),
                 "p_value")]
  hdr <- hdr[, names(out)]
  rbind(hdr, out)
}

# ------------------------------------------- Step 3b: PRISMA flow chart ------

cohort_flow <- function(df) {
  n_all <- nrow(df)
  n_time <- sum(!is.na(df$survival_time))
  n_bm_ever <- sum(df$brain_mets_flag == 1, na.rm = TRUE)
  n_bm_at_dx <- sum(df$brain_met_at_dx == 1, na.rm = TRUE)
  n_no_cns_dx <- sum(is.na(df$brain_met_at_dx) | df$brain_met_at_dx == 0)
  n_tt_bm <- sum(!is.na(df$time_to_brain_mets))
  data.frame(
    step = c(
      "Records in analytic frame",
      "Excluded: no survival time",
      "Analysis cohort (survival time available)",
      "Brain metastasis ever",
      "Brain metastasis at diagnosis",
      "Time-to-brain-met cohort (no CNS at dx)",
      "With non-missing time to brain met"
    ),
    n = c(n_all, n_all - n_time, n_time, n_bm_ever, n_bm_at_dx,
          n_no_cns_dx, n_tt_bm),
    stringsAsFactors = FALSE
  )
}

draw_prisma <- function(flow, path, cohort) {
  boxes <- list(
    list(y = 0.90, txt = sprintf("Records in analytic frame\nn = %d",
                                 flow$n[flow$step == "Records in analytic frame"])),
    list(y = 0.66, txt = sprintf("Analysis cohort\n(survival time available)\nn = %d",
                                 flow$n[flow$step == "Analysis cohort (survival time available)"])),
    list(y = 0.40, txt = sprintf("Time-to-brain-met cohort\n(no CNS at diagnosis)\nn = %d",
                                 flow$n[flow$step == "Time-to-brain-met cohort (no CNS at dx)"])),
    list(y = 0.14, txt = sprintf("Brain metastasis ever\nn = %d",
                                 flow$n[flow$step == "Brain metastasis ever"]))
  )
  excl <- list(
    list(y = 0.78, txt = sprintf("Excluded: no survival time\nn = %d",
                                 flow$n[flow$step == "Excluded: no survival time"])),
    list(y = 0.53, txt = sprintf("Brain metastasis at diagnosis\nn = %d",
                                 flow$n[flow$step == "Brain metastasis at diagnosis"]))
  )

  grDevices::png(path, width = 1100, height = 1200, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(0, 0, 2, 0))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::title(main = sprintf("Figure 1. Cohort flow - %s", cohort))

  bx <- 0.34
  for (b in boxes) {
    graphics::rect(bx - 0.22, b$y - 0.065, bx + 0.22, b$y + 0.065,
                   border = "black", col = "grey97", lwd = 1.4)
    graphics::text(bx, b$y, b$txt, cex = 0.8)
  }
  for (e in excl) {
    graphics::rect(0.66, e$y - 0.055, 0.99, e$y + 0.055,
                   border = "grey40", col = "white", lty = 2)
    graphics::text(0.825, e$y, e$txt, cex = 0.72, col = "grey20")
  }
  # Vertical arrows down the main spine.
  for (i in seq_len(length(boxes) - 1)) {
    graphics::arrows(bx, boxes[[i]]$y - 0.065, bx, boxes[[i + 1]]$y + 0.065,
                     length = 0.10, lwd = 1.4)
  }
  # Horizontal connectors to the exclusion boxes.
  for (e in excl) {
    graphics::segments(bx, e$y, 0.66, e$y, lty = 2, col = "grey40")
  }
  invisible(path)
}

# ------------------------------------------------------------- main ----------

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE), list(
    data = NULL, sheet = 1, cohort = "genie", outdir = NULL, figdir = NULL
  ))
  root <- project_root(depth = 2)
  data_path <- args$data
  if (is.null(data_path)) data_path <- default_data_path(root, args$cohort)
  outdir <- args$outdir
  if (is.null(outdir)) outdir <- file.path(root, "src", "modeling", args$cohort, "descriptive")
  figdir <- args$figdir
  if (is.null(figdir)) figdir <- file.path(root, "manuscript components", args$cohort, "descriptive")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("[1/3] import: %s\n", data_path))
  df <- read_any(data_path, args$sheet)
  cat(sprintf("      %d rows x %d cols\n", nrow(df), ncol(df)))

  cat("[2/3] preprocess: factors, time-to-event, event flags\n")
  df <- preprocess(df)
  audit <- data.frame(
    column = c("survival_time", "time_to_brain_mets", "event_dead", "brain_mets_flag"),
    n_non_missing = c(sum(!is.na(df$survival_time)), sum(!is.na(df$time_to_brain_mets)),
                      sum(!is.na(df$event_dead)), sum(!is.na(df$brain_mets_flag))),
    n_events_or_median = c(
      sprintf("%.1f", stats::median(df$survival_time, na.rm = TRUE)),
      sprintf("%.1f", stats::median(df$time_to_brain_mets, na.rm = TRUE)),
      as.character(sum(df$event_dead == 1, na.rm = TRUE)),
      as.character(sum(df$brain_mets_flag == 1, na.rm = TRUE))
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(audit, file.path(outdir, "preprocessed_summary.csv"), row.names = FALSE)
  print(audit)

  cat("[3/3] descriptive: Table 1 + PRISMA-style flow chart\n")
  t1 <- build_table1(df, group = NULL)
  utils::write.csv(t1, file.path(outdir, "table1_overall.csv"), row.names = FALSE)
  df$brain_met_group <- factor(ifelse(df$brain_mets_flag == 1, "Brain met", "No brain met"),
                               levels = c("No brain met", "Brain met"))
  t1g <- build_table1(df, group = "brain_met_group")
  utils::write.csv(t1g, file.path(outdir, "table1_by_brain_met.csv"), row.names = FALSE)

  flow <- cohort_flow(df)
  utils::write.csv(flow, file.path(outdir, "cohort_flow_counts.csv"), row.names = FALSE)
  fig <- file.path(figdir, "figure1_prisma_flowchart.png")
  draw_prisma(flow, fig, args$cohort)

  cat(sprintf("      wrote %s\n", file.path(outdir, "table1_overall.csv")))
  cat(sprintf("      wrote %s\n", file.path(outdir, "table1_by_brain_met.csv")))
  cat(sprintf("      wrote %s\n", file.path(outdir, "cohort_flow_counts.csv")))
  cat(sprintf("      wrote %s\n", fig))
}

if (sys.nframe() == 0L) main()
