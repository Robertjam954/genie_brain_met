#!/usr/bin/env Rscript
# Step 4 of the analysis order of operations: risk modeling.
#
#   4A Brain metastasis (competing risk = death)
#      - cumulative incidence (Aalen-Johansen), overall and by receptor subtype
#      - multivariable Fine-Gray subdistribution hazards model (clinical covariates)
#      - per-gene subdistribution hazards, each adjusted for the clinical covariates,
#        with Benjamini-Hochberg q-values
#   4B Overall survival
#      - Kaplan-Meier by receptor subtype + log-rank test
#      - multivariable Cox proportional hazards model on the same covariates
#      - per-gene adjusted Cox hazards with BH q-values
#
# The Fine-Gray fit uses survival::finegray() to build the risk-weighted data set
# and coxph(..., weights = fgwt) for the subdistribution hazards, so no package
# beyond the base-recommended `survival` is required.
#
# Usage
#   Rscript "src/modeling/finegray_cox_risk_models.R" [--cohort genie] [--data PATH]
#           [--sheet 1] [--outdir DIR] [--figdir DIR] [--min-events 25]
#           [--strata receptor_primary_cat]
#
# Outputs (in --outdir, default src/modeling/<cohort>/risk_models/)
#   cif_brain_met_overall.csv        cumulative incidence, brain met vs death
#   cif_brain_met_by_<strata>.csv    cumulative incidence by stratum
#   finegray_multivariable.csv       clinical-covariate subhazard ratios + 95% CI
#   finegray_gene_subhazards.csv     per-gene sHR, 95% CI, p, BH q
#   km_os_by_<strata>.csv            KM estimates by stratum
#   logrank_os_by_<strata>.txt       log-rank chi-square / df / p
#   cox_os_multivariable.csv         clinical-covariate HR + 95% CI (+ PH test)
#   cox_os_gene_hazards.csv          per-gene HR, 95% CI, p, BH q
#   risk_model_summary.txt           cohort sizes, event counts, model status
# Figures (in --figdir, default "manuscript components/<cohort>/risk_models/")
#   figure2_cif_brain_met.png, figure3_km_os.png

suppressWarnings(suppressPackageStartupMessages(library(survival)))

.self_dir <- local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a) == 1) dirname(gsub("~\\+~", " ", sub("^--file=", "", a))) else getwd()
})
source(file.path(.self_dir, "_lib.R"))

CLINICAL_COVARIATES <- c(
  "age_dx_num", "receptor_primary_cat", "stage_dx_cat", "grade_ord",
  "smoking_status", "insurance", "race_clean", "ethnicity_clean", "sex"
)

# ------------------------------------------------ competing-risk endpoint -----

#' Build the multi-state brain-metastasis endpoint with death as the competing
#' risk: time = first of (brain met, death, censoring);
#' status = 0 censored / 1 brain met / 2 death without brain met.
build_competing_endpoint <- function(df) {
  t_bm <- df$time_to_brain_mets
  t_os <- df$survival_time
  bm <- df$brain_met_event
  dead <- df$event_dead

  time <- ifelse(bm == 1, t_bm, t_os)
  # Fall back to whichever duration is available.
  time <- ifelse(is.na(time), ifelse(is.na(t_os), t_bm, t_os), time)
  status <- ifelse(bm == 1, 1L, ifelse(dead == 1, 2L, 0L))

  df$cr_time <- as.numeric(time)
  df$cr_status <- factor(status, levels = c(0, 1, 2),
                         labels = c("censored", "brain_met", "death"))
  df
}

#' Rows usable for a competing-risk / survival fit: positive finite time.
usable_rows <- function(df, time_col) {
  ok <- !is.na(df[[time_col]]) & df[[time_col]] > 0
  df[ok, , drop = FALSE]
}

# ------------------------------------------------------ 4A cumulative inc -----

write_cif <- function(fit, path, stratum_label = NA_character_) {
  s <- summary(fit)
  # survfit on a multi-state factor returns pstate columns per event type.
  out <- data.frame(time = fit$time, stringsAsFactors = FALSE)
  states <- fit$states
  pstate <- fit$pstate
  if (is.null(dim(pstate))) pstate <- matrix(pstate, ncol = 1)
  if (length(dim(pstate)) == 3) {
    # stratified fits are flattened by the caller, so this should not happen
    stop("stratified survfit passed to write_cif")
  }
  for (j in seq_along(states)) {
    if (states[[j]] == "(s0)") next
    out[[paste0("cif_", states[[j]])]] <- pstate[, j]
  }
  if (!is.na(stratum_label)) out$stratum <- stratum_label
  utils::write.csv(out, path, row.names = FALSE)
  invisible(out)
}

cif_overall <- function(df, outdir) {
  fit <- survfit(Surv(cr_time, cr_status) ~ 1, data = df)
  write_cif(fit, file.path(outdir, "cif_brain_met_overall.csv"))
  fit
}

cif_by_strata <- function(df, strata, outdir) {
  parts <- list()
  for (lev in levels(droplevels(factor(df[[strata]])))) {
    sub <- df[!is.na(df[[strata]]) & df[[strata]] == lev, , drop = FALSE]
    if (nrow(sub) < 10 || sum(sub$cr_status == "brain_met") < 3) {
      message(sprintf("  skip CIF stratum %s (n=%d, events=%d)", lev, nrow(sub),
                      sum(sub$cr_status == "brain_met")))
      next
    }
    fit <- survfit(Surv(cr_time, cr_status) ~ 1, data = sub)
    out <- data.frame(time = fit$time, stringsAsFactors = FALSE)
    for (j in seq_along(fit$states)) {
      if (fit$states[[j]] == "(s0)") next
      out[[paste0("cif_", fit$states[[j]])]] <- fit$pstate[, j]
    }
    out$stratum <- lev
    out$n <- nrow(sub)
    parts[[lev]] <- out
  }
  if (!length(parts)) return(NULL)
  res <- do.call(rbind, parts)
  utils::write.csv(res, file.path(outdir, sprintf("cif_brain_met_by_%s.csv", strata)),
                   row.names = FALSE)
  res
}

plot_cif <- function(df, strata, path, cohort) {
  grDevices::png(path, width = 1200, height = 800, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  levs <- if (!is.null(strata) && strata %in% names(df)) {
    levels(droplevels(factor(df[[strata]])))
  } else {
    character(0)
  }
  cols <- grDevices::hcl.colors(max(length(levs), 1), "Dark 3")
  graphics::par(mar = c(4.5, 4.5, 3, 1))
  drawn <- FALSE
  for (i in seq_along(levs)) {
    sub <- df[!is.na(df[[strata]]) & df[[strata]] == levs[[i]], , drop = FALSE]
    if (nrow(sub) < 10 || sum(sub$cr_status == "brain_met") < 3) next
    fit <- survfit(Surv(cr_time, cr_status) ~ 1, data = sub)
    j <- which(fit$states == "brain_met")
    if (!length(j)) next
    if (!drawn) {
      graphics::plot(fit$time, fit$pstate[, j], type = "s", col = cols[[i]], lwd = 2,
                     xlab = "Months from diagnosis",
                     ylab = "Cumulative incidence of brain metastasis",
                     ylim = c(0, 1),
                     main = sprintf("Figure 2. Cumulative incidence of brain metastasis\n(death as competing risk) - %s", cohort))
      drawn <- TRUE
    } else {
      graphics::lines(fit$time, fit$pstate[, j], type = "s", col = cols[[i]], lwd = 2)
    }
  }
  if (!drawn) {
    fit <- survfit(Surv(cr_time, cr_status) ~ 1, data = df)
    j <- which(fit$states == "brain_met")
    graphics::plot(fit$time, fit$pstate[, j], type = "s", lwd = 2,
                   xlab = "Months from diagnosis",
                   ylab = "Cumulative incidence of brain metastasis", ylim = c(0, 1),
                   main = sprintf("Figure 2. Cumulative incidence of brain metastasis\n(death as competing risk) - %s", cohort))
  } else {
    graphics::legend("topleft", legend = levs, col = cols, lwd = 2, bty = "n", cex = 0.85)
  }
  invisible(path)
}

# ------------------------------------------------------ 4A Fine-Gray ----------

#' Tidy a coxph fit into a term/estimate table. `label` names the ratio column
#' ("subhazard_ratio" for Fine-Gray, "hazard_ratio" for Cox).
tidy_cox <- function(fit, label) {
  s <- summary(fit)
  co <- s$coefficients
  ci <- s$conf.int
  out <- data.frame(
    term = rownames(co),
    ratio = as.numeric(co[, "exp(coef)"]),
    ci_low = as.numeric(ci[, 3]),
    ci_high = as.numeric(ci[, 4]),
    se_log_ratio = as.numeric(co[, "se(coef)"]),
    z = as.numeric(co[, ncol(co) - 1]),
    p_value = as.numeric(co[, ncol(co)]),
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "ratio"] <- label
  names(out)[names(out) == "ci_low"] <- paste0(label, "_ci_low")
  names(out)[names(out) == "ci_high"] <- paste0(label, "_ci_high")
  out
}

#' Fine-Gray subdistribution hazards for the brain-met cause.
#' Returns NULL when the model cannot be fit.
fit_finegray <- function(df, covariates, etype = "brain_met") {
  keep <- c("cr_time", "cr_status", covariates)
  sub <- stats::na.omit(df[, keep, drop = FALSE])
  if (!nrow(sub) || sum(sub$cr_status == etype) < 5) return(NULL)
  sub$cr_status <- droplevels(sub$cr_status)
  if (!etype %in% levels(sub$cr_status)) return(NULL)
  fg <- tryCatch(
    survival::finegray(Surv(cr_time, cr_status) ~ ., data = sub, etype = etype),
    error = function(e) { message("  finegray() failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(fg)) return(NULL)
  form <- stats::as.formula(paste(
    "Surv(fgstart, fgstop, fgstatus) ~", paste(sprintf("`%s`", covariates), collapse = " + ")
  ))
  tryCatch(
    coxph(form, data = fg, weights = fgwt),
    error = function(e) { message("  Fine-Gray coxph failed: ", conditionMessage(e)); NULL }
  )
}

#' Per-gene models: each gene entered alone alongside `covariates`.
#' `fitter` returns a coxph fit (or NULL); `label` names the ratio column.
per_gene_table <- function(df, genes, covariates, fitter, label) {
  rows <- list()
  for (g in genes) {
    v <- suppressWarnings(as.numeric(df[[g]]))
    n_alt <- sum(v == 1, na.rm = TRUE)
    if (n_alt < 5 || n_alt > sum(!is.na(v)) - 5) {
      rows[[g]] <- data.frame(gene = g, n_altered = n_alt, ratio = NA_real_,
                              ci_low = NA_real_, ci_high = NA_real_,
                              p_value = NA_real_, note = "too few altered/unaltered",
                              stringsAsFactors = FALSE)
      next
    }
    d <- df
    d[[g]] <- v
    fit <- fitter(d, c(covariates, g))
    if (is.null(fit)) {
      rows[[g]] <- data.frame(gene = g, n_altered = n_alt, ratio = NA_real_,
                              ci_low = NA_real_, ci_high = NA_real_,
                              p_value = NA_real_, note = "model did not converge",
                              stringsAsFactors = FALSE)
      next
    }
    tid <- tidy_cox(fit, "ratio")
    hit <- tid[tid$term %in% c(g, sprintf("`%s`", g)), , drop = FALSE]
    if (!nrow(hit)) {
      rows[[g]] <- data.frame(gene = g, n_altered = n_alt, ratio = NA_real_,
                              ci_low = NA_real_, ci_high = NA_real_,
                              p_value = NA_real_, note = "gene term dropped",
                              stringsAsFactors = FALSE)
      next
    }
    rows[[g]] <- data.frame(gene = g, n_altered = n_alt, ratio = hit$ratio[[1]],
                            ci_low = hit$ratio_ci_low[[1]], ci_high = hit$ratio_ci_high[[1]],
                            p_value = hit$p_value[[1]], note = "",
                            stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$q_value <- NA_real_
  ok <- !is.na(out$p_value)
  if (any(ok)) out$q_value[ok] <- stats::p.adjust(out$p_value[ok], method = "BH")
  names(out)[names(out) == "ratio"] <- label
  names(out)[names(out) == "ci_low"] <- paste0(label, "_ci_low")
  names(out)[names(out) == "ci_high"] <- paste0(label, "_ci_high")
  out[order(out$p_value, na.last = TRUE), , drop = FALSE]
}

# ------------------------------------------------------ 4B OS models ----------

km_by_strata <- function(df, strata, outdir, figdir, cohort) {
  d <- usable_rows(df, "survival_time")
  d <- d[!is.na(d[[strata]]), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d[[strata]] <- droplevels(factor(d[[strata]]))
  form <- stats::as.formula(sprintf("Surv(survival_time, event_dead) ~ `%s`", strata))
  fit <- survfit(form, data = d)

  s <- summary(fit)
  est <- data.frame(
    stratum = if (is.null(s$strata)) strata else as.character(s$strata),
    time = s$time, n_risk = s$n.risk, n_event = s$n.event,
    survival = s$surv, ci_low = s$lower, ci_high = s$upper,
    stringsAsFactors = FALSE
  )
  utils::write.csv(est, file.path(outdir, sprintf("km_os_by_%s.csv", strata)),
                   row.names = FALSE)

  lr_txt <- "log-rank not computed (fewer than 2 strata with events)"
  if (nlevels(d[[strata]]) >= 2) {
    lr <- tryCatch(survdiff(form, data = d), error = function(e) NULL)
    if (!is.null(lr)) {
      p <- stats::pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
      lr_txt <- sprintf("log-rank chisq=%.4f  df=%d  p=%.4g", lr$chisq,
                        length(lr$n) - 1, p)
    }
  }
  writeLines(lr_txt, file.path(outdir, sprintf("logrank_os_by_%s.txt", strata)))

  path <- file.path(figdir, "figure3_km_os.png")
  grDevices::png(path, width = 1200, height = 800, res = 140)
  cols <- grDevices::hcl.colors(nlevels(d[[strata]]), "Dark 3")
  graphics::par(mar = c(4.5, 4.5, 3.5, 1))
  graphics::plot(fit, col = cols, lwd = 2, xlab = "Months from diagnosis",
                 ylab = "Overall survival probability",
                 main = sprintf("Figure 3. Overall survival by %s - %s\n%s",
                                strata, cohort, lr_txt))
  graphics::legend("bottomleft", legend = levels(d[[strata]]), col = cols, lwd = 2,
                   bty = "n", cex = 0.85)
  grDevices::dev.off()
  list(fit = fit, logrank = lr_txt, figure = path)
}

fit_cox_os <- function(df, covariates) {
  keep <- c("survival_time", "event_dead", covariates)
  sub <- stats::na.omit(df[, keep, drop = FALSE])
  sub <- sub[sub$survival_time > 0, , drop = FALSE]
  if (!nrow(sub) || sum(sub$event_dead) < 5) return(NULL)
  form <- stats::as.formula(paste(
    "Surv(survival_time, event_dead) ~",
    paste(sprintf("`%s`", covariates), collapse = " + ")
  ))
  tryCatch(coxph(form, data = sub), error = function(e) {
    message("  Cox PH failed: ", conditionMessage(e)); NULL
  })
}

# ------------------------------------------------------------- main ----------

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE), list(
    data = NULL, sheet = 1, cohort = "genie", outdir = NULL, figdir = NULL,
    min_events = 25, strata = "receptor_primary_cat"
  ))
  min_events <- as.integer(args$min_events)
  root <- project_root(depth = 2)
  data_path <- if (is.null(args$data)) default_data_path(root, args$cohort) else args$data
  outdir <- if (is.null(args$outdir)) {
    file.path(root, "src", "modeling", args$cohort, "risk_models")
  } else args$outdir
  figdir <- if (is.null(args$figdir)) {
    file.path(root, "manuscript components", args$cohort, "risk_models")
  } else args$figdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("[import]  %s\n", data_path))
  df <- preprocess(read_any(data_path, args$sheet))
  cat(sprintf("          %d rows\n", nrow(df)))

  # Aim-2 style cohort: exclude CNS involvement at diagnosis when known.
  incident <- df[is.na(df$brain_met_at_dx) | df$brain_met_at_dx == 0, , drop = FALSE]
  incident <- build_competing_endpoint(incident)
  incident <- usable_rows(incident, "cr_time")

  uc <- usable_covariates(incident, CLINICAL_COVARIATES)
  incident <- uc$data
  covariates <- uc$covariates
  genes <- gene_columns(incident)
  n_bm <- sum(incident$cr_status == "brain_met")
  n_death <- sum(incident$cr_status == "death")
  cat(sprintf("[4A]      competing-risk cohort n=%d  brain-met=%d  death=%d\n",
              nrow(incident), n_bm, n_death))
  cat(sprintf("          covariates: %s\n", paste(covariates, collapse = ", ")))
  cat(sprintf("          genes: %d\n", length(genes)))

  summary_lines <- c(
    sprintf("cohort: %s", args$cohort),
    sprintf("input: %s", data_path),
    sprintf("competing-risk cohort n=%d (brain-met=%d, death=%d)", nrow(incident),
            n_bm, n_death),
    sprintf("clinical covariates: %s", paste(covariates, collapse = ", ")),
    sprintf("gene columns: %d", length(genes))
  )

  # --- 4A cumulative incidence ---
  cif_overall(incident, outdir)
  strata <- args$strata
  if (strata %in% names(incident)) {
    cif_by_strata(incident, strata, outdir)
    plot_cif(incident, strata, file.path(figdir, "figure2_cif_brain_met.png"), args$cohort)
  } else {
    plot_cif(incident, NULL, file.path(figdir, "figure2_cif_brain_met.png"), args$cohort)
  }
  cat("          wrote cumulative-incidence tables + figure 2\n")

  # --- 4A Fine-Gray ---
  if (n_bm < min_events) {
    msg <- sprintf("SKIP Fine-Gray: %d brain-met events (< %d)", n_bm, min_events)
    cat("         ", msg, "\n")
    summary_lines <- c(summary_lines, msg)
  } else {
    fg_fit <- fit_finegray(incident, covariates)
    if (is.null(fg_fit)) {
      summary_lines <- c(summary_lines, "Fine-Gray multivariable model did not fit")
    } else {
      tid <- tidy_cox(fg_fit, "subhazard_ratio")
      utils::write.csv(tid, file.path(outdir, "finegray_multivariable.csv"),
                       row.names = FALSE)
      cat(sprintf("          wrote finegray_multivariable.csv (%d terms)\n", nrow(tid)))
    }
    if (length(genes)) {
      gene_tab <- per_gene_table(
        incident, genes, covariates,
        fitter = function(d, cov) fit_finegray(d, cov),
        label = "subhazard_ratio"
      )
      utils::write.csv(gene_tab, file.path(outdir, "finegray_gene_subhazards.csv"),
                       row.names = FALSE)
      cat(sprintf("          wrote finegray_gene_subhazards.csv (%d genes, %d with q<0.05)\n",
                  nrow(gene_tab), sum(gene_tab$q_value < 0.05, na.rm = TRUE)))
    }
  }

  # --- 4B overall survival ---
  os <- usable_rows(df, "survival_time")
  uc_os <- usable_covariates(os, CLINICAL_COVARIATES)
  os <- uc_os$data
  cov_os <- uc_os$covariates
  n_dead <- sum(os$event_dead == 1, na.rm = TRUE)
  cat(sprintf("[4B]      OS cohort n=%d  deaths=%d\n", nrow(os), n_dead))
  summary_lines <- c(summary_lines,
                     sprintf("OS cohort n=%d (deaths=%d)", nrow(os), n_dead))

  if (strata %in% names(os)) {
    km <- km_by_strata(os, strata, outdir, figdir, args$cohort)
    if (!is.null(km)) {
      cat(sprintf("          %s\n", km$logrank))
      summary_lines <- c(summary_lines, km$logrank)
    }
  }

  if (n_dead < min_events) {
    msg <- sprintf("SKIP Cox PH: %d deaths (< %d)", n_dead, min_events)
    cat("         ", msg, "\n")
    summary_lines <- c(summary_lines, msg)
  } else {
    cox_fit <- fit_cox_os(os, cov_os)
    if (is.null(cox_fit)) {
      summary_lines <- c(summary_lines, "Cox multivariable model did not fit")
    } else {
      tid <- tidy_cox(cox_fit, "hazard_ratio")
      # Scaled Schoenfeld residual test for the PH assumption.
      zph <- tryCatch(cox.zph(cox_fit), error = function(e) NULL)
      if (!is.null(zph)) {
        ztab <- as.data.frame(zph$table)
        ztab$term <- rownames(ztab)
        tid$ph_test_p <- ztab$p[match(tid$term, ztab$term)]
        summary_lines <- c(summary_lines, sprintf(
          "global PH test p=%.4g", ztab$p[rownames(ztab) == "GLOBAL"][[1]]
        ))
      }
      utils::write.csv(tid, file.path(outdir, "cox_os_multivariable.csv"),
                       row.names = FALSE)
      cat(sprintf("          wrote cox_os_multivariable.csv (%d terms)\n", nrow(tid)))
    }
    genes_os <- gene_columns(os)
    if (length(genes_os)) {
      gene_tab <- per_gene_table(
        os, genes_os, cov_os,
        fitter = function(d, cov) fit_cox_os(d, cov),
        label = "hazard_ratio"
      )
      utils::write.csv(gene_tab, file.path(outdir, "cox_os_gene_hazards.csv"),
                       row.names = FALSE)
      cat(sprintf("          wrote cox_os_gene_hazards.csv (%d genes, %d with p<0.05 & HR>1)\n",
                  nrow(gene_tab),
                  sum(gene_tab$p_value < 0.05 & gene_tab$hazard_ratio > 1, na.rm = TRUE)))
    }
  }

  writeLines(summary_lines, file.path(outdir, "risk_model_summary.txt"))
  cat(sprintf("[done]    outputs in %s\n", outdir))
}

if (sys.nframe() == 0L) main()
