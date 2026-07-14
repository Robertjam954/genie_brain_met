#!/usr/bin/env Rscript
## Generate Table 1 (demographics) and Table 2 (classical + non-classical predictors)
## Writes CSVs and basic stats to the Output folder

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(broom)
  library(purrr)
  library(rlang)
})

outdir <- file.path("C:", "Users", "jamesr4", "OneDrive - Memorial Sloan Kettering Cancer Center", "Documents", "Research", "Projects", "impact", "Output")
infile <- file.path(outdir, 'merged_impact_patient_level.csv')
if (!file.exists(infile)) stop('Merged input not found: ', infile)
message('Reading: ', infile)
df <- read_csv(infile, show_col_types = FALSE)

# Determine DFS binary column
if ('DFS_EVENT' %in% names(df) && is.numeric(df$DFS_EVENT)) {
  df <- df %>% mutate(.dfs = as.integer(DFS_EVENT))
} else if ('DFS_EVENT_f' %in% names(df)) {
  df <- df %>% mutate(.dfs = as.integer(DFS_EVENT_f == 'Event' | DFS_EVENT_f == 'Event' ))
} else if ('DFS_EVENT' %in% names(df)) {
  # attempt coercion
  df <- df %>% mutate(.dfs = as.integer(as.character(DFS_EVENT) %in% c('1','Event','Eventual')))
} else stop('No DFS_EVENT or DFS_EVENT_f found in data')

df <- df %>% mutate(.dfs = ifelse(is.na(.dfs), NA_integer_, .dfs))

# Variables requested
dem_vars <- c('age_cat','Race_clean','Ethnicity_clean','sample_site','sample_type_bin')
classical <- c('stage_diag_group','overall_tumor_grade_ord','receptor_primary_cat')
possible_mut_vars <- c('tmb_quartile','fga_quartile','MSI_type','mutation_count_q','tmb_nonsynonmous_q','t_alt_count_q','mutation_count_all_sites_sum','tmb_val','t_alt_count')
gene_flags <- grep('^G__', names(df), value = TRUE)

mut_vars_present <- intersect(possible_mut_vars, names(df))
non_classical <- c(mut_vars_present, gene_flags)

message('Detected gene flags: ', paste(head(gene_flags,10), collapse=', '))

# create top_5_gene_status: 1 if any G__ flag is 1
if (length(gene_flags) > 0) {
  df <- df %>% mutate(top_5_gene_status = as.integer((rowSums(across(all_of(gene_flags), ~ as.numeric(coalesce(.x, 0)))) ) > 0))
  message('Created top_5_gene_status (1 if any top5 gene mutated)')
} else {
  df <- df %>% mutate(top_5_gene_status = NA_integer_)
}

make_table <- function(data, vars, dfs_col = '.dfs'){
  out_rows <- list()
  mode_fun <- function(x){
    x2 <- x[!is.na(x)]
    if(length(x2)==0) return(NA_character_)
    ux <- unique(x2)
    ux[which.max(tabulate(match(x2, ux)))]
  }
  for (v in vars){
    if (!v %in% names(data)) next
    vec <- data[[v]]
    # treat ordered factors and numeric as continuous
    if (is.numeric(vec) || is.integer(vec)){
      vals_all_raw <- as.numeric(data %>% pull(!!sym(v)))
      vals_all <- vals_all_raw[!is.na(vals_all_raw)]
      is_normal <- FALSE
      if (length(vals_all) >= 3 && length(vals_all) <= 5000) {
        st <- tryCatch(shapiro.test(vals_all), error = function(e) NULL)
        if (!is.null(st) && !is.na(st$p.value) && st$p.value > 0.05) is_normal <- TRUE
      }
      for (g in c(0,1)){
        vals <- data %>% filter(!is.na(.data[[dfs_col]]), .data[[dfs_col]]==g) %>% pull(!!sym(v))
        vals <- as.numeric(vals)
        vals <- vals[!is.na(vals)]
        if (length(vals)==0) {
          val <- ''
        } else {
          if (is_normal) {
            val <- sprintf('%0.2f (SD %0.2f)', mean(vals), sd(vals))
          } else {
            # non-normal: median (IQR) and mode
            md <- median(vals)
            iqrv <- IQR(vals)
            mv <- mode_fun(round(vals,2))
            val <- sprintf('%0.2f (IQR %0.2f); mode=%s', md, iqrv, mv)
          }
        }
        out_rows[[length(out_rows)+1]] <- tibble::tibble(variable=v, level='', group=as.character(g), value=val)
      }
      vals_all <- vals_all
      allval <- ''
      if (length(vals_all)>0) {
        if (is_normal) allval <- sprintf('%0.2f (SD %0.2f)', mean(vals_all), sd(vals_all)) else allval <- sprintf('%0.2f (IQR %0.2f); mode=%s', median(vals_all), IQR(vals_all), mode_fun(round(vals_all,2)))
      }
      out_rows[[length(out_rows)+1]] <- tibble::tibble(variable=v, level='', group='Total', value=allval)
    } else {
      # categorical
      levs <- sort(unique(as.character(vec[!is.na(vec) & vec != ''])))
      if (length(levs)==0) next
      for (lev in levs){
        for (g in c(0,1)){
          sub <- data %>% filter(!is.na(.data[[dfs_col]]), .data[[dfs_col]]==g) %>% pull(!!sym(v)) %>% as.character()
          sub <- sub[!is.na(sub) & sub != '']
          cnt <- sum(sub==lev)
          pct <- if (length(sub)==0) 0 else 100*cnt/length(sub)
          out_rows[[length(out_rows)+1]] <- tibble::tibble(variable=v, level=lev, group=as.character(g), value=sprintf('%d (%.1f%%)', cnt, pct))
        }
        # total
        allsub <- as.character(vec)
        allsub <- allsub[!is.na(allsub) & allsub != '']
        cnt_all <- sum(allsub==lev)
        pct_all <- if (length(allsub)==0) 0 else 100*cnt_all/length(allsub)
        out_rows[[length(out_rows)+1]] <- tibble::tibble(variable=v, level=lev, group='Total', value=sprintf('%d (%.1f%%)', cnt_all, pct_all))
      }
    }
  }
  bind_rows(out_rows)
}

# Create Table 1
## Wrap the main analysis in a warning collector so we can save warnings to a file
warnings_collector <- character()
res <- withCallingHandlers({
  table1_df <- make_table(df, dem_vars)
  readr::write_csv(table1_df, file.path(outdir, 'table1_demographics_by_dfs.csv'))
  message('Wrote table1_demographics_by_dfs.csv')

  # Create Table 2: classical + non-classical under separate stubs
  table2_classical <- make_table(df, classical)
  table2_nonclassical <- make_table(df, non_classical)
  readr::write_csv(table2_classical, file.path(outdir, 'table2_classical_by_dfs.csv'))
  readr::write_csv(table2_nonclassical, file.path(outdir, 'table2_nonclassical_by_dfs.csv'))
  message('Wrote table2 CSVs')

  # Association tests and correlations (predictors in table2 vs DFS)
  all_preds <- c(classical, non_classical)
  all_preds <- all_preds[all_preds %in% names(df)]

  assoc_rows <- list()
  univar_rows <- list()
  for (v in all_preds){
  vec <- df[[v]]
  # skip if almost all NA
  if (all(is.na(vec))) next
  # prepare a clean df with predictor and dfs
  dsub <- df %>% select(.dfs, !!sym(v)) %>% filter(!is.na(.dfs))
  # handle categorical vs numeric
  if (is.numeric(vec) || is.integer(vec)){
    # numeric predictor: compute t-test (compare means by DFS), Cohen's d, and Pearson r vs DFS (point-biserial)
    g0 <- dsub %>% filter(.dfs==0) %>% pull(!!sym(v)) %>% as.numeric() %>% na.omit()
    g1 <- dsub %>% filter(.dfs==1) %>% pull(!!sym(v)) %>% as.numeric() %>% na.omit()
    t_stat <- NA_real_
    t_p <- NA_real_
    cohen_d <- NA_real_
    if (length(g0)>1 & length(g1)>1) {
      tt <- tryCatch(t.test(g1, g0), error=function(e) NULL)
      if (!is.null(tt)) {
        t_stat <- as.numeric(tt$statistic)
        t_p <- as.numeric(tt$p.value)
      }
      # Cohen's d (mean1 - mean0)/pooled_sd
      n0 <- length(g0); n1 <- length(g1)
      sd0 <- sd(g0); sd1 <- sd(g1)
      pooled_sd <- tryCatch(sqrt(((n0-1)*sd0^2 + (n1-1)*sd1^2)/(n0+n1-2)), error=function(e) NA_real_)
      if (!is.na(pooled_sd) && pooled_sd>0) cohen_d <- (mean(g1)-mean(g0))/pooled_sd
    }
    # Pearson correlation (point-biserial)
    numpred <- as.numeric(dsub %>% pull(!!sym(v)))
    dfsnum <- as.numeric(dsub$.dfs)
    cor_r <- NA_real_
    cor_p <- NA_real_
    cor_res <- tryCatch(cor.test(numpred, dfsnum, method='pearson', use='complete.obs'), error=function(e) NULL)
    if (!is.null(cor_res)){
      cor_r <- as.numeric(cor_res$estimate)
      cor_p <- as.numeric(cor_res$p.value)
    }
    assoc_rows[[length(assoc_rows)+1]] <- tibble::tibble(variable=v, type='numeric', t_stat=t_stat, t_p=t_p, cohen_d=cohen_d, cor_r=cor_r, cor_p=cor_p)
    # univariate logistic (as before)
    form <- as.formula(paste('.dfs ~', paste0('`', v, '`')))
    m <- tryCatch(glm(form, data=dsub, family=binomial()), error=function(e) NULL)
    if (!is.null(m)){
      s <- broom::tidy(m, conf.int=TRUE, exponentiate=TRUE)
      # take the predictor row (second row)
      if (nrow(s)>=2) {
        r <- s[2,]
        univar_rows[[length(univar_rows)+1]] <- tibble::tibble(variable=v, term=r$term, estimate= r$estimate, conf.low=r$conf.low, conf.high=r$conf.high, p.value=r$p.value)
      }
    }
  } else {
    # categorical: chi-square or fisher
    tab <- table(dsub$.dfs, as.character(dsub[[v]]))
    use_fisher <- any(tab < 5)
    pval <- NA_real_
    if (nrow(tab)>0 & ncol(tab)>0){
      if (use_fisher) pval <- tryCatch(fisher.test(tab)$p.value, error=function(e) NA_real_) else pval <- tryCatch(chisq.test(tab)$p.value, error=function(e) NA_real_)
    }
    assoc_rows[[length(assoc_rows)+1]] <- tibble::tibble(variable=v, type='categorical', p_chisq_or_fisher=pval)
    # univariate logistic: convert to factor and use glm with one-hot as necessary
    # we'll run glm using the variable as factor (R will handle contrasts)
    form <- as.formula(paste('.dfs ~', paste0('`', v, '`')))
    m <- tryCatch(glm(form, data=dsub, family=binomial()), error=function(e) NULL)
    if (!is.null(m)){
      s <- broom::tidy(m, conf.int=TRUE, exponentiate=TRUE)
      # store all non-intercept terms
      if (nrow(s)>=2) {
        s2 <- s[-1, , drop=FALSE]
        s2 <- s2 %>% transmute(variable=v, term=term, estimate=estimate, conf.low=conf.low, conf.high=conf.high, p.value=p.value)
        univar_rows <- c(univar_rows, split(s2, seq(nrow(s2))))
      }
    }
  }
}

assoc_df <- bind_rows(assoc_rows)
univar_df <- bind_rows(univar_rows)

# add q-values (BH) to univariate regression terms (adjust across all terms)
if (nrow(univar_df)>0){
  univar_df <- univar_df %>% mutate(p.value = as.numeric(p.value)) %>%
    mutate(q.value = p.adjust(p.value, method = 'BH'))
  readr::write_csv(univar_df, file.path(outdir, 'univariate_logistic_table2_with_q.csv'))
} else {
  readr::write_csv(tibble::tibble(), file.path(outdir, 'univariate_logistic_table2_with_q.csv'))
}

if (nrow(assoc_df)>0) readr::write_csv(assoc_df, file.path(outdir, 'associations_table2_predictors_vs_dfs.csv'))
if (nrow(univar_df)>0) readr::write_csv(univar_df, file.path(outdir, 'univariate_logistic_table2.csv'))
message('Wrote association and univariate logistic results')

# Multivariable logistic models: classical-only and full
mv_rows <- list()
if (all(classical %in% names(df)) ){
  dsub <- df %>% filter(!is.na(.dfs)) %>% select(.dfs, all_of(classical))
  form <- as.formula(paste('.dfs ~', paste(paste0('`', classical, '`'), collapse = ' + ')))
  m <- tryCatch(glm(form, data=dsub, family=binomial()), error=function(e) NULL)
  if (!is.null(m)){
    tc <- broom::tidy(m, conf.int=TRUE, exponentiate=TRUE)
    if (nrow(tc)>0) tc <- tc %>% mutate(q.value = p.adjust(p.value, method='BH'))
    mv_rows[['classical']] <- tc
  }
}
## Build multivariable 'full' model using predictors with univariate p < 0.5
if (exists('univar_df') && nrow(univar_df)>0){
  # For categorical variables univar_df may contain multiple terms per variable.
  # Include the predictor if any term for that variable has q < 0.5 (BH-adjusted)
  selected_preds <- univar_df %>% dplyr::filter(!is.na(q.value) & q.value < 0.5) %>% dplyr::pull(variable) %>% unique()
  # Ensure predictors exist in df
  preds_here <- selected_preds[selected_preds %in% names(df)]
  readr::write_csv(tibble::tibble(selected_predictors = preds_here), file.path(outdir, 'multivariable_selected_predictors.csv'))
  message('Selected predictors for multivariable (q<0.5 univar BH): ', paste(head(preds_here, 50), collapse = ', '))
  if (length(preds_here)>0){
    dsub <- df %>% filter(!is.na(.dfs)) %>% select(.dfs, all_of(preds_here))
    form <- as.formula(paste('.dfs ~', paste(paste0('`', preds_here, '`'), collapse = ' + ')))
    mfull <- tryCatch(glm(form, data=dsub, family=binomial()), error=function(e) e)
    if (!inherits(mfull, 'error')){
      tf <- broom::tidy(mfull, conf.int=TRUE, exponentiate=TRUE)
      if (nrow(tf)>0) tf <- tf %>% mutate(q.value = p.adjust(p.value, method='BH'))
      mv_rows[['full']] <- tf
    } else {
      mv_rows[['full_error']] <- tibble::tibble(message = conditionMessage(mfull))
    }
  }
} else {
  # Fallback: use all_preds if univar_df not available
  preds_here <- all_preds[all_preds %in% names(df)]
  if (length(preds_here)>0){
    dsub <- df %>% filter(!is.na(.dfs)) %>% select(.dfs, all_of(preds_here))
    form <- as.formula(paste('.dfs ~', paste(paste0('`', preds_here, '`'), collapse = ' + ')))
    mfull <- tryCatch(glm(form, data=dsub, family=binomial()), error=function(e) NULL)
    if (!is.null(mfull)) mv_rows[['full']] <- broom::tidy(mfull, conf.int=TRUE, exponentiate=TRUE)
  }
}
if (length(mv_rows)>0){
  for (n in names(mv_rows)) readr::write_csv(mv_rows[[n]], file.path(outdir, paste0('multivariable_logistic_', n, '.csv')))
  message('Wrote multivariable logistic results')
}

}, warning = function(w){
  # collect warning message and then continue
  warnings_collector <<- c(warnings_collector, conditionMessage(w))
  invokeRestart('muffleWarning')
})

# write warnings to file
if (length(warnings_collector)>0){
  writeLines(unique(warnings_collector), con = file.path(outdir, 'generate_tables_1_2_warnings.txt'))
  message('Wrote warnings to: ', file.path(outdir, 'generate_tables_1_2_warnings.txt'))
}

message('All done. Outputs in: ', outdir)
