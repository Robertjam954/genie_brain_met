#!/usr/bin/env python3
"""
Survival and ML analysis pipeline

Features:
- Univariate and multivariate Cox PH
- Kaplan–Meier curves (overall and stratified by SUBTYPE and MANTIS_BIN)
- Fine–Gray competing risks model (if competing event available)
- XGBoost AFT survival model

Inputs:
- datasets_analysis_dictionary/merged_genie.xlsx (expects time + event columns)
- output/fine_gray_python/per_gene_cs_cox_sksurv_bootstrap.csv (for top genes; optional)

Outputs:
- output/survival/... (CSVs, PNGs, JSON model)

Notes:
- Column names are normalized/auto-detected when possible.
- If required columns are missing, the step is skipped with a message.
"""
import os
import re
import json
import math
import argparse
import warnings
from io import StringIO
from contextlib import redirect_stdout

import numpy as np
import pandas as pd

# local normalize helper if present
import sys
sys.path.insert(0, os.path.dirname(__file__))
try:
    from col_normalize import normalize_columns
except Exception:
    def normalize_columns(df):
        return df

# plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# survival
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
try:
    from lifelines import FineAndGrayFitter
    HAS_FINEGRAY = True
except Exception:
    HAS_FINEGRAY = False

# ML
try:
    import xgboost as xgb
    HAS_XGB = True
except Exception:
    HAS_XGB = False

from lifelines.utils import concordance_index


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def detect_col(df: pd.DataFrame, candidates, regex=False, case_insensitive=True):
    cols = list(df.columns)
    if case_insensitive:
        lut = {c.lower(): c for c in cols}
        if not regex:
            for cand in candidates:
                key = cand.lower()
                if key in lut:
                    return lut[key]
        else:
            for c in cols:
                for pat in candidates:
                    flags = re.IGNORECASE if case_insensitive else 0
                    if re.search(pat, c, flags):
                        return c
    else:
        for cand in candidates:
            if cand in cols:
                return cand
    return None


def to_numeric_safe(s):
    return pd.to_numeric(s, errors='coerce')


def load_top_genes(csv_path: str, topn=10):
    if not os.path.exists(csv_path):
        return []
    try:
        pg = pd.read_csv(csv_path)
        if 'n' in pg.columns:
            pg = pg.sort_values('n', ascending=False)
        if 'gene' not in pg.columns:
            return []
        return pg['gene'].head(topn).dropna().astype(str).tolist()
    except Exception:
        return []


def build_gene_indicators(df: pd.DataFrame, top_genes: list):
    # Try to parse HUGO gene list if present
    hugo_col = None
    for c in df.columns:
        if isinstance(c, str) and c.strip().upper().startswith('HUGO'):
            hugo_col = c
            break
    if hugo_col is not None:
        parsed = [set([g.strip() for g in re.split(r'[;|,]', str(v)) if g.strip()]) for v in df[hugo_col].fillna('')]
        for g in top_genes:
            gcol = 'G__' + re.sub(r"\s+","_", g)
            if gcol not in df.columns:
                df[gcol] = [1 if g in s else 0 for s in parsed]
    return df


def prepare_common_covariates(df: pd.DataFrame):
    # Basic numeric covariates with defensive casting
    for col in ['BUFFA_HYPOXIA_SCORE','MUT_COUNT','TBL_LOW','TBL_HIGH','AGE','AJCC_STAGE_NUM','ETHNICITY_BIN']:
        if col not in df.columns:
            df[col] = 0
        df[col] = to_numeric_safe(df[col]).fillna(0)

    # If MUT_COUNT missing, derive from gene indicators
    if 'MUT_COUNT' not in df.columns:
        gene_cols = [c for c in df.columns if isinstance(c,str) and c.startswith('G__')]
        df['MUT_COUNT'] = df[gene_cols].sum(axis=1) if gene_cols else 0

    # dummies for MANTIS_BIN and SUBTYPE
    if 'MANTIS_BIN' in df.columns and df['MANTIS_BIN'].dtype.name == 'category':
        df['MANTIS_BIN'] = df['MANTIS_BIN'].astype(str)
    if 'SUBTYPE' in df.columns and df['SUBTYPE'].dtype.name == 'category':
        df['SUBTYPE'] = df['SUBTYPE'].astype(str)

    mantis = pd.get_dummies(df.get('MANTIS_BIN', 'Unknown'), prefix='MANTIS', drop_first=True)
    subtype = pd.get_dummies(df.get('SUBTYPE', 'Unknown'), prefix='SUBTYPE', drop_first=True)

    cov_base = df[['BUFFA_HYPOXIA_SCORE','MUT_COUNT','TBL_LOW','TBL_HIGH','AGE','AJCC_STAGE_NUM','ETHNICITY_BIN']].copy()
    X_common = pd.concat([cov_base.reset_index(drop=True), mantis.reset_index(drop=True), subtype.reset_index(drop=True)], axis=1)
    return X_common


def save_km_plot(df, time_col, event_col, group_col, out_png, title):
    kmf = KaplanMeierFitter()
    plt.figure(figsize=(8,6))
    for level, gdf in df.groupby(group_col):
        if gdf[time_col].notna().sum() == 0:
            continue
        kmf.fit(durations=gdf[time_col], event_observed=gdf[event_col], label=f"{group_col}={level}")
        kmf.plot(ci_show=True)
    plt.title(title)
    plt.xlabel('Time')
    plt.ylabel('Survival probability')
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def run_univariate_cox(df, time_col, event_col, covariate_cols, out_csv):
    rows = []
    for col in covariate_cols:
        try:
            sub = df[[time_col, event_col, col]].copy()
            # ensure numeric
            if sub[col].dtype.name not in ('int64','float64'):
                # if binary-like strings, coerce
                sub[col] = to_numeric_safe(sub[col])
            # drop rows with NA in covariate
            sub = sub.dropna(subset=[time_col, event_col, col])
            if sub[col].nunique() < 2:
                continue
            cph = CoxPHFitter()
            cph.fit(sub, duration_col=time_col, event_col=event_col)
            summ = cph.summary.loc[col]
            rows.append({
                'variable': col,
                'coef': float(summ['coef']),
                'HR': float(summ['exp(coef)']),
                'CI_lower': float(summ['exp(coef) lower 95%']),
                'CI_upper': float(summ['exp(coef) upper 95%']),
                'p': float(summ['p'])
            })
        except Exception as e:
            # skip problematic covariates silently
            continue
    out = pd.DataFrame(rows).sort_values('p')
    out.to_csv(out_csv, index=False)
    return out


def run_multivariate_cox(df, time_col, event_col, covariate_cols, out_csv, out_diag_txt):
    sub = df[[time_col, event_col] + covariate_cols].copy()
    # numeric conversion where possible
    for c in covariate_cols:
        if sub[c].dtype.name not in ('int64','float64'):
            sub[c] = to_numeric_safe(sub[c])
    # drop rows with any NA
    sub = sub.dropna()
    # filter degenerate columns
    keep = [c for c in covariate_cols if sub[c].nunique() >= 2]
    sub = sub[[time_col, event_col] + keep]
    cph = CoxPHFitter()
    cph.fit(sub, duration_col=time_col, event_col=event_col)
    cph.summary.to_csv(out_csv)

    # PH diagnostics
    with open(out_diag_txt, 'w') as f:
        try:
            # global PH test
            ph_test = proportional_hazard_test(cph, sub, time_transform='rank')
            f.write('Global Schoenfeld test and per-variable tests:\n')
            f.write(str(ph_test.summary))
            f.write('\n\nModel summary:\n')
            f.write(str(cph.summary))
        except Exception as e:
            f.write('PH diagnostics failed: ' + repr(e))
    return cph


def run_fine_gray(df, time_col, event_of_interest_col, competing_event_col, covariate_cols, out_csv):
    if not HAS_FINEGRAY:
        return None
    # Fine & Gray expects cause (1) vs competing (1) mutually exclusive; others censored
    # lifelines FineAndGrayFitter: event_col is event of interest indicator, competing_event_col is competing
    sub = df[[time_col, event_of_interest_col, competing_event_col] + covariate_cols].copy().dropna()
    # ensure binary integers
    for c in [event_of_interest_col, competing_event_col]:
        sub[c] = (to_numeric_safe(sub[c]) > 0).astype(int)
    # filter degenerate
    keep = [c for c in covariate_cols if sub[c].nunique() >= 2]
    sub = sub[[time_col, event_of_interest_col, competing_event_col] + keep]

    fg = FineAndGrayFitter()
    fg.fit(sub, duration_col=time_col, event_col=event_of_interest_col, competing_event_col=competing_event_col)
    # collect summary if available
    try:
        summ = fg.summary
        summ.to_csv(out_csv)
    except Exception:
        pd.DataFrame({'var': keep}).to_csv(out_csv, index=False)
    return fg


def prepare_aft_labels(durations, events):
    # exact events: lower=upper=duration; censored: lower=duration, upper=+inf
    lower = np.asarray(durations, dtype=float)
    upper = np.asarray(durations, dtype=float)
    events = np.asarray(events).astype(bool)
    upper[~events] = np.inf
    return lower, upper


def run_xgb_aft(df, time_col, event_col, feature_cols, out_dir):
    if not HAS_XGB:
        return None
    sub = df[[time_col, event_col] + feature_cols].copy().dropna()
    # cast to float
    for c in feature_cols:
        sub[c] = to_numeric_safe(sub[c]).fillna(0)
    durations = sub[time_col].values.astype(float)
    events = (to_numeric_safe(sub[event_col]).fillna(0).values > 0).astype(int)
    X = sub[feature_cols].values.astype(float)

    # train/val split
    rng = np.random.RandomState(42)
    idx = np.arange(len(sub))
    rng.shuffle(idx)
    split = int(0.8 * len(sub))
    tr, va = idx[:split], idx[split:]

    dtrain = xgb.DMatrix(X[tr, :])
    dvalid = xgb.DMatrix(X[va, :])
    l_lower, l_upper = prepare_aft_labels(durations, events)
    dtrain.set_float_info('label_lower_bound', l_lower[tr])
    dtrain.set_float_info('label_upper_bound', l_upper[tr])
    dvalid.set_float_info('label_lower_bound', l_lower[va])
    dvalid.set_float_info('label_upper_bound', l_upper[va])

    params = {
        'objective': 'survival:aft',
        'eval_metric': 'aft-nloglik',
        'tree_method': 'hist',
        'learning_rate': 0.05,
        'max_depth': 4,
        'min_child_weight': 1.0,
        'lambda': 1.0,
        'alpha': 0.0,
        'aft_loss_distribution': 'normal',
        'aft_loss_distribution_scale': 1.5,
        'verbosity': 0
    }

    evals = [(dtrain, 'train'), (dvalid, 'valid')]
    bst = xgb.train(params, dtrain, num_boost_round=1000, evals=evals,
                    early_stopping_rounds=50, verbose_eval=False)

    # predict median survival time (quantile 0.5)
    dall = xgb.DMatrix(X)
    pred = bst.predict(dall, output_margin=False, pred_contribs=False, training=False, 
                       iteration_range=(0, bst.best_iteration+1),
                       base_margin=None, 
                       pred_interactions=False,
                       validate_features=True,
                       missing=np.nan,
                       strict_shape=False,
                       iteration_end=None,
                       approx_contribs=False,
                       pred_leaf=False,
                       pred_model_shap_contrib=False,
                       pred_margin=False,
                       ntree_limit=0,
                       pred_side_outputs=False,
                       pred_aft_quantile=0.5)

    # Evaluate concordance on validation
    cidx = concordance_index(durations[va], -pred[va], events[va])  # negative for risk ordering

    # Save model
    ensure_dir(out_dir)
    model_path = os.path.join(out_dir, 'xgb_aft_model.json')
    bst.save_model(model_path)
    with open(os.path.join(out_dir, 'xgb_aft_metrics.json'), 'w') as f:
        json.dump({'concordance_index_valid': float(cidx), 'best_iteration': int(bst.best_iteration)}, f, indent=2)

    # Feature importance
    imp = bst.get_score(importance_type='gain')
    if imp:
        fi = pd.DataFrame({'feature': list(imp.keys()), 'gain': list(imp.values())}).sort_values('gain', ascending=False)
        fi.to_csv(os.path.join(out_dir, 'xgb_feature_importance.csv'), index=False)
    return bst


def main():
    p = argparse.ArgumentParser(description='Survival and ML analysis pipeline')
    p.add_argument('--xlsx', default=os.path.join('datasets_analysis_dictionary','merged_genie.xlsx'))
    p.add_argument('--per-gene-csv', default=os.path.join('output','fine_gray_python','per_gene_cs_cox_sksurv_bootstrap.csv'))
    p.add_argument('--outdir', default=os.path.join('output','survival'))
    p.add_argument('--time-col', default=None, help='Override time-to-event column name')
    p.add_argument('--event-col', default=None, help='Override event indicator column name (1=event,0=censored)')
    p.add_argument('--competing-col', default=None, help='Competing event column for Fine–Gray')
    p.add_argument('--strata', default='SUBTYPE,MANTIS_BIN', help='Comma-separated columns to stratify KM and Cox by')
    p.add_argument('--topn-genes', type=int, default=10)
    args = p.parse_args()

    ensure_dir(args.outdir)
    subdirs = {k: os.path.join(args.outdir, k) for k in ['km','cox','finegray','xgb','logs']}
    for d in subdirs.values():
        ensure_dir(d)

    if not os.path.exists(args.xlsx):
        raise FileNotFoundError(args.xlsx)

    df = pd.read_excel(args.xlsx, engine='openpyxl')
    df = normalize_columns(df)

    # columns detection
    time_col = args.time_col or detect_col(df, ['DFS_MONTHS','DFS_TIME','TIME_TO_DFS','TIME_TO_EVENT','MONTHS_TO_DFS','PFS_MONTHS','OS_MONTHS'], regex=False)
    event_col = args.event_col or detect_col(df, ['DFS_STATUS','DFS_EVENT','EVENT','PFS_EVENT','OS_EVENT','OS_STATUS'], regex=False)

    if time_col is None or event_col is None:
        raise RuntimeError(f"Could not detect time/event columns. time_col={time_col}, event_col={event_col}. Use --time-col/--event-col.")

    # clean time/event
    df[time_col] = to_numeric_safe(df[time_col])
    df[event_col] = (to_numeric_safe(df[event_col]) > 0).astype(int)

    # top genes
    top_genes = load_top_genes(args.per_gene_csv, topn=args.topn_genes)
    df = build_gene_indicators(df, top_genes)

    # prepare covariates
    X_common = prepare_common_covariates(df)

    # assemble covariates list for Cox
    gene_cols = [c for c in df.columns if isinstance(c,str) and c.startswith('G__')]
    base_cov = list(X_common.columns)

    # KM overall by top-5 gene status if available, otherwise by any mutation
    if 'top_5_gene_status' in df.columns:
        km_group = 'top_5_gene_status'
    elif gene_cols:
        df['any_top_gene'] = (df[gene_cols].sum(axis=1) > 0).astype(int)
        km_group = 'any_top_gene'
    else:
        km_group = None

    if km_group is not None:
        try:
            save_km_plot(df[[time_col,event_col,km_group]].dropna(), time_col, event_col, km_group,
                         os.path.join(subdirs['km'], f'KM_by_{km_group}.png'), title=f'KM by {km_group}')
        except Exception:
            pass

    # Stratified KM by strata
    strata_cols = [s.strip() for s in args.strata.split(',') if s.strip()]
    for s in strata_cols:
        if s not in df.columns or km_group is None:
            continue
        for level, sdf in df.dropna(subset=[s]).groupby(s):
            try:
                outpng = os.path.join(subdirs['km'], f'KM_by_{km_group}_strat_{s}_{str(level)}.png')
                save_km_plot(sdf[[time_col,event_col,km_group]].dropna(), time_col, event_col, km_group,
                             outpng, title=f'KM by {km_group} stratified by {s} = {level}')
            except Exception:
                continue

    # Univariate Cox: run over gene indicators and base covariates
    uni_covs = gene_cols + base_cov
    uni = run_univariate_cox(df, time_col, event_col, uni_covs, os.path.join(subdirs['cox'], 'cox_univariate.csv'))

    # Multivariate Cox: include genes with p<0.1 + base covariates (cap number)
    selected_genes = uni.query('p < 0.1 and variable.str.startswith("G__")', engine='python')['variable'].tolist()
    selected_genes = selected_genes[: min(10, len(selected_genes))]
    multi_covs = list(dict.fromkeys(base_cov + selected_genes))  # preserve order unique
    if multi_covs:
        cph = run_multivariate_cox(df, time_col, event_col, multi_covs,
                                   os.path.join(subdirs['cox'], 'cox_multivariate_summary.csv'),
                                   os.path.join(subdirs['logs'], 'cox_ph_diagnostics.txt'))

    # Fine–Gray if competing available
    comp_col = args.competing_col or detect_col(df, ['DEATH_EVENT','OS_STATUS','DECEASED','COMPETING_EVENT'], regex=False)
    if HAS_FINEGRAY and comp_col is not None and comp_col != event_col:
        try:
            df[comp_col] = (to_numeric_safe(df[comp_col]) > 0).astype(int)
            # ensure mutually exclusive: if event==1, set competing=0
            df.loc[df[event_col] == 1, comp_col] = 0
            fg_covs = multi_covs if multi_covs else base_cov
            run_fine_gray(df, time_col, event_col, comp_col, fg_covs,
                          os.path.join(subdirs['finegray'], 'finegray_summary.csv'))
        except Exception:
            warnings.warn('Fine–Gray step failed; skipping.')

    # Stratified Cox by subtype and MANTIS_BIN (if present)
    for s in strata_cols:
        if s not in df.columns:
            continue
        for level, sdf in df.dropna(subset=[s]).groupby(s):
            try:
                sub_out = os.path.join(subdirs['cox'], f'cox_multivariate_{s}_{str(level)}.csv')
                # reuse same covariates
                if multi_covs:
                    run_multivariate_cox(sdf, time_col, event_col, multi_covs, sub_out,
                                         os.path.join(subdirs['logs'], f'cox_ph_diag_{s}_{str(level)}.txt'))
            except Exception:
                continue

    # XGBoost AFT model
    if HAS_XGB:
        # Build feature cols: base + genes
        feat_cols = list(dict.fromkeys(base_cov + gene_cols))
        xgb_dir = subdirs['xgb']
        try:
            run_xgb_aft(df, time_col, event_col, feat_cols, xgb_dir)
        except Exception:
            warnings.warn('XGB AFT step failed; skipping.')

    # Write a small run log
    with open(os.path.join(args.outdir, 'run_summary.json'), 'w') as f:
        json.dump({
            'time_col': time_col,
            'event_col': event_col,
            'km_group': km_group,
            'n_rows': int(df.shape[0]),
            'has_finegray': bool(HAS_FINEGRAY),
            'has_xgb': bool(HAS_XGB),
            'strata_used': [s for s in strata_cols if s in df.columns]
        }, f, indent=2)

    print('Analysis complete. Outputs written to', args.outdir)


if __name__ == '__main__':
    main()
