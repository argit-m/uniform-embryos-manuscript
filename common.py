"""Shared code for the manuscript figure notebooks.

Run order:  fig2 -> fig3 -> fig5 -> fig6
  fig3 writes the MF time-course cluster assignments and fig5 writes the larval
  tissue-enriched protein list; fig6 reads both (from outputs/intermediate/).

Sections
  1. Paths and saving
  2. TMT loading + normalisation
  3. FlyBase gene groups
  4. MF lines at a single stage (Fig 2, Fig 3E, Fig S2, Tables S1/S4)
  5. Larval tissue colours (Fig 5, Fig 6)
"""
import gzip
import re
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
from scipy.stats import pearsonr, ttest_ind
from statsmodels.stats.multitest import multipletests

mpl.rcParams['pdf.fonttype'] = 42          # editable text in Illustrator
mpl.rcParams['ps.fonttype'] = 42


# ══════════════════════════════════════════════════════════════════════════════
# 1. Paths and saving
# ══════════════════════════════════════════════════════════════════════════════
FIG_DIR          = Path(__file__).resolve().parent      # manuscript_figures/
DATA_DIR         = FIG_DIR / 'raw data'                  # TMT protein_quant TSVs
HELPER_DIR       = FIG_DIR / 'helper files'              # annotation / reference files
OUT_DIR          = FIG_DIR / 'outputs'
TABLE_DIR        = OUT_DIR / 'supplementary_tables'
INTERMEDIATE_DIR = OUT_DIR / 'intermediate'              # files passed between notebooks

GENETABLE       = HELPER_DIR / 'genetable.xlsx'
GENE_GROUP_FILE = HELPER_DIR / 'gene_group_data_fb_2025_05.tsv.gz'

for _d in (OUT_DIR, TABLE_DIR, INTERMEDIATE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


def fig_dir(name):
    """outputs/<name>/, created on first use."""
    d = OUT_DIR / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def save_fig(fig, out_dir, stem, dpi=300):
    """Save fig as <out_dir>/<stem>.pdf and .png."""
    for ext in ('pdf', 'png'):
        fig.savefig(Path(out_dir) / f'{stem}.{ext}', dpi=dpi, bbox_inches='tight', facecolor='white')
    print(f'Saved {Path(out_dir).name}/{stem}.pdf / .png')


def save_table(df, name):
    """Write a supplementary table to outputs/supplementary_tables/<name>.csv."""
    path = TABLE_DIR / f'{name}.csv'
    df.to_csv(path, index=False)
    print(f'Saved {path.name}  ({len(df):,} rows x {df.shape[1]} columns)')


def require(*paths):
    """Raise a clear error if an input file is missing."""
    missing = [str(p) for p in paths if not Path(p).exists()]
    if missing:
        raise FileNotFoundError('Missing input file(s):\n  ' + '\n  '.join(missing))


# ══════════════════════════════════════════════════════════════════════════════
# 2. TMT loading + normalisation
#    Per protein over the selected channels:
#      row / row mean  ->  column / column median  ->  row / row mean
# ══════════════════════════════════════════════════════════════════════════════
def normalize_tmt(filename, channel_names, variable_names, min_sumSN, genetable_path=GENETABLE):
    """Select channels by NAME ('126_sn', ...); output columns are <variable_name>_norm.
    Used for the MF time course (Fig 3), larval tissues (Fig 5) and CanS series (Fig 6)."""
    if len(channel_names) != len(variable_names):
        raise ValueError('channel_names and variable_names must have the same length.')

    raw = pd.read_csv(filename, sep='\t')
    raw.columns = raw.columns.str.replace(r'^\w+~rq_(.*) sum$', r'\1', regex=True)
    missing = [c for c in channel_names if c not in raw.columns]
    if missing:
        raise KeyError(f'Columns not found in {filename}: {missing}')

    raw['sumSN'] = raw[channel_names].sum(axis=1)
    kept = raw[raw['sumSN'] > min_sumSN]
    kept = kept[~kept['Protein Id'].str.startswith('##', na=False)
                & ~kept['Protein Id'].str.endswith('contaminant', na=False)].reset_index(drop=True)

    vals   = kept[channel_names].to_numpy(float)
    norm   = vals / vals.mean(axis=1, keepdims=True)
    norm   = norm / np.nanmedian(norm, axis=0)
    scaled = norm / norm.mean(axis=1, keepdims=True)
    normed = pd.DataFrame(scaled, columns=[f'{v}_norm' for v in variable_names], index=kept.index)

    data = pd.concat([kept, normed], axis=1).rename(columns={'Protein Id': 'ProteinId'})
    out = data.merge(pd.read_excel(genetable_path), on='ProteinId', how='inner')
    print(f'{Path(filename).name}: {len(raw)} rows -> {len(kept)} pass S/N + decoy filter '
          f'-> {len(out)} with gene annotation')
    return out


def normalize_tmt_by_index(filename, selected_indices, variable_names, min_sumSN,
                           genetable_path=GENETABLE):
    """Select channels by 0-based TSV column INDEX; output columns are named exactly
    variable_names (e.g. NE_a). Also adds 'parent_value' (FBgn parsed from Description).
    Used for the MF lines at a single stage (Fig 2, Fig 3E)."""
    rep = pd.read_csv(filename, sep='\t')
    rep['parent_value'] = rep['Description'].str.extract(r'parent=(FBgn\d+)')

    rep['sumSN'] = rep.iloc[:, selected_indices].sum(axis=1)
    filtered = rep[rep['sumSN'] > min_sumSN].copy()
    pid = filtered['Protein Id'].astype(str)
    filtered = filtered[~pid.str.startswith('##') & ~pid.str.endswith('contaminant')].reset_index(drop=True)

    raw = filtered.iloc[:, selected_indices].to_numpy(dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        normalized        = raw / raw.mean(axis=1, keepdims=True)
        median_normalized = normalized / np.nanmedian(normalized, axis=0)
        scaled            = median_normalized / median_normalized.mean(axis=1, keepdims=True)

    processed = pd.concat([filtered, pd.DataFrame(scaled, columns=variable_names)],
                          axis=1).rename(columns={'Protein Id': 'ProteinId'})
    return processed.merge(pd.read_excel(genetable_path), on='ProteinId', how='inner')


def detect_fbgn_col(df, fallback='gene_ID'):
    """First column whose values look like FlyBase gene IDs (FBgn...)."""
    pat = re.compile(r'^FBgn\d+$')
    for col in df.columns:
        sample = df[col].dropna().astype(str).head(20)
        if len(sample) and sample.str.match(pat).mean() > 0.5:
            return col
    return fallback


# ══════════════════════════════════════════════════════════════════════════════
# 3. FlyBase gene groups
# ══════════════════════════════════════════════════════════════════════════════
def read_gene_groups(path=GENE_GROUP_FILE):
    """FlyBase gene_group_data .tsv.gz, skipping the ## metadata header."""
    skip = 0
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('##') or not line.strip() or (line.startswith('#') and '\t' not in line):
                skip += 1
            else:
                break
    df = pd.read_csv(path, sep='\t', compression='gzip', skiprows=skip, dtype=str, low_memory=False)
    if len(df.columns) and df.columns[0].startswith('#'):
        df = df.rename(columns={df.columns[0]: df.columns[0].lstrip('#')})
    df.columns = df.columns.str.strip()
    return df


# ══════════════════════════════════════════════════════════════════════════════
# 4. MF lines at a single stage (Fig 2, Fig 3E, Fig S2, Tables S1 and S4)
# ══════════════════════════════════════════════════════════════════════════════
MF_FILES = {
    'b0': DATA_DIR / 'protein_quant_2228_cns15_b0.tsv',     # blastoderm, 3 hpf
    'b4': DATA_DIR / 'protein_quant_2229_cns16_b4.tsv',     # blastoderm + 4 h, 7 hpf
}
MF_STAGE_NAMES = {'b0': 'blastoderm (3 hpf)', 'b4': 'blastoderm +4 (7 hpf)'}
MF_STAGE_HPF   = {'b0': '3 hpf', 'b4': '7 hpf'}

MF_MIN_SUMSN_PAIR = 160     # summed S/N over the 8 channels of a pairwise subset
MF_MIN_SUMSN_ALL  = 320     # summed S/N over all 16 channels
Q_THRESHOLD   = 0.05        # BH-FDR cutoff
LOG2FC_CUTOFF = 1.0         # 2-fold

# 4 lines x 4 replicates; 'start' = TSV column index of replicate a
MF_CONDITIONS = [
    dict(short='NE', full='neuroectoderm',      start=7,  dark='#1b5e8c', light='#a9cce3'),
    dict(short='LP', full='lateral proctodeum', start=11, dark='#2e7d32', light='#a5d6a7'),
    dict(short='DP', full='dorsal proctodeum',  start=15, dark='#c65911', light='#f5c39b'),
    dict(short='DE', full='dorsal ectoderm',    start=19, dark='#6a1b9a', light='#cba6dd'),
]
MF_REPS     = 'abcd'
MF_ORDER    = [c['short'] for c in MF_CONDITIONS]
MF_COND     = {c['short']: c for c in MF_CONDITIONS}
MF_COLS     = {c['short']: [f"{c['short']}_{r}" for r in MF_REPS] for c in MF_CONDITIONS}
MF_IDX      = {c['short']: list(range(c['start'], c['start'] + len(MF_REPS))) for c in MF_CONDITIONS}
# 6 comparisons, +log2FC = higher in g1
MF_COMPARISONS = [dict(key=f"{a['short']}_{b['short']}", g1=a['short'], g2=b['short'])
                  for a, b in combinations(MF_CONDITIONS, 2)]


def build_mf_datasets(stage):
    """{'all': 16-channel dataset, 'NE_LP': 8-channel pairwise dataset, ...} for one stage."""
    fname = MF_FILES[stage]
    require(fname, GENETABLE)
    datasets = {'all': normalize_tmt_by_index(
        fname, [i for s in MF_ORDER for i in MF_IDX[s]],
        [n for s in MF_ORDER for n in MF_COLS[s]], MF_MIN_SUMSN_ALL)}
    for comp in MF_COMPARISONS:
        g1, g2 = comp['g1'], comp['g2']
        datasets[comp['key']] = normalize_tmt_by_index(
            fname, MF_IDX[g1] + MF_IDX[g2], MF_COLS[g1] + MF_COLS[g2], MF_MIN_SUMSN_PAIR)
    for key, df in datasets.items():
        print(f'  {stage} {key}: {len(df):,} proteins')
    return datasets


def run_ttest(data, group1_cols, group2_cols):
    """Welch t-test per protein + BH-FDR; log2FC = log2(median g1 / median g2).
    A protein needs >= 2 finite values per group to be tested."""
    g1m = data[group1_cols].to_numpy(dtype=float)
    g2m = data[group2_cols].to_numpy(dtype=float)
    n = len(data)

    p_values = np.full(n, np.nan)
    for i in range(n):
        a = g1m[i][np.isfinite(g1m[i])]
        b = g2m[i][np.isfinite(g2m[i])]
        if len(a) >= 2 and len(b) >= 2:
            p_values[i] = ttest_ind(a, b, equal_var=False).pvalue

    with np.errstate(divide='ignore', invalid='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        log2fc = np.log2(np.nanmedian(g1m, axis=1) / np.nanmedian(g2m, axis=1))

    q_values = np.full(n, np.nan)
    valid = np.isfinite(p_values)
    if valid.any():
        q_values[valid] = multipletests(p_values[valid], method='fdr_bh')[1]
    sig = np.isfinite(q_values) & (q_values < Q_THRESHOLD) & (np.abs(log2fc) >= LOG2FC_CUTOFF)

    out = data.copy()
    out['p_value'], out['q_value'], out['log2FC'], out['significant'] = p_values, q_values, log2fc, sig
    return out


def run_mf_comparisons(datasets, dataset_key=None):
    """Run the 6 comparisons for one stage.
      dataset_key=None  -> each comparison on its own pairwise-normalised subset (used for the % matrix)
      dataset_key='all' -> every comparison on the 16-channel dataset (used for Tables S1/S4)
    Returns (results {key: df}, counts {key: (n_sig, n_total)})."""
    results, counts = {}, {}
    for comp in MF_COMPARISONS:
        df = run_ttest(datasets[dataset_key or comp['key']], MF_COLS[comp['g1']], MF_COLS[comp['g2']])
        results[comp['key']] = df
        counts[comp['key']] = (int(df['significant'].sum()), len(df))
    return results, counts


def mf_stats_table(datasets, results_full):
    """Supplementary Table S1 / S4: 16 normalised channels + log2FC/q/sig per comparison."""
    table = datasets['all'].copy()
    for key, df in results_full.items():
        table[f'log2FC_{key}'] = df['log2FC'].to_numpy()
        table[f'q_{key}']      = df['q_value'].to_numpy()
        table[f'sig_{key}']    = df['significant'].to_numpy()
    drop = ['Gene Symbol', 'Description', 'Group ID', 'Collapsed', 'Parsimony', 'parent_value']
    drop += [c for c in table.columns if c.startswith('default~rq_')]       # raw S/N channels
    return table.drop(columns=drop)


def plot_pct_matrix(counts, stage):
    """Lower-triangular grid: % of proteins differentially abundant per pair of lines (Fig 2B, 3E)."""
    conds = MF_CONDITIONS
    n = len(conds)

    def pct(code_a, code_b):
        key = f'{code_a}_{code_b}' if any(c['key'] == f'{code_a}_{code_b}' for c in MF_COMPARISONS) \
              else f'{code_b}_{code_a}'
        n_sig, total = counts[key]
        return 100 * n_sig / total if total else np.nan

    fig, ax = plt.subplots(figsize=(4.6, 4.6))
    box, gap = 1.0, 0.06
    for i in range(1, n):                 # rows: LP, DP, DE
        for j in range(i):                # columns: lines before the row line
            x, y = j, -(i - 1)
            ax.add_patch(mpatches.FancyBboxPatch(
                (x + gap / 2, y + gap / 2), box - gap, box - gap,
                boxstyle='square,pad=0', facecolor='#f0efed', edgecolor='#dcdcd8', linewidth=1.2))
            ax.text(x + 0.5, y + 0.5, f"{pct(conds[j]['short'], conds[i]['short']):.1f}%",
                    ha='center', va='center', fontsize=20, fontweight='bold', color='#2b2b2b')
    for i in range(1, n):
        ax.text(-0.12, -(i - 1) + 0.5, conds[i]['short'], ha='right', va='center',
                fontsize=17, fontweight='bold', color=conds[i]['dark'])
    for j in range(n - 1):
        ax.text(j + 0.5, -(n - 2) - 0.2, conds[j]['short'], ha='center', va='top',
                fontsize=17, fontweight='bold', color=conds[j]['dark'])

    ax.set_xlim(-0.7, n - 1 + 0.05)
    ax.set_ylim(-(n - 2) - 0.65, 1.05)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(MF_STAGE_HPF[stage], fontsize=15, pad=10)
    fig.suptitle('Proteome differentially abundant\nbetween uniform-fate lines', fontsize=16, y=1.02)
    fig.text(0.5, -0.03, f'% of proteins with FDR q < {Q_THRESHOLD} & ≥ {2 ** LOG2FC_CUTOFF:g}-fold change',
             ha='center', fontsize=10, color='#555555')
    fig.tight_layout()
    return fig


def plot_best_replicate_density(all_by_stage, out_dir, gridsize=60):
    """Fig S2: per line, the replicate pair with the highest Pearson r, as a hexbin.
    all_by_stage = {stage: 16-channel dataset}. Axis range is shared across all stages given."""
    def xy(df, a, b):
        x, y = df[a].to_numpy(float), df[b].to_numpy(float)
        m = np.isfinite(x) & np.isfinite(y)
        return x[m], y[m]

    panels, gmin, gmax = {}, np.inf, -np.inf
    for stage, df in all_by_stage.items():
        print(f'=== {MF_STAGE_NAMES[stage]} ===')
        for short in MF_ORDER:
            pair_r = [((a, b), pearsonr(*xy(df, a, b))[0]) for a, b in combinations(MF_COLS[short], 2)]
            (a, b), best_r = pair_r[int(np.argmax([r for _, r in pair_r]))]
            xl, yl = xy(df, a, b)
            panels[(stage, short)] = (a, b, xl, yl, best_r)
            gmin, gmax = min(gmin, xl.min(), yl.min()), max(gmax, xl.max(), yl.max())
            print(f'  {short}: ' + ', '.join(f'{p[0][-1]}~{p[1][-1]} {r:.3f}' for p, r in pair_r)
                  + f'  -> shown: {a} vs {b} (r = {best_r:.3f})')
    pad = 0.03 * (gmax - gmin)
    lims = [gmin - pad, gmax + pad]

    for stage in all_by_stage:
        fig, axes = plt.subplots(1, len(MF_ORDER), figsize=(4.2 * len(MF_ORDER), 4.4), constrained_layout=True)
        hexes, counts_max = [], 0
        for ax, short in zip(axes, MF_ORDER):
            a, b, xl, yl, best_r = panels[(stage, short)]
            accent = MF_COND[short]['dark']
            hb = ax.hexbin(xl, yl, gridsize=gridsize, cmap='viridis', mincnt=1, linewidths=0.0,
                           extent=(lims[0], lims[1], lims[0], lims[1]), rasterized=True)
            hexes.append(hb)
            counts_max = max(counts_max, hb.get_array().max())
            ax.plot(lims, lims, ls='--', color='#333333', lw=1.0, alpha=0.7, zorder=3)
            ax.set_xlim(lims); ax.set_ylim(lims); ax.set_aspect('equal')
            ax.text(0.05, 0.95, f'r = {best_r:.3f}\nn = {xl.size:,}', transform=ax.transAxes,
                    fontsize=10, fontweight='bold', va='top', ha='left',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=accent,
                              linewidth=1.2, alpha=0.95))
            ax.set_xlabel(f'{a}\n(normalized abundance)', fontsize=10)
            ax.set_ylabel(f'{b}\n(normalized abundance)', fontsize=10)
            ax.set_title(short, fontsize=12, fontweight='bold', color=accent, pad=6)
            ax.spines[['top', 'right']].set_visible(False)
            ax.tick_params(labelsize=8)
        norm = LogNorm(vmin=1, vmax=counts_max)
        for hb in hexes:
            hb.set_norm(norm)
        cbar = fig.colorbar(hexes[-1], ax=axes, fraction=0.025, pad=0.01)
        cbar.set_label('proteins per bin', fontsize=9)
        cbar.ax.tick_params(labelsize=8)
        fig.suptitle(f'Replicate concordance (best pair per condition) — {MF_STAGE_NAMES[stage]}',
                     fontsize=12, fontweight='bold')
        save_fig(fig, out_dir, f'figS2_best_replicate_pair_density_{stage}')
        plt.show()


# ══════════════════════════════════════════════════════════════════════════════
# 5. Larval tissue colours (Fig 5, Fig 6, Fig S3, S5)
# ══════════════════════════════════════════════════════════════════════════════
TISSUE_COLORS_DARK  = {'gut': '#C0392B', 'brain': '#2471A3', 'mt': '#27AE60', 'discs': '#8E44AD'}
TISSUE_COLORS_LIGHT = {'gut': '#F1948A', 'brain': '#85C1E9', 'mt': '#82E0AA', 'discs': '#C39BD3'}
TISSUE_LABELS       = {'brain': 'BRAIN', 'mt': 'MT', 'discs': 'DISCS', 'gut': 'INTESTINE'}
