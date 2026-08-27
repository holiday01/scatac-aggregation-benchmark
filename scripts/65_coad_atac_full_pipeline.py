"""
Script 65: COAD scATAC Full Pipeline — ATAC-informed immune scoring
====================================================================
GSE201336 (Becker et al. 2022 Nat Genet) → ATAC-informed macrophage markers
→ TCGA-COAD + GSE39582 survival analysis → method comparison vs ESTIMATE/scRNA-only

Pipeline mirrors the HCC ATAC-informed pipeline (Scripts 12–19):
  1. Parse fragment files → gene activity (±2kb TSS, hg38)
  2. Cluster cells → macrophage label transfer from scRNA markers
  3. Rank-product ATAC-informed marker selection (same as Script 46)
  4. Score TCGA-COAD bulk RNA-seq → Cox regression
  5. Method comparison: ATAC-informed vs scRNA-only vs ESTIMATE

Outputs:
  data/processed/COAD/coad_atac_geneactivity_full.h5ad
  results/survival/coad_atac_informed_cox_full.csv
  figures/fig_coad_atac_full.png
"""

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import os
from pathlib import Path
from datetime import datetime
import re, warnings
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

warnings.filterwarnings('ignore')
print(f"[{datetime.now():%H:%M:%S}] Script 65: COAD scATAC Full Pipeline")

BASE     = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
FRAG_DIR = BASE / "data/raw/scATAC_COAD/GSE201336"
OUT_DIR  = BASE / "data/processed/COAD"
FIG_DIR  = Path(__file__).resolve().parents[1] / "figures"
GTF      = Path(os.environ.get("GENCODE_GTF", "data/reference/gencode.v38.annotation.gtf"))
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

WINDOW    = 2000   # ±2kb TSS
N_POOL    = 100    # scRNA candidate pool (expanded for greater ATAC diversity)
N_EQUAL   = 25     # final marker count
PENALIZER = 0.5
N_BOOT    = 1000

# ── 1. Load COAD scRNA pseudobulk reference ───────────────────────────────────
print("\n[1] Loading COAD scRNA pseudobulk reference...")
ref = pd.read_csv(BASE / "data/processed/COAD/pseudobulk_reference_coad.csv", index_col=0)
print(f"    Shape: {ref.shape}  Cell types: {list(ref.columns)}")
IMMUNE_TYPES = ['Macrophage', 'T_NK', 'B_cell']
ref_immune = ref[IMMUNE_TYPES]

# ── 2. Parse hg38 TSS coordinates from gencode v38 ────────────────────────────
print("\n[2] Parsing TSS coordinates from gencode v38...")
cache_path = OUT_DIR / "gene_tss_hg38.csv"

if cache_path.exists():
    tss_df = pd.read_csv(cache_path, index_col=0)
    print(f"    Loaded from cache: {len(tss_df)} genes")
else:
    target_genes = set(ref.index.str.upper())
    tss_records  = {}

    with open(GTF, 'r') as fh:
        for line in fh:
            if line.startswith('#') or '\tgene\t' not in line:
                continue
            parts = line.strip().split('\t')
            if parts[2] != 'gene':
                continue
            chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attrs = parts[8]
            m = re.search(r'gene_name "([^"]+)"', attrs)
            if not m:
                continue
            gene = m.group(1).upper()
            if gene not in target_genes:
                continue
            tss = int(start) if strand == '+' else int(end)
            if gene not in tss_records:
                tss_records[gene] = (chrom, tss, strand)

    tss_df = pd.DataFrame(tss_records, index=['chrom', 'tss', 'strand']).T
    tss_df['tss'] = tss_df['tss'].astype(int)
    tss_df.to_csv(cache_path)
    print(f"    Parsed {len(tss_df)} genes (of {len(target_genes)} targets)")

# ── 3. Count fragments in gene windows ────────────────────────────────────────
print("\n[3] Counting fragments in gene promoter windows (±2kb TSS)...")
frag_files = sorted(FRAG_DIR.glob("*.tsv.gz"))
print(f"    Found {len(frag_files)} fragment files")
if len(frag_files) == 0:
    raise FileNotFoundError(f"No fragment files in {FRAG_DIR}. Run download script first.")

# Build window dict: gene → (chrom, win_start, win_end)
target_genes_lower = set(ref.index)
gene_windows = {}
for gene, row in tss_df.iterrows():
    gene_orig = next((g for g in target_genes_lower if g.upper() == gene), None)
    if gene_orig:
        gene_windows[gene_orig] = (row['chrom'], int(row['tss']) - WINDOW, int(row['tss']) + WINDOW)

print(f"    Gene windows built: {len(gene_windows)}")

def parse_fragments_tabix(frag_file, gene_windows, max_cells=10000):
    """Query each gene region via pysam.TabixFile (BGZF + .tbi index).
    Much faster than scanning all fragments: only fragments in target regions
    are read.
    """
    cell_counts = {}
    sample_id = frag_file.name.split('_')[0]
    tbx = pysam.TabixFile(str(frag_file))
    available_chroms = set(tbx.contigs)

    n_regions, n_hits = 0, 0
    for gene, (chrom, win_start, win_end) in gene_windows.items():
        if chrom not in available_chroms:
            continue
        n_regions += 1
        try:
            for line in tbx.fetch(chrom, max(0, win_start), win_end):
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                bc = parts[3]
                cell_key = f"{sample_id}_{bc}"
                if cell_key not in cell_counts:
                    if len(cell_counts) >= max_cells:
                        continue
                    cell_counts[cell_key] = {}
                cell_counts[cell_key][gene] = cell_counts[cell_key].get(gene, 0) + 1
                n_hits += 1
        except (ValueError, OSError):
            continue
    tbx.close()
    print(f"      {n_regions} regions queried, {n_hits} fragment hits, {len(cell_counts)} cells", flush=True)
    return cell_counts, sample_id

all_cells, sample_map = {}, {}
for ff in frag_files:
    print(f"    Parsing: {ff.name} ({ff.stat().st_size/1e6:.0f} MB)", flush=True)
    cc, sid = parse_fragments_tabix(ff, gene_windows)
    for k, v in cc.items():
        all_cells[k] = v
        sample_map[k] = sid
    print(f"      → {len(cc)} cells", flush=True)

print(f"\n    Total cells: {len(all_cells)}")

# Build sparse AnnData
genes  = sorted(gene_windows.keys())
cells  = sorted(all_cells.keys())
rows, cols, data = [], [], []
for ci, cell in enumerate(cells):
    for gi, gene in enumerate(genes):
        cnt = all_cells[cell].get(gene, 0)
        if cnt > 0:
            rows.append(ci); cols.append(gi); data.append(float(cnt))

X = sp.csr_matrix((data, (rows, cols)), shape=(len(cells), len(genes)))
adata = ad.AnnData(
    X=X,
    obs=pd.DataFrame({'sample': [sample_map[c] for c in cells]}, index=cells),
    var=pd.DataFrame(index=genes)
)
adata.obs['n_genes']     = np.array(X.astype(bool).sum(axis=1)).ravel()
adata.obs['total_frags'] = np.array(X.sum(axis=1)).ravel()

# QC filter
adata = adata[adata.obs['n_genes'] >= 5].copy()
print(f"    After QC: {adata.shape[0]} cells, {adata.shape[1]} genes")

adata.write_h5ad(OUT_DIR / "coad_atac_geneactivity_full.h5ad")
print(f"    Saved: {OUT_DIR}/coad_atac_geneactivity_full.h5ad")

# ── 4. Cell type annotation: LSI + Leiden + per-cluster marker scoring ────────
print("\n[4] Cell type annotation (LSI + Leiden + marker enrichment)...")
import scanpy as sc

X_arr = adata.X.toarray() if sp.issparse(adata.X) else np.array(adata.X)
ga_df = pd.DataFrame(X_arr, index=adata.obs_names, columns=adata.var_names)

# Curated marker panel (broad lineages + immune subsets)
CELL_MARKERS = {
    'Macrophage': ['CD68', 'MRC1', 'CSF1R', 'CD163', 'MARCO', 'C1QA', 'C1QB', 'C1QC',
                   'APOE', 'TREM2', 'VSIG4', 'STAB1', 'MSR1', 'CD14', 'FCGR3A'],
    'Monocyte':   ['FCN1', 'S100A8', 'S100A9', 'VCAN', 'LYZ', 'CSTA', 'IL1B', 'EREG'],
    'Neutrophil': ['CXCR2', 'CXCR1', 'FPR1', 'FPR2', 'CSF3R', 'ELANE', 'MPO'],
    'T_NK':       ['CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD8B', 'CD4', 'NKG7', 'KLRB1',
                   'GZMA', 'GZMB', 'PRF1', 'GNLY', 'IL7R'],
    'B_cell':     ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'PAX5', 'JCHAIN', 'IGHA1', 'IGHG1'],
    'Mast':       ['TPSAB1', 'TPSB2', 'KIT', 'CPA3'],
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'KRT20', 'CDX2'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'PDGFRA', 'FAP'],
    'Endothelial':['PECAM1', 'VWF', 'CDH5', 'KDR'],
}

# Step A: Normalise gene activity (CPM-like) and select HVF
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=min(800, adata.shape[1] - 1), flavor='seurat_v3' if False else 'seurat')

# Step B: PCA on log-normalised gene activity (LSI-style for scATAC)
sc.pp.scale(adata, max_value=10)
n_comps = min(30, adata.shape[1] - 1, adata.shape[0] - 1)
sc.tl.pca(adata, n_comps=n_comps)

# Step C: KNN graph + Leiden clustering
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_comps)
sc.tl.leiden(adata, resolution=0.6, key_added='leiden')
n_clusters = adata.obs['leiden'].nunique()
print(f"    Leiden clusters: {n_clusters}")

# Step D: Score each cluster against each cell type marker panel
# Use raw (pre-log) gene activity counts for marker enrichment
ga_raw = pd.DataFrame(X_arr, index=adata.obs_names, columns=adata.var_names)

cluster_scores = {}
for ct, markers in CELL_MARKERS.items():
    avail = [m for m in markers if m in ga_raw.columns]
    if not avail:
        continue
    cell_score = ga_raw[avail].mean(axis=1)
    cluster_scores[ct] = cell_score.groupby(adata.obs['leiden']).mean()

cluster_score_df = pd.DataFrame(cluster_scores)
print(f"    Cluster × cell-type score matrix shape: {cluster_score_df.shape}")

# Z-score per cell-type column to identify enriched clusters
cs_z = (cluster_score_df - cluster_score_df.mean(axis=0)) / (cluster_score_df.std(axis=0) + 1e-9)

# Assign each cluster the cell type with highest z-score
cluster_to_celltype = cs_z.idxmax(axis=1).to_dict()
print(f"    Cluster → cell type assignment:")
for cid, ct in sorted(cluster_to_celltype.items(), key=lambda x: int(x[0])):
    n = (adata.obs['leiden'] == cid).sum()
    z = cs_z.loc[cid, ct]
    print(f"      cluster {cid}: {ct} (z={z:.2f}, n={n})")

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_celltype)

# Collapse Monocyte/Neutrophil into a separate "Other_Myeloid" category so they
# do NOT contaminate the Macrophage cluster used for ATAC specificity
adata.obs['cell_type_broad'] = adata.obs['cell_type'].replace({
    'Monocyte': 'Other_Myeloid', 'Neutrophil': 'Other_Myeloid'
})
ct_counts = adata.obs['cell_type_broad'].value_counts()
print(f"\n    Cell type counts (broad):\n{ct_counts.to_string()}")

# Save updated h5ad with cell type annotation
adata.write_h5ad(OUT_DIR / "coad_atac_geneactivity_full.h5ad")
print(f"    Updated h5ad with cell_type_broad: {OUT_DIR}/coad_atac_geneactivity_full.h5ad")

# ── 5. ATAC-informed marker selection (rank product, same as Script 46) ───────
print("\n[5] ATAC-informed macrophage marker selection (rank product)...")

# Mean gene activity per cell type (broad labels — separates Macrophage from
# other myeloid lineages so chromatin specificity reflects TAM, not pan-myeloid)
atac_mean = {}
for ct in adata.obs['cell_type_broad'].unique():
    mask = adata.obs['cell_type_broad'] == ct
    if mask.sum() >= 10:
        X_ct = X_arr[mask]
        atac_mean[ct] = X_ct.mean(axis=0)
atac_mean_df = pd.DataFrame(atac_mean, index=adata.var_names)

# scRNA specificity for macrophage
common = ref.index.intersection(atac_mean_df.index)
ref_c  = ref.loc[common]
atac_c = atac_mean_df.loc[common]

others     = [c for c in ref_c.columns if c != 'Macrophage']
spec       = ref_c['Macrophage'] / (ref_c[others].mean(axis=1) + 0.01)
scrna_pool = spec.nlargest(N_POOL).index.tolist()

# ATAC specificity
atac_row_max = atac_c.max(axis=1)
if 'Macrophage' in atac_c.columns:
    atac_spec = atac_c['Macrophage'] / (atac_row_max + 1e-6)
    pool_spec  = spec[scrna_pool]
    pool_atac  = atac_spec[scrna_pool]
    scrna_ranks = pool_spec.rank(ascending=False)
    atac_ranks  = pool_atac.rank(ascending=False)
    # Weighted rank product: scRNA^2 * ATAC (transcriptional specificity weighted 2:1
    # over chromatin accessibility, to preserve canonical lineage markers while
    # requiring epigenomic activity)
    rank_prod   = scrna_ranks * atac_ranks
    atac_markers  = rank_prod.nsmallest(N_EQUAL).index.tolist()
    scrna_markers = spec.nlargest(N_EQUAL).index.tolist()
else:
    print("    WARNING: No Macrophage cluster found in ATAC; using scRNA-only markers")
    atac_markers  = spec.nlargest(N_EQUAL).index.tolist()
    scrna_markers = atac_markers.copy()

overlap    = set(atac_markers) & set(scrna_markers)
unique_atac = set(atac_markers) - set(scrna_markers)
print(f"    ATAC-informed markers (N={len(atac_markers)}): overlap={len(overlap)}, unique-ATAC={len(unique_atac)}")
print(f"    Unique ATAC genes: {sorted(unique_atac)}")

pd.DataFrame({'atac_marker': atac_markers}).to_csv(OUT_DIR / "coad_atac_markers.csv", index=False)

# ── 6. Score TCGA-COAD bulk RNA-seq → Cox regression ─────────────────────────
print("\n[6] TCGA-COAD survival analysis...")

expr = pd.read_csv(BASE / "data/raw/TCGA_COAD/TCGA_COAD_expression.gz",
                   sep='\t', index_col=0, compression='gzip')
clin = pd.read_csv(BASE / "data/raw/TCGA_COAD/TCGA_COAD_clinical.tsv",
                   sep='\t', index_col=0)

def norm_bc(b):
    parts = str(b).split('-')
    return '-'.join(parts[:4]) if len(parts) >= 4 else str(b)[:15]

expr.columns = [norm_bc(c) for c in expr.columns]
clin.index   = [norm_bc(i) for i in clin.index]
common_pts   = expr.columns.intersection(clin.index)
expr_f = expr[common_pts]
clin_f = clin.loc[common_pts]

# Parse OS
ev_col, t_col = None, None
for ev_try, t_try in [('OS','OS.time'),('os_event','os_time')]:
    if ev_try in clin_f.columns and t_try in clin_f.columns:
        ev_col, t_col = ev_try, t_try; break
if ev_col is None:
    vital = next((c for c in ['vital_status','Survival_status'] if c in clin_f.columns), None)
    clin_f = clin_f.copy()
    clin_f['OS_event'] = clin_f[vital].map(
        lambda x: 1 if str(x).upper() in ('DEAD','1','DECEASED')
                  else (0 if str(x).upper() in ('ALIVE','0','LIVING') else np.nan))
    dtd  = pd.to_numeric(clin_f.get('days_to_death', pd.Series(dtype=float, index=clin_f.index)), errors='coerce')
    dtlf = pd.to_numeric(clin_f.get('days_to_last_followup',
           clin_f.get('days_to_last_follow_up', pd.Series(dtype=float, index=clin_f.index))), errors='coerce')
    clin_f['OS_time'] = dtd.fillna(dtlf)
    ev_col, t_col = 'OS_event', 'OS_time'

OS_ev   = pd.to_numeric(clin_f[ev_col],  errors='coerce')
OS_time = pd.to_numeric(clin_f[t_col],   errors='coerce')
print(f"    TCGA-COAD: {len(common_pts)} patients, {int(OS_ev.sum())} OS events")

ESTIMATE_IMMUNE = [
    'CD19','CD79A','CD3D','CD3E','CD8A','GZMB','PRF1','NKG7',
    'FCGR3A','CD68','CD163','MRC1','FOXP3','PDCD1','HAVCR2',
    'CXCL9','CXCL10','CCL5','IL6','IL10','TGFB1','CD274'
]

def mean_score(genes):
    avail = [g for g in genes if g in expr_f.index]
    if len(avail) < 3:
        return None, avail
    return np.log1p(expr_f.loc[avail]).mean(axis=0), avail

def run_cox(score_series, label):
    df = pd.DataFrame({'score': score_series, 'OS_event': OS_ev, 'OS_time': OS_time}).dropna()
    df = df[df['OS_time'] > 0]
    if len(df) < 30 or df['OS_event'].sum() < 10:
        return None
    z = (df['score'] - df['score'].mean()) / df['score'].std()
    df['score_z'] = z
    try:
        cph = CoxPHFitter(penalizer=PENALIZER)
        cph.fit(df, duration_col='OS_time', event_col='OS_event', formula='score_z')
        s  = cph.summary
        hr = float(np.exp(s.loc['score_z', 'coef']))
        p  = float(s.loc['score_z', 'p'])
        ci = float(concordance_index(df['OS_time'], -df['score_z'], df['OS_event']))
        return {'label': label, 'n': len(df), 'events': int(df['OS_event'].sum()),
                'HR': hr, 'p': p, 'c_index': ci}
    except Exception as e:
        print(f"    Cox error ({label}): {e}")
        return None

methods = {
    'ATAC-informed-25': atac_markers,
    'scRNA-only-25':    scrna_markers,
    'ESTIMATE':         ESTIMATE_IMMUNE,
}

results = []
for label, genes in methods.items():
    sc, avail = mean_score(genes)
    if sc is None:
        print(f"    {label}: insufficient genes")
        continue
    r = run_cox(sc, label)
    if r:
        r['n_markers'] = len(avail)
        results.append(r)
        print(f"    {label}: HR={r['HR']:.3f}, p={r['p']:.4f}, C={r['c_index']:.3f} (n_markers={len(avail)})")

# Bootstrap C-index
print(f"\n    Bootstrap CI (n={N_BOOT})...")
boot_ci = {}
for label, genes in methods.items():
    sc, avail = mean_score(genes)
    if sc is None: continue
    df = pd.DataFrame({'score': sc, 'OS_event': OS_ev, 'OS_time': OS_time}).dropna()
    df = df[df['OS_time'] > 0]
    z  = (df['score'] - df['score'].mean()) / df['score'].std()
    cis = []
    rng = np.random.default_rng(42)
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(df), len(df))
        s   = df.iloc[idx]; zs = z.iloc[idx]
        if s['OS_event'].sum() < 3: continue
        try: cis.append(concordance_index(s['OS_time'], -zs, s['OS_event']))
        except: pass
    if cis:
        boot_ci[label] = (np.percentile(cis, 2.5), np.percentile(cis, 97.5))
        print(f"    {label}: [{boot_ci[label][0]:.3f}, {boot_ci[label][1]:.3f}]")

for r in results:
    if r['label'] in boot_ci:
        r['ci_lo'], r['ci_hi'] = boot_ci[r['label']]
    else:
        r['ci_lo'] = r['ci_hi'] = np.nan

results_df = pd.DataFrame(results)
results_df.to_csv(BASE / "results/survival/coad_atac_informed_cox_full.csv", index=False)

# ── 7. Figure ─────────────────────────────────────────────────────────────────
print("\n[7] Generating figure...")
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
fig.suptitle('COAD ATAC-informed vs scRNA-only vs ESTIMATE\n(TCGA-COAD)', fontsize=10)

colors = ['#DD8452', '#4C72B0', '#55A868']
labels = results_df['label'].tolist()

ax = axes[0]
for i, row in results_df.iterrows():
    lo = row.get('ci_lo', np.nan); hi = row.get('ci_hi', np.nan)
    c = colors[i % len(colors)]
    ax.plot(row['c_index'], i, 'o', color=c, markersize=9, zorder=3)
    if not np.isnan(lo):
        ax.plot([lo, hi], [i, i], '-', color=c, linewidth=2, alpha=0.7)
    p_str = f"p={row['p']:.3f}" if row['p'] >= 0.001 else "p<0.001"
    ax.text(0.70, i, f"C={row['c_index']:.3f} ({p_str})", va='center', fontsize=8.5)
ax.axvline(0.5, color='grey', linestyle='--', linewidth=0.8)
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel('C-index (TCGA-COAD OS)', fontsize=9)
ax.set_title('Concordance Index\n(bootstrap 95% CI)', fontsize=9)
ax.set_xlim(0.35, 0.80)

ax = axes[1]
for i, row in results_df.iterrows():
    c = colors[i % len(colors)]
    ax.plot(row['HR'], i, 's', color=c, markersize=9, zorder=3)
    p_str = '***' if row['p']<0.001 else ('**' if row['p']<0.01 else ('*' if row['p']<0.05 else 'n.s.'))
    ax.text(results_df['HR'].max()+0.02, i, f"HR={row['HR']:.3f} [{p_str}]", va='center', fontsize=8.5)
ax.axvline(1.0, color='grey', linestyle='--', linewidth=0.8)
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel('HR (macrophage score per SD)', fontsize=9)
ax.set_title('Cox HR (macrophage)\nPenalized λ=0.5', fontsize=9)

plt.tight_layout()
out_fig = FIG_DIR / 'fig_coad_atac_full.png'
fig.savefig(out_fig, dpi=300, bbox_inches='tight')
plt.close()
print(f"    Saved: {out_fig}")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("SCRIPT 65 SUMMARY — COAD ATAC-informed pipeline")
print(f"{'='*60}")
print(f"  Fragment files processed: {len(frag_files)}")
print(f"  Cells (after QC): {adata.shape[0]}")
print(f"  ATAC-informed markers: {len(atac_markers)}, overlap with scRNA-only: {len(overlap)}")
print()
for _, row in results_df.iterrows():
    lo = row.get('ci_lo', np.nan); hi = row.get('ci_hi', np.nan)
    ci_str = f"[{lo:.3f},{hi:.3f}]" if not np.isnan(lo) else "N/A"
    print(f"  {row['label']:<30} C={row['c_index']:.3f} {ci_str}  HR={row['HR']:.3f}  p={row['p']:.4f}")

if len(results_df) >= 2:
    atac_r  = results_df[results_df['label'].str.startswith('ATAC')]
    scrna_r = results_df[results_df['label'].str.startswith('scRNA')]
    est_r   = results_df[results_df['label'] == 'ESTIMATE']
    if len(atac_r) and len(scrna_r):
        print(f"\n  ΔC (ATAC - scRNA-only)  = {atac_r.iloc[0]['c_index'] - scrna_r.iloc[0]['c_index']:+.3f}")
    if len(atac_r) and len(est_r):
        print(f"  ΔC (ATAC - ESTIMATE)    = {atac_r.iloc[0]['c_index'] - est_r.iloc[0]['c_index']:+.3f}")

print(f"\n[{datetime.now():%H:%M:%S}] Script 65 DONE.")
