"""
Script 94: Build the LUAD chromatin prior + scRNA reference (COAD replacement)
=============================================================================
COAD is void (Script 92: 0/1434 macrophages, 0/1759 T_NK survive >=1000-frag QC).
LUAD GSE270148 is the replacement (Script 93: 7,474 macrophage / 8,080 B /
15,351 T-NK REAL cells; 0 of 8 patients show the "B>70% artefact" the manuscript
used to reject it).

Builds, genome-wide, from each sample's 10x peak matrix (peak coords in the h5;
no fragments needed):
  proximal = peaks overlapping gene body + 2 kb upstream   (Signac model)
  enhancer = distance-weighted over peaks within +-100 kb  (ArchR-style, tau=5kb)
and the matching scRNA pseudobulk reference from GSE131907.

Excluded: GSM8069377 only (median 304 counts per cell). NOTE 2026-08-27:
GSM8069383 was previously excluded on the basis of a local download that was a
duplicate copy of GSM8069382; the file was re-fetched from GEO and the sample is
valid, so it is now included.
          GSM8069377 (median 304 counts/cell).

Outputs: data/processed/LUAD/luad_geneactivity_proximal.h5ad
         data/processed/LUAD/luad_geneactivity_enhancer.h5ad
         data/processed/LUAD/pseudobulk_reference_luad.csv
"""
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp
from pathlib import Path
from datetime import datetime
import glob, os, warnings
warnings.filterwarnings('ignore')

BASE = Path("/mnt/10t/scrna_atac")
GTF  = "/mnt/10t/holiday/hnsc_analysis/gencode.v38.annotation.gtf"
OUT  = BASE / "data/processed/LUAD"
OUT.mkdir(parents=True, exist_ok=True)

MIN_COUNTS = 1000
EXCLUDE    = {'GSM8069377'}                 # low coverage (median 304 counts/cell)
EXTEND, DECAY = 100_000, 5000

MARKERS = {
    'Macrophage': ['CD68', 'CSF1R', 'MRC1', 'MSR1', 'AIF1', 'ITGAM', 'MARCO'],
    'B_cell':     ['MS4A1', 'CD79A', 'CD79B', 'BANK1', 'CD19', 'PAX5'],
    'T_NK':       ['CD3D', 'CD3E', 'CD2', 'CD8A', 'IL7R', 'THEMIS', 'NKG7', 'GNLY', 'KLRD1'],
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'CDH1', 'NKX2-1', 'SFTPC'],
    'Endothelial':['PECAM1', 'VWF', 'CDH5'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'],
}


def log(m):
    print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)


# ── gene models ──────────────────────────────────────────────────────────────
log("Parsing GENCODE v38 (protein-coding)")
rows = []
with open(GTF) as f:
    for line in f:
        if line[0] == '#':
            continue
        p = line.split('\t')
        if p[2] != 'gene' or 'gene_type "protein_coding"' not in p[8]:
            continue
        chrom, s, e, strand = p[0], int(p[3]), int(p[4]), p[6]
        ps, pe = (s - 2000, e) if strand == '+' else (s, e + 2000)   # +2kb upstream
        rows.append((chrom, s, e, ps, pe, p[8].split('gene_name "')[1].split('"')[0]))
genes = pd.DataFrame(rows, columns=['chrom', 'gs', 'ge', 'ps', 'pe', 'gene']).drop_duplicates('gene').reset_index(drop=True)
log(f"  genes: {len(genes)}")

prox_parts, enh_parts = [], []

for f in sorted(glob.glob(str(BASE / "data/raw/scATAC_LUAD/GSE270148/*.h5"))):
    gsm = os.path.basename(f).split('_')[0]
    if gsm in EXCLUDE:
        log(f"{gsm}: excluded"); continue
    a = sc.read_10x_h5(f, gex_only=False); a.var_names_make_unique()
    d = np.asarray(a.X.sum(1)).ravel()
    a = a[d >= MIN_COUNTS].copy()

    n = pd.Series(a.var_names)
    if {'chrom', 'chromStart', 'chromEnd'} <= set(a.var.columns):
        pc = a.var[['chrom', 'chromStart', 'chromEnd']].copy(); pc.columns = ['chrom', 'start', 'end']
    else:
        pc = pd.DataFrame({'chrom': n.str.split(':').str[0]})
        se = n.str.split(':').str[1].str.split('-')
        pc['start'] = se.str[0].astype(int); pc['end'] = se.str[1].astype(int)
    pc = pc.reset_index(drop=True)
    pc['start'] = pc.start.astype(int); pc['end'] = pc.end.astype(int)
    pc['mid'] = (pc.start + pc.end) // 2
    pc['idx'] = np.arange(len(pc))
    by = {c: g for c, g in pc.groupby('chrom')}

    pi, pj = [], []                       # proximal overlap
    ei, ej, ev = [], [], []               # enhancer distance-weighted
    for gi, g in enumerate(genes.itertuples(index=False)):
        ch = by.get(g.chrom)
        if ch is None:
            continue
        m = ch[(ch.end >= g.ps) & (ch.start <= g.pe)]
        pi.extend(m.idx.values); pj.extend([gi] * len(m))

        m2 = ch[(ch.mid >= g.gs - EXTEND) & (ch.mid <= g.ge + EXTEND)]
        if len(m2):
            dist = np.where(m2.mid.values < g.gs, g.gs - m2.mid.values,
                   np.where(m2.mid.values > g.ge, m2.mid.values - g.ge, 0))
            ei.extend(m2.idx.values); ej.extend([gi] * len(m2))
            ev.extend(np.exp(-dist / DECAY))

    Wp = sp.csr_matrix((np.ones(len(pi)), (pi, pj)), shape=(len(pc), len(genes)))
    We = sp.csr_matrix((ev, (ei, ej)), shape=(len(pc), len(genes)))
    obs = [f"{gsm}_{b}" for b in a.obs_names]

    prox_parts.append(ad.AnnData(X=(a.X @ Wp).tocsr(),
                                 obs=pd.DataFrame(index=obs),
                                 var=pd.DataFrame(index=genes.gene.values)))
    enh_parts.append(ad.AnnData(X=(a.X @ We).tocsr(),
                                obs=pd.DataFrame(index=obs),
                                var=pd.DataFrame(index=genes.gene.values)))
    log(f"  {gsm}: {a.shape[0]} cells, {len(pc)} peaks")


def annotate(A):
    X = A.X
    d = np.asarray(X.sum(1)).ravel(); d[d == 0] = 1
    sub = {}
    for lin, gl in MARKERS.items():
        cols = [g for g in gl if g in A.var_names]
        idx = [A.var_names.get_loc(g) for g in cols]
        v = np.asarray(X[:, idx].todense()) / d[:, None] * 1e4
        v = (v - v.mean(0)) / (v.std(0) + 1e-9)
        sub[lin] = v.mean(1)
    return pd.DataFrame(sub, index=A.obs_names).idxmax(axis=1)


for tag, parts in [('proximal', prox_parts), ('enhancer', enh_parts)]:
    A = ad.concat(parts)
    A.obs['sample'] = [i.split('_')[0] for i in A.obs_names]
    A.obs['cell_type'] = annotate(A).values
    A.var['gene_name'] = A.var_names
    A.write_h5ad(OUT / f"luad_geneactivity_{tag}.h5ad")
    log(f"  saved {tag}: {A.shape}  {A.obs.cell_type.value_counts().to_dict()}")

# ── scRNA pseudobulk reference (GSE131907) ───────────────────────────────────
log("Building LUAD scRNA pseudobulk reference")
r = ad.read_h5ad(BASE / "data/processed/LUAD/luad_scrna_annotated.h5ad")
MAP = {'Macrophage': 'Macrophage', 'B lymphocytes': 'B_cell', 'T_NK': 'T_NK',
       'Epithelial': 'Epithelial', 'Fibroblast': 'Fibroblast', 'Endothelial': 'Endothelial'}
r.obs['ct'] = r.obs['celltype_std'].map(MAP)
r = r[~r.obs['ct'].isna()].copy()
if 'log1p' not in r.uns_keys():
    sc.pp.normalize_total(r, target_sum=1e4); sc.pp.log1p(r)
X = r.X.toarray() if sp.issparse(r.X) else np.asarray(r.X)
pb = pd.DataFrame({ct: X[(r.obs.ct == ct).values].mean(0) for ct in r.obs.ct.unique()},
                  index=r.var_names)
pb.to_csv(OUT / "pseudobulk_reference_luad.csv")
log(f"  saved reference: {pb.shape}  cols={list(pb.columns)}")
log("Script 94 done.")
