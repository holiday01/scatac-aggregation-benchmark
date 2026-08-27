"""
Script 92: Rebuild COAD scATAC gene activity from fragments — correctly.
========================================================================
The canonical COAD prior `coad_atac_geneactivity_full.h5ad` (Script 65) is
broken three ways:
  1. Z-SCORED. Script 65 line 209 runs `sc.pp.scale(adata, max_value=10)` for
     clustering and then SAVES that matrix. 91.7% of entries are negative.
     Script 76 then computes atac_spec = mean[target]/(max+1e-6) on z-scores —
     meaningless. 88.4% of all genes end up with argmax = Mast.
  2. TRUNCATED GENE UNIVERSE. Only 1,896 genes, not genome-wide.
  3. Consequently the COAD arm of the paper (3 of 9 combinations) carries no
     valid chromatin signal.

This rebuilds from the 6 GSE201336 fragment files with snapatac2, genome-wide,
NO scaling, using the same two models as Script 90 (HCC):
     proximal = gene body + 2 kb upstream (Signac model)
     enhancer = ArchR-style distance-weighted, 5 kb tiles, +-100 kb, exp decay

Cell set + labels are taken from the existing COAD h5ad so the comparison is
like-for-like (same cells, same annotations — only the prior changes).

Outputs: data/processed/COAD/coad_geneactivity_proximal.h5ad
         data/processed/COAD/coad_geneactivity_enhancer.h5ad
"""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
import snapatac2 as snap
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

BASE    = Path("/mnt/10t/scrna_atac")
FRAGDIR = BASE / "data/raw/scATAC_COAD/GSE201336"
GTF     = Path("/mnt/10t/holiday/hnsc_analysis/gencode.v38.annotation.gtf")
OLD     = BASE / "data/processed/COAD/coad_atac_geneactivity_full.h5ad"
OUT_DIR = BASE / "data/processed/COAD"
WORK    = OUT_DIR / "_s92_work"
WORK.mkdir(parents=True, exist_ok=True)

BIN, EXTEND, DECAY = 5000, 100_000, 5000


def log(m):
    print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)


old = ad.read_h5ad(OLD)
meta = old.obs[['sample', 'cell_type_broad']].copy()
meta['bc'] = [n.split('_', 1)[1] for n in old.obs_names]   # 'GSM..._AAAC-1' -> 'AAAC-1'
log(f"COAD reference cells: {len(meta)} across {meta['sample'].nunique()} samples")

frag_files = sorted(FRAGDIR.glob("*fragments.tsv.gz"))
log(f"fragment files: {len(frag_files)}")

prox_parts, enh_parts = [], []

for f in frag_files:
    gsm = f.name.split('_')[0]
    sub = meta[meta['sample'] == gsm]
    if len(sub) == 0:
        log(f"  {gsm}: no cells in reference — skipped"); continue
    wl = list(sub['bc'])
    log(f"  {gsm}: {len(wl)} cells")

    frag_h5 = WORK / f"{gsm}.h5ad"
    if not frag_h5.exists():
        d = snap.pp.import_fragments(f, chrom_sizes=snap.genome.hg38, file=str(frag_h5),
                                     whitelist=wl, min_num_fragments=0,
                                     sorted_by_barcode=False, n_jobs=8)
        d.close()

    # proximal (Signac model)
    d = snap.read(str(frag_h5), backed='r')
    gm = snap.pp.make_gene_matrix(d, gene_anno=GTF, upstream=2000, downstream=0,
                                  include_gene_body=True, file=str(WORK / f"{gsm}_gm.h5ad"))
    Xp = sp.csr_matrix(gm.X[:]); pv = list(gm.var_names); po = list(gm.obs_names)
    gm.close(); d.close()
    ap = ad.AnnData(X=Xp, obs=pd.DataFrame(index=[f"{gsm}_{b}" for b in po]),
                    var=pd.DataFrame(index=pv))
    prox_parts.append(ap)

    # tiles -> enhancer-aware
    d = snap.read(str(frag_h5), backed='r+')
    if d.X is None or d.shape[1] == 0:
        snap.pp.add_tile_matrix(d, bin_size=BIN, exclude_chroms=['chrM', 'chrY', 'M', 'Y'], n_jobs=8)
    T = sp.csr_matrix(d.X[:]); tn = list(d.var_names); to = list(d.obs_names)
    d.close()
    enh_parts.append((T, tn, [f"{gsm}_{b}" for b in to]))
    (WORK / f"{gsm}_gm.h5ad").unlink(missing_ok=True)

# ── proximal: concat ─────────────────────────────────────────────────────────
log("Concatenating proximal parts")
prox = ad.concat(prox_parts, join='outer', fill_value=0)
prox.obs['cell_type'] = old.obs['cell_type_broad'].reindex(prox.obs_names).values
prox.obs['sample']    = [n.split('_', 1)[0] for n in prox.obs_names]
prox = prox[~prox.obs['cell_type'].isna()].copy()
prox.var['gene_name'] = prox.var_names
prox.write_h5ad(OUT_DIR / "coad_geneactivity_proximal.h5ad")
log(f"  saved proximal: {prox.shape}")

# ── enhancer: build W once on the union tile space, then per-part multiply ───
log("Parsing gene bodies (protein-coding, GENCODE v38)")
rows = []
with open(GTF) as fh:
    for line in fh:
        if line[0] == '#':
            continue
        p = line.split('\t')
        if p[2] != 'gene' or 'gene_type "protein_coding"' not in p[8]:
            continue
        rows.append((p[0], int(p[3]), int(p[4]), p[8].split('gene_name "')[1].split('"')[0]))
genes = pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'gene']).drop_duplicates('gene')
log(f"  genes: {len(genes)}")

enh_mats, enh_obs = [], []
for T, tn, obs in enh_parts:
    tc = pd.DataFrame({'name': tn})
    tc['chrom'] = tc.name.str.split(':').str[0]
    se = tc.name.str.split(':').str[1].str.split('-')
    tc['mid'] = (se.str[0].astype(int) + se.str[1].astype(int)) // 2
    tc['idx'] = np.arange(len(tc))
    by = {c: g.sort_values('mid') for c, g in tc.groupby('chrom')}

    wi, wj, wv = [], [], []
    for gi, g in enumerate(genes.itertuples(index=False)):
        ch = by.get(g.chrom)
        if ch is None:
            continue
        m = ch[(ch.mid >= g.start - EXTEND) & (ch.mid <= g.end + EXTEND)]
        if len(m) == 0:
            continue
        d = np.where(m.mid.values < g.start, g.start - m.mid.values,
            np.where(m.mid.values > g.end,   m.mid.values - g.end, 0))
        w = np.exp(-d / DECAY)
        wi.extend(m.idx.values); wj.extend([gi] * len(m)); wv.extend(w)
    W = sp.csr_matrix((wv, (wi, wj)), shape=(len(tc), len(genes)))
    enh_mats.append((T @ W).tocsr())
    enh_obs.extend(obs)

E = sp.vstack(enh_mats).tocsr()
enh = ad.AnnData(X=E, obs=pd.DataFrame(index=enh_obs),
                 var=pd.DataFrame(index=genes.gene.values))
enh.obs['cell_type'] = old.obs['cell_type_broad'].reindex(enh.obs_names).values
enh.obs['sample']    = [n.split('_', 1)[0] for n in enh.obs_names]
enh = enh[~enh.obs['cell_type'].isna()].copy()
enh.var['gene_name'] = enh.var_names
enh.var_names_make_unique()
enh.write_h5ad(OUT_DIR / "coad_geneactivity_enhancer.h5ad")
log(f"  saved enhancer: {enh.shape}")

# ── validate ─────────────────────────────────────────────────────────────────
LIN = {'CD3D': 'T_NK', 'CD3E': 'T_NK', 'CD2': 'T_NK',
       'MS4A1': 'B_cell', 'CD79A': 'B_cell', 'CD79B': 'B_cell', 'BANK1': 'B_cell',
       'CSF1R': 'Macrophage', 'CD68': 'Macrophage', 'MSR1': 'Macrophage',
       'CD163': 'Macrophage', 'AIF1': 'Macrophage',
       'EPCAM': 'Epithelial', 'KRT8': 'Epithelial', 'CDH1': 'Epithelial'}

for tag, a in [('old_zscored', old), ('new_proximal', prox), ('new_enhancer', enh)]:
    col = 'cell_type_broad' if 'cell_type_broad' in a.obs else 'cell_type'
    X = a.X.tocsc() if sp.issparse(a.X) else a.X
    ct = a.obs[col].astype(str).values
    hit = 0; tot = 0
    for g, exp in LIN.items():
        if g not in a.var_names:
            continue
        v = X[:, a.var_names.get_loc(g)]
        v = np.asarray(v.todense()).ravel() if sp.issparse(v) else np.asarray(v).ravel()
        am = pd.Series(v).groupby(ct).mean().idxmax()
        hit += int(am == exp); tot += 1
    print(f"  {tag:13} markers in correct lineage: {hit}/{tot}")

log("Script 92 done.")
