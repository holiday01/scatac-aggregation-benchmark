"""Rebuild the NSCLC gene-activity matrices after the GSM8069383 download was
found to be a duplicate copy of GSM8069382, and re-run every analysis that
depends on them.

The corrected file (MD5 826ccbbc53226c88a6469e17a2a2e6ff, 36,768,155 bytes) was
re-fetched from GEO; the other eight samples are byte-identical to fresh
downloads. GSM8069383 (Lu934) is a valid sample (4,818 cells, median 10,316
counts) and is now included; only GSM8069377 (Lu883, median 304 counts, 371
cells at the 1000-count threshold) remains excluded.

Window models replicate scripts/analysis/94_luad_prior_and_reference.py exactly.
"""
import glob, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp

BASE = Path("/mnt/10t/scrna_atac")
GTF = "/mnt/10t/holiday/hnsc_analysis/gencode.v38.annotation.gtf"
OUT = BASE / "data/processed/LUAD"
MIN_COUNTS = 1000
EXCLUDE = {"GSM8069377"}          # low coverage only; GSM8069383 is now included
EXTEND, DECAY = 100_000, 5000
MARKERS = {
    "Macrophage": ["CD68", "CSF1R", "MRC1", "MSR1", "AIF1", "ITGAM", "MARCO"],
    "B_cell":     ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19", "PAX5"],
    "T_NK":       ["CD3D", "CD3E", "CD2", "CD8A", "IL7R", "THEMIS", "NKG7", "GNLY", "KLRD1"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "CDH1", "NKX2-1", "SFTPC"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
}
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

log("parsing GENCODE v38 (protein-coding)")
rows = []
with open(GTF) as f:
    for line in f:
        if line[0] == "#": continue
        p = line.split("\t")
        if p[2] != "gene" or 'gene_type "protein_coding"' not in p[8]: continue
        chrom, s, e, strand = p[0], int(p[3]), int(p[4]), p[6]
        ps, pe = (s - 2000, e) if strand == "+" else (s, e + 2000)
        rows.append((chrom, s, e, ps, pe, p[8].split('gene_name "')[1].split('"')[0]))
genes = pd.DataFrame(rows, columns=["chrom", "gs", "ge", "ps", "pe", "gene"]).drop_duplicates("gene").reset_index(drop=True)
log(f"  genes: {len(genes)}")

prox_parts, enh_parts = [], []
for f in sorted(glob.glob(str(BASE / "data/raw/scATAC_LUAD/GSE270148/*.h5"))):
    gsm = os.path.basename(f).split("_")[0]
    if gsm in EXCLUDE:
        log(f"{gsm}: excluded (low coverage)"); continue
    a = sc.read_10x_h5(f, gex_only=False); a.var_names_make_unique()
    d = np.asarray(a.X.sum(1)).ravel()
    a = a[d >= MIN_COUNTS].copy()
    n = pd.Series(a.var_names)
    if {"chrom", "chromStart", "chromEnd"} <= set(a.var.columns):
        pc = a.var[["chrom", "chromStart", "chromEnd"]].copy(); pc.columns = ["chrom", "start", "end"]
    else:
        pc = pd.DataFrame({"chrom": n.str.split(":").str[0]})
        se = n.str.split(":").str[1].str.split("-")
        pc["start"] = se.str[0].astype(int); pc["end"] = se.str[1].astype(int)
    pc = pc.reset_index(drop=True)
    pc["start"] = pc.start.astype(int); pc["end"] = pc.end.astype(int)
    pc["mid"] = (pc.start + pc.end) // 2
    pc["idx"] = np.arange(len(pc))
    by = {c: g for c, g in pc.groupby("chrom")}
    pi, pj, ei, ej, ev = [], [], [], [], []
    for gi, g in enumerate(genes.itertuples(index=False)):
        ch = by.get(g.chrom)
        if ch is None: continue
        m = ch[(ch.end >= g.ps) & (ch.start <= g.pe)]
        pi.extend(m.idx.values); pj.extend([gi] * len(m))
        m2 = ch[(ch.mid >= g.gs - EXTEND) & (ch.mid <= g.ge + EXTEND)]
        if len(m2):
            dist = np.where(m2.mid.values < g.gs, g.gs - m2.mid.values,
                   np.where(m2.mid.values > g.ge, m2.mid.values - g.ge, 0))
            ei.extend(m2.idx.values); ej.extend([gi] * len(m2)); ev.extend(np.exp(-dist / DECAY))
    Wp = sp.csr_matrix((np.ones(len(pi)), (pi, pj)), shape=(len(pc), len(genes)))
    We = sp.csr_matrix((ev, (ei, ej)), shape=(len(pc), len(genes)))
    obs = [f"{gsm}_{b}" for b in a.obs_names]
    prox_parts.append(ad.AnnData(X=(a.X @ Wp).tocsr(), obs=pd.DataFrame(index=obs), var=pd.DataFrame(index=genes.gene.values)))
    enh_parts.append(ad.AnnData(X=(a.X @ We).tocsr(), obs=pd.DataFrame(index=obs), var=pd.DataFrame(index=genes.gene.values)))
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

for tag, parts in [("proximal", prox_parts), ("enhancer", enh_parts)]:
    A = ad.concat(parts)
    A.obs["sample"] = [i.split("_")[0] for i in A.obs_names]
    A.obs["cell_type"] = annotate(A).values
    A.var["gene_name"] = A.var_names
    A.write_h5ad(OUT / f"luad_geneactivity_{tag}.h5ad")
    log(f"  saved {tag}: {A.shape}  {A.obs.cell_type.value_counts().to_dict()}")
log("REBUILD DONE")
