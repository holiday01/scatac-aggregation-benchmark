"""
01: Build ATAC gene-activity matrices for the 10x PBMC multiome (pbmc_granulocyte_sorted_10k)
    from the fragment file, restricted to the cells of the filtered feature-barcode matrix.

  proximal  : snapatac2 make_gene_matrix(upstream=2000, downstream=0, include_gene_body=True)
              == Signac GeneActivity model (same as scripts/analysis/90 stage B)
  enhancer  : ArchR-style distance-weighted gene score over 5 kb tiles,
              w = exp(-d/5kb) within 100 kb of the gene body, full weight inside the body
              (same as scripts/analysis/90 stage C)

True per-cell depth = n_fragment recorded by snapatac2 at import (never a matrix row sum).

Out: data/processed/multiome_pbmc10k/{fragments.h5ad, pbmc10k_geneactivity_proximal.h5ad,
                                       pbmc10k_geneactivity_enhancer.h5ad}
"""
import sys, time, warnings
import os
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, hdf5plugin  # noqa
import scanpy as sc
import snapatac2 as snap
warnings.filterwarnings("ignore")

BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
RAW  = BASE / "data/raw/multiome_pbmc10k"
H5   = RAW / "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
FRAG = RAW / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
GTF  = Path(os.environ.get("GENCODE_GTF", "data/reference/gencode.v38.annotation.gtf"))
OUT  = BASE / "data/processed/multiome_pbmc10k"; OUT.mkdir(parents=True, exist_ok=True)
FRAG_H5 = OUT / "fragments.h5ad"
BIN, EXTEND, DECAY = 5000, 100_000, 5000
STAGE = sys.argv[1] if len(sys.argv) > 1 else "all"

def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

mm = sc.read_10x_h5(H5, gex_only=False)
barcodes = list(mm.obs_names)
log(f"filtered matrix: {mm.shape}; {len(barcodes)} cell barcodes")

if STAGE in ("all", "import") and not FRAG_H5.exists():
    log("importing fragments (coordinate-sorted) with barcode whitelist")
    data = snap.pp.import_fragments(FRAG, chrom_sizes=snap.genome.hg38, file=str(FRAG_H5),
                                    whitelist=barcodes, min_num_fragments=0,
                                    sorted_by_barcode=False, n_jobs=8)
    log(f"  imported {data.shape}; n_fragment median {np.median(data.obs['n_fragment'][:]):.0f}")
    data.close()

if STAGE in ("all", "proximal"):
    log("proximal gene matrix (gene body + 2 kb upstream)")
    data = snap.read(str(FRAG_H5), backed="r")
    nfrag = pd.Series(np.asarray(data.obs["n_fragment"][:]), index=list(data.obs_names))
    gm = snap.pp.make_gene_matrix(data, gene_anno=GTF, upstream=2000, downstream=0,
                                  include_gene_body=True, file=str(OUT / "_gm_proximal.h5ad"))
    X = gm.X[:]
    a = ad.AnnData(X=sp.csr_matrix(X).astype(np.float32),
                   obs=pd.DataFrame(index=list(gm.obs_names)), var=pd.DataFrame(index=list(gm.var_names)))
    gm.close(); data.close()
    a.var_names_make_unique()
    a.obs["n_fragment"] = nfrag.reindex(a.obs_names).values
    a = a[[b for b in barcodes if b in a.obs_names]].copy()
    a.write_h5ad(OUT / "pbmc10k_geneactivity_proximal.h5ad")
    log(f"  saved proximal {a.shape}")

if STAGE in ("all", "enhancer"):
    log(f"tile matrix bin={BIN}")
    data = snap.read(str(FRAG_H5), backed="r+")
    if data.X is None or data.shape[1] == 0:
        snap.pp.add_tile_matrix(data, bin_size=BIN, exclude_chroms=["chrM", "chrY", "M", "Y"], n_jobs=8)
    T = sp.csr_matrix(data.X[:]); tile_names = list(data.var_names); tobs = list(data.obs_names)
    nfrag = pd.Series(np.asarray(data.obs["n_fragment"][:]), index=tobs)
    data.close()
    log(f"  tiles {T.shape}")
    tc = pd.DataFrame({"name": tile_names})
    tc["chrom"] = tc.name.str.split(":").str[0]
    se = tc.name.str.split(":").str[1].str.split("-")
    tc["start"] = se.str[0].astype(int); tc["end"] = se.str[1].astype(int)
    tc["mid"] = (tc.start + tc.end) // 2; tc["idx"] = np.arange(len(tc))
    rows = []
    with open(GTF) as f:
        for line in f:
            if line[0] == "#": continue
            p = line.split("\t")
            if p[2] != "gene" or 'gene_type "protein_coding"' not in p[8]: continue
            rows.append((p[0], int(p[3]), int(p[4]), p[6], p[8].split('gene_name "')[1].split('"')[0]))
    genes = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gene"]).drop_duplicates("gene")
    genes = genes[genes.chrom.isin(tc.chrom.unique())].reset_index(drop=True)
    log(f"  protein-coding genes {len(genes)}")
    by_chrom = {c: g.sort_values("mid") for c, g in tc.groupby("chrom")}
    wi, wj, wv = [], [], []
    for gi, g in enumerate(genes.itertuples(index=False)):
        ch = by_chrom.get(g.chrom)
        if ch is None: continue
        m = ch[(ch.mid >= g.start - EXTEND) & (ch.mid <= g.end + EXTEND)]
        if len(m) == 0: continue
        d = np.where(m.mid.values < g.start, g.start - m.mid.values,
                     np.where(m.mid.values > g.end, m.mid.values - g.end, 0))
        wi.extend(m.idx.values); wj.extend([gi] * len(m)); wv.extend(np.exp(-d / DECAY))
    W = sp.csr_matrix((wv, (wi, wj)), shape=(len(tc), len(genes)))
    G = (T @ W).tocsr().astype(np.float32)
    a = ad.AnnData(X=G, obs=pd.DataFrame(index=tobs), var=pd.DataFrame(index=genes.gene.values))
    a.var_names_make_unique()
    a.obs["n_fragment"] = nfrag.reindex(a.obs_names).values
    a = a[[b for b in barcodes if b in a.obs_names]].copy()
    a.write_h5ad(OUT / "pbmc10k_geneactivity_enhancer.h5ad")
    log(f"  saved enhancer {a.shape}")
log("done")
