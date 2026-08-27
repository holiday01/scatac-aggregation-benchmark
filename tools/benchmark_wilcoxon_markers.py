"""
The same question asked of a real tool: scanpy's rank_genes_groups (Wilcoxon,
the FindAllMarkers-equivalent that Seurat, Signac and scanpy users run to get
per-cell-type markers) on the HCC proximal gene-activity matrix, run on the raw
counts and on the per-cell normalised matrix.

Read-outs per run:
  - of the top-50 markers called for each cell type, how many are the type's own
    canonical markers, and how many of the canonical panel land in the wrong type;
  - how the number of significant "up" genes per type (padj<0.05, log2FC>0.5)
    tracks each type's median true depth.

Run:  python3 tools/benchmark_wilcoxon_markers.py
Out:  results/benchmark_wilcoxon.csv
"""
import os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, h5py, hdf5plugin  # noqa
import scipy.sparse as sp, anndata as ad, scanpy as sc
from scipy.stats import spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_aggregation import PANEL, true_depth, BASE, OUT

PATH = BASE / "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad"
TOPN = 50

def run(a, tag):
    sc.tl.rank_genes_groups(a, "cell_type", method="wilcoxon", n_genes=a.n_vars, use_raw=False)
    types = list(a.obs.cell_type.cat.categories)
    rows = []
    for t in types:
        df = sc.get.rank_genes_groups_df(a, group=t)
        top = set(df.head(TOPN).names)
        n_sig = int(((df.pvals_adj < 0.05) & (df.logfoldchanges > 0.5)).sum())
        own = [g for g, tg in PANEL["HCC"].items() if tg == t]
        other = [g for g, tg in PANEL["HCC"].items() if tg != t]
        rows.append(dict(norm=tag, cell_type=t, n_sig_up=n_sig,
                         own_markers_in_top50=sum(g in top for g in own), own_markers=len(own),
                         foreign_markers_in_top50=sum(g in top for g in other),
                         top10=",".join(df.head(10).names)))
    return pd.DataFrame(rows)

def main():
    a = ad.read_h5ad(PATH); a.obs["cell_type"] = a.obs.cell_type.astype("category")
    dep = true_depth("HCC", a.obs_names.values); a.obs["depth"] = dep
    med = a.obs.groupby("cell_type", observed=True).depth.median()
    print("median true depth per type:\n" + med.round(0).to_string())
    a.X = a.X.astype(np.float32)
    out = []
    raw = a.copy(); out.append(run(raw, "raw_counts"))
    nrm = a.copy(); sc.pp.normalize_total(nrm, target_sum=1e4); sc.pp.log1p(nrm); out.append(run(nrm, "normalize_total_log1p"))
    R = pd.concat(out); R["median_depth"] = R.cell_type.map(med)
    for tag, s in R.groupby("norm"):
        rho, p = spearmanr(s.median_depth, s.n_sig_up)
        print(f"\n[{tag}]  Spearman(depth, n significant up-genes) = {rho:+.2f} (p={p:.3f})")
        print(s[["cell_type", "median_depth", "n_sig_up", "own_markers_in_top50", "own_markers", "foreign_markers_in_top50", "top10"]].to_string(index=False))
    R.to_csv(OUT / "benchmark_wilcoxon.csv", index=False)
    print(f"\nwrote {OUT/'benchmark_wilcoxon.csv'}")

if __name__ == "__main__":
    main()
