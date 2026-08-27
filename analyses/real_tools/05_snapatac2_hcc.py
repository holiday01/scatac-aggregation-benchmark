"""SnapATAC2 v2.9: snap.pp.make_gene_matrix with ALL defaults on the HCC fragment import
(upstream=2000, downstream=0, include_gene_body=True, id_type='gene', counting_strategy='paired-insertion'),
GENCODE v38 GTF as in scripts/analysis/90. Then the aggregations an analyst would apply to the tool's own output:
raw mean, mean of log1p(normalize_total) [the SnapATAC2 tutorial's sc.pp.normalize_total() + sc.pp.log1p(), i.e.
target_sum=None -> median total], the same with target_sum=1e4, and the CPM(1e4) mean. Also checks whether the
default output equals the paper's 'proximal' matrix."""
import time, sys, numpy as np, pandas as pd, hdf5plugin, anndata as ad, scipy.sparse as sp, scanpy as sc, snapatac2 as snap
import os
from pathlib import Path
from score_profiles import SCRATCH, OUTDIR, BASE, relabel, score, true_depth
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
GTF = Path(os.environ.get("GENCODE_GTF", "data/reference/gencode.v38.annotation.gtf"))
FRAG_H5 = BASE / "data/processed/integrated/_s90_work/hcc_fragments.h5ad"
RES = OUTDIR / "results"; RES.mkdir(exist_ok=True)

data = snap.read(str(FRAG_H5), backed="r")
log(f"fragment AnnData: {data.shape}; obs n_fragment median {np.median(data.obs['n_fragment'][:])}")
t0 = time.time()
gm = snap.pp.make_gene_matrix(data, gene_anno=GTF)       # all defaults
data.close()
log(f"make_gene_matrix defaults: {gm.shape}, {time.time()-t0:.0f}s; X dtype {gm.X.dtype}")
gm = ad.AnnData(X=sp.csr_matrix(gm.X), obs=pd.DataFrame(index=list(gm.obs_names)), var=pd.DataFrame(index=list(gm.var_names)))
ndup = gm.var_names.duplicated().sum(); gm.var_names_make_unique()
log(f"duplicated gene names in tool output: {ndup}")

# compare with the paper's proximal matrix
ref = ad.read_h5ad(BASE / "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad")
gm = gm[ref.obs_names].copy()
same_var = list(gm.var_names) == list(ref.var_names)
Xg = sp.csr_matrix(gm.X).astype(np.int64); Xr = sp.csr_matrix(ref.X).astype(np.int64)
if same_var:
    diff = (Xg - Xr); ndiff = diff.nnz
    log(f"var names identical: True; entries differing from paper's proximal matrix: {ndiff} of {Xg.nnz} nonzeros (sum|diff|={abs(diff).sum()})")
else:
    common = ref.var_names.intersection(gm.var_names); log(f"var names differ; common {len(common)}")
    Xg2 = sp.csr_matrix(gm[:, common].X).astype(np.int64); Xr2 = sp.csr_matrix(ref[:, common].X).astype(np.int64)
    ndiff = (Xg2 - Xr2).nnz; log(f"entries differing on common genes: {ndiff}")
lab = relabel(ref.obs["cell_type"].values); types = sorted(pd.unique(lab))
genes = np.asarray(gm.var_names); X = sp.csr_matrix(gm.X).astype(np.float32)
expressed = np.asarray(X.sum(0)).ravel() > 0
tot = np.asarray(X.sum(1)).ravel(); log(f"median matrix row sum {np.median(tot):.0f}")

def tmean(M): return np.vstack([np.asarray(M[lab == t].mean(0)).ravel() for t in types])
rows, mrows = [], []
def add(name, what, P):
    r, mt = score(P, types, genes, expressed); r = dict(tool="SnapATAC2", path=name, computes=what, **r); rows.append(r)
    mt.insert(0, "path", name); mrows.append(mt)
    log(f"{name:40} top {r['top_type']:>20} {r['top_share']:5.1%}  markers {r['markers_correct']}/{r['markers_total']} perm p={r['perm_p']:.4f}")
add("raw_mean", "mean of make_gene_matrix counts per type", tmean(X))
a = ad.AnnData(X=X.copy()); sc.pp.normalize_total(a); sc.pp.log1p(a)         # tutorial call: target_sum=None -> median
add("tutorial_normalize_total_median_log1p_mean", "sc.pp.normalize_total() [target_sum=None=median]; sc.pp.log1p; mean per type", tmean(a.X))
a = ad.AnnData(X=X.copy()); sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
add("normalize_total_1e4_log1p_mean", "sc.pp.normalize_total(target_sum=1e4); sc.pp.log1p; mean per type", tmean(a.X))
a = ad.AnnData(X=X.copy()); sc.pp.normalize_total(a, target_sum=1e4)
add("cpm_mean_1e4", "sc.pp.normalize_total(target_sum=1e4) then mean per type (no log)", tmean(a.X))
a = ad.AnnData(X=X.copy()); sc.pp.normalize_total(a)
add("cpm_mean_median", "sc.pp.normalize_total() [median target] then mean per type (no log)", tmean(a.X))
S = np.vstack([np.asarray(X[lab == t].sum(0)).ravel() for t in types]).astype(np.float64)
add("pseudobulk_cpm", "sum per type then CPM", S / S.sum(1, keepdims=True) * 1e6)
R = pd.DataFrame(rows); R["identical_to_paper_proximal"] = bool(same_var and ndiff == 0); R["n_entries_differing"] = ndiff
R.to_csv(RES / "snapatac2_results.csv", index=False); pd.concat(mrows).to_csv(RES / "snapatac2_markers.csv", index=False)
log("wrote snapatac2_results.csv")
