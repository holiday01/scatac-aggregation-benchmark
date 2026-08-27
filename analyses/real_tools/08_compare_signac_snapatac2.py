"""How different is Signac GeneActivity() (gene body + 2 kb upstream, fragments overlapping, protein-coding <=500 kb,
GENCODE v38 gene spans) from SnapATAC2 make_gene_matrix() defaults (same window, paired-insertion counting, all
GENCODE genes) on the same 12,029 cells?"""
import numpy as np, pandas as pd, hdf5plugin, anndata as ad, scipy.sparse as sp, scipy.io
from scipy.stats import spearmanr, pearsonr
from score_profiles import SCRATCH, OUTDIR, BASE, PANEL, relabel, true_depth
RES = OUTDIR / "results"
ref = ad.read_h5ad(BASE / "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad")
M = scipy.io.mmread(str(SCRATCH / "signac_ga/matrix.mtx")).tocsr()            # genes x cells
feats = np.array(open(SCRATCH / "signac_ga/features.txt").read().split("\n")[:M.shape[0]])
cells = np.array(open(SCRATCH / "signac_ga/cells.txt").read().split("\n")[:M.shape[1]])
sig = ad.AnnData(X=M.T.tocsr().astype(np.float32), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=feats))
ndup = int(sig.var_names.duplicated().sum()); sig.var_names_make_unique(); sig = sig[ref.obs_names].copy()
print(f'duplicated gene names in Signac GeneActivity output: {ndup}')
common = ref.var_names.intersection(sig.var_names)
A = sp.csr_matrix(ref[:, common].X).astype(np.float64); B = sp.csr_matrix(sig[:, common].X).astype(np.float64)
ta, tb = np.asarray(A.sum(1)).ravel(), np.asarray(B.sum(1)).ravel()
ga, gb = np.asarray(A.sum(0)).ravel(), np.asarray(B.sum(0)).ravel()
dep = true_depth("HCC", ref.obs_names.values)
print(f"Signac genes {sig.n_vars}, SnapATAC2 genes {ref.n_vars}, common names {len(common)}")
print(f"per-cell totals on common genes: Signac/SnapATAC2 median ratio {np.median(tb/ta):.3f}; Spearman {spearmanr(ta, tb)[0]:.4f}")
print(f"per-cell totals vs true fragment depth: SnapATAC2 rho {spearmanr(ta, dep)[0]:.3f}, Signac rho {spearmanr(tb, dep)[0]:.3f}")
print(f"per-gene totals: Pearson(log1p) {pearsonr(np.log1p(ga), np.log1p(gb))[0]:.4f}, Spearman {spearmanr(ga, gb)[0]:.4f}")
print(f"identical entries fraction (common block): {(A != B).nnz / max(A.nnz, 1):.3f} of SnapATAC2 nonzeros differ")
panel = [g for g in PANEL['HCC'] if g in set(common)]
print("22-marker panel present in Signac output:", len(panel), "/ 22; missing:", [g for g in PANEL['HCC'] if g not in set(sig.var_names)])
lab = relabel(ref.obs['cell_type'].values); types = sorted(pd.unique(lab))
rows = []
for g in panel:
    i = list(common).index(g)
    rows.append(dict(gene=g, snapatac2_total=int(ga[i]), signac_total=int(gb[i]), ratio=gb[i] / max(ga[i], 1)))
print(pd.DataFrame(rows).to_string(index=False))
pd.DataFrame(dict(metric=["n_genes_signac", "n_genes_snapatac2", "n_common", "median_percell_ratio_signac_over_snapatac2",
                          "spearman_percell_totals", "spearman_pergene_totals", "frac_entries_differ_common_block"],
                  value=[sig.n_vars, ref.n_vars, len(common), np.median(tb/ta), spearmanr(ta, tb)[0], spearmanr(ga, gb)[0], (A != B).nnz / A.nnz])
             ).to_csv(RES / "signac_vs_snapatac2_matrix.csv", index=False)
