"""
05: Is the residual argmax concentration in pDC a depth effect or a small-n effect?
    (a) RNA truth itself: share of argmaxima per type in the RNA CPM-mean profile (same genes).
    (b) Equal-n control: subsample every RNA-defined type to n = min type size (112 cells),
        20 draws, recompute top share / rho(depth, share) / RNA-argmax agreement per method.
Out: multiome_equal_n.csv, multiome_rna_truth_share.csv
"""
import sys, time, warnings
import os
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, hdf5plugin  # noqa
import scanpy as sc
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from benchmark_aggregation import m_raw, m_zscale, m_pseudobulk_cpm, m_depth_matched, m_cpm, m_lognorm, m_archr  # noqa

HERE = Path(__file__).resolve().parent
BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
PROC = BASE / "data/processed/multiome_pbmc10k"
H5 = BASE / "data/raw/multiome_pbmc10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
N_DRAW = 20
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

lab_df = pd.read_csv(HERE / "rna_labels.csv", index_col=0); lab_df = lab_df[lab_df.cell_type != "dropped"]
types = sorted(lab_df.cell_type.unique())
rna = sc.read_10x_h5(H5, gex_only=True); rna.var_names_make_unique(); rna = rna[lab_df.index].copy()
Xr = sp.csr_matrix(rna.X).astype(np.float32); lab_r = lab_df.cell_type.values.astype(str)
R = m_cpm(Xr, lab_r, types); rna_det = np.asarray((Xr > 0).mean(0)).ravel()
rna_genes = rna.var_names.astype(str).values; rmap = {g: i for i, g in enumerate(rna_genes)}
srt = np.sort(R, 0); rna_inf = (srt[-1] >= 2 * np.maximum(srt[-2], 1e-9)) & (srt[-1] > 0)
rng = np.random.default_rng(0)
rows, trows = [], []
for model in ["proximal", "enhancer"]:
    a = ad.read_h5ad(PROC / f"pbmc10k_geneactivity_{model}.h5ad"); a = a[lab_df.index.intersection(a.obs_names)].copy()
    lab = lab_df.loc[a.obs_names, "cell_type"].values.astype(str); dep = a.obs["n_fragment"].values.astype(float)
    X = sp.csr_matrix(a.X).astype(np.float32); genes = a.var_names.astype(str).values
    expressed = np.asarray(X.sum(0)).ravel() > 0
    ai = np.array([i for i, g in enumerate(genes) if g in rmap and expressed[i] and rna_det[rmap[g]] >= 0.01])
    ri = np.array([rmap[genes[i]] for i in ai]); r_arg = R[:, ri].argmax(0); inf = rna_inf[ri]
    # (a) RNA truth concentration on the same gene sets
    sh_all = np.bincount(r_arg, minlength=len(types)) / len(r_arg)
    sh_inf = np.bincount(r_arg[inf], minlength=len(types)) / inf.sum()
    med_depth = np.array([np.median(dep[lab == t]) for t in types])
    for t, s1, s2, md in zip(types, sh_all, sh_inf, med_depth):
        trows.append(dict(model=model, cell_type=t, n_cells=int((lab == t).sum()), median_true_depth=md,
                          rna_argmax_share_all=s1, rna_argmax_share_informative=s2))
    log(f"{model}: RNA-truth argmax share (all/informative): " + ", ".join(f"{t} {s1:.2f}/{s2:.2f}" for t, s1, s2 in zip(types, sh_all, sh_inf)))
    # (b) equal-n subsampling
    n_min = min((lab == t).sum() for t in types)
    methods = [("raw_mean", lambda M, l: m_raw(M, l, types)), ("zscale_mean", lambda M, l: m_zscale(M, l, types)),
               ("lognorm_mean", lambda M, l: m_lognorm(M, l, types)), ("pseudobulk_cpm", lambda M, l: m_pseudobulk_cpm(M, l, types)),
               ("depth_matched", lambda M, l: m_depth_matched(M, l, types, rng)), ("cpm_mean", lambda M, l: m_cpm(M, l, types)),
               ("archr_style", lambda M, l: m_archr(M, l, types))]
    for d in range(N_DRAW):
        idx = np.concatenate([rng.choice(np.where(lab == t)[0], n_min, replace=False) for t in types])
        Xs, ls, ds = X[idx], lab[idx], dep[idx]
        md = np.array([np.median(ds[ls == t]) for t in types])
        for name, fn in methods:
            res = fn(Xs, ls); P = res[0] if isinstance(res, tuple) else res
            arg = P[:, expressed].argmax(0); share = np.bincount(arg, minlength=len(types)) / arg.size
            top = int(share.argmax()); rho = spearmanr(md, share)[0]
            a_arg = P[:, ai].argmax(0)
            rows.append(dict(model=model, draw=d, method=name, n_per_type=int(n_min), top_share=float(share[top]),
                             top_type=types[top], deepest_type=types[int(md.argmax())], rho_depth_share=float(rho),
                             share_pDC=float(share[types.index("pDC")]),
                             argmax_agree_all=float((a_arg == r_arg).mean()), argmax_agree_informative=float((a_arg[inf] == r_arg[inf]).mean())))
        if d % 5 == 0: log(f"  {model} draw {d} done")
E = pd.DataFrame(rows); E.to_csv(HERE / "multiome_equal_n.csv", index=False)
pd.DataFrame(trows).to_csv(HERE / "multiome_rna_truth_share.csv", index=False)
print(E.groupby(["model", "method"]).agg(top_share=("top_share", "median"), share_pDC=("share_pDC", "median"),
      rho=("rho_depth_share", "median"), top_is_deepest=("top_type", lambda s: (s == E.loc[s.index, "deepest_type"]).mean()),
      agree_inf=("argmax_agree_informative", "median")).round(3).to_string())
log("done")
