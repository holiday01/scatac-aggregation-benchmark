"""Analyses required by the pre-submission audit (2026-08-27).

R1  integrality of every matrix, and the depth-matched control restricted to
    integer count matrices (the enhancer-aware matrices hold fractional
    distance-weighted scores, on which binomial thinning is not defined).
R2  donor-aware sensitivity: leave-one-sample-out recomputation of the top
    argmax share, so that uncertainty reflects between-sample variation rather
    than cell-level resampling.
R3  HCC marker concordance under the ORIGINAL motif-derived labels (which are
    independent of every gene-activity matrix) as well as under the two
    corrected labels, to expose the circularity of the correction.
R4  the realised depth spread of the "equal-depth" baseline on the mean as well
    as the median.

Out: revision/*.csv
"""
import sys, time
import os
from pathlib import Path
import numpy as np, pandas as pd, h5py, hdf5plugin  # noqa
import scipy.sparse as sp, anndata as ad
from scipy.stats import spearmanr

BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
OUT = Path(__file__).resolve().parents[1] / "results/revision"; OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_aggregation import PANEL, true_depth, perm_p  # noqa

MAT = [("HCC","proximal","data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad","sample_id"),
       ("HCC","enhancer","data/processed/integrated/atac_hcc_geneactivity_enhancer.h5ad",None),
       ("COAD","proximal","data/processed/COAD/coad_geneactivity_proximal.h5ad","sample"),
       ("COAD","enhancer","data/processed/COAD/coad_geneactivity_enhancer.h5ad",None),
       ("LUAD","proximal","data/processed/LUAD/luad_geneactivity_proximal.h5ad","sample"),
       ("LUAD","enhancer","data/processed/LUAD/luad_geneactivity_enhancer.h5ad",None)]
RELABEL = {"B_cell":"Endothelial_stromal","DC":"B_cell"}

def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
def scale(X, t=1e4):
    d = np.asarray(X.sum(1)).ravel().astype(np.float64); d[d==0]=1
    return (sp.diags((t/d).astype(np.float32)) @ X).tocsr()
def tmean(X, lab, types):
    return np.vstack([np.asarray(X[lab==t].mean(0)).ravel() for t in types])
def share(P, expressed, types):
    arg = P[:, expressed].argmax(0)
    s = np.bincount(arg, minlength=len(types))/arg.size
    return float(s.max()), types[int(s.argmax())], s

def methods(X, lab, types, integer, rng):
    Xc = scale(X); Xl = Xc.copy(); Xl.data = np.log1p(Xl.data)
    S = np.vstack([np.asarray(X[lab==t].sum(0)).ravel() for t in types]).astype(np.float64)
    out = {"raw_mean": tmean(X,lab,types), "cpm_mean": tmean(Xc,lab,types),
           "lognorm_mean": tmean(Xl,lab,types), "pseudobulk_cpm": S/S.sum(1,keepdims=True)*1e6}
    if integer:
        tot = np.asarray(X.sum(1)).ravel().astype(np.float64)
        target = min(np.median(tot[lab==t]) for t in types)
        p = np.minimum(1.0, target/np.maximum(tot,1.0))
        Xd = X.tocsr(copy=True); Xd.data = rng.binomial(np.rint(Xd.data).astype(np.int64),
                                                        np.repeat(p, np.diff(Xd.indptr))).astype(np.float32)
        Xd.eliminate_zeros(); out["depth_matched"] = tmean(Xd,lab,types)
    return out

rows_int, rows_loso, rows_rel = [], [], []
for coh, model, rel, samplecol in MAT:
    log(f"{coh}/{model}")
    a = ad.read_h5ad(BASE/rel)
    lab = a.obs["cell_type"].astype(str).values
    if coh == "HCC": lab_rel = np.array([RELABEL.get(x,x) for x in lab])
    else: lab_rel = lab
    X = (a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)).tocsr().astype(np.float32)
    d = np.asarray(X.data, dtype=np.float64)
    integer = bool(np.allclose(d, np.rint(d)))
    rows_int.append(dict(cohort=coh, model=model, integer=integer, min_value=float(d.min()),
                         frac_below_0_5=float(np.mean(d<0.5)), frac_below_1=float(np.mean(d<1)),
                         n_nonzero=int(d.size)))
    log(f"  integer={integer} frac<0.5={np.mean(d<0.5):.3f}")
    types = sorted(pd.unique(lab_rel)); expressed = np.asarray(X.sum(0)).ravel() > 0
    rng = np.random.default_rng(42)
    P = methods(X, lab_rel, types, integer, rng)
    dep = true_depth(coh, a.obs_names.values)
    med_depth = np.array([np.nanmedian(dep[lab_rel==t]) for t in types])
    for m, prof in P.items():
        ts, tt, s = share(prof, expressed, types)
        rho, pv = spearmanr(med_depth, s)
        rows_int[-1][f"{m}_top_share"] = ts
        rows_int[-1][f"{m}_top_type"] = tt

    # R3: HCC, marker concordance under original vs corrected labels
    if coh == "HCC":
        gidx = {g:i for i,g in enumerate(a.var_names.astype(str))}
        for tag, L in [("motif_original", lab), ("corrected", lab_rel)]:
            T = sorted(pd.unique(L))
            Pm = methods(X, L, T, integer, np.random.default_rng(42))
            panel = {g:t for g,t in PANEL["HCC"].items() if g in gidx and t in T}
            for m, prof in Pm.items():
                marg = np.array([T[int(prof[:, gidx[g]].argmax())] for g in panel])
                mtgt = np.array(list(panel.values())); obs = int((marg==mtgt).sum())
                pp,_,_ = perm_p(marg, mtgt, obs, np.random.default_rng(7))
                rows_rel.append(dict(model=model, labels=tag, method=m, n_types=len(T),
                                     markers_correct=obs, markers_total=len(panel), perm_p=pp))
                log(f"   {tag:15} {m:15} {obs}/{len(panel)} p={pp:.4f}")

    # R2: leave-one-sample-out on the proximal (integer) matrices
    if samplecol is not None:
        smp = a.obs[samplecol].astype(str).values
        for s_out in sorted(pd.unique(smp)):
            k = smp != s_out
            Xk, Lk = X[k], lab_rel[k]
            Tk = [t for t in types if (Lk==t).sum() >= 20]
            if len(Tk) < len(types): log(f"    LOSO {s_out}: {len(types)-len(Tk)} type(s) dropped (<20 cells)")
            ek = np.asarray(Xk.sum(0)).ravel() > 0
            Pk = methods(Xk, Lk, Tk, integer, np.random.default_rng(42))
            for m, prof in Pk.items():
                ts, tt, _ = share(prof, ek, Tk)
                rows_loso.append(dict(cohort=coh, model=model, sample_left_out=s_out,
                                      n_cells=int(k.sum()), n_types=len(Tk), method=m,
                                      top_share=ts, top_type=tt))
        log(f"  LOSO done over {len(pd.unique(smp))} samples")
    del a, X

pd.DataFrame(rows_int).to_csv(OUT/"matrix_integrality_and_shares.csv", index=False)
pd.DataFrame(rows_loso).to_csv(OUT/"loso_top_share.csv", index=False)
pd.DataFrame(rows_rel).to_csv(OUT/"hcc_label_sensitivity.csv", index=False)
print("\n=== integrality ===")
print(pd.DataFrame(rows_int)[["cohort","model","integer","min_value","frac_below_0_5"]].to_string(index=False))
L = pd.DataFrame(rows_loso)
print("\n=== leave-one-sample-out top share (min-max over samples) ===")
g = L.groupby(["cohort","model","method"]).top_share.agg(["min","max","median","count"])
print((g*1).round(3).to_string())
print("\n=== HCC label sensitivity ===")
print(pd.DataFrame(rows_rel).to_string(index=False))
