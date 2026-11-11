"""
01: RNA-side analyses on the three matched RNA-ATAC datasets (PBMC, brain, bone marrow).
Nothing here modifies earlier results.

A. Sensitivity of the RNA reference (the CP10K RNA reference itself equalizes cell totals,
   so agreement with it alone cannot show that normalizing ATAC is right). The RNA reference is
   rebuilt four ways -- CP10K mean (as in the paper), raw mean (no normalization), pooled CPM and
   log-normalized mean -- and every ATAC strategy is scored against each:
     agreement of the per-gene maximum over RNA-informative genes (informative set defined by the
     CP10K reference, as in the paper, and also by each reference itself),
     and the L1 distance between the ATAC and RNA distributions of maxima over types.
B. The same aggregation question asked of RNA in the same nuclei (is the log-vs-linear
   difference specific to chromatin?). For each modality on the same cells and labels: per-cell
   detection rate, top share and coverage association under raw / log-normalized / CP10K / pooled
   aggregation, and the spread across types of the per-type mean of log(1+CP10K) over all genes
   (the residual factor of Eq. 1), with coverage = UMIs for RNA and fragments for ATAC.

Strategies and definitions are imported from tools/benchmark_aggregation.py and the loading
follows multiome/03, normal_tissues/05 and normal_tissues/06 exactly.
Out: rna_reference_sensitivity.csv, rna_vs_atac_aggregation.csv, run.log
"""
import os, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, hdf5plugin  # noqa
import scanpy as sc
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
ADD = HERE.parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
from benchmark_aggregation import m_raw, m_pseudobulk_cpm, m_depth_matched, m_cpm, m_lognorm, type_mean  # noqa

BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "/mnt/10t/scrna_atac"))
SEED = 42
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
def P_of(res): return res[0] if isinstance(res, tuple) else res


def informative(R):
    srt = np.sort(R, 0)
    return (srt[-1] >= 2 * np.maximum(srt[-2], 1e-9)) & (srt[-1] > 0)


def shares(P, cols, k):
    arg = P[:, cols].argmax(0)
    return np.bincount(arg, minlength=k) / arg.size


# ── datasets: (name, RNA counts cells x genes, rna genes, labels, list of ATAC (model, X, genes, coverage)) ──
def load_pbmc_like(name, lab_csv, h5, models):
    lab_df = pd.read_csv(lab_csv, index_col=0); lab_df = lab_df[lab_df.cell_type != "dropped"]
    rna = sc.read_10x_h5(h5, gex_only=True); rna.var_names_make_unique()
    atacs = []
    for model, path in models:
        a = ad.read_h5ad(path)
        cells = lab_df.index.intersection(a.obs_names)
        a = a[cells].copy()
        atacs.append((model, sp.csr_matrix(a.X).astype(np.float32), a.var_names.astype(str).values,
                      a.obs["n_fragment"].values.astype(float), list(cells)))
    cells0 = atacs[0][4]
    assert all(list(c[4]) == cells0 for c in atacs), "models differ in cells"
    rna = rna[cells0].copy()
    Xr = sp.csr_matrix(rna.X).astype(np.float32)
    lab = lab_df.loc[cells0, "cell_type"].values.astype(str)
    return name, Xr, rna.var_names.astype(str).values, lab, [(m, X, g, d) for m, X, g, d, _ in atacs]


def load_bmmc():
    a = ad.read_h5ad(str(BASE / "data/raw/multiome_bmmc_neurips/bmmc_multiome.h5ad"))
    lab_all = a.obs["cell_type"].astype(str).values
    vc = pd.Series(lab_all).value_counts(); types = sorted(vc.index[vc >= 100]); keep = np.isin(lab_all, types)
    dep = a.obs["ATAC_atac_fragments"].values.astype(float)[keep]
    D = sp.csr_matrix(a.obsm["ATAC_gene_activity"]).astype(np.float64)
    genes = np.asarray(a.uns["ATAC_gene_activity_var_names"]).astype(str)
    E = D.copy(); E.data = np.expm1(E.data); n = np.diff(E.indptr)
    unit = np.minimum.reduceat(E.data, E.indptr[:-1]); C = E.copy(); C.data = C.data / np.repeat(unit, n)
    assert np.abs(C.data - np.rint(C.data)).max() < 1e-3
    C.data = np.rint(C.data); X = C.astype(np.float32).tocsr()[keep]
    gex = (a.var["feature_types"] == "GEX").values
    Xr = sp.csr_matrix(a.layers["counts"][:, gex]).astype(np.float32)[keep]
    return "bone marrow", Xr, a.var_names[gex].astype(str).values, lab_all[keep], [("deposited", X, genes, dep)]


DATASETS = [
    lambda: load_pbmc_like("PBMC", ADD / "multiome" / "rna_labels.csv",
                           str(BASE / "data/raw/multiome_pbmc10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"),
                           [("count", str(BASE / "data/processed/multiome_pbmc10k/pbmc10k_geneactivity_proximal.h5ad")),
                            ("weighted", str(BASE / "data/processed/multiome_pbmc10k/pbmc10k_geneactivity_enhancer.h5ad"))]),
    lambda: load_pbmc_like("brain", ADD / "normal_tissues" / "brain3k_rna_labels.csv",
                           str(BASE / "data/raw/multiome_brain3k/human_brain_3k_filtered_feature_bc_matrix.h5"),
                           [("count", str(BASE / "data/processed/multiome_brain3k/brain3k_geneactivity_proximal.h5ad")),
                            ("weighted", str(BASE / "data/processed/multiome_brain3k/brain3k_geneactivity_enhancer.h5ad"))]),
    load_bmmc,
]

rng = np.random.default_rng(SEED)
sens_rows, cmp_rows = [], []
for loader in DATASETS:
    name, Xr, rgenes, lab, atacs = loader()
    types = sorted(set(lab)); k = len(types)
    umi = np.asarray(Xr.sum(1)).ravel()
    log(f"== {name}: {Xr.shape[0]} cells, {k} types, RNA genes {Xr.shape[1]}")
    rna_det = np.asarray((Xr > 0).mean(0)).ravel()
    RNA_REFS = {"cpm_mean": P_of(m_cpm(Xr, lab, types)), "raw_mean": P_of(m_raw(Xr, lab, types)),
                "pseudobulk_cpm": P_of(m_pseudobulk_cpm(Xr, lab, types)), "lognorm_mean": P_of(m_lognorm(Xr, lab, types))}
    inf_cpm = informative(RNA_REFS["cpm_mean"])
    # ── B (RNA side): aggregation of RNA itself ──
    rna_expr = np.asarray(Xr.sum(0)).ravel() > 0
    med_umi = np.array([np.median(umi[lab == t]) for t in types])
    nnz_r = np.diff(Xr.indptr) / Xr.shape[1]
    det_r = np.array([nnz_r[lab == t].mean() for t in types])
    for sname, P in RNA_REFS.items():
        sh = shares(P, rna_expr, k)
        cmp_rows.append(dict(dataset=name, modality="RNA", model="", strategy=sname, n_genes=int(rna_expr.sum()),
                             coverage_spread=float(med_umi.max() / med_umi.min()),
                             detection_min=float(det_r.min()), detection_max=float(det_r.max()),
                             top_share=float(sh.max()), top_type=types[int(sh.argmax())],
                             rho_coverage_share=float(spearmanr(med_umi, sh)[0])))
    Xn = Xr.multiply(1e4 / np.maximum(umi, 1)[:, None]).tocsr(); Xn.data = np.log1p(Xn.data)
    lm = np.array([np.asarray(Xn[lab == t][:, rna_expr].mean(1)).ravel().mean() for t in types])
    cmp_rows.append(dict(dataset=name, modality="RNA", model="", strategy="logmean_spread", top_share=float(lm.max() / lm.min()),
                         rho_coverage_share=float(spearmanr(det_r, lm)[0])))
    for model, X, agenes, dep in atacs:
        expressed = np.asarray(X.sum(0)).ravel() > 0
        med = np.array([np.median(dep[lab == t]) for t in types])
        nnz_a = np.diff(X.indptr) / X.shape[1]
        det_a = np.array([nnz_a[lab == t].mean() for t in types])
        rmap = {g: i for i, g in enumerate(rgenes)}
        pairs = [(i, rmap[g]) for i, g in enumerate(agenes) if g in rmap and expressed[i] and rna_det[rmap[g]] >= 0.01]
        ai = np.array([p[0] for p in pairs]); ri = np.array([p[1] for p in pairs])
        ATAC = {"raw_mean": P_of(m_raw(X, lab, types)), "lognorm_mean": P_of(m_lognorm(X, lab, types)),
                "pseudobulk_cpm": P_of(m_pseudobulk_cpm(X, lab, types)), "cpm_mean": P_of(m_cpm(X, lab, types))}
        if model != "weighted":      # thinning is defined for integer counts only
            ATAC["depth_matched"] = P_of(m_depth_matched(X, lab, types, rng))
        # ── B (ATAC side), same cells and labels ──
        for sname in ("raw_mean", "lognorm_mean", "cpm_mean", "pseudobulk_cpm"):
            sh = shares(ATAC[sname], expressed, k)
            cmp_rows.append(dict(dataset=name, modality="ATAC", model=model, strategy=sname, n_genes=int(expressed.sum()),
                                 coverage_spread=float(med.max() / med.min()),
                                 detection_min=float(det_a.min()), detection_max=float(det_a.max()),
                                 top_share=float(sh.max()), top_type=types[int(sh.argmax())],
                                 rho_coverage_share=float(spearmanr(med, sh)[0])))
        rs = np.asarray(X.sum(1)).ravel()
        Xa = X.multiply(1e4 / np.maximum(rs, 1)[:, None]).tocsr(); Xa.data = np.log1p(Xa.data)
        lma = np.array([np.asarray(Xa[lab == t][:, expressed].mean(1)).ravel().mean() for t in types])
        cmp_rows.append(dict(dataset=name, modality="ATAC", model=model, strategy="logmean_spread",
                             top_share=float(lma.max() / lma.min()), rho_coverage_share=float(spearmanr(det_a, lma)[0])))
        # ── A: every ATAC strategy against every RNA reference ──
        for rname, R in RNA_REFS.items():
            Rs = R[:, ri]; r_arg = Rs.argmax(0); inf_fixed = inf_cpm[ri]; inf_own = informative(R)[ri]
            rshare = np.bincount(r_arg, minlength=k) / r_arg.size
            for aname, P in ATAC.items():
                a_arg = P[:, ai].argmax(0)
                ashare = np.bincount(a_arg, minlength=k) / a_arg.size
                sens_rows.append(dict(dataset=name, model=model, rna_reference=rname, atac_strategy=aname,
                                      n_shared=len(ai), n_inf_cpm=int(inf_fixed.sum()), n_inf_own=int(inf_own.sum()),
                                      agree_inf_cpmset=float((a_arg[inf_fixed] == r_arg[inf_fixed]).mean()),
                                      agree_inf_ownset=float((a_arg[inf_own] == r_arg[inf_own]).mean()),
                                      agree_all=float((a_arg == r_arg).mean()),
                                      l1_to_rna=float(np.abs(ashare - rshare).sum()),
                                      rna_top_share=float(rshare.max()), rna_top_type=types[int(rshare.argmax())]))
        log(f"  {model}: done ({len(ai)} shared genes)")
    del Xr, atacs
    pd.DataFrame(sens_rows).to_csv(HERE / "rna_reference_sensitivity.csv", index=False)
    pd.DataFrame(cmp_rows).to_csv(HERE / "rna_vs_atac_aggregation.csv", index=False)
log("done")
