"""
Downstream consequence of the aggregation choice: NNLS deconvolution of pseudo-bulk
mixtures against cell-type reference profiles built from the SAME scATAC gene-activity
matrix by each aggregation method of tools/benchmark_aggregation.py.

Design (per cohort: HCC, COAD, LUAD proximal matrices; HCC relabel B_cell->Endothelial_stromal,
DC->B_cell applied after loading, as in the benchmark):
  1. cells split 50/50 stratified by type, seeds 0,1,2;
  2. training half -> reference (genes x types) by raw mean, z-scale mean, log1p(CPM) mean,
     pseudobulk CPM, depth-matched raw mean, CPM mean (functions imported from the benchmark
     so the definitions are identical); plus expm1(log1p-CPM mean) as a variant of the log ref;
  3. held-out half -> 100 mixtures per seed: Dirichlet(alpha=1) proportions over types,
     500 cells drawn multinomially, mixture = SUM OF RAW COUNTS of the sampled cells,
     then CPM on the linear scale (primary, "rawsum").  Secondary construction "cellnorm":
     each sampled cell scaled to 1e4 before summation (every cell contributes equally).
  4. NNLS (scipy.optimize.nnls) of each mixture on each reference, fractions normalised to 1.
Gene sets: "common" = union of the top-200 genes per type by specificity ratio
     (type mean / max other-type mean, with a small pseudocount) computed ONCE per split from
     the CPM-mean reference and used for every method (only the profile VALUES differ);
     "own" = each method selects its own top-200 by its own ratio (z-scale: difference,
     because z-scaled means can be negative).
Two truths are recorded for every mixture: the realised CELL fraction (n_t / 500) and the
FRAGMENT share (fraction of the mixture's total counts contributed by type t).  A reference
whose columns are per-cell-normalised models fragment share; a raw-mean reference models
cell fraction (see REPORT.md).

Run:  nohup python3 deconvolution_benchmark.py > run.log 2>&1 &
Out:  estimates.csv, metrics_by_seed.csv, metrics_summary.csv, bias_by_type.csv,
      signature_sizes.csv, common_signature_genes.csv, split_depth.csv
"""
import os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, h5py, hdf5plugin  # noqa: F401
import scipy.sparse as sp
import anndata as ad
from scipy.optimize import nnls

BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "/mnt/10t/scrna_atac"))
HERE = Path(__file__).resolve().parent
APBC = HERE.parent.parent.parent                      # .../apbc2026
sys.path.insert(0, str(APBC / "tools"))
from benchmark_aggregation import (m_raw, m_zscale, m_pseudobulk_cpm, m_depth_matched,  # noqa: E402
                                   m_cpm, m_lognorm, percell_scale)

SEEDS = [0, 1, 2]
N_MIX, N_CELLS, TOP_N, ALPHA = 100, 500, 200, 1.0
MATRICES = [("HCC",  "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad"),
            ("COAD", "data/processed/COAD/coad_geneactivity_proximal.h5ad"),
            ("LUAD", "data/processed/LUAD/luad_geneactivity_proximal.h5ad")]
METHODS = ["raw_mean", "zscale_mean", "lognorm_mean", "lognorm_mean_expm1",
           "pseudobulk_cpm", "depth_matched", "cpm_mean"]


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load(cohort, rel):
    a = ad.read_h5ad(BASE / rel)
    lab = a.obs["cell_type"].astype(str).values
    keep = ~pd.isna(a.obs["cell_type"].values) & (lab != "nan")
    a = a[keep].copy(); lab = a.obs["cell_type"].astype(str).values
    if cohort == "HCC":   # same relabel as tools/benchmark_aggregation.py (tools/check_hcc_labels.py)
        lab = np.where(lab == "B_cell", "Endothelial_stromal", np.where(lab == "DC", "B_cell", lab))
    X = (a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)).tocsr().astype(np.float32)
    genes = a.var_names.astype(str).values
    obs = a.obs_names.values
    del a
    return X, lab, genes, obs


def build_refs(Xtr, labtr, types, rng):
    """types x genes for each method (benchmark definitions)."""
    R = {}
    R["raw_mean"] = m_raw(Xtr, labtr, types)
    R["zscale_mean"] = m_zscale(Xtr, labtr, types)
    R["lognorm_mean"] = m_lognorm(Xtr, labtr, types)
    R["lognorm_mean_expm1"] = np.expm1(R["lognorm_mean"])
    R["pseudobulk_cpm"] = m_pseudobulk_cpm(Xtr, labtr, types)
    R["depth_matched"], target = m_depth_matched(Xtr, labtr, types, rng)
    R["cpm_mean"] = m_cpm(Xtr, labtr, types)
    return {k: np.asarray(v, dtype=np.float64) for k, v in R.items()}, target


def top_genes(P, expressed, mode="ratio", top_n=TOP_N):
    """Union over types of the top_n genes by specificity.  ratio: (m_t + c)/(max_other + c),
    c = 1% of the mean positive entry of P, and m_t must be > 0.  diff: m_t - max_other."""
    T, G = P.shape
    sel, per_type = set(), {}
    if mode == "ratio":
        c = 0.01 * P[P > 0].mean() if (P > 0).any() else 1e-9
    for i in range(T):
        other = np.max(np.delete(P, i, axis=0), axis=0)
        if mode == "ratio":
            score = (P[i] + c) / (other + c)
            score[P[i] <= 0] = -np.inf
        else:
            score = P[i] - other
        score[~expressed] = -np.inf
        order = np.argsort(-score, kind="stable")[:top_n]
        order = order[np.isfinite(score[order])]
        per_type[i] = order
        sel.update(order.tolist())
    return np.array(sorted(sel)), per_type


def deconvolve(Rt, mix, gidx):
    """Rt: types x genes reference; mix: (n_mix x genes) CPM mixtures; returns n_mix x types fractions."""
    A = Rt[:, gidx].T                    # genes x types
    out = np.zeros((mix.shape[0], Rt.shape[0])); degenerate = 0
    for k in range(mix.shape[0]):
        w, _ = nnls(A, mix[k, gidx])
        s = w.sum()
        if s <= 0:
            degenerate += 1; out[k] = 1.0 / Rt.shape[0]
        else:
            out[k] = w / s
    return out, degenerate


def metrics(est, true):
    d = est - true
    rmse = float(np.sqrt(np.mean(d ** 2)))
    r = float(np.corrcoef(est.ravel(), true.ravel())[0, 1])
    return rmse, r, d.mean(0)


def main():
    depth_tab = pd.read_csv(APBC / "results/benchmark_depth_by_type.csv")
    depth_tab = depth_tab[depth_tab.model == "proximal"]
    est_rows, met_rows, sig_rows, gene_rows, split_rows = [], [], [], [], []
    for cohort, rel in MATRICES:
        t0 = time.time()
        X, lab, genes, obs = load(cohort, rel)
        types = sorted(pd.unique(lab)); T = len(types)
        rowsum = np.asarray(X.sum(1)).ravel().astype(np.float64)
        dtab = depth_tab[depth_tab.cohort == cohort].set_index("cell_type")
        deepest_true = dtab.median_true_depth.idxmax()
        log(f"{cohort}: {X.shape[0]} cells x {X.shape[1]} genes, {T} types; deepest by true depth = {deepest_true}")
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            tr_idx, te_idx = [], []
            for t in types:
                idx = rng.permutation(np.where(lab == t)[0])
                h = len(idx) // 2
                tr_idx.append(idx[:h]); te_idx.append(idx[h:])
            tr = np.sort(np.concatenate(tr_idx)); te = np.sort(np.concatenate(te_idx))
            Xtr, labtr = X[tr], lab[tr]
            for t in types:
                split_rows.append(dict(cohort=cohort, seed=seed, cell_type=t,
                                       n_train=int((labtr == t).sum()), n_test=int((lab[te] == t).sum()),
                                       median_rowsum_train=float(np.median(rowsum[tr][labtr == t])),
                                       median_rowsum_test=float(np.median(rowsum[te][lab[te] == t])),
                                       median_true_depth=float(dtab.loc[t, "median_true_depth"]),
                                       deepest_by_true_depth=int(t == deepest_true)))
            t1 = time.time()
            refs, dm_target = build_refs(Xtr, labtr, types, rng)
            expressed = np.asarray(Xtr.sum(0)).ravel() > 0
            log(f"  seed {seed}: refs built in {time.time()-t1:.0f}s (depth-matched target {dm_target:.0f}); "
                f"{int(expressed.sum())} genes expressed in training half")
            # gene sets
            common, per_type = top_genes(refs["cpm_mean"], expressed, "ratio")
            for i, t in enumerate(types):
                for rank, g in enumerate(per_type[i]):
                    gene_rows.append(dict(cohort=cohort, seed=seed, cell_type=t, rank=rank + 1, gene=genes[g]))
            gsets = {("common", m): common for m in METHODS}
            for m in METHODS:
                own, _ = top_genes(refs[m], expressed, "diff" if m == "zscale_mean" else "ratio")
                gsets[("own", m)] = own
                sig_rows.append(dict(cohort=cohort, seed=seed, geneset="own", method=m, n_genes=len(own),
                                     overlap_with_common=int(len(np.intersect1d(own, common)))))
            sig_rows.append(dict(cohort=cohort, seed=seed, geneset="common", method="cpm_mean(all)",
                                 n_genes=len(common), overlap_with_common=len(common)))
            # mixtures from the held-out half
            te_by_type = [np.where(lab[te] == t)[0] for t in types]
            true_cell = np.zeros((N_MIX, T)); true_frag = np.zeros((N_MIX, T)); dir_p = np.zeros((N_MIX, T))
            mix_raw = np.zeros((N_MIX, X.shape[1])); mix_cn = np.zeros((N_MIX, X.shape[1]))
            Xte = X[te]; rs_te = rowsum[te]
            for k in range(N_MIX):
                p = rng.dirichlet(np.full(T, ALPHA)); dir_p[k] = p
                n = rng.multinomial(N_CELLS, p)
                rows = []
                for i in range(T):
                    pool = te_by_type[i]
                    if n[i] == 0: continue
                    pick = rng.choice(pool, n[i], replace=(n[i] > len(pool)))
                    rows.append(pick)
                    true_frag[k, i] = rs_te[pick].sum()
                rows = np.concatenate(rows)
                true_cell[k] = n / N_CELLS
                true_frag[k] /= true_frag[k].sum()
                sub = Xte[rows]
                v = np.asarray(sub.sum(0)).ravel().astype(np.float64)
                mix_raw[k] = v / v.sum() * 1e6
                vc = np.asarray(percell_scale(sub).sum(0)).ravel().astype(np.float64)
                mix_cn[k] = vc / vc.sum() * 1e6
            log(f"  seed {seed}: {N_MIX} mixtures built; common signature {len(common)} genes")
            # deconvolution
            for mixname, mix in [("rawsum", mix_raw), ("cellnorm", mix_cn)]:
                for gs in ["common", "own"]:
                    for m in METHODS:
                        gidx = gsets[(gs, m)]
                        est, degen = deconvolve(refs[m], mix, gidx)
                        for k in range(N_MIX):
                            for i, t in enumerate(types):
                                est_rows.append(dict(cohort=cohort, seed=seed, mixture=k, mixture_variant=mixname,
                                                     geneset=gs, method=m, cell_type=t, est=est[k, i],
                                                     true_cellfrac=true_cell[k, i], true_fragfrac=true_frag[k, i],
                                                     dirichlet_p=dir_p[k, i]))
                        for truth, tv in [("cellfrac", true_cell), ("fragfrac", true_frag)]:
                            rmse, r, bias = metrics(est, tv)
                            row = dict(cohort=cohort, seed=seed, mixture_variant=mixname, geneset=gs, method=m,
                                       truth=truth, n_genes=len(gidx), rmse=rmse, pearson_r=r,
                                       n_degenerate=degen, deepest_type=deepest_true,
                                       bias_deepest=float(bias[types.index(deepest_true)]),
                                       max_abs_bias_type=types[int(np.argmax(np.abs(bias)))],
                                       max_abs_bias=float(np.max(np.abs(bias))))
                            for i, t in enumerate(types):
                                row[f"bias__{t}"] = float(bias[i])
                            met_rows.append(row)
                        if mixname == "rawsum" and gs == "common":
                            mc = [r for r in met_rows if r["cohort"] == cohort and r["seed"] == seed and
                                  r["mixture_variant"] == "rawsum" and r["geneset"] == "common" and r["method"] == m]
                            log(f"    {m:19} rawsum/common  cellfrac RMSE {mc[0]['rmse']:.3f} r {mc[0]['pearson_r']:.3f} "
                                f"bias[{deepest_true}] {mc[0]['bias_deepest']:+.3f} | fragfrac RMSE {mc[1]['rmse']:.3f} "
                                f"bias {mc[1]['bias_deepest']:+.3f}")
            del Xtr, Xte, refs
        del X
        log(f"{cohort} done in {time.time()-t0:.0f}s")

    E = pd.DataFrame(est_rows); E.to_csv(HERE / "estimates.csv", index=False)
    M = pd.DataFrame(met_rows); M.to_csv(HERE / "metrics_by_seed.csv", index=False)
    pd.DataFrame(sig_rows).to_csv(HERE / "signature_sizes.csv", index=False)
    pd.DataFrame(gene_rows).to_csv(HERE / "common_signature_genes.csv", index=False)
    pd.DataFrame(split_rows).to_csv(HERE / "split_depth.csv", index=False)
    key = ["cohort", "mixture_variant", "geneset", "method", "truth"]
    S = (M.groupby(key).agg(rmse_mean=("rmse", "mean"), rmse_sd=("rmse", "std"),
                           r_mean=("pearson_r", "mean"), r_sd=("pearson_r", "std"),
                           bias_deepest_mean=("bias_deepest", "mean"), bias_deepest_sd=("bias_deepest", "std"),
                           deepest_type=("deepest_type", "first"), n_genes_mean=("n_genes", "mean"),
                           n_degenerate=("n_degenerate", "sum")).reset_index())
    S.to_csv(HERE / "metrics_summary.csv", index=False)
    # per-type bias, averaged over seeds and mixtures
    B = (E.assign(bias_cell=E.est - E.true_cellfrac, bias_frag=E.est - E.true_fragfrac)
          .groupby(["cohort", "mixture_variant", "geneset", "method", "cell_type"])
          .agg(bias_cellfrac=("bias_cell", "mean"), bias_fragfrac=("bias_frag", "mean"),
               mean_true_cellfrac=("true_cellfrac", "mean"), mean_true_fragfrac=("true_fragfrac", "mean"),
               mean_est=("est", "mean")).reset_index())
    dmap = {c: depth_tab[depth_tab.cohort == c].set_index("cell_type").median_true_depth for c in depth_tab.cohort.unique()}
    B["median_true_depth"] = [dmap[c][t] for c, t in zip(B.cohort, B.cell_type)]
    B["is_deepest"] = [int(t == dmap[c].idxmax()) for c, t in zip(B.cohort, B.cell_type)]
    B.to_csv(HERE / "bias_by_type.csv", index=False)
    print("\n=== summary: rawsum mixtures, common gene set ===")
    for truth in ["cellfrac", "fragfrac"]:
        print(f"--- truth = {truth}")
        s = S[(S.mixture_variant == "rawsum") & (S.geneset == "common") & (S.truth == truth)]
        print(s.pivot(index="method", columns="cohort", values="rmse_mean").loc[METHODS].round(3).to_string())
        print("bias of deepest type:")
        print(s.pivot(index="method", columns="cohort", values="bias_deepest_mean").loc[METHODS].round(3).to_string())
    log("all done")


if __name__ == "__main__":
    main()
