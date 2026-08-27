"""
03: Aggregation benchmark on the 10x PBMC multiome with a GENOME-WIDE RNA ground truth.

Cells are labelled from RNA only (02_label_rna.py). For each ATAC gene-activity model
(proximal = Signac model; enhancer = ArchR-style distance-weighted) and each aggregation
method (functions imported from tools/benchmark_aggregation.py) we compute:

  top_share / top_type           share of ATAC-expressed genes whose per-type argmax falls in one type
  rho_depth_share                Spearman(median true per-cell fragments, argmax share) across types
  argmax_agree_all               (a) fraction of genes expressed in BOTH modalities whose ATAC argmax
                                 type == RNA argmax type (RNA profile = CPM mean of RNA counts)
  argmax_agree_informative       same, restricted to genes whose RNA argmax is well defined
                                 (top type >= 2x the second type in RNA CPM mean)
  gene_rho_mean                  (b) mean per-gene Spearman between the ATAC and RNA per-type profiles
                                 (across the 8 types), over genes expressed in both
  markers_correct/total, perm_p  (c) held-out PBMC canonical panel, disjoint from the labelling genes
  annot_acc_*                    annotation task: reference profiles from a 50% training split, each
                                 held-out cell (log1p CPM) assigned to the best-correlated reference

Out: multiome_benchmark.csv, multiome_markers.csv, multiome_depth_by_type.csv, multiome_annotation.csv
"""
import sys, time, warnings
import os
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, hdf5plugin  # noqa
import scanpy as sc
from scipy.stats import spearmanr, rankdata
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from benchmark_aggregation import (percell_scale, type_mean, m_raw, m_zscale, m_pseudobulk_cpm,
                                   m_depth_matched, m_cpm, m_lognorm, m_archr, perm_p)  # noqa

HERE = Path(__file__).resolve().parent
BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
PROC = BASE / "data/processed/multiome_pbmc10k"
H5   = BASE / "data/raw/multiome_pbmc10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
SEED = 42
MODELS = [("proximal", PROC / "pbmc10k_geneactivity_proximal.h5ad"),
          ("enhancer", PROC / "pbmc10k_geneactivity_enhancer.h5ad")]

# held-out PBMC panel: gene -> set of acceptable RNA-defined types. None of these genes is in
# labelling_genes.txt (asserted below). Pan-lineage genes accept every type of that lineage.
T = {"T_CD4", "T_CD8"}; MYE = {"Mono_CD14", "Mono_CD16", "DC"}
PANEL = {
    **{g: T for g in ["CD3G", "CD2", "CD5", "CD6", "LCK", "ITK", "SKAP1", "CD247", "LEF1", "TCF7", "CD28", "TRAT1", "BCL11B"]},
    **{g: {"NK"} for g in ["KLRF1", "NCR1", "KLRC1", "SH2D1B", "NCR3"]},
    **{g: {"Mono_CD14"} for g in ["VCAN", "S100A12", "CSF3R", "CD36", "CCR2"]},
    **{g: {"Mono_CD16"} for g in ["CDKN1C", "TCF7L2", "CX3CR1", "LYPD2", "HES4"]},
    **{g: {"B"} for g in ["PAX5", "EBF1", "BANK1", "CD22", "FCRL1", "TNFRSF13C", "BLK", "VPREB3"]},
    **{g: {"DC"} for g in ["CD1E", "PKIB", "ENHO", "CLEC9A"]}, "FLT3": {"DC", "pDC"},
    **{g: {"pDC"} for g in ["PLD4", "SERPINF1", "PTPRS", "DERL3", "TCF4", "IRF7", "ITM2C"]},
    "CD68": MYE, "CSF1R": {"Mono_CD14", "Mono_CD16"}, "ITGAX": MYE,
}

def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

labelling = set(Path(HERE / "labelling_genes.txt").read_text().split())
assert not (set(PANEL) & labelling), set(PANEL) & labelling
lab_df = pd.read_csv(HERE / "rna_labels.csv", index_col=0)
lab_df = lab_df[lab_df.cell_type != "dropped"]
types = sorted(lab_df.cell_type.unique())
log(f"{len(lab_df)} labelled cells, {len(types)} types: {types}")

# ── RNA truth: per-type CPM mean of RNA counts ───────────────────────────────
rna = sc.read_10x_h5(H5, gex_only=True); rna.var_names_make_unique()
rna = rna[lab_df.index].copy()
Xr = sp.csr_matrix(rna.X).astype(np.float32)
lab_r = lab_df.cell_type.values.astype(str)
R = m_cpm(Xr, lab_r, types)                              # types x genes (RNA CPM-mean profile)
rna_det = np.asarray((Xr > 0).mean(0)).ravel()           # detection fraction per gene
rna_genes = rna.var_names.astype(str).values
rna_arg = R.argmax(0)
srt = np.sort(R, 0)
rna_informative = (srt[-1] >= 2 * np.maximum(srt[-2], 1e-9)) & (srt[-1] > 0)
log(f"RNA: {Xr.shape}, genes detected in >=1% cells: {(rna_det >= 0.01).sum()}, informative: {rna_informative.sum()}")

def annotate(P, genes_mask, Q_log, Q_rank, ylab, types):
    """assign each query cell to the reference row with the highest correlation.
    Pearson on values (query = log1p CPM) and Spearman (both ranked)."""
    Pm = P[:, genes_mask]
    # Pearson: centre rows
    def corr_assign(Qc, Pc):
        Qc = Qc - Qc.mean(1, keepdims=True); Pc = Pc - Pc.mean(1, keepdims=True)
        Qc /= (np.linalg.norm(Qc, axis=1, keepdims=True) + 1e-12); Pc /= (np.linalg.norm(Pc, axis=1, keepdims=True) + 1e-12)
        return (Qc @ Pc.T).argmax(1)
    pred_p = corr_assign(Q_log.copy(), Pm.astype(np.float64).copy())
    Pr = np.vstack([rankdata(r) for r in Pm]).astype(np.float64)
    pred_s = corr_assign(Q_rank.copy(), Pr)
    out = {}
    for tag, pred in [("pearson", pred_p), ("spearman", pred_s)]:
        y = np.array(types)[pred]
        out[f"annot_acc_{tag}"] = float((y == ylab).mean())
        out[f"annot_balacc_{tag}"] = float(np.mean([(y[ylab == t] == t).mean() for t in types]))
        out[f"annot_pred_{tag}"] = y
    return out

rng = np.random.default_rng(SEED)
rows, mrows, drows, arows = [], [], [], []
for model, path in MODELS:
    log(f"== {model}: loading {path.name}")
    a = ad.read_h5ad(path)
    a = a[lab_df.index.intersection(a.obs_names)].copy()
    lab = lab_df.loc[a.obs_names, "cell_type"].values.astype(str)
    dep = a.obs["n_fragment"].values.astype(float)
    X = sp.csr_matrix(a.X).astype(np.float32).tocsr()
    genes = a.var_names.astype(str).values
    gidx = {g: i for i, g in enumerate(genes)}
    med_depth = np.array([np.median(dep[lab == t]) for t in types])
    nnz = np.diff(X.indptr) / X.shape[1]
    det_rate = np.array([float(nnz[lab == t].mean()) for t in types])
    rowsum = np.asarray(X.sum(1)).ravel()
    for t, md, dr in zip(types, med_depth, det_rate):
        drows.append(dict(model=model, cell_type=t, n_cells=int((lab == t).sum()), median_true_depth=md,
                          median_matrix_rowsum=float(np.median(rowsum[lab == t])), detection_rate=dr,
                          median_rna_umi=float(lab_df.loc[a.obs_names[lab == t], "rna_umi"].median())))
    log(f"  {X.shape}, depth spread {med_depth.max()/med_depth.min():.2f}x "
        f"(min {types[int(med_depth.argmin())]} {med_depth.min():.0f}, max {types[int(med_depth.argmax())]} {med_depth.max():.0f})")
    expressed = np.asarray(X.sum(0)).ravel() > 0
    # genes shared with RNA truth
    rmap = {g: i for i, g in enumerate(rna_genes)}
    shared = [(i, rmap[g]) for i, g in enumerate(genes) if g in rmap and expressed[i] and rna_det[rmap[g]] >= 0.01]
    ai = np.array([s[0] for s in shared]); ri = np.array([s[1] for s in shared])
    inf_mask = rna_informative[ri]
    Rs = R[:, ri]
    log(f"  shared expressed genes: {len(ai)}; RNA-informative among them: {inf_mask.sum()}")
    panel = {g: s for g, s in PANEL.items() if g in gidx}
    # annotation split (stratified 50/50)
    train = np.zeros(X.shape[0], bool)
    for t in types:
        idx = np.where(lab == t)[0]; rng.shuffle(idx); train[idx[: len(idx) // 2]] = True
    test = ~train
    Xn_test = percell_scale(X[test]).tocsr(); Xn_test.data = np.log1p(Xn_test.data)
    ann_genes = np.asarray(X[train].sum(0)).ravel() > 0
    Q_log = Xn_test[:, ann_genes].toarray().astype(np.float64)
    Q_rank = np.vstack([rankdata(r) for r in Q_log]).astype(np.float64)
    ylab = lab[test]
    Xtr, labtr = X[train], lab[train]

    methods = [("raw_mean", lambda M, l: m_raw(M, l, types)),
               ("zscale_mean", lambda M, l: m_zscale(M, l, types)),
               ("lognorm_mean", lambda M, l: m_lognorm(M, l, types)),
               ("pseudobulk_cpm", lambda M, l: m_pseudobulk_cpm(M, l, types)),
               ("depth_matched", lambda M, l: m_depth_matched(M, l, types, rng)),
               ("cpm_mean", lambda M, l: m_cpm(M, l, types)),
               ("archr_style", lambda M, l: m_archr(M, l, types))]
    for name, fn in methods:
        t0 = time.time()
        res = fn(X, lab); P, extra = (res if isinstance(res, tuple) else (res, None))
        arg = P[:, expressed].argmax(0)
        share = np.bincount(arg, minlength=len(types)) / arg.size
        top = int(share.argmax())
        rho, pv = spearmanr(med_depth, share); rho_det, _ = spearmanr(det_rate, share)
        # (a) argmax agreement with RNA
        Ps = P[:, ai]; a_arg = Ps.argmax(0); r_arg = rna_arg[ri]
        agree_all = float((a_arg == r_arg).mean()); agree_inf = float((a_arg[inf_mask] == r_arg[inf_mask]).mean())
        # (b) per-gene Spearman across types (vectorised: rank columns, then Pearson)
        def colrank(M): return np.apply_along_axis(rankdata, 0, M)
        Pa, Rr = colrank(Ps), colrank(Rs)
        Pa = Pa - Pa.mean(0); Rr = Rr - Rr.mean(0)
        g_rho = (Pa * Rr).sum(0) / (np.sqrt((Pa**2).sum(0) * (Rr**2).sum(0)) + 1e-12)
        g_rho_mean = float(np.nanmean(g_rho)); g_rho_inf = float(np.nanmean(g_rho[inf_mask]))
        # per-type profile correlation across genes (log1p of both, Spearman)
        type_rho = float(np.mean([spearmanr(Ps[k], Rs[k])[0] for k in range(len(types))]))
        # (c) held-out markers
        marg = np.array([types[int(P[:, gidx[g]].argmax())] for g in panel])
        correct = np.array([m in panel[g] for g, m in zip(panel, marg)])
        obs = int(correct.sum())
        null = np.array([sum(m in panel[g] for g, m in zip(panel, rng.permutation(marg))) for _ in range(10_000)])
        pp = float(np.mean(null >= obs)); nmean = float(null.mean())
        # annotation task
        res_tr = fn(Xtr, labtr); Ptr = res_tr[0] if isinstance(res_tr, tuple) else res_tr
        ann = annotate(Ptr, ann_genes, Q_log, Q_rank, ylab, types)
        for tag in ("pearson", "spearman"):
            y = ann.pop(f"annot_pred_{tag}")
            for t in types:
                arows.append(dict(model=model, method=name, corr=tag, cell_type=t, n_test=int((ylab == t).sum()),
                                  recall=float((y[ylab == t] == t).mean()),
                                  confusion="; ".join(f"{u}:{int((y[ylab == t] == u).sum())}" for u in types)))
        rows.append(dict(model=model, method=name, n_cells=X.shape[0], n_genes_expressed=int(expressed.sum()),
                         n_types=len(types), top_share=float(share[top]), top_type=types[top],
                         deepest_type=types[int(med_depth.argmax())], depth_spread=float(med_depth.max() / med_depth.min()),
                         rho_depth_share=float(rho), p_depth_share=float(pv), rho_detection_share=float(rho_det),
                         n_shared_genes=len(ai), n_informative_genes=int(inf_mask.sum()),
                         argmax_agree_all=agree_all, argmax_agree_informative=agree_inf,
                         gene_rho_mean=g_rho_mean, gene_rho_mean_informative=g_rho_inf, type_profile_rho=type_rho,
                         markers_correct=obs, markers_total=len(panel), perm_p=pp, null_mean=nmean,
                         **ann, n_train=int(train.sum()), n_test=int(test.sum()), seconds=round(time.time() - t0, 1)))
        for g, m, c in zip(panel, marg, correct):
            mrows.append(dict(model=model, method=name, gene=g, target="|".join(sorted(panel[g])), argmax=m, correct=int(c)))
        log(f"  {name:15} top {types[top]:>9} {share[top]:5.1%} rho={rho:+.2f} | RNA argmax agree {agree_all:.3f} "
            f"(inf {agree_inf:.3f}) gene-rho {g_rho_mean:+.3f} | markers {obs}/{len(panel)} p={pp:.3f} | "
            f"annot P {ann['annot_acc_pearson']:.3f} S {ann['annot_acc_spearman']:.3f} ({time.time()-t0:.0f}s)")
    del a, X, Xn_test, Q_log, Q_rank

pd.DataFrame(rows).to_csv(HERE / "multiome_benchmark.csv", index=False)
pd.DataFrame(mrows).to_csv(HERE / "multiome_markers.csv", index=False)
pd.DataFrame(drows).to_csv(HERE / "multiome_depth_by_type.csv", index=False)
pd.DataFrame(arows).to_csv(HERE / "multiome_annotation.csv", index=False)
log("done")
