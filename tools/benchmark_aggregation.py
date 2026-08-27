"""
Benchmark: does the way a single-cell gene-activity matrix is aggregated into
cell-type profiles make the per-gene argmax report lineage, or sequencing depth?

Seven aggregation methods, as they are actually used to turn a gene-activity
matrix into a per-cell-type profile (for marker scoring, annotation, specificity
ratios, or a deconvolution reference):

  raw_mean        mean of raw counts per cell type            (naive; what a consumer of a
                                                               deposited matrix does first)
  zscale_mean     scanpy.pp.scale(max_value=10) per gene, then mean   (scaled matrix stored as data)
  pseudobulk_cpm  sum counts per cell type, then CPM          (deconvolution-reference style)
  depth_matched   binomially downsample every cell to the shallowest cell type's median
                  depth, then raw mean                        (control: equalises depth, no normaliser)
  cpm_mean        counts-per-10k per cell, then mean          (Seurat AverageExpression on the
                                                               normalised slot; the step we recommend)
  lognorm_mean    counts-per-10k, log1p, then mean            (scanpy normalize_total+log1p / Seurat
                                                               LogNormalize data slot)
  archr_style     scale to 10k per cell, log2(x+1), then mean (ArchR GeneScoreMatrix convention)

Three read-outs per (matrix, method):
  1. argmax concentration: the share of ALL expressed genes whose per-type argmax
     falls in a single cell type, and which type that is;
  2. Spearman correlation between each type's median TRUE per-cell depth
     (fragments, recomputed from fragment files or peak matrices, never the
     matrix row sum) and its share of the argmaxima;
  3. marker--lineage concordance on a held-out canonical marker panel: markers
     that were NOT used to assign the cell labels (COAD: Script 65 panel excluded;
     LUAD: Script 94 panel excluded; HCC labels come from TF-motif activity, so the
     standard 22-marker panel is already independent), with a 10,000-shuffle
     permutation null on the marker-to-argmax assignment.

Seven count matrices, three cancers, two gene-activity models, one negative
control (the HCC matrix whose rows were named by coordinate order).

HCC labels: the clusters labelled "B_cell" and "DC" by the motif annotation are
relabelled Endothelial_stromal and B_cell respectively (see tools/check_hcc_labels.py
and the comment at the relabel site).

Run:  python3 tools/benchmark_aggregation.py
Out:  results/benchmark_aggregation.csv      one row per matrix x method
      results/benchmark_markers.csv          one row per matrix x method x marker
      results/benchmark_depth_by_type.csv    per-type n cells and median true depth
"""
import os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, h5py, hdf5plugin  # noqa: F401  (plugin needed to read the SnapATAC2 files)
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
from scipy.stats import spearmanr

BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "/mnt/10t/scrna_atac"))
OUT = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
N_PERM = 10_000
SEED = 42

MATRICES = [
    ("HCC",  "misannotated", "data/processed/integrated/atac_hcc_geneactivity.h5ad",          "cell_type"),
    ("HCC",  "proximal",     "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad", "cell_type"),
    ("HCC",  "enhancer",     "data/processed/integrated/atac_hcc_geneactivity_enhancer.h5ad", "cell_type"),
    ("COAD", "proximal",     "data/processed/COAD/coad_geneactivity_proximal.h5ad",           "cell_type"),
    ("COAD", "enhancer",     "data/processed/COAD/coad_geneactivity_enhancer.h5ad",           "cell_type"),
    ("LUAD", "proximal",     "data/processed/LUAD/luad_geneactivity_proximal.h5ad",           "cell_type"),
    ("LUAD", "enhancer",     "data/processed/LUAD/luad_geneactivity_enhancer.h5ad",           "cell_type"),
]

# ── held-out canonical marker panels: marker -> the ONE annotated type it must peak in ──
# HCC: labels are from the deposited TF-motif matrix, so the standard panel is independent.
PANEL = {
    "HCC": {
        "CD68": "Macrophage", "CSF1R": "Macrophage", "MRC1": "Macrophage", "MSR1": "Macrophage",
        "CD163": "Macrophage", "ITGAM": "Macrophage",
        "MS4A1": "B_cell", "CD79A": "B_cell", "CD79B": "B_cell", "CD19": "B_cell", "BANK1": "B_cell",
        "CD3D": "NK_cytotoxic_T", "CD3E": "NK_cytotoxic_T", "CD8A": "NK_cytotoxic_T",
        "IL7R": "NK_cytotoxic_T", "NKG7": "NK_cytotoxic_T", "GNLY": "NK_cytotoxic_T",
        "ALB": "Hepatocyte", "APOB": "Hepatocyte", "HNF4A": "Hepatocyte", "TTR": "Hepatocyte",
        "APOA1": "Hepatocyte",
    },
    # COAD: labels were assigned by scripts/analysis/65 from a 60-gene panel; none of these is in it.
    "COAD": {
        "MS4A7": "Macrophage", "FOLR2": "Macrophage", "LYVE1": "Macrophage", "CD86": "Macrophage",
        "MPEG1": "Macrophage", "SIGLEC1": "Macrophage",
        "BANK1": "B_cell", "CD22": "B_cell", "FCRL1": "B_cell", "BLK": "B_cell",
        "TNFRSF13C": "B_cell", "FCRLA": "B_cell",
        "CD2": "T_NK", "CD247": "T_NK", "SKAP1": "T_NK", "THEMIS": "T_NK", "CD96": "T_NK", "ITK": "T_NK",
        "HDC": "Mast", "MS4A2": "Mast", "GATA2": "Mast", "HPGDS": "Mast", "SLC18A2": "Mast",
        "CDH1": "Epithelial", "CLDN4": "Epithelial", "ELF3": "Epithelial", "VIL1": "Epithelial",
        "LGALS4": "Epithelial", "CDX1": "Epithelial",
        "PDGFRB": "Fibroblast", "FBLN1": "Fibroblast", "COL6A3": "Fibroblast", "COL5A1": "Fibroblast",
        "THY1": "Fibroblast", "MMP2": "Fibroblast",
        "CLDN5": "Endothelial", "FLT1": "Endothelial", "ERG": "Endothelial", "ENG": "Endothelial",
        "TEK": "Endothelial", "PLVAP": "Endothelial",
    },
    # LUAD: labels were assigned by scripts/analysis/94 from a 35-gene panel; none of these is in it.
    "LUAD": {
        "CD163": "Macrophage", "C1QA": "Macrophage", "C1QB": "Macrophage", "CD14": "Macrophage",
        "MS4A7": "Macrophage", "FOLR2": "Macrophage",
        "CD22": "B_cell", "FCRL1": "B_cell", "BLK": "B_cell", "TNFRSF13C": "B_cell",
        "VPREB3": "B_cell", "FCRLA": "B_cell",
        "CD3G": "T_NK", "CD247": "T_NK", "SKAP1": "T_NK", "CD96": "T_NK", "ITK": "T_NK", "LCK": "T_NK",
        "KRT19": "Epithelial", "KRT7": "Epithelial", "CLDN4": "Epithelial", "ELF3": "Epithelial",
        "SFTPB": "Epithelial", "NAPSA": "Epithelial",
        "CLDN5": "Endothelial", "FLT1": "Endothelial", "ERG": "Endothelial", "ENG": "Endothelial",
        "TEK": "Endothelial", "PLVAP": "Endothelial",
        "COL3A1": "Fibroblast", "PDGFRA": "Fibroblast", "PDGFRB": "Fibroblast", "FBLN1": "Fibroblast",
        "COL6A3": "Fibroblast", "COL5A1": "Fibroblast",
    },
}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ── true per-cell depth, keyed by the gene-activity obs_names ────────────────
def true_depth(cohort, obs_names):
    d = {}
    if cohort == "HCC":
        with h5py.File(BASE / "data/processed/integrated/_s90_work/hcc_fragments.h5ad") as f:
            idx = [x.decode() for x in f["obs/index"][:]]
            d = dict(zip(idx, f["obs/n_fragment"][:].astype(float)))
    elif cohort == "COAD":
        for p in sorted((BASE / "data/processed/COAD/_s92_work").glob("GSM*.h5ad")):
            with h5py.File(p) as f:
                idx = [f"{p.stem}_{x.decode()}" for x in f["obs/index"][:]]
                d.update(zip(idx, f["obs/n_fragment"][:].astype(float)))
    elif cohort == "LUAD":
        gsms = sorted({o.split("_")[0] for o in obs_names})
        for gsm in gsms:
            h5 = next((BASE / "data/raw/scATAC_LUAD/GSE270148").glob(f"{gsm}_*filtered_peak_bc_matrix.h5"))
            a = sc.read_10x_h5(h5, gex_only=False)
            tot = np.asarray(a.X.sum(1)).ravel()
            d.update(zip([f"{gsm}_{b}" for b in a.obs_names], tot.astype(float)))
    dep = np.array([d.get(o, np.nan) for o in obs_names])
    return dep


# ── aggregation methods: X (cells x genes, csr float32), lab (str array), types ──
def type_mean(M, lab, types):
    return np.vstack([np.asarray(M[lab == t].mean(0)).ravel() for t in types])


def percell_scale(X, target=1e4):
    tot = np.asarray(X.sum(1)).ravel().astype(np.float64)
    tot[tot == 0] = 1.0
    return sp.diags((target / tot).astype(np.float32)) @ X


def m_raw(X, lab, types):
    return type_mean(X, lab, types)


def m_pseudobulk_cpm(X, lab, types):
    S = np.vstack([np.asarray(X[lab == t].sum(0)).ravel() for t in types]).astype(np.float64)
    return S / S.sum(1, keepdims=True) * 1e6


def m_cpm(X, lab, types):
    return type_mean(percell_scale(X), lab, types)


def m_lognorm(X, lab, types):
    Xn = percell_scale(X).tocsr(); Xn.data = np.log1p(Xn.data)
    return type_mean(Xn, lab, types)


def m_archr(X, lab, types):
    Xn = percell_scale(X).tocsr(); Xn.data = np.log2(Xn.data + 1)
    return type_mean(Xn, lab, types)


def m_zscale(X, lab, types, max_value=10.0, block=2000):
    """scanpy.pp.scale(zero_center=True, max_value=10) per gene, then the per-type mean.
    Done in dense gene blocks because the clip is entry-wise and breaks sparsity."""
    n, g = X.shape
    masks = [(lab == t) for t in types]
    out = np.zeros((len(types), g), dtype=np.float64)
    Xc = X.tocsc()
    for j in range(0, g, block):
        D = Xc[:, j:j + block].toarray().astype(np.float64)
        mu = D.mean(0); sd = D.std(0, ddof=1); sd[sd == 0] = 1.0
        D = (D - mu) / sd
        D[D > max_value] = max_value
        for i, m in enumerate(masks):
            out[i, j:j + block] = D[m].mean(0)
    return out


def m_depth_matched(X, lab, types, rng):
    """Binomially thin every cell to the shallowest type's median matrix depth; raw mean.
    No normaliser is applied: this is the control that shows the effect is depth."""
    tot = np.asarray(X.sum(1)).ravel().astype(np.float64)
    target = min(np.median(tot[lab == t]) for t in types)
    p = np.minimum(1.0, target / np.maximum(tot, 1.0))
    Xd = X.tocsr(copy=True)
    per_nz = np.repeat(p, np.diff(Xd.indptr))
    Xd.data = rng.binomial(np.rint(Xd.data).astype(np.int64), per_nz).astype(np.float32)
    Xd.eliminate_zeros()
    return type_mean(Xd, lab, types), target


# ── permutation null for marker concordance ─────────────────────────────────
def perm_p(argmax, target, obs, rng, n=N_PERM):
    """One-sided permutation p with the (1 + k) / (B + 1) estimator, so that the
    smallest reportable value is 1 / (B + 1) rather than exactly zero."""
    null = np.array([(rng.permutation(argmax) == target).sum() for _ in range(n)])
    p = (1.0 + float((null >= obs).sum())) / (n + 1.0)
    return p, float(null.mean()), float(np.percentile(null, 95))


def main():
    rng = np.random.default_rng(SEED)
    rows, mrows, drows = [], [], []
    for cohort, model, rel, labcol in MATRICES:
        path = BASE / rel
        if not path.exists():
            log(f"skip {cohort}/{model}: {path} missing"); continue
        log(f"{cohort}/{model}: loading")
        a = ad.read_h5ad(path)
        lab_all = a.obs[labcol].astype(str).values
        keep = ~pd.isna(a.obs[labcol].values) & (lab_all != "nan")
        a = a[keep].copy(); lab = a.obs[labcol].astype(str).values
        if cohort == "HCC":
            # The motif-derived labels of GSE227265 name two clusters wrongly: the cluster
            # labelled "B_cell" carries the endothelial/stromal programme (PECAM1, VWF, CDH5,
            # KDR, CLDN5, FLT1, ENG, TEK, ERG, COL1A2, COL3A1, DCN, PDGFRB all peak there and
            # no B-cell gene does), and the cluster labelled "DC" carries the B/plasma programme
            # (PAX5, EBF1, BLK, CD22, FCRL1, TNFRSF13C, JCHAIN, MZB1, IGHG1) while bona fide
            # DC markers (CLEC9A, XCR1, CD1C, CLEC10A, LAMP3, ITGAX, ZBTB46) peak in the
            # macrophage compartment. Evidence: tools/check_hcc_labels.py -> results/hcc_label_check.csv,
            # using markers disjoint from the test panel. Relabelled accordingly.
            lab = np.where(lab == "B_cell", "Endothelial_stromal", np.where(lab == "DC", "B_cell", lab))
        X = (a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)).tocsr().astype(np.float32)
        genes = a.var_names.astype(str).values
        is_counts = bool(X.min() >= 0 and np.allclose(X.data[:100000], np.rint(X.data[:100000])))
        types = sorted(pd.unique(lab))
        dep = true_depth(cohort, a.obs_names.values)
        log(f"  {X.shape[0]} cells x {X.shape[1]} genes, {len(types)} types, counts={is_counts}, "
            f"depth mapped for {np.isfinite(dep).mean():.1%} of cells")
        med_depth = np.array([np.nanmedian(dep[lab == t]) for t in types])
        med_rowsum = np.array([np.median(np.asarray(X[lab == t].sum(1)).ravel()) for t in types])
        nnz = np.diff(X.indptr) / X.shape[1]
        det_rate = np.array([float(nnz[lab == t].mean()) for t in types])
        for t, md, mr, dr in zip(types, med_depth, med_rowsum, det_rate):
            drows.append(dict(cohort=cohort, model=model, cell_type=t, n_cells=int((lab == t).sum()),
                              median_true_depth=md, median_matrix_rowsum=mr, detection_rate=dr))
        expressed = np.asarray(X.sum(0)).ravel() > 0
        panel = {g: t for g, t in PANEL[cohort].items() if g in set(genes) and t in types}
        gidx = {g: i for i, g in enumerate(genes)}

        methods = [("raw_mean", lambda: m_raw(X, lab, types)),
                   ("zscale_mean", lambda: m_zscale(X, lab, types)),
                   ("pseudobulk_cpm", lambda: m_pseudobulk_cpm(X, lab, types)),
                   ("depth_matched", lambda: m_depth_matched(X, lab, types, rng)),
                   ("cpm_mean", lambda: m_cpm(X, lab, types)),
                   ("lognorm_mean", lambda: m_lognorm(X, lab, types)),
                   ("archr_style", lambda: m_archr(X, lab, types))]
        for name, fn in methods:
            t0 = time.time()
            res = fn()
            P, extra = (res if isinstance(res, tuple) else (res, None))
            arg = P[:, expressed].argmax(0)
            share = np.bincount(arg, minlength=len(types)) / arg.size
            top = int(share.argmax())
            rho, pv = spearmanr(med_depth, share)
            rho_det, _ = spearmanr(det_rate, share)
            # held-out marker concordance
            marg = np.array([types[int(P[:, gidx[g]].argmax())] for g in panel])
            mtgt = np.array(list(panel.values()))
            correct = (marg == mtgt)
            obs = int(correct.sum())
            pp, nmean, n95 = perm_p(marg, mtgt, obs, rng) if len(panel) else (np.nan, np.nan, np.nan)
            per_lin = {t: f"{int(correct[mtgt == t].sum())}/{int((mtgt == t).sum())}" for t in types if (mtgt == t).any()}
            rows.append(dict(cohort=cohort, model=model, method=name, n_cells=X.shape[0],
                             n_genes_expressed=int(expressed.sum()), n_types=len(types),
                             top_share=float(share[top]), top_type=types[top],
                             top_type_median_depth=float(med_depth[top]),
                             deepest_type=types[int(np.nanargmax(med_depth))],
                             rho_depth_share=float(rho), p_depth_share=float(pv), rho_detection_share=float(rho_det),
                             markers_correct=obs, markers_total=len(panel),
                             perm_p=pp, null_mean=nmean, null_95=n95,
                             per_lineage="; ".join(f"{k} {v}" for k, v in per_lin.items()),
                             depth_matched_target=(float(extra) if extra is not None else np.nan),
                             seconds=round(time.time() - t0, 1)))
            for g, tg, ag, c in zip(panel, mtgt, marg, correct):
                mrows.append(dict(cohort=cohort, model=model, method=name, gene=g,
                                  target=tg, argmax=ag, correct=int(c)))
            log(f"  {name:15} top {types[top]:>14} {share[top]:5.1%}  rho={rho:+.2f} p={pv:.3f}  "
                f"markers {obs}/{len(panel)} perm p={pp:.4f}   ({time.time()-t0:.0f}s)")
        del a, X

    R = pd.DataFrame(rows); R.to_csv(OUT / "benchmark_aggregation.csv", index=False)
    pd.DataFrame(mrows).to_csv(OUT / "benchmark_markers.csv", index=False)
    pd.DataFrame(drows).to_csv(OUT / "benchmark_depth_by_type.csv", index=False)
    print("\n=== summary across matrices (excluding the mis-annotated negative control) ===")
    S = R[R.model != "misannotated"]
    for name in [m for m, _ in methods]:
        s = S[S.method == name]
        print(f"  {name:15} top share median {s.top_share.median():5.1%} ({s.top_share.min():.0%}-{s.top_share.max():.0%})"
              f"   argmax in deepest type {(s.top_type == s.deepest_type).sum()}/{len(s)}"
              f"   rho median {s.rho_depth_share.median():+.2f}, sig {(s.p_depth_share < 0.05).sum()}/{len(s)}"
              f"   markers {s.markers_correct.sum()}/{s.markers_total.sum()}"
              f"   perm p<0.05 {(s.perm_p < 0.05).sum()}/{len(s)}")
    print("\n=== negative control (HCC, rows named by coordinate order) ===")
    for r in R[R.model == "misannotated"].itertuples():
        print(f"  {r.method:15} markers {r.markers_correct}/{r.markers_total} perm p={r.perm_p:.3f}   top share {r.top_share:.1%}")
    print(f"\nwrote {OUT/'benchmark_aggregation.csv'}")


if __name__ == "__main__":
    main()
