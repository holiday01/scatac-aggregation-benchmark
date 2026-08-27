"""
Additional per-cell normalisers for the aggregation benchmark (journal addition).

Same six real gene-activity matrices, same labels (HCC relabel rule applied after
loading), same held-out marker panels, same true-depth source and the same
permutation null as tools/benchmark_aggregation.py, from which PANEL, MATRICES,
true_depth(), perm_p() and the original six methods are imported unchanged.
The mis-annotated HCC negative control is excluded.

New per-type profiles (each is "transform per cell, then mean per type" unless stated):

  tfidf1_mean       Signac RunTFIDF method 1: log1p(TF x IDF x 1e4),
                    TF = count / cell total, IDF = n_cells / gene total       (Stuart 2021)
  tfidf2_mean       Signac method 2: TF x log(1 + IDF)   -- linear in TF, so per gene it is
                    the CPM mean times a constant: argmax provably identical to cpm_mean
  tfidf3_mean       Signac method 3: log1p(TF x 1e4) x log(1 + IDF) -- per gene the
                    log1p(CPM) mean times a constant: argmax provably identical to lognorm_mean
  pearson_mean      analytic Pearson residuals (Lause 2021; scanpy.experimental.pp.
                    normalize_pearson_residuals, theta = 100, clip = sqrt(n_cells)),
                    computed on the FULL matrix in dense gene blocks (exact, because the
                    residual of entry (i,g) needs only the cell total, the gene total and
                    the grand total); verified against scanpy on a sub-matrix at run time
  medscale_mean     each cell scaled linearly to the median cell total (scanpy
                    normalize_total(target_sum=None)); per gene = cpm_mean x constant
  rank_mean         each cell's non-zero entries replaced by their within-cell rank
                    (average ties; zeros stay 0), then mean
  rankfrac_mean     as rank_mean but ranks divided by the cell's number of non-zeros
                    (so every cell's top gene is 1)
  detection_mean    fraction of cells in the type with count > 0 (binarised mean)
  lognorm_maxnorm   log1p(CPM) mean, then each gene divided by its maximum over types
                    (per-gene monotone rescaling AFTER aggregation: argmax unchanged)
  pseudobulk_logcpm pseudobulk sum per type, CPM, then log1p (log AFTER aggregation:
                    monotone per gene, argmax unchanged)

The original six (raw_mean, zscale_mean, pseudobulk_cpm, depth_matched, cpm_mean,
lognorm_mean) are recomputed here only so that the per-gene argmax vectors of all
methods are available for the identity checks in --merge; the published numbers remain
those of results/benchmark_aggregation.csv.

Run:  python3 benchmark_extra_normalisers.py --matrix i      (i = 0..5, one real matrix)
      python3 benchmark_extra_normalisers.py --merge
Out:  benchmark_extra_normalisers.csv   new methods, same columns as benchmark_aggregation.csv
      benchmark_all_methods.csv         original six (recomputed) + new, plus argmax_identical_to
      benchmark_extra_markers.csv       one row per matrix x method x marker
"""
import argparse, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, h5py, hdf5plugin  # noqa: F401
import scipy.sparse as sp
import anndata as ad
from scipy.stats import spearmanr, rankdata

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parents[2] / "tools"               # <repo>/tools
sys.path.insert(0, str(TOOLS))
import benchmark_aggregation as B               # noqa: E402  (PANEL, MATRICES, true_depth, perm_p, methods)

PARTS = HERE / "parts"; PARTS.mkdir(exist_ok=True)
REAL = [m for m in B.MATRICES if m[1] != "misannotated"]   # the six real matrices
SEED = 42
COLS = ["cohort", "model", "method", "n_cells", "n_genes_expressed", "n_types", "top_share", "top_type",
        "top_type_median_depth", "deepest_type", "rho_depth_share", "p_depth_share", "rho_detection_share",
        "markers_correct", "markers_total", "perm_p", "null_mean", "null_95", "per_lineage",
        "depth_matched_target", "seconds"]


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ── loading: identical to tools/benchmark_aggregation.py main(), incl. the HCC relabel ──
def _read_index(g):
    key = g.attrs.get("_index", "_index")
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[key][:]])


def _read_labels(g, col):
    """anndata categorical / string column -> str array ('nan' for missing), as
    a.obs[col].astype(str) would give."""
    node = g[col]
    if isinstance(node, h5py.Group):                       # categorical encoding
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in node["categories"][:]])
        codes = node["codes"][:]
        out = np.where(codes >= 0, cats[np.maximum(codes, 0)], "nan")
        return out.astype(str)
    vals = node[:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in vals])


def read_csr_float32(path, chunk=50_000_000):
    """CSR gene-activity matrix straight from the h5ad with h5py: data is read in chunks
    into a preallocated float32 array (no float64 copy, no anndata subsetting copy)."""
    with h5py.File(path) as f:
        X = f["X"]; assert X.attrs["encoding-type"] == "csr_matrix", X.attrs["encoding-type"]
        shape = tuple(int(v) for v in X.attrs["shape"])
        nnz = X["data"].shape[0]
        data = np.empty(nnz, dtype=np.float32)
        for s_ in range(0, nnz, chunk):
            data[s_:s_ + chunk] = X["data"][s_:s_ + chunk]
        indices = X["indices"][:].astype(np.int32, copy=False)
        indptr = X["indptr"][:].astype(np.int64, copy=False)
        obs_names = _read_index(f["obs"]); var_names = _read_index(f["var"])
        return sp.csr_matrix((data, indices, indptr), shape=shape), obs_names, var_names


def load_matrix(i):
    cohort, model, rel, labcol = REAL[i]
    path = B.BASE / rel
    log(f"{cohort}/{model}: loading {path}")
    with h5py.File(path) as f:
        lab = _read_labels(f["obs"], labcol)
    X, obs_names, genes = read_csr_float32(path)
    keep = lab != "nan"                                    # same filter as the benchmark
    if not keep.all():
        X = X[keep]; obs_names = obs_names[keep]; lab = lab[keep]
    if cohort == "HCC":
        # same relabel rule as tools/benchmark_aggregation.py (tools/check_hcc_labels.py)
        lab = np.where(lab == "B_cell", "Endothelial_stromal", np.where(lab == "DC", "B_cell", lab))
    if not X.has_sorted_indices:
        X.sort_indices()
    types = sorted(pd.unique(lab))
    dep = B.true_depth(cohort, obs_names)
    med_depth = np.array([np.nanmedian(dep[lab == t]) for t in types])
    nnz = np.diff(X.indptr) / X.shape[1]
    det_rate = np.array([float(nnz[lab == t].mean()) for t in types])
    expressed = np.asarray(X.sum(0)).ravel() > 0
    panel = {g: t for g, t in B.PANEL[cohort].items() if g in set(genes) and t in types}
    gidx = {g: k for k, g in enumerate(genes)}
    is_counts = bool(np.allclose(X.data[:100000], np.rint(X.data[:100000])))
    log(f"  {X.shape[0]} cells x {X.shape[1]} genes, nnz {X.nnz:,}, {len(types)} types, counts={is_counts}, "
        f"depth mapped {np.isfinite(dep).mean():.1%}, panel {len(panel)}")
    return dict(cohort=cohort, model=model, X=X, lab=lab, types=types, genes=genes, dep=dep,
                med_depth=med_depth, det_rate=det_rate, expressed=expressed, panel=panel, gidx=gidx)


# ── read-outs: identical to the loop body of tools/benchmark_aggregation.py ──────────
def evaluate(P, ctx, rng, name, seconds=np.nan, extra=None):
    types, expressed, panel, gidx, med_depth, det_rate = (ctx[k] for k in
        ("types", "expressed", "panel", "gidx", "med_depth", "det_rate"))
    arg = P[:, expressed].argmax(0)
    share = np.bincount(arg, minlength=len(types)) / arg.size
    top = int(share.argmax())
    rho, pv = spearmanr(med_depth, share)
    rho_det, _ = spearmanr(det_rate, share)
    marg = np.array([types[int(P[:, gidx[g]].argmax())] for g in panel])
    mtgt = np.array(list(panel.values()))
    correct = (marg == mtgt); obs = int(correct.sum())
    pp, nmean, n95 = B.perm_p(marg, mtgt, obs, rng)
    per_lin = {t: f"{int(correct[mtgt == t].sum())}/{int((mtgt == t).sum())}" for t in types if (mtgt == t).any()}
    row = dict(cohort=ctx["cohort"], model=ctx["model"], method=name, n_cells=ctx["X"].shape[0],
               n_genes_expressed=int(expressed.sum()), n_types=len(types),
               top_share=float(share[top]), top_type=types[top], top_type_median_depth=float(med_depth[top]),
               deepest_type=types[int(np.nanargmax(med_depth))],
               rho_depth_share=float(rho), p_depth_share=float(pv), rho_detection_share=float(rho_det),
               markers_correct=obs, markers_total=len(panel), perm_p=pp, null_mean=nmean, null_95=n95,
               per_lineage="; ".join(f"{k} {v}" for k, v in per_lin.items()),
               depth_matched_target=(float(extra) if extra is not None else np.nan), seconds=round(seconds, 1))
    mrows = [dict(cohort=ctx["cohort"], model=ctx["model"], method=name, gene=g, target=tg, argmax=ag, correct=int(c))
             for g, tg, ag, c in zip(panel, mtgt, marg, correct)]
    full_arg = P.argmax(0).astype(np.int16)          # per-gene argmax over ALL genes (identity checks)
    log(f"  {name:18} top {types[top]:>20} {share[top]:5.1%}  rho={rho:+.2f}  markers {obs}/{len(panel)} "
        f"perm p={pp:.4f}  ({seconds:.0f}s)")
    return row, mrows, full_arg


# ── new normalisers ───────────────────────────────────────────────────────────
def tfidf_matrix(X, method=1, scale=1e4, chunk=20_000_000):
    """Signac RunTFIDF.default (v1.17.1), transposed to cells x genes. The result shares
    X's index arrays (only a new data array is allocated, filled in chunks)."""
    tot = np.asarray(X.sum(1)).ravel().astype(np.float64); tot[tot == 0] = 1.0
    gsum = np.asarray(X.sum(0)).ravel().astype(np.float64)
    idf = np.zeros_like(gsum); nz = gsum > 0
    idf[nz] = X.shape[0] / gsum[nz]
    gfac = (idf if method == 1 else np.log1p(idf)).astype(np.float32)
    inv_tot = (1.0 / tot).astype(np.float32)
    row = np.repeat(np.arange(X.shape[0], dtype=np.int32), np.diff(X.indptr))
    d = np.empty(X.nnz, dtype=np.float32)
    for s in range(0, X.nnz, chunk):
        e = min(s + chunk, X.nnz)
        tf = X.data[s:e] * inv_tot[row[s:e]]
        if method == 1:
            d[s:e] = np.log1p(tf * gfac[X.indices[s:e]] * np.float32(scale))
        elif method == 2:
            d[s:e] = tf * gfac[X.indices[s:e]]
        elif method == 3:
            d[s:e] = np.log1p(tf * np.float32(scale)) * gfac[X.indices[s:e]]
    d[~np.isfinite(d)] = 0
    return sp.csr_matrix((d, X.indices, X.indptr), shape=X.shape)


def pearson_residual_blocks(X, theta=100.0, clip=None, block=1000):
    """Analytic Pearson residuals exactly as scanpy.experimental.pp.normalize_pearson_residuals
    (theta=100, clip=None -> sqrt(n_cells)), yielded in dense gene blocks. mu[i,g] =
    cell_total[i] * gene_total[g] / grand_total uses the FULL-matrix totals, so blocking over
    genes is exact."""
    n, g = X.shape
    sc_ = np.asarray(X.sum(1)).ravel().astype(np.float64)
    sg = np.asarray(X.sum(0)).ravel().astype(np.float64)
    tot = sg.sum()
    clip = np.sqrt(n) if clip is None else clip
    Xc = X.tocsc()
    for j in range(0, g, block):
        D = Xc[:, j:j + block].toarray().astype(np.float64)
        mu = np.outer(sc_, sg[j:j + block]) / tot
        with np.errstate(invalid="ignore", divide="ignore"):
            R = (D - mu) / np.sqrt(mu + mu ** 2 / theta)
        R[~np.isfinite(R)] = 0.0            # genes with zero total (mu = 0)
        np.clip(R, -clip, clip, out=R)
        yield j, R


def m_pearson(X, lab, types, theta=100.0, block=1000):
    masks = [(lab == t) for t in types]
    out = np.zeros((len(types), X.shape[1]), dtype=np.float64)
    for j, R in pearson_residual_blocks(X, theta=theta, block=block):
        for i, m in enumerate(masks):
            out[i, j:j + R.shape[1]] = R[m].mean(0)
    return out


def check_pearson_against_scanpy(X, n_cells=2000, n_genes=1500, theta=100.0):
    """Self-check: my block implementation vs scanpy on a dense-able sub-matrix."""
    import scanpy as sc
    sub = X[:n_cells].tocsc()
    keep = np.flatnonzero(np.asarray(sub.sum(0)).ravel() > 0)[:n_genes]
    sub = sub[:, keep].tocsr()
    a = ad.AnnData(sub.copy())
    sc.experimental.pp.normalize_pearson_residuals(a, theta=theta, clip=None, check_values=False)
    ref = np.asarray(a.X, dtype=np.float64)
    mine = np.zeros_like(ref)
    for j, R in pearson_residual_blocks(sub, theta=theta, block=400):
        mine[:, j:j + R.shape[1]] = R
    return float(np.abs(mine - ref).max()), float(np.abs(ref).max())


def rank_within_cell(X, normalise=False):
    d = X.data.copy(); ip = X.indptr
    Xr = sp.csr_matrix((d, X.indices, ip), shape=X.shape)   # shares X's index arrays
    for i in range(Xr.shape[0]):
        s, e = ip[i], ip[i + 1]
        if e > s:
            r = rankdata(d[s:e]).astype(np.float32)
            d[s:e] = r / (e - s) if normalise else r
    return Xr


def m_detection(X, lab, types):
    Xb = sp.csr_matrix(((X.data > 0).astype(np.float32), X.indices, X.indptr), shape=X.shape)
    return B.type_mean(Xb, lab, types)


def run_matrix(i):
    ctx = load_matrix(i)
    X, lab, types = ctx["X"], ctx["lab"], ctx["types"]
    rng = np.random.default_rng(SEED + i)
    rows, mrows, argm = [], [], {}

    def do(name, fn):
        t0 = time.time(); res = fn()
        P, extra = (res if isinstance(res, tuple) else (res, None))
        row, mr, arg = evaluate(P, ctx, rng, name, time.time() - t0, extra)
        rows.append(row); mrows.extend(mr); argm[name] = arg

    # original six (recomputed; identity checks only)
    do("raw_mean", lambda: B.m_raw(X, lab, types))
    do("zscale_mean", lambda: B.m_zscale(X, lab, types))
    do("pseudobulk_cpm", lambda: B.m_pseudobulk_cpm(X, lab, types))
    do("depth_matched", lambda: B.m_depth_matched(X, lab, types, rng))
    do("cpm_mean", lambda: B.m_cpm(X, lab, types))
    do("lognorm_mean", lambda: B.m_lognorm(X, lab, types))
    # new normalisers
    do("tfidf1_mean", lambda: B.type_mean(tfidf_matrix(X, 1), lab, types))
    do("tfidf2_mean", lambda: B.type_mean(tfidf_matrix(X, 2), lab, types))
    do("tfidf3_mean", lambda: B.type_mean(tfidf_matrix(X, 3), lab, types))
    err, scale = check_pearson_against_scanpy(X)
    log(f"  pearson self-check vs scanpy on {2000}x{1500} sub-matrix: max|diff| = {err:.2e} (max|resid| {scale:.1f})")
    do("pearson_mean", lambda: m_pearson(X, lab, types))
    med = float(np.median(np.asarray(X.sum(1)).ravel()))
    do("medscale_mean", lambda: (B.type_mean(B.percell_scale(X, target=med), lab, types), med))
    do("rank_mean", lambda: B.type_mean(rank_within_cell(X, False), lab, types))
    do("rankfrac_mean", lambda: B.type_mean(rank_within_cell(X, True), lab, types))
    do("detection_mean", lambda: m_detection(X, lab, types))

    def lognorm_maxnorm():
        P = B.m_lognorm(X, lab, types); mx = P.max(0); mx[mx == 0] = 1.0
        return P / mx
    do("lognorm_maxnorm", lognorm_maxnorm)
    do("pseudobulk_logcpm", lambda: np.log1p(B.m_pseudobulk_cpm(X, lab, types)))

    tag = f"{ctx['cohort']}_{ctx['model']}"
    pd.DataFrame(rows)[COLS].to_csv(PARTS / f"extra_{tag}.csv", index=False)
    pd.DataFrame(mrows).to_csv(PARTS / f"extra_markers_{tag}.csv", index=False)
    np.savez_compressed(PARTS / f"extra_argmax_{tag}.npz", expressed=ctx["expressed"], types=np.array(types), **argm)
    with open(PARTS / f"extra_{tag}.pearson_check.txt", "w") as f:
        f.write(f"max_abs_diff_vs_scanpy\t{err:.3e}\nmax_abs_residual\t{scale:.3f}\n")
    log(f"done {tag}")


ORIG = ["raw_mean", "zscale_mean", "pseudobulk_cpm", "depth_matched", "cpm_mean", "lognorm_mean"]
NEW = ["tfidf1_mean", "tfidf2_mean", "tfidf3_mean", "pearson_mean", "medscale_mean", "rank_mean",
       "rankfrac_mean", "detection_mean", "lognorm_maxnorm", "pseudobulk_logcpm"]


def merge():
    R = pd.concat([pd.read_csv(p) for p in sorted(PARTS.glob("extra_*_*.csv")) if "markers" not in p.name])
    M = pd.concat([pd.read_csv(p) for p in sorted(PARTS.glob("extra_markers_*.csv"))])
    order = {f"{c}_{m}": k for k, (c, m, _, _) in enumerate(REAL)}
    R["_o"] = (R.cohort + "_" + R.model).map(order); R = R.sort_values(["_o", "method"]).drop(columns="_o")
    R["method"] = pd.Categorical(R.method, ORIG + NEW, ordered=True)
    R = R.sort_values(["cohort", "model", "method"], key=lambda s: s.map(order) if s.name == "cohort" else s)
    # identity of the per-gene argmax vectors (over expressed genes) with the original methods
    ident = {}
    for c, m, _, _ in REAL:
        if not (PARTS / f"extra_argmax_{c}_{m}.npz").exists():
            print(f"WARNING: {c}/{m} part missing, merged output is partial"); continue
        z = np.load(PARTS / f"extra_argmax_{c}_{m}.npz"); ex = z["expressed"]
        for name in NEW:
            same = [o for o in ORIG if np.array_equal(z[name][ex], z[o][ex])]
            frac = {o: float(np.mean(z[name][ex] == z[o][ex])) for o in ORIG}
            best = max(frac, key=frac.get)
            ident[(c, m, name)] = (same[0] if same else "none", best, frac[best])
    R["argmax_identical_to"] = [ident.get((r.cohort, r.model, r.method), ("", "", np.nan))[0] if r.method in NEW else ""
                                for r in R.itertuples()]
    R["closest_original"] = [f"{ident[(r.cohort, r.model, r.method)][1]} ({ident[(r.cohort, r.model, r.method)][2]:.3f})"
                             if r.method in NEW else "" for r in R.itertuples()]
    R.to_csv(HERE / "benchmark_all_methods.csv", index=False)
    R[R.method.isin(NEW)][COLS].to_csv(HERE / "benchmark_extra_normalisers.csv", index=False)
    M.to_csv(HERE / "benchmark_extra_markers.csv", index=False)
    print("\n=== summary across the six real matrices ===")
    for name in ORIG + NEW:
        s = R[R.method == name]
        print(f"  {name:18} top share median {s.top_share.median():5.1%} ({s.top_share.min():.0%}-{s.top_share.max():.0%})"
              f"   argmax in deepest {(s.top_type == s.deepest_type).sum()}/{len(s)}"
              f"   rho median {s.rho_depth_share.median():+.2f}, p<0.05 {(s.p_depth_share < 0.05).sum()}/{len(s)}"
              f"   markers {s.markers_correct.sum()}/{s.markers_total.sum()}"
              f"   perm p<0.05 {(s.perm_p < 0.05).sum()}/{len(s)}"
              + (f"   argmax identical to: {sorted(set(s.argmax_identical_to))}" if name in NEW else ""))
    for p in sorted(PARTS.glob("extra_*.pearson_check.txt")):
        print(p.name, open(p).read().replace("\n", "  "))
    print(f"wrote {HERE/'benchmark_extra_normalisers.csv'}, {HERE/'benchmark_all_methods.csv'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=int, default=None)
    ap.add_argument("--merge", action="store_true")
    args = ap.parse_args()
    if args.merge:
        merge()
    elif args.matrix is not None:
        run_matrix(args.matrix)
    else:
        for i in range(len(REAL)):
            run_matrix(i)
        merge()
