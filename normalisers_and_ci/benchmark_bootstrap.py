"""
Cell-level bootstrap intervals and paired tests for the six aggregation methods of
tools/benchmark_aggregation.py (raw_mean, zscale_mean, pseudobulk_cpm, depth_matched,
cpm_mean, lognorm_mean; archr_style is log2 instead of ln and has the same argmax as
lognorm_mean, so it is not repeated).

Bootstrap (per matrix, B = 200): cells are resampled WITH replacement WITHIN each cell
type (stratified, type sizes fixed), and every method is re-aggregated on the same
resample (paired across methods). A resample is a weight vector w (how many times each
cell was drawn), and every method here is a weighted mean of per-cell rows:

  raw_mean        (W X) / n_t
  pseudobulk_cpm  row-CPM of (W X)
  zscale_mean     ((W X)/n_t - mu)/sd - (W E)/n_t, where E holds the clip excess
                  max(0, (x-mu)/sd - 10) of the non-zero entries (exact identity for
                  scanpy.pp.scale(max_value=10); mu, sd per gene are held at their
                  full-data values, n >= 12k cells)
  cpm_mean        (W X_cpm)/n_t
  lognorm_mean    (W log1p X_cpm)/n_t
  depth_matched   (W X_thin)/n_t, X_thin = ONE binomial thinning to the shallowest type's
                  median depth (the bootstrap is conditional on that thinning; the extra
                  thinning noise is measured separately from 10 independent thinnings
                  without resampling and reported as thinning_sd)

All five data variants share one CSR index structure (explicit zeros kept after thinning),
so the largest matrix (LUAD enhancer, 488 M non-zeros) needs ~12 GB.

Read-outs per replicate and method: top argmax share over expressed genes, its type, and
the held-out marker count (same panels as the benchmark). Percentile 95 % intervals are
taken over the 200 replicates; paired differences vs cpm_mean use the same replicates.

Across the six matrices: paired Wilcoxon signed-rank tests (exact, n = 6; smallest
attainable two-sided p = 2/64 = 0.031) on the point-estimate top share and marker
fraction, method vs cpm_mean, plus the exact sign test. New normalisers from
benchmark_extra_normalisers.csv are tested the same way when that file exists.

Run:  python3 benchmark_bootstrap.py --matrix i    (i = 0..5)
      python3 benchmark_bootstrap.py --merge
Out:  benchmark_bootstrap_ci.csv        one row per matrix x method
      benchmark_bootstrap_replicates.npz
      benchmark_paired_tests.csv        Wilcoxon / sign tests across the six matrices
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd
import scipy.sparse as sp
from scipy.stats import wilcoxon, binomtest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from benchmark_extra_normalisers import load_matrix, log, REAL, PARTS, ORIG, NEW  # noqa: E402
import benchmark_aggregation as B  # noqa: E402

N_BOOT = 200
N_THIN = 10
SEED = 42
METHODS = ORIG   # raw_mean, zscale_mean, pseudobulk_cpm, depth_matched, cpm_mean, lognorm_mean


CHUNK = 20_000_000   # entries per chunk for the element-wise passes over X.data


def _row_of_entry(X):
    return np.repeat(np.arange(X.shape[0], dtype=np.int32), np.diff(X.indptr))


def thin_data(X, lab, types, rng, row=None):
    """Binomial thinning exactly as B.m_depth_matched, but returning the data array only
    (explicit zeros kept so the CSR structure stays shared); chunked to bound memory."""
    tot = np.asarray(X.sum(1)).ravel().astype(np.float64)
    target = min(np.median(tot[lab == t]) for t in types)
    p = np.minimum(1.0, target / np.maximum(tot, 1.0))
    row = _row_of_entry(X) if row is None else row
    d = np.empty(X.nnz, dtype=np.float32)
    for s in range(0, X.nnz, CHUNK):
        e = min(s + CHUNK, X.nnz)
        d[s:e] = rng.binomial(np.rint(X.data[s:e]).astype(np.int64), p[row[s:e]])
    return d, float(target)


def build_variants(X, lab, types, rng):
    n = X.shape[0]
    tot = np.asarray(X.sum(1)).ravel().astype(np.float64); tot[tot == 0] = 1.0
    inv = (1e4 / tot).astype(np.float32)
    row = _row_of_entry(X)
    d_cpm = np.empty(X.nnz, dtype=np.float32); d_log = np.empty(X.nnz, dtype=np.float32)
    d_exc = np.empty(X.nnz, dtype=np.float32)
    # per-gene mean and sd (ddof=1) of the raw matrix, as scanpy.pp.scale / B.m_zscale
    mu = np.asarray(X.sum(0)).ravel().astype(np.float64) / n
    sq = np.zeros(X.shape[1], dtype=np.float64)
    for s in range(0, X.nnz, CHUNK):
        e = min(s + CHUNK, X.nnz)
        sq += np.bincount(X.indices[s:e], weights=X.data[s:e].astype(np.float64) ** 2, minlength=X.shape[1])
    var = (sq - n * mu ** 2) / (n - 1); var[var < 0] = 0
    sd = np.sqrt(var); sd[sd == 0] = 1.0
    for s in range(0, X.nnz, CHUNK):
        e = min(s + CHUNK, X.nnz)
        d_cpm[s:e] = X.data[s:e] * inv[row[s:e]]
        d_log[s:e] = np.log1p(d_cpm[s:e])
        z = (X.data[s:e].astype(np.float64) - mu[X.indices[s:e]]) / sd[X.indices[s:e]] - 10.0
        d_exc[s:e] = np.maximum(z, 0.0)
    d_thin, target = thin_data(X, lab, types, rng, row)
    return dict(raw=X.data, cpm=d_cpm, log=d_log, exc=d_exc, thin=d_thin), mu, sd, target


def weighted_sums(X, data, W):
    """W (types x cells, dense float32) @ X with X.data replaced by `data` -> types x genes."""
    M = sp.csr_matrix((data, X.indices, X.indptr), shape=X.shape)
    return np.asarray((M.T @ W.T).T, dtype=np.float64)


def profiles(X, V, mu, sd, W):
    n_t = W.sum(1)[:, None]
    S = weighted_sums(X, V["raw"], W)
    P = {"raw_mean": S / n_t,
         "pseudobulk_cpm": S / S.sum(1, keepdims=True) * 1e6,
         "zscale_mean": (S / n_t - mu) / sd - weighted_sums(X, V["exc"], W) / n_t,
         "cpm_mean": weighted_sums(X, V["cpm"], W) / n_t,
         "lognorm_mean": weighted_sums(X, V["log"], W) / n_t,
         "depth_matched": weighted_sums(X, V["thin"], W) / n_t}
    return P


def readouts(P, ctx):
    ex, types, panel, gidx = ctx["expressed"], ctx["types"], ctx["panel"], ctx["gidx"]
    pidx = np.array([gidx[g] for g in panel]); ptgt = np.array([types.index(t) for t in panel.values()])
    out = {}
    for m, Pm in P.items():
        arg = Pm[:, ex].argmax(0)
        share = np.bincount(arg, minlength=len(types)) / arg.size
        top = int(share.argmax())
        mk = int((Pm[:, pidx].argmax(0) == ptgt).sum())
        out[m] = (float(share[top]), top, mk)
    return out


def run_matrix(i):
    ctx = load_matrix(i)
    X, lab, types = ctx["X"], ctx["lab"], ctx["types"]
    n = X.shape[0]; rng = np.random.default_rng(SEED + i)
    t0 = time.time()
    V, mu, sd, target = build_variants(X, lab, types, rng)
    log(f"  variants built ({time.time()-t0:.0f}s), depth-matched target {target:.0f}")
    masks = [np.flatnonzero(lab == t) for t in types]
    W_full = np.zeros((len(types), n), dtype=np.float32)
    for k, m in enumerate(masks): W_full[k, m] = 1.0
    point = readouts(profiles(X, V, mu, sd, W_full), ctx)
    for m in METHODS:
        log(f"  point {m:15} top {types[point[m][1]]:>20} {point[m][0]:5.1%}  markers {point[m][2]}/{len(ctx['panel'])}")
    # thinning-only noise: N_THIN independent thinnings, no resampling
    thin_share, thin_mk = [], []
    for r in range(N_THIN):
        d, _ = thin_data(X, lab, types, rng)
        Pm = weighted_sums(X, d, W_full) / W_full.sum(1)[:, None]
        s, _, k = readouts({"depth_matched": Pm}, ctx)["depth_matched"]
        thin_share.append(s); thin_mk.append(k)
    log(f"  thinning-only: top share sd {np.std(thin_share, ddof=1):.4f}, markers {np.min(thin_mk)}-{np.max(thin_mk)}")
    # stratified cell bootstrap
    share = np.zeros((N_BOOT, len(METHODS))); topt = np.zeros((N_BOOT, len(METHODS)), dtype=np.int16)
    mk = np.zeros((N_BOOT, len(METHODS)), dtype=np.int16)
    t0 = time.time()
    for b in range(N_BOOT):
        W = np.zeros((len(types), n), dtype=np.float32)
        for k, m in enumerate(masks):
            draw = rng.integers(0, m.size, m.size)
            W[k, m] = np.bincount(draw, minlength=m.size)
        R = readouts(profiles(X, V, mu, sd, W), ctx)
        for j, mth in enumerate(METHODS):
            share[b, j], topt[b, j], mk[b, j] = R[mth]
        if b % 20 == 19:
            log(f"  bootstrap {b+1}/{N_BOOT} ({time.time()-t0:.0f}s)")
    tag = f"{ctx['cohort']}_{ctx['model']}"
    np.savez_compressed(PARTS / f"boot_{tag}.npz", share=share, top=topt, markers=mk, methods=np.array(METHODS),
                        types=np.array(types), point_share=np.array([point[m][0] for m in METHODS]),
                        point_top=np.array([point[m][1] for m in METHODS]),
                        point_markers=np.array([point[m][2] for m in METHODS]),
                        thin_share=np.array(thin_share), thin_markers=np.array(thin_mk),
                        n_markers=len(ctx["panel"]), deepest=int(np.nanargmax(ctx["med_depth"])), target=target)
    log(f"done {tag}")


def merge():
    rows = []
    for c, m, _, _ in REAL:
        if not (PARTS / f"boot_{c}_{m}.npz").exists():
            print(f"WARNING: {c}/{m} bootstrap part missing, merged output is partial"); continue
        z = np.load(PARTS / f"boot_{c}_{m}.npz")
        types = list(z["types"]); meth = list(z["methods"]); j_cpm = meth.index("cpm_mean")
        for j, name in enumerate(meth):
            s, t, k = z["share"][:, j], z["top"][:, j], z["markers"][:, j]
            d = s - z["share"][:, j_cpm]; dk = k.astype(int) - z["markers"][:, j_cpm].astype(int)
            top_counts = np.bincount(t, minlength=len(types))
            ps, pk = float(z["point_share"][j]), int(z["point_markers"][j])
            rows.append(dict(cohort=c, model=m, method=name, n_boot=len(s),
                             top_share=ps, top_share_ci_lo=float(np.percentile(s, 2.5)),
                             top_share_ci_hi=float(np.percentile(s, 97.5)), top_share_boot_sd=float(s.std(ddof=1)),
                             top_share_boot_mean=float(s.mean()), top_share_boot_bias=float(s.mean() - ps),
                             top_share_basic_lo=float(2 * ps - np.percentile(s, 97.5)),
                             top_share_basic_hi=float(2 * ps - np.percentile(s, 2.5)),
                             top_type=types[int(z["point_top"][j])],
                             top_type_mode=types[int(top_counts.argmax())],
                             top_type_stability=float(top_counts.max() / len(t)),
                             frac_top_is_deepest=float(np.mean(t == int(z["deepest"]))),
                             markers_correct=int(z["point_markers"][j]), markers_total=int(z["n_markers"]),
                             markers_ci_lo=float(np.percentile(k, 2.5)), markers_ci_hi=float(np.percentile(k, 97.5)),
                             markers_boot_mean=float(k.mean()),
                             markers_basic_lo=float(2 * pk - np.percentile(k, 97.5)), markers_basic_hi=float(2 * pk - np.percentile(k, 2.5)),
                             diff_share_vs_cpm=float(z["point_share"][j] - z["point_share"][j_cpm]),
                             diff_share_ci_lo=float(np.percentile(d, 2.5)), diff_share_ci_hi=float(np.percentile(d, 97.5)),
                             diff_markers_vs_cpm=int(z["point_markers"][j] - z["point_markers"][j_cpm]),
                             diff_markers_ci_lo=float(np.percentile(dk, 2.5)), diff_markers_ci_hi=float(np.percentile(dk, 97.5)),
                             thinning_sd=(float(np.std(z["thin_share"], ddof=1)) if name == "depth_matched" else np.nan),
                             thinning_markers_range=(f"{z['thin_markers'].min()}-{z['thin_markers'].max()}" if name == "depth_matched" else ""),
                             depth_matched_target=(float(z["target"]) if name == "depth_matched" else np.nan)))
    R = pd.DataFrame(rows); R.to_csv(HERE / "benchmark_bootstrap_ci.csv", index=False)
    # paired tests across the six matrices, on the point estimates
    pub = pd.read_csv(HERE.parents[2] / "results" / "benchmark_aggregation.csv")
    pub = pub[pub.model != "misannotated"]
    src = [("published benchmark_aggregation.csv", pub)]
    fx = HERE / "benchmark_extra_normalisers.csv"
    if fx.exists():
        src.append(("benchmark_extra_normalisers.csv", pd.read_csv(fx)))
    trows = []
    key = ["cohort", "model"]
    ref = pub[pub.method == "cpm_mean"].set_index(key)
    for label, df in src:
        for name in sorted(set(df.method)):
            if name in ("cpm_mean", "archr_style"): continue
            a = df[df.method == name].set_index(key)
            common = ref.index.intersection(a.index); a = a.loc[common]; ref_ = ref.loc[common]
            if len(common) < 2: continue
            for metric, xa, xb in [("top_share", a.top_share.values, ref_.top_share.values),
                                   ("marker_fraction", (a.markers_correct / a.markers_total).values,
                                    (ref_.markers_correct / ref_.markers_total).values)]:
                diff = xa - xb; nz = diff != 0
                if nz.sum() == 0:
                    p2 = p_gt = p_sign = 1.0
                else:
                    p2 = float(wilcoxon(xa, xb, alternative="two-sided", method="exact").pvalue)
                    p_gt = float(wilcoxon(xa, xb, alternative="greater", method="exact").pvalue)
                    p_sign = float(binomtest(int((diff > 0).sum()), int(nz.sum()), 0.5).pvalue)
                trows.append(dict(source=label, method=name, vs="cpm_mean", metric=metric, n_pairs=len(diff),
                                  n_nonzero_diffs=int(nz.sum()), n_positive=int((diff > 0).sum()),
                                  median_diff=float(np.median(diff)), mean_diff=float(diff.mean()),
                                  wilcoxon_p_two_sided=p2, wilcoxon_p_greater=p_gt, sign_test_p=p_sign,
                                  min_attainable_two_sided_p=2 / 2 ** int(nz.sum()) if nz.sum() else 1.0))
    T = pd.DataFrame(trows); T.to_csv(HERE / "benchmark_paired_tests.csv", index=False)
    print(R[["cohort", "model", "method", "top_share", "top_share_ci_lo", "top_share_ci_hi", "top_type", "top_type_stability",
             "markers_correct", "markers_ci_lo", "markers_ci_hi"]].to_string(index=False))
    print(T[T.metric == "top_share"][["source", "method", "n_positive", "median_diff", "wilcoxon_p_two_sided",
                                      "wilcoxon_p_greater", "sign_test_p"]].to_string(index=False))
    print(f"wrote {HERE/'benchmark_bootstrap_ci.csv'}, {HERE/'benchmark_paired_tests.csv'}")


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
