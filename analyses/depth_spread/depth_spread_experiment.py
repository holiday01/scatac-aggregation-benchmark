"""
Controlled depth-spread experiment for the aggregation benchmark.

The benchmark (tools/benchmark_aggregation.py) shows the effect of the aggregation
method at the three NATURAL depth spreads of the cohorts (HCC ~2x, LUAD ~3.5x,
COAD ~37x; spread = median per-cell depth of the deepest type / shallowest type).
Here the spread is IMPOSED on two count matrices so the whole curve is measured:

  1. equalise: every cell is binomially thinned to a common target depth
     (the 20th percentile of matrix row sums; cells already below it are kept
     as they are), so every cell type starts at spread ~1x;
  2. impose:   one cell type ("deep") is kept at the target, the others are thinned
     further by factors spaced evenly on a log scale between 1 and 1/spread
     (their order is a random permutation per seed), for spread = 1, 2, 5, 10,
     20, 40x, two choices of the deep type, three seeds.  Both thinnings compose
     into one binomial draw from the original counts with p = min(1, target/rowsum)
     * factor(type).
  3. score:    five aggregation methods (raw_mean, zscale_mean, lognorm_mean,
     pseudobulk_cpm, cpm_mean; definitions identical to the benchmark script) ->
     top argmax share, winning type, whether the deep type wins, Spearman of
     per-type share vs realised median depth and vs detection rate, held-out
     marker concordance (PANEL of the benchmark script) with a permutation null;
  4. mechanism: per type and condition, the detection rate, the mean over
     expressed genes of the log1p(CPM) profile, of the CPM profile and of the raw
     profile, and the argmax share under each method.

zscale_mean is computed with a sparse-exact shortcut: scanpy.pp.scale(max_value=10)
clips only entries with (x - mu)/sd > 10, and x = 0 never reaches -10 (min over
genes of -mu/sd is -1.75 on HCC), so clipping x at mu + 10 sd keeps sparsity and the
per-type mean of the clipped z-score is (mean of clipped x - mu)/sd.  Checked
against benchmark_aggregation.m_zscale on a 3000 x 6000 HCC block: max |diff| 4e-5.

Usage:
  python3 depth_spread_experiment.py run HCC        # -> depth_spread_results_HCC.csv etc.
  python3 depth_spread_experiment.py run LUAD
  python3 depth_spread_experiment.py summarise      # -> merged CSVs, summary, figure
Outputs are written next to this script.  No file outside this directory is written.
"""
import os, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd, hdf5plugin  # noqa: F401  (plugin needed to read the files)
import scipy.sparse as sp
import anndata as ad
from scipy.stats import spearmanr

sys.path.insert(0, "/mnt/10t/scrna_atac/apbc2026/tools")
from benchmark_aggregation import PANEL, percell_scale  # noqa: E402

HERE = Path(__file__).resolve().parent
BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "/mnt/10t/scrna_atac"))
COHORTS = {
    "HCC":  dict(path=BASE / "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad",
                 deep=["B_cell", "NK_cytotoxic_T"]),
    "LUAD": dict(path=BASE / "data/processed/LUAD/luad_geneactivity_proximal.h5ad",
                 deep=["B_cell", "T_NK"]),
}
SPREADS = [1, 2, 5, 10, 20, 40]
SEEDS = [0, 1, 2]
METHODS = ["raw_mean", "zscale_mean", "lognorm_mean", "pseudobulk_cpm", "cpm_mean"]
TARGET_PCT = 20          # common target depth = this percentile of matrix row sums
N_PERM = 2000
SEED_BASE = 20260827


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ── aggregation via a one-hot indicator matrix (one pass over the matrix per method) ──
def indicator(lab, types):
    n = len(lab)
    row = np.searchsorted(types, lab)
    return sp.csr_matrix((np.ones(n, np.float32), (row, np.arange(n))), shape=(len(types), n))


def type_sum(M, X):
    return np.asarray((M @ X).todense(), dtype=np.float64)


def profiles(X, lab, types):
    """All five per-type profiles (types x genes) from one csr float32 count matrix."""
    M = indicator(lab, types)
    cnt = np.asarray(M.sum(1)).ravel()
    n = X.shape[0]
    out = {}
    S = type_sum(M, X)
    out["raw_mean"] = S / cnt[:, None]
    out["pseudobulk_cpm"] = S / np.maximum(S.sum(1, keepdims=True), 1.0) * 1e6
    # z-scale (scanpy.pp.scale zero_center, max_value=10), sparse-exact
    mu = S.sum(0) / n
    ex2 = np.asarray(X.multiply(X).sum(0)).ravel() / n
    sd = np.sqrt(np.maximum(ex2 - mu ** 2, 0.0) * n / (n - 1)); sd[sd == 0] = 1.0
    thr = (mu + 10.0 * sd).astype(np.float32)
    Xc = X.copy(); Xc.data = np.minimum(Xc.data, thr[Xc.indices])
    out["zscale_mean"] = (type_sum(M, Xc) / cnt[:, None] - mu) / sd
    del Xc
    Xn = percell_scale(X).tocsr()
    out["cpm_mean"] = type_sum(M, Xn) / cnt[:, None]
    Xn.data = np.log1p(Xn.data)
    out["lognorm_mean"] = type_sum(M, Xn) / cnt[:, None]
    del Xn
    return out


def thin(X, p, rng):
    """Binomial thinning of every count with a per-cell keep probability p."""
    Xd = X.copy()
    per_nz = np.repeat(p, np.diff(Xd.indptr))
    Xd.data = rng.binomial(np.rint(Xd.data).astype(np.int64), per_nz).astype(np.float32)
    Xd.eliminate_zeros()
    return Xd


def perm_p(argmax, target, obs, rng, n=N_PERM):
    null = np.array([(rng.permutation(argmax) == target).sum() for _ in range(n)])
    return float(np.mean(null >= obs)), float(null.mean())


def score(X, lab, types, cohort, panel, gidx, rng, cond, trows, rows, mrows):
    """Score one thinned matrix: per-type rows, per-method rows, per-marker rows."""
    k = len(types)
    rs = np.asarray(X.sum(1)).ravel()
    nnz = np.diff(X.indptr) / X.shape[1]
    med = np.array([np.median(rs[lab == t]) for t in types])
    det = np.array([float(nnz[lab == t].mean()) for t in types])
    expressed = np.asarray(X.sum(0)).ravel() > 0
    P = profiles(X, lab, types)
    shares = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for m in METHODS:
            arg = P[m][:, expressed].argmax(0)
            share = np.bincount(arg, minlength=k) / arg.size
            shares[m] = share
            top = int(share.argmax())
            rho, pv = (spearmanr(med, share) if cond["spread"] != 1 else (np.nan, np.nan))
            rho_det, _ = spearmanr(det, share)
            marg = np.array([types[int(P[m][:, gidx[g]].argmax())] for g in panel])
            mtgt = np.array(list(panel.values()))
            correct = marg == mtgt
            obs = int(correct.sum())
            pp, nmean = perm_p(marg, mtgt, obs, rng)
            per_lin = {t: f"{int(correct[mtgt == t].sum())}/{int((mtgt == t).sum())}" for t in types if (mtgt == t).any()}
            rows.append(dict(**cond, method=m, n_cells=X.shape[0], n_genes_expressed=int(expressed.sum()),
                             realised_spread=float(med.max() / med.min()),
                             top_share=float(share[top]), top_type=types[top],
                             deep_wins=int(types[top] == cond["deep_type"]) if cond["deep_type"] != "none" else np.nan,
                             deep_type_share=float(share[types.index(cond["deep_type"])]) if cond["deep_type"] != "none" else np.nan,
                             shallowest_type=types[int(med.argmin())],
                             rho_depth_share=float(rho), p_depth_share=float(pv), rho_detection_share=float(rho_det),
                             markers_correct=obs, markers_total=len(panel), perm_p=pp, null_mean=nmean,
                             per_lineage="; ".join(f"{a} {b}" for a, b in per_lin.items())))
            for g, tg, ag, c in zip(panel, mtgt, marg, correct):
                mrows.append(dict(**cond, method=m, gene=g, target=tg, argmax=ag, correct=int(c)))
    for i, t in enumerate(types):
        trows.append(dict(**cond, cell_type=t, n_cells=int((lab == t).sum()),
                          imposed_factor=float(cond["_factors"][t]) if cond["_factors"] else np.nan,
                          realised_median_rowsum=float(med[i]), detection_rate=float(det[i]),
                          mean_lognorm_profile=float(P["lognorm_mean"][i, expressed].mean()),
                          mean_cpm_profile=float(P["cpm_mean"][i, expressed].mean()),
                          mean_raw_profile=float(P["raw_mean"][i, expressed].mean()),
                          mean_zscale_profile=float(P["zscale_mean"][i, expressed].mean()),
                          **{f"share_{m}": float(shares[m][i]) for m in METHODS}))
    return {m: types[int(shares[m].argmax())] for m in METHODS}, {m: float(shares[m].max()) for m in METHODS}


def run(cohort):
    cfg = COHORTS[cohort]
    log(f"{cohort}: loading {cfg['path']}")
    a = ad.read_h5ad(cfg["path"])
    lab = a.obs["cell_type"].astype(str).values
    keep = lab != "nan"
    a = a[keep].copy(); lab = a.obs["cell_type"].astype(str).values
    if cohort == "HCC":   # relabel rule of the benchmark (tools/check_hcc_labels.py)
        lab = np.where(lab == "B_cell", "Endothelial_stromal", np.where(lab == "DC", "B_cell", lab))
    X = a.X.tocsr().astype(np.float32); X.sort_indices()
    genes = a.var_names.astype(str).values
    del a
    types = sorted(pd.unique(lab)); k = len(types)
    panel = {g: t for g, t in PANEL[cohort].items() if g in set(genes) and t in types}
    gidx = {g: i for i, g in enumerate(genes)}
    rs0 = np.asarray(X.sum(1)).ravel().astype(np.float64)
    target = float(np.percentile(rs0, TARGET_PCT))
    p_base = np.minimum(1.0, target / np.maximum(rs0, 1.0))
    nat = {t: float(np.median(rs0[lab == t])) for t in types}
    log(f"  {X.shape[0]} cells x {X.shape[1]} genes, {k} types; target depth (p{TARGET_PCT}) = {target:.0f}; "
        f"natural spread {max(nat.values())/min(nat.values()):.2f}x; panel {len(panel)} markers")
    rows, trows, mrows = [], [], []
    rng_nat = np.random.default_rng([SEED_BASE, 999])
    base = dict(cohort=cohort, target_depth=target)
    # natural (unthinned) reference condition
    cond = dict(**base, spread="natural", deep_type="none", seed=-1, _factors=None)
    t0 = time.time()
    win, sh = score(X, lab, types, cohort, panel, gidx, rng_nat, cond, trows, rows, mrows)
    log(f"  natural           " + "  ".join(f"{m} {sh[m]:.0%}({win[m]})" for m in METHODS) + f"  [{time.time()-t0:.0f}s]")
    for di, deep in enumerate(cfg["deep"]):
        others = [t for t in types if t != deep]
        for seed in SEEDS:
            for si, S in enumerate(SPREADS):
                rng = np.random.default_rng([SEED_BASE, list(COHORTS).index(cohort), di, seed, si])
                order = list(rng.permutation(others))
                factors = {deep: 1.0}
                for i, t in enumerate(order, start=1):
                    factors[t] = float(S ** (-i / (k - 1)))
                f_cell = np.array([factors[t] for t in lab])
                t0 = time.time()
                Xd = thin(X, p_base * f_cell, rng)
                cond = dict(**base, spread=S, deep_type=deep, seed=seed, _factors=factors)
                win, sh = score(Xd, lab, types, cohort, panel, gidx, rng, cond, trows, rows, mrows)
                del Xd
                log(f"  deep={deep:>14} seed={seed} spread={S:>2}x  " +
                    "  ".join(f"{m.split('_')[0]} {sh[m]:.0%}({'D' if win[m]==deep else win[m][:6]})" for m in METHODS) +
                    f"  [{time.time()-t0:.0f}s]")
    R = pd.DataFrame(rows).drop(columns=["_factors"])
    T = pd.DataFrame(trows).drop(columns=["_factors"])
    Mk = pd.DataFrame(mrows).drop(columns=["_factors"])
    R.to_csv(HERE / f"depth_spread_results_{cohort}.csv", index=False)
    T.to_csv(HERE / f"depth_spread_by_type_{cohort}.csv", index=False)
    Mk.to_csv(HERE / f"depth_spread_markers_{cohort}.csv", index=False)
    log(f"{cohort}: wrote {len(R)} method rows, {len(T)} type rows, {len(Mk)} marker rows")


# ── summary + figure ─────────────────────────────────────────────────────────
LABEL = {"raw_mean": "raw mean", "zscale_mean": "z-scale mean", "lognorm_mean": "log1p(CPM) mean",
         "pseudobulk_cpm": "pseudobulk CPM", "cpm_mean": "CPM mean"}
COLOR = {"raw_mean": "#2a78d6", "zscale_mean": "#eb6834", "lognorm_mean": "#1baf7a",
         "pseudobulk_cpm": "#eda100", "cpm_mean": "#e87ba4"}
MARK = {"raw_mean": "o", "zscale_mean": "s", "lognorm_mean": "^", "pseudobulk_cpm": "D", "cpm_mean": "v"}


def summarise():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R = pd.concat([pd.read_csv(HERE / f"depth_spread_results_{c}.csv") for c in COHORTS], ignore_index=True)
    T = pd.concat([pd.read_csv(HERE / f"depth_spread_by_type_{c}.csv") for c in COHORTS], ignore_index=True)
    Mk = pd.concat([pd.read_csv(HERE / f"depth_spread_markers_{c}.csv") for c in COHORTS], ignore_index=True)
    R.to_csv(HERE / "depth_spread_results.csv", index=False)
    T.to_csv(HERE / "depth_spread_by_type.csv", index=False)
    Mk.to_csv(HERE / "depth_spread_markers.csv", index=False)
    A = R[R.spread != "natural"].copy(); A["spread"] = A.spread.astype(int)
    g = A.groupby(["cohort", "spread", "method"])
    S = g.agg(n=("top_share", "size"), top_share_mean=("top_share", "mean"), top_share_min=("top_share", "min"),
              top_share_max=("top_share", "max"), deep_wins=("deep_wins", "sum"),
              deep_type_share_mean=("deep_type_share", "mean"),
              rho_depth_share_mean=("rho_depth_share", "mean"), rho_detection_share_mean=("rho_detection_share", "mean"),
              markers_correct_mean=("markers_correct", "mean"), markers_correct_min=("markers_correct", "min"),
              markers_correct_max=("markers_correct", "max"), markers_total=("markers_total", "first"),
              perm_p_max=("perm_p", "max"), realised_spread_mean=("realised_spread", "mean")).reset_index()
    S["method"] = pd.Categorical(S.method, METHODS); S = S.sort_values(["cohort", "spread", "method"])
    S.to_csv(HERE / "depth_spread_summary.csv", index=False)
    # per deep type as well
    S2 = A.groupby(["cohort", "deep_type", "spread", "method"]).agg(
        top_share_mean=("top_share", "mean"), deep_wins=("deep_wins", "sum"), n=("deep_wins", "size"),
        markers_correct_mean=("markers_correct", "mean")).reset_index()
    S2.to_csv(HERE / "depth_spread_summary_by_deep_type.csv", index=False)
    # mechanism table: per (cohort, spread) Spearman across types of share vs detection / lognorm mean
    Tm = T[T.spread != "natural"].copy(); Tm["spread"] = Tm.spread.astype(int)
    mech = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for (c, S_, d, s), grp in Tm.groupby(["cohort", "spread", "deep_type", "seed"]):
            mech.append(dict(cohort=c, spread=S_, deep_type=d, seed=s,
                             lognorm_mean_range=grp.mean_lognorm_profile.max() / grp.mean_lognorm_profile.min(),
                             cpm_mean_range=grp.mean_cpm_profile.max() / grp.mean_cpm_profile.min(),
                             raw_mean_range=grp.mean_raw_profile.max() / grp.mean_raw_profile.min(),
                             detection_range=grp.detection_rate.max() / grp.detection_rate.min(),
                             rho_detection_vs_lognormmean=spearmanr(grp.detection_rate, grp.mean_lognorm_profile)[0],
                             rho_depth_vs_detection=spearmanr(grp.realised_median_rowsum, grp.detection_rate)[0],
                             rho_lognormmean_vs_share_lognorm=spearmanr(grp.mean_lognorm_profile, grp.share_lognorm_mean)[0],
                             rho_rawmean_vs_share_raw=spearmanr(grp.mean_raw_profile, grp.share_raw_mean)[0]))
    Mech = pd.DataFrame(mech)
    Mech.to_csv(HERE / "depth_spread_mechanism.csv", index=False)
    MechS = Mech.groupby(["cohort", "spread"]).mean(numeric_only=True).drop(columns=["seed"]).reset_index()
    MechS.to_csv(HERE / "depth_spread_mechanism_summary.csv", index=False)

    # ── figure ──
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b",
                         "xtick.color": "#898781", "ytick.color": "#898781", "axes.titlecolor": "#0b0b0b",
                         "grid.color": "#e1e0d9", "grid.linewidth": 0.6})
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), dpi=600, sharex=True)
    letters = "ABCD"
    for j, c in enumerate(COHORTS):
        Sc = S[S.cohort == c]; nat = R[(R.cohort == c) & (R.spread == "natural")]
        Tn = T[(T.cohort == c) & (T.spread == "natural")]
        nat_spread = Tn.realised_median_rowsum.max() / Tn.realised_median_rowsum.min()
        k = Tn.shape[0]; ntot = int(Sc.markers_total.iloc[0])
        for i, (col, ylab) in enumerate([("top_share_mean", "top argmax share"),
                                         ("markers_correct_mean", f"held-out markers correct (of {ntot})")]):
            ax = axes[i, j]
            ax.set_xscale("log"); ax.grid(True, axis="y", zorder=0); ax.set_axisbelow(True)
            lo, hi = ("top_share_min", "top_share_max") if i == 0 else ("markers_correct_min", "markers_correct_max")
            for m in METHODS:
                s = Sc[Sc.method == m]
                ax.fill_between(s.spread, s[lo], s[hi], color=COLOR[m], alpha=0.12, lw=0, zorder=1)
                ax.plot(s.spread, s[col], color=COLOR[m], lw=2, marker=MARK[m], ms=4.5, mec="white", mew=0.6,
                        label=LABEL[m], zorder=3)
                nv = nat[nat.method == m]
                yv = nv.top_share.iloc[0] if i == 0 else nv.markers_correct.iloc[0]
                ax.plot([nat_spread], [yv], marker=MARK[m], ms=6, mfc="none", mec=COLOR[m], mew=1.2, ls="none", zorder=4)
            if i == 0:
                ax.axhline(1 / k, color="#898781", lw=0.8, ls=(0, (4, 3)), zorder=2)
                ax.text(1.05, 1 / k + 0.015, f"uniform (1/{k})", color="#52514e", fontsize=6.5, va="bottom")
                ax.set_ylim(0, 1); ax.set_yticks(np.arange(0, 1.01, 0.2))
                ax.set_yticklabels([f"{v:.0%}" for v in np.arange(0, 1.01, 0.2)])
                ax.set_title(f"{ {'LUAD':'NSCLC'}.get(c, c) } count matrix ({int(Tn.n_cells.sum()):,} cells, {k} types)", fontsize=8.5, loc="left")
            else:
                ax.set_ylim(0, ntot * 1.05); ax.set_xlabel("imposed depth spread (deepest / shallowest type)")
                ax.set_xticks(SPREADS); ax.set_xticklabels([f"{s}x" for s in SPREADS])
                ax.tick_params(axis="x", which="minor", bottom=False)
            ax.set_ylabel(ylab)
            ax.text(-0.14, 1.04, letters[i * 2 + j], transform=ax.transAxes, fontsize=11, fontweight="bold")
            ax.axvline(nat_spread, color="#c3c2b7", lw=0.6, ls=":", zorder=1)
            ax.set_xlim(0.85, 75)
            if i == 0:
                ax.text(nat_spread * 1.04, 0.97, f"natural {nat_spread:.1f}x\n(open markers)", color="#52514e",
                        fontsize=6, va="top")
                # direct labels at the right end (pairs of near-identical curves share one label)
                end = {m: Sc[(Sc.method == m) & (Sc.spread == 40)][col].iloc[0] for m in METHODS}
                for txt, ms_, dy in [("raw / z-scale", ("raw_mean", "zscale_mean"), 0.0),
                                     ("log1p(CPM)", ("lognorm_mean",), 0.0),
                                     ("pseudobulk / CPM mean", ("pseudobulk_cpm", "cpm_mean"), 0.0)]:
                    y = np.mean([end[m] for m in ms_]) + dy
                    ax.text(44, y, txt, color="#52514e", fontsize=6, va="center", ha="left")
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.0),
               handlelength=2.4, columnspacing=1.4)
    fig.text(0.005, 0.005, "Lines: mean over 2 deep-type choices x 3 seeds; band: min-max. "
             "Cells thinned to a common p20 depth, then types thinned to the imposed spread.",
             fontsize=6, color="#52514e")
    fig.tight_layout(rect=(0, 0.02, 1, 0.965))
    fig.savefig(HERE / "fig_depth_spread.png", dpi=600)
    fig.savefig(HERE / "fig_depth_spread.pdf")
    log("wrote summary CSVs and fig_depth_spread.png")
    # console summary
    for c in COHORTS:
        print(f"\n=== {c}: top argmax share, mean over 2 deep types x 3 seeds (deep type wins / 6) ===")
        print(f"{'method':16}" + "".join(f"{s:>5}x        " for s in SPREADS))
        for m in METHODS:
            s = S[(S.cohort == c) & (S.method == m)].set_index("spread")
            print(f"{m:16}" + "".join(f"{s.loc[sp_,'top_share_mean']:5.0%} ({int(s.loc[sp_,'deep_wins'])}/6)   " for sp_ in SPREADS))
        print(f"  markers correct (of {int(S[S.cohort==c].markers_total.iloc[0])}), mean:")
        for m in METHODS:
            s = S[(S.cohort == c) & (S.method == m)].set_index("spread")
            print(f"{m:16}" + "".join(f"{s.loc[sp_,'markers_correct_mean']:5.1f}        " for sp_ in SPREADS))


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2])
    elif sys.argv[1] == "summarise":
        summarise()
