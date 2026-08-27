"""Figure: top argmax share (A) and held-out marker concordance (B) for all aggregation
methods on the six real matrices; bootstrap 95 % intervals (B = 200 stratified cell
resamples) for the six original methods. Reads benchmark_all_methods.csv and
benchmark_bootstrap_ci.csv; writes fig_normalisers_ci.png / .pdf."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
A = pd.read_csv(HERE / "benchmark_all_methods.csv")
C = pd.read_csv(HERE / "benchmark_bootstrap_ci.csv")
ORDER = ["raw_mean", "zscale_mean", "rank_mean", "rankfrac_mean", "detection_mean", "tfidf1_mean",
         "lognorm_mean", "tfidf3_mean", "lognorm_maxnorm", "pearson_mean", "depth_matched",
         "pseudobulk_cpm", "pseudobulk_logcpm", "cpm_mean", "medscale_mean", "tfidf2_mean"]
LABEL = {"raw_mean": "raw mean", "zscale_mean": "z-scale mean (scanpy scale)", "rank_mean": "within-cell rank mean",
         "rankfrac_mean": "within-cell rank / nnz mean", "detection_mean": "detection rate (binarised mean)",
         "tfidf1_mean": "TF-IDF, Signac m1: log1p(TF·IDF·1e4) mean", "lognorm_mean": "log1p(CPM) mean (scanpy/Seurat/ArchR)",
         "tfidf3_mean": "TF-IDF, Signac m3: log1p(TF·1e4)·log(1+IDF)", "lognorm_maxnorm": "log1p(CPM) mean, then per-gene max-norm",
         "pearson_mean": "Pearson residual mean (θ = 100)", "depth_matched": "depth-matched thinning, raw mean",
         "pseudobulk_cpm": "pseudobulk sum → CPM", "pseudobulk_logcpm": "pseudobulk sum → log1p(CPM)",
         "cpm_mean": "CPM mean (ours)", "medscale_mean": "median-total scaling mean", "tfidf2_mean": "TF-IDF, Signac m2: TF·log(1+IDF) mean"}
COHORT_COL = {"HCC": "#2a78d6", "COAD": "#eb6834", "LUAD": "#1baf7a"}
INK, MUTED, GRID = "#0b0b0b", "#898781", "#e1e0d9"
MATS = [("HCC", "proximal"), ("HCC", "enhancer"), ("COAD", "proximal"), ("COAD", "enhancer"), ("LUAD", "proximal"), ("LUAD", "enhancer")]
OFF = np.linspace(-0.3, 0.3, len(MATS))

plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(10.5, 6.2), sharey=True, gridspec_kw=dict(wspace=0.08, left=0.30, right=0.98, top=0.90, bottom=0.12))
y = {m: len(ORDER) - 1 - k for k, m in enumerate(ORDER)}
for ax, metric in zip(axes, ["top_share", "marker_frac"]):
    ax.set_facecolor("#fcfcfb"); ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for k in range(len(ORDER)):
        if k % 2: ax.axhspan(k - 0.5, k + 0.5, color="#f3f2ee", zorder=0, lw=0)
    for m in ORDER:
        for (c, mod), off in zip(MATS, OFF):
            r = A[(A.cohort == c) & (A.model == mod) & (A.method == m)]
            if r.empty: continue
            r = r.iloc[0]
            val = r.top_share if metric == "top_share" else r.markers_correct / r.markers_total
            ci = C[(C.cohort == c) & (C.model == mod) & (C.method == m)]
            if not ci.empty:
                ci = ci.iloc[0]
                lo, hi = ((ci.top_share_ci_lo, ci.top_share_ci_hi) if metric == "top_share"
                          else (ci.markers_ci_lo / ci.markers_total, ci.markers_ci_hi / ci.markers_total))
                ax.plot([lo, hi], [y[m] + off] * 2, color=COHORT_COL[c], lw=1.6, solid_capstyle="round", alpha=0.9, zorder=2)
            ax.plot(val, y[m] + off, marker="o", ms=5.5, mfc=(COHORT_COL[c] if mod == "proximal" else "#fcfcfb"),
                    mec=COHORT_COL[c], mew=1.4, zorder=3)
    ax.set_ylim(-0.6, len(ORDER) - 0.4)
    ax.set_xlim(0, 1.0); ax.set_xticks(np.arange(0, 1.01, 0.2)); ax.set_xticklabels([f"{v:.0%}" for v in np.arange(0, 1.01, 0.2)])
axes[0].set_yticks([y[m] for m in ORDER]); axes[0].set_yticklabels([LABEL[m] for m in ORDER])
for lab_ in axes[0].get_yticklabels():
    if lab_.get_text() == LABEL["cpm_mean"]: lab_.set_fontweight("bold")
axes[0].set_xlabel("share of expressed genes whose argmax falls in one cell type")
axes[1].set_xlabel("held-out marker concordance (fraction of panel)")
axes[0].set_title("A  top argmax share", loc="left", fontsize=9.5, color=INK)
axes[1].set_title("B  held-out marker concordance", loc="left", fontsize=9.5, color=INK)
axes[1].tick_params(axis="y", length=0)
handles = [Line2D([], [], marker="o", ls="", ms=5.5, mfc=COHORT_COL[c], mec=COHORT_COL[c], label=f"{c} proximal") for c in COHORT_COL]
handles += [Line2D([], [], marker="o", ls="", ms=5.5, mfc="#fcfcfb", mec=COHORT_COL[c], mew=1.4, label=f"{c} enhancer") for c in COHORT_COL]
handles += [Line2D([], [], color=MUTED, lw=1.6, label="95 % bootstrap interval")]
fig.legend(handles=handles, loc="upper center", ncol=7, frameon=False, fontsize=7, bbox_to_anchor=(0.62, 0.995), handletextpad=0.3, columnspacing=0.9)
fig.savefig(HERE / "fig_normalisers_ci.png", dpi=300); fig.savefig(HERE / "fig_normalisers_ci.pdf")
print("wrote", HERE / "fig_normalisers_ci.png")
