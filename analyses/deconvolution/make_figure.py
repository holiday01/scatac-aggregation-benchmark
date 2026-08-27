"""Figure: NNLS deconvolution error per reference-aggregation method (row 1: RMSE; row 2: bias of the
deepest cell type), one panel per cohort, rawsum mixtures.  Two series = the two definitions of the
true proportion (cell fraction; fragment share).  One file per gene-set variant (common / own)."""
import sys
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
S = pd.read_csv(HERE / "metrics_summary.csv")
METHODS = ["raw_mean", "zscale_mean", "lognorm_mean", "lognorm_mean_expm1", "pseudobulk_cpm", "depth_matched", "cpm_mean"]
LABELS = ["raw mean", "z-scale mean", "log1p(CPM) mean", "expm1(log mean)", "pseudobulk CPM", "depth-matched", "CPM mean"]
COHORTS = ["HCC", "COAD", "LUAD"]
SERIES = [("cellfrac", "true = cell fraction", "#2a78d6"), ("fragfrac", "true = fragment share", "#eb6834")]
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlecolor": INK, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF})


def draw(geneset, mixture="rawsum"):
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.6), sharex="col")
    x = np.arange(len(METHODS)); w = 0.38
    for j, co in enumerate(COHORTS):
        sub = S[(S.cohort == co) & (S.geneset == geneset) & (S.mixture_variant == mixture)]
        deepest = sub.deepest_type.iloc[0]
        for i, (metric, ylabel) in enumerate([("rmse", "RMSE (estimated - true proportion)"),
                                              ("bias_deepest", f"bias of deepest type")]):
            ax = axes[i, j]
            for k, (truth, name, col) in enumerate(SERIES):
                d = sub[sub.truth == truth].set_index("method").loc[METHODS]
                ax.bar(x + (k - 0.5) * w, d[f"{metric}_mean"], w * 0.92, color=col, label=name, zorder=3)
                ax.errorbar(x + (k - 0.5) * w, d[f"{metric}_mean"], yerr=d[f"{metric}_sd"], fmt="none",
                            ecolor=INK2, elinewidth=0.8, capsize=2, zorder=4)
            ax.axhline(0, color=AXIS, lw=0.8, zorder=2)
            ax.grid(axis="y", color=GRID, lw=0.6, zorder=0); ax.set_axisbelow(True)
            for s in ["top", "right"]: ax.spines[s].set_visible(False)
            ax.tick_params(length=0)
            if i == 0:
                ax.set_title(f"{co}  (deepest type: {deepest.replace('_', ' ')})", fontsize=9, loc="left")
                ax.set_ylim(0, max(0.2, sub.rmse_mean.max() * 1.15))
            else:
                lim = max(0.12, np.abs(sub.bias_deepest_mean).max() * 1.3); ax.set_ylim(-lim, lim)
                ax.set_xticks(x); ax.set_xticklabels(LABELS, fontsize=7.5, rotation=35, ha="right", rotation_mode="anchor")
            if j == 0:
                ax.set_ylabel(ylabel if i == 0 else f"bias of deepest type\n(mean estimated - true)")
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper left", ncol=2, frameon=False, fontsize=8, bbox_to_anchor=(0.01, 0.955))
    fig.suptitle(f"NNLS deconvolution of 300 pseudo-bulk mixtures (3 seeds x 100 cells-drawn mixtures) against references "
                 f"built by each aggregation method\n{'shared' if geneset == 'common' else 'per-method'} signature genes; "
                 f"mixtures = {'sum of raw counts of the sampled cells' if mixture == 'rawsum' else 'sum of per-cell-normalised (10k) cells'}; "
                 f"error bars = sd over seeds",
                 fontsize=8.5, x=0.01, y=0.995, ha="left", va="top", color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out = HERE / f"fig_deconvolution_{geneset}_{mixture}.png"
    fig.savefig(out, dpi=200); plt.close(fig); print("wrote", out)


if __name__ == "__main__":
    for gs in ["common", "own"]:
        draw(gs, "rawsum")
    draw("common", "cellnorm")
