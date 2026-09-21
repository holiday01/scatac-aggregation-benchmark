"""04: figure for the multiome ground-truth benchmark (reads multiome_benchmark.csv)."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'axes.labelsize': 10, 'ytick.labelsize': 9.5})

HERE = Path(__file__).resolve().parent
B = pd.read_csv(HERE / "multiome_benchmark.csv")
AR = pd.read_csv(HERE / "multiome_annotation_repeats.csv")
AGG = AR.groupby(["model", "method"]).agg(["mean", "min", "max"])
ORDER = ["raw_mean", "zscale_mean", "lognorm_mean", "archr_style", "pseudobulk_cpm", "depth_matched", "cpm_mean"]
LABEL = {"raw_mean": "raw mean", "zscale_mean": "z-scaled mean", "lognorm_mean": "log-normalized mean",
         "archr_style": "ArchR log2 mean", "pseudobulk_cpm": "pooled CPM", "depth_matched": "coverage-matched",
         "cpm_mean": "CP10K mean"}
COL = {"proximal": "#4C78A8", "enhancer": "#F58518"}
panels = [("top_share", "Top argmax share\n(all ATAC-expressed genes)", (0, 1)),
          ("argmax_agree_informative", "ATAC argmax = RNA argmax\n(RNA-informative genes)", (0, 1)),
          ("gene_rho_mean_informative", "Mean per-gene Spearman\nATAC vs RNA profile", None),
          ("markers_correct", "Held-out markers correct", None),
          ("acc_spearman", "Annotation accuracy\n(rank correlation, 5 splits)", (0, 1)),
          ("balacc_spearman", "Balanced annotation accuracy\n(rank correlation, 5 splits)", (0, 1))]
fig, axes = plt.subplots(2, 3, figsize=(12.6, 7.6))
x = np.arange(len(ORDER)); w = 0.38
for ax, (col, title, ylim) in zip(axes.ravel(), panels):
    repeat = col in ("acc_spearman", "balacc_spearman")
    for k, model in enumerate(["proximal", "enhancer"]):
        if repeat:
            s = AGG.loc[model].reindex(ORDER)
            vals = s[(col, "mean")].values
            err = np.vstack([vals - s[(col, "min")].values, s[(col, "max")].values - vals])
        else:
            s = B[B.model == model].set_index("method").reindex(ORDER)
            vals = s[col].values; err = None
        ax.bar(x + (k - 0.5) * w, vals, w, color=COL[model], label=f"{model} model",
               yerr=err, error_kw=dict(lw=0.8, capsize=1.6, ecolor="#444444"))
        for xi, v in zip(x + (k - 0.5) * w, vals):
            if np.isfinite(v):
                ax.text(xi, v, f"{v:.2f}" if col != "markers_correct" else f"{int(v)}", ha="center", va="bottom", fontsize=8.5)
    if col == "markers_correct":
        ax.axhline(B.null_mean.iloc[0], ls="--", color="grey", lw=1); ax.text(6.4, B.null_mean.iloc[0], "perm null", fontsize=8, va="bottom", ha="right", color="grey")
        ax.set_ylim(0, B.markers_total.iloc[0] + 3)
    elif ylim: ax.set_ylim(*ylim)
    ax.set_xticks(x); ax.set_xticklabels([LABEL[m] for m in ORDER], rotation=40, ha="right", fontsize=9.5)
    ax.set_title(title, fontsize=11); ax.spines[["top", "right"]].set_visible(False)
axes[0, 0].legend(fontsize=9.5, frameon=False)
d = B.iloc[0]
fig.suptitle(f"10x PBMC multiome (n={int(d.n_cells):,} RNA-labeled cells, {int(d.n_types)} types, "
             f"coverage spread {d.depth_spread:.1f}x): ATAC aggregation against the genome-wide RNA reference", fontsize=12)
fig.tight_layout()
fig.savefig(HERE / "fig_multiome_benchmark.png", dpi=400); fig.savefig(HERE / "fig_multiome_benchmark.pdf")
print("saved", HERE / "fig_multiome_benchmark.png")
