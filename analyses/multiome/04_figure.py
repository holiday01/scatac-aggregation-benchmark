"""04: figure for the multiome ground-truth benchmark (reads multiome_benchmark.csv)."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
B = pd.read_csv(HERE / "multiome_benchmark.csv")
ORDER = ["raw_mean", "zscale_mean", "lognorm_mean", "archr_style", "pseudobulk_cpm", "depth_matched", "cpm_mean"]
LABEL = {"raw_mean": "raw mean", "zscale_mean": "z-scale mean", "lognorm_mean": "log1p(CPM) mean",
         "archr_style": "ArchR log2 mean", "pseudobulk_cpm": "pseudobulk CPM", "depth_matched": "depth-matched",
         "cpm_mean": "CPM mean (ours)"}
COL = {"proximal": "#4C78A8", "enhancer": "#F58518"}
panels = [("top_share", "Top argmax share\n(all ATAC-expressed genes)", (0, 1)),
          ("argmax_agree_informative", "ATAC argmax = RNA argmax\n(RNA-informative genes)", (0, 1)),
          ("gene_rho_mean_informative", "Mean per-gene Spearman\nATAC vs RNA profile", None),
          ("markers_correct", "Held-out markers correct", None),
          ("annot_acc_pearson", "Annotation accuracy\n(held-out 50 %, Pearson)", (0, 1)),
          ("annot_balacc_pearson", "Balanced annotation accuracy\n(held-out 50 %, Pearson)", (0, 1))]
fig, axes = plt.subplots(2, 3, figsize=(13, 7.2))
x = np.arange(len(ORDER)); w = 0.38
for ax, (col, title, ylim) in zip(axes.ravel(), panels):
    for k, model in enumerate(["proximal", "enhancer"]):
        s = B[B.model == model].set_index("method").reindex(ORDER)
        vals = s[col].values
        ax.bar(x + (k - 0.5) * w, vals, w, color=COL[model], label=f"{model} model")
        for xi, v in zip(x + (k - 0.5) * w, vals):
            if np.isfinite(v):
                ax.text(xi, v, f"{v:.2f}" if col != "markers_correct" else f"{int(v)}", ha="center", va="bottom", fontsize=7)
    if col == "markers_correct":
        ax.axhline(B.null_mean.iloc[0], ls="--", color="grey", lw=1); ax.text(6.4, B.null_mean.iloc[0], "perm null", fontsize=7, va="bottom", ha="right", color="grey")
        ax.set_ylim(0, B.markers_total.iloc[0] + 3)
    elif ylim: ax.set_ylim(*ylim)
    ax.set_xticks(x); ax.set_xticklabels([LABEL[m] for m in ORDER], rotation=40, ha="right", fontsize=8)
    ax.set_title(title, fontsize=10); ax.spines[["top", "right"]].set_visible(False)
axes[0, 0].legend(fontsize=8, frameon=False)
d = B.iloc[0]
fig.suptitle(f"10x PBMC multiome (n={int(d.n_cells):,} RNA-labelled cells, {int(d.n_types)} types, "
             f"true-depth spread {d.depth_spread:.1f}x): ATAC aggregation vs genome-wide RNA truth", fontsize=11)
fig.tight_layout()
fig.savefig(HERE / "fig_multiome_benchmark.png", dpi=200); fig.savefig(HERE / "fig_multiome_benchmark.pdf")
print("saved", HERE / "fig_multiome_benchmark.png")
