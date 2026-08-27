"""
Figure 1 of the benchmark paper.

  (A) share of all expressed genes whose cell-type argmax falls in a single
      cell type, for six gene-activity matrices x seven aggregation methods;
  (B) held-out marker--lineage concordance for the same grid, plus the
      mis-annotated negative-control matrix;
  (C) the mechanism, on the colorectal proximal matrix: the per-type mean of
      log1p(CPM) rises with the type's detection rate (which tracks depth),
      whereas the per-type mean of CPM is constant by construction.

Drawn at exactly 16 cm / 600 dpi and saved without a tight bounding box, so the
canvas LaTeX places at \\textwidth is the canvas matplotlib drew.

Inputs : results/benchmark_aggregation.csv, results/benchmark_depth_by_type.csv,
         the COAD proximal matrix (for panel C)
Output : figures/fig1_benchmark.png
"""
import os
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
OUT = ROOT / "figures"; OUT.mkdir(exist_ok=True)
TEXTW, DPI = 6.3, 600
plt.rcParams.update({"font.size": 7.5, "font.family": "DejaVu Sans", "axes.linewidth": 0.7,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7})

METHODS = ["raw_mean", "zscale_mean", "lognorm_mean", "pseudobulk_cpm", "depth_matched", "cpm_mean"]
MLAB = ["raw mean", "z-scaled mean", "log-normalised mean", "pseudobulk CPM", "depth-matched", "CPM mean"]
MATS = [("COAD", "proximal"), ("COAD", "enhancer"), ("LUAD", "proximal"), ("LUAD", "enhancer"),
        ("HCC", "proximal"), ("HCC", "enhancer")]
NEG = ("HCC", "misannotated")


def mlabel(c, m):
    return f"{c} {'prox.' if m == 'proximal' else 'enh.' if m == 'enhancer' else 'mis-annot.'}"


def save(fig, name):
    fig.savefig(OUT / name, dpi=DPI, facecolor="white", bbox_inches=None, pad_inches=0)
    plt.close(fig)
    got = Image.open(OUT / name).size
    assert got[0] == round(TEXTW * DPI), got
    print(f"saved {name}  {got[0]}x{got[1]} px")


def panel_c_data():
    import anndata as ad, scipy.sparse as sp, h5py, hdf5plugin  # noqa
    a = ad.read_h5ad(BASE / "data/processed/COAD/coad_geneactivity_proximal.h5ad")
    X = a.X.tocsr().astype(np.float32); lab = a.obs.cell_type.astype(str).values
    tot = np.asarray(X.sum(1)).ravel(); Xc = (sp.diags(1e4 / np.maximum(tot, 1)) @ X).tocsr()
    Xl = Xc.copy(); Xl.data = np.log1p(Xl.data)
    nnz = np.diff(X.indptr) / X.shape[1]
    rows = []
    for t in sorted(set(lab)):
        m = lab == t
        rows.append(dict(type=t, det=float(nnz[m].mean()), cpm=float(Xc[m].mean()), log=float(Xl[m].mean())))
    d = pd.DataFrame(rows)
    d.to_csv(ROOT / "results" / "benchmark_mechanism_coad.csv", index=False)
    return d


def main():
    R = pd.read_csv(ROOT / "results" / "benchmark_aggregation.csv")
    R["key"] = list(zip(R.cohort, R.model))
    share = np.array([[R[(R.key == k) & (R.method == m)].top_share.iloc[0] for m in METHODS] for k in MATS])
    mats_b = MATS + [NEG]
    conc = np.array([[R[(R.key == k) & (R.method == m)].markers_correct.iloc[0] for m in METHODS] for k in mats_b])
    totm = np.array([R[(R.key == k)].markers_total.iloc[0] for k in mats_b])
    frac = conc / totm[:, None]
    C = panel_c_data()

    fig = plt.figure(figsize=(TEXTW, 5.6), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 0.85], width_ratios=[1, 1])
    cm_share = LinearSegmentedColormap.from_list("s", ["#f7f7f7", "#fdae61", "#b2182b"])
    cm_conc = LinearSegmentedColormap.from_list("c", ["#f7f7f7", "#92c5de", "#2166ac"])

    # ── A ──
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(share, cmap=cm_share, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(METHODS))); ax.set_xticklabels(MLAB, fontsize=6.2, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(MATS))); ax.set_yticklabels([mlabel(*k) for k in MATS], fontsize=6.5)
    for i in range(share.shape[0]):
        for j in range(share.shape[1]):
            v = share[i, j]
            ax.text(j, i, f"{v:.0%}", ha="center", va="center", fontsize=6,
                    color="white" if v > 0.55 else "#222222")
    ax.set_title("A  Share of argmaxima in one cell type", loc="left", fontsize=8, fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.ax.tick_params(labelsize=6)
    ax.tick_params(length=2)

    # ── B ──
    ax = fig.add_subplot(gs[0, 1])
    im = ax.imshow(frac, cmap=cm_conc, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(METHODS))); ax.set_xticklabels(MLAB, fontsize=6.2, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(mats_b)))
    ax.set_yticklabels([f"{mlabel(*k)} (n={t})" for k, t in zip(mats_b, totm)], fontsize=6.5)
    for i in range(frac.shape[0]):
        for j in range(frac.shape[1]):
            ax.text(j, i, f"{conc[i, j]}", ha="center", va="center", fontsize=6,
                    color="white" if frac[i, j] > 0.6 else "#222222")
    ax.axhline(len(MATS) - 0.5, color="#333333", lw=0.8, ls="--")
    ax.set_title("B  Held-out markers correct", loc="left", fontsize=8, fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.ax.tick_params(labelsize=6)
    cb.set_label("fraction of panel", fontsize=6)
    ax.tick_params(length=2)

    # ── C ──
    ax = fig.add_subplot(gs[1, :])
    C = C.sort_values("det")
    ax.plot(C.det * 100, C.log / C.log.max(), "o-", color="#b2182b", ms=4, lw=1, label="mean of log1p(CPM), all genes")
    ax.plot(C.det * 100, C.cpm / C.cpm.max(), "s-", color="#2166ac", ms=4, lw=1, label="mean of CPM, all genes")
    shallow = list(C.type[:3])
    for i, r in enumerate(C.itertuples()):
        if r.type in shallow:
            continue
        ax.annotate(r.type.replace("_", " "), (r.det * 100, r.log / C.log.max()), xytext=(4, -9 if r.type == "B_cell" else 4),
                    textcoords="offset points", fontsize=6, color="#b2182b")
    ax.annotate(", ".join(t.replace("_", " ") for t in shallow), (C.det.iloc[2] * 100, C.log.iloc[2] / C.log.max()),
                xytext=(10, -14), textcoords="offset points", fontsize=6, color="#b2182b",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#b2182b"))
    ax.set_xlabel("genes detected per cell (%), colorectal proximal matrix, per cell type", fontsize=7)
    ax.set_ylabel("per-type mean / max over types", fontsize=7)
    ax.set_ylim(0, 1.08); ax.legend(fontsize=6.5, frameon=False, loc="lower right")
    ax.set_title("C  Why the log-normalised mean still follows depth", loc="left", fontsize=8, fontweight="bold")
    ax.tick_params(labelsize=6.5, length=2)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

    save(fig, "fig1_benchmark.png")


if __name__ == "__main__":
    main()
