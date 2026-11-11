"""Earlier three-figure layout of the display figures (tools/make_figures_v2.py draws the current four).

fig1_benchmark.png   A share of argmaxima in one cell type, B held-out marker
                     concordance, for the five aggregation methods that are
                     defined on both count and weighted-score matrices, over the
                     six matrices plus the mis-annotated control; C the
                     detection-rate mechanism on the colorectal count matrix.
                     The depth-matched control is omitted here because binomial
                     thinning is defined only for the three integer count
                     matrices; it is reported in the text.
fig2_depth_spread.png  copied from the controlled-thinning analysis.
fig3_multiome.png      copied from the matched-modality analysis.
"""
import os
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "/mnt/10t/scrna_atac"))
OUT = ROOT / "figures"; OUT.mkdir(exist_ok=True)
TEXTW, DPI = 6.3, 600
plt.rcParams.update({"font.size": 7.5, "font.family": "DejaVu Sans", "axes.linewidth": 0.7,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7})

METHODS = ["raw_mean", "zscale_mean", "lognorm_mean", "pseudobulk_cpm", "cpm_mean"]
MLAB = ["raw mean", "z-scaled mean", "log-normalized mean", "pooled CPM", "CP10K mean"]
MATS = [("COAD", "proximal"), ("COAD", "enhancer"), ("LUAD", "proximal"), ("LUAD", "enhancer"),
        ("HCC", "proximal"), ("HCC", "enhancer")]
NEG = ("HCC", "misannotated")
DISPLAY = {'LUAD': 'NSCLC'}
SPREAD = {"COAD": "36.8x", "LUAD": "3.2x", "HCC": "2.0x"}
lbl = lambda c, m: (f"{DISPLAY.get(c, c)} {'count' if m=='proximal' else 'weighted'} ({SPREAD[c]})"
                    if m in ("proximal", "enhancer") else f"{DISPLAY.get(c, c)} mislabeled")


def save(fig, name):
    fig.savefig(OUT / name, dpi=DPI, facecolor="white", bbox_inches=None, pad_inches=0)
    plt.close(fig)
    got = Image.open(OUT / name).size
    assert got[0] == round(TEXTW * DPI), got
    print(f"saved {name}  {got[0]}x{got[1]} px")


def main():
    R = pd.read_csv(ROOT / "results" / "benchmark_aggregation.csv")
    R["key"] = list(zip(R.cohort, R.model))
    C = pd.read_csv(ROOT / "results" / "benchmark_mechanism_coad.csv").sort_values("det")
    LO = pd.read_csv(ROOT / "results" / "revision" / "loso_top_share.csv")
    share = np.array([[R[(R.key == k) & (R.method == m)].top_share.iloc[0] for m in METHODS] for k in MATS])
    mats_b = MATS + [NEG]
    conc = np.array([[R[(R.key == k) & (R.method == m)].markers_correct.iloc[0] for m in METHODS] for k in mats_b])
    totm = np.array([R[R.key == k].markers_total.iloc[0] for k in mats_b])
    frac = conc / totm[:, None]

    fig = plt.figure(figsize=(TEXTW, 5.9), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 0.85])
    cm_s = LinearSegmentedColormap.from_list("s", ["#f7f7f7", "#fdae61", "#b2182b"])
    cm_c = LinearSegmentedColormap.from_list("c", ["#f7f7f7", "#92c5de", "#2166ac"])

    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(share, cmap=cm_s, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(METHODS))); ax.set_xticklabels(MLAB, fontsize=6.2, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(MATS))); ax.set_yticklabels([lbl(*k) for k in MATS], fontsize=6.5)
    for i in range(share.shape[0]):
        for j in range(share.shape[1]):
            ax.text(j, i, f"{share[i,j]:.0%}", ha="center", va="center", fontsize=6,
                    color="white" if share[i, j] > 0.55 else "#222222")
    ax.set_title("A  Largest cell-type share of per-gene maxima", loc="left", fontsize=8, fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.ax.tick_params(labelsize=6)
    ax.tick_params(length=2)

    ax = fig.add_subplot(gs[0, 1])
    im = ax.imshow(frac, cmap=cm_c, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(METHODS))); ax.set_xticklabels(MLAB, fontsize=6.2, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(mats_b)))
    ax.set_yticklabels([f"{lbl(*k)} (n={t})" for k, t in zip(mats_b, totm)], fontsize=6.5)
    for i in range(frac.shape[0]):
        for j in range(frac.shape[1]):
            ax.text(j, i, f"{conc[i,j]}", ha="center", va="center", fontsize=6,
                    color="white" if frac[i, j] > 0.6 else "#222222")
    ax.axhline(len(MATS) - 0.5, color="#333333", lw=0.8, ls="--")
    ax.set_title("B  Held-out markers correct", loc="left", fontsize=8, fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.ax.tick_params(labelsize=6)
    cb.set_label("fraction of panel", fontsize=6)
    ax.tick_params(length=2)

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(C.det * 100, C.log / C.log.max(), "o-", color="#b2182b", ms=4, lw=1, label=None)
    ax.plot(C.det * 100, C.cpm / C.cpm.max(), "s-", color="#2166ac", ms=4, lw=1, label=None)
    for r, dx, dy, ha in [(C.iloc[0], 3, 6, "left"), (C.iloc[-1], -3, -2, "right")]:
        ax.annotate(r.type.replace("_", " "), (r.det * 100, r.log / C.log.max()),
                    xytext=(dx, dy), textcoords="offset points", fontsize=6.2, color="#b2182b", ha=ha)
    ax.set_xlabel("genes detected per cell (%), colorectal count matrix", fontsize=7)
    ax.set_ylabel("per-type mean / max over types", fontsize=7)
    ax.text(0.96, 0.30, r"$\rho$ = +0.95", transform=ax.transAxes, ha="right", va="top", fontsize=7, color="#b2182b")
    ax.set_ylim(0, 1.08); ax.text(C.det.iloc[3]*100, C.log.iloc[3]/C.log.max()-0.16, "mean of log(1+CP10K)", fontsize=6.4, color="#b2182b")
    ax.text(C.det.iloc[-1]*100, 1.03, "mean of CP10K", fontsize=6.4, color="#2166ac", ha="right")
    ax.set_title("C  Log-mean tracks the detection rate", loc="left", fontsize=8, fontweight="bold")
    ax.tick_params(labelsize=6.5, length=2)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    # ── D: leave-one-sample-out ranges ──
    ax = fig.add_subplot(gs[1, 1])
    ORD = [("COAD", "COAD"), ("LUAD", "NSCLC"), ("HCC", "HCC")]
    MSET = [("raw_mean", "raw mean", "#b2182b"), ("lognorm_mean", "log-normalized", "#2f8f5b"),
            ("cpm_mean", "CP10K mean", "#2166ac")]
    ypos, ylab_ = [], []
    for i, (coh, disp) in enumerate(ORD):
        for j, (m, lab_, col) in enumerate(MSET):
            s_ = LO[(LO.cohort == coh) & (LO.method == m)].top_share
            full = R[(R.cohort == coh) & (R.model == "proximal") & (R.method == m)].top_share.iloc[0]
            y = i * 1.0 + (1 - j) * 0.22
            ax.plot([s_.min() * 100, s_.max() * 100], [y, y], color=col, lw=2.6, solid_capstyle="round",
                    zorder=3, label=lab_ if i == 0 else None)
            ax.plot([full * 100], [y], "o", color=col, ms=3.4, mec="white", mew=0.7, zorder=4)
        ypos.append(i); ylab_.append(f"{disp}\n(n={LO[LO.cohort==coh].sample_left_out.nunique()})")
    ax.set_yticks(ypos); ax.set_yticklabels(ylab_, fontsize=6.5)
    ax.set_ylim(-0.5, len(ORD) - 0.5); ax.set_xlim(0, 100)
    ax.set_xlabel("top share when one sample is left out (%)", fontsize=7)
    ax.set_title("D  Leave-one-sample-out", loc="left", fontsize=8, fontweight="bold")
    ax.legend(fontsize=6.2, frameon=False, loc="lower right", handlelength=1.6)
    ax.tick_params(labelsize=6.5, length=2); ax.grid(True, axis="x", lw=0.4, color="#e6e6e6")
    ax.set_axisbelow(True)
    for sp_ in ("top", "right", "left"): ax.spines[sp_].set_visible(False)

    save(fig, "fig1_benchmark.png")

    # Figures 2 and 3 are produced by their own analyses and copied in by the
    # driver; they are not written here.


if __name__ == "__main__":
    main()
