"""Display figures 1-4, redrawn 2026-11-11 for legibility (reviewer comment: low
resolution, small labels). Reads only existing result tables; recomputes nothing.

Shared style: 7.5 in wide (PLOS maximum), 600 dpi, no text below 7 pt at final
size, one color and marker per aggregation strategy across all figures.

  fig1_benchmark.png       results/benchmark_aggregation.csv, results/benchmark_mechanism_coad.csv,
                           results/revision/loso_top_share.csv
  fig2_depth_spread.png    analyses/depth_spread/depth_spread_{summary,results,by_type}.csv
  fig3_multiome.png        analyses/multiome/multiome_{benchmark,annotation_repeats,rna_truth_share}.csv
  fig4_normal_tissues.png  analyses/normal_tissues/*.csv and the above

Supersedes tools/make_revision_figures.py and the per-analysis figure scripts for the display figures.
"""
import os
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent            # .../apbc2026
ADD = ROOT / "analyses"
# FIG_PROFILE=plos (default): 7.5 in wide -> figures/ (read by plos_compbiol/make_ploscb.py)
# FIG_PROFILE=bmc: 170 mm full-page width of BMC / Springer journals -> figures_ec/, PNG + vector PDF
PROFILE = os.environ.get("FIG_PROFILE", "plos")
DPI = 600
if PROFILE == "bmc":
    W, OUT, HEIGHT = 170 / 25.4, ROOT / "figures_ec", {1: 7.3, 2: 5.9, 3: 5.6, 4: 7.7}
else:
    W, OUT, HEIGHT = 7.5, ROOT / "figures", {1: 7.6, 2: 6.2, 3: 5.6, 4: 8.7}
OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5, "axes.titleweight": "bold",
    "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlelocation": "left", "axes.titlepad": 6, "legend.frameon": False,
})

LAB = {"raw_mean": "raw mean", "zscale_mean": "z-scaled mean", "lognorm_mean": "log-normalized mean",
       "pseudobulk_cpm": "pooled CPM", "cpm_mean": "CP10K mean", "depth_matched": "coverage-matched"}
COL = {"raw_mean": "#b2182b", "zscale_mean": "#e08214", "lognorm_mean": "#2f8f5b",
       "pseudobulk_cpm": "#7b5ea7", "cpm_mean": "#2166ac", "depth_matched": "#7f7f7f"}
MRK = {"raw_mean": "o", "zscale_mean": "s", "lognorm_mean": "^", "pseudobulk_cpm": "D", "cpm_mean": "v",
       "depth_matched": "P"}
DISPLAY = {"LUAD": "NSCLC"}
GRID = dict(lw=0.4, color="#e3e3e3")


def save(fig, name, h):
    fig.savefig(OUT / name, dpi=DPI, facecolor="white")
    if PROFILE == "bmc":
        fig.savefig(OUT / name.replace(".png", ".pdf"), facecolor="white")
    plt.close(fig)
    from PIL import Image
    sz = Image.open(OUT / name).size
    assert abs(sz[0] - W * DPI) <= 1 and abs(sz[1] - h * DPI) <= 1, sz
    print(f"saved {name}: {sz[0]}x{sz[1]} px = {sz[0]/DPI:.2f}x{sz[1]/DPI:.2f} in")


# ═════════════════════════════════════════════════════════════════════ Fig 1
def fig1():
    R = pd.read_csv(ROOT / "results" / "benchmark_aggregation.csv")
    R["key"] = list(zip(R.cohort, R.model))
    C = pd.read_csv(ROOT / "results" / "benchmark_mechanism_coad.csv").sort_values("det")
    LO = pd.read_csv(ROOT / "results" / "revision" / "loso_top_share.csv")
    METHODS = ["raw_mean", "zscale_mean", "lognorm_mean", "pseudobulk_cpm", "cpm_mean"]
    MATS = [("COAD", "proximal"), ("COAD", "enhancer"), ("LUAD", "proximal"), ("LUAD", "enhancer"),
            ("HCC", "proximal"), ("HCC", "enhancer")]
    NEG = ("HCC", "misannotated")
    TD = pd.read_csv(ROOT / "results" / "benchmark_depth_by_type.csv")
    SPREAD = {(c, m): f"{g.median_true_depth.max() / g.median_true_depth.min():.1f}×" for (c, m), g in TD.groupby(["cohort", "model"])}
    lbl = lambda c, m: (f"{DISPLAY.get(c, c)} {'count' if m == 'proximal' else 'weighted'} ({SPREAD[(c, m)]})"
                        if m in ("proximal", "enhancer") else f"{DISPLAY.get(c, c)} mislabeled")
    share = np.array([[R[(R.key == k) & (R.method == m)].top_share.iloc[0] for m in METHODS] for k in MATS])
    mats_b = MATS + [NEG]
    conc = np.array([[R[(R.key == k) & (R.method == m)].markers_correct.iloc[0] for m in METHODS] for k in mats_b])
    totm = np.array([R[R.key == k].markers_total.iloc[0] for k in mats_b])
    frac = conc / totm[:, None]

    H = HEIGHT[1]
    fig = plt.figure(figsize=(W, H), layout="constrained")
    fig.get_layout_engine().set(w_pad=0.04, h_pad=0.06, wspace=0.04, hspace=0.04)
    top, bot = fig.subfigures(2, 1, height_ratios=[1.25, 1.0])
    axA, axB = top.subplots(1, 2, width_ratios=[1.0, 0.9])
    axC, axD = bot.subplots(1, 2)
    cm_s = LinearSegmentedColormap.from_list("s", ["#f7f7f7", "#fdae61", "#b2182b"])
    cm_c = LinearSegmentedColormap.from_list("c", ["#f7f7f7", "#92c5de", "#2166ac"])
    shortlbl = lambda c, m: (f"{DISPLAY.get(c, c)} {'count' if m == 'proximal' else 'weighted'}"
                             if m in ("proximal", "enhancer") else f"{DISPLAY.get(c, c)} mislabeled")

    for col, (ax, mat, rows, title, txt, thr, cmap) in enumerate([
            (axA, share, [lbl(*k) for k in MATS], "A  Top share (%)",
             lambda i, j: f"{share[i, j] * 100:.0f}", 0.55, cm_s),
            (axB, frac, [f"{shortlbl(*k)} (of {t})" for k, t in zip(mats_b, totm)], "B  Held-out markers correct",
             lambda i, j: f"{conc[i, j]}", 0.6, cm_c)]):
        im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(METHODS)))
        ax.set_xticklabels([LAB[m] for m in METHODS], rotation=40, ha="right", rotation_mode="anchor")
        ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                ax.text(j, i, txt(i, j), ha="center", va="center", fontsize=7.5,
                        color="white" if mat[i, j] > thr else "#1a1a1a")
        if col == 1:
            ax.axhline(len(MATS) - 0.5, color="#333333", lw=0.9, ls="--")
        for s in ("top", "right"): ax.spines[s].set_visible(True)
        ax.set_title(title)
        ax.tick_params(length=0)

    # C: detection-rate mechanism
    ax = axC
    x = C.det * 100
    ax.plot(x, C.cpm / C.cpm.max(), "-", color=COL["cpm_mean"], lw=1.4, marker="v", ms=5, mec="white", mew=0.5)
    ax.plot(x, C.log / C.log.max(), "-", color=COL["lognorm_mean"], lw=1.4, marker="^", ms=5, mec="white", mew=0.5)
    lo, hi = C.iloc[0], C.iloc[-1]
    low3 = ", ".join(t.replace("_", "/") if t == "T_NK" else t.lower() for t in C.type.iloc[:3])
    ax.annotate(low3, (lo.det * 100, lo.log / C.log.max()), xytext=(8, -4),
                textcoords="offset points", fontsize=7.5, va="top", color="#333333")
    ax.annotate(hi.type.replace("_", "/").lower(), (hi.det * 100, hi.log / C.log.max()), xytext=(-6, 5),
                textcoords="offset points", fontsize=7.5, ha="right", va="bottom", color="#333333")
    ax.text(x.min(), 1.07, "mean of CP10K (equal in every type)", color=COL["cpm_mean"], fontsize=7.5, va="bottom")
    ax.text(8.5, 0.42, "mean of log(1+CP10K)", color=COL["lognorm_mean"], fontsize=7.5)
    ax.text(0.5, 0.93, "log mean vs detection rate:\nSpearman ρ = +0.95", ha="left", va="top",
            fontsize=7.5, color="#333333")
    ax.set_ylim(0, 1.2); ax.set_xlim(0, x.max() * 1.06)
    ax.set_xlabel("genes detected per cell (%), colorectal count matrix")
    ax.set_ylabel("per-type mean / maximum over types")
    ax.set_title("C  The log mean follows the detection rate")
    ax.grid(True, **GRID); ax.set_axisbelow(True)

    # D: leave-one-sample-out
    ax = axD
    ORD = [("COAD", "COAD"), ("LUAD", "NSCLC"), ("HCC", "HCC")]
    MSET = ["raw_mean", "lognorm_mean", "cpm_mean"]
    yt, yl = [], []
    for i, (coh, disp) in enumerate(ORD):
        for j, m in enumerate(MSET):
            s_ = LO[(LO.cohort == coh) & (LO.method == m)].top_share * 100
            full = R[(R.cohort == coh) & (R.model == "proximal") & (R.method == m)].top_share.iloc[0] * 100
            y = (len(ORD) - 1 - i) + (1 - j) * 0.24
            ax.plot([s_.min(), s_.max()], [y, y], color=COL[m], lw=3.2, solid_capstyle="round", zorder=3,
                    label=LAB[m] if i == 0 else None)
            ax.plot([full], [y], "o", color="white", ms=3.2, mec=COL[m], mew=1.0, zorder=4)
        yt.append(len(ORD) - 1 - i); yl.append(f"{disp}\n{LO[LO.cohort == coh].sample_left_out.nunique()} samples")
    ax.set_yticks(yt); ax.set_yticklabels(yl)
    ax.set_ylim(-1.55, 2.5); ax.set_xlim(0, 100)
    ax.set_xlabel("top share with one sample left out (%)")
    ax.set_title("D  Leave-one-sample-out range")
    h, l = ax.get_legend_handles_labels()
    h.append(Line2D([], [], marker="o", color="white", mec="#333333", mew=1.0, ls="", ms=4)); l.append("all samples")
    ax.legend(h, l, loc="lower left", ncol=2, handlelength=1.8, columnspacing=1.2, borderaxespad=0.2)
    ax.grid(True, axis="x", **GRID); ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    save(fig, "fig1_benchmark.png", H)


# ═════════════════════════════════════════════════════════════════════ Fig 2
def fig2():
    D = ADD / "depth_spread"
    S = pd.read_csv(D / "depth_spread_summary.csv")
    Rn = pd.read_csv(D / "depth_spread_results.csv"); Rn = Rn[Rn.spread == "natural"]
    T = pd.read_csv(D / "depth_spread_by_type.csv"); T = T[T.spread == "natural"]
    METHODS = ["raw_mean", "zscale_mean", "lognorm_mean", "pseudobulk_cpm", "cpm_mean"]
    SPREADS = [1, 2, 5, 10, 20, 40]
    H = HEIGHT[2]
    fig, axes = plt.subplots(2, 2, figsize=(W, H), sharex=True, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.04, h_pad=0.05, wspace=0.06, hspace=0.04)
    for j, c in enumerate(["HCC", "LUAD"]):
        Sc = S[S.cohort == c]; nat = Rn[Rn.cohort == c]; Tn = T[T.cohort == c]
        nat_spread = Tn.realised_median_rowsum.max() / Tn.realised_median_rowsum.min()
        k = Tn.shape[0]; ntot = int(Sc.markers_total.iloc[0])
        for i, (col, ylab) in enumerate([("top_share_mean", "top share (%)"),
                                         ("markers_correct_mean", f"held-out markers in own lineage (of {ntot})")]):
            ax = axes[i, j]
            ax.set_xscale("log")
            ax.axvspan(14, 60, color="#8a8a8a", alpha=0.10, lw=0, zorder=0)
            ax.axvline(nat_spread, color="#9a9a9a", lw=0.7, ls=":", zorder=1)
            lo, hi = ("top_share_min", "top_share_max") if i == 0 else ("markers_correct_min", "markers_correct_max")
            sc = 100 if i == 0 else 1
            for m in METHODS:
                s = Sc[Sc.method == m]
                ax.fill_between(s.spread, s[lo] * sc, s[hi] * sc, color=COL[m], alpha=0.13, lw=0, zorder=1)
                ax.plot(s.spread, s[col] * sc, color=COL[m], lw=1.6, marker=MRK[m], ms=4.5, mec="white", mew=0.5,
                        label=LAB[m], zorder=3)
                nv = nat[nat.method == m]
                yv = (nv.top_share.iloc[0] * 100) if i == 0 else nv.markers_correct.iloc[0]
                ax.plot([nat_spread], [yv], marker=MRK[m], ms=6.5, mfc="none", mec=COL[m], mew=1.2, ls="none", zorder=4)
            if i == 0:
                ax.axhline(100 / k, color="#666666", lw=0.8, ls=(0, (4, 3)), zorder=2)
                ax.text(5.5, 100 / k - 1.5, f"uniform share, 1/{k}", color="#444444", fontsize=7.5, va="top")
                ax.set_ylim(0, 100)
                ax.set_title(f"{'AB'[j]}  {DISPLAY.get(c, c)} count matrix, top share")
            else:
                ax.set_ylim(0, ntot * 1.05)
                ax.set_xlabel("imposed row-sum ratio (highest / lowest type)")
                ax.set_title(f"{'CD'[j]}  {DISPLAY.get(c, c)} held-out markers")
            ax.set_xticks(SPREADS); ax.set_xticklabels([f"{s}×" for s in SPREADS])
            ax.tick_params(axis="x", which="minor", bottom=False)
            ax.set_xlim(0.85, 60); ax.set_ylabel(ylab)
            ax.grid(True, axis="y", **GRID); ax.set_axisbelow(True)
    h, l = axes[0, 0].get_legend_handles_labels()
    h += [Line2D([], [], marker="o", ls="none", mfc="none", mec="#333333", mew=1.2, ms=6.5),
          Patch(color="#8a8a8a", alpha=0.25, lw=0)]
    l += ["unthinned matrix", "below usual quality range"]
    fig.legend(h, l, loc="outside lower center", ncol=4, handlelength=2.0, columnspacing=1.4)
    save(fig, "fig2_depth_spread.png", H)


# ═════════════════════════════════════════════════════════════════════ Fig 3
def fig3():
    D = ADD / "multiome"
    B = pd.read_csv(D / "multiome_benchmark.csv")
    AR = pd.read_csv(D / "multiome_annotation_repeats.csv")
    RT = pd.read_csv(D / "multiome_rna_truth_share.csv")
    AGG = AR.groupby(["model", "method"]).agg(["mean", "min", "max"])
    ORDER = ["raw_mean", "zscale_mean", "lognorm_mean", "pseudobulk_cpm", "cpm_mean", "depth_matched"]
    rna_top = RT[RT.model == "proximal"].rna_argmax_share_all.max() * 100
    print(f"fig3: RNA reference top share {rna_top:.1f} %")
    panels = [("top_share", "A  Top share", "top share (%)", 100, (0, 50)),
              ("argmax_agree_informative", "B  Same maximum as RNA", "RNA-informative genes (%)", 100, (0, 100)),
              ("gene_rho_mean_informative", "C  Correlation with RNA", "mean per-gene Spearman ρ", 1, (0, 0.8)),
              ("markers_correct", "D  Held-out markers", f"in own lineage (of {int(B.markers_total.iloc[0])})", 1,
               (0, int(B.markers_total.iloc[0]) + 2)),
              ("acc_spearman", "E  Annotation accuracy", "accuracy", 1, (0, 1)),
              ("balacc_spearman", "F  Balanced accuracy", "balanced accuracy", 1, (0, 1))]
    H = HEIGHT[3]
    fig, axes = plt.subplots(2, 3, figsize=(W, H), layout="constrained")
    fig.get_layout_engine().set(w_pad=0.04, h_pad=0.05, wspace=0.06, hspace=0.06)
    x = np.arange(len(ORDER)); w = 0.4
    for ax, (col, title, ylab, sc, ylim) in zip(axes.ravel(), panels):
        rep = col in ("acc_spearman", "balacc_spearman")
        for kk, model in enumerate(["proximal", "enhancer"]):
            if rep:
                s = AGG.loc[model].reindex(ORDER)
                v = s[(col, "mean")].values
                err = np.vstack([v - s[(col, "min")].values, s[(col, "max")].values - v])
            else:
                v = B[B.model == model].set_index("method").reindex(ORDER)[col].values * sc; err = None
            if model == "enhancer":     # binomial thinning is defined for integer counts: no weighted-model control
                v = np.where(np.array(ORDER) == "depth_matched", np.nan, v)
                if err is not None:
                    err = np.where(np.array(ORDER) == "depth_matched", np.nan, err)
            fc = [COL[m] if kk == 0 else to_rgba(COL[m], 0.38) for m in ORDER]
            ax.bar(x + (kk - 0.5) * w, v, w * 0.94, color=fc, edgecolor=[COL[m] for m in ORDER], linewidth=0.7,
                   yerr=err, error_kw=dict(lw=0.8, capsize=1.5, ecolor="#333333"), zorder=2)
        if col == "top_share":
            ax.axhline(rna_top, color="#333333", lw=0.9, ls="--", zorder=3)
            ax.text(5.6, rna_top + 1.2, f"RNA reference {rna_top:.1f}%", ha="right", va="bottom", fontsize=7.5, zorder=4, bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))
        if col == "markers_correct":
            nm = B.null_mean.iloc[0]
            ax.axhline(nm, color="#333333", lw=0.9, ls="--", zorder=3)
            ax.text(5.6, nm + 1, "permutation null", ha="right", va="bottom", fontsize=7.5, zorder=4, bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))
        ax.set_ylim(*ylim); ax.set_ylabel(ylab); ax.set_title(title)
        ax.set_xticks(x); ax.set_xticklabels([LAB[m] for m in ORDER], rotation=45, ha="right", rotation_mode="anchor")
        ax.grid(True, axis="y", **GRID); ax.set_axisbelow(True)
    h = [Patch(facecolor="#555555", edgecolor="#555555"), Patch(facecolor=to_rgba("#555555", 0.38), edgecolor="#555555")]
    fig.legend(h, ["count model (solid)", "weighted model (light)"], loc="outside lower center", ncol=2)
    save(fig, "fig3_multiome.png", H)


# ═════════════════════════════════════════════════════════════════════ Fig 4
def fig4():
    D = ADD / "normal_tissues"
    NICE = {"colon_transverse": "colon, transverse", "colon_sigmoid": "colon, sigmoid", "heart_lv": "heart, left ventricle",
            "heart_atrial_appendage": "heart, atrial appendage", "artery_aorta": "aorta", "artery_tibial": "tibial artery",
            "nerve_tibial": "tibial nerve", "adipose_omentum": "omental adipose", "adrenal_gland": "adrenal gland",
            "esophagus_ge_junction": "gastroesophageal junction", "esophagus_mucosa": "esophagus, mucosa",
            "esophagus_muscularis": "esophagus, muscularis", "mammary_tissue": "mammary tissue", "muscle": "skeletal muscle",
            "skin_sun_exposed": "skin, sun-exposed", "skin": "skin, unexposed", "small_intestine": "small intestine",
            "frontal_cortex": "frontal cortex"}
    R = pd.read_csv(D / "catlas_benchmark.csv")
    T = pd.read_csv(ROOT / "results" / "benchmark_aggregation.csv"); T = T[T.model.isin(["proximal", "enhancer"])]
    TD = pd.read_csv(ROOT / "results" / "benchmark_depth_by_type.csv")
    PB = pd.read_csv(ADD / "multiome" / "multiome_benchmark.csv")
    BR = pd.read_csv(D / "brain3k_benchmark.csv")
    BM = pd.read_csv(D / "bmmc_benchmark.csv"); LD = pd.read_csv(D / "bmmc_lodo.csv")
    MSET = ["raw_mean", "lognorm_mean", "cpm_mean"]

    def get(df, method, col, **kw):
        s = df[df.method == method]
        for kk, vv in kw.items(): s = s[s[kk] == vv]
        return float(s[col].iloc[0])

    H = HEIGHT[4]
    fig = plt.figure(figsize=(W, H), layout="constrained")
    fig.get_layout_engine().set(w_pad=0.04, h_pad=0.05, wspace=0.05, hspace=0.05)
    gs = fig.add_gridspec(2, 2, height_ratios=[2.25, 1.0], width_ratios=[1.12, 1.0])
    order = R[(R.model == "proximal") & (R.method == "raw_mean")].sort_values("coverage_spread").organ.tolist()
    for col, (model, title) in enumerate([("proximal", f"A  Count matrix, {R.organ.nunique()} non-tumor tissues"),
                                          ("enhancer", "B  Weighted matrix")]):
        ax = fig.add_subplot(gs[0, col])
        S = R[R.model == model].set_index(["organ", "method"])
        for y, organ in enumerate(order):
            vals = [S.loc[(organ, m), "top_share"] * 100 for m in MSET]
            ax.plot([min(vals), max(vals)], [y, y], color="#cfcfcf", lw=1.1, zorder=1)
            ax.plot([S.loc[(organ, "raw_mean"), "uniform_share"] * 100], [y], marker="|", color="#444444",
                    ms=6, mew=1.0, zorder=2)
            for m, v in zip(MSET, vals):
                ax.plot([v], [y], MRK[m], color=COL[m], ms=5, mec="white", mew=0.5, zorder=3,
                        label=LAB[m] if y == 0 else None)
        ax.set_yticks(range(len(order)))
        if col == 0:
            ax.set_yticklabels([f"{NICE.get(o, o.replace('_', ' '))} ({S.loc[(o, 'raw_mean'), 'coverage_spread']:.1f}×, "
                                f"{int(S.loc[(o, 'raw_mean'), 'n_types'])})" for o in order], fontsize=7.2)
        else:
            ax.set_yticklabels([])
        ax.set_ylim(-0.7, len(order) - 0.3)
        ax.set_xlabel("top share (%)")
        ax.set_title(title)
        ax.grid(True, axis="x", **GRID); ax.set_axisbelow(True)
        if col == 1:
            h, l = ax.get_legend_handles_labels()
            h.append(Line2D([], [], marker="|", color="#444444", ls="", ms=6, mew=1.0)); l.append("uniform share, 1/k")
            fig.legend(h, l, loc="outside upper center", ncol=4, handletextpad=0.3, columnspacing=1.6)

    # C
    ax = fig.add_subplot(gs[1, 0])
    pts = []
    for model in ("proximal", "enhancer"):
        s = R[R.model == model].pivot(index="organ", columns="method", values="top_share")
        sp_ = R[(R.model == model) & (R.method == "raw_mean")].set_index("organ").coverage_spread
        for o in s.index:
            pts.append(dict(group="non-tumor tissue", model=model, spread=sp_[o], ratio=s.loc[o, "raw_mean"] / s.loc[o, "cpm_mean"]))
        for c in T.cohort.unique():
            d = TD[(TD.cohort == c) & (TD.model == model)]
            pts.append(dict(group="cancer-associated", model=model, spread=d.median_true_depth.max() / d.median_true_depth.min(),
                            ratio=get(T, "raw_mean", "top_share", cohort=c, model=model) / get(T, "cpm_mean", "top_share", cohort=c, model=model)))
        for df in (PB, BR):
            pts.append(dict(group="matched RNA–ATAC", model=model, spread=get(df, "raw_mean", "depth_spread", model=model),
                            ratio=get(df, "raw_mean", "top_share", model=model) / get(df, "cpm_mean", "top_share", model=model)))
    pts.append(dict(group="matched RNA–ATAC", model="deposited", spread=get(BM, "raw_mean", "coverage_spread"),
                    ratio=get(BM, "raw_mean", "top_share") / get(BM, "cpm_mean", "top_share")))
    P = pd.DataFrame(pts)
    GC = {"non-tumor tissue": "#7a7a7a", "cancer-associated": "#b2182b", "matched RNA–ATAC": "#2166ac"}
    for g, c in GC.items():
        for model, mk in (("proximal", "o"), ("enhancer", "s"), ("deposited", "D")):
            d = P[(P.group == g) & (P.model == model)]
            if len(d):
                ax.scatter(d.spread, d.ratio, s=20, marker=mk, color=c, edgecolor="white", linewidth=0.4, zorder=3, alpha=0.9)
    ax.axhline(1, color="#333333", lw=0.7, ls="--")
    ax.set_xscale("log"); ax.set_xticks([1, 2, 5, 10, 20, 40]); ax.set_xticklabels(["1", "2", "5", "10", "20", "40"])
    ax.tick_params(axis="x", which="minor", bottom=False)
    ax.set_ylim(0.6, 7.4)
    ax.set_xlabel("coverage spread between cell types (fold)")
    ax.set_ylabel("raw-mean / CP10K top share")
    rho = spearmanr(P.spread, P.ratio)[0]
    ax.set_title("C  Raw-mean excess")
    ax.text(0.97, 0.12, f"Spearman ρ = {rho:+.2f}\n{len(P)} matrices", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7.5)
    h = [Line2D([], [], marker="o", color=c, ls="", ms=5) for c in GC.values()]
    l = ["atlas tissue", "HCC/COAD/NSCLC", "RNA–ATAC"]
    for mk, nm in (("o", "count"), ("s", "weighted"), ("D", "deposited")):
        h.append(Line2D([], [], marker=mk, color="#333333", ls="", ms=4.5, mfc="none", mew=0.8)); l.append(nm)
    ax.legend(h, l, loc="upper left", ncol=2, handletextpad=0.2, columnspacing=0.8, labelspacing=0.3,
              borderaxespad=0.1, fontsize=7.2)
    ax.grid(True, **GRID); ax.set_axisbelow(True)

    # D
    ax = fig.add_subplot(gs[1, 1])
    sets = [("PBMC\n1 donor", lambda m: np.mean([get(PB, m, "argmax_agree_informative", model=k) for k in ("proximal", "enhancer")]), None),
            ("brain\n1 donor", lambda m: np.mean([get(BR, m, "argmax_agree_informative", model=k) for k in ("proximal", "enhancer")]), None),
            ("bone marrow\n10 donors", lambda m: get(BM, m, "argmax_agree_informative"), LD)]
    w = 0.27
    for i, (name, f, lodo) in enumerate(sets):
        for j, m in enumerate(MSET):
            xx = i + (j - 1) * w
            ax.bar(xx, f(m) * 100, width=w * 0.92, color=COL[m], label=LAB[m] if i == 0 else None, zorder=2)
            if lodo is not None:
                r = lodo[lodo.method == m].argmax_agree_informative * 100
                ax.plot([xx, xx], [r.min(), r.max()], color="#111111", lw=1.0, zorder=3)
    ax.set_xticks(range(len(sets))); ax.set_xticklabels([s[0] for s in sets])
    ax.set_ylabel("per-gene maximum agrees with RNA (%)")
    ax.set_title("D  Agreement with RNA")
    ax.legend(loc="upper left", handletextpad=0.4, ncol=1, handlelength=1.2, fontsize=7.2, borderaxespad=0.2)
    ax.grid(True, axis="y", **GRID); ax.set_axisbelow(True); ax.set_ylim(0, 100)
    save(fig, "fig4_normal_tissues.png", H)
    print(f"fig4: panel C Spearman rho {rho:+.3f} over {len(P)} matrices")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4()
