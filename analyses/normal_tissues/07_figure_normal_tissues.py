"""
07: Figure for the non-tumor extension (style constants follow tools/make_revision_figures.py).

  A  top share per tissue, count matrix: raw / log-normalized / CP10K, tissues ordered by their
     natural coverage spread, with the uniform share 1/k marked
  B  the same for the distance-weighted matrix
  C  inflation of the raw mean (raw top share / CP10K top share) against coverage spread for every
     matrix in the study: the non-tumor tissues, three tumor datasets, three matched RNA-ATAC samples
  D  agreement of the per-gene maximum with the RNA reference (RNA-informative genes) in the three
     matched RNA-ATAC samples; bone-marrow bars carry the leave-one-donor-out range

Out: fig_normal_tissues.png / .pdf, fig_normal_tissues_data.csv (the points of panel C)
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TEXTW, DPI = 6.3, 600
plt.rcParams.update({"font.size": 7.5, "font.family": "DejaVu Sans", "axes.linewidth": 0.7,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7})
MSET = [("raw_mean", "raw mean", "#b2182b"), ("lognorm_mean", "log-normalized mean", "#2f8f5b"),
        ("cpm_mean", "CP10K mean", "#2166ac")]
NICE = {"colon_transverse": "colon, transverse", "colon_sigmoid": "colon, sigmoid", "heart_lv": "heart, left ventricle",
        "heart_atrial_appendage": "heart, atrial appendage", "artery_aorta": "aorta", "artery_tibial": "tibial artery",
        "nerve_tibial": "tibial nerve", "adipose_omentum": "omental adipose", "adrenal_gland": "adrenal gland",
        "esophagus_ge_junction": "gastroesophageal junction", "esophagus_mucosa": "esophagus, mucosa",
        "esophagus_muscularis": "esophagus, muscularis", "mammary_tissue": "mammary tissue", "muscle": "skeletal muscle",
        "skin_sun_exposed": "skin, sun-exposed", "skin": "skin, unexposed", "small_intestine": "small intestine",
        "frontal_cortex": "frontal cortex"}

R = pd.read_csv(HERE / "catlas_benchmark.csv")
T = pd.read_csv(ROOT / "results" / "benchmark_aggregation.csv"); T = T[T.model.isin(["proximal", "enhancer"])]
TD = pd.read_csv(ROOT / "results" / "benchmark_depth_by_type.csv")
PB = pd.read_csv(HERE.parent / "multiome" / "multiome_benchmark.csv")
BR = pd.read_csv(HERE / "brain3k_benchmark.csv")
BM = pd.read_csv(HERE / "bmmc_benchmark.csv"); LD = pd.read_csv(HERE / "bmmc_lodo.csv")


def get(df, method, col, **kw):
    s = df[df.method == method]
    for k, v in kw.items(): s = s[s[k] == v]
    return float(s[col].iloc[0])


fig = plt.figure(figsize=(TEXTW, 7.4), layout="constrained")
gs = fig.add_gridspec(2, 2, height_ratios=[2.15, 1.0])

# ── A, B: per-tissue top share ───────────────────────────────────────────────
order = (R[(R.model == "proximal") & (R.method == "raw_mean")].sort_values("coverage_spread").organ.tolist())
for col, (model, title) in enumerate([("proximal", f"Count matrix, {R.organ.nunique()} non-tumor tissues"), ("enhancer", "Distance-weighted matrix")]):
    ax = fig.add_subplot(gs[0, col])
    S = R[R.model == model].set_index(["organ", "method"])
    for y, organ in enumerate(order):
        vals = [S.loc[(organ, m), "top_share"] * 100 for m, _, _ in MSET]
        ax.plot([min(vals), max(vals)], [y, y], color="#cfcfcf", lw=1.0, zorder=1)
        ax.plot([S.loc[(organ, "raw_mean"), "uniform_share"] * 100], [y], marker="|", color="#555555", ms=5, mew=0.9, zorder=2)
        for (m, lab_, c), v in zip(MSET, vals):
            ax.plot([v], [y], "o", color=c, ms=3.6, mec="white", mew=0.5, zorder=3, label=lab_ if y == 0 else None)
    ax.set_yticks(range(len(order)))
    if col == 0:
        ax.set_yticklabels([f"{NICE.get(o, o.replace('_', ' '))} ({S.loc[(o, 'raw_mean'), 'coverage_spread']:.1f}x, "
                            f"{int(S.loc[(o, 'raw_mean'), 'n_types'])})" for o in order], fontsize=6.0)
    else:
        ax.set_yticklabels([])
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel("largest share of per-gene maxima (%)", fontsize=7)
    ax.set_title(f"{'AB'[col]}  {title}", fontsize=7.6, loc="left")
    ax.tick_params(labelsize=6.3, length=2); ax.grid(True, axis="x", lw=0.4, color="#e6e6e6"); ax.set_axisbelow(True)
    if col == 1:
        h, l = ax.get_legend_handles_labels()
        h.append(plt.Line2D([], [], marker="|", color="#555555", ls="", ms=5, mew=0.9)); l.append("uniform share, 1/k")
        fig.legend(h, l, fontsize=6.4, frameon=False, loc="outside upper center", ncol=4, handletextpad=0.3, columnspacing=1.6)

# ── C: inflation of the raw mean against coverage spread ─────────────────────
ax = fig.add_subplot(gs[1, 0])
pts = []
for model, mk in (("proximal", "o"), ("enhancer", "s")):
    s = R[R.model == model].pivot(index="organ", columns="method", values="top_share")
    sp_ = R[(R.model == model) & (R.method == "raw_mean")].set_index("organ").coverage_spread
    for o in s.index:
        pts.append(dict(group="non-tumor tissue", name=o, model=model, spread=sp_[o], ratio=s.loc[o, "raw_mean"] / s.loc[o, "cpm_mean"]))
    for c in T.cohort.unique():
        d = TD[(TD.cohort == c) & (TD.model == model)]
        pts.append(dict(group="tumor dataset", name=c, model=model, spread=d.median_true_depth.max() / d.median_true_depth.min(),
                        ratio=get(T, "raw_mean", "top_share", cohort=c, model=model) / get(T, "cpm_mean", "top_share", cohort=c, model=model)))
    for nm, df, scol in (("PBMC", PB, "depth_spread"), ("brain", BR, "depth_spread")):
        pts.append(dict(group="matched RNA-ATAC", name=nm, model=model, spread=get(df, "raw_mean", scol, model=model),
                        ratio=get(df, "raw_mean", "top_share", model=model) / get(df, "cpm_mean", "top_share", model=model)))
pts.append(dict(group="matched RNA-ATAC", name="bone marrow", model="deposited", spread=get(BM, "raw_mean", "coverage_spread"),
                ratio=get(BM, "raw_mean", "top_share") / get(BM, "cpm_mean", "top_share")))
P = pd.DataFrame(pts); P.to_csv(HERE / "fig_normal_tissues_data.csv", index=False)
GC = {"non-tumor tissue": "#7a7a7a", "tumor dataset": "#b2182b", "matched RNA-ATAC": "#2166ac"}
for g, c in GC.items():
    for model, mk in (("proximal", "o"), ("enhancer", "s"), ("deposited", "D")):
        d = P[(P.group == g) & (P.model == model)]
        if len(d): ax.scatter(d.spread, d.ratio, s=13, marker=mk, color=c, edgecolor="white", linewidth=0.4, zorder=3,
                              alpha=0.9, label=g if model == "proximal" else None)
ax.axhline(1, color="#333333", lw=0.7, ls="--")
ax.set_xscale("log"); ax.set_xticks([1, 2, 5, 10, 20, 40]); ax.set_xticklabels(["1", "2", "5", "10", "20", "40"])
ax.set_xlabel("coverage spread between cell types (fold)", fontsize=7)
ax.set_ylabel("raw-mean top share / CP10K top share", fontsize=7)
rho = spearmanr(P.spread, P.ratio)[0]
ax.set_title("C  Raw-mean inflation vs coverage spread", fontsize=7.6, loc="left")
ax.text(0.98, 0.20, rf"Spearman $\rho$ = {rho:+.2f}" + f"\n{len(P)} matrices", transform=ax.transAxes, fontsize=6.3, ha="right", va="bottom", linespacing=1.25)
h, l = ax.get_legend_handles_labels()
for mk, nm in (("o", "count matrix"), ("s", "distance-weighted"), ("D", "deposited matrix")):
    h.append(plt.Line2D([], [], marker=mk, color="#333333", ls="", ms=3.2, mfc="none", mew=0.7)); l.append(nm)
ax.legend(h, l, fontsize=5.9, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.84), handletextpad=0.2, labelspacing=0.35, borderaxespad=0.1)
ax.tick_params(labelsize=6.5, length=2); ax.grid(True, lw=0.4, color="#e6e6e6"); ax.set_axisbelow(True)

# ── D: agreement with the RNA reference ──────────────────────────────────────
ax = fig.add_subplot(gs[1, 1])
sets = [("PBMC\n1 donor", lambda m: np.mean([get(PB, m, "argmax_agree_informative", model=k) for k in ("proximal", "enhancer")]), None),
        ("brain\n1 donor", lambda m: np.mean([get(BR, m, "argmax_agree_informative", model=k) for k in ("proximal", "enhancer")]), None),
        ("bone marrow\n10 donors", lambda m: get(BM, m, "argmax_agree_informative"), LD)]
w = 0.26
for i, (name, f, lodo) in enumerate(sets):
    for j, (m, lab_, c) in enumerate(MSET):
        x = i + (j - 1) * w; v = f(m) * 100
        ax.bar(x, v, width=w * 0.92, color=c, label=lab_ if i == 0 else None, zorder=2)
        if lodo is not None:
            r = lodo[lodo.method == m].argmax_agree_informative * 100
            ax.plot([x, x], [r.min(), r.max()], color="#222222", lw=0.9, zorder=3)
ax.set_xticks(range(len(sets))); ax.set_xticklabels([s[0] for s in sets], fontsize=6.4)
ax.set_ylabel("per-gene maximum agrees with RNA (%)", fontsize=7)
ax.set_title("D  Agreement with RNA reference", fontsize=7.6, loc="left")
ax.legend(fontsize=6.0, frameon=False, loc="upper left", handletextpad=0.3, ncol=1)
ax.tick_params(labelsize=6.5, length=2); ax.grid(True, axis="y", lw=0.4, color="#e6e6e6"); ax.set_axisbelow(True)
ax.set_ylim(0, 100)

fig.savefig(HERE / "fig_normal_tissues.png", dpi=DPI, facecolor="white")
fig.savefig(HERE / "fig_normal_tissues.pdf", facecolor="white")
print("saved fig_normal_tissues.png/.pdf;  panel C Spearman rho = %+.3f over %d matrices" % (rho, len(P)))
print(P.groupby("group").ratio.describe()[["count", "min", "50%", "max"]].round(2).to_string())
