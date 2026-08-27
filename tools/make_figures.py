"""
Figure for the APBC2026 paper on silent failure modes of chromatin gene-activity
priors.

Figure 1  (A) marker-lineage concordance for all 22 canonical markers under the
              positionally-named prior and the two priors rebuilt with correct
              gene identities;
          (B) what that costs downstream: per-combination concordance gain on the
              six HCC/HNSC combinations with everything held constant except the
              prior, 1/6 -> 4/6, with the sign reversing in 5 of 6.

Drawn at exactly the printed size (6.3 in = 16 cm, the text width of the APBC
a4paper/2.5 cm template) and saved WITHOUT bbox_inches="tight", so the canvas
LaTeX places at \\textwidth is the canvas matplotlib drew: the point sizes below
are the point sizes on the page, and 600 dpi is 600 dpi. A tight bounding box
crops the canvas and lets \\includegraphics silently upscale it.

Panel B's x-limits are derived from the data and asserted, never hard-coded, so
no confidence interval can be clipped.

Inputs : results/marker_concordance.csv                       (Test 1, normalised)
         results/validation/corrected_prior_headtohead.csv        (Script 91)
Output : figures/fig1_concordance.png (600 dpi)
"""
import os
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from PIL import Image
from scipy import stats

BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
OUT = Path(__file__).resolve().parent.parent / "figures"
OUT.mkdir(parents=True, exist_ok=True)

TEXTW = 6.3   # inches == 16 cm, the template's text width
DPI = 600

plt.rcParams.update({
    "font.size": 8,
    "font.family": "DejaVu Sans",
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})

C_BAD = "#c0392b"
C_PROX = "#2c7fb8"
C_ENH = "#238b45"

# Colour-blind- and greyscale-safe pair: pale grey versus deep blue separates on
# luminance rather than hue, so it survives every dichromacy and a monochrome
# printer. The glyphs remain the primary encoding; the fill is redundant.
CMAP_OK = ListedColormap(["#eeeeee", "#3f6fb0"])

PRIORS = ["old_guessed_names", "new_proximal", "new_enhancer"]
PRIOR_LABEL = ["positional\nnames", "rebuilt\nproximal", "rebuilt\nenhancer"]
HEAD = {"old_guessed": "positional names", "new_proximal": "rebuilt proximal",
        "new_enhancer": "rebuilt enhancer"}


def _read(name, subdir):
    for q in (BASE / subdir / name,
              Path(__file__).resolve().parent.parent / "data" / name):
        if q.exists():
            return pd.read_csv(q)
    raise FileNotFoundError(f"{name} not found under {BASE / subdir} or ../data")


def save(fig, name):
    """Save at the delivered size, and prove it."""
    fig.savefig(OUT / name, dpi=DPI, facecolor="white",
                bbox_inches=None, pad_inches=0)
    plt.close(fig)
    want = round(TEXTW * DPI)
    got = Image.open(OUT / name).size
    assert got[0] == want, f"{name}: width {got[0]} px, expected {want}"
    print(f"saved {name}  {got[0]}x{got[1]} px  ({got[0] / TEXTW:.0f} dpi on the page)")


def figure1():
    # Panel A comes from this paper's own test script, which depth-normalises each
    # cell before the cell-type mean. Script 90's gene_activity_rebuild_validation
    # .csv is the same schema computed WITHOUT that step and is superseded; see
    # STATUS_NUMBERS_SUPERSEDED.md.
    val = pd.read_csv(Path(__file__).resolve().parent.parent
                      / "results" / "marker_concordance.csv")
    head = _read("corrected_prior_headtohead.csv", "results/validation")

    fig = plt.figure(figsize=(TEXTW, 4.35), layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=[0.86, 1.24], wspace=0.06)

    # ── A: every canonical marker under every prior ──────────────────────
    axA = fig.add_subplot(gs[0, 0])
    piv = val.pivot_table(index="gene", columns="prior", values="correct",
                          aggfunc="first")
    lin = val.drop_duplicates("gene").set_index("gene").lineage
    groups = ["Mac", "B", "T", "Hep"]
    gnames = {"Mac": "macrophage", "B": "B cell", "T": "T / NK",
              "Hep": "hepatocyte"}
    genes, sep = [], []
    for g in groups:
        genes += sorted(x for x in piv.index if lin.get(x) == g)
        sep.append(len(genes))
    mat = piv.loc[genes, PRIORS].values.astype(float)

    axA.imshow(mat, cmap=CMAP_OK, vmin=0, vmax=1, aspect="auto")
    axA.set_xticks(range(3))
    axA.set_xticklabels(PRIOR_LABEL, fontsize=6.8)
    axA.set_yticks(range(len(genes)))
    axA.set_yticklabels(genes, fontsize=6.5)
    for i in range(len(genes)):
        for j in range(3):
            hit = mat[i, j] == 1
            axA.text(j, i, "✓" if hit else "✗", ha="center", va="center",
                     fontsize=7.5, color="white" if hit else "#333333")
    for s in sep[:-1]:
        axA.axhline(s - 0.5, color="white", lw=2.2)
    start = 0
    for g, end in zip(groups, sep):
        axA.text(-1.45, (start + end - 1) / 2, gnames[g], rotation=90,
                 ha="center", va="center", fontsize=6.3, color="#555555")
        start = end
    score = [int(val[val.prior == p].correct.sum()) for p in PRIORS]
    tot = int(val[val.prior == PRIORS[0]].shape[0])
    for j, s in enumerate(score):
        axA.text(j, len(genes) - 0.25, f"{s}/{tot}", ha="center", va="top",
                 fontsize=7.8, weight="bold", color=[C_BAD, C_PROX, C_ENH][j])
    axA.set_ylim(len(genes) + 0.45, -0.5)
    axA.set_title("A   Is each marker's accessibility\n"
                  "      maximal in its own lineage?",
                  fontsize=7.9, loc="left", weight="bold")
    axA.tick_params(length=0)
    for s in axA.spines.values():
        s.set_visible(False)

    # ── B: what the naming error costs downstream ────────────────────────
    axB = fig.add_subplot(gs[0, 1])
    head = head.copy()
    head["label"] = (head.cancer + " "
                     + head.cell_type.replace({"Macrophage": "Mac",
                                               "B_cell": "B",
                                               "B_Plasma": "B/plasma",
                                               "CD8_T": "CD8 T"}))
    lo = head.ci_lo.min()
    hi = head.ci_hi.max()
    span = hi - lo
    xmin, xmax = lo - 0.07 * span, hi + 0.04 * span

    y, yticks, ylabs, rows, ann_top = 0, [], [], [], 0.0
    for prior, colour in [("old_guessed", C_BAD), ("new_proximal", C_PROX),
                          ("new_enhancer", C_ENH)]:
        blk = head[head.prior == prior].sort_values("dAbsC")
        y0 = y
        for _, r in blk.iterrows():
            axB.plot([r.ci_lo, r.ci_hi], [y, y], color=colour, lw=1.2,
                     alpha=0.75, zorder=1, solid_capstyle="butt")
            axB.scatter(r.dAbsC, y, s=17, color=colour, zorder=3,
                        edgecolor="white", linewidth=0.4)
            yticks.append(y)
            ylabs.append(r.label)
            y += 1
        rows.append((y0 - 0.6, y - 0.4))
        d = blk.dAbsC.values
        w = stats.wilcoxon(d, alternative="greater", zero_method="wilcox").pvalue
        axB.axhspan(y - 0.35, y + 0.95, color=colour, alpha=0.07)
        axB.text(xmin + 0.002, y + 0.30,
                 f"{HEAD[prior]}: {int((d > 0).sum())}/{len(d)} favour, "
                 f"mean {d.mean():+.4f}, $p$={w:.2f}",
                 fontsize=6.3, weight="bold", color=colour, va="center", zorder=5,
                 bbox=dict(facecolor="white", edgecolor="none", pad=1.2))
        ann_top = max(ann_top, y + 0.95)
        y += 2.1

    for y0, y1 in rows:
        axB.vlines(0, y0, y1, color="black", lw=0.8, zorder=0)

    axB.set_yticks(yticks)
    axB.set_yticklabels(ylabs, fontsize=6.4)
    axB.set_xlabel(r"$\Delta|C-0.5|$    chromatin-informed $-$ transcription-only",
                   fontsize=7.4)
    axB.set_xlim(xmin, xmax)
    axB.tick_params(axis="x", labelsize=6.8)
    axB.set_ylim(-1.0, ann_top + 0.4)
    axB.set_title("B   Swapping only the prior reverses the\n"
                  "      sign in 5 of 6 combinations",
                  fontsize=7.9, loc="left", weight="bold")
    axB.spines[["top", "right"]].set_visible(False)

    assert xmin <= lo and xmax >= hi, "a confidence interval would be clipped"

    save(fig, "fig1_concordance.png")
    print(f"  concordance {dict(zip(PRIORS, score))} of {tot}; "
          f"panel B x-range [{xmin:.3f}, {xmax:.3f}] covers all CIs "
          f"[{lo:.3f}, {hi:.3f}]")


if __name__ == "__main__":
    figure1()
