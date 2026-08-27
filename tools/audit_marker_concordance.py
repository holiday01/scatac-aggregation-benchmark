"""
Test 1: does each canonical lineage marker have its accessibility maximal in the
cell type it is a marker of?

This is the cheapest check on a gene-activity prior, and the one that catches a
guessed row annotation. It needs only the prior, its cell labels, and a marker
panel with an externally known answer.

Two properties of the scoring rule are worth knowing before reading the output.
The panel targets four lineages, but the HCC cohort carries six annotated cell
types (NK_cytotoxic_T, CD4_T, B_cell, DC, Hepatocyte, Macrophage), and a marker
counts as concordant only if its argmax is the ONE designated type. IL7R is the
case where that bites: its argmax is CD4_T, which is lineage-consistent but is
scored as a failure because the T/NK panel is keyed to NK_cytotoxic_T. The rule
is deliberately strict and identical across the priors being compared.

Run:  python3 tools/audit_marker_concordance.py
Out:  results/marker_concordance.csv + a console summary, per prior and per lineage
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
OUT = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(parents=True, exist_ok=True)

PRIORS = {
    "old_guessed_names": "data/processed/integrated/atac_hcc_geneactivity.h5ad",
    "new_proximal": "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad",
    "new_enhancer": "data/processed/integrated/atac_hcc_geneactivity_enhancer.h5ad",
}

# marker -> (lineage label, the one cell type its accessibility must peak in)
PANEL = {
    "CD68": ("Mac", "Macrophage"), "CSF1R": ("Mac", "Macrophage"),
    "MRC1": ("Mac", "Macrophage"), "MSR1": ("Mac", "Macrophage"),
    "CD163": ("Mac", "Macrophage"), "ITGAM": ("Mac", "Macrophage"),
    "MS4A1": ("B", "B_cell"), "CD79A": ("B", "B_cell"),
    "CD79B": ("B", "B_cell"), "CD19": ("B", "B_cell"),
    "BANK1": ("B", "B_cell"),
    "CD3D": ("T", "NK_cytotoxic_T"), "CD3E": ("T", "NK_cytotoxic_T"),
    "CD8A": ("T", "NK_cytotoxic_T"), "IL7R": ("T", "NK_cytotoxic_T"),
    "NKG7": ("T", "NK_cytotoxic_T"), "GNLY": ("T", "NK_cytotoxic_T"),
    "ALB": ("Hep", "Hepatocyte"), "APOB": ("Hep", "Hepatocyte"),
    "HNF4A": ("Hep", "Hepatocyte"), "TTR": ("Hep", "Hepatocyte"),
    "APOA1": ("Hep", "Hepatocyte"),
}


def celltype_means(path, normalise=True):
    """Mean gene activity per cell type.

    Cells are depth-normalised BEFORE the group mean. Without this the mean is a
    sum over raw counts, so a cell type sequenced 1.8x deeper contributes 1.8x the
    signal to its own column and wins the argmax for genes that have nothing to do
    with it: on the HCC priors the B-cell compartment (the deepest, median 13,139
    counts per cell) takes 75.6 % of all genome-wide argmaxima un-normalised
    against 15.4 % normalised. That is Mode 3 of Table 1 operating on this test's
    own intermediate, and it inverts the per-lineage result.

    The choice of normaliser does not matter; CPM, log1p(CPM), median-scaling and a
    depth-matched subsample with no normalisation at all agree to within one marker.
    `normalise=False` reproduces the superseded numbers for comparison.
    """
    a = ad.read_h5ad(path)
    X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
    X = X.astype(np.float64)
    if normalise:
        tot = X.sum(1, keepdims=True)
        tot[tot == 0] = 1.0
        X = X / tot * 1e4
    ct = a.obs["cell_type"].astype(str).values
    return pd.DataFrame({c: X[ct == c].mean(0) for c in pd.unique(ct)},
                        index=a.var_names.astype(str))


def main(normalise=True):
    out_csv = OUT / ("marker_concordance.csv" if normalise
                     else "marker_concordance_RAW_superseded.csv")
    print(f"cell-type profiles: {'depth-normalised (CPM)' if normalise else 'RAW COUNTS -- superseded, for comparison only'}\n")
    rows = []
    for tag, rel in PRIORS.items():
        path = BASE / rel
        if not path.exists():
            print(f"  {tag}: {path} not found, skipped")
            continue
        m = celltype_means(path, normalise=normalise)
        # Mode 3 check on this test's own intermediate: if one cell type takes a
        # large share of ALL argmaxima, the panel result is reporting depth.
        share = m.idxmax(axis=1).value_counts(normalize=True)
        print(f"  {tag}: top argmax share " +
              ", ".join(f"{k} {v:.1%}" for k, v in share.head(3).items()))
        print(f"{tag}: {m.shape[0]} genes x {m.shape[1]} cell types "
              f"({', '.join(sorted(m.columns))})")
        for gene, (lineage, target) in PANEL.items():
            if gene not in m.index or target not in m.columns:
                continue
            argmax = m.loc[gene].idxmax()
            rows.append({"prior": tag, "gene": gene, "lineage": lineage,
                         "target_cell_type": target, "atac_argmax": argmax,
                         "correct": int(argmax == target)})
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)

    print("\nconcordant markers, overall and per lineage:")
    tot = df.groupby("prior").correct.agg(["sum", "count"])
    per = df.pivot_table(index="lineage", columns="prior", values="correct",
                         aggfunc=["sum", "count"])
    order = [p for p in PRIORS if p in tot.index]
    print("  overall   " + "   ".join(
        f"{p}: {tot.loc[p, 'sum']}/{tot.loc[p, 'count']}" for p in order))
    for lin in ["Mac", "B", "T", "Hep"]:
        if lin not in per.index:
            continue
        print(f"  {lin:9} " + "   ".join(
            f"{p}: {int(per.loc[lin, ('sum', p)])}/{int(per.loc[lin, ('count', p)])}"
            for p in order))

    # what a single aggregate number hides: which markers move, and which way
    piv = df.pivot_table(index="gene", columns="prior", values="correct",
                         aggfunc="first")
    base = "old_guessed_names"
    for p in order:
        if p == base or base not in piv.columns:
            continue
        gained = sorted(g for g in piv.index if piv.loc[g, base] == 0 and piv.loc[g, p] == 1)
        lost = sorted(g for g in piv.index if piv.loc[g, base] == 1 and piv.loc[g, p] == 0)
        print(f"\n  {base} -> {p}")
        print(f"    gained ({len(gained)}): {', '.join(gained) or 'none'}")
        print(f"    lost   ({len(lost)}): {', '.join(lost) or 'none'}")

    # A permutation null turns the count into an instrument: shuffle which
    # marker each observed argmax belongs to, and see how often chance alone
    # reaches the observed concordance.
    rng = np.random.default_rng(42)
    print("\n  permutation null (10,000 shuffles of the marker-to-argmax assignment):")
    for pr in order:
        s = df[df.prior == pr]
        arg = s.atac_argmax.values
        tgt = s.target_cell_type.values
        null = np.array([(rng.permutation(arg) == tgt).sum() for _ in range(10000)])
        obs = int(s.correct.sum())
        print(f"    {pr:18} observed {obs}/{len(s)}   null mean {null.mean():.2f}"
              f"   95th pct {np.percentile(null, 95):.0f}   p={np.mean(null >= obs):.4f}")

    # A looser rule, in which a lineage may map to more than one annotated type,
    # as a sensitivity check on the strict one used above.
    loose_map = {"Mac": {"Macrophage"}, "B": {"B_cell"},
                 "T": {"NK_cytotoxic_T", "CD4_T"}, "Hep": {"Hepatocyte"}}
    df["loose"] = [int(r.atac_argmax in loose_map[r.lineage]) for r in df.itertuples()]
    print("\n  strict vs lineage-consistent scoring:")
    for pr in order:
        s = df[df.prior == pr]
        print(f"    {pr:18} strict {int(s.correct.sum())}/{len(s)}   "
              f"loose {int(s.loose.sum())}/{len(s)}")

    df.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")


if __name__ == "__main__":
    import sys
    main(normalise="--raw" not in sys.argv)
