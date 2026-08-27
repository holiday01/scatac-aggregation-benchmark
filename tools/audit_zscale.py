"""
Test 2: is this stored matrix data, or a transform?

Our own pipeline (65_coad_atac_full_pipeline.py, line 209) called
`sc.pp.scale(adata, max_value=10)` — a z-scaling step meant to prepare a matrix
for PCA and clustering — and then wrote the result to disk as the cohort's
gene-activity matrix. Everything downstream read it as accessibility.

Two fingerprints separate a stored transform from stored counts, and both are
computed here:

  * the share of negative entries — a count matrix cannot be mostly negative;
  * the per-gene mean and standard deviation — z-scores have mean 0 and unit
    variance. Note that `max_value=10` clips the upper tail, so the per-gene sd
    is 1 only up to that clipping: the median is 1.000 but the mean is lower.

The script also settles what the z-scaling does NOT explain. Per-gene z-scaling
is an affine map within each gene, so the argmax of a gene's per-cell-type means
is invariant to it. The heavy concentration of per-gene argmaxima on one cell
type is therefore a property of the underlying counts, not an artefact of the
transform — it tracks per-cell-type detection rate, which is Test 3's territory.
This script demonstrates that invariance rather than asserting it.

Run:  python3 tools/audit_zscale.py
Out:  results/zscale_diagnostic.csv, results/zscale_gene_argmax.csv + a summary
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
MATRIX = BASE / "data/processed/COAD/coad_atac_geneactivity_full.h5ad"
OUT = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(parents=True, exist_ok=True)


def _decode(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])


def main():
    with h5py.File(MATRIX, "r") as f:
        X = f["X"][:]                                   # dense, as written
        ct = f["obs/cell_type_broad"]
        cats = _decode(ct["categories"][:])
        labels = cats[ct["codes"][:]]
        vkey = f["var"].attrs.get("_index", "_index")
        vkey = vkey.decode() if isinstance(vkey, bytes) else vkey
        genes = _decode(f["var"][vkey][:])

    n_cells, n_genes = X.shape
    neg = float((X < 0).mean())
    sd = X.std(axis=0)
    mu = X.mean(axis=0)
    clipped = float((X >= 10 - 1e-9).mean())

    print(f"matrix {n_cells} cells x {n_genes} genes")
    print(f"  negative entries : {100 * neg:.2f} %  (a count matrix cannot be)")
    print(f"  min / mean / max : {X.min():.3f} / {X.mean():.3f} / {X.max():.3f}")
    print(f"  per-gene mean    : max |mean| {np.abs(mu).max():.4f}  (z-scores: 0)")
    print(f"  per-gene sd      : min {sd.min():.6f}  median {np.median(sd):.6f}  "
          f"mean {sd.mean():.6f}  (z-scores: 1)")
    print(f"                     {100 * np.mean(np.abs(sd - 1) < 1e-3):.1f} % within "
          f"1e-3 of 1; the rest are pulled down by the max_value=10 clip")
    print(f"  entries clipped at +10 : {100 * clipped:.4f} %")

    # Per-gene argmax over cell types, before and after undoing the per-gene
    # shift. If they agree, the transform cannot be what concentrates them.
    means_scaled = pd.DataFrame(X, columns=genes).groupby(labels).mean()
    means_shift = pd.DataFrame(X - X.min(axis=0), columns=genes).groupby(labels).mean()
    a_scaled = means_scaled.idxmax(axis=0)
    a_shift = means_shift.idxmax(axis=0)
    agree = bool((a_scaled == a_shift).all())
    share = a_scaled.value_counts(normalize=True).sort_values(ascending=False)

    print(f"\nper-gene argmax over cell types is identical before and after "
          f"undoing the per-gene shift: {agree}")
    print("  share of the panel taking its maximum in each cell type:")
    for cty, frac in share.items():
        print(f"    {cty:15s} {100 * frac:5.1f} %")

    detect = pd.Series((X > X.min(axis=0)).mean(axis=1)).groupby(labels).mean()
    print("\n  mean fraction of the panel detected per cell (the actual driver):")
    for cty, frac in detect.sort_values(ascending=False).items():
        print(f"    {cty:15s} {100 * frac:5.2f} %")

    pd.DataFrame({
        "metric": ["n_cells", "n_genes", "pct_entries_negative",
                   "per_gene_mean_max_abs", "per_gene_sd_min", "per_gene_sd_median",
                   "per_gene_sd_mean", "pct_genes_sd_within_1e-3_of_1",
                   "pct_entries_clipped_at_10", "argmax_invariant_to_shift",
                   "top_argmax_cell_type", "pct_genes_at_top_argmax"],
        "value": [n_cells, n_genes, round(100 * neg, 2),
                  round(float(np.abs(mu).max()), 4), round(float(sd.min()), 6),
                  round(float(np.median(sd)), 6), round(float(sd.mean()), 6),
                  round(100 * float(np.mean(np.abs(sd - 1) < 1e-3)), 1),
                  round(100 * clipped, 4), agree,
                  share.index[0], round(100 * float(share.iloc[0]), 2)],
    }).to_csv(OUT / "zscale_diagnostic.csv", index=False)
    pd.DataFrame({"argmax_cell_type": a_scaled,
                  "argmax_after_unshift": a_shift}).to_csv(
        OUT / "zscale_gene_argmax.csv", index_label="gene")
    print(f"\nwrote {OUT / 'zscale_diagnostic.csv'} and zscale_gene_argmax.csv")


if __name__ == "__main__":
    main()
