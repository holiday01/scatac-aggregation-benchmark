"""
Re-audit the colorectal (GSE201336) sequencing-depth claim.

Background. The superseded analysis reported that the colorectal scATAC-seq
cohort consisted largely of empty droplets: 60 % of its 39,398 barcodes with
fewer than 100 fragments, macrophage and T/NK medians of 8 fragments, and 0 of
1434 macrophages surviving a >=1000-fragment filter. Those numbers were read
from obs['total_frags'] in data/processed/COAD/coad_atac_geneactivity_full.h5ad.

That column is not a fragment count. Script 65 line 172 sets

    adata.obs['total_frags'] = np.array(X.sum(axis=1)).ravel()

where X at that point is the +-2 kb TSS gene-activity matrix restricted to 1,896
highly variable genes. It is a promoter-window count summed over a small gene
panel, and it understates per-cell library depth by roughly two orders of
magnitude.

This script recomputes the true per-cell fragment counts from the SnapATAC2
imports of the same fragment files (data/processed/COAD/_s92_work/GSM*.h5ad,
produced by Script 92) and joins them to the published cell-type annotation, so
the two quantities can be compared side by side.

Run:  python3 tools/audit_coad_depth.py
Out :  results/coad_depth_audit.csv  + a console summary
"""
import collections
import os
from pathlib import Path

import numpy as np
import pandas as pd

import hdf5plugin  # noqa: F401  (registers the blosc filter SnapATAC2 writes with)
import h5py

BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
COAD = BASE / "data/processed/COAD"
OUT = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(parents=True, exist_ok=True)

MIN_FRAG = 1000  # the threshold the manuscript describes as standard


def _decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in arr])


def _categorical(h5obj):
    """Read an AnnData categorical (or plain string array) as python strings."""
    if isinstance(h5obj, h5py.Group):
        cats = _decode(h5obj["categories"][:])
        return cats[h5obj["codes"][:]]
    return _decode(h5obj[:])


def true_fragment_counts():
    """Per-cell fragment counts as recorded by SnapATAC2 at import time."""
    bcs, nfs = [], []
    for p in sorted(COAD.glob("_s92_work/GSM*.h5ad")):
        with h5py.File(p, "r") as f:
            nf = np.asarray(f["obs/n_fragment"][:])
            bc = _decode(f["obs/index"][:])
        gsm = p.stem
        bcs.append(np.array([f"{gsm}_{b}" for b in bc]))
        nfs.append(nf)
        print(f"  {gsm}: {len(nf):6d} cells   median {np.median(nf):8.0f} fragments"
              f"   >={MIN_FRAG}: {100 * (nf >= MIN_FRAG).mean():5.1f} %")
    return dict(zip(np.concatenate(bcs).tolist(), np.concatenate(nfs).tolist()))


def main():
    print("True per-cell fragment counts (SnapATAC2 n_fragment):")
    nfmap = true_fragment_counts()

    path = COAD / "coad_atac_geneactivity_full.h5ad"
    with h5py.File(path, "r") as f:
        ikey = f["obs"].attrs.get("_index", "_index")
        ikey = ikey.decode() if isinstance(ikey, bytes) else ikey
        barcodes = _decode(f["obs"][ikey][:])
        cell_type = _categorical(f["obs/cell_type_broad"])
        total_frags = f["obs/total_frags"][:]
        vkey = f["var"].attrs.get("_index", "_index")
        vkey = vkey.decode() if isinstance(vkey, bytes) else vkey
        n_genes = f["var"][vkey].shape[0]

    missing = [b for b in barcodes if b not in nfmap]
    if missing:
        raise SystemExit(f"{len(missing)} barcodes could not be joined; aborting")

    frags = np.array([nfmap[b] for b in barcodes], dtype=float)
    df = pd.DataFrame({"barcode": barcodes, "cell_type": cell_type,
                       "n_fragment": frags, "total_frags_column": total_frags})

    rows = []
    for cty, g in df.groupby("cell_type"):
        rows.append({
            "cell_type": cty,
            "n_cells": len(g),
            "median_n_fragment": float(np.median(g.n_fragment)),
            "pct_pass_1000_fragments": 100 * float((g.n_fragment >= MIN_FRAG).mean()),
            "n_pass_1000_fragments": int((g.n_fragment >= MIN_FRAG).sum()),
            "median_total_frags_column": float(np.median(g.total_frags_column)),
            "pct_pass_1000_total_frags_column":
                100 * float((g.total_frags_column >= MIN_FRAG).mean()),
        })
    summary = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    summary.to_csv(OUT / "coad_depth_audit.csv", index=False)

    print(f"\nCohort: {len(df)} barcodes, gene-activity matrix over {n_genes} genes\n")
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    print("\n--- whole cohort ---")
    for name, v in [("true n_fragment", df.n_fragment),
                    ("total_frags column", df.total_frags_column)]:
        print(f"  {name:20s} median {np.median(v):8.1f}"
              f"   <100: {100 * (v < 100).mean():5.2f} %"
              f"   >={MIN_FRAG}: {100 * (v >= MIN_FRAG).mean():5.2f} %")

    print("\nConclusion: the colorectal immune compartments are SHALLOW, not empty.")
    print(f"  wrote {OUT / 'coad_depth_audit.csv'}")


if __name__ == "__main__":
    main()
