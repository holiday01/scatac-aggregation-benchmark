"""
Are the HCC cell labels right? The labels of GSE227265 were assigned by us from
Leiden clusters of the deposited TF-motif matrix (EBF1 -> "B_cell", IRF8 -> "DC").
The marker test flagged that every canonical B marker peaks in the "DC" cluster.
This script asks the question with marker panels DISJOINT from the paper's
22-marker test panel, on per-cell-normalised (CPM) cell-type means:

  * endothelial / stromal markers  -> which cluster?
  * B-cell markers not in the test panel (PAX5, EBF1, BLK, ...) and plasma markers -> which cluster?
  * bona fide dendritic-cell markers (CLEC9A, XCR1, CD1C, ...) -> which cluster?

Run:  python3 tools/check_hcc_labels.py
Out:  results/hcc_label_check.csv
"""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, h5py, hdf5plugin  # noqa
import os
from pathlib import Path
BASE = Path(os.environ.get("SCRNA_ATAC_BASE", ".")); OUT = Path(__file__).resolve().parent.parent / "results"
PANELS = {
    "endothelial": ["PECAM1", "VWF", "CDH5", "KDR", "CLDN5", "FLT1", "ENG", "TEK", "ERG", "EMCN", "PLVAP", "ESAM"],
    "stromal":     ["COL1A2", "COL3A1", "DCN", "PDGFRB", "RGS5", "COL6A3", "FBLN1", "THY1"],
    "B_heldout":   ["PAX5", "EBF1", "BLK", "CD22", "FCRL1", "TNFRSF13C", "VPREB3", "FCRLA"],
    "plasma":      ["JCHAIN", "MZB1", "IGHG1", "TNFRSF17", "IGKC", "DERL3"],
    "DC_bonafide": ["CLEC9A", "XCR1", "CD1C", "FCER1A", "CLEC10A", "LAMP3", "ITGAX", "ZBTB46", "FLT3", "IRF8"],
    "hepatocyte":  ["ALB", "APOB", "HNF4A", "TTR", "APOA1"],
    "macrophage":  ["CD68", "CSF1R", "MRC1", "MSR1", "CD163", "ITGAM"],
    "T_NK":        ["CD3D", "CD3E", "CD8A", "IL7R", "NKG7", "GNLY"],
}
rows = []
for model in ["proximal", "enhancer"]:
    a = ad.read_h5ad(BASE / f"data/processed/integrated/atac_hcc_geneactivity_{model}.h5ad")
    X = a.X.tocsr().astype(np.float32); lab = a.obs.cell_type.astype(str).values; types = sorted(set(lab))
    tot = np.asarray(X.sum(1)).ravel(); Xc = sp.diags(1e4 / np.maximum(tot, 1)) @ X
    P = pd.DataFrame(np.vstack([np.asarray(Xc[lab == t].mean(0)).ravel() for t in types]), index=types, columns=a.var_names)
    print(f"== {model} ==")
    for k, gl in PANELS.items():
        gl = [g for g in gl if g in P.columns]; am = P[gl].idxmax(0)
        vc = am.value_counts()
        print(f"  {k:12} {len(gl):2d} markers -> " + ", ".join(f"{t} {n}" for t, n in vc.items()))
        for g, t in am.items():
            rows.append(dict(model=model, panel=k, gene=g, argmax_label=t))
pd.DataFrame(rows).to_csv(OUT / "hcc_label_check.csv", index=False)
print(f"wrote {OUT/'hcc_label_check.csv'}")
