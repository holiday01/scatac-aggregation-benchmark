"""Export the HCC proximal gene-activity counts to Matrix Market (genes x cells) + cell/feature tables
so that Seurat can read it with ReadMtx. Labels: original and relabelled (B_cell->Endothelial_stromal, DC->B_cell)."""
import time, numpy as np, pandas as pd, hdf5plugin, anndata as ad, scipy.sparse as sp, scipy.io
from score_profiles import SCRATCH, BASE, true_depth, relabel
t0 = time.time()
a = ad.read_h5ad(BASE / "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad")
X = sp.csr_matrix(a.X)
print(a.shape, "nnz", X.nnz, "dtype", X.dtype, "min", X.min(), "max", X.max(), flush=True)
out = SCRATCH / "hcc_mtx"; out.mkdir(exist_ok=True, parents=True)
cells = pd.DataFrame({"barcode": a.obs_names, "cell_type_original": a.obs["cell_type"].astype(str).values,
                      "cell_type": relabel(a.obs["cell_type"].values),
                      "true_depth": true_depth("HCC", a.obs_names.values)})
cells.to_csv(out / "cells.tsv", sep="\t", index=False)
pd.DataFrame({"gene": a.var_names}).to_csv(out / "features.tsv", sep="\t", index=False, header=False)
print(cells.cell_type.value_counts(), flush=True)
scipy.io.mmwrite(str(out / "matrix.mtx"), X.T.tocsc().astype(np.int32), field="integer")
print("wrote", out, f"{time.time()-t0:.0f}s", flush=True)
