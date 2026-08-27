"""Reduced problem for FindAllMarkers (the full 12,029 x 59,385 default-Wilcoxon run without presto exceeds the
time budget): at most 500 cells per type (seed 1), expressed protein-coding genes (GENCODE v38 gene_type)."""
import numpy as np, pandas as pd, hdf5plugin, anndata as ad, scipy.sparse as sp, scipy.io, re
from score_profiles import SCRATCH, BASE, relabel, true_depth
rng = np.random.default_rng(1)
a = ad.read_h5ad(BASE / "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad")
lab = relabel(a.obs["cell_type"].values)
keep = np.concatenate([rng.choice(np.where(lab == t)[0], size=min(500, (lab == t).sum()), replace=False) for t in sorted(set(lab))])
keep.sort()
pc = set()
for line in open(SCRATCH / "hcc_frag/gencode.v38.genes.gtf"):
    if 'gene_type "protein_coding"' in line:
        pc.add(re.search(r'gene_name "([^"]+)"', line).group(1))
X = sp.csr_matrix(a.X)[keep]
gmask = np.array([g in pc for g in a.var_names]) & (np.asarray(X.sum(0)).ravel() > 0)
out = SCRATCH / "hcc_sub"; out.mkdir(exist_ok=True)
pd.DataFrame({"barcode": a.obs_names[keep], "cell_type": lab[keep], "true_depth": true_depth("HCC", a.obs_names[keep].values)}).to_csv(out / "cells.tsv", sep="\t", index=False)
pd.DataFrame({"gene": a.var_names[gmask]}).to_csv(out / "features.tsv", sep="\t", index=False, header=False)
scipy.io.mmwrite(str(out / "matrix.mtx"), X[:, gmask].T.tocsc().astype(np.int32), field="integer")
print(f"subsample: {len(keep)} cells, {int(gmask.sum())} expressed protein-coding genes; per type:", pd.Series(lab[keep]).value_counts().to_dict())
