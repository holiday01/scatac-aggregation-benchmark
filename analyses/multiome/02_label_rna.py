"""
02: Label PBMC multiome cells from the RNA modality ONLY.
    scanpy normalize_total -> log1p -> HVG -> PCA -> neighbours -> leiden; clusters are named
    by the mean scaled expression of a canonical PBMC labelling panel. Those labelling genes
    are recorded so the benchmark can exclude them from the held-out marker panel.
Out: rna_labels.csv (barcode, leiden, cell_type, rna_umi), rna_cluster_scores.csv, labelling_genes.txt
"""
import warnings; warnings.filterwarnings("ignore")
import os
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc

BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
H5  = BASE / "data/raw/multiome_pbmc10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
OUT = Path(__file__).resolve().parent
MIN_CELLS = 100

# labelling panel: gene -> broad type (these genes are EXCLUDED from the held-out test panel)
LABEL = {
    "CD3D": "T_CD4", "CD3E": "T_CD4", "CD4": "T_CD4", "IL7R": "T_CD4", "CCR7": "T_CD4",
    "CD8A": "T_CD8", "CD8B": "T_CD8",
    "NKG7": "NK", "GNLY": "NK", "KLRD1": "NK", "NCAM1": "NK",
    "CD14": "Mono_CD14", "LYZ": "Mono_CD14", "S100A8": "Mono_CD14", "S100A9": "Mono_CD14",
    "FCGR3A": "Mono_CD16", "MS4A7": "Mono_CD16", "LST1": "Mono_CD16",
    "MS4A1": "B", "CD79A": "B", "CD79B": "B", "CD19": "B",
    "FCER1A": "DC", "CD1C": "DC", "CLEC10A": "DC",
    "LILRA4": "pDC", "IL3RA": "pDC", "CLEC4C": "pDC",
    "MZB1": "Plasma", "JCHAIN": "Plasma", "XBP1": "Plasma",
}

a = sc.read_10x_h5(H5, gex_only=True); a.var_names_make_unique()
a.obs["rna_umi"] = np.asarray(a.X.sum(1)).ravel()
a.layers["counts"] = a.X.copy()
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat")
sc.pp.pca(a, n_comps=30, mask_var="highly_variable")
sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, random_state=0)
sc.tl.leiden(a, resolution=1.0, random_state=0, flavor="igraph", n_iterations=2, directed=False)
print("leiden clusters:", a.obs.leiden.nunique())

# per-cluster score for each type: mean over that type's panel genes of the z-scored cluster mean
genes = [g for g in LABEL if g in a.var_names]
M = pd.DataFrame(np.vstack([np.asarray(a[a.obs.leiden == c, genes].X.mean(0)).ravel() for c in a.obs.leiden.cat.categories]),
                 index=a.obs.leiden.cat.categories, columns=genes)
Z = (M - M.mean(0)) / (M.std(0) + 1e-9)
types = sorted(set(LABEL.values()))
S = pd.DataFrame({t: Z[[g for g in genes if LABEL[g] == t]].mean(1) for t in types})
# T_CD8 vs T_CD4 vs NK need care: require the T score to be positive for T calls
# lineage groups for the doublet rule: a cluster whose best score sits in one group but which
# also scores > DOUBLET_Z in another group is a cross-lineage mixture (doublets) and is dropped.
GROUP = {"T_CD4": "TNK", "T_CD8": "TNK", "NK": "TNK", "B": "B", "Plasma": "B",
         "Mono_CD14": "MYE", "Mono_CD16": "MYE", "DC": "MYE", "pDC": "pDC"}
DOUBLET_Z = 0.4
call = {}
for c in S.index:
    s = S.loc[c]
    best = s.idxmax()
    if best == "T_CD8" and s["NK"] > s["T_CD8"] and Z.loc[c, "CD3E"] < 0: best = "NK"
    # pDC genuinely express MZB1/JCHAIN, so the Plasma score is not evidence of a pDC+plasma doublet
    other = max(s[t] for t in types if GROUP[t] != GROUP[best] and not (best == "pDC" and t == "Plasma"))
    if other > DOUBLET_Z: best = "mixed_doublet"
    call[c] = best
a.obs["cell_type"] = a.obs.leiden.map(call).astype(str)
S["call"] = pd.Series(call); S["n_cells"] = a.obs.leiden.value_counts()
S.to_csv(OUT / "rna_cluster_scores.csv")
vc = a.obs.cell_type.value_counts()
print(vc.to_string())
keep = vc[(vc >= MIN_CELLS) & (vc.index != "mixed_doublet")].index
a.obs["cell_type"] = np.where(a.obs.cell_type.isin(keep), a.obs.cell_type, "dropped")
a.obs[["leiden", "cell_type", "rna_umi"]].to_csv(OUT / "rna_labels.csv")
Path(OUT / "labelling_genes.txt").write_text("\n".join(genes) + "\n")
print("kept types:", list(keep))
print(a.obs.cell_type.value_counts().to_string())
# sanity: key marker means per called type (log-normalised)
chk = ["CD3E","CD4","CD8A","NKG7","GNLY","CD14","LYZ","FCGR3A","MS4A1","FCER1A","LILRA4","MZB1"]
print(pd.DataFrame(np.vstack([np.asarray(a[a.obs.cell_type == t, chk].X.mean(0)).ravel() for t in a.obs.cell_type.unique()]),
                   index=a.obs.cell_type.unique(), columns=chk).round(2).to_string())
