"""
04: Label the 10x healthy-brain multiome cells (human_brain_3k) from the RNA modality ONLY.
    Same procedure as multiome/02_label_rna.py: scanpy normalize_total -> log1p -> HVG -> PCA ->
    neighbours -> leiden; clusters are named by the mean scaled expression of a canonical brain
    labelling panel. The labelling genes are recorded so the benchmark can exclude them from
    the held-out marker panel.
Out: brain3k_rna_labels.csv (barcode, leiden, cell_type, rna_umi), brain3k_rna_cluster_scores.csv,
     brain3k_labelling_genes.txt
"""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc

H5  = Path("/mnt/10t/scrna_atac/data/raw/multiome_brain3k/human_brain_3k_filtered_feature_bc_matrix.h5")
OUT = Path(__file__).resolve().parent
MIN_CELLS = 100

# labelling panel: gene -> broad type (these genes are EXCLUDED from the held-out test panel)
LABEL = {
    "MOBP": "Oligodendrocyte", "MBP": "Oligodendrocyte", "PLP1": "Oligodendrocyte", "MOG": "Oligodendrocyte",
    "PDGFRA": "OPC", "VCAN": "OPC", "CSPG4": "OPC",
    "AQP4": "Astrocyte", "GFAP": "Astrocyte", "SLC1A2": "Astrocyte", "SLC1A3": "Astrocyte",
    "CSF1R": "Microglia", "P2RY12": "Microglia", "CX3CR1": "Microglia",
    "SLC17A7": "Excitatory_neuron", "CUX2": "Excitatory_neuron", "RORB": "Excitatory_neuron", "SATB2": "Excitatory_neuron",
    "GAD1": "Inhibitory_neuron", "GAD2": "Inhibitory_neuron", "SLC32A1": "Inhibitory_neuron",
    # PDGFRB is not used for the vascular score: in these nuclei it is higher in astrocytes (0.7-0.8)
    # than in the vascular cluster (0.4), which made the doublet rule discard both astrocyte clusters.
    "CLDN5": "Vascular", "FLT1": "Vascular", "VWF": "Vascular",
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
# lineage groups for the doublet rule: a cluster whose best score sits in one group but which
# also scores > DOUBLET_Z in another group is a cross-lineage mixture (doublets) and is dropped.
GROUP = {"Oligodendrocyte": "OL", "OPC": "OL", "Astrocyte": "AST", "Microglia": "MIC",
         "Excitatory_neuron": "NEU", "Inhibitory_neuron": "NEU", "Vascular": "VAS"}
DOUBLET_Z = 0.4
call = {}
for c in S.index:
    s = S.loc[c]
    best = s.idxmax()
    other = max(s[t] for t in types if GROUP[t] != GROUP[best])
    if other > DOUBLET_Z: best = "mixed_doublet"
    call[c] = best
a.obs["cell_type"] = a.obs.leiden.map(call).astype(str)
S["call"] = pd.Series(call); S["n_cells"] = a.obs.leiden.value_counts()
S.to_csv(OUT / "brain3k_rna_cluster_scores.csv")
print(S.round(2).to_string())
vc = a.obs.cell_type.value_counts()
print(vc.to_string())
keep = vc[(vc >= MIN_CELLS) & (vc.index != "mixed_doublet")].index
a.obs["cell_type"] = np.where(a.obs.cell_type.isin(keep), a.obs.cell_type, "dropped")
a.obs[["leiden", "cell_type", "rna_umi"]].to_csv(OUT / "brain3k_rna_labels.csv")
Path(OUT / "brain3k_labelling_genes.txt").write_text("\n".join(genes) + "\n")
print("kept types:", list(keep))
print(a.obs.cell_type.value_counts().to_string())
