"""Score every per-type profile produced by the real tools (Seurat profile CSVs, Signac profile CSVs) with the
paper's read-outs (top argmax share over expressed genes; held-out 22-marker concordance + permutation null),
and summarise the FindAllMarkers tables (markers per type vs true depth; canonical markers in own list)."""
import re, glob, numpy as np, pandas as pd, hdf5plugin, anndata as ad, scipy.sparse as sp, scipy.io
from scipy.stats import spearmanr
from pathlib import Path
from score_profiles import OUTDIR, SCRATCH, BASE, PANEL, relabel, score, true_depth
RES = OUTDIR / "results"
a = ad.read_h5ad(BASE / "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad")
lab = relabel(a.obs["cell_type"].values); types = sorted(pd.unique(lab)); genes = a.var_names.values.astype(str)
X = sp.csr_matrix(a.X); expressed = np.asarray(X.sum(0)).ravel() > 0
dep = true_depth("HCC", a.obs_names.values); med_depth = {t: float(np.nanmedian(dep[lab == t])) for t in types}
seurat_names = np.array([g.replace("_", "-") for g in genes])          # Seurat replaces '_' by '-' in feature names

WHAT = {  # exact call -> what it computes
    "seurat_profile_AverageExpression_default": "NormalizeData() [LogNormalize, 1e4] -> AverageExpression(group.by='cell_type') [layer='data' default: expm1 then mean = CP10k mean]",
    "seurat_profile_AverageExpression_counts": "AverageExpression(group.by='cell_type', layer='counts') [mean of raw counts]",
    "seurat_profile_AverageExpression_logdata_layer": "AverageExpression(layer=<copy of data layer under another name>) [mean of log1p(CP10k); no expm1 because layer name != 'data']",
    "seurat_profile_AverageExpression_returnseurat_data": "AverageExpression(return.seurat=TRUE) data layer = log1p(mean(expm1(data))) [monotone transform of the CP10k mean]",
    "seurat_profile_manual_mean_of_data_layer": "rowMeans(LayerData(obj, layer='data')) per type [what DotPlot/DoHeatmap-style summaries use: mean of log values]",
    "seurat_profile_AggregateExpression_sum": "AggregateExpression(group.by='cell_type') [sum of counts per type]",
    "seurat_profile_AggregateExpression_returnseurat_data": "AggregateExpression(return.seurat=TRUE) data layer = LogNormalize(sums, 1e4) [pseudobulk CPM, log]",
    "signac_profile_AverageExpression_default_scalefactor_median": "Signac vignette: GeneActivity() -> NormalizeData(LogNormalize, scale.factor=median(nCount_RNA)) -> AverageExpression() [expm1 then mean]",
    "signac_profile_mean_of_data_layer_scalefactor_median": "Signac vignette normalisation, then mean of the log data layer per type",
    "signac_profile_mean_of_counts": "GeneActivity() raw counts, mean per type",
    "signac_profile_AverageExpression_default_scalefactor_1e4": "GeneActivity() -> NormalizeData(scale.factor=1e4) -> AverageExpression() [expm1 then mean]",
    "signac_profile_mean_of_data_layer_scalefactor_1e4": "GeneActivity() -> NormalizeData(scale.factor=1e4), mean of the log data layer per type",
}
rows, mrows = [], []
for f in sorted(glob.glob(str(RES / "seurat_profile_*.csv"))) + sorted(glob.glob(str(RES / "signac_profile_*.csv"))):
    key = Path(f).stem; df = pd.read_csv(f)
    if key.startswith("seurat_profile"):
        assert len(df) == len(genes) and (df.gene.values == seurat_names).all(), f"{key}: gene order mismatch"
        P = df[types].values.T; g_use, expr = genes, expressed
    else:  # Signac GeneActivity has its own gene set
        gsig = df.gene.values.astype(str)
        M = scipy.io.mmread(str(SCRATCH / "signac_ga/matrix.mtx")).tocsc()
        feats = np.array(open(SCRATCH / "signac_ga/features.txt").read().split("\n")[:M.shape[0]])
        # CreateAssayObject() applied make.unique() to duplicated gene names; rows keep the FeatureMatrix order
        assert len(feats) == len(gsig) and (feats == gsig).mean() > 0.99, f"{key}: gene order mismatch"
        expr = np.asarray(M.sum(1)).ravel() > 0; P = df[types].values.T; g_use = gsig
    r, mt = score(P, types, g_use, expr)
    tool = "Signac" if key.startswith("signac") else "Seurat"
    rows.append(dict(tool=tool, path=key.replace("seurat_profile_", "").replace("signac_profile_", ""), computes=WHAT.get(key, ""), **r))
    mt.insert(0, "path", key); mrows.append(mt)
    print(f"{key:60} top {r['top_type']:>20} {r['top_share']:5.1%}  markers {r['markers_correct']}/{r['markers_total']} perm p={r['perm_p']:.4f}  {r['per_lineage']}")
pd.DataFrame(rows).to_csv(RES / "seurat_signac_profile_results.csv", index=False)
pd.concat(mrows).to_csv(RES / "seurat_signac_profile_markers.csv", index=False)

# ── FindAllMarkers summaries ──────────────────────────────────────────────────
panel = PANEL["HCC"]; out = []
for f in sorted(f for f in glob.glob(str(RES / "seurat_findallmarkers_*.csv")) if "summary" not in f):
    tag = Path(f).stem.replace("seurat_findallmarkers_", ""); fm = pd.read_csv(f)
    fm["gene_orig"] = fm.gene.map(dict(zip(seurat_names, genes)))
    fm["cluster"] = fm.cluster.astype(str)
    n_used = {t: int((lab == t).sum()) for t in types}
    if tag.startswith("reduced"):
        sub = pd.read_csv(SCRATCH / "hcc_sub/cells.tsv", sep="\t"); n_used = sub.cell_type.value_counts().to_dict()
    for t in types:
        s = fm[fm.cluster == t].reset_index(drop=True)           # already ordered by p_val within cluster
        own = [g for g, tg in panel.items() if tg == t]; foreign = [g for g, tg in panel.items() if tg != t]
        pos = {g: int(i) + 1 for i, g in enumerate(s.gene_orig)}
        out.append(dict(layer=tag, cell_type=t, median_true_depth=med_depth[t], n_cells=int(n_used[t]),
                        n_returned=len(s), n_padj_lt_0_05=int((s.p_val_adj < 0.05).sum()),
                        n_padj_lt_0_05_log2fc_gt_0_5=int(((s.p_val_adj < 0.05) & (s.avg_log2FC > 0.5)).sum()),
                        own_markers=len(own), own_in_list=sum(g in pos for g in own),
                        own_in_list_padj=sum((g in pos) and (s.p_val_adj[pos[g]-1] < 0.05) for g in own),
                        own_in_top50=sum(pos.get(g, 10**9) <= 50 for g in own),
                        own_ranks="; ".join(f"{g}:{pos.get(g, 'absent')}" for g in own),
                        foreign_in_list=sum(g in pos for g in foreign), foreign_in_top50=sum(pos.get(g, 10**9) <= 50 for g in foreign),
                        top10=",".join(s.gene_orig.head(10))))
if out:
    F = pd.DataFrame(out)
    for tag, s in F.groupby("layer"):
        for col in ["n_returned", "n_padj_lt_0_05", "n_padj_lt_0_05_log2fc_gt_0_5"]:
            rho, p = spearmanr(s.median_true_depth, s[col]); F.loc[s.index, f"rho_depth_{col}"] = rho
            print(f"[{tag}] Spearman(median true depth, {col}) = {rho:+.2f} (p={p:.3f})")
        print(s[["cell_type", "median_true_depth", "n_returned", "n_padj_lt_0_05", "n_padj_lt_0_05_log2fc_gt_0_5", "own_in_list", "own_markers", "own_in_top50", "foreign_in_top50", "own_ranks"]].to_string(index=False))
    F.to_csv(RES / "seurat_findallmarkers_summary.csv", index=False)
print("done")
