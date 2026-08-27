"""Merge all real-tool profile read-outs into one table next to the paper's re-implementation numbers
(results/benchmark_aggregation.csv, HCC proximal, corrected labels). Writes results/all_tool_paths.csv and prints markdown."""
import pandas as pd, numpy as np
from pathlib import Path
RES = Path(__file__).resolve().parent / "results"
ROOT = Path(__file__).resolve().parents[2]
paper = pd.read_csv(ROOT / "results/benchmark_aggregation.csv")
paper = paper[(paper.cohort == "HCC") & (paper.model == "proximal")].set_index("method")
EQUIV = {  # tool path -> paper method it should reproduce
    ("Seurat", "AverageExpression_default"): "cpm_mean", ("Seurat", "AverageExpression_returnseurat_data"): "cpm_mean",
    ("Seurat", "AverageExpression_counts"): "raw_mean", ("Seurat", "AverageExpression_logdata_layer"): "lognorm_mean",
    ("Seurat", "manual_mean_of_data_layer"): "lognorm_mean", ("Seurat", "AggregateExpression_sum"): "(sum; no paper equivalent)",
    ("Seurat", "AggregateExpression_returnseurat_data"): "pseudobulk_cpm",
    ("Signac", "mean_of_counts"): "raw_mean", ("Signac", "AverageExpression_default_scalefactor_median"): "cpm_mean",
    ("Signac", "AverageExpression_default_scalefactor_1e4"): "cpm_mean", ("Signac", "mean_of_data_layer_scalefactor_median"): "lognorm_mean",
    ("Signac", "mean_of_data_layer_scalefactor_1e4"): "lognorm_mean",
    ("SnapATAC2", "raw_mean"): "raw_mean", ("SnapATAC2", "tutorial_normalize_total_median_log1p_mean"): "lognorm_mean",
    ("SnapATAC2", "normalize_total_1e4_log1p_mean"): "lognorm_mean", ("SnapATAC2", "cpm_mean_1e4"): "cpm_mean",
    ("SnapATAC2", "cpm_mean_median"): "cpm_mean", ("SnapATAC2", "pseudobulk_cpm"): "pseudobulk_cpm",
    ("scanpy", "sc.get.aggregate_mean"): "lognorm_mean", ("scanpy", "obs_df_groupby_mean"): "lognorm_mean",
    ("scanpy", "paper_lognorm_mean"): "lognorm_mean", ("scanpy", "paper_cpm_mean"): "cpm_mean", ("scanpy", "paper_raw_mean"): "raw_mean",
}
parts = []
for f in ["seurat_signac_profile_results.csv", "snapatac2_results.csv", "scanpy_results.csv"]:
    if (RES / f).exists(): parts.append(pd.read_csv(RES / f))
T = pd.concat(parts, ignore_index=True)
cols = ["tool", "path", "computes", "top_type", "top_share", "markers_correct", "markers_total", "perm_p", "per_lineage", "n_genes_expressed"]
T = T[cols]
T["paper_equivalent"] = [EQUIV.get((a, b), "") for a, b in zip(T.tool, T.path)]
T["paper_top_share"] = [paper.top_share.get(m, np.nan) for m in T.paper_equivalent]
T["paper_markers"] = [f"{int(paper.markers_correct[m])}/{int(paper.markers_total[m])}" if m in paper.index else "" for m in T.paper_equivalent]
T.to_csv(RES / "all_tool_paths.csv", index=False)
print("| tool | path | top type | top share | held-out markers | paper equivalent | paper top share | paper markers |")
print("|---|---|---|---|---|---|---|---|")
for r in T.itertuples():
    print(f"| {r.tool} | {r.path} | {r.top_type} | {r.top_share:.1%} | {r.markers_correct}/{r.markers_total} (perm p={r.perm_p:.4f}) | {r.paper_equivalent} | "
          f"{'' if np.isnan(r.paper_top_share) else f'{r.paper_top_share:.1%}'} | {r.paper_markers} |")
