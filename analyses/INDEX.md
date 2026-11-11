# Supporting analyses

Each folder holds its scripts, result tables and a REPORT.md (or README.md). Aggregation
functions are imported from `tools/benchmark_aggregation.py` throughout.

| folder | question |
|---|---|
| `depth_spread/` | how each strategy behaves when a coverage spread of 1- to 40-fold between cell types is imposed |
| `normalisers_and_ci/` | ten further transforms; cell-level bootstrap intervals |
| `multiome/` | agreement with RNA from the same nuclei and annotation accuracy (10x PBMC multiome) |
| `normal_tissues/` | 26 non-tumor tissues (GSE184462), brain multiome, ten-donor bone-marrow multiome |
| `deconvolution/` | NNLS deconvolution with references built by each strategy |
| `real_tools/` | what Seurat, Signac, scanpy and SnapATAC2 return on their documented paths |
| `practice_survey/` | how published studies summarize gene activity per cell type |
| `independent_markers/` | marker concordance with markers from independent atlases |
| `rna_sensitivity/` | RNA reference variants; the same comparison on RNA from the same nuclei |
| `real_bulk/` | deconvolution of real bulk ATAC-seq with flow-cytometry proportions |
