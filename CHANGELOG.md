# Changelog

## 1.0.1 — 2026-11-11

Documentation, paths and run logs tidied; analysis code and result tables unchanged.

- Scripts import `tools/` from this repository; `tools/revision_analyses.py` writes to
  `results/revision/`; the real-tool comparison writes its intermediate files to
  `$REAL_TOOLS_WORK` (default `$SCRNA_ATAC_BASE/work/real_tools`).
- Script and report headers describe what each script does in this benchmark;
  `analyses/INDEX.md` lists the supporting analyses.
- `CITATION.cff` and README: Zenodo DOIs.

## 1.0.0 — 2026-11-11

First archived release: code and result tables for the benchmark of six aggregation
strategies on three cancer-associated datasets (GSE227265, GSE201336, GSE270148), 26
non-tumor tissues of GSE184462, and the matched RNA-ATAC datasets (10x PBMC and brain
multiome, GSE194122 bone marrow), with the supporting analyses in `analyses/`
(controlled coverage series, further transforms, deconvolution, real-tool paths,
independent marker panels, RNA-side analyses, real bulk samples).
