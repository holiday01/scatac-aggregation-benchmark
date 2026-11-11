# Re-run of 2026-11-11

`../deconvolution_benchmark.py`, unchanged, re-run so that NSCLC (LUAD in the tables) uses the
current matrix (58,112 cells, eight samples); the tables one level up were produced on an earlier
NSCLC build (53,308 cells). HCC and COAD reproduce the earlier tables exactly (largest absolute
difference about 1e-16). `estimates.csv` (18 MB, every mixture × type estimate) is not stored.

The manuscript's deconvolution table uses `metrics_summary.csv` from this folder
(`geneset == "common"`; `mixture_variant` rawsum = raw-count mixture, cellnorm =
equal-contribution mixture; `truth` cellfrac / fragfrac).
