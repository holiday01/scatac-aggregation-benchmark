# Changelog

## 2026-08-27 — NSCLC sample correction (affects every NSCLC-dependent result)

The local copy of `GSM8069383` (patient Lu934) used in the first analysis was a
duplicate of `GSM8069382` (Lu931): same MD5, same file size and the same
3,825 × 145,746 barcode matrix. All nine deposited peak matrices of GSE270148
were re-downloaded and compared with GEO; **only that one file was affected**.

`GSM8069383` is a valid sample (4,818 cells, median 10,316 counts, 4,804
barcodes above the 1000-count threshold) and is now included. `GSM8069377`
(median 304 counts, 371 barcodes above threshold) remains the only exclusion,
re-verified on the correct file. The NSCLC dataset is therefore **eight samples
and 58,112 barcodes**, not seven and 53,308.

Consequences:

- `scripts/93_luad_scatac_evaluation.py` no longer claims a GEO duplicate.
- `scripts/94_luad_prior_and_reference.py` excludes only `GSM8069377`.
- The gene-activity matrices, the benchmark, the leave-one-sample-out ranges, the
  controlled depth series, the additional normalisers and the marker analyses
  were all recomputed; every table in `results/` is from the corrected data.

## 2026-08-27 — other corrections

- `tools/make_revision_figures.py` ended by copying a stale figure over the one
  it had just generated, so re-runs silently restored an old panel. The copy was
  removed.
- The permutation p estimator was changed from `mean(null >= observed)` to
  `(1 + k) / (B + 1)`, so p is never exactly zero. Re-running changed no
  reported value.
- `rank_genes_groups` on unnormalised counts was dropped from the analysis: the
  function documents log-transformed input, so that comparison was outside its
  stated assumptions.
