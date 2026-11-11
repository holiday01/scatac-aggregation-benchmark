# Changelog

## 1.0.0 — 2026-11-11 — first archived release

Release archived on Zenodo for the manuscript submission; no existing result table changed.

- `analyses/normal_tissues/08_brain_share_by_type.py` → `brain3k_share_by_type.csv`: share of
  per-gene maxima held by every cell type of the brain multiome under each strategy (the benchmark
  table stores only the largest share; the manuscript also quotes the excitatory-neuron share
  where it is not the largest). Checked against `brain3k_benchmark.csv` when it runs.
- `CITATION.cff` and `.zenodo.json` (citation and archive metadata).
- README: citation section; permutation counts; row-sum/coverage range over all three
  cancer-associated count matrices (0.82–1.16).

## 2026-11-11 — independent markers, RNA-side analyses, real bulk, figures

Added in response to reviewer comments; no existing result table changed.

- `analyses/independent_markers/`: marker concordance with markers from CellTypist models trained
  on independent atlases, label-assignment genes excluded.
- `analyses/rna_sensitivity/`: RNA reference rebuilt without per-cell normalization (raw, pooled,
  log mean), and the same aggregation comparison on RNA from the same nuclei.
- `analyses/real_bulk/`: deconvolution of five PBMC bulk ATAC-seq samples with flow-cytometry
  proportions (EPIC-ATAC, Zenodo 10.5281/zenodo.12571947), plus a post-hoc sensitivity analysis.
- `analyses/deconvolution/rerun_2026-11-11/`: the deconvolution benchmark re-run so that NSCLC
  uses the current matrix; HCC and COAD unchanged. `deconvolution_benchmark.py` now resolves
  `tools/` relative to the repository root.
- `tools/make_figures_v2.py`: the four display figures from the result tables (shared style,
  7.5 in or 170 mm wide).
- GSE201336 is described as one familial-adenomatous-polyposis donor with one carcinoma, three
  polyps and two unaffected mucosae, as annotated in GEO. `celltypist` added to requirements.

## 2026-09-21 — non-tumor tissues and multi-donor RNA references

Added `analyses/normal_tissues/`: the six strategies, imported unchanged, applied
to 26 non-tumor tissues of the human chromatin accessibility atlas (GSE184462,
72 samples, 487,097 cells), to the 10x `human_brain_3k` multiome with RNA-only
labels, and to the ten-donor bone-marrow multiome of GSE194122 using the
gene-activity matrix deposited by its authors. Result tables are included.

Repository housekeeping in the same update:

- `normalisers_and_ci/` moved back under `analyses/`, with
  `benchmark_all_methods.csv` and `benchmark_extra_normalisers.csv` replaced by
  the re-run that uses the `(1 + k) / (B + 1)` permutation estimator. No other
  existing table changed.
- The figure generators (`tools/make_revision_figures.py`,
  `analyses/depth_spread/depth_spread_experiment.py`,
  `analyses/multiome/04_figure.py`) were brought up to date, and the inputs of
  the annotation panel were added (`analyses/multiome/06_annotation_repeats.py`,
  `multiome_annotation_repeats.csv`).
- `tools/report_numbers.py` now also reports the non-tumor, brain and
  bone-marrow quantities.
- Removed files that no current analysis uses: the duplicate top-level
  `practice_survey/`, superseded figure and diagnostic scripts in `tools/`
  (`make_figures.py`, `make_benchmark_figure.py`, `audit_*.py`,
  `benchmark_wilcoxon_markers.py`), `scripts/91_corrected_prior_headtohead.py`,
  `results/validation/`, and the seven output tables of those scripts in
  `results/` (`benchmark_wilcoxon.csv`, `coad_depth_audit.csv`,
  `depth_generality.csv`, `marker_concordance.csv`,
  `marker_concordance_RAW_superseded.csv`, `zscale_diagnostic.csv`,
  `zscale_gene_argmax.csv`). They remain available in the repository history.

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
