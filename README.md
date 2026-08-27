# Aggregation scale and library-size-associated variation in scATAC-seq gene-activity profiles

Analysis code and derived result tables for a comparison of the scales on which
single cells are aggregated into per-cell-type gene-activity profiles.

This repository holds **code and result tables only**. It contains no manuscript
text and no figure image files; the figures are regenerated from the scripts and
the tables below.

## What is here

```
tools/      matrix construction, benchmark, revision analyses, figure generation,
            claim and reference verification
scripts/    gene-activity matrix construction from the deposited GEO files
results/    derived tables that every reported number is computed from,
            including results/revision/
analyses/   supporting analyses, each with its scripts, tables and a REPORT.md
data/       accessions and checksums (no data files are stored here)
```

## Datasets

| accession | tissue | note |
|---|---|---|
| GSE227265 | primary liver cancer | 13 hepatocellular tumours of 16 patients used, one per patient |
| GSE201336 | colorectal | six colon specimens from one donor, unsorted tissue nuclei |
| GSE270148 | non-small-cell lung cancer | CD45-positive bead-enriched nuclei; eight of nine tumour samples used |
| 10x `pbmc_granulocyte_sorted_10k` | PBMC multiome | RNA and ATAC in the same nuclei, one donor |

No new data were generated.

## Reproducing

```bash
pip install -r requirements.txt            # or see environment.txt
python3 scripts/94_luad_prior_and_reference.py   # NSCLC gene-activity matrices
python3 tools/rebuild_nsclc_and_rerun.py         # equivalent standalone rebuild
python3 tools/benchmark_aggregation.py           # -> results/benchmark_*.csv
python3 tools/revision_analyses.py               # integrality, leave-one-sample-out, label sensitivity
python3 tools/check_hcc_labels.py                # evidence for the two renamed HCC clusters
python3 tools/make_revision_figures.py           # regenerates the figures locally
python3 tools/report_numbers.py                  # recomputes and prints every reported quantity from the tables
```

Scripts resolve paths from `$SCRNA_ATAC_BASE` (default `/mnt/10t/scrna_atac`).
Random seeds are fixed in the code: 42 for the benchmark and the permutation
nulls, 0–2 for the controlled depth series.

## Notes on the analysis

- The permutation null uses the `(1 + k) / (B + 1)` estimator, so the smallest
  reportable p with B = 10,000 is 9.999 × 10⁻⁵ rather than zero.
- Binomial thinning is applied only to the three integer count matrices. The
  distance-weighted matrices hold fractional scores (52–56 % of their non-zero
  entries are below 0.5), on which count-based operations are not defined.
- Per-cell coverage is the fragment count recorded at import (liver, colorectal)
  or the total peak count (lung), never the gene-window row sum; the two differ
  by 0.93–1.16-fold across cell types because gene windows overlap.
- Uncertainty is assessed by leaving out one sample at a time, not by resampling
  cells.
