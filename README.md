# Aggregation scale and library-size-associated variation in scATAC-seq gene-activity profiles

Analysis code and derived result tables for a comparison of the scales on which
single cells are aggregated into per-cell-type gene-activity profiles.

This repository holds **code and result tables only**. It contains no manuscript
text and no figure image files; the figures are regenerated from the scripts and
the tables below.

## What is here

```
tools/      matrix construction, benchmark, revision analyses, figure generation,
            and report_numbers.py, which recomputes every reported quantity
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
| GSE201336 | colorectal | six colon specimens from one donor with familial adenomatous polyposis (one carcinoma, three polyps, two unaffected mucosae per GEO), unsorted tissue nuclei |
| GSE270148 | non-small-cell lung cancer | CD45-positive bead-enriched nuclei; eight of nine tumour samples used |
| 10x `pbmc_granulocyte_sorted_10k` | PBMC multiome | RNA and ATAC in the same nuclei, one donor |
| GSE184462 (adult samples) | 26 non-tumor tissues | human chromatin accessibility atlas, 72 fragment files, depositors' cell-type labels |
| 10x `human_brain_3k` | brain multiome | RNA and ATAC in the same nuclei, one donor, labels from RNA only |
| GSE194122 (multiome) | bone marrow | ten donors, 13 batches; gene-activity matrix as deposited by its authors |
| Zenodo [10.5281/zenodo.12571947](https://doi.org/10.5281/zenodo.12571947) | PBMC bulk ATAC-seq | five samples with flow-cytometry proportions (EPIC-ATAC study), used for real-bulk deconvolution |

Marker genes for `analyses/independent_markers/` come from CellTypist models (celltypist 1.7.1)
trained on independent atlases; see that folder's REPORT.md.

No new data were generated.

## Citation

Each release of this repository is archived on Zenodo; please cite the archived
version (metadata in `CITATION.cff` and `.zenodo.json`).

## Reproducing

```bash
pip install -r requirements.txt            # or see environment.txt
python3 scripts/94_luad_prior_and_reference.py   # NSCLC gene-activity matrices
python3 tools/rebuild_nsclc_and_rerun.py         # equivalent standalone rebuild
python3 tools/benchmark_aggregation.py           # -> results/benchmark_*.csv
python3 tools/revision_analyses.py               # integrality, leave-one-sample-out, label sensitivity
python3 tools/check_hcc_labels.py                # evidence for the two renamed HCC clusters
python3 tools/make_revision_figures.py           # regenerates the figures locally (earlier version)
python3 tools/make_figures_v2.py                 # current display figures (FIG_PROFILE=bmc for 170 mm)
python3 tools/report_numbers.py                  # recomputes and prints every reported quantity from the tables
# non-tumor tissues, brain and bone-marrow multiome: run order in analyses/normal_tissues/README.md
python3 analyses/independent_markers/01_independent_markers.py      # markers from independent atlases
python3 analyses/rna_sensitivity/01_rna_reference_sensitivity.py   # RNA reference variants; RNA vs ATAC
python3 analyses/real_bulk/01_real_bulk_deconvolution.py           # real bulk with flow cytometry
```

Scripts resolve paths from `$SCRNA_ATAC_BASE` (default `/mnt/10t/scrna_atac`).
Random seeds are fixed in the code: 42 for the benchmark and the permutation
nulls, 0–2 for the controlled depth series.

## Notes on the analysis

- The permutation null uses the one-sided `(1 + k) / (B + 1)` estimator, so the
  smallest reportable p with B = 10,000 is 9.999 × 10⁻⁵ rather than zero (B = 2,000
  in the controlled depth series).
- Binomial thinning is applied only to the three integer count matrices. The
  distance-weighted matrices hold fractional scores (52–56 % of their non-zero
  entries are below 0.5), on which count-based operations are not defined.
- Per-cell coverage is the fragment count recorded at import (liver, colorectal)
  or the total peak count (lung), never the gene-window row sum; across the cell
  types of the three cancer-associated count matrices the two differ by
  0.82–1.16-fold because gene windows overlap.
- Uncertainty is assessed by leaving out one sample at a time, not by resampling
  cells; in the bone-marrow multiome one donor is left out at a time.
- The non-tumor tissues are aggregated one tissue at a time with the depositors'
  cell-type labels; cell types with fewer than 100 cells in a tissue are dropped.
- The bone-marrow gene-activity matrix is the one deposited by its authors. It is
  stored log-normalized, and the integer counts are recovered exactly from it
  (see `analyses/normal_tissues/README.md`).
