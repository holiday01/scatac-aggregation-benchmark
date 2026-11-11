# Deconvolution of real bulk ATAC-seq with flow-cytometry proportions

**Data.** Five PBMC bulk ATAC-seq samples from the EPIC-ATAC study (Gabriel et al., eLife 2024),
Zenodo [10.5281/zenodo.12571947](https://doi.org/10.5281/zenodo.12571947) (CC BY 4.0):
`PBMC_counts.txt` (peak × sample counts, hg38) and `PBMC_cell_fractions.txt` (flow cytometry).
Not stored here; download both files unchanged to `$SCRNA_ATAC_BASE/data/raw/epic_atac_pbmc_bulk/`
(checksums in `../../data/checksums.md`).

**Method.** Bulk gene activity = counts of all peaks overlapping a gene's count-model window
(gene body + 2 kb upstream, GENCODE v38; path via `$GENCODE_GTF`). References from the 10x PBMC
multiome count matrix with RNA-derived labels, one per aggregation strategy. Signature genes:
200 per type by specificity on the CP10K-mean reference (shared by all strategies), restricted to
genes with bulk counts. NNLS per sample on CPM over signature genes; estimates summed into the
flow groups (B; monocytes; CD4 T incl. Treg; CD8 T; NK; DC incl. pDC).

**Result** (`real_bulk_metrics.csv`): every strategy recovers the flow proportions poorly
(RMSE 0.15–0.18, Pearson r 0.16–0.33), mainly through over-estimating dendritic cells
(+0.24 to +0.31); the bulk-peak vs single-nucleus difference dominates.

**Post-hoc sensitivity** (`02_real_bulk_noDC_posthoc.py`, decided after seeing the primary
result): DC and pDC removed from the reference and the truth renormalized. RMSE 0.13–0.15 for
every strategy except the z-scaled mean (0.21); the reference cell types differ in median row sum
by only about 1.2-fold. These samples do not rank the linear, log and raw references.

Run: `python3 01_real_bulk_deconvolution.py` and `python3 02_real_bulk_noDC_posthoc.py`.
