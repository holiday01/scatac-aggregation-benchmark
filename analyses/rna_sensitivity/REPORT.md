# RNA-side analyses on the matched RNA–ATAC datasets

Datasets: PBMC (10x `pbmc_granulocyte_sorted_10k`), brain (10x `human_brain_3k`), bone marrow
(GSE194122, deposited gene activity). Same cells and labels as `../multiome` and
`../normal_tissues`; aggregation functions imported from `tools/benchmark_aggregation.py`.

## A. Does the conclusion depend on the RNA reference being normalized per cell?

The RNA reference was rebuilt as the CP10K mean (as in the benchmark), raw mean, pooled CPM and
log-normalized mean, and every ATAC strategy was scored against each
(`rna_reference_sensitivity.csv`).

- **Pooled RNA reference** (no per-cell rescaling): normalized ATAC strategies agree with it for
  64–69 % (PBMC), 74–79 % (brain) and 50–56 % (bone marrow) of its RNA-informative genes, against
  63–65 %, 42–43 % and 23 % for the raw ATAC mean.
- **Raw RNA mean** is itself coverage-associated: 23–73 % of maxima in one type (12–26 % for the
  CP10K reference), share of maxima vs UMI coverage ρ = +0.88 to +0.95. On the genes it calls
  informative it agrees best with the raw ATAC mean in PBMC and brain (60–82 % vs 51–60 %), not in
  bone marrow.
- On the genes informative under the CP10K reference, every normalized ATAC strategy agrees
  better than the raw ATAC mean with every RNA reference in every dataset.

## B. The same aggregation question asked of RNA in the same nuclei

`rna_vs_atac_aggregation.csv`:

- Detection: 13–42 % of genes per cell in the count/deposited gene-activity matrices, 4–17 % in RNA
  (gene activity is not sparser than RNA here).
- Log-normalized mean of RNA: share of maxima follows UMI coverage (ρ = +0.70 to +0.83; CP10K mean
  −0.36 to +0.12).
- Spread across types of the per-type mean of log(1+CP10K) over all genes (the residual of the
  log mean): RNA 1.30 / 1.47 / 2.36 (PBMC / brain / bone marrow); gene activity 1.08–1.09 /
  1.14–1.16 / 1.60. The log-scale residual is not specific to chromatin.

Run: `python3 01_rna_reference_sensitivity.py`.
