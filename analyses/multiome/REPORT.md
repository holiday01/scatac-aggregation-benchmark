# Multiome ground truth for the aggregation benchmark (10x PBMC, RNA + ATAC in the same cell)

Date: 2026-08-27. Directory: `apbc2026/journal_version/additions/multiome/`.
Nothing in `apbc2026/` outside this directory, and no existing data directory, was modified.

## Why

The main benchmark (`apbc2026/tools/benchmark_aggregation.py`) scores each aggregation method
against 22-41 held-out canonical markers. That is a thin ground truth and it depends on the
cell labels being right. In a 10x Multiome experiment RNA and ATAC are measured in the same
cell, so (i) labels can come from RNA alone, independent of the ATAC matrix being scored, and
(ii) the ATAC per-type profile of every method can be compared genome-wide with the RNA
per-type profile of the same cells.

## Data and processing

| item | value |
|---|---|
| dataset | 10x Genomics `pbmc_granulocyte_sorted_10k` (Cell Ranger ARC 2.0.0), downloaded 2026-08-27 to `data/raw/multiome_pbmc10k/` (filtered_feature_bc_matrix.h5 192 MB, atac_fragments.tsv.gz 2.4 GB + .tbi) |
| cells in filtered matrix | 11,898 (36,601 genes + 143,887 peaks); median RNA UMI 3,778 |
| ATAC gene activity | built from fragments with SnapATAC2 2.9 for the 11,898 filtered barcodes (`01_build_gene_activity.py`). **proximal** = `make_gene_matrix(upstream=2000, downstream=0, include_gene_body=True)` on GENCODE v38 (Signac model; 59,385 genes). **enhancer** = 5 kb tiles x exp(-d/5 kb) distance weights within 100 kb of the gene body, protein-coding genes only (ArchR-style; 19,871 genes). Same code paths as `scripts/analysis/90`. |
| true depth | `n_fragment` recorded by SnapATAC2 at import (median 17,298 per cell), never a matrix row sum |
| RNA labels | `02_label_rna.py`: scanpy normalize_total/log1p/2000 HVG/30 PC/leiden(res 1.0) on Gene Expression only -> 18 clusters, named by z-scored cluster means of a 31-gene labelling panel (`labelling_genes.txt`). Clusters scoring > 0.4 in a second lineage group (myeloid+T, B+myeloid, B+T) were dropped as doublets (648 cells); Plasma (41 cells) < 100-cell floor. Platelets are absent from this granulocyte-sorted sample (PPBP/PF4 ~ 0), so no platelet type. |
| kept types (n) | T_CD4 3,542; Mono_CD14 3,143; T_CD8 2,273; B 794; NK 682; Mono_CD16 442; DC 221; pDC 112 -> **11,209 cells, 8 types** |

Outputs: `multiome_benchmark.csv` (one row per model x method), `multiome_markers.csv`,
`multiome_depth_by_type.csv`, `multiome_annotation.csv` (per-type recall + confusion),
`multiome_equal_n.csv`, `multiome_rna_truth_share.csv`, `fig_multiome_benchmark.png/.pdf`.
Logs: `03.log`, `05.log`. Processed matrices: `data/processed/multiome_pbmc10k/` (new directory, 3.5 GB).

## Depth spread between RNA-defined types (say it up front: it is small)

| type | n | median fragments | median RNA UMI |
|---|---|---|---|
| T_CD4 | 3,542 | 16,001 | 3,241 |
| T_CD8 | 2,273 | 16,434 | 3,296 |
| NK | 682 | 16,579 | 3,529 |
| Mono_CD16 | 442 | 17,216 | 6,699 |
| Mono_CD14 | 3,143 | 18,515 | 4,551 |
| B | 794 | 19,497 | 3,160 |
| DC | 221 | 26,038 | 8,262 |
| pDC | 112 | 26,573 | 5,607 |

**Depth spread = 1.66x** (pDC / T_CD4), against 3-35x in the tumour cohorts of the main benchmark.
The two deepest types (DC, pDC) are also the two rarest and the two most transcriptionally
distinct, so depth and biology are confounded across these 8 points: the RNA truth's own
argmax share correlates with ATAC depth at rho = 0.43 (all genes) / 0.67 (informative genes).
A depth-driven failure mode cannot be large here; the value of this dataset is the
genome-wide RNA truth, not the depth stress.

## Read-outs (7 methods x 2 models; `multiome_benchmark.csv`)

Genes "expressed in both modalities": ATAC column sum > 0 and RNA detected in >= 1 % of cells
(11,920 proximal / 11,134 enhancer). "RNA-informative": RNA argmax type >= 2x the second type
in the RNA CPM-mean profile (1,338 / 1,215 genes) - on the rest the RNA argmax is itself noise.
Held-out panel: 51 PBMC markers disjoint from the 31 labelling genes (pan-T and pan-myeloid
genes accept any type of that lineage); permutation null 8.6/51.
Annotation: references from a stratified 50 % training split (5,603 cells); each held-out cell
(5,606) is log1p(CPM)-normalised and assigned to the reference row with the highest
correlation (Pearson on values; Spearman on ranks).

### Proximal model (Signac gene body + 2 kb)

| method | top share (type) | rho(depth, share) | RNA argmax agree, informative (all) | per-gene rho, informative (all) | markers /51 | annot acc Pearson (balanced) | annot acc Spearman (balanced) |
|---|---|---|---|---|---|---|---|
| raw mean | **34.1 % (pDC)** | +0.93 (p=0.001) | 0.621 (0.329) | 0.545 (0.394) | 41 | 0.880 (0.901) | 0.663 (0.822) |
| z-scale mean | **34.5 % (pDC)** | +0.93 | 0.622 (0.328) | 0.545 (0.393) | 41 | **0.148 (0.433)** | **0.082 (0.256)** |
| log1p(CPM) mean | 22.6 % (pDC) | +0.88 | 0.649 (0.350) | 0.553 (0.429) | 41 | 0.880 (0.903) | 0.709 (0.841) |
| ArchR log2 mean | 22.6 % (pDC) | +0.88 | 0.649 (0.350) | 0.553 (0.429) | 41 | 0.880 (0.903) | 0.709 (0.841) |
| pseudobulk CPM | 20.9 % (Mono_CD14) | +0.83 | 0.632 (0.340) | 0.552 (0.425) | 41 | 0.880 (0.901) | 0.663 (0.823) |
| depth-matched | 22.0 % (pDC) | +0.88 | 0.649 (0.353) | 0.558 (0.429) | 40 | 0.882 (0.894) | 0.658 (0.827) |
| **CPM mean (ours)** | **19.3 % (pDC)** | +0.90 | 0.632 (0.332) | 0.548 (0.409) | 41 | 0.880 (0.870) | 0.729 (0.843) |

### Enhancer model (ArchR-style distance-weighted tiles)

| method | top share (type) | rho(depth, share) | RNA argmax agree, informative (all) | per-gene rho, informative (all) | markers /51 | annot acc Pearson (balanced) | annot acc Spearman (balanced) |
|---|---|---|---|---|---|---|---|
| raw mean | **44.9 % (pDC)** | +0.71 (p=0.047) | 0.650 (0.335) | 0.559 (0.392) | 40 | 0.884 (0.894) | 0.725 (0.807) |
| z-scale mean | **45.2 % (pDC)** | +0.71 | 0.653 (0.335) | 0.560 (0.392) | 40 | **0.330 (0.541)** | **0.159 (0.464)** |
| log1p(CPM) mean | 22.5 % (pDC) | +0.24 | 0.681 (0.359) | 0.578 (0.438) | 41 | 0.884 (0.896) | 0.869 (0.869) |
| ArchR log2 mean | 22.5 % (pDC) | +0.24 | 0.681 (0.359) | 0.578 (0.438) | 41 | 0.884 (0.896) | 0.869 (0.869) |
| pseudobulk CPM | 17.8 % (pDC) | +0.19 | 0.662 (0.353) | 0.577 (0.433) | 40 | 0.884 (0.894) | 0.725 (0.807) |
| depth-matched | 23.3 % (pDC) | +0.48 | 0.684 (0.361) | 0.584 (0.439) | 40 | 0.885 (0.883) | 0.841 (0.855) |
| **CPM mean (ours)** | **18.1 % (pDC)** | +0.38 | 0.662 (0.339) | 0.573 (0.420) | 41 | 0.881 (0.855) | 0.862 (0.845) |

### The RNA truth's own argmax concentration (`multiome_rna_truth_share.csv`)

On the same shared genes the RNA CPM-mean profile puts **22.5 %** of all argmaxima (33 % of
informative-gene argmaxima) in pDC, 16 % in NK, 8-12 % in each remaining type. So a pDC share
near 20 % is what a faithful ATAC profile *should* show.

### Equal-n control (`multiome_equal_n.csv`; every type subsampled to 112 cells, 20 draws, medians)

| method | proximal: pDC share / rho / top = deepest | enhancer: pDC share / rho / top = deepest |
|---|---|---|
| raw mean | 26.2 % / +0.93 / 19 of 20 | 42.1 % / +0.81 / 20 of 20 |
| z-scale mean | 28.3 % / +0.90 / 19 of 20 | 42.1 % / +0.81 / 20 of 20 |
| log1p(CPM) mean | 17.8 % / +0.71 / 7 of 20 | 20.0 % / +0.20 / 20 of 20 |
| pseudobulk CPM | 14.6 % / +0.30 / 0 of 20 | 16.0 % / +0.06 / 20 of 20 |
| depth-matched | 15.3 % / +0.80 / 0 of 20 | 19.3 % / +0.48 / 17 of 20 |
| CPM mean | 15.6 % / +0.69 / 0 of 20 | 16.0 % / +0.12 / 20 of 20 |

The raw-mean inflation of the deepest type survives equalising cell numbers, so it is not a
small-n artefact; under the normalised methods the pDC share drops to the 15-20 % that the
RNA truth predicts.

## What the multiome data establish

1. **Raw and z-scaled means over-concentrate argmaxima in the deepest type even at a 1.7x
   depth spread**: 34 % (proximal) and 45 % (enhancer) of all genes peak in pDC, versus the
   22.5 % that the RNA truth of the same cells shows. Every per-cell-normalised method
   (log1p(CPM), ArchR log2, pseudobulk CPM, depth-matched, CPM mean) returns 18-23 %, i.e.
   the RNA-truth level. This is the first place in the project where the "expected" argmax
   share is measured rather than assumed to be uniform.
2. **Genome-wide fidelity to RNA is essentially method-independent** once the profile is
   normalised: ATAC argmax = RNA argmax on 62-68 % of RNA-informative genes (33-36 % of all
   shared genes), mean per-gene Spearman across the 8 types 0.55-0.58 (informative) /
   0.39-0.44 (all). Raw mean is consistently the lowest (0.621/0.650 agreement; 0.394/0.392
   rho on all genes) and log-type means or depth matching the highest (0.649/0.684;
   0.429/0.439), but the gap is 0.03 in agreement and 0.04 in rho. Enhancer-aware activity is
   slightly closer to RNA than proximal activity on every read-out (+0.03 agreement).
3. **The held-out marker test does not separate methods here**: 40-41/51 under all 14
   combinations (null 8.6). The consistent misses are panel imperfections (SKAP1, CD247 -> NK,
   which does express them; CD36 -> Mono_CD16; CLEC9A is a cDC1 gene inside a mostly-cDC2
   cluster; DERL3/ITM2C -> B). What differs is *which* markers the raw mean loses:
   CSF3R, CD36, CCR2, CX3CR1, BLK, CD68, CSF1R all go to pDC (deepest type) under raw/z-scale
   means, whereas the normalised methods lose ENHO/CLEC9A/DERL3 instead. Same count, different
   failure signature - consistent with the depth mechanism, but PBMC cannot make it bite.
4. **Annotation of held-out cells**: 0.88 accuracy (0.85-0.90 balanced) for every method
   except the z-scaled mean (0.15 proximal / 0.33 enhancer with Pearson; 0.08 / 0.16 with
   Spearman), which is the one clearly broken reference: z-scaling per gene across cells makes
   the reference rows uncorrelated with any single cell's profile (Mono_CD14, Mono_CD16, T_CD8
   recall 0.00-0.01). Pearson accuracy is identical for raw mean and pseudobulk CPM because
   the two references differ only by a per-row constant, and near-identical for all other
   unscaled methods because Pearson is dominated by the same high-count genes. With Spearman
   (rank) correlation the CPM mean and the log-type means annotate best (0.73/0.86 and
   0.71/0.87 for proximal/enhancer) and the raw mean worst (0.66/0.73). Residual errors are
   biological: T_CD8 -> NK (296/1,137 with the proximal CPM reference) and T_CD4 <-> T_CD8.
   The enhancer-aware reference reduces the T_CD8/NK confusion under Spearman
   (T_CD8 recall 0.67 -> 0.80).

## Caveats

- **Depth spread is 1.66x.** This dataset cannot reproduce the 60-98 % argmax capture seen in
  COAD/HCC; it shows the direction (raw mean 1.5-2x the RNA-truth share) at a small spread. The
  claim that depth-normalisation is necessary rests on the tumour cohorts; the multiome result
  is that it costs nothing genome-wide and restores the RNA-truth argmax distribution.
- pDC (n = 112) is simultaneously the deepest, the rarest, and the most distinct type, so the
  8-point Spearman between depth and share is confounded (the RNA truth itself has rho 0.43-0.67
  with depth). The equal-n control addresses n; nothing here can separate depth from biology
  across only 8 types. Do not headline rho(depth, share) from this dataset.
- Labels are leiden clusters named by a marker panel, not an expert annotation; 689 cells were
  dropped (doublet-scoring clusters and 41 plasma cells). T_CD4 vs T_CD8 vs NK boundaries in the
  RNA labels are themselves imperfect and cap the annotation accuracy.
- The RNA truth is a CPM mean of RNA counts, so the argmax-agreement read-out compares ATAC
  argmax to an RNA argmax computed with the very normalisation we recommend. The per-gene
  Spearman and the annotation task do not share this dependence.
- "Expressed in both" uses a 1 %-of-cells RNA detection floor; the "informative" subset (2x
  margin) is where the RNA argmax is meaningful. Both thresholds are arbitrary and stated in
  `03_benchmark_multiome.py`.
- The enhancer model covers protein-coding genes only (19,871) while the proximal model
  covers all GENCODE genes (59,385); shared-gene sets therefore differ (11,134 vs 11,920).
- Depth-matched thinning uses the matrix row sum as its target (as in the main benchmark),
  not `n_fragment`.
- Only one multiome dataset (healthy PBMC, one donor). A tumour multiome from GEO was not
  attempted: a quick search did not surface a small fragment-file deposit, and multi-GB GEO
  downloads were out of the time budget. Candidates to try next: any 10x Multiome tumour
  deposit with `atac_fragments.tsv.gz` (search GEO for "Chromium Single Cell Multiome" +
  tumour type); the same five scripts run unchanged given an h5 + fragments pair.

## Reproduce

```
python3 01_build_gene_activity.py all      # ~5 min, writes data/processed/multiome_pbmc10k/
python3 02_label_rna.py                    # rna_labels.csv, labelling_genes.txt
python3 03_benchmark_multiome.py           # ~3 min, the four multiome_*.csv
python3 04_figure.py                       # fig_multiome_benchmark.png/.pdf
python3 05_equal_n_check.py                # multiome_equal_n.csv, multiome_rna_truth_share.csv
```
Aggregation functions are imported from `apbc2026/tools/benchmark_aggregation.py` (unchanged).
Peak RAM < 8 GB.
