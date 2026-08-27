# Matched RNA–ATAC ground truth (10x PBMC Multiome)

The main benchmark scores each aggregation strategy against 22–41 held-out
canonical markers. That is a thin ground truth and it depends on the cell labels
being right. In a 10x Multiome experiment RNA and ATAC are measured in the same
cell, so labels can come from RNA alone — independent of the ATAC matrix being
scored — and the ATAC per-type profile of every strategy can be compared
genome-wide with the RNA per-type profile of the same cells.

> Naming: this directory predates the manuscript's final terminology. In file and
> column names, `depth` = coverage, `proximal` = the count matrix, `enhancer` =
> the weighted matrix, and `pseudobulk_cpm` = pooled CPM.

## Data

| item | value |
|---|---|
| dataset | 10x Genomics `pbmc_granulocyte_sorted_10k` (Cell Ranger ARC 2.0.0), one donor |
| cells in filtered matrix | 11,898 (36,601 genes + 143,887 peaks); median RNA UMI 3,778 |
| count matrix | SnapATAC2 2.9 `make_gene_matrix(upstream=2000, downstream=0, include_gene_body=True)` on GENCODE v38 — the Signac window; 59,385 genes |
| weighted matrix | 5 kb tiles × exp(−d/5 kb) within 100 kb of the gene body, protein-coding only — ArchR-style; 19,871 genes |
| coverage | `n_fragment` recorded by SnapATAC2 at import (median 17,298 per cell), never a matrix row sum |
| RNA labels | scanpy normalize_total / log1p / 2000 HVG / 30 PC / leiden(res 1.0) on Gene Expression only → 18 clusters, named by z-scored cluster means of a 31-gene panel (`labelling_genes.txt`). Clusters scoring > 0.4 in a second lineage group were dropped as doublets (648 cells); Plasma (41 cells) fell below a 100-cell floor. Platelets are absent from this granulocyte-sorted sample. |
| kept types (n) | T_CD4 3,542; Mono_CD14 3,143; T_CD8 2,273; B 794; NK 682; Mono_CD16 442; DC 221; pDC 112 → **11,209 cells, 8 types** |

## Coverage spread between RNA-defined types — it is small

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

**Coverage spread = 1.66×** (pDC / T_CD4), against 3–35× in the tumour cohorts.
The two highest-coverage types (DC, pDC) are also the two rarest and the two most
transcriptionally distinct, so coverage and biology are confounded across these
eight points: the RNA truth's own argmax share correlates with ATAC coverage at
ρ = 0.43 (all genes) / 0.67 (informative genes). A coverage-driven failure cannot
be large here; the value of this dataset is the genome-wide RNA truth, not the
coverage stress.

## Read-outs

Genes "expressed in both modalities": ATAC column sum > 0 and RNA detected in ≥ 1%
of cells (11,920 count / 11,134 weighted). "RNA-informative": RNA argmax type ≥ 2×
the second type in the RNA CPM-mean profile (1,338 / 1,215 genes) — on the rest the
RNA argmax is itself noise. Held-out panel: 51 PBMC markers disjoint from the 31
labelling genes; permutation null 8.6/51. Annotation: references from a stratified
50% training split (5,603 cells); each held-out cell (5,606) is log-normalised and
assigned to the reference row with the highest correlation.

### Count matrix

| strategy | top share (type) | ρ(coverage, share) | RNA argmax agree, informative (all) | per-gene ρ, informative (all) | markers /51 | annot acc Pearson (balanced) | annot acc Spearman (balanced) |
|---|---|---|---|---|---|---|---|
| raw mean | **34.1% (pDC)** | +0.93 (p=0.001) | 0.621 (0.329) | 0.545 (0.394) | 41 | 0.880 (0.901) | 0.663 (0.822) |
| z-scaled mean | **34.5% (pDC)** | +0.93 | 0.622 (0.328) | 0.545 (0.393) | 41 | **0.148 (0.433)** | **0.082 (0.256)** |
| log-normalised mean | 22.6% (pDC) | +0.88 | 0.649 (0.350) | 0.553 (0.429) | 41 | 0.880 (0.903) | 0.709 (0.841) |
| ArchR log2 mean | 22.6% (pDC) | +0.88 | 0.649 (0.350) | 0.553 (0.429) | 41 | 0.880 (0.903) | 0.709 (0.841) |
| pooled CPM | 20.9% (Mono_CD14) | +0.83 | 0.632 (0.340) | 0.552 (0.425) | 41 | 0.880 (0.901) | 0.663 (0.823) |
| coverage-matched | 22.0% (pDC) | +0.88 | 0.649 (0.353) | 0.558 (0.429) | 40 | 0.882 (0.894) | 0.658 (0.827) |
| **CPM mean** | **19.3% (pDC)** | +0.90 | 0.632 (0.332) | 0.548 (0.409) | 41 | 0.880 (0.870) | 0.729 (0.843) |

### Weighted matrix

| strategy | top share (type) | ρ(coverage, share) | RNA argmax agree, informative (all) | per-gene ρ, informative (all) | markers /51 | annot acc Pearson (balanced) | annot acc Spearman (balanced) |
|---|---|---|---|---|---|---|---|
| raw mean | **44.9% (pDC)** | +0.71 (p=0.047) | 0.650 (0.335) | 0.559 (0.392) | 40 | 0.884 (0.894) | 0.725 (0.807) |
| z-scaled mean | **45.2% (pDC)** | +0.71 | 0.653 (0.335) | 0.560 (0.392) | 40 | **0.330 (0.541)** | **0.159 (0.464)** |
| log-normalised mean | 22.5% (pDC) | +0.24 | 0.681 (0.359) | 0.578 (0.438) | 41 | 0.884 (0.896) | 0.869 (0.869) |
| ArchR log2 mean | 22.5% (pDC) | +0.24 | 0.681 (0.359) | 0.578 (0.438) | 41 | 0.884 (0.896) | 0.869 (0.869) |
| pooled CPM | 17.8% (pDC) | +0.19 | 0.662 (0.353) | 0.577 (0.433) | 40 | 0.884 (0.894) | 0.725 (0.807) |
| coverage-matched | 23.3% (pDC) | +0.48 | 0.684 (0.361) | 0.584 (0.439) | 40 | 0.885 (0.883) | 0.841 (0.855) |
| **CPM mean** | **18.1% (pDC)** | +0.38 | 0.662 (0.339) | 0.573 (0.420) | 41 | 0.881 (0.855) | 0.862 (0.845) |

### The RNA truth's own argmax concentration

On the same shared genes the RNA CPM-mean profile puts **22.5%** of all argmaxima
(33% of informative-gene argmaxima) in pDC, 16% in NK, and 8–12% in each remaining
type. A pDC share near 20% is therefore what a faithful ATAC profile *should* show.

### Equal-n control

Every type subsampled to 112 cells, 20 draws, medians.

| strategy | count: pDC share / ρ / top = highest-coverage | weighted: pDC share / ρ / top = highest-coverage |
|---|---|---|
| raw mean | 26.2% / +0.93 / 19 of 20 | 42.1% / +0.81 / 20 of 20 |
| z-scaled mean | 28.3% / +0.90 / 19 of 20 | 42.1% / +0.81 / 20 of 20 |
| log-normalised mean | 17.8% / +0.71 / 7 of 20 | 20.0% / +0.20 / 20 of 20 |
| pooled CPM | 14.6% / +0.30 / 0 of 20 | 16.0% / +0.06 / 20 of 20 |
| coverage-matched | 15.3% / +0.80 / 0 of 20 | 19.3% / +0.48 / 17 of 20 |
| CPM mean | 15.6% / +0.69 / 0 of 20 | 16.0% / +0.12 / 20 of 20 |

The raw-mean inflation of the highest-coverage type survives equalising cell
numbers, so it is not a small-n artefact; under the normalised strategies the pDC
share drops to the 15–20% the RNA truth predicts.

## What these data establish

1. **Raw and z-scaled means over-concentrate argmaxima in the highest-coverage type
   even at a 1.7× spread**: 34% (count) and 45% (weighted) of all genes peak in pDC,
   against the 22.5% the RNA truth of the same cells shows. Every per-cell-normalised
   strategy returns 18–23%, i.e. the RNA-truth level. This is the one place in the
   project where the expected argmax share is measured rather than assumed uniform.
2. **Genome-wide fidelity to RNA is essentially strategy-independent** once the
   profile is normalised: ATAC argmax = RNA argmax on 62–68% of RNA-informative genes
   (33–36% of all shared genes), mean per-gene Spearman across the eight types
   0.55–0.58 (informative) / 0.39–0.44 (all). The raw mean is consistently lowest and
   log-type means or coverage matching highest, but the gap is 0.03 in agreement and
   0.04 in ρ. The weighted matrix is slightly closer to RNA than the count matrix on
   every read-out (+0.03 agreement).
3. **The held-out marker test does not separate strategies here**: 40–41/51 under all
   14 combinations (null 8.6). The consistent misses are panel imperfections. What
   differs is *which* markers the raw mean loses — CSF3R, CD36, CCR2, CX3CR1, BLK,
   CD68 and CSF1R all go to pDC under raw and z-scaled means, whereas the normalised
   strategies lose ENHO/CLEC9A/DERL3 instead. Same count, different failure signature:
   consistent with the coverage mechanism, but PBMC cannot make it bite.
4. **Annotation of held-out cells**: 0.88 accuracy (0.85–0.90 balanced) for every
   strategy except the z-scaled mean (0.15 count / 0.33 weighted with Pearson; 0.08 /
   0.16 with Spearman), the one clearly broken reference — z-scaling per gene across
   cells makes the reference rows uncorrelated with any single cell's profile
   (Mono_CD14, Mono_CD16, T_CD8 recall 0.00–0.01). Pearson accuracy is identical for
   the raw mean and pooled CPM because the two references differ only by a per-row
   constant. With Spearman the CPM mean and the log-type means annotate best
   (0.73/0.86 and 0.71/0.87) and the raw mean worst (0.66/0.73). Residual errors are
   biological: T_CD8 → NK and T_CD4 ↔ T_CD8. The weighted reference reduces the
   T_CD8/NK confusion under Spearman (T_CD8 recall 0.67 → 0.80).

## Caveats

- **Coverage spread is 1.66×.** This dataset cannot reproduce the 60–98% argmax
  capture seen in COAD and HCC; it shows the direction (raw mean 1.5–2× the RNA-truth
  share) at a small spread. The claim that per-cell normalisation is necessary rests
  on the tumour cohorts; the multiome result is that it costs nothing genome-wide and
  restores the RNA-truth argmax distribution.
- pDC (n = 112) is simultaneously the highest-coverage, the rarest and the most
  distinct type, so the eight-point Spearman between coverage and share is confounded
  — the RNA truth itself has ρ = 0.43–0.67 with coverage. The equal-n control
  addresses n; nothing here can separate coverage from biology across only eight
  types. Do not headline ρ(coverage, share) from this dataset.
- Labels are leiden clusters named by a marker panel, not an expert annotation; 689
  cells were dropped. T_CD4 / T_CD8 / NK boundaries in the RNA labels are themselves
  imperfect and cap the annotation accuracy.
- The RNA truth is a CPM mean of RNA counts, so the argmax-agreement read-out compares
  the ATAC argmax to an RNA argmax computed with the very normalisation recommended
  here. The per-gene Spearman and the annotation task do not share this dependence.
- "Expressed in both" uses a 1%-of-cells RNA detection floor; the "informative" subset
  uses a 2× margin. Both thresholds are arbitrary and stated in `03_benchmark_multiome.py`.
- The weighted model covers protein-coding genes only (19,871) while the count model
  covers all GENCODE genes (59,385), so the shared-gene sets differ (11,134 vs 11,920).
- Coverage-matched thinning uses the matrix row sum as its target, as in the main
  benchmark, not `n_fragment`.
- One multiome dataset only (healthy PBMC, one donor).

## Files

`multiome_benchmark.csv` (one row per model × strategy), `multiome_markers.csv`,
`multiome_depth_by_type.csv`, `multiome_annotation.csv` (per-type recall and
confusion), `multiome_equal_n.csv`, `multiome_rna_truth_share.csv`,
`rna_labels.csv`, `rna_cluster_scores.csv`, `labelling_genes.txt`.
