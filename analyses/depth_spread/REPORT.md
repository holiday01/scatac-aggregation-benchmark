# Controlled depth-spread experiment — 2026-08-27

Addition to the aggregation benchmark (`tools/benchmark_aggregation.py`). The benchmark shows
the effect of the aggregation method at the three natural depth spreads of the cohorts
(HCC ~2x, LUAD ~3.5x, COAD ~37x). Here the spread is *imposed* on two count matrices so the
whole curve is measured under one design.

Script: `depth_spread_experiment.py` (this directory). Reads only the two h5ad files and
imports `PANEL`/`percell_scale` from the benchmark script; writes only into this directory.
Runtime: HCC 2 min, LUAD 5 min (single process each, peak RSS < 8 GB).

## Design

| | HCC proximal | LUAD proximal |
|---|---|---|
| matrix | 12,029 cells x 59,385 genes, 6 types (relabel rule of the benchmark applied: "B_cell" -> Endothelial_stromal, "DC" -> B_cell) | 58,112 cells x 19,930 genes, 6 types |
| common target depth (p20 of row sums) | 5,630 | 2,679 |
| natural spread (median row sum, deepest / shallowest) | 1.84x | 3.57x |
| deep-type choices | B_cell (n=416), NK_cytotoxic_T (n=3,677) | B_cell (n=9,484), T_NK (n=12,221) |
| held-out panel | 22 markers, 4 lineages | 36 markers, 6 lineages |

1. **Equalise.** Every cell is binomially thinned to the target (cells already below it are kept),
   so every type starts at a spread of 1.00x (realised 1.003x HCC, 1.015x LUAD).
2. **Impose.** The deep type stays at the target; the other five are thinned by factors
   `S^(-i/5)`, i = 1..5, i.e. evenly log-spaced between 1 and 1/S, their order a random
   permutation per seed. S = 1, 2, 5, 10, 20, 40; 2 deep types x 3 seeds = 6 replicates per
   spread, 36 conditions per cohort (+1 natural, unthinned). Both thinnings compose into one
   binomial draw from the original counts. Realised spreads: 2.0, 5.0-5.1, 10.1-10.3, 20.3-20.5,
   40.9-41.7x.
3. **Score.** Five methods, definitions identical to the benchmark: raw_mean, zscale_mean
   (scanpy.pp.scale max_value=10; computed with a sparse-exact shortcut, max |diff| 4e-5 vs the
   benchmark's dense code), lognorm_mean (normalize_total + log1p), pseudobulk_cpm, cpm_mean.
   Read-outs: share of expressed genes whose argmax falls in one type ("top share") and which
   type; whether the deep type wins; Spearman of per-type share vs realised median depth and vs
   detection rate; held-out marker concordance with a 2,000-shuffle permutation null.
4. **Mechanism.** Per type and condition: detection rate (fraction of non-zero genes per cell),
   the mean over expressed genes of the log1p(CPM), CPM and raw profiles, and the argmax share
   under each method.

## Headline: top argmax share vs imposed spread (mean of 6 replicates; in brackets: replicates in which the deep type wins, /6)

**HCC**

| method | 1x | 2x | 5x | 10x | 20x | 40x | natural 1.8x |
|---|---|---|---|---|---|---|---|
| raw mean | 30% (0) | 36% (3) | 48% (6) | 59% (6) | 62% (6) | 68% (6) | 67% (Endothelial_stromal) |
| z-scale mean | 31% (0) | 35% (3) | 47% (6) | 57% (6) | 61% (6) | 67% (6) | 65% |
| log1p(CPM) mean | 31% (0) | 28% (2) | 32% (6) | 38% (4) | 43% (5) | 50% (4) | 34% (Hepatocyte) |
| pseudobulk CPM | 30% (0) | 29% (0) | 28% (0) | 26% (0) | 24% (0) | 25% (0) | 32% (Hepatocyte) |
| **CPM mean** | 29% (0) | 29% (0) | 27% (0) | 26% (0) | 24% (0) | 24% (0) | 30% (Hepatocyte) |

**LUAD**

| method | 1x | 2x | 5x | 10x | 20x | 40x | natural 3.6x |
|---|---|---|---|---|---|---|---|
| raw mean | 28% (3*) | 69% (6) | 85% (6) | 87% (6) | 92% (6) | 93% (6) | 66% (B_cell) |
| z-scale mean | 28% (3*) | 68% (6) | 83% (6) | 85% (6) | 90% (6) | 92% (6) | 65% |
| log1p(CPM) mean | 31% (3*) | 54% (6) | 72% (6) | 75% (6) | 81% (6) | 83% (6) | 37% (T_NK) |
| pseudobulk CPM | 25% (3*) | 24% (3*) | 23% (3*) | 22% (1) | 23% (0) | 25% (0) | 24% (B_cell) |
| **CPM mean** | 24% (3*) | 24% (3*) | 22% (3*) | 22% (0) | 23% (0) | 25% (0) | 25% (T_NK) |

\* At 1x (and under the depth-invariant methods at 2-5x) the winner on LUAD is T_NK under every
method, i.e. the residual, non-depth structure of the matrix; T_NK is also one of the two
"deep" choices, hence 3/6 "deep wins" that have nothing to do with depth. Under CPM mean the
winner on HCC is Hepatocyte in 36/36 conditions regardless of which type is deep; on LUAD it is
T_NK at 1-5x and drifts to Macrophage/Epithelial at 20-40x, where the shallowest types have
median row sums of 65-140 counts (see caveats).

Min-max over the 6 replicates (see `depth_spread_summary.csv`): raw mean at 40x spans 59-79%
(HCC) and 89-96% (LUAD); CPM mean at 40x spans 22-27% in both.

**Spearman of per-type argmax share vs realised median depth** (mean over replicates; NaN at 1x):
raw / z-scale +0.71 (2x) -> +0.89-0.96 (5-40x) on both cohorts; log1p(CPM) +0.54 -> +0.89 (HCC),
+0.69 -> +0.98 (LUAD); pseudobulk CPM and CPM mean -0.06 -> -0.41 (HCC) and +0.17 -> -0.69
(LUAD) -- zero-to-negative, i.e. the residual share drifts *away* from the deep type as its
competitors lose their detection floor.

## Held-out marker concordance (mean of 6 replicates)

| cohort / method | 1x | 2x | 5x | 10x | 20x | 40x |
|---|---|---|---|---|---|---|
| HCC (of 22) raw mean | 19.2 | 19.0 | 16.2 | 10.5 | 10.3 | 11.2 |
| HCC z-scale mean | 19.2 | 19.3 | 16.7 | 10.8 | 10.7 | 11.5 |
| HCC log1p(CPM) mean | 19.2 | 19.7 | 18.8 | 15.2 | 14.7 | 16.0 |
| HCC pseudobulk CPM | 19.2 | 19.5 | 18.7 | 18.7 | 19.0 | 18.3 |
| HCC CPM mean | 19.2 | 19.0 | 18.8 | 18.7 | 18.8 | 18.5 |
| LUAD (of 36) raw mean | 24.0 | 21.2 | 15.0 | 13.0 | 11.5 | 10.3 |
| LUAD z-scale mean | 23.7 | 21.5 | 15.5 | 14.2 | 12.2 | 11.0 |
| LUAD log1p(CPM) mean | 23.7 | 23.3 | 19.8 | 17.3 | 17.3 | 14.2 |
| LUAD pseudobulk CPM | 23.8 | 23.5 | 23.3 | 24.2 | 23.2 | 23.3 |
| LUAD CPM mean | 23.8 | 23.5 | 23.0 | 24.5 | 22.8 | 22.8 |

At 40x the raw mean is at or near the permutation null (perm p up to 0.04 HCC, 0.29 LUAD; null
mean 5.5/22 and 6/36) while CPM mean stays at p < 0.001 in every condition. Per lineage, the
loss under raw mean is concentrated in whichever lineages were thinned: HCC macrophage markers
fall from 6/6 to 0.4/6 (mean over replicates), LUAD macrophage 6/6 -> 1/6, epithelial
4/6 -> 1.5/6, T_NK 6/6 -> 3.2/6. Under CPM mean the same lineages keep 5.5/6, 4.7/6, 4.7/6,
6/6. (The LUAD fibroblast and endothelial panels are poor under *every* method already at
1x -- 0/6 and 2/6 -- a property of the panel or the labels, not of the method.)

## Mechanism (`depth_spread_mechanism_summary.csv`, `depth_spread_by_type.csv`)

The per-type detection rate scales with the imposed depth: its range across types is 1.07
(1x) -> 4.2 (5x) -> 8.1 (10x) -> 32 (40x) on HCC, 1.3 -> 30 on LUAD (Spearman detection vs
depth = 1.00 at every spread >= 2x). The mean of the log1p(CPM) profile follows the detection
rate (rho = 0.95-1.00 in every condition) with a compressed range: 1.04 -> 8.6 (HCC) and
1.19 -> 10.5 (LUAD), because a cell contributes 0 to the mean of a gene it does not detect,
however it is normalised. The mean of the CPM profile is 1.000 across types at every spread
(each cell sums to 1e4), and the mean of the raw profile scales 1:1 with depth (range 40.6 at
40x). The argmax share under log1p(CPM) tracks its own profile mean (rho 0.54 -> 0.91 HCC,
0.88 -> 0.98 LUAD) and under raw mean tracks depth (0.64 -> 0.96).

Worked example, HCC deep = B_cell, seed 0, 40x: the B-cell cluster keeps a detection rate of
0.070 and the others fall to 0.002-0.038; raw share 75%, log1p(CPM) share 57%, CPM share 14%.

## Interpretation

1. The effect is monotone in the spread and saturates: raw and z-scaled means put the argmax
   in the deep type in 36/36 replicates at >= 5x, reaching 68% (HCC) and 93% (LUAD) at 40x. The
   z-scale is not a fix: it is within 2 points of the raw mean at every spread, because the
   per-gene centring and clipping do not change which type has the highest mean.
2. log1p(CPM), the scanpy/Seurat/ArchR convention, halves the effect but does not remove it:
   50% (HCC) and 83% (LUAD) at 40x, deep type winning in 25/36 (HCC) and 36/36 (LUAD) replicates
   at >= 2x. Its residual depth dependence is the detection-rate term.
3. CPM mean and pseudobulk CPM are flat at 22-30% across the whole range, indistinguishable
   from each other, with the deep type never winning at >= 10x, and their marker concordance is
   unchanged from 1x to 40x. The two methods behave identically here because equal-depth cells
   contribute equally to both; they differ only when depth varies *within* a type (the
   pseudobulk weights cells by depth), which the thinning removes.
4. The cohort sensitivity differs: at 2x LUAD raw mean is already at 69% while HCC is at 36%.
   The difference is the amount of genuine lineage structure in the matrix: on HCC the
   Hepatocyte profile wins ~30% of genes at 1x under every method (a real, strongly specific
   programme that a 13% depth edge does not overturn), whereas the LUAD gene-activity profiles
   are nearly tied at equal depth, so a small depth edge decides most genes.
5. The natural matrices (open markers on the figure) do not sit on the imposed curve at their
   median spread: HCC raw mean 67% at 1.8x (imposed 2x: 36%), LUAD log1p(CPM) 37% at 3.6x
   (imposed: 65%). See caveat 1.

## Caveats

1. **The median spread understates what the raw mean sees.** The raw mean responds to the
   per-type *mean* depth and its tail. Natural HCC: median spread 1.84x but mean spread 2.86x
   (Endothelial_stromal mean 22,752, p90 48,462, vs CD4_T mean 7,963); natural LUAD: B_cell has
   the highest *mean* row sum (16,728 vs T_NK 14,341) although T_NK has the highest median --
   which is why the natural raw-mean argmax goes to B_cell. The imposed design caps every cell at
   the p20 target, so mean = median and the tails are removed; natural matrices with heavy
   within-type tails will behave like a larger spread than their median suggests. Conversely the
   natural LUAD log1p(CPM) value (37%) is lower than the imposed 3.6x curve (65%) because in the
   natural matrix five of six types sit within 1.5x of each other and only Fibroblast is far
   below, whereas the imposed profile spreads all five competitors evenly. A single scalar
   "spread" is therefore a coarse descriptor; the direction of every conclusion above is
   unaffected, the exact share at a given spread is design-dependent.
2. Thinning is binomial on the gene-activity counts, not on fragments; it is exact for a
   Poisson/binomial sampling model but does not reproduce depth-dependent changes in the
   fragment-to-gene assignment (duplicates, peak calling). Depth is measured as the matrix row
   sum, which the benchmark deliberately avoids; here it is the imposed quantity, so it is the
   right scale.
3. At 20-40x the shallowest types have median row sums of 65-300 counts (LUAD Macrophage 65 at
   40x, i.e. 0.3% detection). These are extreme but not unrealistic: the COAD cohort's
   macrophages have a median of 582 promoter-window counts, a 37x spread.
4. Two deep-type choices and three seeds give 6 replicates per point; the min-max bands are
   wide for the depth-sensitive methods (HCC raw at 40x: 59-79%) mainly because the random
   ordering of the intermediate types changes which real programme sits second. Deep-type-
   specific means are in `depth_spread_summary_by_deep_type.csv`.
5. The held-out panels are the benchmark's; the HCC labels are motif-derived and two clusters
   were relabelled (benchmark rule). The 1x concordance (19/22, 24/36) is the ceiling any method
   can reach on these labels, not 100%.
6. Only the proximal matrices were run; the enhancer matrices are expected to behave the same
   (same normalisers), but this was not tested here. COAD was not used because its natural
   spread is already 37x.

## Files

- `depth_spread_experiment.py` -- the script (`run HCC`, `run LUAD`, `summarise`).
- `depth_spread_results.csv` -- one row per condition x method (370 rows): top share, top type,
  deep_wins, deep_type_share, Spearman vs depth / detection, markers correct, perm p, per lineage.
- `depth_spread_by_type.csv` -- one row per condition x type (444 rows): imposed factor,
  realised median row sum, detection rate, mean of the log1p(CPM) / CPM / raw / z profiles,
  share under each method.
- `depth_spread_markers.csv` -- one row per condition x method x marker.
- `depth_spread_summary.csv`, `depth_spread_summary_by_deep_type.csv` -- means / min / max
  over replicates; `depth_spread_mechanism.csv`, `_summary.csv` -- per-condition ranges and
  correlations of the mechanism variables.
- `fig_depth_spread.png` (600 dpi) / `.pdf` -- A/B top share vs spread, C/D held-out marker
  concordance; lines = mean of 6 replicates, bands = min-max, open markers = natural matrix.
- `*_HCC.csv`, `*_LUAD.csv` -- per-cohort intermediates read by `summarise`.
