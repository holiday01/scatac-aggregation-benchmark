# Controlled coverage-spread experiment

The main benchmark measures the effect of the aggregation strategy at the three
natural coverage spreads of the cohorts (HCC ~2×, NSCLC ~3.5×, COAD ~37×). Here
the spread is *imposed* on two count matrices, so the whole curve is measured
under one design.

> Naming: this directory predates the manuscript's final terminology. In file and
> column names, `depth` = coverage and `pseudobulk_cpm` = pooled CPM.

## Design

|  | HCC count matrix | NSCLC count matrix |
|---|---|---|
| matrix | 12,029 cells × 59,385 genes, 6 types | 53,308 cells × 19,930 genes, 6 types |
| common target coverage (p20 of row sums) | 5,630 | 2,679 |
| natural spread (median row sum, highest / lowest) | 1.84× | 3.57× |
| high-coverage type choices | B_cell (n=416), NK_cytotoxic_T (n=3,677) | B_cell (n=9,484), T_NK (n=12,221) |
| held-out panel | 22 markers, 4 lineages | 36 markers, 6 lineages |

Every cell is binomially thinned to the target, so every type starts at a spread
of 1.00×. The high-coverage type then stays at the target while the other five
are thinned by `S^(-i/5)`, i = 1..5 — evenly log-spaced between 1 and 1/S, their
order a random permutation per seed. S = 1, 2, 5, 10, 20, 40; two high-coverage
choices × three seeds = 6 replicates per spread, plus the natural unthinned
matrix. Five strategies are scored with definitions identical to the benchmark.
The HCC relabel rule of the benchmark is applied.

## Top argmax share vs imposed spread

Mean of 6 replicates; in brackets, replicates in which the high-coverage type wins (/6).

**HCC**

| strategy | 1× | 2× | 5× | 10× | 20× | 40× | natural 1.8× |
|---|---|---|---|---|---|---|---|
| raw mean | 30% (0) | 36% (3) | 48% (6) | 59% (6) | 62% (6) | 68% (6) | 67% (Endothelial_stromal) |
| z-scaled mean | 31% (0) | 35% (3) | 47% (6) | 57% (6) | 61% (6) | 67% (6) | 65% |
| log-normalised mean | 31% (0) | 28% (2) | 32% (6) | 38% (4) | 43% (5) | 50% (4) | 34% (Hepatocyte) |
| pooled CPM | 30% (0) | 29% (0) | 28% (0) | 26% (0) | 24% (0) | 25% (0) | 32% (Hepatocyte) |
| **CPM mean** | 29% (0) | 29% (0) | 27% (0) | 26% (0) | 24% (0) | 24% (0) | 30% (Hepatocyte) |

**NSCLC**

| strategy | 1× | 2× | 5× | 10× | 20× | 40× | natural 3.6× |
|---|---|---|---|---|---|---|---|
| raw mean | 28% (3*) | 69% (6) | 85% (6) | 87% (6) | 92% (6) | 93% (6) | 66% (B_cell) |
| z-scaled mean | 28% (3*) | 68% (6) | 83% (6) | 85% (6) | 90% (6) | 92% (6) | 65% |
| log-normalised mean | 31% (3*) | 54% (6) | 72% (6) | 75% (6) | 81% (6) | 83% (6) | 37% (T_NK) |
| pooled CPM | 25% (3*) | 24% (3*) | 23% (3*) | 22% (1) | 23% (0) | 25% (0) | 24% (B_cell) |
| **CPM mean** | 24% (3*) | 24% (3*) | 22% (3*) | 22% (0) | 23% (0) | 25% (0) | 25% (T_NK) |

\* At 1× (and under the coverage-invariant strategies at 2–5×) the winner on NSCLC
is T_NK under every strategy — the residual, non-coverage structure of the matrix.
T_NK is also one of the two high-coverage choices, hence 3/6 "wins" that have
nothing to do with coverage. Under the CPM mean the winner on HCC is Hepatocyte in
36/36 conditions regardless of which type is high-coverage; on NSCLC it is T_NK at
1–5× and drifts to Macrophage/Epithelial at 20–40×, where the lowest-coverage types
have median row sums of 65–140 counts (caveat 3).

Min–max over the 6 replicates: raw mean at 40× spans 59–79% (HCC) and 89–96%
(NSCLC); CPM mean at 40× spans 22–27% in both.

**Spearman of per-type argmax share vs realised median coverage** (mean over
replicates; NaN at 1×): raw and z-scaled +0.71 (2×) → +0.89–0.96 (5–40×) on both
cohorts; log-normalised +0.54 → +0.89 (HCC), +0.69 → +0.98 (NSCLC); pooled CPM and
CPM mean −0.06 → −0.41 (HCC) and +0.17 → −0.69 (NSCLC) — zero to negative, i.e. the
residual share drifts *away* from the high-coverage type as its competitors lose
their detection floor.

## Held-out marker concordance

Mean of 6 replicates.

| cohort / strategy | 1× | 2× | 5× | 10× | 20× | 40× |
|---|---|---|---|---|---|---|
| HCC (of 22) raw mean | 19.2 | 19.0 | 16.2 | 10.5 | 10.3 | 11.2 |
| HCC z-scaled mean | 19.2 | 19.3 | 16.7 | 10.8 | 10.7 | 11.5 |
| HCC log-normalised mean | 19.2 | 19.7 | 18.8 | 15.2 | 14.7 | 16.0 |
| HCC pooled CPM | 19.2 | 19.5 | 18.7 | 18.7 | 19.0 | 18.3 |
| HCC CPM mean | 19.2 | 19.0 | 18.8 | 18.7 | 18.8 | 18.5 |
| NSCLC (of 36) raw mean | 24.0 | 21.2 | 15.0 | 13.0 | 11.5 | 10.3 |
| NSCLC z-scaled mean | 23.7 | 21.5 | 15.5 | 14.2 | 12.2 | 11.0 |
| NSCLC log-normalised mean | 23.7 | 23.3 | 19.8 | 17.3 | 17.3 | 14.2 |
| NSCLC pooled CPM | 23.8 | 23.5 | 23.3 | 24.2 | 23.2 | 23.3 |
| NSCLC CPM mean | 23.8 | 23.5 | 23.0 | 24.5 | 22.8 | 22.8 |

At 40× the raw mean is at or near the permutation null (perm p up to 0.04 HCC,
0.29 NSCLC; null mean 5.5/22 and 6/36) while the CPM mean stays at p < 0.001 in
every condition. Per lineage, the loss under the raw mean is concentrated in
whichever lineages were thinned: HCC macrophage markers fall from 6/6 to 0.4/6,
NSCLC macrophage 6/6 → 1/6, epithelial 4/6 → 1.5/6, T_NK 6/6 → 3.2/6. Under the
CPM mean the same lineages keep 5.5/6, 4.7/6, 4.7/6 and 6/6. The NSCLC fibroblast
and endothelial panels are poor under *every* strategy already at 1× (0/6 and 2/6)
— a property of the panel or the labels, not of the strategy.

## Mechanism

The per-type detection rate scales with the imposed coverage: its range across
types is 1.07 (1×) → 4.2 (5×) → 8.1 (10×) → 32 (40×) on HCC, 1.3 → 30 on NSCLC
(Spearman detection vs coverage = 1.00 at every spread ≥ 2×). The mean of the
log-normalised profile follows the detection rate (ρ = 0.95–1.00 in every
condition) with a compressed range: 1.04 → 8.6 (HCC) and 1.19 → 10.5 (NSCLC),
because a cell contributes 0 to the mean of a gene it does not detect, however it
is normalised. The mean of the CPM profile is 1.000 across types at every spread
(each cell sums to 1e4), and the mean of the raw profile scales 1:1 with coverage
(range 40.6 at 40×). The argmax share under the log-normalised mean tracks its own
profile mean (ρ 0.54 → 0.91 HCC, 0.88 → 0.98 NSCLC); under the raw mean it tracks
coverage (0.64 → 0.96).

Worked example — HCC, high-coverage type B_cell, seed 0, 40×: the B-cell cluster
keeps a detection rate of 0.070 while the others fall to 0.002–0.038; raw share
75%, log-normalised share 57%, CPM share 14%.

## Interpretation

1. The effect is monotone in the spread and saturates: raw and z-scaled means put
   the argmax in the high-coverage type in 36/36 replicates at ≥ 5×, reaching 68%
   (HCC) and 93% (NSCLC) at 40×. The z-scale is not a fix — it is within 2 points
   of the raw mean at every spread, because per-gene centring and clipping do not
   change which type has the highest mean.
2. The log-normalised mean, the scanpy/Seurat/ArchR convention, halves the effect
   but does not remove it: 50% (HCC) and 83% (NSCLC) at 40×. Its residual coverage
   dependence is the detection-rate term.
3. The CPM mean and pooled CPM are flat at 22–30% across the whole range and
   indistinguishable from each other, the high-coverage type never winning at
   ≥ 10×, with marker concordance unchanged from 1× to 40×. The two behave
   identically here because equal-coverage cells contribute equally to both; they
   differ only when coverage varies *within* a type, which the thinning removes.
4. Cohort sensitivity differs: at 2× the NSCLC raw mean is already at 69% while HCC
   is at 36%. The difference is the amount of genuine lineage structure in the
   matrix — on HCC the Hepatocyte profile wins ~30% of genes at 1× under every
   strategy, a strongly specific programme that a 13% coverage edge does not
   overturn, whereas the NSCLC profiles are nearly tied at equal coverage, so a
   small edge decides most genes.
5. The natural matrices do not sit on the imposed curve at their median spread
   (caveat 1).

## Caveats

1. **The median spread understates what the raw mean sees.** The raw mean responds
   to the per-type *mean* coverage and its tail. Natural HCC: median spread 1.84×
   but mean spread 2.86× (Endothelial_stromal mean 22,752, p90 48,462, vs CD4_T
   7,963). Natural NSCLC: B_cell has the highest *mean* row sum (16,728 vs T_NK
   14,341) although T_NK has the highest median — which is why the natural raw-mean
   argmax goes to B_cell. The imposed design caps every cell at the p20 target, so
   mean = median and the tails are removed. A single scalar "spread" is a coarse
   descriptor; the direction of every conclusion above is unaffected, the exact
   share at a given spread is design-dependent.
2. Thinning is binomial on the gene-activity counts, not on fragments. It is exact
   for a Poisson/binomial sampling model but does not reproduce coverage-dependent
   changes in fragment-to-gene assignment (duplicates, peak calling). Coverage is
   measured here as the matrix row sum — the quantity the main benchmark
   deliberately avoids, but here it is the imposed quantity, so it is the right scale.
3. At 20–40× the lowest-coverage types have median row sums of 65–300 counts
   (NSCLC Macrophage 65 at 40×, i.e. 0.3% detection). Extreme but not unrealistic:
   the COAD cohort's macrophages have a median of 582 gene-window counts at a
   natural 37× spread.
4. Six replicates per point; min–max bands are wide for the coverage-sensitive
   strategies (HCC raw at 40×: 59–79%) mainly because the random ordering of the
   intermediate types changes which real programme sits second.
5. The held-out panels are the benchmark's, and the 1× concordance (19/22, 24/36)
   is the ceiling any strategy can reach on these labels, not 100%.
6. Only the count matrices were run. The weighted matrices are expected to behave
   the same (same normalisers) but this was not tested. COAD was not used because
   its natural spread is already 37×.

## Files

- `depth_spread_experiment.py` — the script (`run HCC`, `run LUAD`, `summarise`).
- `depth_spread_results.csv` — one row per condition × strategy (370 rows).
- `depth_spread_by_type.csv` — one row per condition × type (444 rows).
- `depth_spread_markers.csv` — one row per condition × strategy × marker.
- `depth_spread_summary.csv`, `depth_spread_summary_by_deep_type.csv` — means,
  min and max over replicates.
- `depth_spread_mechanism.csv`, `depth_spread_mechanism_summary.csv` — per-condition
  ranges and correlations of the mechanism variables.
- `*_HCC.csv`, `*_LUAD.csv` — per-cohort intermediates read by `summarise`.
