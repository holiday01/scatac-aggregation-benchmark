# Additional normalisers and bootstrap uncertainty

Extends the main benchmark with (A) ten further per-type profiling strategies on
the same six real matrices and (B) cell-level bootstrap intervals and paired tests
for the original six.

> Naming: this directory predates the manuscript's final terminology. Identifiers
> in scripts and CSV columns keep the older names: `depth` = coverage,
> `pseudobulk` = pooled CPM, `proximal` = count matrix, `enhancer` = weighted matrix.

## Headline

| method (per-cell transform → mean per type unless stated) | top argmax share, median (range) | argmax in highest-coverage type | held-out markers | coverage-safe? |
|---|---|---|---|---|
| raw mean | 73 % (66–98) | 4/6 | 92/198 | no |
| z-scaled mean (scanpy `scale`) | 73 % (65–98) | 4/6 | 93/198 | no |
| within-cell rank mean | 66 % (40–99) | 5/6 | 93/198 | no (worse than raw on 4/6) |
| within-cell rank / nnz mean | 43 % (33–98) | 5/6 | 100/198 | no |
| detection-rate mean (binarised) | 52 % (33–100) | 6/6 | 87/198 | no — this *is* the mechanism |
| TF-IDF mean, Signac method 1 | 37 % (28–95) | 5/6 | 107/198 | no (91–95 % on COAD) |
| log-normalised mean (scanpy / Seurat / ArchR) | 36 % (28–83) | 4/6 | 105/198 | no |
| TF-IDF mean, Signac method 3 | = log1p(CPM) mean | | | no (argmax identical) |
| log1p(CPM) mean → per-gene max-norm | = log1p(CPM) mean | | | no (argmax identical) |
| Pearson-residual mean (θ = 100) | 27 % (22–35) | 2/6 | 110/198 | nearly (0–4 points above CPM on 5/6, +14 on COAD weighted) |
| coverage-matched thinning, raw mean | 28 % (17–31) | 3/6 | 107/198 | yes (control) |
| pooled CPM (sum per type → CPM) | 25 % (24–32) | 0/6 | 117/198 | yes |
| pooled sum → log-normalised | = pooled CPM | | | yes (argmax identical) |
| **CPM mean** | **25 % (19–30)** | 3/6 | 106/198 | yes (by construction) |
| median-total scaling mean | = CPM mean | | | yes (argmax identical) |
| TF-IDF mean, Signac method 2 | = CPM mean | | | yes (argmax identical) |

Bootstrap (B = 200, cells resampled within type): every top-share ordering in the tables above is
reproduced in essentially every replicate — the paired difference vs CPM mean excludes zero on all
six matrices for raw, z-scaled and log-normalised means (e.g. raw − CPM +37 to +77 points), the top type
is the same in ≥ 90 % of replicates (100 % in 29/36 matrix × method cells), and the marker-count
intervals are ±1–2 genes on HCC/NSCLC but 5–13 of 41 on COAD. Paired exact Wilcoxon tests across the
six matrices reach the smallest attainable p (0.031 two-sided) for the seven coverage-driven methods vs
CPM on top share; none reaches 0.05 on marker fraction. With n = 6 (three cohorts × two priors that
share cells) these tests only encode a consistent sign pattern.

Figure: `fig_normalisers_ci.png` (A top argmax share, B marker concordance, all 16 methods, percentile
bootstrap intervals for the original six).

## Scope

All computations use the six real gene-activity matrices of the main benchmark
(HCC, COAD and NSCLC × count and weighted; the mislabelled negative control is
excluded), with the same labels, the same true per-cell coverage, the same
held-out marker panels (22 HCC / 41 COAD / 36 NSCLC genes present) and the same
10,000-shuffle permutation null. `PANEL`, `MATRICES`, `true_depth()`, `perm_p()`
and the original six strategy functions are imported unchanged from
`tools/benchmark_aggregation.py`.

### Task A — additional normalisers (`benchmark_extra_normalisers.py`)

| method | definition | note |
|---|---|---|
| `tfidf1_mean` | Signac `RunTFIDF` method 1 (v1.17.1 source checked): log1p(TF × IDF × 10⁴), TF = count / cell total, IDF = n cells / gene total; then mean per type | the default Signac normaliser |
| `tfidf2_mean` | Signac method 2: TF × log(1 + IDF); mean | linear in TF: per gene it is the CPM mean times a constant, so the argmax over types is provably that of `cpm_mean` |
| `tfidf3_mean` | Signac method 3: log1p(TF × 10⁴) × log(1 + IDF); mean | per gene the log1p(CPM) mean times a constant: argmax provably that of `lognorm_mean` |
| `pearson_mean` | analytic Pearson residuals, `scanpy.experimental.pp.normalize_pearson_residuals` formula, θ = 100, clip = √n cells (the default); mean | computed on the **full** matrix in dense blocks of 1,000 genes (exact: residual (i,g) needs only the cell total, gene total and grand total); checked against scanpy on a 2,000 × 1,500 sub-matrix at run time, max |difference| ≤ 6×10⁻⁶ |
| `medscale_mean` | each cell scaled linearly to the median cell total (`scanpy.pp.normalize_total(target_sum=None)`); mean | CPM with a different constant: argmax identical to `cpm_mean` |
| `rank_mean` | each cell's non-zero entries replaced by their within-cell rank (average ties; zeros stay 0); mean | |
| `rankfrac_mean` | as above, ranks divided by the cell's number of non-zeros (top gene of every cell = 1); mean | |
| `detection_mean` | fraction of cells of the type with count > 0 (binarised mean) | the mechanism itself |
| `lognorm_maxnorm` | log1p(CPM) mean, then each gene divided by its maximum over types | per-gene rescaling after aggregation: argmax unchanged by construction |
| `pseudobulk_logcpm` | pseudobulk sum per type → CPM → log1p | log after aggregation: monotone per gene, argmax unchanged |

Read-outs, columns and CSV layout are those of `results/benchmark_aggregation.csv`. The original six
methods were recomputed in the same run only to obtain their per-gene argmax vectors for the identity
checks (their recomputed top shares equal the published ones to all printed digits, except
`depth_matched`, whose binomial thinning uses a different random stream).

### Task B — bootstrap and paired tests (`benchmark_bootstrap.py`)

Cell-level bootstrap, B = 200 per matrix: cells are resampled with replacement
**within each cell type** (type sizes fixed), and all six strategies are
re-aggregated on the *same* resample, so differences between strategies are
paired. The z-scaled mean uses the exact identity
mean_t(min(z, 10)) = (mean_t(x) − μ)/σ − mean_t(max(0, z − 10)), with per-gene
μ and σ held at their full-data values. The coverage-matched bootstrap is
conditional on one binomial thinning; thinning-only variability is measured
separately from 10 independent thinnings.

Per replicate and strategy: top argmax share over expressed genes, the top type,
and the held-out marker count. Intervals are the percentile 2.5–97.5% of the 200
replicates plus the basic (reverse-percentile) interval, because the top share is
a maximum over types and the percentile interval is biased downward. Across the
six matrices: exact paired Wilcoxon signed-rank tests and the exact sign test on
the point estimates, each strategy vs the CPM mean. With n = 6 pairs the smallest
attainable two-sided Wilcoxon p is 2/2⁶ = 0.031, so these tests can at most say
"consistent across all six matrices".

## Interpretation

### Which normalisers are coverage-safe, which are not, and why

Write the per-type profile of gene *g* in type *t* as the mean over that type's cells of a per-cell
transform *f*(x<sub>ig</sub>, d<sub>i</sub>) (count and coverage). Splitting the cells into detected
(x > 0) and undetected,

  mean<sub>t</sub> f = P<sub>t</sub>(detected) × E<sub>t</sub>[f | detected]  (+ f(0), zero for every method here).

P<sub>t</sub>(detected) rises with coverage (`results/benchmark_depth_by_type.csv`: COAD count Mast 0.178 vs
T/NK 0.009, a 35× coverage spread). A method is coverage-safe only if the second factor falls with coverage fast
enough to cancel the first. This happens exactly in one case: when *f* is linear in the per-cell
fraction x/d. Then P(det) × E[x/d | det] = E[x/d], whose expectation is the gene's fraction of the
cell's fragments and does not depend on d. That is the CPM mean, and it is why median-total scaling
(same thing with a different constant) and Signac's TF-IDF method 2 (TF × a per-gene constant) give the
*identical* per-gene argmax on all six matrices (0–2 of ~20–58k genes differ, float32 ties).
Any transform that is non-linear in x — log1p, rank, binarisation, clipping, log1p(TF × IDF) — breaks
the cancellation and leaves a residual P<sub>t</sub>(detected) term whose size depends on how strongly
the transform compresses the detected values:

- **Detection rate** (binarised mean) is the limiting case, *f* = 1{x > 0}: the profile *is*
  P<sub>t</sub>(detected). 6/6 argmax in the highest-coverage type; 100 % of COAD weighted genes in Mast;
  ρ(coverage, share) = +0.89.
- **Within-cell rank** is worse than the raw count on 4/6 matrices (HCC count 74 % vs 67 %;
  COAD 97–99 %): E[rank | detected] scales with the cell's number of non-zeros, i.e. with coverage, so
  *both* factors grow with coverage. Dividing by nnz (rank / nnz, bounded in (0, 1]) removes the second
  factor but not the first: still 43 % median, 92–98 % on COAD.
- **log1p(CPM)** compresses the detected values (a single fragment in a shallow cell is a *large* CPM,
  so E[log1p(CPM) | detected] actually decreases with coverage), which is why it is much better than the
  raw mean (36 % vs 73 %), but the compression is not enough when the coverage spread is large:
  60 % / 83 % of COAD genes still peak in Mast. Because ln, log2 and a subsequent per-gene
  max-normalisation are monotone per gene, ArchR style, Signac method 3 and "log then max-norm" all
  have the identical argmax.
- **TF-IDF method 1** (Signac default) is *worse* than log1p(CPM) on the deep-spread matrices
  (COAD 91 % / 95 %; 5/6 argmax in the highest-coverage type): the IDF factor n<sub>cells</sub>/gene-total is
  large for sparse genes, so log1p(TF × IDF × 10⁴) saturates and every detected entry maps to a
  nearly constant value — the profile degenerates towards the detection rate. On HCC/NSCLC (≤ 2×
  spread) it is close to log1p(CPM). Its held-out concordance (107/198) is nevertheless as good as
  CPM's, because the marker genes are the well-detected ones.
- **z-scaling** is affine per gene and therefore cannot change the argmax except through the clip at
  10: 73 %, indistinguishable from the raw mean.
- **Pearson residuals** are a different case: the residual (x − μ)/√(μ + μ²/θ) has zero expectation
  in every cell under the null model regardless of coverage, so the per-type mean is coverage-neutral
  *under the model*. Empirically it is nearly safe — 27 % (22–35 %), the highest-coverage type wins on only 2/6
  matrices, ρ = +0.60, and the best held-out concordance of any per-cell method (110/198; COAD
  count 12/41, p = 0.006, vs 10/41 for CPM). It sits above the CPM mean on all six matrices
  (p = 0.031) — by 0–4 points on five of them, but by 14 points on COAD weighted (35 % vs 21 %), where
  it moves the argmax mass to B cells (a mid-coverage type with 3,860 cells), i.e. the residual model's
  variance term and the √n clip leave a type-size/coverage footprint that the argmax-share diagnostic
  still catches. It requires dense gene blocks and a full pass over the matrix (5–45 s per matrix here)
  versus a single sparse scaling for CPM.
- **Pooled CPM** is a coverage-*weighted* average of per-cell fractions (deep cells weigh more
  within a type, but the per-type CPM removes the between-type coverage difference), so it is safe
  (25 %, 0/6 in the highest-coverage type) and has the best held-out total (117/198), mostly from COAD
  (16/41, 14/41, CI 12–18), where the intra-type coverage weighting favours the well-covered cells.
  Taking log1p after the sum does not change a single argmax.
- **Coverage-matched thinning** (no normaliser) is the control that isolates coverage: 28 %, with the
  thinning noise negligible (top-share SD ≤ 0.002 over 10 thinnings).

Practical rule for gene-activity matrices: normalise each cell *linearly* (CPM or any constant), take
the mean, and apply any log, scaling or max-normalisation *after* aggregation, where it is monotone per
gene and cannot move the argmax. Means of log-, rank-, TF-IDF-m1- or binarised values across cells
report P<sub>t</sub>(detected), i.e. coverage, whenever coverage differs between types; z-scaling before
averaging is equivalent to not normalising at all. Pearson residuals are an acceptable alternative
if their cost is affordable, with a small residual bias that the argmax-share diagnostic detects.

### Bootstrap intervals: what they say and their bias

The percentile intervals of the top share lie 0.1–4.3 points *below* the point estimate on every
matrix × method (column "bias" in Table 4). This is a property of the statistic, not a bug: the top
share is a maximum over types, and re-sampling cells adds a second layer of sampling noise, so genes
whose type means are near-tied flip argmax at random and the share regresses towards 1/n<sub>types</sub>.
(A statistic with the opposite curvature would be biased upward.) Table 4 therefore also gives the
basic (reverse-percentile) interval, which corrects first-order bias; the percentile interval is what
is plotted. The bias is shared across methods on the same resample, so the paired differences are
much less affected: raw − CPM, z-scaled − CPM and log1p(CPM) − CPM exclude zero on all six matrices;
coverage-matched − CPM excludes zero on 5/6 (both signs, |Δ| ≤ 4 points); pooled CPM − CPM mean on 4/6 (both
signs, |Δ| ≤ 11 points, the largest on COAD weighted where pooled CPM moves the argmax to B cells).
The top type is unchanged in ≥ 90 % of replicates in every cell of the table (100 % in 29/36), so the
"argmax in highest-coverage type" counts are not sampling artefacts. Interval widths for the top share are
1–4 points for the safe methods and 9–19 points for the raw/z-scaled means on HCC and NSCLC, where the
highest-coverage types are small (HCC Endothelial_stromal 1,137 cells) and their coverage profile moves under
resampling; on COAD (n = 39,398) all widths are ≤ 3 points.

Marker-count intervals are ±1–2 genes on HCC (22-gene panel) and NSCLC (36) but 5–13 of 41 on COAD
under every per-cell method (CPM 10/41, CI 6–13, permutation-null mean 5.3): the COAD held-out test is
marginal (perm p 0.03–0.09) for all per-cell methods and clearly positive only for the pseudobulk
methods (16/41, CI 12–17, p < 0.001) and Pearson residuals (12/41, p = 0.006). This should be stated
in the paper wherever COAD concordance is quoted.

### Paired tests: n = 6 caveat

Exact Wilcoxon signed-rank tests on six pairs cannot go below p = 0.031 (two-sided) or 0.016
(one-sided); they are reached whenever all six differences have the same sign, and that is what
happens for raw, z-scaled, rank, rank/nnz, detection-rate, TF-IDF m1 and log1p(CPM) vs CPM on the top
share (median Δ from +0.11 to +0.48). Pearson also reaches it, with a median Δ of only +0.023 (but +0.137 on COAD weighted).
Coverage-matched (5/6, p = 0.094) and pooled CPM (4/6, p = 0.16) are not distinguishable from CPM. On the
marker fraction no method differs from CPM at p < 0.05 (raw/z-scale 0/6 better, p = 0.125;
detection-rate 0/6, p = 0.062). Two further caveats: the six matrices are three cohorts × two
gene-activity priors computed on the *same* cells, so the effective number of independent units is
three, and the test is on point estimates, which carry the max-statistic bias discussed above. The
p-values are descriptive; the evidence for the ranking is the consistency of the sign pattern together
with the per-matrix bootstrap differences.

## Files

| file | content |
|---|---|
| `benchmark_extra_normalisers.py` | Task A; `--matrix i` per matrix, `--merge`; imports `PANEL`, `MATRICES`, `true_depth`, `perm_p` and the original methods from `tools/benchmark_aggregation.py` |
| `benchmark_bootstrap.py` | Task B; same interface; paired Wilcoxon / sign tests in `--merge` |
| `make_report_tables.py`, `make_figure.py` | tables below and the figure |
| `benchmark_extra_normalisers.csv` | new methods, same columns as `results/benchmark_aggregation.csv` |
| `benchmark_all_methods.csv` | original six (recomputed) + new, plus `argmax_identical_to` / `closest_original` |
| `benchmark_extra_markers.csv` | per marker × method × matrix |
| `benchmark_bootstrap_ci.csv` | per matrix × method: point, percentile and basic 95 % CIs, bias, top-type stability, marker CI, paired Δ vs CPM with CI, thinning SD |
| `benchmark_paired_tests.csv` | Wilcoxon (two-sided, one-sided) and sign test, every method vs CPM, top share and marker fraction |

## Tables

The five result tables are in [`report_tables.md`](report_tables.md):
all strategies on the six matrices; per-matrix top argmax share; per-matrix
held-out marker concordance; the cell-level bootstrap; and the paired tests.
