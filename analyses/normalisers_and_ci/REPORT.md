# Additional normalisers and bootstrap uncertainty for the aggregation benchmark

Journal-version addition to `apbc2026` (2026-08-27). Extends `tools/benchmark_aggregation.py` /
`results/benchmark_aggregation.csv` with (A) ten further per-type profiling methods on the same six
real matrices and (B) cell-level bootstrap intervals and paired tests for the original six methods.
Nothing outside this directory was modified.

## Headline

| method (per-cell transform → mean per type unless stated) | top argmax share, median (range) | argmax in deepest type | held-out markers | depth-safe? |
|---|---|---|---|---|
| raw mean | 73 % (66–98) | 4/6 | 92/198 | no |
| z-scale mean (scanpy `scale`) | 73 % (65–98) | 4/6 | 93/198 | no |
| within-cell rank mean | 66 % (40–99) | 5/6 | 93/198 | no (worse than raw on 4/6) |
| within-cell rank / nnz mean | 43 % (33–98) | 5/6 | 100/198 | no |
| detection-rate mean (binarised) | 52 % (33–100) | 6/6 | 87/198 | no — this *is* the mechanism |
| TF-IDF mean, Signac method 1 | 37 % (28–95) | 5/6 | 107/198 | no (91–95 % on COAD) |
| log1p(CPM) mean (scanpy / Seurat / ArchR) | 36 % (28–83) | 4/6 | 105/198 | no |
| TF-IDF mean, Signac method 3 | = log1p(CPM) mean | | | no (argmax identical) |
| log1p(CPM) mean → per-gene max-norm | = log1p(CPM) mean | | | no (argmax identical) |
| Pearson-residual mean (θ = 100) | 27 % (22–35) | 2/6 | 110/198 | nearly (0–4 points above CPM on 5/6, +14 on COAD enhancer) |
| depth-matched thinning, raw mean | 28 % (17–31) | 3/6 | 107/198 | yes (control) |
| pseudobulk sum → CPM | 25 % (24–32) | 0/6 | 117/198 | yes |
| pseudobulk sum → log1p(CPM) | = pseudobulk CPM | | | yes (argmax identical) |
| **CPM mean (ours)** | **25 % (19–30)** | 3/6 | 106/198 | yes (by construction) |
| median-total scaling mean | = CPM mean | | | yes (argmax identical) |
| TF-IDF mean, Signac method 2 | = CPM mean | | | yes (argmax identical) |

Bootstrap (B = 200, cells resampled within type): every top-share ordering in the tables above is
reproduced in essentially every replicate — the paired difference vs CPM mean excludes zero on all
six matrices for raw, z-scale and log1p(CPM) means (e.g. raw − CPM +37 to +77 points), the top type
is the same in ≥ 90 % of replicates (100 % in 29/36 matrix × method cells), and the marker-count
intervals are ±1–2 genes on HCC/LUAD but 5–13 of 41 on COAD. Paired exact Wilcoxon tests across the
six matrices reach the smallest attainable p (0.031 two-sided) for the seven depth-driven methods vs
CPM on top share; none reaches 0.05 on marker fraction. With n = 6 (three cohorts × two priors that
share cells) these tests only encode a consistent sign pattern.

Figure: `fig_normalisers_ci.png` (A top argmax share, B marker concordance, all 16 methods, percentile
bootstrap intervals for the original six).

## What was run

All computations use the six real gene-activity matrices of `tools/benchmark_aggregation.py`
(HCC proximal/enhancer, COAD proximal/enhancer, LUAD proximal/enhancer; the mis-annotated HCC
negative control is excluded), with the same labels (HCC clusters "B_cell" → Endothelial_stromal and
"DC" → B_cell, applied after loading), the same true per-cell depth (fragment files for HCC/COAD,
peak matrices for LUAD), the same held-out marker panels (22 HCC / 41 COAD / 36 LUAD genes present)
and the same 10,000-shuffle permutation null. `PANEL`, `MATRICES`, `true_depth()`, `perm_p()` and the
original six method functions are imported unchanged from `tools/benchmark_aggregation.py`; nothing
under `apbc2026/` outside this directory or under the data directories was modified.

Matrices are read directly from the h5ad files with h5py into float32 CSR (verified identical to the
`anndata` path on HCC proximal: labels, cell names, genes and every matrix entry). All per-cell
transforms are applied to the CSR data array in chunks of 2×10⁷ entries; dense operations (Pearson
residuals, z-scaling) are done in gene blocks (1,000 and 2,000 genes). The bootstrap shares one CSR
index structure between the five transformed data arrays. Measured resident memory (sum over
concurrently running jobs, sampled every 10 s): see "Resources" below.

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

Cell-level bootstrap, B = 200 per matrix: cells are resampled with replacement **within each
cell type** (type sizes fixed), and all six methods are re-aggregated on the *same* resample, so
differences between methods are paired. Every method is a weighted mean of per-cell rows, so a
resample is a weight vector and re-aggregation is one sparse product per data variant:
raw (→ raw mean and pseudobulk CPM), CPM, log1p(CPM), the binomially thinned matrix, and the
z-scale clip excess. The z-scale mean uses the exact identity
mean_t(min(z, 10)) = (mean_t(x) − μ)/σ − mean_t(max(0, z − 10)), with per-gene μ, σ held at their
full-data values (n ≥ 12,029 cells). The depth-matched bootstrap is conditional on one binomial
thinning; the thinning-only variability is measured separately from 10 independent thinnings.

Per replicate and method: top argmax share over expressed genes, the top type, and the held-out
marker count. Intervals: percentile 2.5–97.5 % of the 200 replicates, plus the basic
(reverse-percentile, 2θ̂ − q) interval, because the top share is a maximum over types and the
percentile interval is biased downward (see Interpretation). Across the six matrices: exact paired
Wilcoxon signed-rank tests (two-sided and one-sided "method > CPM") and the exact sign test on the
point estimates, method vs `cpm_mean`, for top share and marker fraction. With n = 6 pairs the
smallest attainable two-sided Wilcoxon p is 2/2⁶ = 0.031 (one-sided 0.016), so these tests can at most
say "consistent across all six matrices"; they carry no more information than the sign pattern.

## Interpretation

### Which normalisers are depth-safe, which are not, and why

Write the per-type profile of gene *g* in type *t* as the mean over that type's cells of a per-cell
transform *f*(x<sub>ig</sub>, d<sub>i</sub>) (count and depth). Splitting the cells into detected
(x > 0) and undetected,

  mean<sub>t</sub> f = P<sub>t</sub>(detected) × E<sub>t</sub>[f | detected]  (+ f(0), zero for every method here).

P<sub>t</sub>(detected) rises with depth (`results/benchmark_depth_by_type.csv`: COAD proximal Mast 0.178 vs
T/NK 0.009, a 35× depth spread). A method is depth-safe only if the second factor falls with depth fast
enough to cancel the first. This happens exactly in one case: when *f* is linear in the per-cell
fraction x/d. Then P(det) × E[x/d | det] = E[x/d], whose expectation is the gene's fraction of the
cell's fragments and does not depend on d. That is the CPM mean, and it is why median-total scaling
(same thing with a different constant) and Signac's TF-IDF method 2 (TF × a per-gene constant) give the
*identical* per-gene argmax on all six matrices (0–2 of ~20–58k genes differ, float32 ties).
Any transform that is non-linear in x — log1p, rank, binarisation, clipping, log1p(TF × IDF) — breaks
the cancellation and leaves a residual P<sub>t</sub>(detected) term whose size depends on how strongly
the transform compresses the detected values:

- **Detection rate** (binarised mean) is the limiting case, *f* = 1{x > 0}: the profile *is*
  P<sub>t</sub>(detected). 6/6 argmax in the deepest type; 100 % of COAD enhancer genes in Mast;
  ρ(depth, share) = +0.89.
- **Within-cell rank** is worse than the raw count on 4/6 matrices (HCC proximal 74 % vs 67 %;
  COAD 97–99 %): E[rank | detected] scales with the cell's number of non-zeros, i.e. with depth, so
  *both* factors grow with depth. Dividing by nnz (rank / nnz, bounded in (0, 1]) removes the second
  factor but not the first: still 43 % median, 92–98 % on COAD.
- **log1p(CPM)** compresses the detected values (a single fragment in a shallow cell is a *large* CPM,
  so E[log1p(CPM) | detected] actually decreases with depth), which is why it is much better than the
  raw mean (36 % vs 73 %), but the compression is not enough when the depth spread is large:
  60 % / 83 % of COAD genes still peak in Mast. Because ln, log2 and a subsequent per-gene
  max-normalisation are monotone per gene, ArchR style, Signac method 3 and "log then max-norm" all
  have the identical argmax.
- **TF-IDF method 1** (Signac default) is *worse* than log1p(CPM) on the deep-spread matrices
  (COAD 91 % / 95 %; 5/6 argmax in the deepest type): the IDF factor n<sub>cells</sub>/gene-total is
  large for sparse genes, so log1p(TF × IDF × 10⁴) saturates and every detected entry maps to a
  nearly constant value — the profile degenerates towards the detection rate. On HCC/LUAD (≤ 2×
  spread) it is close to log1p(CPM). Its held-out concordance (107/198) is nevertheless as good as
  CPM's, because the marker genes are the well-detected ones.
- **z-scaling** is affine per gene and therefore cannot change the argmax except through the clip at
  10: 73 %, indistinguishable from the raw mean.
- **Pearson residuals** are a different case: the residual (x − μ)/√(μ + μ²/θ) has zero expectation
  in every cell under the null model regardless of depth, so the per-type mean is depth-neutral
  *under the model*. Empirically it is nearly safe — 27 % (22–35 %), the deepest type wins on only 2/6
  matrices, ρ = +0.60, and the best held-out concordance of any per-cell method (110/198; COAD
  proximal 12/41, p = 0.006, vs 10/41 for CPM). It sits above the CPM mean on all six matrices
  (p = 0.031) — by 0–4 points on five of them, but by 14 points on COAD enhancer (35 % vs 21 %), where
  it moves the argmax mass to B cells (a mid-depth type with 3,860 cells), i.e. the residual model's
  variance term and the √n clip leave a type-size/depth footprint that the argmax-share diagnostic
  still catches. It requires dense gene blocks and a full pass over the matrix (5–45 s per matrix here)
  versus a single sparse scaling for CPM.
- **Pseudobulk sum → CPM** is a depth-*weighted* average of per-cell fractions (deep cells weigh more
  within a type, but the per-type CPM removes the between-type depth difference), so it is safe
  (25 %, 0/6 in the deepest type) and has the best held-out total (117/198), mostly from COAD
  (16/41, 14/41, CI 12–18), where the intra-type depth weighting favours the well-covered cells.
  Taking log1p after the sum does not change a single argmax.
- **Depth-matched thinning** (no normaliser) is the control that isolates depth: 28 %, with the
  thinning noise negligible (top-share SD ≤ 0.002 over 10 thinnings).

Practical rule for gene-activity matrices: normalise each cell *linearly* (CPM or any constant), take
the mean, and apply any log, scaling or max-normalisation *after* aggregation, where it is monotone per
gene and cannot move the argmax. Means of log-, rank-, TF-IDF-m1- or binarised values across cells
report P<sub>t</sub>(detected), i.e. depth, whenever depth differs between types; z-scaling before
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
much less affected: raw − CPM, z-scale − CPM and log1p(CPM) − CPM exclude zero on all six matrices;
depth-matched − CPM excludes zero on 5/6 (both signs, |Δ| ≤ 4 points); pseudobulk − CPM on 4/6 (both
signs, |Δ| ≤ 11 points, the largest on COAD enhancer where pseudobulk moves the argmax to B cells).
The top type is unchanged in ≥ 90 % of replicates in every cell of the table (100 % in 29/36), so the
"argmax in deepest type" counts are not sampling artefacts. Interval widths for the top share are
1–4 points for the safe methods and 9–19 points for the raw/z-scale means on HCC and LUAD, where the
deepest types are small (HCC Endothelial_stromal 1,137 cells) and their depth profile moves under
resampling; on COAD (n = 39,398) all widths are ≤ 3 points.

Marker-count intervals are ±1–2 genes on HCC (22-gene panel) and LUAD (36) but 5–13 of 41 on COAD
under every per-cell method (CPM 10/41, CI 6–13, permutation-null mean 5.3): the COAD held-out test is
marginal (perm p 0.03–0.09) for all per-cell methods and clearly positive only for the pseudobulk
methods (16/41, CI 12–17, p < 0.001) and Pearson residuals (12/41, p = 0.006). This should be stated
in the paper wherever COAD concordance is quoted.

### Paired tests: n = 6 caveat

Exact Wilcoxon signed-rank tests on six pairs cannot go below p = 0.031 (two-sided) or 0.016
(one-sided); they are reached whenever all six differences have the same sign, and that is what
happens for raw, z-scale, rank, rank/nnz, detection-rate, TF-IDF m1 and log1p(CPM) vs CPM on the top
share (median Δ from +0.11 to +0.48). Pearson also reaches it, with a median Δ of only +0.023 (but +0.137 on COAD enhancer).
Depth-matched (5/6, p = 0.094) and pseudobulk (4/6, p = 0.16) are not distinguishable from CPM. On the
marker fraction no method differs from CPM at p < 0.05 (raw/z-scale 0/6 better, p = 0.125;
detection-rate 0/6, p = 0.062). Two further caveats: the six matrices are three cohorts × two
gene-activity priors computed on the *same* cells, so the effective number of independent units is
three, and the test is on point estimates, which carry the max-statistic bias discussed above. The
p-values are descriptive; the evidence for the ranking is the consistency of the sign pattern together
with the per-matrix bootstrap differences.

## Resources and provenance

- Hardware: 24 cores, 125 GB host; jobs run with 4 BLAS threads each, ≤ 2–3 concurrently under a
  memory-budgeted scheduler (`run_all.py`). Measured resident memory of all concurrently running
  benchmark jobs, sampled every 10 s (`logs/rss_trace.log`): **peak 24.0 GB** (LUAD enhancer bootstrap
  16.9 GB + LUAD proximal bootstrap ≈ 7 GB). A first launch that read the matrices through
  `anndata.read_h5ad` (float64 + a subsetting copy + a float32 copy) reached 32 GB and was aborted at
  12:57; the loader was replaced by the direct h5py reader before the run reported here (restart
  12:58:24). Wall-clock for all twelve jobs: 22 min; the longest was the LUAD enhancer bootstrap
  (488 M non-zeros, 606 s for 200 replicates).
- Verification: the direct loader reproduces the anndata path exactly (labels, cell names, genes,
  every entry; HCC proximal); the recomputed original methods equal the published top shares and
  marker counts on every matrix (depth-matched differs by ≤ 0.4 points and ≤ 3 marker genes, different random stream); the
  bootstrap point estimates (weights = 1, including the exact clip-excess identity for z-scaling)
  equal the published values; the Pearson block implementation matches
  `scanpy.experimental.pp.normalize_pearson_residuals` to ≤ 5.4 × 10⁻⁶ on every matrix.
- Not done: Seurat/Signac in R were not run (the Signac 1.17.1 `RunTFIDF.default` source was read to
  fix the exact formulas); the bootstrap holds per-gene μ, σ for z-scaling at their full-data values;
  the depth-matched bootstrap is conditional on one thinning.

## Files

| file | content |
|---|---|
| `benchmark_extra_normalisers.py` | Task A; `--matrix i` per matrix, `--merge`; imports `PANEL`, `MATRICES`, `true_depth`, `perm_p` and the original methods from `tools/benchmark_aggregation.py` |
| `benchmark_bootstrap.py` | Task B; same interface; paired Wilcoxon / sign tests in `--merge` |
| `run_all.py` | memory-budgeted scheduler used for the run |
| `make_report_tables.py`, `make_figure.py` | tables below and the figure |
| `benchmark_extra_normalisers.csv` | new methods, same columns as `results/benchmark_aggregation.csv` |
| `benchmark_all_methods.csv` | original six (recomputed) + new, plus `argmax_identical_to` / `closest_original` |
| `benchmark_extra_markers.csv` | per marker × method × matrix |
| `benchmark_bootstrap_ci.csv` | per matrix × method: point, percentile and basic 95 % CIs, bias, top-type stability, marker CI, paired Δ vs CPM with CI, thinning SD |
| `benchmark_paired_tests.csv` | Wilcoxon (two-sided, one-sided) and sign test, every method vs CPM, top share and marker fraction |
| `parts/boot_*.npz` | the 200 × 6 replicate arrays (share, top type, marker count) per matrix |
| `parts/extra_argmax_*.npz` | per-gene argmax vectors of all 16 methods per matrix |
| `fig_normalisers_ci.png/.pdf` | figure |
| `logs/` | per-job logs, scheduler log, RSS trace |

## Tables

### Table 1. All methods, six real matrices (published values for the original six; this run for the new ten)

| method | top argmax share, median (range) | argmax in deepest type | ρ(depth, share) median; p<0.05 | held-out markers | perm p<0.05 | argmax identical to |
|---|---|---|---|---|---|---|
| raw mean | 73% (66%–98%) | 4/6 | +0.89; 5/6 | 92/198 | 6/6 | — |
| z-scale mean (scanpy `scale`, clip 10) | 73% (65%–98%) | 4/6 | +0.86; 5/6 | 93/198 | 6/6 | — |
| within-cell rank mean | 66% (40%–99%) | 5/6 | +0.88; 6/6 | 93/198 | 5/6 | none |
| within-cell rank / nnz mean | 43% (33%–98%) | 5/6 | +0.89; 5/6 | 100/198 | 5/6 | none |
| detection-rate mean (binarised) | 52% (33%–100%) | 6/6 | +0.89; 5/6 | 87/198 | 5/6 | none |
| TF-IDF mean, Signac method 1 (log1p(TF·IDF·1e4)) | 37% (28%–95%) | 5/6 | +0.86; 4/6 | 107/198 | 6/6 | none |
| log1p(CPM) mean (scanpy / Seurat / ArchR) | 36% (28%–83%) | 4/6 | +0.70; 2/6 | 105/198 | 6/6 | — |
| TF-IDF mean, Signac method 3 | 36% (28%–83%) | 4/6 | +0.70; 2/6 | 105/198 | 6/6 | lognorm_mean |
| log1p(CPM) mean, then per-gene max-normalisation | 36% (28%–83%) | 4/6 | +0.70; 2/6 | 105/198 | 6/6 | lognorm_mean |
| Pearson-residual mean (θ = 100) | 27% (22%–35%) | 2/6 | +0.60; 1/6 | 110/198 | 5/6 | none |
| depth-matched thinning, raw mean | 28% (17%–31%) | 3/6 | +0.75; 3/6 | 107/198 | 5/6 | — |
| pseudobulk sum → CPM | 25% (24%–32%) | 0/6 | +0.57; 1/6 | 117/198 | 6/6 | — |
| pseudobulk sum → log1p(CPM) | 25% (24%–32%) | 0/6 | +0.57; 1/6 | 117/198 | 6/6 | pseudobulk_cpm |
| **CPM mean (ours)** | 25% (19%–30%) | 3/6 | +0.50; 1/6 | 106/198 | 5/6 | — |
| median-total scaling mean | 25% (19%–30%) | 3/6 | +0.50; 1/6 | 106/198 | 5/6 | cpm_mean |
| TF-IDF mean, Signac method 2 (TF·log(1+IDF)) | 25% (19%–30%) | 3/6 | +0.50; 1/6 | 106/198 | 5/6 | cpm_mean |

### Table 2. Per-matrix top argmax share (top type) for every method

| method | HCC prox | HCC enha | COAD prox | COAD enha | LUAD prox | LUAD enha |
|---|---|---|---|---|---|---|
| raw mean | 67% Endothelial * | 76% Endothelial * | 95% Mast* | 98% Mast* | 66% B cell | 71% B cell |
| z-scale mean (scanpy `scale`, clip 10) | 65% Endothelial * | 75% Endothelial * | 95% Mast* | 98% Mast* | 65% B cell | 70% B cell |
| within-cell rank mean | 74% Endothelial * | 47% Endothelial * | 97% Mast* | 99% Mast* | 59% B cell | 40% T NK* |
| within-cell rank / nnz mean | 46% Endothelial * | 33% Hepatocyte | 92% Mast* | 98% Mast* | 34% T NK* | 41% T NK* |
| detection-rate mean (binarised) | 55% Endothelial * | 33% Endothelial * | 95% Mast* | 100% Mast* | 41% T NK* | 49% T NK* |
| TF-IDF mean, Signac method 1 (log1p(TF·IDF·1e4)) | 34% Endothelial * | 28% Hepatocyte | 91% Mast* | 95% Mast* | 37% T NK* | 36% T NK* |
| log1p(CPM) mean (scanpy / Seurat / ArchR) | 34% Hepatocyte | 28% Hepatocyte | 60% Mast* | 83% Mast* | 37% T NK* | 36% T NK* |
| TF-IDF mean, Signac method 3 | 34% Hepatocyte | 28% Hepatocyte | 60% Mast* | 83% Mast* | 37% T NK* | 36% T NK* |
| log1p(CPM) mean, then per-gene max-normalisation | 34% Hepatocyte | 28% Hepatocyte | 60% Mast* | 83% Mast* | 37% T NK* | 36% T NK* |
| Pearson-residual mean (θ = 100) | 31% Hepatocyte | 25% Hepatocyte | 22% B cell | 35% B cell | 28% T NK* | 27% T NK* |
| depth-matched thinning, raw mean | 31% Hepatocyte | 29% NK cytotoxic | 17% Macrophage | 23% Mast* | 28% T NK* | 28% T NK* |
| pseudobulk sum → CPM | 32% Hepatocyte | 26% Hepatocyte | 24% Macrophage | 32% B cell | 24% B cell | 25% B cell |
| pseudobulk sum → log1p(CPM) | 32% Hepatocyte | 26% Hepatocyte | 24% Macrophage | 32% B cell | 24% B cell | 25% B cell |
| **CPM mean (ours)** | 30% Hepatocyte | 25% Hepatocyte | 19% Macrophage | 21% Mast* | 25% T NK* | 25% T NK* |
| median-total scaling mean | 30% Hepatocyte | 25% Hepatocyte | 19% Macrophage | 21% Mast* | 25% T NK* | 25% T NK* |
| TF-IDF mean, Signac method 2 (TF·log(1+IDF)) | 30% Hepatocyte | 25% Hepatocyte | 19% Macrophage | 21% Mast* | 25% T NK* | 25% T NK* |

\* = the deepest type (highest median true per-cell depth).

### Table 3. Per-matrix held-out marker concordance (permutation p) for every method

| method | HCC prox | HCC enha | COAD prox | COAD enha | LUAD prox | LUAD enha | total |
|---|---|---|---|---|---|---|---|
| raw mean | 15/22 (p=0.000) | 14/22 (p=0.000) | 9/41 (p=0.010) | 9/41 (p=0.001) | 23/36 (p=0.000) | 22/36 (p=0.000) | 92/198 |
| z-scale mean (scanpy `scale`, clip 10) | 16/22 (p=0.000) | 14/22 (p=0.000) | 9/41 (p=0.009) | 9/41 (p=0.001) | 23/36 (p=0.000) | 22/36 (p=0.000) | 93/198 |
| within-cell rank mean | 15/22 (p=0.000) | 18/22 (p=0.000) | 8/41 (p=0.018) | 5/41 (p=1.000) | 24/36 (p=0.000) | 23/36 (p=0.000) | 93/198 |
| within-cell rank / nnz mean | 18/22 (p=0.000) | 20/22 (p=0.000) | 9/41 (p=0.011) | 5/41 (p=0.874) | 25/36 (p=0.000) | 23/36 (p=0.000) | 100/198 |
| detection-rate mean (binarised) | 18/22 (p=0.000) | 14/22 (p=0.000) | 8/41 (p=0.040) | 5/41 (p=1.000) | 23/36 (p=0.000) | 19/36 (p=0.000) | 87/198 |
| TF-IDF mean, Signac method 1 (log1p(TF·IDF·1e4)) | 20/22 (p=0.000) | 21/22 (p=0.000) | 9/41 (p=0.017) | 10/41 (p=0.000) | 23/36 (p=0.000) | 24/36 (p=0.000) | 107/198 |
| log1p(CPM) mean (scanpy / Seurat / ArchR) | 19/22 (p=0.000) | 21/22 (p=0.000) | 9/41 (p=0.025) | 9/41 (p=0.018) | 23/36 (p=0.000) | 24/36 (p=0.000) | 105/198 |
| TF-IDF mean, Signac method 3 | 19/22 (p=0.000) | 21/22 (p=0.000) | 9/41 (p=0.026) | 9/41 (p=0.022) | 23/36 (p=0.000) | 24/36 (p=0.000) | 105/198 |
| log1p(CPM) mean, then per-gene max-normalisation | 19/22 (p=0.000) | 21/22 (p=0.000) | 9/41 (p=0.028) | 9/41 (p=0.021) | 23/36 (p=0.000) | 24/36 (p=0.000) | 105/198 |
| Pearson-residual mean (θ = 100) | 20/22 (p=0.000) | 21/22 (p=0.000) | 12/41 (p=0.006) | 10/41 (p=0.052) | 23/36 (p=0.000) | 24/36 (p=0.000) | 110/198 |
| depth-matched thinning, raw mean | 19/22 (p=0.000) | 21/22 (p=0.000) | 10/41 (p=0.030) | 9/41 (p=0.063) | 24/36 (p=0.000) | 24/36 (p=0.000) | 107/198 |
| pseudobulk sum → CPM | 20/22 (p=0.000) | 21/22 (p=0.000) | 16/41 (p=0.000) | 14/41 (p=0.001) | 23/36 (p=0.000) | 23/36 (p=0.000) | 117/198 |
| pseudobulk sum → log1p(CPM) | 20/22 (p=0.000) | 21/22 (p=0.000) | 16/41 (p=0.000) | 14/41 (p=0.001) | 23/36 (p=0.000) | 23/36 (p=0.000) | 117/198 |
| **CPM mean (ours)** | 19/22 (p=0.000) | 21/22 (p=0.000) | 10/41 (p=0.032) | 9/41 (p=0.081) | 23/36 (p=0.000) | 24/36 (p=0.000) | 106/198 |
| median-total scaling mean | 19/22 (p=0.000) | 21/22 (p=0.000) | 10/41 (p=0.032) | 9/41 (p=0.087) | 23/36 (p=0.000) | 24/36 (p=0.000) | 106/198 |
| TF-IDF mean, Signac method 2 (TF·log(1+IDF)) | 19/22 (p=0.000) | 21/22 (p=0.000) | 10/41 (p=0.029) | 9/41 (p=0.080) | 23/36 (p=0.000) | 24/36 (p=0.000) | 106/198 |

### Table 4. Cell-level bootstrap (B = 200, cells resampled with replacement within each type) for the original six methods

Top argmax share: point estimate, percentile 95 % interval, basic (bias-corrected) 95 % interval, bootstrap bias; held-out marker count: point estimate and percentile 95 % interval; paired difference vs CPM mean (same resamples).

| matrix | method | top share | percentile 95 % CI | basic 95 % CI | bias | top type (stability) | markers | 95 % CI | Δshare vs CPM [95 % CI] | Δmarkers vs CPM [95 % CI] |
|---|---|---|---|---|---|---|---|---|---|---|
| HCC proximal | raw mean | 67.1% | 55.8%–68.7% | 65.4%–78.3% | -4.3% | Endothelial stromal (100%) | 15/22 | 14–17 | +37.0% [+27.4%, +40.5%] | -4 [-5, -2] |
| HCC proximal | z-scale mean (scanpy `scale`, clip 10) | 65.3% | 54.2%–66.8% | 63.9%–76.5% | -4.0% | Endothelial stromal (100%) | 16/22 | 14–17 | +35.2% [+26.1%, +38.9%] | -3 [-5, -2] |
| HCC proximal | log1p(CPM) mean (scanpy / Seurat / ArchR) | 33.6% | 30.4%–32.4% | 34.8%–36.9% | -2.2% | Hepatocyte (100%) | 19/22 | 19–20 | +3.5% [+2.8%, +4.1%] | +0 [+0, +1] |
| HCC proximal | depth-matched thinning, raw mean | 30.3% | 27.5%–29.2% | 31.5%–33.2% | -2.0% | Hepatocyte (100%) | 19/22 | 18–20 | +0.2% [-0.4%, +0.9%] | +0 [-1, +1] |
| HCC proximal | pseudobulk sum → CPM | 31.6% | 28.1%–30.6% | 32.7%–35.2% | -2.2% | Hepatocyte (100%) | 20/22 | 18–20 | +1.6% [+0.3%, +2.3%] | +1 [-1, +1] |
| HCC proximal | **CPM mean (ours)** | 30.1% | 27.3%–28.8% | 31.4%–32.8% | -2.0% | Hepatocyte (100%) | 19/22 | 19–20 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |
| HCC enhancer | raw mean | 75.6% | 63.0%–81.9% | 69.3%–88.2% | -1.6% | Endothelial stromal (100%) | 14/22 | 13–16 | +50.4% [+38.4%, +57.6%] | -7 [-8, -4] |
| HCC enhancer | z-scale mean (scanpy `scale`, clip 10) | 75.1% | 63.0%–81.0% | 69.1%–87.1% | -1.2% | Endothelial stromal (100%) | 14/22 | 13–16 | +49.8% [+39.1%, +56.6%] | -7 [-8, -5] |
| HCC enhancer | log1p(CPM) mean (scanpy / Seurat / ArchR) | 27.7% | 26.3%–27.9% | 27.6%–29.2% | -0.6% | Hepatocyte (100%) | 21/22 | 20–21 | +2.5% [+2.0%, +3.6%] | +0 [-1, +1] |
| HCC enhancer | depth-matched thinning, raw mean | 29.3% | 26.0%–29.3% | 29.3%–32.6% | -1.6% | NK cytotoxic T (100%) | 21/22 | 20–21 | +4.1% [+1.4%, +5.0%] | +0 [-1, +0] |
| HCC enhancer | pseudobulk sum → CPM | 26.4% | 24.8%–26.4% | 26.4%–28.0% | -0.8% | Hepatocyte (100%) | 21/22 | 19–21 | +1.2% [+0.5%, +1.8%] | +0 [-2, +1] |
| HCC enhancer | **CPM mean (ours)** | 25.2% | 23.8%–25.0% | 25.4%–26.6% | -0.8% | Hepatocyte (100%) | 21/22 | 20–21 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |
| COAD proximal | raw mean | 94.8% | 91.0%–93.2% | 96.3%–98.6% | -2.5% | Mast (100%) | 9/41 | 9–10 | +76.0% [+73.5%, +75.8%] | -1 [-3, +4] |
| COAD proximal | z-scale mean (scanpy `scale`, clip 10) | 94.9% | 91.3%–93.5% | 96.3%–98.5% | -2.3% | Mast (100%) | 9/41 | 9–10 | +76.1% [+74.0%, +76.1%] | -1 [-4, +4] |
| COAD proximal | log1p(CPM) mean (scanpy / Seurat / ArchR) | 60.1% | 55.7%–57.5% | 62.7%–64.4% | -3.5% | Mast (100%) | 9/41 | 8–11 | +41.3% [+38.3%, +40.2%] | -1 [-4, +4] |
| COAD proximal | depth-matched thinning, raw mean | 17.4% | 15.0%–15.8% | 19.0%–19.8% | -2.0% | Macrophage (100%) | 10/41 | 6–13 | -1.4% [-2.5%, -1.7%] | +0 [-4, +5] |
| COAD proximal | pseudobulk sum → CPM | 23.9% | 21.0%–22.3% | 25.5%–26.8% | -2.3% | Macrophage (90%) | 16/41 | 12–17 | +5.1% [+3.5%, +4.9%] | +6 [+2, +10] |
| COAD proximal | **CPM mean (ours)** | 18.8% | 17.1%–17.7% | 19.9%–20.5% | -1.4% | Macrophage (98%) | 10/41 | 6–13 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |
| COAD enhancer | raw mean | 97.6% | 97.2%–97.8% | 97.4%–98.1% | -0.1% | Mast (100%) | 9/41 | 7–10 | +76.7% [+78.5%, +80.2%] | +0 [-3, +4] |
| COAD enhancer | z-scale mean (scanpy `scale`, clip 10) | 97.7% | 97.3%–97.9% | 97.5%–98.1% | -0.1% | Mast (100%) | 9/41 | 7–10 | +76.8% [+78.6%, +80.2%] | +0 [-3, +3] |
| COAD enhancer | log1p(CPM) mean (scanpy / Seurat / ArchR) | 82.8% | 80.7%–82.2% | 83.5%–84.9% | -1.3% | Mast (100%) | 9/41 | 8–10 | +61.9% [+62.2%, +64.3%] | +0 [-2, +4] |
| COAD enhancer | depth-matched thinning, raw mean | 22.8% | 19.2%–20.4% | 25.2%–26.4% | -3.0% | Mast (100%) | 11/41 | 6–13 | +1.9% [+1.0%, +2.2%] | +2 [-2, +5] |
| COAD enhancer | pseudobulk sum → CPM | 31.7% | 28.7%–30.0% | 33.3%–34.6% | -2.3% | B cell (100%) | 14/41 | 12–18 | +10.7% [+10.2%, +12.0%] | +5 [+2, +10] |
| COAD enhancer | **CPM mean (ours)** | 20.9% | 17.5%–18.9% | 23.0%–24.4% | -2.7% | Mast (99%) | 9/41 | 5–12 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |
| LUAD proximal | raw mean | 66.1% | 59.9%–70.4% | 61.8%–72.3% | -0.8% | B cell (100%) | 23/36 | 22–24 | +40.8% [+36.0%, +46.3%] | +0 [-3, +1] |
| LUAD proximal | z-scale mean (scanpy `scale`, clip 10) | 64.6% | 58.4%–69.0% | 60.2%–70.9% | -0.8% | B cell (100%) | 23/36 | 21–24 | +39.4% [+34.3%, +45.0%] | +0 [-3, +1] |
| LUAD proximal | log1p(CPM) mean (scanpy / Seurat / ArchR) | 36.9% | 33.7%–38.1% | 35.7%–40.1% | -1.0% | T NK (100%) | 23/36 | 23–24 | +11.7% [+9.7%, +14.0%] | +0 [-2, +1] |
| LUAD proximal | depth-matched thinning, raw mean | 28.2% | 25.8%–27.6% | 28.8%–30.6% | -1.5% | T NK (100%) | 24/36 | 22–25 | +2.9% [+1.8%, +3.6%] | +1 [-2, +1] |
| LUAD proximal | pseudobulk sum → CPM | 24.1% | 22.9%–23.9% | 24.4%–25.4% | -0.8% | B cell (98%) | 23/36 | 23–25 | -1.1% [-1.3%, +0.0%] | +0 [-2, +1] |
| LUAD proximal | **CPM mean (ours)** | 25.3% | 23.5%–24.5% | 26.0%–27.0% | -1.3% | T NK (100%) | 23/36 | 23–25 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |
| LUAD enhancer | raw mean | 71.3% | 65.8%–74.4% | 68.1%–76.7% | -0.5% | B cell (100%) | 22/36 | 21–23 | +46.1% [+42.2%, +51.0%] | -2 [-4, +1] |
| LUAD enhancer | z-scale mean (scanpy `scale`, clip 10) | 70.4% | 64.7%–73.9% | 66.8%–76.1% | -0.6% | B cell (100%) | 22/36 | 21–23 | +45.2% [+40.9%, +50.3%] | -2 [-4, +1] |
| LUAD enhancer | log1p(CPM) mean (scanpy / Seurat / ArchR) | 35.5% | 32.5%–36.3% | 34.7%–38.5% | -1.0% | T NK (100%) | 24/36 | 23–25 | +10.4% [+9.3%, +12.9%] | +0 [-2, +2] |
| LUAD enhancer | depth-matched thinning, raw mean | 27.7% | 25.3%–27.9% | 27.4%–30.0% | -1.2% | T NK (100%) | 25/36 | 23–26 | +2.6% [+1.8%, +4.5%] | +1 [-2, +3] |
| LUAD enhancer | pseudobulk sum → CPM | 24.6% | 23.3%–24.5% | 24.7%–25.8% | -0.7% | B cell (100%) | 23/36 | 22–25 | -0.5% [-0.3%, +1.3%] | -1 [-3, +1] |
| LUAD enhancer | **CPM mean (ours)** | 25.1% | 22.9%–24.0% | 26.3%–27.3% | -1.7% | T NK (100%) | 24/36 | 22–26 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |

Thinning-only noise (10 independent binomial thinnings, no cell resampling), depth-matched method: top-share SD HCC proximal 0.0009 (markers 19-20), HCC enhancer 0.0016 (markers 21-21), COAD proximal 0.0009 (markers 9-15), COAD enhancer 0.0020 (markers 9-11), LUAD proximal 0.0014 (markers 23-24), LUAD enhancer 0.0015 (markers 23-25).

### Table 5. Paired tests across the six matrices, each method vs CPM mean (point estimates; exact Wilcoxon signed-rank and exact sign test; n = 6 pairs, smallest attainable two-sided Wilcoxon p = 0.031, one-sided 0.016)

| method | metric | matrices with method > CPM | median Δ | Wilcoxon p (two-sided) | Wilcoxon p (method > CPM) | sign-test p |
|---|---|---|---|---|---|---|
| raw mean | top share | 6/6 | +0.483 | 0.031 | 0.016 | 0.031 |
| raw mean | marker fraction | 0/6 | -0.040 | 0.125 | 1.000 | 0.125 |
| z-scale mean (scanpy `scale`, clip 10) | top share | 6/6 | +0.475 | 0.031 | 0.016 | 0.031 |
| z-scale mean (scanpy `scale`, clip 10) | marker fraction | 0/6 | -0.040 | 0.125 | 1.000 | 0.125 |
| within-cell rank mean | top share | 6/6 | +0.386 | 0.031 | 0.016 | 0.031 |
| within-cell rank mean | marker fraction | 1/6 | -0.073 | 0.094 | 0.984 | 0.219 |
| within-cell rank / nnz mean | top share | 6/6 | +0.157 | 0.031 | 0.016 | 0.031 |
| within-cell rank / nnz mean | marker fraction | 1/6 | -0.037 | 0.312 | 0.891 | 0.219 |
| detection-rate mean (binarised) | top share | 6/6 | +0.247 | 0.031 | 0.016 | 0.031 |
| detection-rate mean (binarised) | marker fraction | 0/6 | -0.073 | 0.062 | 1.000 | 0.062 |
| TF-IDF mean, Signac method 1 (log1p(TF·IDF·1e4)) | top share | 6/6 | +0.115 | 0.031 | 0.016 | 0.031 |
| TF-IDF mean, Signac method 1 (log1p(TF·IDF·1e4)) | marker fraction | 2/6 | +0.000 | 0.750 | 0.375 | 1.000 |
| log1p(CPM) mean (scanpy / Seurat / ArchR) | top share | 6/6 | +0.110 | 0.031 | 0.016 | 0.031 |
| log1p(CPM) mean (scanpy / Seurat / ArchR) | marker fraction | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
| TF-IDF mean, Signac method 3 | top share | 6/6 | +0.110 | 0.031 | 0.016 | 0.031 |
| TF-IDF mean, Signac method 3 | marker fraction | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
| log1p(CPM) mean, then per-gene max-normalisation | top share | 6/6 | +0.110 | 0.031 | 0.016 | 0.031 |
| log1p(CPM) mean, then per-gene max-normalisation | marker fraction | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
| Pearson-residual mean (θ = 100) | top share | 6/6 | +0.023 | 0.031 | 0.016 | 0.031 |
| Pearson-residual mean (θ = 100) | marker fraction | 3/6 | +0.012 | 0.250 | 0.125 | 0.250 |
| depth-matched thinning, raw mean | top share | 5/6 | +0.023 | 0.094 | 0.047 | 0.219 |
| depth-matched thinning, raw mean | marker fraction | 1/6 | +0.000 | 1.000 | 0.500 | 1.000 |
| pseudobulk sum → CPM | top share | 4/6 | +0.014 | 0.156 | 0.078 | 0.688 |
| pseudobulk sum → CPM | marker fraction | 3/6 | +0.023 | 0.250 | 0.125 | 0.625 |
| pseudobulk sum → log1p(CPM) | top share | 4/6 | +0.014 | 0.156 | 0.078 | 0.688 |
| pseudobulk sum → log1p(CPM) | marker fraction | 3/6 | +0.023 | 0.250 | 0.125 | 0.625 |
| median-total scaling mean | top share | 1/6 | +0.000 | 1.000 | 0.500 | 1.000 |
| median-total scaling mean | marker fraction | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
| TF-IDF mean, Signac method 2 (TF·log(1+IDF)) | top share | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
| TF-IDF mean, Signac method 2 (TF·log(1+IDF)) | marker fraction | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
