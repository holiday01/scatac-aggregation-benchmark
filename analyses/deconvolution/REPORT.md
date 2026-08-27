# Deconvolution with references built by each aggregation method — REPORT (2026-08-27)

**Question.** If an analyst builds a cell-type reference from a scATAC gene-activity matrix
with each of the aggregation methods benchmarked in `tools/benchmark_aggregation.py`, and then
deconvolves bulk-like mixtures against it with NNLS, how large is the error, and is the deepest
cell type systematically over-estimated with the raw-mean or log-normalised reference?

**Short answer.** The error is set by whether the reference's per-cell weighting matches the
weighting implicit in the bulk signal, not by "raw vs normalised" as such.

* With mixtures built **as specified — the sum of raw counts of the sampled cells — the raw-mean
  reference is the correctly specified linear model for cell fractions** and gives the lowest
  RMSE in all three cohorts (0.039 / 0.088 / 0.046 for HCC / COAD / LUAD), with no deepest-type
  bias (−0.016 / −0.013 / −0.005). Every depth-normalised reference (CPM mean, pseudobulk CPM,
  depth-matched) **over-estimates the deepest type as a cell fraction**: +0.07 (HCC),
  +0.17 to +0.27 (COAD), +0.03 to +0.14 (LUAD). This is not a failure of normalisation: those
  references estimate the **fragment share** of each type (n_t × depth_t), and against that truth
  pseudobulk CPM is nearly exact (RMSE 0.032 / 0.050 / 0.037, r 0.98 / 0.96 / 0.97, deepest bias
  ≤ 0.02) while the raw-mean reference under-estimates the deepest type (−0.10 / −0.30 / −0.03).
* The **log1p(CPM)-mean reference (scanpy / Seurat / ArchR default) is biased under BOTH truths**:
  it over-estimates the deepest type as a cell fraction by +0.07 / +0.08 / +0.15 and is also
  biased on fragment share (COAD −0.21, LUAD +0.12), with r 0.39–0.89 (cell fraction). expm1 of
  the log mean does not repair it (bias +0.06 / +0.06 / +0.07; RMSE 0.059 / 0.157 / 0.103).
* The **z-scaled mean reference is unusable** for deconvolution (r 0.10–0.35; RMSE ≈ 0.17,
  i.e. no better than guessing 1/T; deepest-type bias +0.31 on COAD).
* When the bulk is instead built so that **every cell contributes equally** (each sampled cell
  scaled to 10k before summation), the ranking flips: the CPM mean is best (RMSE 0.042 / 0.064 /
  0.028, deepest bias −0.008 / +0.038 / +0.004) and the raw-mean reference under-estimates the
  deepest type (−0.08 / −0.13 / −0.12).

So the hypothesis in the task ("deepest type over-estimated with raw/log references") holds for
the **log** reference and, in the equal-contribution regime, not for the raw reference; in the
raw-count-sum regime it is the **normalised** references that over-estimate the deepest type as a
cell fraction. The safe statement for the paper is: *the reference must be aggregated with the
same per-cell weighting as the bulk it will be fitted to; the log-mean and z-scaled references are
wrong under every weighting.*

## What was run

Script: `deconvolution_benchmark.py` (this directory); figures by `make_figure.py`.
Inputs: the three proximal gene-activity matrices from `MATRICES` in
`tools/benchmark_aggregation.py` (HCC 12,029 × 59,385; COAD 39,398 × 59,385; LUAD 53,308 × 19,930;
all integer counts). HCC relabel `B_cell→Endothelial_stromal`, `DC→B_cell` applied after loading,
as in the benchmark. The six aggregation functions are **imported from the benchmark script**, so
the reference definitions are identical to the paper's.

1. **Split.** Cells split 50/50 stratified by type; seeds 0, 1, 2 (`split_depth.csv` has per-type
   n and median depth for each half; halves are depth-matched to within a few percent).
2. **References** (genes × types) from the training half: raw mean; z-scale mean
   (`scanpy.pp.scale`, `max_value=10`); log1p(CPM) mean; expm1 of the log1p(CPM) mean (variant);
   pseudobulk CPM; depth-matched raw mean (binomial thinning to the shallowest type's median row
   sum: 7.1–7.2k HCC, 434–439 COAD, 3.0k LUAD); CPM mean.
3. **Signature genes.** *Shared ("common") set*: for each type, the top-200 genes by specificity
   ratio (type mean + c) / (max other-type mean + c), c = 1 % of the mean positive entry, computed
   **once per split from the CPM-mean reference** and used for every method, so that only the
   profile values differ between methods (1,200 genes HCC/LUAD, 1,600 COAD; the top-200 lists are
   disjoint because a gene's ratio exceeds 1 in at most one type; lists in
   `common_signature_genes.csv`). *Per-method ("own") set*: each method selects its own top-200 by
   its own ratio (z-scale: difference, since z-scaled means can be negative). Overlap of the own
   sets with the shared set (`signature_sizes.csv`): raw mean 622 / 219 / 651 of 1,200 / 1,600 /
   1,200 (HCC / COAD / LUAD); pseudobulk 755 / 258 / 680; lognorm 943 / 585 / 852; z-scale
   459 / 136 / 466. On COAD the raw-mean own set has only 1,543 genes and the z-scale set 1,314:
   some types do not have 200 genes in which they are the raw-mean argmax, so their lists overlap
   with other types' lists (consistent with the argmax concentration reported in the paper).
4. **Mixtures** from the held-out half: 100 per seed; proportions ~ Dirichlet(α = 1) over types;
   500 cells drawn multinomially; cells sampled without replacement unless the draw exceeds the
   type's held-out pool (then with replacement; relevant only for HCC B_cell, 208 cells, and COAD
   Endothelial, 359). Primary construction **rawsum**: mixture = sum of raw counts of the sampled
   cells, then CPM over all genes (linear scale). Secondary **cellnorm**: each sampled cell scaled
   to 10k before summation, then CPM. Two truths recorded per mixture: the realised **cell fraction**
   n_t / 500 and the **fragment share** (the fraction of the mixture's total counts contributed by
   type t, from the sampled cells' row sums).
5. **NNLS** (`scipy.optimize.nnls`) of each CPM mixture on each reference restricted to the
   signature genes; fractions normalised to sum 1 (no degenerate all-zero solutions occurred).
   Metrics over all 100 mixtures × T types per seed: RMSE, Pearson r, and per-type bias
   (mean estimated − true); means ± sd over the three seeds in `metrics_summary.csv`; per-type
   bias in `bias_by_type.csv`; every estimate in `estimates.csv` (18 MB).

Deepest type by true depth (from `results/benchmark_depth_by_type.csv`): HCC Endothelial_stromal
(13.4k), COAD Mast (19.2k; 35× the shallowest), LUAD T_NK (10.6k; 3.5×).

## Results — primary setting: rawsum mixtures, shared signature genes (`fig_deconvolution_common_rawsum.png`)

RMSE (mean over 3 seeds; sd ≤ 0.01 except where noted) and bias of the deepest type, truth = **cell fraction**:

| reference | HCC RMSE | COAD RMSE | LUAD RMSE | HCC r | COAD r | LUAD r | deepest-type bias HCC / COAD / LUAD |
|---|---|---|---|---|---|---|---|
| raw mean | **0.039** | **0.088** | **0.046** | 0.96 | 0.72 | 0.95 | −0.016 / −0.013 / −0.005 |
| z-scale mean | 0.174 | 0.165 | 0.168 | 0.18 | 0.35 | 0.10 | −0.001 / **+0.305** / +0.161 |
| log1p(CPM) mean | 0.069 | 0.143 (sd 0.03) | 0.097 | 0.89 | 0.39 | 0.79 | +0.069 / +0.080 / **+0.146** |
| expm1(log mean) | 0.059 | 0.157 | 0.103 | 0.91 | 0.30 | 0.74 | +0.062 / +0.064 / +0.070 |
| pseudobulk CPM | 0.063 | 0.150 | 0.063 | 0.91 | 0.51 | 0.91 | +0.078 / **+0.266** / +0.025 |
| depth-matched | 0.058 | 0.158 | 0.114 | 0.91 | 0.22 | 0.71 | +0.065 / +0.174 / +0.127 |
| CPM mean | 0.067 | 0.165 | 0.125 | 0.89 | 0.20 | 0.65 | +0.074 / +0.224 / +0.141 |

Same mixtures, truth = **fragment share**:

| reference | HCC RMSE | COAD RMSE | LUAD RMSE | HCC r | COAD r | LUAD r | deepest-type bias HCC / COAD / LUAD |
|---|---|---|---|---|---|---|---|
| raw mean | 0.076 | 0.170 | 0.066 | 0.87 | 0.41 | 0.90 | **−0.102 / −0.300 / −0.033** |
| z-scale mean | 0.183 | 0.098 | 0.152 | 0.16 | 0.84 | 0.32 | −0.087 / +0.018 / +0.133 |
| log1p(CPM) mean | 0.047 | 0.167 | 0.082 | 0.95 | 0.48 | 0.85 | −0.017 / −0.207 / +0.117 |
| expm1(log mean) | 0.047 | 0.185 | 0.067 | 0.95 | 0.37 | 0.89 | −0.025 / −0.224 / +0.041 |
| pseudobulk CPM | **0.032** | **0.050** | **0.037** | 0.98 | 0.96 | 0.97 | −0.009 / −0.021 / −0.004 |
| depth-matched | 0.049 | 0.155 | 0.075 | 0.95 | 0.54 | 0.88 | −0.021 / −0.113 / +0.099 |
| CPM mean | 0.056 | 0.150 | 0.086 | 0.93 | 0.58 | 0.84 | −0.012 / −0.063 / +0.112 |

Where the mis-attributed mass goes (per-type bias, `bias_by_type.csv`): on COAD with the CPM-mean
reference, Mast +0.22 and **Fibroblast +0.19** as cell fractions, taken from Macrophage −0.12, T_NK
−0.12, Endothelial −0.10; on LUAD, T_NK +0.14 and B_cell +0.10 from Fibroblast −0.17 and
Macrophage −0.07. The shallow types are the ones that lose mass. On COAD, fragment shares are very
skewed (mean over mixtures: Mast 0.42, Macrophage / T_NK / Fibroblast 0.02–0.03 each, although each
has mean cell fraction 0.11–0.13), so the shallow types are nearly invisible in a raw-count mixture,
which is why COAD is the hardest cohort for every reference.

## Results — per-method signature genes ("own", `fig_deconvolution_own_rawsum.png`)

Letting each method choose its own top-200 changes little (truth = cell fraction; RMSE
HCC / COAD / LUAD): raw mean 0.028 / 0.076 / 0.048; z-scale 0.207 / 0.228 / 0.167; log1p(CPM)
0.056 / 0.152 / 0.068; expm1 0.053 / 0.170 / 0.086; pseudobulk 0.062 / 0.157 / 0.064;
depth-matched 0.058 / 0.131 / 0.107; CPM mean 0.067 / 0.165 / 0.125 (identical to the shared set
by construction). Deepest-type biases keep their sign and size (log +0.07 / +0.04 / +0.07;
pseudobulk +0.08 / +0.30 / +0.03; CPM mean +0.07 / +0.22 / +0.14; raw −0.02 / −0.01 / −0.00).
Truth = fragment share: pseudobulk CPM 0.027 / 0.041 / 0.036. So the gene set is not what
separates the methods here; the profile scaling is.

## Results — equal-contribution mixtures ("cellnorm", `fig_deconvolution_common_cellnorm.png`)

Each sampled cell scaled to 10k before summation; cell fraction and fragment share coincide.
RMSE HCC / COAD / LUAD and deepest-type bias: CPM mean **0.042 / 0.064 / 0.028** (−0.008 / +0.038 /
+0.004); depth-matched 0.035 / 0.065 / 0.043 (−0.014 / +0.007 / −0.022); pseudobulk CPM
0.036 / 0.113 / 0.100 (−0.009 / −0.103 / −0.105); expm1(log mean) 0.033 / 0.099 / 0.068;
log1p(CPM) mean 0.040 / 0.102 / 0.115 (−0.013 / −0.101 / −0.024); raw mean 0.063 / 0.164 / 0.153
(**−0.075 / −0.130 / −0.124**); z-scale 0.176 / 0.144 / 0.166. Here the raw-mean reference
under-estimates the deepest type and the equal-weight references (CPM mean, depth-matched) are best;
pseudobulk CPM now under-estimates the deepest type on COAD/LUAD because it weights training
cells by depth while the bulk does not.

## Interpretation

For a linear mixture m = Σ_t w_t r_t, NNLS recovers w_t up to the scale of r_t. A raw-mean profile
is r_t ≈ (depth_t / 10k) × CPM-mean profile, so its coefficient is w_t ∝ n_t (cell fraction) when
the mixture is a raw-count sum, and w_t ∝ n_t / depth_t when every cell contributes equally.
A depth-normalised profile has coefficient w_t ∝ n_t × depth_t (fragment share) for a raw-count
sum and w_t ∝ n_t for an equal-contribution bulk. The simulation reproduces exactly this: the
"right" reference is the one whose per-cell weighting matches the bulk. Among the depth-normalised
references, pseudobulk CPM (the depth-weighted mean of per-cell CPM) is markedly better than the
CPM mean and the depth-matched mean on raw-count mixtures (COAD fragment-share RMSE 0.050 vs 0.150
and 0.155), because within a type the deep cells dominate a raw-count mixture and pseudobulk
weights the training cells the same way; the equal-weight CPM mean does not, and on COAD (35×
depth spread, and many very shallow macrophage / T / fibroblast cells) that mismatch alone costs
r 0.96 → 0.58.

The log1p(CPM) mean is not a consistent estimator of either weighting: the mean of log values
tracks the per-type detection rate (paper §3.2), so the deepest type's profile is inflated
non-linearly across all genes. This is the reference for which "the deepest type is
over-estimated" holds in the primary setting (+0.07 to +0.15 as a cell fraction), and it remains
biased on fragment share (COAD −0.21, LUAD +0.12), i.e. no re-definition of the truth rescues it;
expm1 of the log mean is not the CPM mean (Jensen) and stays biased. The z-scaled reference has no
linear relation to the mixture at all.

## Caveats (read before citing any number)

1. **The direction of the deepest-type bias depends on the bulk model, and this simulation cannot
   say which model a real bulk ATAC library follows.** The raw-count-sum construction treats the
   35× per-cell depth difference between COAD Mast and Fibroblast as real fragment yield per cell;
   in the source data those differences are at least partly technical (matrix row sums are
   promoter-window counts; types come from different samples/patients at different depths). The
   equal-contribution construction assumes the opposite. Real bulk lies in between (nuclear
   chromatin content differs between types, but not 35×). Report both regimes, or report fragment
   share, which is what a depth-normalised reference estimates.
2. This favours no method a priori but it does mean that **"CPM mean" — the paper's recommended
   aggregation for argmax / marker profiles — is not the best deconvolution reference for a
   raw-count bulk** (COAD r 0.20 on cell fraction, 0.58 on fragment share, vs pseudobulk CPM
   0.51 / 0.96). The paper's claim is about per-gene argmax and marker tests, not about
   deconvolution references; do not extend it.
3. In-sample split: training and held-out halves are the same cells, samples and labels, so
   there are no batch, patient or platform effects; NNLS has no intercept, weights or gene
   filtering by detection. This is an optimistic setting for every method; absolute RMSEs
   will be larger in practice.
4. HCC labels are TF-motif-derived and two clusters are relabelled (`tools/check_hcc_labels.py`);
   HCC B_cell has only 208 held-out cells, so draws above that were with replacement.
5. Three seeds; sds are in `metrics_summary.csv` (RMSE sd ≤ 0.01 in the primary setting except
   log1p(CPM) mean on COAD, 0.03–0.04, and depth-matched / CPM mean deepest-type bias on COAD,
   0.09–0.11, which reflects the binomial thinning to 438 fragments and the shallow-cell noise).
6. Signature genes were chosen from the CPM-mean reference (shared set); the per-method variant
   shows this choice does not drive the ranking. A fully independent gene set (e.g. from RNA) was
   not tested. Dirichlet(1) mixtures have mean fraction 1/T per type; rare-type performance
   (fractions < 0.02) is not specifically assessed.
7. Nothing outside this directory was modified. Runtime 3 min; RSS 3.3 GB observed during the COAD pass (well under the 25 GB cap).

## Files

* `deconvolution_benchmark.py` — the analysis (imports the aggregation functions from `tools/benchmark_aggregation.py`).
* `make_figure.py` — figures.
* `metrics_summary.csv` — mean ± sd over seeds per cohort × mixture construction × gene set × method × truth (RMSE, r, deepest-type bias).
* `metrics_by_seed.csv` — per seed, with per-type bias columns `bias__<type>`.
* `bias_by_type.csv` — per-type bias under both truths, with median true depth and deepest flag.
* `estimates.csv` — every mixture × type estimate with both truths and the Dirichlet parameter.
* `signature_sizes.csv`, `common_signature_genes.csv`, `split_depth.csv`, `run.log`.
* `fig_deconvolution_common_rawsum.png` (primary), `fig_deconvolution_own_rawsum.png`, `fig_deconvolution_common_cellnorm.png`.
