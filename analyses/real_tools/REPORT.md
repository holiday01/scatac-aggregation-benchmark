# Real-tool confirmation of the benchmark (HCC count matrix)

The benchmark's aggregation strategies are Python re-implementations of what
Seurat, Signac, scanpy, ArchR and SnapATAC2 do. What does an analyst who runs the
*real* tools with their documented defaults obtain on the same matrix, scored with
the same two read-outs — top argmax share over expressed genes, and held-out
22-marker concordance against a 10,000-shuffle permutation null?

> Naming: this directory predates the manuscript's final terminology. In file and
> column names, `depth` = coverage, `proximal` = the count matrix, and
> `pseudobulk_cpm` = pooled CPM.

**Input.** `atac_hcc_geneactivity_proximal.h5ad` (12,029 cells × 59,385 genes,
integer counts) with the benchmark's relabel rule; 6 types (B_cell 416, CD4_T
3,305, Endothelial_stromal 1,137, Hepatocyte 1,415, Macrophage 2,079,
NK_cytotoxic_T 3,677). True coverage = `obs['n_fragment']` from the SnapATAC2
fragment import (median per type: Endothelial_stromal 13,435; Hepatocyte 12,522;
Macrophage 8,737; B_cell 8,040; NK_cytotoxic_T 7,835; CD4_T 6,643).

**Software run.** R 4.5.2, Seurat 5.5.0, SeuratObject 5.4.0, Signac 1.17.1,
Rsamtools 2.26.0, GenomicRanges 1.62.1, rtracklayer 1.70.1. Python 3.11: scanpy
1.11.5, anndata 0.12.10, SnapATAC2 2.9.0, pysam 0.24.0. EnsDb.Hsapiens.v86,
ensembldb, BSgenome, presto and ArchR are **not** installed; ArchR was not run.

## 1. Headline table (every tool path, same read-outs as the paper)

| tool | call | what it computes | top type | top share | held-out markers | paper method it reproduces | paper (top share, markers) |
|---|---|---|---|---|---|---|---|
| Seurat | `NormalizeData()` -> `AverageExpression(group.by='cell_type')` (defaults) | expm1 of the LogNormalize `data` layer, then mean per type = CP10k mean | Hepatocyte | **30.1%** | **19/22** (p<1e-4) | cpm_mean | 30.1%, 19/22 |
| Seurat | `AverageExpression(..., return.seurat=TRUE)` `data` layer | log1p(mean(expm1(data))) — same argmax as above | Hepatocyte | 30.1% | 19/22 | cpm_mean | 30.1%, 19/22 |
| Seurat | `AverageExpression(..., layer='counts')` | mean of raw counts per type | Endothelial_stromal (highest coverage) | **67.1%** | **15/22** | raw_mean | 67.1%, 15/22 |
| Seurat | `AverageExpression(..., layer=<copy of data under another name>)`; also `rowMeans(LayerData(layer='data'))` per type | mean of log1p(CP10k) values (no expm1, because the expm1 is triggered by the layer *name* 'data') | Hepatocyte | **33.6%** | 19/22 | lognorm_mean | 33.6%, 19/22 |
| Seurat | `AggregateExpression(group.by='cell_type')` (default) | sum of counts per type (pooled sum, un-normalised) | NK_cytotoxic_T (largest n) | **41.0%** | **17/22** | (sum: no paper method; tracks n cells x depth) | — |
| Seurat | `AggregateExpression(..., return.seurat=TRUE)` `data` layer | LogNormalize(sums, 1e4) = log pooled CPM | Hepatocyte | 31.6% | 20/22 | pooled CPM | 31.6%, 20/22 |
| Signac | `GeneActivity()` (defaults) -> `NormalizeData(LogNormalize, scale.factor=median(nCount_RNA))` -> `AverageExpression()` (the vignette) | CP-median mean of Signac's own 19,622-gene matrix | Hepatocyte | **25.7%** | **18/21** | cpm_mean | 30.1%, 19/22 |
| Signac | same with `scale.factor=1e4` | identical argmax (global scale factor) | Hepatocyte | 25.7% | 18/21 | cpm_mean | 30.1%, 19/22 |
| Signac | `GeneActivity()` counts, mean per type | raw mean of Signac's matrix | Endothelial_stromal (highest coverage) | **72.4%** | **15/21** | raw_mean | 67.1%, 15/22 |
| Signac | vignette normalisation, mean of the `data` layer | mean of log values | Hepatocyte | 27.8% (27.9% at 1e4) | 18/21 | lognorm_mean | 33.6%, 19/22 |
| SnapATAC2 | `snap.pp.make_gene_matrix(data, gene_anno=GTF)` (all defaults), raw mean | mean of the tool's default counts per type | Endothelial_stromal (highest coverage) | **67.1%** | **15/22** | raw_mean | 67.1%, 15/22 |
| SnapATAC2 | tutorial: `sc.pp.normalize_total()` (target = median) + `sc.pp.log1p()`, mean | mean of log1p(CP-median) | Hepatocyte | **33.4%** | 19/22 | lognorm_mean | 33.6%, 19/22 |
| SnapATAC2 | `sc.pp.normalize_total(target_sum=1e4)` + `log1p`, mean | mean of log1p(CP10k) | Hepatocyte | 33.6% | 19/22 | lognorm_mean | 33.6%, 19/22 |
| SnapATAC2 | `sc.pp.normalize_total(1e4)` (no log), mean | CP10k mean | Hepatocyte | **30.1%** | 19/22 | cpm_mean | 30.1%, 19/22 |
| SnapATAC2 | sum per type, CPM | pooled CPM | Hepatocyte | 31.6% | 20/22 | pooled CPM | 31.6%, 20/22 |
| scanpy | `sc.pp.normalize_total(1e4)`+`log1p` -> `sc.get.aggregate(by, func='mean')` | mean of log1p(CP10k) | Hepatocyte | **33.6%** | 19/22 | lognorm_mean | 33.6%, 19/22 |
| scanpy | `sc.get.obs_df(...).groupby('cell_type').mean()` | same | Hepatocyte | 33.6% | 19/22 | lognorm_mean | 33.6%, 19/22 |
| scanpy | `sc.pl.DotPlot(...).dot_color_df` (22 markers) | same (dot colour); dot size = fraction expressing | — | — | colour 19/22, size 18/22 | lognorm_mean | 19/22 |

Permutation p for every path above is < 1e-4 (10,000 shuffles; null mean ~4.2/22) — as in the paper,
the marker test does not discriminate the normalisers on the correctly named HCC matrix; the argmax
share and the identity of the misassigned markers do.

Which markers are wrong (identical across tools for a given aggregation):
- every per-cell-normalised mean (Seurat AverageExpression default, Signac vignette, SnapATAC2 tutorial,
  scanpy aggregate/obs_df/dotplot): CD19 -> Macrophage, IL7R -> CD4_T, APOA1 -> Endothelial_stromal (19/22);
- every raw-count mean (Seurat layer='counts', Signac counts, SnapATAC2 default output): the macrophage
  panel collapses (MRC1, MSR1 -> Endothelial_stromal; CD163 -> Hepatocyte; ITGAM -> B_cell) plus BANK1 ->
  Endothelial_stromal, IL7R, APOA1 (15/22; macrophage 2/6);
- Seurat `AggregateExpression()` sums: CD79A, CD79B -> NK_cytotoxic_T, CD19 -> Macrophage, BANK1 ->
  Endothelial_stromal, APOA1 -> NK_cytotoxic_T (17/22; B cell 1/5, the smallest population loses its
  markers to the largest);
- pooled CPM (Seurat AggregateExpression return.seurat, SnapATAC2 sum/CPM): IL7R, APOA1 only (20/22).

## 2. Seurat (`02_seurat_hcc.R`)

Object: `ReadMtx` -> `CreateSeuratObject(counts, meta.data)` with defaults (`min.cells = 0`,
`min.features = 0`; all 12,029 cells and 59,385 features kept; Seurat rewrote the 27 underscore-containing
feature names with dashes — handled by position). `NormalizeData()` defaults = `LogNormalize`,
`scale.factor = 1e4`.

Verified numerically inside R (max abs difference against explicit `Matrix::rowMeans`):
- `AverageExpression(group.by='cell_type')` [default `layer = 'data'`] = mean(expm1(data)) to 2.9e-14, and
  differs from mean(data) by up to 4.78. **So Seurat's documented default is the paper's CPM mean
  (cpm_mean), not the log-mean.** Seurat prints "As of Seurat v5, we recommend using AggregateExpression to
  perform pseudo-bulk analysis" — the documented alternative is pooled sums.
- `AverageExpression(layer='counts')` = mean of raw counts to 2e-13 (raw_mean).
- There is no documented argument for "mean of the log values"; the expm1 is applied whenever the layer is
  called `data`. Copying the data layer under another name (`obj[["RNA"]]$logdata <- data`) and averaging
  that layer gives mean(data) to 7.5e-15, i.e. the paper's lognorm_mean. This is also what
  `DotPlot`/`DoHeatmap`-style summaries and any `rowMeans` on the data slot compute.
- `AverageExpression(return.seurat=TRUE)` stores log1p(mean(expm1(data))) (to 4.7e-15): a monotone
  transform of the CPM mean, same argmax.
- `AggregateExpression()` = per-type column sums exactly; with `return.seurat=TRUE` the data layer is
  `LogNormalize(sums, 1e4)` = log pooled CPM.

Read-outs (scored by `07_score_tool_profiles.py`, table `results/seurat_signac_profile_results.csv`):
see Section 1. Agreement with the paper's re-implementation is exact for every path (30.1%/19,
67.1%/15, 33.6%/19, 31.6%/20).

### FindAllMarkers
Full-size run: `FindAllMarkers(obj, only.pos = TRUE)` on all 12,029 cells × 59,385
genes with Seurat's default test (`test.use = 'wilcox'`, which without the optional
`presto` package falls back to base-R `wilcox.test` gene by gene;
`logfc.threshold = 0.1`, `min.pct = 0.01`, `return.thresh = 0.01`) returned 82,336
rows. Per-type counts are in `results/seurat_findallmarkers_summary.csv`; the
full per-gene tables are large and are not committed — re-run
`11_seurat_findallmarkers_reduced.R` to regenerate them.

Data layer (`slot = 'data'`, the default), Wilcoxon on log1p(CP10k), all cells:

| type | n cells | median true coverage | genes returned | padj < 0.05 | padj < 0.05 & log2FC > 0.5 | own canonical markers in list | in top 50 | ranks of own markers |
|---|---|---|---|---|---|---|---|---|
| Endothelial_stromal | 1,137 | 13,435 | 25,285 | 18,244 | 12,142 | 0/0 | 0 | — |
| Hepatocyte | 1,415 | 12,522 | 23,356 | 16,053 | 13,076 | 4/5 | 2 | ALB:442; APOB:41; HNF4A:33; TTR:1614; APOA1:absent |
| Macrophage | 2,079 | 8,737 | 8,345 | 5,133 | 4,098 | 6/6 | 1 | CD68:170; CSF1R:8; MRC1:354; MSR1:140; CD163:1340; ITGAM:960 |
| B_cell | 416 | 8,040 | 8,187 | 2,532 | 2,247 | 5/5 | 1 | MS4A1:32; CD79A:228; CD79B:631; CD19:1269; BANK1:273 |
| NK_cytotoxic_T | 3,677 | 7,835 | 12,124 | 9,416 | 3,504 | 6/6 | 0 | CD3D:122; CD3E:571; CD8A:192; IL7R:1608; NKG7:210; GNLY:1388 |
| CD4_T | 3,305 | 6,643 | 5,039 | 1,074 | 531 | 0/0 | 0 | — |

Spearman(median true coverage, genes returned) = +0.83 (p = 0.042); with padj < 0.05:
+0.83 (p = 0.042); with padj < 0.05 and log2FC > 0.5: +0.89 (p = 0.019).
The two highest-coverage types receive 18,244 and 16,053 significant positive markers, the shallowest
(CD4_T, which has 3,305 cells — more than either) 1,074: a 17-fold spread that follows coverage, not
population size (B_cell, n = 416, still gets 2,532). The scanpy `rank_genes_groups` run on the same
matrix (`results/benchmark_wilcoxon.csv`) gave rho 0.77.

Canonical markers: 21 of the 22 held-out markers appear in their own type's positive-marker list, all
21 with padj < 0.05; the missing one is APOA1 (Hepatocyte). But only **4/22 are in their type's top 50**
(the order Seurat returns, by p): MS4A1 (B), APOB and HNF4A (Hepatocyte), CSF1R (Macrophage); the rest sit at
ranks 122–1,614 (CD163 1,340, GNLY 1,388, CD19 1,269, IL7R 1,608, TTR 1,614). The top-10 lists are
coverage/accessibility-driven rather than lineage-driven — Endothelial_stromal: COL4A1,ADAMTS9,PCDH17,ZNF521,FAT4,RP11-613D13.4,KIRREL1,ARHGAP29,IGFBP7,ADGRF5;
Hepatocyte: MET,SYBU,EGFR,AGT,SLC30A10,HNF4G,KHDRBS3,RP11-99J16__A.2,SDC2,PROX1-AS1; Macrophage: RAB31,DAPK1,RBM47,FPR1,EPB41L3,FMNL2,FPR3,CSF1R,MS4A6E,SDS; B_cell: RP11-534L6.2,MIR339,RP11-1134I14.6,AC104024.1,U62631.5,DTNB-AS1,LL22NC03-30E12.13,LA16c-329F2.1,CH17-84K15.2,RP11-136K14.2 (a
list of lncRNA/clone-named loci); NK_cytotoxic_T: BCL11B,FYN,IKZF1,PTPRC,PRKCH,CD247,RUNX3,RP11-553L6.2,SCML4,PRKCQ. One foreign marker enters a
top-50 list (CD4_T). CD4_T and Endothelial_stromal have no held-out markers in the panel.

Reduced problem (equal cell numbers): to separate coverage from population size, and to obtain the
`slot = 'counts'` comparison within the time budget, the same calls were repeated on at most 500 cells per
type (seed 1; B_cell keeps all 416 -> 2,916 cells) and the 19,730 expressed protein-coding genes
(GENCODE v38 `gene_type`) (`10_subsample_export.py`, `11_seurat_findallmarkers_reduced.R`).

Reduced, data layer:

| type | n cells used | median true coverage | genes returned | padj < 0.05 | padj < 0.05 & log2FC > 0.5 | own canonical markers in list | in top 50 | ranks of own markers |
|---|---|---|---|---|---|---|---|---|
| Endothelial_stromal | 500 | 13,435 | 6,401 | 4,634 | 2,644 | 0/0 | 0 | — |
| Hepatocyte | 500 | 12,522 | 5,853 | 3,701 | 3,196 | 4/5 | 2 | ALB:405; APOB:21; HNF4A:34; TTR:967; APOA1:absent |
| Macrophage | 500 | 8,737 | 1,946 | 918 | 813 | 6/6 | 1 | CD68:54; CSF1R:9; MRC1:253; MSR1:173; CD163:813; ITGAM:415 |
| B_cell | 416 | 8,040 | 3,883 | 1,280 | 1,121 | 5/5 | 2 | MS4A1:22; CD79A:41; CD79B:321; CD19:883; BANK1:146 |
| NK_cytotoxic_T | 500 | 7,835 | 4,460 | 1,760 | 1,370 | 6/6 | 1 | CD3D:40; CD3E:230; CD8A:1329; IL7R:100; NKG7:123; GNLY:1403 |
| CD4_T | 500 | 6,643 | 978 | 274 | 238 | 0/0 | 0 | — |

Spearman(median true coverage, genes returned) = +0.77 (p = 0.072);
with padj < 0.05: +0.77 (p = 0.072);
with padj < 0.05 and log2FC > 0.5: +0.71 (p = 0.111).
The two highest-coverage types (Endothelial_stromal, Hepatocyte) receive 4,634 and 3,701
significant positive markers, the shallowest (CD4_T) 274 — a 17-fold spread at equal
cell numbers, ordered by coverage (the one exception to the ordering is Macrophage, whose 918 sits below
NK_cytotoxic_T and B_cell despite slightly higher coverage) — the full-size result is reproduced with
population size removed.

Canonical markers: 21 of the 22 held-out markers appear in their own type's positive-marker list
(21 of them with padj < 0.05); the missing one is APOA1 (Hepatocyte). Only 6/22 are in
their type's top 50 (ordered by p, as Seurat returns them): MS4A1, CD79A (B); APOB, HNF4A (Hepatocyte);
CSF1R (Macrophage); CD3D (NK/T). The rest sit at ranks ~100–1,400 (CD163 813, CD8A 1,329, GNLY 1,403),
i.e. an analyst reading the top of each list would not see them; 1 foreign marker enters a top-50
list (CD4_T). CD4_T and Endothelial_stromal have no held-out markers in the panel.

Reduced, raw counts (`slot = 'counts'`; Seurat accepts it and switches the fold-change to log2 of mean counts):

| type | n cells used | median true coverage | genes returned | padj < 0.05 | padj < 0.05 & log2FC > 0.5 | own canonical markers in list | in top 50 | ranks of own markers |
|---|---|---|---|---|---|---|---|---|
| Endothelial_stromal | 500 | 13,435 | 6,490 | 5,550 | 2,668 | 0/0 | 0 | — |
| Hepatocyte | 500 | 12,522 | 6,022 | 3,871 | 3,201 | 4/5 | 2 | ALB:411; APOB:24; HNF4A:33; TTR:1040; APOA1:absent |
| Macrophage | 500 | 8,737 | 1,313 | 584 | 579 | 6/6 | 2 | CD68:42; CSF1R:9; MRC1:210; MSR1:191; CD163:546; ITGAM:326 |
| B_cell | 416 | 8,040 | 4,121 | 1,265 | 1,133 | 5/5 | 2 | MS4A1:14; CD79A:27; CD79B:318; CD19:757; BANK1:181 |
| NK_cytotoxic_T | 500 | 7,835 | 2,979 | 963 | 920 | 6/6 | 1 | CD3D:26; CD3E:138; CD8A:990; IL7R:66; NKG7:73; GNLY:824 |
| CD4_T | 500 | 6,643 | 2,051 | 141 | 85 | 0/0 | 0 | — |

Spearman(coverage, padj < 0.05 count) = +0.83 (p = 0.042); 21/22 own markers in list,
7/22 in the top 50. The Wilcoxon test is rank-based, so per-cell scaling changes the numbers only
through ties and the fold-change filter; the depth ordering of the marker counts is the same under both
layers. Note the contrast with the profile read-outs of Section 1: normalisation fixes the *argmax* (the
per-type mean), but the *number* of significant markers per type keeps tracking coverage because the test
counts genes whose detection rate differs, and detection rate is a coverage quantity.

## 3. Signac (`03_filter_fragments_and_index.sh`, `04_signac_hcc.R`)

What was feasible and how:
- Fragment file: `GSE227265_fragments_AllSamples.tsv.gz` is coordinate-sorted, bgzipped, and a `.tbi`
  index exists (`...tsv.gz.tbi`, 4.3 MB) — Signac/Rsamtools `scanTabix` uses it (0.5–0.7 s per 500 kb
  query). But the file holds every unfiltered barcode of all 13 samples (~1.5e9 lines, 62 GB
  uncompressed; the first 2e6 lines alone contain 358,676 distinct barcodes), so `FeatureMatrix` over
  ~20k gene bodies would parse ~40% of it as R character vectors. The documented remedy is
  `Signac::FilterCells()` (subset a fragment file to chosen cells); it is single-threaded R over the whole
  file, so the same operation was done with `grep -F -w` + `mawk` on the barcode column, then
  `pysam.tabix_compress` + `pysam.tabix_index(preset='bed')` (11 min; 1.56 GB bgzipped, 12,029 cells,
  all cells retained — no 3,000-cell subsampling was needed).
- Annotation: `EnsDb.Hsapiens.v86` (the vignette's `GetGRangesFromEnsDb`) is not installed and could not be
  installed, so the GENCODE v38 gene records (`type == "gene"`, 60,649 genes) were imported with
  `rtracklayer` and given the columns Signac reads (`gene_id`, `gene_name`, `gene_biotype`, `type`).
  Signac's `CollapseToLongestTranscript()` takes min-start/max-end per `gene_id`, which for EnsDb
  transcripts is the gene span, so a gene record is the same input. 19,955 protein-coding genes (of 60,649) on standard
  chromosomes.
- Object: `CreateFragmentObject(path, cells)` (validated all 12,029 cells against the index),
  `FeatureMatrix` over chr22 protein-coding genes as the container counts, `CreateChromatinAssay(counts,
  fragments, annotation)`, `CreateSeuratObject(assay='peaks')`, then **`GeneActivity(sobj)` with all
  defaults** (`extend.upstream = 2000`, `extend.downstream = 0`, `biotypes = 'protein_coding'`,
  `max.width = 5e5`) → **19,622 genes × 12,029 cells, 54.2e6 non-zeros**
  (13 genes on non-standard contigs dropped; 1 duplicated gene name made unique by `CreateAssayObject`).
  Then the vignette's `NormalizeData(assay='RNA', normalization.method='LogNormalize',
  scale.factor = median(nCount_RNA))` and `AverageExpression(group.by='cell_type')`.

How Signac's matrix relates to the SnapATAC2/paper matrix (`08_compare_signac_snapatac2.py`,
`results/signac_vs_snapatac2_matrix.csv`): 19,594 gene names in common; per-cell totals on the common
genes have median ratio 1.005 and Spearman 1.000; per-gene totals Pearson(log1p) 0.9999; only 1.0% of
non-zero entries differ (Signac counts fragment overlaps, SnapATAC2 paired insertions; both use gene body
+ 2 kb upstream and the same GENCODE spans). 18 of the 22 panel genes have identical totals; MSR1, IL7R,
GNLY, ALB, HNF4A, CD79A differ by 0.1–8.6%. **BANK1 is absent** from Signac's output because its GENCODE
v38 span is 663,527 bp > `max.width = 5e5` — hence 21, not 22, panel genes are scorable.

Read-outs: vignette path (AverageExpression default = expm1 mean): top share **25.7%** (Hepatocyte),
markers **18/21** (CD19, IL7R, APOA1 wrong, exactly the paper's cpm_mean misassignments); scale factor
median vs 1e4 makes no difference to the argmax. Raw counts mean: **72.4%** in Endothelial_stromal (the
highest-coverage type), 15/21, macrophage 2/6 — the same collapse as raw_mean on the full matrix. Mean of the log
`data` layer: 27.8% (median) / 27.9% (1e4), 18/21. The lower top shares relative to the paper (25.7 vs
30.1%; 27.8 vs 33.6%) come from the denominator: Signac's matrix contains only protein-coding genes
<= 500 kb (19,622 expressed genes vs 58,450), and the non-coding / lncRNA rows that the full matrix carries
are the ones most concentrated in the hepatocyte compartment. The raw-mean share is *higher* (72.4 vs
67.1%) for the same reason in reverse. Conclusion: the Signac vignette workflow yields the CPM-mean
behaviour of the paper; a Signac user only reaches the coverage-dominated profile by averaging the counts
assay (or the pooled sums) directly.

## 4. SnapATAC2 (`05_snapatac2_hcc.py`, `results/snapatac2_results.csv`)

`snap.read(hcc_fragments.h5ad, backed='r')` (12,029 cells; `obs['n_fragment']` median 8,133), then
`snap.pp.make_gene_matrix(data, gene_anno=GENCODE v38 GTF)` with **every argument at its default**
(`upstream=2000`, `downstream=0`, `include_gene_body=True`, `id_type='gene'`,
`counting_strategy='paired-insertion'`): 12,029 × 59,385, uint32, 0 duplicated gene names.
**The tool's default output is byte-identical to the manuscript's count matrix: 0 of 79,413,387 non-zero
entries differ** (the paper's Script 90 called the same function with the same, default, arguments).
So every HCC count-matrix number in the paper is a statement about SnapATAC2's own default gene matrix.

Aggregating that default output:
- raw mean: 67.1% of argmaxima in Endothelial_stromal (highest-coverage type), 15/22 (macrophage 2/6);
- the SnapATAC2 tutorial's normalisation (`sc.pp.normalize_total()` with the default
  `target_sum=None`, i.e. the median total 8,766, then `sc.pp.log1p()`), mean per type: 33.4%, 19/22
  (33.6% with `target_sum=1e4` — the paper's lognorm_mean);
- CPM mean (`normalize_total` without log, 1e4 or median target): 30.1%, 19/22 — the paper's cpm_mean;
- pooled CPM: 31.6%, 20/22.
(The tutorial continues with `sc.external.pp.magic` imputation before plotting marker genes; MAGIC was
not run — it is a smoothing of the log-normalised matrix and does not change the aggregation question.)

## 5. scanpy (`06_scanpy_groupby_check.py`, `results/scanpy_results.csv`)

On `sc.pp.normalize_total(target_sum=1e4)` + `sc.pp.log1p` data:
- `sc.get.aggregate(adata, by='cell_type', func='mean')` vs the paper's `m_lognorm`: max |diff| 3.3e-6
  (float32 accumulation), top share 33.6%, 19/22;
- `sc.get.obs_df(keys=genes+['cell_type']).groupby('cell_type').mean()` (in 4,000-gene blocks): max |diff|
  3.2e-6, 33.6%, 19/22;
- `sc.pl.DotPlot(adata, var_names=<22 markers>, groupby='cell_type').dot_color_df`: max |diff| 4.2e-7
  vs `m_lognorm`; colour argmax correct 19/22 (CD19, IL7R, APOA1 wrong — same three); the dot-size channel
  (fraction of cells with a non-zero value) is correct 18/22 and additionally sends MRC1 to
  Endothelial_stromal and CD163 to Hepatocyte, i.e. the detection-rate channel of the dotplot carries the
  depth signal even after normalisation (`results/scanpy_dotplot_markers.csv`).
So the paper's lognorm_mean is exactly what scanpy's group summaries and dotplots display. The
`rank_genes_groups` Wilcoxon result already in `results/benchmark_wilcoxon.csv` was not re-run.

## 6. What this establishes

1. The re-implementations are faithful: for every documented path the real tool reproduces the paper's
   number to numerical precision (Seurat AverageExpression default = cpm_mean 30.1%/19; layer='counts' =
   raw_mean 67.1%/15; scanpy aggregate/obs_df/dotplot = lognorm_mean 33.6%/19; SnapATAC2 default matrix =
   the paper's matrix, 0 differing entries).
2. Attribution correction for Table 2 of the paper: Seurat `AverageExpression()` averages expm1 of the
   log-normalised layer (CPM mean); it is the *mean of the log layer* that corresponds to
   scanpy/Signac dotplot-style summaries. Seurat v5's recommended `AggregateExpression()` returns raw sums,
   whose argmax follows population size × coverage (41.0% in the largest type, B-cell markers 1/5) until the
   user normalises the pooled sum (`return.seurat=TRUE`: 31.6%/20).
3. The Signac vignette (GeneActivity -> NormalizeData with median scale factor -> AverageExpression) lands
   on the CPM-mean behaviour (25.7%, 18/21); its matrix is the SnapATAC2 default restricted to
   protein-coding genes <= 500 kb (99% identical entries), which drops BANK1.
4. The raw-count collapse is reproduced identically by all three real tools when their counts are averaged
   (Seurat 67.1%, Signac 72.4%, SnapATAC2 67.1%; macrophage panel 2/6 in each).
5. `FindAllMarkers(only.pos = TRUE)` with defaults, all cells (42 min): the number of significant positive
   markers per type follows true coverage (rho +0.83, p = 0.042; +0.89 with the log2FC > 0.5 filter),
   independent of population size (confirmed at 500 cells/type: rho +0.77 data, +0.83 counts);
   21/22 canonical markers are *in* their own list but only 4/22 are in the top 50, so a marker list
   read from the top does not recover the lineage panel under either layer.
6. ArchR could not be run (not installed); its `GeneScoreMatrix` convention (per-cell scaling to 1e4,
   log2(x+1)) remains a re-implementation only (`archr_style`, 33.6%/19 — identical argmax to lognorm_mean
   because log2 and log1p differ by a constant factor after the same per-cell scaling).
