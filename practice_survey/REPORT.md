# What pipelines actually do: a survey of the aggregation scale for scATAC gene activity

**Purpose.** Rebuttal evidence for the anticipated reviewer objection that *"depth
normalisation is common knowledge."* The benchmark's claim is not that depth
normalisation is unknown, but that **the aggregation step — the scale on which
single cells are combined into a per-cell-type profile — is unstandardised in
practice, and a large share of real pipelines perform it on a scale that reports
sequencing depth or detection rate.** This document surveys 28 items (14 toolkit
functions, 6 deconvolution/reference methods, 8 primary papers) and records, for
each, the exact code line or Methods sentence that establishes the scale.

Compiled 2026-08-27. Companion machine-readable file: `survey_table.csv`.

**Evidence standard.** Every row carries a verbatim quote and a confidence level.
Where a quote could not be obtained the row is marked `unstated` and the
"Could not verify" section below lists exactly what was checked. Nothing in the
table is inferred from a tool's reputation; where a classification rests on a
documented default rather than the paper's own words, the confidence column says so.

---

## 1. Headline counts

Classification is by the scale on which the **per-cell-type / per-cluster summary**
is formed in that pipeline's documented default workflow.

| aggregation class | n | depth-invariant? |
|---|---:|---|
| linear-CPM-mean (per-cell linear normalisation, then mean) | 10 | yes |
| **log-mean** (mean or median of log-transformed values) | **8** | **no — tracks detection rate** |
| pseudobulk-CPM (sum cells, then normalise) | 4 | yes |
| **unstated** (normalisation named but not specified, or absent) | **3** | **unknown** |
| **raw mean** (mean of counts, no per-cell normalisation) | **1** | **no — tracks depth** |
| **z-scaled average** (of an unstated underlying scale) | **1** | **no** |
| bulk-sorted reference (quantile-normalised; not scATAC) | 1 | n/a |
| **total** | **28** | |

**13 of 28 surveyed items (46%) aggregate on a scale the benchmark identifies as
depth-reporting, or leave the scale unverifiable.** Restricting to the eight
primary research papers: 3 are depth-safe only because they inherit an ArchR
default they never state, 2 are on the log scale, 1 publishes a *z-scaled average
of an unstated scale*, 1 is unstated outright, and 1 (Zhang 2021) is explicitly
pseudobulk-CPM.

The split is not random. It falls along a clean line:

- **Toolkits that normalise per cell before storing the matrix are safe**: ArchR
  (`scaleTo = 10000`), Cicero (`normalize_gene_activities`), and Seurat's
  `AverageExpression`/`AggregateExpression`/`DotPlot`/`FoldChange`, all of which
  either exponentiate the log slot before averaging or sum counts first.
- **Toolkits that hand the user a raw count matrix and a log-normalised slot are
  where the failure enters**: Signac's `GeneActivity`, SnapATAC2's
  `make_gene_matrix`, and epiScanpy's `geneactivity` all return counts; the
  vignettes then apply `LogNormalize`/`normalize_total`+`log1p`, and every scanpy
  group summary downstream (`rank_genes_groups`, `pl.dotplot`) averages the logged
  values. scanpy's own docstring concedes the point: the fold change is *"an
  approximation calculated from mean-log values."*

This is precisely the boundary the paper's Discussion identifies — *"The failure
modes arise between the tools."* The survey substantiates it with quotes.

### The single sharpest contrast for the rebuttal

Two dot-plot functions, same conceptual operation, opposite scales:

- **Seurat** `DotPlot`: `avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) { return(mean(x = expm1(x = x))) })` — mean of **exponentiated** values (CPM mean).
- **scanpy** `sc.pl.dotplot`: `dot_color_df = self.obs_tidy.groupby(level=0, observed=True).mean()` — mean of **log** values, no `expm1` anywhere.

A reviewer claiming the choice is common knowledge has to explain why the two most
used single-cell toolkits make it differently and neither warns the user.

### A trap inside a "safe" toolkit

ArchR is depth-safe for gene scores only because `addGeneScoreMatrix` pre-normalises
each cell. Its own group-summary export does not:

```r
# getGroupSE() signature:  divideN = TRUE, scaleTo = NULL
if(divideN){ groupMat <- t(t(groupMat) / as.vector(nCells)) }          # always
if(!is.null(scaleTo)){ groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo }  # opt-in
```

With `useMatrix = "PeakMatrix"` or `"TileMatrix"` and the default `scaleTo = NULL`,
`getGroupSE` returns a per-group mean of **raw counts**. The depth normalisation is
an opt-in argument, not a default.

---

## 2. Evidence table

Full quotes, DOIs and per-row notes are in `survey_table.csv` (columns:
`id, category, tool_or_paper, version_year, doi_or_url, evidence_quote,
matrix_scale, aggregation_class, confidence, notes`). Abridged here:

### (a) Toolkits

| # | Tool / function | Version | Key evidence (abridged; full quote in CSV) | Class | Confidence |
|---|---|---|---|---|---|
| 1 | Signac `GeneActivity()` + PBMC/brain vignette | 1.17.1 / 2021 | Rd: "Compute counts per cell in gene body and promoter region"; source ends `return(counts)`. Vignette: `NormalizeData(..., normalization.method = 'LogNormalize', scale.factor = median(pbmc$nCount_RNA))` then `FeaturePlot(...)` | log-mean (slot handed on) | verified from code (local) + docs |
| 2 | Seurat `AverageExpression()` | 5.5.0 | `if (layer == "data") { data.use <- expm1(x = data.use) }`; Rd: "feature values are exponentiated prior to averaging so that averaging is done in non-log space" | linear-CPM-mean | verified from code (local) |
| 3 | Seurat `AggregateExpression()` | 4.0.0+/5.x | NEWS 4.0.0: "Add `AggregateExpression()` for summation based pseudobulk calculations"; `method == "aggregate"` → sum then `NormalizeData` | pseudobulk-CPM | verified from code + docs |
| 4 | Seurat `DotPlot()` | 5.5.0 | `return(mean(x = expm1(x = x)))`; colour then `scale(x = log1p(data.use))` | linear-CPM-mean (z for colour) | verified from code (local) |
| 5 | Seurat `FindMarkers()`/`FoldChange()` | 5.5.0 | `log1pdata.mean.fxn <- function(x) log((rowSums(expm1(x)) + pseudocount.use)/NCOL(x), base)`; default `slot = "data"` | linear-CPM-mean | verified from code (local) |
| 6 | scanpy `sc.tl.rank_genes_groups()` | 1.11.5 | means via `_get_mean_var(x_mask)` on `adata.X`; `foldchanges = (expm1_func(mean_group)+1e-9)/(expm1_func(mean_rest)+1e-9)`; docstring: "this is an approximation calculated from mean-log values" | **log-mean** | verified from code (local) |
| 7 | scanpy `sc.pl.dotplot()` | 1.11.5 | `dot_color_df = self.obs_tidy.groupby(level=0, observed=True).mean()` | **log-mean** | verified from code (local) |
| 8 | ArchR `addGeneScoreMatrix()` | 1.0.x / 2021 | `scaleTo = 10000`; "Each column ... will be normalized to a column sum designated by scaleTo"; `matGS@x <- as.numeric(scaleTo * matGS@x/rep.int(totalGS, ...))` | linear-CPM (per cell) | verified from code + docs |
| 9 | ArchR `getMarkerFeatures()` | 1.0.x | `normFactors <- scaleTo / mColSums`; `log2FC <- log2((m1 + offset) / (m2 + offset))` — log **of** the linear mean | linear-CPM-mean | verified from code |
| 10 | ArchR `getGroupSE()` / `plotEmbedding` / `plotBrowserTrack` | 1.0.x | `divideN = TRUE, scaleTo = NULL`; depth rescale is opt-in. `plotBrowserTrack` default `normMethod="ReadsInTSS"` | **raw mean by default** (tracks: pseudobulk-CPM) | verified from code |
| 11 | SnapATAC2 `pp.make_gene_matrix()` + PBMC tutorial | 2.9.0 / 2024 | Docstring: "counting the TN5 insertions in each gene's regulatory domain"; tutorial `normalize_total` → `log1p` → `magic` → `sc.pl.umap` | **log-mean** | verified from code (local) + docs |
| 12 | epiScanpy `epi.tl.geneactivity()` | 0.3.x / 2021 | "the number of open features ... overlapping genes"; `np.sum(matrix[:, indice].todense(), axis=1)`; no normalisation in function | **log-mean** (raw out, scanpy logs it) | verified from code (GitHub) |
| 13 | Cicero `build_gene_activity_matrix` + `normalize_gene_activities` | 1.x / 2018 | `lm(log(total_activity) ~ log(total_sites))`; `scores@x <- pmin(1e9, exp(scores@x) - 1)`; `scale_factors <- Matrix::Diagonal(x=1/sum_activity_scores)` | linear per-cell (columns sum to 1) | verified from code |
| 14 | pycisTopic `get_gene_activity()` (SCENIC+) | 2.x / 2023 | weighted dot product on `imputed_acc_object.mtx`; `gene_act = gene_act * scale_factor` (default 1); no log1p, no per-cell norm | **unstated** | verified from code |

### (b) Deconvolution / reference-profile methods

| # | Method | Year | Key evidence (abridged) | Class | Confidence |
|---|---|---|---|---|---|
| 15 | EPIC-ATAC | 2024, eLife | "dividing counts by peak length, correcting sample counts for depth and rescaling counts so that the counts of each sample sum to 10^6"; "median of the TPM-like counts across all samples from each cell type" | linear TPM | from Methods |
| 16 | DeconPeaker | 2020, Front Genet | "quantile normalization on the count matrix across all reference samples"; "the average of the peaks in each cell type's samples" | bulk-sorted reference | from Methods |
| 17 | CIBERSORTx (single-cell signature) | 2019, Nat Biotech | "aggregated the profiles by summation in **non-log linear space** and normalized each population-level GEP into TPM" | pseudobulk-CPM | from Methods |
| 18 | EPISCORE | 2020, Genome Biol | per-cell rescale to a common total, "+1", "log2 transformation"; reference = "the median over all cells within each cluster (cell type)" | **log-mean (median of logs)** | from Methods |
| 19 | DECA | 2025, Brief Bioinform | "min-max normalization strategy to adjust for discrepancies in sequencing depth"; per-cell aggregation scale never stated | **unstated** | from Methods |
| 20 | scATAcat | 2024, NAR GAB | "adding the read coverage for each feature across the cluster member cells"; then "library size normalization and log2 transformation to the pseudobulk matrix" | pseudobulk-CPM | from Methods |

### (c) Primary papers, 2022–2026

| # | Paper | Year / journal | Key evidence (abridged) | Class | Confidence |
|---|---|---|---|---|---|
| 21 | Muto et al., ADPKD | 2022, Nat Commun | "The gene activity matrix was log-normalized." "Each cluster was annotated based on gene activities." | **log-mean** | from Methods |
| 22 | Zhang et al., human chromatin accessibility atlas | 2021, Cell | "We aggregated open chromatin fragments from each cluster and utilized the promoter accessibility, defined as RPM of +/- 1kb around TSS, as the proxy for gene activity" | pseudobulk-CPM | from Methods |
| 23 | Braband et al., tissue CD4+ T cells | 2023, Front Immunol | `addGeneScoreMat = TRUE`; `getMarkerFeatures(..., useMatrix = 'GeneScoreMatrix', ..., testMethod = 'wilcoxon')` | linear-CPM-mean (ArchR default) | Methods + documented default |
| 24 | Suen et al., gonadal lineage | 2024, BMC Genomics | "Gene scores were calculated using the addGeneScoreMatrix function"; getMarkerFeatures with `bias` = TSSEnrichment, log10(nFrags) | linear-CPM-mean (ArchR default) | Methods + documented default |
| 25 | Farghli et al., fibrolamellar carcinoma | 2026, Sci Rep | "violin plots that depict ArchR calculated gene activity scores"; "`getMarkerFeatures` with default parameters" | linear-CPM-mean (ArchR default) | Methods + documented default |
| 26 | Camp et al., ccRCC | 2026, Sci Adv | "we calculated gene activity scores in Signac"; Fig 1D: "**scaled average gene activity scores** across broad nontumor cell types" | **z-scaled average of an unstated scale** | from Methods / figure legend |
| 27 | Ma et al., retina & brain microglia | 2025, Genome Res | "gene activity matrix was calculated using the GeneActivity function"; "normalized using NormalizeData with default parameters"; "module scores ... visualized using Seurat's DotPlot" | **log-mean** (module score on log slot) | from Methods |
| 28 | dbscATAC resource | 2025 | "After normalizing this gene matrix, the 'FindAllMarkers' function in Seurat was applied" — method never named | **unstated** | from Methods |

---

## 3. Is there a prior benchmark of this step?

Searched: Europe PMC (four structured queries over 2019–2026 combining
*gene activity / gene score* × *sequencing depth / library size* × *pseudobulk /
normalization* × *benchmark / bias / confound*), plus targeted web searches on
"benchmarking gene activity scores scATAC", "gene score depth bias", and the
scRNA-seq mean-of-logs literature. Findings:

**No prior work benchmarks the per-cell-type aggregation scale for scATAC gene
activity.** The four nearest candidates each stop short:

1. **Luo et al. 2024, *Genome Biology*** (10.1186/s13059-024-03356-x),
   "Benchmarking computational methods for single-cell chromatin data analysis."
   Benchmarks QC, feature engineering, dimension reduction, graph construction and
   clustering. It *documents* the per-tool convention — "Signac employs a basic
   technique, counting fragments in the gene body and promoter regions for each
   gene, followed by log-normalization of these counts to derive a gene activity
   score"; SnapATAC "normalizes them using log-transformed count-per-million reads"
   — and it quantifies depth bias in embeddings ("LSI-based methods (Signac, ArchR)
   showed a strong library size bias across all datasets"). It does **not** evaluate
   aggregation scale. **Verdict: does not pre-empt** — and it is a useful citation,
   because it confirms in a benchmark paper's own words that the toolkits' gene
   activity convention is log-normalisation.
2. **Kwok et al. 2025, *Genome Biology*** (10.1186/s13059-025-03735-y), hierarchical
   count model for scATAC. Directly attacks depth normalisation in scATAC:
   "TF-IDF transformation, being a rehash of CPM, suffers from the same issues";
   "we recommend against treating TF-IDF as a depth normalization method";
   "the resulting counts are not 'depth-normalized.' In many cases, the sequencing
   depth effect is even exaggerated." But this concerns TF-IDF on the **peak**
   matrix at the **cell** level; the paper confirms it does not discuss averaging
   cells within a cell type on raw vs normalised vs log scale.
   **Verdict: partly relevant, does not pre-empt** — cite it as evidence that depth
   handling in scATAC is contested rather than settled.
3. **Teo et al. 2024, *Nat Commun*** (10.1038/s41467-024-53089-5), "Best practices
   for differential accessibility analysis in single-cell epigenomics." Evaluates
   gene-level activity scores and tests sequencing depth as a performance factor,
   but compares statistical *tests*, not aggregation scales; the aggregation scale
   for its own pseudobulks is not specified in the accessible text.
   **Verdict: does not pre-empt.**
4. **Benchmarking automated cell type annotation for scATAC-seq** (Genes 2022,
   PMC9792779) reports that "all methods were sensitive to lower sequencing depth"
   but does not describe reference-profile construction. The gene-set-scoring
   benchmark (GPB 2024, PMC11423854) explicitly does neither: it does not evaluate
   depth/library-size effects on scores and does not benchmark per-cell-type
   aggregation.

**The mechanism, however, is established in scRNA-seq** and the paper must cite it
rather than claim the phenomenon is new:

- **Lun 2018** (bioRxiv 10.1101/404962), "Overcoming systematic errors caused by
  log-transformation of normalized single-cell RNA sequencing data": the mean of
  log-counts is not the log of the mean count, and the discrepancy varies with
  coverage, producing spurious differences between groups of cells.
- **Booeshaghi & Pachter 2021**, *Bioinformatics* (10.1093/bioinformatics/btab085).
- **Ahlmann-Eltze & Huber 2023**, *Nat Methods* 20:665–672 (10.1038/s41592-023-01814-1),
  comparison of transformations for scRNA-seq.
- **Crowell et al. 2020** *muscat* (10.1038/s41467-020-19894-4) and **Squair et al.
  2021** (10.1038/s41467-021-25960-2) for pseudobulk over cell-level aggregation.
- **CIBERSORTx** (row 17) already *requires* summation in non-log linear space.

### Novelty verdict

**The claim holds in its narrow form and must be stated narrowly.**

- ✅ Defensible: *no prior work benchmarks the aggregation scale of scATAC
  gene-activity matrices into per-cell-type profiles against sequencing depth, and
  none quantifies the resulting damage with an argmax-share read-out.* Nothing
  found contradicts this.
- ❌ Not defensible: *the mean-of-logs problem is unknown.* It is well established
  in scRNA-seq (Lun 2018 onward) and enforced as an input rule by CIBERSORTx.
- ✅ The genuinely new contributions are (i) that scATAC gene activity is sparse
  enough (1–18% detection per cell in the COAD matrix) for the effect to dominate
  rather than perturb, (ii) that the log-normalised mean — the aggregation most
  practitioners believe is safe — is *not* safe at realistic depth spreads, and
  (iii) the survey below, which shows the choice is genuinely unstandardised in
  the scATAC literature.

**Recommended framing for the rebuttal:** concede that depth normalisation is
common knowledge and that mean-of-logs is a known scRNA-seq pathology; then show
this table. The claim is about *practice*, not about *principle* — 8 of 28
surveyed pipelines average on the log scale, one on a raw mean by default, one
publishes a z-scaled average of an unstated scale, and three do not state the
scale at all. Common knowledge that half the field does not apply is worth a paper.

---

## 4. Could not verify (checked, not resolved)

These were investigated and are recorded as unverified rather than guessed:

- **Terekhanova et al. 2023, *Nature*** (10.1038/s41586-023-06682-5, PMC10632147),
  pan-cancer snATAC. Fetched the Nature page (auth redirect), the PMC article, and
  the Europe PMC full-text XML; in every case the retrieved Methods were truncated
  before the computational section, and the terms "gene activity", "GeneActivity",
  "ArchR", "Signac", "NormalizeData" did not appear in the retrieved text. The
  `ding-lab/PanCan_snATAC_publication` README describes no scripts or normalisation.
  **Not included in the table.**
- **Sussman et al. 2026, *Cell Rep Med*** (PMC13198283) — Methods truncated;
  Signac v1.14.0 appears in the software table but no gene-activity aggregation
  sentence was retrievable. **Not included.**
- **Li et al. 2026, *Nat Genet*** retina atlas (PMC13050529) — annotation is done
  by GLUE co-embedding + logistic regression, not by averaging gene activity; no
  aggregation-scale sentence found. **Not included** (not an instance of the
  operation surveyed).
- **Wu et al. 2023, *Nat Commun*** ccRCC (PMC10042888) — a figure legend refers to
  "gene activity of the epithelial and mesenchymal marker genes for tumor
  clusters", but the retrieved Methods excerpt cut off before the computation was
  described. **Not included.**
- **ArchR `projectBulkATAC`** — `R/ProjectBulkATAC.R` returned 404 on the master
  branch; the function may have been renamed or relocated. Not classified.
- **CIBERSORTx web documentation** (cibersortx.stanford.edu/tutorial.php) returned
  only navigation/ToS text; the "non-log linear space" requirement in row 17 is
  quoted from the Newman et al. 2019 paper Methods instead, which is the stronger
  source anyway.
- **SnapATAC2 notebook source** (`docs/tutorials/pbmc.ipynb`) is a Git-LFS pointer
  on GitHub; row 11's tutorial code was taken from the rendered documentation site
  (snapatac2.scverse.org) and cross-checked against the Galaxy Training Network
  rendering of the same pipeline, which lists the identical step order
  (make_gene_matrix → filter_genes → normalize_total → log1p → MAGIC → plot).
- **Rows 23–25** (Braband, Suen, Farghli): the papers state the ArchR functions but
  never state an aggregation scale of their own. Their classification rests on
  ArchR's documented defaults (rows 8–9), and the confidence column says so. This
  is itself a finding — the depth-safety of three of the eight primary papers is
  inherited from a toolkit default the authors never mention.

## 5. Reproduction

Sources are the installed packages (Seurat 5.5.0, Signac 1.17.1, scanpy 1.11.5,
snapatac2 2.9.0) inspected via `deparse()`/`inspect.getsource()`, GitHub raw source
for ArchR / Cicero / epiScanpy / pycisTopic, the tools' own rendered documentation,
and journal or PMC full text for the papers. Every quote in `survey_table.csv`
carries the URL or the local package it came from in the `doi_or_url` and
`confidence` columns.
