# Supporting analyses

Five independent analyses, each in its own directory with its scripts, summary
tables and a `REPORT.md`. None depends on the results of another, and each can be
re-run on its own. All import the aggregation functions, marker panels and
permutation null from `tools/benchmark_aggregation.py`, so the strategy
definitions are identical to the main benchmark's.

| directory | question | headline |
|---|---|---|
| [`depth_spread/`](depth_spread/REPORT.md) | How does each strategy's failure scale with an *imposed* coverage spread of 1–40×? | Raw and z-scaled means reach 68% (HCC) and 93% (NSCLC) top share at 40×; the log-normalised mean 50% and 83%; the CPM mean and pooled CPM stay flat at 22–30% across the whole range. |
| [`real_tools/`](real_tools/REPORT.md) | What do Seurat, Signac, scanpy and SnapATAC2 actually produce on their documented default paths? | Every documented path reproduces the corresponding re-implementation to numerical precision. SnapATAC2's default output is byte-identical to the manuscript's HCC count matrix. Seurat's `AverageExpression` default is the CPM mean, not the log mean. |
| [`multiome/`](multiome/REPORT.md) | Genome-wide RNA ground truth in the same cells (10x PBMC Multiome): ATAC-vs-RNA argmax agreement and annotation accuracy. | Even at a 1.66× coverage spread, raw and z-scaled means put 34% (count) and 45% (weighted) of argmaxima in the highest-coverage type against the 22.5% the RNA truth shows; every normalised strategy returns 18–23%. |
| [`normalisers_and_ci/`](normalisers_and_ci/REPORT.md) | Ten further transforms (TF-IDF, Pearson residuals, median scaling, rank, detection rate) plus bootstrap intervals and paired tests. | Only transforms linear in the per-cell fraction are coverage-safe. Signac's default TF-IDF is *worse* than the log-normalised mean on wide-spread matrices (91–95% on COAD); Pearson residuals are the only near-safe non-linear transform. |
| [`practice_survey/`](practice_survey/REPORT.md) | Is the aggregation scale standardised in published pipelines and toolkits? | Of 28 surveyed items, 8 average on the log scale, 1 uses a raw mean by default, 1 publishes a z-scaled average of an unstated scale, and 3 do not state the scale at all. No prior work benchmarks this step. |

> Naming: these directories predate the manuscript's final terminology.
> Identifiers in scripts and CSV columns keep the older names: `depth` = coverage,
> `pseudobulk` = pooled CPM, `proximal` = count matrix, `enhancer` = weighted
> matrix, `LUAD` = NSCLC. Each `REPORT.md` repeats this note.
