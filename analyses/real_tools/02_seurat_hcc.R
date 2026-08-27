#!/usr/bin/Rscript
# Real Seurat (v5) run on the HCC proximal gene-activity counts (SnapATAC2 make_gene_matrix output,
# exported by 01_export_hcc_mtx.py). Everything below is the documented default unless stated.
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(future); library(data.table)})
options(future.globals.maxSize = 16 * 1024^3)
SCR <- "/tmp/claude-1000/-mnt-10t-scrna-atac/b849b398-690b-4bff-b242-236b530a2a05/scratchpad"
OUT <- Sys.getenv("REAL_TOOLS_OUT", "analyses/real_tools/results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
save_profile <- function(M, name) {                       # genes x types -> CSV, group names back to underscores
  M <- as.matrix(M); colnames(M) <- gsub("-", "_", colnames(M))
  fwrite(data.table(gene = rownames(M), M), file.path(OUT, paste0("seurat_profile_", name, ".csv")))
  lg("saved profile", name, "dims", paste(dim(M), collapse = "x"), "cols:", paste(colnames(M), collapse = ","))
}
maxdiff <- function(A, B) { A <- as.matrix(A); B <- as.matrix(B); colnames(A) <- gsub("-", "_", colnames(A)); colnames(B) <- gsub("-", "_", colnames(B)); max(abs(A[rownames(B), colnames(B)] - B)) }

# ── load ─────────────────────────────────────────────────────────────────────
m <- ReadMtx(mtx = file.path(SCR, "hcc_mtx/matrix.mtx"), cells = file.path(SCR, "hcc_mtx/cells.tsv"),
             features = file.path(SCR, "hcc_mtx/features.tsv"), cell.column = 1, feature.column = 1,
             skip.cell = 1, skip.feature = 0, unique.features = TRUE)
meta <- read.delim(file.path(SCR, "hcc_mtx/cells.tsv"), stringsAsFactors = FALSE); rownames(meta) <- meta$barcode
lg("matrix", nrow(m), "genes x", ncol(m), "cells; nnz", length(m@x), "; max", max(m@x))
obj <- CreateSeuratObject(counts = m, meta.data = meta[colnames(m), ])   # defaults: min.cells = 0, min.features = 0
Idents(obj) <- "cell_type"
lg("cells per type:"); print(table(obj$cell_type))
lg("median nCount_RNA (matrix row sum) per type:"); print(round(tapply(obj$nCount_RNA, obj$cell_type, median)))
lg("median TRUE depth (fragments) per type:"); print(round(tapply(obj$true_depth, obj$cell_type, median)))

# ── NormalizeData: LogNormalize, scale.factor = 1e4 (defaults) ───────────────
obj <- NormalizeData(obj)
D <- LayerData(obj, layer = "data"); C <- LayerData(obj, layer = "counts")
grp <- split(colnames(obj), obj$cell_type)
manual_expm1_mean <- sapply(grp, function(cl) Matrix::rowMeans(expm1(D[, cl])))   # = CPM(1e4) mean
manual_log_mean   <- sapply(grp, function(cl) Matrix::rowMeans(D[, cl]))          # = mean of log1p(CP10k)
manual_raw_mean   <- sapply(grp, function(cl) Matrix::rowMeans(C[, cl]))
manual_sum        <- sapply(grp, function(cl) Matrix::rowSums(C[, cl]))

# ── (a) AverageExpression / AggregateExpression ───────────────────────────────
lg("AverageExpression(group.by='cell_type')  [default layer='data']")
ae_default <- AverageExpression(obj, group.by = "cell_type")$RNA
save_profile(ae_default, "AverageExpression_default")
cat("  max|AE_default - mean(expm1(data))| =", maxdiff(ae_default, manual_expm1_mean),
    "   max|AE_default - mean(data)| =", maxdiff(ae_default, manual_log_mean), "\n")

lg("AverageExpression(layer='counts')")
ae_counts <- AverageExpression(obj, group.by = "cell_type", layer = "counts")$RNA
save_profile(ae_counts, "AverageExpression_counts")
cat("  max|AE_counts - mean(counts)| =", maxdiff(ae_counts, manual_raw_mean), "\n")

lg("AverageExpression on a copy of the data layer stored under another name (mean of log values, no expm1)")
obj[["RNA"]]$logdata <- D
ae_logmean <- tryCatch(AverageExpression(obj, group.by = "cell_type", layer = "logdata")$RNA,
                       error = function(e) { lg("  layer='logdata' failed:", conditionMessage(e)); NULL })
if (!is.null(ae_logmean)) {
  save_profile(ae_logmean, "AverageExpression_logdata_layer")
  cat("  max|AE_logdata - mean(data)| =", maxdiff(ae_logmean, manual_log_mean),
      "   max|AE_logdata - mean(expm1(data))| =", maxdiff(ae_logmean, manual_expm1_mean), "\n")
}
# the only documented way to average log values without expm1: return.seurat=TRUE gives log1p(mean(expm1)), NOT mean(log)
ae_seu <- AverageExpression(obj, group.by = "cell_type", return.seurat = TRUE)
save_profile(LayerData(ae_seu, layer = "data"), "AverageExpression_returnseurat_data")
cat("  max|AE(return.seurat) data - log1p(mean(expm1(data)))| =", maxdiff(LayerData(ae_seu, layer = "data"), log1p(manual_expm1_mean)), "\n")
save_profile(manual_log_mean, "manual_mean_of_data_layer")

lg("AggregateExpression(group.by='cell_type')  [sum of counts]")
agg <- AggregateExpression(obj, group.by = "cell_type")$RNA
save_profile(agg, "AggregateExpression_sum")
cat("  max|Agg - colSums by type| =", maxdiff(agg, manual_sum), "\n")
agg_seu <- AggregateExpression(obj, group.by = "cell_type", return.seurat = TRUE)   # data = LogNormalize(sums, 1e4)
save_profile(LayerData(agg_seu, layer = "data"), "AggregateExpression_returnseurat_data")
rm(D, C); gc()

# ── (b) FindAllMarkers ────────────────────────────────────────────────────────
plan("multicore", workers = 6)
lg("FindAllMarkers(only.pos=TRUE)  [test.use='wilcox', slot='data', logfc.threshold=0.1, min.pct=0.01]")
t1 <- Sys.time()
fm_data <- FindAllMarkers(obj, only.pos = TRUE, verbose = FALSE)
lg("  done in", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min; rows", nrow(fm_data))
fwrite(fm_data, file.path(OUT, "seurat_findallmarkers_data.csv"))
print(table(fm_data$cluster))

lg("FindAllMarkers(only.pos=TRUE, slot='counts')  [Wilcoxon on raw counts]")
t1 <- Sys.time()
fm_counts <- tryCatch(FindAllMarkers(obj, only.pos = TRUE, slot = "counts", verbose = FALSE),
                      error = function(e) { lg("  slot='counts' failed:", conditionMessage(e)); NULL })
if (!is.null(fm_counts)) {
  lg("  done in", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min; rows", nrow(fm_counts))
  fwrite(fm_counts, file.path(OUT, "seurat_findallmarkers_counts.csv"))
  print(table(fm_counts$cluster))
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "seurat_sessionInfo.txt"))
lg("ALL DONE")
