#!/usr/bin/Rscript
# FindAllMarkers with Seurat defaults on the reduced problem (<=500 cells/type, expressed protein-coding genes).
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(future); library(data.table)})
options(future.globals.maxSize = 8 * 1024^3)
SCR <- "/tmp/claude-1000/-mnt-10t-scrna-atac/b849b398-690b-4bff-b242-236b530a2a05/scratchpad/hcc_sub"
OUT <- Sys.getenv("REAL_TOOLS_OUT", "analyses/real_tools/results")
lg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
m <- ReadMtx(mtx = file.path(SCR, "matrix.mtx"), cells = file.path(SCR, "cells.tsv"), features = file.path(SCR, "features.tsv"),
             cell.column = 1, feature.column = 1, skip.cell = 1)
meta <- read.delim(file.path(SCR, "cells.tsv"), stringsAsFactors = FALSE); rownames(meta) <- meta$barcode
obj <- CreateSeuratObject(counts = m, meta.data = meta[colnames(m), ]); Idents(obj) <- "cell_type"
lg("reduced object:", nrow(obj), "genes x", ncol(obj), "cells"); print(table(Idents(obj)))
obj <- NormalizeData(obj, verbose = FALSE)
plan("multicore", workers = 6)
t1 <- Sys.time()
fm <- FindAllMarkers(obj, only.pos = TRUE, verbose = FALSE)      # wilcox, slot='data', logfc.threshold=0.1, min.pct=0.01
lg("data layer done in", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min; rows", nrow(fm)); print(table(fm$cluster))
fwrite(fm, file.path(OUT, "seurat_findallmarkers_reduced_data.csv"))
t1 <- Sys.time()
fmc <- tryCatch(FindAllMarkers(obj, only.pos = TRUE, slot = "counts", verbose = FALSE),
                error = function(e) { lg("slot='counts' failed:", conditionMessage(e)); NULL })
if (!is.null(fmc)) { lg("counts layer done in", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min; rows", nrow(fmc)); print(table(fmc$cluster))
  fwrite(fmc, file.path(OUT, "seurat_findallmarkers_reduced_counts.csv")) }
lg("ALL DONE")
