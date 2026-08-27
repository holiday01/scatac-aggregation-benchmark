#!/usr/bin/Rscript
# Signac vignette workflow on the HCC fragment file (subset to the 12,029 annotated cells):
#   CreateFragmentObject -> CreateChromatinAssay(fragments, annotation) -> GeneActivity() [defaults:
#   gene body + 2 kb upstream, protein_coding, max.width 5e5] -> CreateAssayObject -> NormalizeData(LogNormalize,
#   scale.factor = median(nCount_RNA), as in the vignette) -> AverageExpression(group.by = cell_type).
# Annotation: EnsDb.Hsapiens.v86 is not installed, so the GENCODE v38 GTF (gene records) is used, with the
# columns Signac expects (gene_id, gene_name, gene_biotype, type). Gene span = min start / max end of the gene,
# which is what CollapseToLongestTranscript() produces from transcript-level EnsDb ranges.
suppressPackageStartupMessages({library(Signac); library(Seurat); library(GenomicRanges); library(rtracklayer)
  library(future); library(Matrix); library(data.table)})
options(future.globals.maxSize = 16 * 1024^3)
plan("multicore", workers = 8)
SCR <- "/tmp/claude-1000/-mnt-10t-scrna-atac/b849b398-690b-4bff-b242-236b530a2a05/scratchpad"
OUT <- Sys.getenv("REAL_TOOLS_OUT", "analyses/real_tools/results")
lg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
FRAG <- file.path(SCR, "hcc_frag/hcc_12029.fragments.tsv.gz")
meta <- read.delim(file.path(SCR, "hcc_mtx/cells.tsv"), stringsAsFactors = FALSE); rownames(meta) <- meta$barcode
cells <- setNames(rownames(meta), rownames(meta))

lg("annotation from GENCODE v38 gene records")
gtf <- rtracklayer::import(file.path(SCR, "hcc_frag/gencode.v38.genes.gtf"))
ann <- gtf[gtf$type == "gene"]
mcols(ann) <- DataFrame(gene_id = ann$gene_id, gene_name = ann$gene_name, gene_biotype = ann$gene_type,
                        type = "gene", tx_id = ann$gene_id)
ann <- GenomeInfoDb::keepStandardChromosomes(ann, pruning.mode = "coarse")
lg("  ", length(ann), "genes;", sum(ann$gene_biotype == "protein_coding"), "protein_coding")

lg("CreateFragmentObject (validates the 12,029 cells against the index)")
frags <- CreateFragmentObject(path = FRAG, cells = cells, validate.fragments = TRUE)
lg("container counts: FeatureMatrix over chr22 protein-coding genes (only needed to build a ChromatinAssay)")
feat22 <- ann[seqnames(ann) == "chr22" & ann$gene_biotype == "protein_coding"]
cnt22 <- FeatureMatrix(fragments = frags, features = feat22, cells = cells, verbose = FALSE)
chrom <- CreateChromatinAssay(counts = cnt22, fragments = frags, annotation = ann)
sobj <- CreateSeuratObject(counts = chrom, assay = "peaks", meta.data = meta)
lg("object:", ncol(sobj), "cells")

lg("GeneActivity() with defaults")
t1 <- Sys.time()
ga <- GeneActivity(sobj)      # extend.upstream = 2000, extend.downstream = 0, biotypes = 'protein_coding', max.width = 5e5
lg("  done in", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min:", nrow(ga), "genes x", ncol(ga), "cells; nnz", length(ga@x))
dir.create(file.path(SCR, "signac_ga"), showWarnings = FALSE)
writeMM(ga, file.path(SCR, "signac_ga/matrix.mtx"))
writeLines(rownames(ga), file.path(SCR, "signac_ga/features.txt")); writeLines(colnames(ga), file.path(SCR, "signac_ga/cells.txt"))

sobj[["RNA"]] <- CreateAssayObject(counts = ga)
DefaultAssay(sobj) <- "RNA"
lg("median nCount_RNA (GeneActivity row sum) per type:"); print(round(tapply(sobj$nCount_RNA, sobj$cell_type, median)))
sobj <- NormalizeData(sobj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = median(sobj$nCount_RNA))  # vignette
save_profile <- function(M, name) { M <- as.matrix(M); colnames(M) <- gsub("-", "_", colnames(M))
  fwrite(data.table(gene = rownames(M), M), file.path(OUT, paste0("signac_profile_", name, ".csv"))); lg("saved", name, paste(dim(M), collapse = "x")) }
ae <- AverageExpression(sobj, group.by = "cell_type")$RNA                       # default layer = 'data' -> expm1 -> mean
save_profile(ae, "AverageExpression_default_scalefactor_median")
D <- GetAssayData(sobj, assay = "RNA", layer = "data"); grp <- split(colnames(sobj), sobj$cell_type)
save_profile(sapply(grp, function(cl) Matrix::rowMeans(D[, cl])), "mean_of_data_layer_scalefactor_median")
save_profile(sapply(grp, function(cl) Matrix::rowMeans(ga[, cl])), "mean_of_counts")
sobj <- NormalizeData(sobj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 1e4)
D <- GetAssayData(sobj, assay = "RNA", layer = "data")
save_profile(AverageExpression(sobj, group.by = "cell_type")$RNA, "AverageExpression_default_scalefactor_1e4")
save_profile(sapply(grp, function(cl) Matrix::rowMeans(D[, cl])), "mean_of_data_layer_scalefactor_1e4")
writeLines(capture.output(sessionInfo()), file.path(OUT, "signac_sessionInfo.txt"))
lg("ALL DONE")
