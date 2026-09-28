# Export per-cell labels from Pal's Fig 4B Seurat object.
#
# The object is SeuratObject_NormB1Total.rds from Chen and Smyth's figshare
# deposit (doi:10.6084/m9.figshare.17058077), built with Seurat 3.1.1.
# This script writes one row per cell. It does not fit the proportion model.
#
# Usage:
#   Rscript scripts/export_pal_norm_b1_labels.R path/to/SeuratObject_NormB1Total.rds path/to/pal_norm_b1_labels.csv

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript scripts/export_pal_norm_b1_labels.R <object.rds> <labels.csv>")
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat is not installed. The figshare object cannot be read without it.")
}
suppressPackageStartupMessages(library(Seurat))

object <- readRDS(args[[1]])
meta <- object@meta.data
if (!"seurat_clusters" %in% colnames(meta)) {
  stop("The object has no seurat_clusters column. Fig 4B clusters are required.")
}
if (!"orig.ident" %in% colnames(meta)) {
  stop("The object has no orig.ident column. Sample names are required.")
}

cell_type <- NULL
cell_type_source <- "ident"
for (name in c("cell_type", "celltype", "CellType", "annotated_type")) {
  if (name %in% colnames(meta)) {
    cell_type <- as.character(meta[[name]])
    cell_type_source <- name
    break
  }
}
if (is.null(cell_type)) {
  cell_type <- as.character(Idents(object))
}

exported <- data.frame(
  barcode = rownames(meta),
  sample = as.character(meta$orig.ident),
  cluster = as.character(meta$seurat_clusters),
  cell_type = cell_type,
  cell_type_source = cell_type_source,
  stringsAsFactors = FALSE
)
write.csv(exported, args[[2]], row.names = FALSE)
message("Wrote ", nrow(exported), " cells to ", args[[2]])
message("Cell types came from ", cell_type_source, ".")
