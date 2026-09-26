### NormBRCA1_faithful.R
### Faithful reproduction of Pal et al. (2021) Fig 4A–C + Appendix S1A–C
### Ported from: https://github.com/yunshun/HumanBreast10X/blob/main/RCode/NormBRCA1.R
###
### Expected layout (created by the companion notebook):
###   <data_root>/
###     SampleStats.txt
###     180808_Homo_sapiens.gene_info.gz   (optional; skips official-symbol remap if absent)
###     N-0019-total/   (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
###     ... (all 12 samples)
###
### Usage:
###   Rscript NormBRCA1_faithful.R <data_root> <output_dir>
###
### Requires: Seurat, edgeR, limma, ggplot2, pheatmap  (scater optional for palettes)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript NormBRCA1_faithful.R <data_root> <output_dir>")
}
data_root <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)
setwd(out_dir)

message("data_root = ", data_root)
message("out_dir   = ", out_dir)

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
})

## Colour palettes (scater optional)
if (requireNamespace("scater", quietly = TRUE)) {
  col.pMedium <- scater:::.get_palette("tableau10medium")
  col.pDark <- scater:::.get_palette("tableau20")[2 * (1:10) - 1]
  col.pLight <- scater:::.get_palette("tableau20")[2 * (1:10)]
  col.p <- c(col.pDark, col.pLight)
} else {
  col.p <- rep(rainbow(20), length.out = 40)
}

#####################################
# Fig 4 A–C, Appendix S1 A–C
Samples <- c(
  "N-0019-total", "N-0233-total", "N-0092-total", "N-0230.17-total",
  "N-0093-total", "N-0123-total", "N-0064-total", "N-0169-total",
  "B1-0894", "B1-0033", "B1-0023", "B1-0090"
)
SamplesComb <- gsub("-", "_", Samples)
SamplesNormal <- SamplesComb[1:8]
SamplesPreB1 <- SamplesComb[9:12]

### Read 10X → DGEList (edgeR), matching author pipeline
DGE <- paste0("dge_", SamplesComb)
DD <- paste0("dd_", SamplesComb)

read_one <- function(sample_dir) {
  ## Prefer edgeR::read10X when available with DGEList=TRUE; else Seurat + DGEList
  if ("read10X" %in% getNamespaceExports("edgeR")) {
    return(edgeR::read10X(path = sample_dir, DGEList = TRUE))
  }
  counts <- Seurat::Read10X(data.dir = sample_dir)
  if (is.list(counts)) counts <- counts[[1]]
  y <- edgeR::DGEList(counts = counts)
  y$genes <- data.frame(Symbol = rownames(counts), stringsAsFactors = FALSE)
  y$samples$Barcode <- colnames(counts)
  y
}

for (i in seq_along(SamplesComb)) {
  sample_dir <- file.path(data_root, Samples[[i]])
  if (!dir.exists(sample_dir)) {
    stop("Missing sample directory: ", sample_dir)
  }
  message("Reading ", Samples[[i]])
  assign(DGE[[i]], read_one(sample_dir))
}

### Gene annotation (official symbols) — skip gracefully if NCBI file absent
gene_info <- file.path(data_root, "180808_Homo_sapiens.gene_info.gz")
for (i in seq_along(SamplesComb)) {
  y <- get(DGE[[i]])
  if (file.exists(gene_info)) {
    ann <- limma::alias2SymbolUsingNCBI(
      y$genes$Symbol,
      required.columns = c("GeneID", "Symbol"),
      gene.info.file = gene_info
    )
    Genes <- cbind(y$genes, Official = ann$Symbol, GeneID = ann$GeneID)
  } else {
    message("NCBI gene_info not found; using Cell Ranger symbols as Official")
    Genes <- cbind(y$genes, Official = y$genes$Symbol, GeneID = NA_character_)
  }
  y$genes <- Genes
  assign(DGE[[i]], y)
}

### QC metrics
for (i in seq_along(SamplesComb)) {
  y <- get(DGE[[i]])
  mito <- grep("^MT-", y$genes$Symbol)
  percent.mito <- colSums(y$counts[mito, , drop = FALSE]) / y$samples$lib.size
  y$samples <- cbind(y$samples, percent.mito = percent.mito, nGenes = colSums(y$counts != 0))
  assign(DGE[[i]], y)
}

### Cell filtering from SampleStats.txt (paper thresholds)
SampleStats <- read.delim(file.path(data_root, "SampleStats.txt"), stringsAsFactors = FALSE)
m <- match(Samples, SampleStats$SampleName)
if (anyNA(m)) stop("SampleStats.txt missing rows for: ", paste(Samples[is.na(m)], collapse = ", "))
SampleStats <- SampleStats[m, ]

mito_upper <- SampleStats$Mito
nGenes_lower <- SampleStats$GeneLower
nGenes_upper <- SampleStats$GeneUpper
lib_upper <- SampleStats$LibSize

for (i in seq_along(SamplesComb)) {
  y <- get(DGE[[i]])
  keep.mito <- y$samples$percent.mito < mito_upper[[i]]
  keep.nGenes <- y$samples$nGenes > nGenes_lower[[i]] & y$samples$nGenes < nGenes_upper[[i]]
  keep.nUMIs <- y$samples$lib.size < lib_upper[[i]]
  keep <- keep.mito & keep.nGenes & keep.nUMIs
  message(Samples[[i]], ": keeping ", sum(keep), " / ", length(keep), " cells")
  assign(DGE[[i]], y[, keep])
}

### Gene filtering (≥1% cells; valid Official; unique Official)
for (i in seq_along(SamplesComb)) {
  y <- get(DGE[[i]])
  colnames(y) <- paste(SamplesComb[[i]], y$samples$Barcode, sep = "_")
  o <- order(rowSums(y$counts), decreasing = TRUE)
  y <- y[o, ]
  keep1 <- rowSums(y$counts > 0) >= ncol(y) * 0.01
  keep2 <- !is.na(y$genes$Official)
  keep3 <- !duplicated(y$genes$Official)
  yall <- y[keep2 & keep3, , keep = FALSE]
  rownames(yall) <- yall$genes$Official
  assign(DD[[i]], yall)
  y <- y[keep1 & keep2 & keep3, , keep = FALSE]
  rownames(y) <- y$genes$Official
  assign(DGE[[i]], y)
}

### Per-sample Seurat objects
for (i in seq_along(SamplesComb)) {
  y <- get(DGE[[i]])
  so <- CreateSeuratObject(counts = y$counts, project = SamplesComb[[i]])
  so <- NormalizeData(so)
  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 1500)
  so <- ScaleData(so)
  so$group <- SamplesComb[[i]]
  assign(SamplesComb[[i]], so)
}

### Integration (Seurat CCA anchors) — paper parameters
CombSeurat <- lapply(SamplesComb, get)
names(CombSeurat) <- SamplesComb

dimUsed <- 30
message("FindIntegrationAnchors (dims=1:30, anchor.features=1000)...")
Anchors <- FindIntegrationAnchors(
  object.list = CombSeurat,
  dims = 1:dimUsed,
  anchor.features = 1000,
  scale = TRUE,
  k.anchor = 5,
  k.filter = 30,
  k.score = 20,
  max.features = 100
)
message("IntegrateData...")
NormB1Total <- IntegrateData(anchorset = Anchors, dims = 1:dimUsed, k.weight = 100)

DefaultAssay(NormB1Total) <- "integrated"
NormB1Total <- ScaleData(NormB1Total, verbose = FALSE)
NormB1Total <- RunPCA(NormB1Total, npcs = dimUsed, verbose = FALSE)
NormB1Total <- RunTSNE(NormB1Total, dims = 1:dimUsed, seed.use = 2018)
tSNE <- Embeddings(NormB1Total, "tsne")

Group <- factor(NormB1Total$group, levels = SamplesComb)
col <- col.p[as.integer(Group)]
plotOrd <- sample(ncol(NormB1Total))

pdf("Fig4A.pdf", height = 9, width = 9)
plot(tSNE[plotOrd, ], pch = 16, col = col[plotOrd], cex = 0.7,
     xlab = "tSNE-1", ylab = "tSNE-2", main = "t-SNE - By sample")
dev.off()

col.p2 <- c("darkblue", "gold2")
Group2 <- factor(c("Normal", "BRCA1")[as.integer(Group %in% SamplesPreB1) + 1L],
                 levels = c("Normal", "BRCA1"))

pdf("FigS1A-left.pdf", height = 7, width = 14)
par(mfrow = c(1, 2))
plot(tSNE[Group2 == "Normal", ], pch = 16, col = col.p2[[1]], cex = 0.7,
     xlab = "tSNE-1", ylab = "tSNE-2", main = "Normal")
plot(tSNE[Group2 == "BRCA1", ], pch = 16, col = col.p2[[2]], cex = 0.7,
     xlab = "tSNE-1", ylab = "tSNE-2", main = "BRCA1")
dev.off()

## Clustering — resolution 0.12 (paper)
resolution <- 0.12
NormB1Total <- FindNeighbors(NormB1Total, dims = 1:dimUsed, verbose = FALSE)
NormB1Total <- FindClusters(NormB1Total, resolution = resolution, verbose = FALSE)
Cluster <- as.integer(NormB1Total$seurat_clusters)
ncls <- length(table(Cluster))
message("All-cell clusters: ", ncls)

pdf("Fig4B.pdf", height = 9, width = 9)
plot(tSNE, pch = 16, col = col.p[Cluster], cex = 0.7,
     xlab = "tSNE-1", ylab = "tSNE-2", main = "t-SNE - By cluster")
dev.off()

cellNum <- table(Cluster, Group2)
cellProp <- t(t(cellNum) / colSums(cellNum))
pdf("FigS1A-right.pdf", height = 6, width = 9)
barplot(t(cellProp * 100), beside = TRUE, xlab = "Cluster", ylab = "Cell percentage",
        col = rep(col.p2, ncls), names.arg = seq_len(ncls), space = c(0, 0.3))
legend("topright", legend = levels(Group2), fill = col.p2)
dev.off()

saveRDS(NormB1Total, file = "SeuratObject_NormB1Total.rds")
save(SamplesComb, SamplesNormal, SamplesPreB1, Group, Group2, ncls, Cluster, tSNE,
     file = "NormB1Total_meta.RData")

####################################################################################################
### Micro-environment — remove epithelial clusters 2, 4, 5, 8 (author hard-coded IDs)
Sub <- !(Cluster %in% c(2, 4, 5, 8))
cellNamesSub <- colnames(NormB1Total)[Sub]
message("Non-epithelial cells: ", sum(Sub))

DGESub <- paste0("dge_sub_", SamplesComb)
DDSub <- paste0("dd_sub_", SamplesComb)
for (i in seq_along(SamplesComb)) {
  d <- get(DD[[i]])
  d <- d[, colnames(d) %in% cellNamesSub]
  keep1 <- rowSums(d$counts > 0) >= ncol(d) * 0.01
  assign(DDSub[[i]], d)
  assign(DGESub[[i]], d[keep1, , keep = FALSE])
}

SamplesCombSub <- paste0(SamplesComb, "_Sub")
CombSeuratSub <- list()
for (i in seq_along(SamplesCombSub)) {
  d <- get(DGESub[[i]])
  so <- CreateSeuratObject(counts = d$counts, project = SamplesCombSub[[i]])
  so <- NormalizeData(so)
  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 1500)
  so <- ScaleData(so)
  so$group <- SamplesComb[[i]]
  CombSeuratSub[[i]] <- so
}
names(CombSeuratSub) <- SamplesCombSub

message("Re-integrating non-epithelial cells...")
AnchorsSub <- FindIntegrationAnchors(
  object.list = CombSeuratSub,
  dims = 1:dimUsed,
  anchor.features = 1000,
  scale = TRUE,
  k.anchor = 5,
  k.filter = 30,
  k.score = 20,
  max.features = 100
)
NormB1TotalSub <- IntegrateData(anchorset = AnchorsSub, dims = 1:dimUsed, k.weight = 100)

DefaultAssay(NormB1TotalSub) <- "integrated"
NormB1TotalSub <- ScaleData(NormB1TotalSub, verbose = FALSE)
NormB1TotalSub <- RunPCA(NormB1TotalSub, npcs = dimUsed, verbose = FALSE)
NormB1TotalSub <- RunTSNE(NormB1TotalSub, dims = 1:dimUsed, seed.use = 2018)
tSNE <- Embeddings(NormB1TotalSub, "tsne")

GroupSub <- factor(NormB1TotalSub$group, levels = SamplesComb)
GroupSub2 <- factor(c("Normal", "BRCA1")[as.integer(GroupSub %in% SamplesPreB1) + 1L],
                    levels = c("Normal", "BRCA1"))

## Clustering — resolution 0.08 (paper)
resolution <- 0.08
NormB1TotalSub <- FindNeighbors(NormB1TotalSub, dims = 1:dimUsed, verbose = FALSE)
NormB1TotalSub <- FindClusters(NormB1TotalSub, resolution = resolution, verbose = FALSE)
ClusterSub <- as.integer(NormB1TotalSub$seurat_clusters)
ncls <- length(table(ClusterSub))
message("Stromal/immune clusters: ", ncls)

pdf("Fig4C.pdf", height = 9, width = 9)
plot(tSNE, pch = 16, col = col.p[ClusterSub], cex = 0.7,
     xlab = "tSNE-1", ylab = "tSNE-2", main = "t-SNE - By cluster (stromal/immune)")
dev.off()

### Combine raw counts for pseudo-bulk (author pattern)
allGenesFilter <- unlist(lapply(DGESub, function(nm) rownames(get(nm))))
allGenesFilter <- names(table(allGenesFilter))[table(allGenesFilter) >= 2]
common_genes <- Reduce(intersect, lapply(DDSub, function(nm) rownames(get(nm))))
common_genes <- intersect(common_genes, allGenesFilter)

y <- get(DDSub[[1]])[common_genes, ]
for (i in 2:length(SamplesComb)) {
  y <- cbind(y, get(DDSub[[i]])[common_genes, ])
}
## sample id is the prefix before the first "_" of the cell barcode? Author used SamplesComb_
## Our barcodes are SamplesComb_Barcode — group = SamplesComb prefix up to last sample token.
y$samples$group <- vapply(
  strsplit(colnames(y), "_", fixed = TRUE),
  function(p) paste(p[-length(p)], collapse = "_"),
  character(1)
)
y$samples$group <- factor(y$samples$group, levels = SamplesComb)

### Pseudo-bulk DE (edgeR QL) — Appendix S1C
cl_map <- setNames(as.integer(ClusterSub), colnames(NormB1TotalSub))
sampClust <- paste(as.character(y$samples$group), cl_map[colnames(y)], sep = "_Clst")
stopifnot(!anyNA(cl_map[colnames(y)]))

counts2 <- t(rowsum(t(as.matrix(y$counts)), group = sampClust))
yComb <- DGEList(counts2)
yComb$samples$Patient <- gsub("_Clst.*$", "", colnames(yComb))
yComb$samples$Cluster <- as.numeric(gsub("^.*_Clst", "", colnames(yComb)))
yComb$samples$group <- yComb$samples$Cluster

N <- sort(unique(yComb$samples$Cluster))
yClstSub <- yComb[, yComb$samples$Cluster %in% N]
keep <- filterByExpr(yClstSub, min.count = 7, min.total.count = 15)
yClstSub <- yClstSub[keep, , keep = FALSE]
yClstSub <- calcNormFactors(yClstSub)

Population <- yClstSub$samples$Patient
Population[Population %in% SamplesNormal] <- "Normal"
Population[Population %in% SamplesPreB1] <- "BRCA1-Pre"
Population <- factor(Population, levels = c("Normal", "BRCA1-Pre"))
yClstSub$samples$Population <- Population

sel <- yClstSub$samples$lib.size > 4e4
yClstSub2 <- yClstSub[, sel]

Cls <- as.factor(yClstSub2$samples$Cluster)
Pat <- factor(yClstSub2$samples$Patient, levels = SamplesComb)
Pop <- yClstSub2$samples$Population
design <- model.matrix(~ Cls + Pat)

yClstSub2 <- estimateDisp(yClstSub2, design = design)
qfit2 <- glmQLFit(yClstSub2, design)

prior.count <- 1
zClstSub2 <- edgeR::cpm(yClstSub2, log = TRUE, prior.count = prior.count)

ncls <- nlevels(Cls)
contr <- rbind(
  matrix(1 / (1 - ncls), ncls, ncls),
  matrix(0, ncol(design) - ncls, ncls)
)
diag(contr) <- 1
contr[1, ] <- 0
rownames(contr) <- colnames(design)
colnames(contr) <- paste0("Cls", seq_len(ncls))

ctest <- lapply(seq_len(ncls), function(i) glmQLFTest(qfit2, contrast = contr[, i]))

top <- 15
pseudoMakers <- lapply(seq_len(ncls), function(i) {
  ord <- order(ctest[[i]]$table$PValue, decreasing = FALSE)
  upreg <- ctest[[i]]$table$logFC > 0
  rownames(yClstSub2)[ord[upreg][seq_len(min(top, sum(upreg)))]]
})
Markers <- unique(unlist(pseudoMakers))

annot2 <- data.frame(
  Cluster = paste0("Cluster ", Cls),
  Patient = Pat,
  Population = Pop,
  row.names = colnames(zClstSub2)
)
ann_colors2 <- list(
  Cluster = setNames(col.p[seq_len(ncls)], paste0("Cluster ", levels(Cls))),
  Patient = setNames(col.p[seq_along(SamplesComb)], SamplesComb),
  Population = setNames(col.p2[1:2], c("Normal", "BRCA1-Pre"))
)

pdf("FigS1C.pdf", height = 12, width = 9)
mat4 <- t(scale(t(zClstSub2[Markers, , drop = FALSE])))
pheatmap(
  mat4,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  border_color = NA,
  breaks = seq(-2, 2, length.out = 101),
  cluster_cols = TRUE,
  scale = "none",
  fontsize_row = 5,
  show_colnames = FALSE,
  treeheight_row = 70,
  treeheight_col = 70,
  clustering_method = "ward.D2",
  annotation_col = annot2,
  annotation_colors = ann_colors2
)
dev.off()

saveRDS(NormB1TotalSub, file = "SeuratObject_NormB1TotalSub.rds")
save(ClusterSub, GroupSub, GroupSub2, Markers, file = "NormB1TotalSub_meta.RData")

## Cluster abundance quasi-Poisson note (Appendix S1B): authors report P = 0.14
## Reproduce with NormTotal.R / NormEpi.R style GLM if desired; write contingency for inspection.
write.csv(as.matrix(cellNum), file = "cluster_by_condition_counts.csv")

message("Done. Figures written to ", out_dir)
