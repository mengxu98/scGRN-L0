

# 它会告诉你哪一行的代码消耗了多少时间、内存，释放多少内存，复制了多少向量
# library(devtools)
# devtools::install_github("hadley/lineprof")
# library(lineprof)
# source("D:/test/test.R")
# prof <- lineprof(test("D:/test/testcsv"))
# shine(prof)


# devtools::install_github("shmohammadi86/ACTIONetExperiment", force = TRUE)
# library(ACTIONetExperiment)
# BiocManager::install("batchelor")
# devtools::install_github("shmohammadi86/ACTIONet", ref = "R-devel")
# library(ACTIONet)
# devtools::install_github("shmohammadi86/SCINET")
# library(SCINET)
# devtools::install_github("netbiolab/scHumanNet")
# library(scHumanNet)

# library(SingleCellExperiment)

# data <- SingleCellExperiment(assays = list(logcounts = counts), colData = meta)
# ace <- reduce.ace(data)
# ace = run.ACTIONet(ace = ace, thread_no = 8)

# remotes::install_github("GfellerLab/SuperCell", force = TRUE, upgrade = FALSE)
library(SuperCell)
library(Seurat)
library(dplyr)

testData <- paste0("../scGRN-L0_data/Data1-Qian et al. Cell Research 2020/", "2096-Lungcancer_counts-export/LC_counts/")
counts <- Seurat::Read10X(testData)
head(counts)
meta <- read.table("../scGRN-L0_data/Data1-Qian et al. Cell Research 2020/2097-Lungcancer_metadata.csv", header = T, sep = ",")
head(meta)
# load single-cell (sc) count matrix and cell metadata
sc.counts <- counts
rownames(meta) <- meta$Cell
sc.meta <- meta
table(sc.meta$CellType)
proj.name <- "LUNG"
.color.cell.type <- c(
  "Alveolar" = "#E69F00", "B_cell" = "#026b02", "Cancer" = "#9e0000", "EC" = "#F0E442",
  "Epithelial" = "#CC79A7", "Erythroblast" = "#09fa96", "Fibroblast" = "#a31e8d",
  "Mast_cell" = "#0d0ab1", "Myeloid" = "#127032e7",
  "T_cell" = "#0299f0"
)
sc <- CreateSeuratObject(
  counts = sc.counts,
  project = proj.name,
  meta.data = sc.meta
)
sc <- NormalizeData(sc, verbose = FALSE)
sc <- FindVariableFeatures(
  sc,
  selection.method = "disp", # "vst" is default
  nfeatures = 1000,
  verbose = FALSE
)

hvg <- VariableFeatures(sc, verbose = FALSE)

# Plot variable features
plot1 <- VariableFeaturePlot(sc)
LabelPoints(plot = plot1, points = hvg[1:20], repel = TRUE)
sc <- ScaleData(sc, verbose = FALSE)
sc <- RunPCA(sc, verbose = FALSE)

# Plot PCA (2D representation of scRNA-seq data) colored by cell line
DimPlot(sc, reduction = "pca", group.by = "CellType", cols = .color.cell.type)

sc <- RunUMAP(sc, dims = 1:10)

# Plot UMAP (2D representation of scRNA-seq data) colored by cell line
DimPlot(sc, reduction = "umap", group.by = "CellType", cols = .color.cell.type)
sc <- FindNeighbors(sc, dims = 1:10)
sc <- FindClusters(sc, resolution = 0.05)

# As it is a toy example with well defined cell types (i.e., cell lines), unsupervised clustering fully recapitulates cell line annotation
table(sc@active.ident, sc$CellType)

# Set idents to cell lines (as clusters are the same as cell lines)
Idents(sc) <- "CellType"
levels(sc) <- sort(levels(sc))

# Compute upregulated genes in each cell line (versus other cells)
sc.all.markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "t")
# saveRDS(sc.all.markers, file = file.path(data.folder, "output", "sc_all_markers.Rds"))

# Load markers (backup)
# sc.all.markers <- readRDS(file = file.path(data.folder, "output", "sc_all_markers.Rds"))

# Top markers (select top markers of each cell line)
sc.top.markers <- sc.all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

sc.top.markers
VlnPlot(sc, features = sc.top.markers$gene[c(seq(1, 9, 2), seq(2, 10, 2))], ncol = 5, pt.size = 0.0, cols = .color.cell.type)

gamma <- 20 # Graining level

# Compute metacells using SuperCell package
MC <- SCimplify(
  X = GetAssayData(sc), # single-cell log-normalized gene expression data
  genes.use = hvg,
  gamma = gamma,
  n.pc = 10
)

# Compute gene expression of metacells by simply averaging gene expression within each metacell
MC.ge <- supercell_GE(
  ge = GetAssayData(sc),
  groups = MC$membership
)

# Alternatively, counts can be averaged (summed up) followed by a lognormalization step (this approach is used in the MetaCell and SEACell algorithms)
if (0) {
  MC.counts <- supercell_GE(
    ge = GetAssayData(sc, slot = "counts"),
    mode = "sum", # summing counts instead of the default averaging
    groups = MC$membership
  )
  MC.ge <- Seurat::LogNormalize(MC.counts, verbose = FALSE)
}

# Annotate metacells to cells line
MC$CellType <- supercell_assign(
  cluster = sc.meta$CellType, # single-cell assignment to cell lines
  supercell_membership = MC$membership, # single-cell assignment to metacells
  method = "absolute" # available methods are c("jaccard", "relative", "absolute"), function's help() for explanation
)

# Compute purity of metacells as :
#  * a proportion of the most abundant cell type withing metacells (`method = `"max_proportion)
#  * an entropy of cell type within metacells (`method = "entropy"`)
method_purity <- c("max_proportion", "entropy")[1]
MC$purity <- supercell_purity(
  clusters = sc.meta$CellType,
  supercell_membership = MC$membership,
  method = method_purity
)

# Metacell purity distribution
summary(MC$purity)

hist(MC$purity, main = paste0("Purity of metacells \nin terms of cell line composition (", method_purity, ")"))

supercell_plot(
  MC$graph.supercells,
  group = MC$CellType,
  color.use = .color.cell.type,
  seed = 1,
  alpha = -pi / 2,
  main = "Metacells colored by cell line assignment"
)
# supercell_plot(
#   MC$graph.singlecell,
#   group = sc.meta$CellType,
#   color.use = .color.cell.type,
#   do.frames = FALSE,
#   lay.method = "components",
#   seed = 1,
#   alpha = -pi / 2,
#   main = "Single cells colored by cell line assignment"
# )

MC.seurat <- supercell_2_Seurat(
  SC.GE = MC.ge,
  SC = MC,
  fields = c("CellType", "purity"),
  var.genes = MC$genes.use,
  N.comp = 10
)

Idents(MC.seurat) <- "CellType"
MC.seurat <- RunUMAP(MC.seurat, dims = 1:10)

DimPlot(MC.seurat, cols = .color.cell.type, reduction = "umap")
MC.seurat <- FindClusters(MC.seurat, resolution = 0.5)

DimPlot(MC.seurat, reduction = "umap")

# Set idents to cell lines (as clusters are the same as cell lines)
Idents(MC.seurat) <- "CellType"
levels(MC.seurat) <- sort(levels(Idents(MC.seurat)))

# Compute upregulated genes in each cell line (versus other cells)
MC.seurat.all.markers <- FindAllMarkers(
  MC.seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "t"
)

# saveRDS(MC.seurat.all.markers, file = file.path(data.folder, "output", paste0("MC_gamma_", gamma, "_all_markers_seurat.Rds")))
# MC.seurat.all.markers <- readRDS(file = file.path(data.folder, "output", "MC_gamma_20_all_markers_seurat.Rds"))

# Top markers (select top markers of each cell line)
MC.seurat.top.markers <- MC.seurat.all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

MC.seurat.top.markers

genes.to.plot <- MC.seurat.top.markers$gene[c(seq(1, 9, 2), seq(2, 10, 2))]
VlnPlot(MC.seurat, features = genes.to.plot, ncol = 5, pt.size = 0.0, cols = .color.cell.type)

gene_x <- MC$genes.use[500:505] # 500
gene_y <- MC$genes.use[550:555] # 600

alpha <- 0.7

p.SC <- supercell_GeneGenePlot(GetAssayData(sc, slot = "data"), gene_x = gene_x, gene_y = gene_y, clusters = sc$CellType, color.use = .color.cell.type, sort.by.corr = F, alpha = alpha)
p.SC$p
p.MC <- supercell_GeneGenePlot(MC.ge, gene_x = gene_x, gene_y = gene_y, supercell_size = MC$supercell_size, clusters = MC$CellType, color.use = .color.cell.type, sort.by.corr = F, alpha = alpha)
p.MC$p

MC$PCA <- supercell_prcomp(
  Matrix::t(MC.ge),
  genes.use = MC$genes.use, # or a new set of HVG can be computed
  supercell_size = MC$supercell_size, # provide this parameter to run sample-weighted version of PCA,
  k = 10
)

MC$UMAP <- supercell_UMAP(
  SC = MC,
  PCA_name = "PCA",
  n_neighbors = 50 # large number to repel cells
)

supercell_DimPlot(
  MC,
  groups = MC$CellType,
  dim.name = "UMAP",
  title = paste0("UMAP of metacells colored by cell line assignment"),
  color.use = .color.cell.type
)
# compute distance among metacells
D <- dist(MC$PCA$x)

# cluster metacells
MC$clustering <- supercell_cluster(D = D, k = 10, supercell_size = MC$supercell_size)

# Plot clustering result
supercell_DimPlot(
  MC,
  groups = factor(MC$clustering$clustering),
  dim.name = "UMAP",
  title = paste0("UMAP of metacells colored by metacell clustering")
)

# Compute upregulated genes in each cell line (versus other cells)
MC.all.markers <- supercell_FindAllMarkers(
  ge = MC.ge,
  clusters = MC$CellType,
  supercell_size = MC$supercell_size,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# saveRDS(MC.all.markers, file = file.path(data.folder, "output", paste0("MC_gamma_", gamma, "_all_markers.Rds")))

# Load markers (backup)
# MC.all.markers <- readRDS(file = file.path(data.folder, "output", "paste0("MC_gamma_", gamma, "_all_markers.Rds")))

# Transform the output of `supercell_FindAllMarkers()` to be in the format of the `Seurat::FindAllMarkers()`
MC.all.markers.df <- data.frame()
for (cl in names(MC.all.markers)) {
  cur <- MC.all.markers[[cl]]
  cur$cluster <- cl
  cur$gene <- rownames(cur)
  cur$avg_log2FC <- cur$logFC
  MC.all.markers.df <- rbind(MC.all.markers.df, cur)
}

# Top markers (select top markers of each cell line)
MC.top.markers <- MC.all.markers.df %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

supercell_VlnPlot(
  ge = MC.ge,
  supercell_size = MC$supercell_size,
  clusters = MC$CellType,
  features = MC.top.markers$gene[c(seq(1, 9, 2), seq(2, 10, 2))],
  color.use = .color.cell.type,
  ncol = 5
)
