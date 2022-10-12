

if (T) {
  rm(list = ls())
  gc()
  source("step_function.R")
  # devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
  
  library(tidyverse)
  library(patchwork)
  library(scuttle)
  library(Seurat)
  library(scran)
  # BiocManager::install("celldex")
  # BiocManager::install("SingleR")
  
  library(SingleR)
  library(celldex)
  library(Seurat)
  library(pheatmap)
  load("/data/mengxu/data/SingleR_data/HumanPrimaryCellAtlas_hpca.se_human.RData")
  load("/data/mengxu/data/SingleR_data/BlueprintEncode_bpe.se_human.RData")
  
}

load("/data/mengxu/data/all/lung_L0_data.Rdata") #seu_obj_list
head(mat_com[1:3,1:3])
mat_com <- t(mat_com)

message("[", Sys.time(), "] -----: data normalization!")

seu_obj_data <- CreateSeuratObject(
  counts = mat_com,
  min.features = 200,
  min.cells = 3
)
rm(mat_com)
gc()
### QC
seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^MT-", col.name = "pMT")
seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^HBA|^HBB", col.name = "pHB")
seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^RPS|^RPL", col.name = "pRP")
dim(seu_obj_data)

# Pre-process Seurat object (standard)
# seu_obj_data <- NormalizeData(seu_obj_data)
# seu_obj_data <- FindVariableFeatures(seu_obj_data)
# seu_obj_data <- ScaleData(seu_obj_data)

# Pre-process Seurat object (sctransform)
seu_obj_data <- SCTransform(seu_obj_data)

# ElbowPlot(seu_obj_data)
pc <- pc_num(seu_obj_data) # Auto infer 'pc' value
pc.num <- 1:pc

seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
seu_obj_data <- RunUMAP(seu_obj_data, dims = pc.num)
seu_obj_data <- FindNeighbors(seu_obj_data, dims = pc.num) %>% FindClusters(resolution = 0.3)

DimPlot(
  seu_obj_data,
  # group.by = "orig.ident",
  # label = T,
  # repel = T
  # pt.size = 0.2
) + theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
  theme_bw()
#+NoLegend()

library(harmony)
scRNA_harmony <- RunHarmony(seu_obj_data,
                            group.by.vars = "orig.ident",
                            assay.use = "SCT",
                            max.iter.harmony = 20
)
ElbowPlot(scRNA_harmony, ndims = 50)
pc <- pc_num(scRNA_harmony)
pc.num <- 1:pc
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num) # %>% RunTSNE(reduction="harmony", dims=pc.num)
scRNA_harmony <- FindNeighbors(scRNA_harmony, dims = pc.num) %>% FindClusters(resolution = 1)
DimPlot(scRNA_harmony,
        reduction = "umap",
        group.by = "orig.ident",
        label = T,
        cols = colP,
        repel = T
        # pt.size = 0.2
) + theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + # NoLegend() +
  scale_colour_viridis_d(option = "inferno")

DimPlot(seu_obj_data,
        reduction = "umap",
        group.by = "orig.ident",
        label = T,
        cols = colP,
        repel = T
        # pt.size = 0.2
) + theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + # NoLegend() +
  scale_colour_viridis_d(option = "inferno")

#-----------------------------------------------------------------------------------#
scRNA_harmony <- annotation_celltype(scRNA_harmony, method = "celltypist") # method = "celltypist" or "singleR"
table(scRNA_harmony$predicted_labels)

scRNA_harmony@meta.data$labels <- scRNA_harmony$predicted_labels
print(DimPlot(seu_obj_data, group.by = c("orig.ident"), reduction = "umap"))
print(DimPlot(seu_obj_data, group.by = c("labels"), reduction = "umap"))

print(DimPlot(scRNA_harmony, group.by = c("orig.ident"), reduction = "umap"))
print(DimPlot(scRNA_harmony, group.by = c("labels"), reduction = "umap"))

if (length(which(scRNA_harmony$labels == c("B_mature","B_naive","B_plasma","B_plasmablast"))) == 0) {
  print("No B cell !!!")
} else {
  seu_obj_data_Plasma <- scRNA_harmony[, (scRNA_harmony$labels == c("B_mature","B_naive","B_plasma","B_plasmablast"))]
}

obj_cells <- c("B_mature","B_naive","B_plasma","B_plasmablast")
seu_obj_data_obj_cells <- list()
for (i in 1:length(obj_cells)) {
  obj_cell <- obj_cells[i]
  seu_obj_data_obj_cell <- scRNA_harmony[, (scRNA_harmony$labels == obj_cell)]
  seu_obj_data_obj_cells[[i]] <- seu_obj_data_obj_cell
}
seu_obj_data_B <- merge(seu_obj_data_obj_cells[[1]],seu_obj_data_obj_cells[2:length(seu_obj_data_obj_cells)])

dim(seu_obj_data_B)
seu_obj_data_B <- SCTransform(seu_obj_data_B)
seu_obj_data_B <- RunPCA(seu_obj_data_B, verbose = T)
pc <- pc_num(seu_obj_data_B)
pc.num <- 1:pc

seu_obj_data_B <- RunUMAP(seu_obj_data_B, dims = pc.num) # %>% RunTSNE(reduction="harmony", dims=pc.num)
seu_obj_data_B <- FindNeighbors(seu_obj_data_B, dims = pc.num) %>% FindClusters(resolution = 1)
DimPlot(
  seu_obj_data_B,
  group.by = "labels",
  label = T,
  repel = T,
  pt.size = 0.2
) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
  theme_bw()

DimPlot(
  seu_obj_data_B,
  group.by = "orig.ident",
  label = T,
  repel = T,
  pt.size = 0.2
) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
  theme_bw()


scRNA_harmony_B <- RunHarmony(seu_obj_data_B,
                            group.by.vars = "orig.ident",
                            assay.use = "SCT",
                            max.iter.harmony = 20
)

scRNA_harmony_B <- RunUMAP(scRNA_harmony_B, reduction = "harmony", dims = pc.num) # %>% RunTSNE(reduction="harmony", dims=pc.num)
scRNA_harmony_B <- FindNeighbors(scRNA_harmony_B, dims = pc.num) %>% FindClusters(resolution = 1)
DimPlot(scRNA_harmony_B,
        reduction = "umap",
        group.by = "orig.ident",
        label = T,
        cols = colP,
        repel = T
        # pt.size = 0.2
) + theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + # NoLegend() +
  scale_colour_viridis_d(option = "inferno")

DimPlot(scRNA_harmony_B,
        reduction = "umap",
        group.by = "labels",
        label = T,
        # cols = colP,
        repel = T
        # pt.size = 0.2
)

DimPlot(seu_obj_data_B,
        reduction = "umap",
        group.by = "orig.ident",
        label = T,
        cols = colP,
        repel = T
        # pt.size = 0.2
) + theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + # NoLegend() +
  scale_colour_viridis_d(option = "inferno")








mat_com <- seu_obj_data_B@assays$RNA@counts

seu_obj_data <- CreateSeuratObject(
  counts = mat_com,
  min.features = 200,
  min.cells = 3
)
rm(mat_com)
gc()
### QC
seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^MT-", col.name = "pMT")
seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^HBA|^HBB", col.name = "pHB")
seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^RPS|^RPL", col.name = "pRP")
dim(seu_obj_data)

# Pre-process Seurat object (standard)
# seu_obj_data <- NormalizeData(seu_obj_data)
# seu_obj_data <- FindVariableFeatures(seu_obj_data)
# seu_obj_data <- ScaleData(seu_obj_data)

# Pre-process Seurat object (sctransform)
seu_obj_data <- SCTransform(seu_obj_data)

seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
pc <- pc_num(seu_obj_data)
pc.num <- 1:pc

seu_obj_data <- RunUMAP(seu_obj_data, dims = pc.num)
seu_obj_data <- FindNeighbors(seu_obj_data, dims = pc.num) %>% FindClusters(resolution = 0.3)

DimPlot(
  seu_obj_data,
  # group.by = "orig.ident",
  # label = T,
  # repel = T
  # pt.size = 0.2
) + theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
  theme_bw()
seu_obj_data <- annotation_celltype(seu_obj_data, method = "celltypist") # method = "celltypist" or "singleR"
table(seu_obj_data$predicted_labels)

seu_obj_data@meta.data$labels <- seu_obj_data$predicted_labels
print(DimPlot(seu_obj_data, group.by = c("orig.ident"), reduction = "umap"))
print(DimPlot(seu_obj_data, group.by = c("labels"), reduction = "umap"))

saveRDS(seu_obj_data_Plasma_list, "/data/mengxu/data/all/lung_B_Plasma.rds")
