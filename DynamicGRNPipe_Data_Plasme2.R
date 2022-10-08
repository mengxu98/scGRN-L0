

if (T) {
  rm(list = ls())
  gc()
  source("step_function.R")
  # devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
  library(DoubletFinder)
  library(tidyverse)
  library(scDblFinder)
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
  
  # Python
  library(reticulate)
  pandas <- import("pandas")
  numpy <- import("numpy")
  scanpy <- import("scanpy")
  celltypist <- import("celltypist")
}

nFeature_lower <- 200
nFeature_upper <- 10000
nCount_lower <- 100
nCount_upper <- 150000
pMT_lower <- 0
pMT_upper <- 20

load("/data/mengxu/data/all/lung_L0_data.Rdata") #seu_obj_list
head(mat_com[1:3,1:3])
mat_com <- t(mat_com)

seu_obj_data_Plasma_list <- list()
for (i in 1:length(seu_obj_list)) {
  message("[", Sys.time(), "] -----: data normalization!")
  seu_obj_data <- seu_obj_list[[i]]
  
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
  # seu_obj_data <- subset(seu_obj_data,
  #     subset =
  #         nFeature_RNA > nFeature_lower &
  #             nFeature_RNA < nFeature_upper &
  #             nCount_RNA > nCount_lower &
  #             nCount_RNA < nCount_upper &
  #             pMT < pMT_upper
  # )
  
  # Pre-process Seurat object (standard)
  # seu_obj_data <- NormalizeData(seu_obj_data)
  # seu_obj_data <- FindVariableFeatures(seu_obj_data)
  # seu_obj_data <- ScaleData(seu_obj_data)
  
  # Pre-process Seurat object (sctransform)
  seu_obj_data <- SCTransform(seu_obj_data)
  # ElbowPlot(seu_obj_data)
  pc.num <- 1:30
  
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
  
  #-----------------------------------------------------------------------------------#
  meta <- seu_obj_data@meta.data
  seu_obj_data <- annotation_celltype(seu_obj_data, method = "celltypist") # method = "celltypist" or "singleR"
  # seu_obj_data <- annotation_celltype(seu_obj_data, method = "singleR")
  table(seu_obj_data$predicted_labels)
  
  seu_obj_data@meta.data$labels <- seu_obj_data$predicted_labels
  print(DimPlot(seu_obj_data, group.by = c("orig.ident"), reduction = "umap"))
  print(DimPlot(seu_obj_data, group.by = c("labels"), reduction = "umap"))
  
  if (length(which(seu_obj_data$labels == c("B_mature","B_naive","B_plasma","B_plasmablast"))) == 0) {
    print("No B cell !!!")
  } else {
    seu_obj_data_Plasma <- seu_obj_data[, (seu_obj_data$labels == c("B_mature","B_naive","B_plasma","B_plasmablast"))]
    DimPlot(
      seu_obj_data_Plasma,
      # group.by = "orig.ident",
      group.by = "labels",
      label = T,
      repel = T,
      pt.size = 0.2
    ) + 
      theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
      theme_bw()
    # seu_obj_data_Plasma_list[[i]] <- seu_obj_data_Plasma
  }
  
}

saveRDS(seu_obj_data_Plasma_list, "/data/mengxu/data/all/lung_stage-III_list_Plasma.rds")
