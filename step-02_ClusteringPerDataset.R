

if (T) {
  rm(list = ls())
  gc()
  library(tidyverse)
  library(patchwork)
  library(scuttle)
  library(Seurat)
  library(scran)
  # BiocManager::install("celldex")
  # BiocManager::install("SingleR")
  library(celldex)
  library(pheatmap)
}

samples <- c(
  "PhilipBischoff2021.P030",
  "PhilipBischoff2021.P033",
  "DietherLambrechts2018.P06.1",
  "DietherLambrechts2018.P06.2",
  "DietherLambrechts2018.P06.3",
  "DietherLambrechts2018.P03.1",
  "DietherLambrechts2018.P03.2",
  "DietherLambrechts2018.P03.3",
  "RapolasZilionis2019.P2",
  "RapolasZilionis2019.P7",
  # 'JustinaXCaushi2021.MD01.010',
  # 'JustinaXCaushi2021.MD01.005',
  # 'JustinaXCaushi2021.NY016.021',
  # 'XinyiGuo2018.P0617',
  # 'XinyiGuo2018.P0619',
  # 'XinyiGuo2018.P0729',
  # 'XinyiGuo2018.P0913',
  "NayoungKim2020.P0028",
  "NayoungKim2020.P0031",
  "AndrewMLeader2021.P377.6",
  "AndrewMLeader2021.P403.11",
  "AndrewMLeader2021.P514.29",
  "AshleyMaynard2020.AZ.01",
  "AshleyMaynard2020.AZ.03",
  "AshleyMaynard2020.AZ.04",
  "AshleyMaynard2020.AZ.05",
  "AshleyMaynard2020.LT.S49",
  "AshleyMaynard2020.LT.S69",
  "AshleyMaynard2020.LT.S74"
)

seu_obj_list_filter <- readRDS("/data/mengxu/data/all/lung_stage-III_list_filter.rds")
seu_obj_data_B_cell_list <- list()
sample_no_B_cell <- list()
for (i in 1:length(seu_obj_list_filter)) {
  message("[", Sys.time(), "] -----: data normalization")

  seu_obj_data <- seu_obj_list_filter[[i]]

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
  meta <- seu_obj_data@meta.data # scRNA的meta文件，包含了seurat的聚类结果
  head(meta)
  seu_obj_data_SingleR <- seu_obj_data@assays$SCT@counts
  # seu_obj_data_SingleR <- GetAssayData(seu_obj_data, slot="data") ##获取标准化矩阵
  seu_obj_data_ann <- SingleR(test = seu_obj_data_SingleR, ref = hpca.se, labels = hpca.se$label.main) #
  seu_obj_data_ann

  # seurat 和 singleR的table表
  table(seu_obj_data_ann$labels, meta$seurat_clusters)

  seu_obj_data@meta.data$labels <- seu_obj_data_ann$labels
  print(DimPlot(seu_obj_data, group.by = c("seurat_clusters", "labels"), reduction = "umap"))

  ### 多个数据库注释
  # scRNA3 <- seu_obj_data_ann <- SingleR(test = seu_obj_data_SingleR, ref = list(BP=bpe.se, HPCA=hpca.se),
  #                                  labels = list(bpe.se$label.main, hpca.se$label.main))
  # table(seu_obj_data_ann$labels,meta$seurat_clusters)
  #
  # scRNA3$labels <- seu_obj_data_ann$labels
  #
  # print(DimPlot(seu_obj_data, group.by = c("seurat_clusters", "labels"),reduction = "umap"))

  ### 查看注释结果

  # 基于scores within cells
  print(plotScoreHeatmap(seu_obj_data_ann))

  # 基于 per-cell “deltas”
  plotDeltaDistribution(seu_obj_data_ann, ncol = 3)

  # 与cluster结果比较

  tab <- table(label = seu_obj_data_ann$labels, cluster = meta$seurat_clusters)

  pheatmap(log10(tab + 10))

  if (length(which(seu_obj_data$labels == "B_cell")) == 0) {
    print("No B cell !!!")
    sample_no_B_cell[i] <- samples[i]
  } else {
    seu_obj_data_B_cell <- seu_obj_data[, (seu_obj_data$labels == "B_cell")]
    seu_obj_data_B_cell_list[[i]] <- seu_obj_data_B_cell
  }
  # seu_obj_data_B_cell_list <- seu_obj_data_B_cell_list[-i]
}

saveRDS(seu_obj_data_B_cell_list, "/data/mengxu/data/all/lung_stage-III_list_B_cell.rds")
