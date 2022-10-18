

if (T) {
  rm(list = ls())
  gc()
  # R
  library(Seurat)
  library(harmony)
  library(tidyverse)
  library(patchwork)
  library(tsne)
  library(cowplot)
  library(ggplot2)
  library(dplyr)
  library(sctransform)
  library(tidyr)
  library(magrittr)
  library(reshape2)
  library(readxl)
  library(progeny)
  library(readr)
  library(stringr)
  library(celldex)
  library(ggsci)
  library(ggpubr)
  # library(sscVis)
  library(data.table)
  library(R.utils)
  library(plyr)
  library(grid)
  library(ggrepel)
  library(viridis)
  library(SingleCellExperiment)
  # BiocManager::install('glmGamPoi')
  library(glmGamPoi) # Boost 'SCTransform' function
  source("Function.R")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
}

# Merge samples --------------------------------------------------
if (T) {
  stage <- c("normal", "1", "2", "3", "4")
  for (s in stage) {
    if (file.exists(paste0("/data/mengxu/data/all/lung_stage-", s, "_list_filter.Rdata")) == T) {
      load(paste0("/data/mengxu/data/all/lung_stage-", s, "_list_filter.Rdata")) # From 'step-01_doubletfinder&normalzation'
      seu_obj_filter <- merge_seu_obj(seu_obj_list_filter, samples, stage = s) # stage = "normal" or "1" or "2" or "3" or "4"
      save(seu_obj_filter, samples, file = paste0("/data/mengxu/data/all/lung_stage-", s, "_seu.Rdata"))
      rm(seu_obj_list_filter)
      rm(samples)
      gc()
    }
  }
}

# Merge all samples --------------------------------------------------
if (T) {
  stage <- c("normal", "1", "2", "3", "4")
  seu_obj_filter_list <- list()
  samples_list <- c()
  i <- 0
  for (s in stage) {
    i <- i + 1
    if (file.exists(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu.Rdata")) == T) {
      load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu.Rdata")) # From 'step-01_doubletfinder&normalzation'
      seu_obj_filter$stage <- paste0("stage-", s)
      seu_obj_filter_list[[i]] <- seu_obj_filter
      samples_list <- rbind.data.frame(samples_list, samples)
      rm(seu_obj_filter)
      rm(samples)
      gc()
    }
  }
  seu_obj_data <- merge(seu_obj_filter_list[[1]],
                        y = c(
                          seu_obj_filter_list[2:length(seu_obj_filter_list)]
                        ),
                        # add.cell.ids = samples_list,
                        project = "NSCLC"
  )
  dim(seu_obj_data)
  table(seu_obj_data$orig.ident)
  table(seu_obj_data$celltype)
  save(seu_obj_data, samples_list, file = paste0("/data/mengxu/data/all/lung_seu.Rdata"))
}

# harmony -----------------------------------------------------------------
if (T) {
  
  load(paste0("/data/mengxu/data/all/lung_seu.Rdata"))
  gc()
  seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^MT-", col.name = "pMT")
  seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^HBA|^HBB", col.name = "pHB")
  seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^RPS|^RPL", col.name = "pRP")
  dim(seu_obj_data)
  obj_cells <- c("B_mature", "B_naive", "B_plasma", "B_plasmablast")
  seu_obj_data_obj_cells <- list()
  for (i in 1:length(obj_cells)) {
    obj_cell <- obj_cells[i]
    seu_obj_data_obj_cell <- seu_obj_data[, (seu_obj_data$celltype == obj_cell)]
    seu_obj_data_obj_cells[[i]] <- seu_obj_data_obj_cell
  }
  seu_obj_data <- merge(seu_obj_data_obj_cells[[1]], seu_obj_data_obj_cells[2:length(seu_obj_data_obj_cells)])
  dim(seu_obj_data)
  
  seu_obj_data <- SCTransform(seu_obj_data,
                              method = "glmGamPoi")
  seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
  pc.num <- 1:pc_num(seu_obj_data)
  seu_obj_data <- RunUMAP(seu_obj_data, dims = pc.num)
  seu_obj_data <- FindNeighbors(seu_obj_data, dims = pc.num) %>% FindClusters(resolution = 1)
  p1 <- DimPlot(seu_obj_data,
                reduction = "umap",
                group.by = "orig.ident"
  ) +
    theme_bw() +
    NoLegend()
  p2 <- DimPlot(seu_obj_data,
                reduction = "umap",
                group.by = "stage"
  ) +
    theme_bw()
  p3 <- DimPlot(seu_obj_data,
                reduction = "umap",
                group.by = "platform"
  ) +
    theme_bw()
  p1 + p2 + p3
  ggsave2("Fig1.raw_umap.png",
          path = paste0("Results/"),
          width = 30, height = 9, units = "cm"
  )
  if (length(table(seu_obj_data$platform)) > 1) {
    # Pre-process Seurat object (sctransform)
    seu_obj_data <- SCTransform(seu_obj_data,
                                method = "glmGamPoi",
                                vars.to.regress = c("nCount_RNA", "pMT", "platform"),
                                # features = rownames(seu_obj_data),
                                conserve.memory = T
    )
  } else {
    seu_obj_data <- SCTransform(seu_obj_data,
                                method = "glmGamPoi",
                                vars.to.regress = c("nCount_RNA", "pMT"),
                                # features = rownames(seu_obj_data),
                                conserve.memory = T
    )
  }
  
  seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
  pc.num <- 1:pc_num(seu_obj_data)
  seu_obj_data <- RunUMAP(seu_obj_data, dims = pc.num)
  seu_obj_data <- FindNeighbors(seu_obj_data, dims = pc.num) %>% FindClusters(resolution = 1)
  p1 <- DimPlot(seu_obj_data,
                reduction = "umap",
                group.by = "orig.ident"
  ) +
    theme_bw() +
    NoLegend()
  p2 <- DimPlot(seu_obj_data,
                reduction = "umap",
                group.by = "stage"
  ) +
    theme_bw()
  p3 <- DimPlot(seu_obj_data,
                reduction = "umap",
                group.by = "platform"
  ) +
    theme_bw()
  p1 + p2 + p3
  ggsave2("Fig2.sct_umap.png",
          path = paste0("Results/"),
          width = 30, height = 9, units = "cm"
  )
  
  seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
  seu_obj_data <- CellCycleScoring(seu_obj_data,
                                   s.features = s.genes,
                                   g2m.features = g2m.genes,
                                   set.ident = TRUE
  )
  DimPlot(seu_obj_data)
  # ggsave2(paste0("SuppFig.2_stage-", s, "_CellCycleScoring_raw.png"),
  #   path = paste0("Results/stage-", s),
  #   width = 10, height = 8, units = "cm"
  # )
  
  seu_obj_data$CC.Difference <- seu_obj_data$S.Score - seu_obj_data$G2M.Score
  seu_obj_data <- SCTransform(seu_obj_data,
                              method = "glmGamPoi",
                              vars.to.regress = "CC.Difference",
                              # features = rownames(seu_obj_data),
                              conserve.memory = T
  )
  seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
  DimPlot(seu_obj_data)
  # ggsave2(paste0("SuppFig.2_stage-", s, "_CellCycleScoring_SCT.png"),
  #   path = paste0("Results/stage-", s),
  #   width = 10, height = 8, units = "cm"
  # )
  # 观察细胞周期基因的表达情况
  # RidgePlot(seu_obj_data, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
  # RidgePlot(seu_obj_data, features = c(s.genes[1:25], g2m.genes[1:25]), ncol = 8)
  # RidgePlot(seu_obj_data, features = c(s.genes[26:length(s.genes)], g2m.genes[26:length(g2m.genes)]), ncol = 8)
  if (F) {
    # Cell cycle scoring
    seu_obj_data <- RunPCA(seu_obj_data,
                           features = c(s.genes, g2m.genes)
    )
    DimPlot(seu_obj_data)
    seu_obj_data <- RunPCA(seu_obj_data,
                           features = VariableFeatures(seu_obj_data),
                           nfeatures.print = 10
    )
    # FeatureScatter(scRNA_harmony, "G2M.Score", "S.Score", group.by = "Phase", pt.size = .1) +
    #   coord_fixed(ratio = 1)
    DimPlot(seu_obj_data, reduction = "Phase")
    if (T) {
      seu_obj_data <- ScaleData(seu_obj_data,
                                vars.to.regress = c("S.Score", "G2M.Score"),
                                features = rownames(seu_obj_data)
      )
    } else {
      seu_obj_data$CC.Difference <- seu_obj_data$S.Score - seu_obj_data$G2M.Score
      seu_obj_data <- ScaleData(seu_obj_data,
                                vars.to.regress = "CC.Difference",
                                features = rownames(seu_obj_data)
      )
    }
  }
  # ElbowPlot(seu_obj_data)
  # save(seu_obj_data,
  #   file = paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_filter_SCT.Rdata")
  # )
  # load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_filter_SCT.Rdata")
  pc.num <- 1:pc_num(seu_obj_data)
  seu_obj_data <- seu_obj_data %>%
    RunUMAP(dims = pc.num) %>%
    FindNeighbors(dims = pc.num) %>%
    FindClusters(resolution = 0.5)
  
  # p1 <- DimPlot(seu_obj_data,
  #               reduction = "umap",
  #               group.by = "celltype"
  # ) + theme_bw() + NoLegend()
  # 
  # p2 <- DimPlot(seu_obj_data,
  #               reduction = "umap",
  #               group.by = "orig.ident"
  # ) + theme_bw() + NoLegend()
  # 
  # p3 <- DimPlot(seu_obj_data,
  #               reduction = "umap",
  #               group.by = "platform"
  # ) + theme_bw() + NoLegend()
  # 
  # p1 + p2 + p3
  # ggsave2("Fig1.raw_umap.png",
  #         path = paste0("Results/stage-", s),
  #         width = 27, height = 9, units = "cm"
  # )
  
  if (F) {
    Anchors <- FindIntegrationAnchors(object.list = dataset.list, dims = 1:30)
    ITG.sbj <- IntegrateData(anchorset = Anchors, dims = 1:30)
    DefaultAssay(ITG.sbj) <- "integrated"
  }
  scRNA_harmony <- RunHarmony(seu_obj_data,
                              group.by.vars = "orig.ident",
                              assay.use = "SCT",
                              # lambda = 1, # [0.5-2] The more smaller lambda value, the bigger integration efforts.
                              max.iter.harmony = 20
  )
  
  pc.num <- 1:pc_num(scRNA_harmony)
  scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num) %>%
    FindNeighbors(dims = pc.num) %>%
    FindClusters(resolution = 1)

  
  p1 <- DimPlot(scRNA_harmony,
                reduction = "umap",
                group.by = "orig.ident"
  ) + theme_bw() + NoLegend()
  
  p2 <- DimPlot(scRNA_harmony,
                reduction = "umap",
                group.by = "stage"
  ) + theme_bw() + NoLegend()
  
  p3 <- DimPlot(scRNA_harmony,
                reduction = "umap",
                group.by = "platform"
  ) + theme_bw() + NoLegend()
  p1 + p2 + p3
  ggsave2("Fig3.harmony_umap.png",
          path = paste0("Results/"),
          width = 30, height = 9, units = "cm"
  )
  
  #----------------------------------------------------------------------------#
  sce.merged <- as.SingleCellExperiment(scRNA_harmony)
  # LISI index
  lisi.pca <- lisi::compute_lisi(
    reducedDim(sce.merged, "PCA"),
    colData(sce.merged), c("orig.ident", "platform") #stage
  ) # "dataset","dataset.tech","ClusterID.pca"
  
  lisi.harmony <- lisi::compute_lisi(
    reducedDim(sce.merged, "HARMONY"),
    colData(sce.merged), c("orig.ident", "platform")
  ) # c("dataset","dataset.tech","ClusterID.harmony")
  
  lisi.pca.tb <- cbind(
    data.table(
      cellID = rownames(lisi.pca),
      rd = "Raw", method = "raw"
    ),
    lisi.pca[, c("orig.ident"), drop = F]
  )
  lisi.harmony.tb <- cbind(
    data.table(
      cellID = rownames(lisi.harmony),
      rd = "Harmony", method = "harmony"
    ),
    lisi.harmony[, c("orig.ident"), drop = F]
  )
  lisi.merge.tb <- rbind(lisi.pca.tb, lisi.harmony.tb)
  lisi.merge.tb[, .(mean(orig.ident)), by = "method"]
  lisi.merge.tb[, .(median(orig.ident)), by = "method"]
  
  p <- ggboxplot(lisi.merge.tb,
                 x = "method", y = "orig.ident",
                 fill = "method", alpha = 0.8
  ) +
    stat_compare_means(comparisons = list(c("raw", "harmony"))) +
    ylab("LISI") +
    theme(legend.position = "right") + scale_fill_tron()
  print(p)
  ggsave2(paste0("Results/Fig4.LISI.png"),
          width = 3.2, height = 4,
          dpi = 600
  )
  
  if (T) {
    scRNA_harmony <- annotation_celltype(scRNA_harmony, method = "celltypist") # method = "celltypist" or "singleR"
    levels(seu_obj_data$celltype)
    table(scRNA_harmony$celltype)
  }
  obj_cells <- c("B_mature", "B_naive", "B_plasma", "B_plasmablast")
  seu_obj_data_obj_cells <- list()
  for (i in 1:length(obj_cells)) {
    obj_cell <- obj_cells[i]
    seu_obj_data_obj_cell <- scRNA_harmony[, (scRNA_harmony$celltype == obj_cell)]
    seu_obj_data_obj_cells[[i]] <- seu_obj_data_obj_cell
  }
  seu_obj_data <- merge(seu_obj_data_obj_cells[[1]], seu_obj_data_obj_cells[2:length(seu_obj_data_obj_cells)])
  dim(seu_obj_data)
  
  
  if (F) {
    mainmarkers <- c(
      # Nature Medicine-Phenotype molding of stromal cells in the lung  tumor microenvironment
      # "ACAT2",
      # "CLDN18", # Alveolar
      # "CLDN5", # Endothelial
      # "CAPS", # Epithelial
      # "ALB", # Hepatocytes
      # "COL1A1", # Fibroblast
      "CD79A", # B cell
      # "LYZ", # Myeloid
      # "CD3D", # T cell
      # "EPCAM", # Cancer cell
      # 2022-Cancer cell-Intratumoral plasma cells predict outcomes to PD-L1 blockade in non-small cell lung cancer
      # Follicular B cells
      "BANK1",
      "CD83",
      "CD69",
      "SELL",
      # "LINC00926",
      # "MARCH1",
      # "FCER2",
      # "GAPT",
      # "HVCN1",
      #Germinal center B cells
      "AICDA",
      "HMGA1",
      "RGS13",
      # "GCSAM",
      # "LRMP",
      # "AC023590.1",
      # "SUSD3",
      #Plasma cell
      "MZB1",
      "DERL3",
      "XBP1",
      "IGHG2",
      "IGHGP",
      "IGHA2"
      # "SDC1",
      # "DERL3",
      # "JSRP1",
      # "TNFRSF17",
      # "SLAMF7",
      # "IGLV3-1",
      # "IGLV6-57",
      # "IGKV4-1",
      # "IGKV1-12",
      # "IGLC7",
      # "IGLL5"
    )
    library(ggplot2)
    seu_obj_data <- SCTransform(seu_obj_data)
    seu_obj_data <- RunPCA(seu_obj_data)
    seu_obj_data <- RunUMAP(seu_obj_data, dims = 1:15) %>%
      FindNeighbors(dims = 1:15) %>%
      FindClusters(resolution = 1)
    FeaturePlot(seu_obj_data,
                features = "CD14",
                reduction = "umap",
                coord.fixed = T,
                order = T,
                cols = viridis(10)
    ) +
      scale_color_viridis(discrete = F, option = "inferno")
    DotPlot(seu_obj_data, features = unique(mainmarkers), group.by = "seurat_clusters") +theme_bw()+
      RotatedAxis() +
      scale_x_discrete("") +
      scale_y_discrete("") +
      # coord_flip() +
      scale_color_viridis(discrete = F, option = "C")
    
    DotPlot(scRNA_harmony, features = unique(mainmarkers), group.by = "seurat_clusters") +theme_bw()+
      RotatedAxis() +
      scale_x_discrete("") +
      scale_y_discrete("") +
      # coord_flip() +
      scale_color_viridis(discrete = F, option = "C")
    
    ggsave2(paste0("Fig5.FeaturePlot_mainmarkers_", s, ".png"),
            path = paste0("Results/"),
            width = 20, height = 12, units = "cm"
    )
    
    for (i in seq_along(mainmarkers)) {
      if (mainmarkers[i] %in% rownames(scRNA_harmony)) {
        FeaturePlot(scRNA_harmony,
                    features = mainmarkers[i],
                    coord.fixed = T,
                    order = T,
                    cols = viridis(10)
        ) +
          scale_color_viridis(discrete = F, option = "inferno")
        ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".png"),
                path = paste0("Results/marker/"),
                width = 10, height = 10, units = "cm"
        )
      }
    }
  }
  
  new.cluster.ids <- c(
    "B", "B","Follicular B cells",
    "Follicular B cells",  "B", "B", "Plasma",
    "Plasma", "Plasma", "Plasma", "B", "B", "B",
    "Cancer", "Plasma", "Myeloid",
    "Plasma", "Myeloid", "B",
   "B", "Plasma", "Plasma", "B", "Follicular B cells"
  )
  
  
  scRNA_harmony@meta.data$new.cluster.ids <- NA
  
  for (i in 1:length(scRNA_harmony@meta.data$new.cluster.ids)) {
    # if (scRNA_harmony@meta.data$seurat_clusters[i]%in%c(5,8,12,13,19,20,22)) {
    #   print('1')
    #   scRNA_harmony@meta.data$new.cluster.ids[i] <- "Plasma_B_cell"
    # }else{
    #   scRNA_harmony@meta.data$new.cluster.ids[i] <- "UnDef"
    # }
    scRNA_harmony@meta.data$new.cluster.ids[i] <- new.cluster.ids[scRNA_harmony@meta.data$seurat_clusters[i]]
  }
  # 获取当前用的Idents
  Idents(object = scRNA_harmony)
  levels(scRNA_harmony)
  Idents(scRNA_harmony) <- "seurat_clusters"
  new.cluster.ids <- as.character(new.cluster.ids)
  names(new.cluster.ids) <- levels(scRNA_harmony)
  scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)
  
  scRNA_harmony@meta.data$new.cluster.ids
  
  scRNA_harmony <- RenameIdents(scRNA_harmony,
                                `0` = "B", `1` = "B", `2` = "Follicular B cells",
                                `3` = "Follicular B cells", `4` = "B", `5` = "B", `6` = "Plasma",
                                `7` = "Plasma", `8` = "Plasma", `9` = "Plasma",`10` = "B", `11` = "B",`12` = "B",
                                `13` = "Cancer", `14` = "Plasma", `15` = "Myeloid",
                                `16` = "Plasma", `17` = "Myeloid", `18` = "B",
                                `19` = "B", `20` = "Plasma", `21` = "Plasma",`22` = "B",`23` = "Follicular B cells"
  )
  
  DimPlot(scRNA_harmony, label = T) + NoLegend()+theme_bw()
  DimPlot(scRNA_harmony, group.by = "orig.ident") +theme_bw()+ NoLegend()
  DotPlot(scRNA_harmony, features = unique(mainmarkers), group.by = "new.cluster.ids") + RotatedAxis() +
    scale_x_discrete("") + scale_y_discrete("")
  
  print(DimPlot(scRNA_harmony, reduction = "umap", group.by = c("new.cluster.ids")))
  
  p <- FeaturePlot(scRNA_harmony, features = markers, ncol = 8)
  p
  
  obj_cells <- c("B", "Plasma", "Follicular B cells")
  seu_obj_data_obj_cells <- list()
  for (i in 1:length(obj_cells)) {
    obj_cell <- obj_cells[i]
    seu_obj_data_obj_cell <- scRNA_harmony[, (scRNA_harmony$new.cluster.ids == obj_cell)]
    seu_obj_data_obj_cells[[i]] <- seu_obj_data_obj_cell
  }
  seu_obj_data <- merge(seu_obj_data_obj_cells[[1]], seu_obj_data_obj_cells[2:length(seu_obj_data_obj_cells)])
  dim(seu_obj_data)
  table(seu_obj_data$new.cluster.ids)
  save(seu_obj_data, file="../scGRN-L0_data/seu_obj_data_B_samples.Rdata")
}

# Select samples and cells --------------------------------------------------
if (T) {
  load(paste0("/data/mengxu/data/all/lung_stage-", 4, "_seu.Rdata"))
  dim(seu_obj_filter)
  seu_obj_data <- seu_obj_filter
  seu_obj_data <- SCTransform(seu_obj_data)
  seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
  pc.num <- 1:pc_num(seu_obj_data)
  seu_obj_data <- RunUMAP(seu_obj_data, dims = pc.num)
  seu_obj_data <- FindNeighbors(seu_obj_data, dims = pc.num) %>% FindClusters(resolution = 1)
  p1 <- DimPlot(seu_obj_data,
                reduction = "umap",
                group.by = "orig.ident"
  ) +
    theme_bw() +
    NoLegend()
  p2 <- DimPlot(seu_obj_data,
                reduction = "umap",
                group.by = "celltype"
  ) +
    theme_bw() +
    NoLegend()
  p1 + p2
  
  obj_cells <- c("B_mature", "B_naive", "B_plasma", "B_plasmablast")
  seu_obj_data_obj_cells <- list()
  for (i in 1:length(obj_cells)) {
    obj_cell <- obj_cells[i]
    seu_obj_data_obj_cell <- seu_obj_data[, (seu_obj_data$celltype == obj_cell)]
    seu_obj_data_obj_cells[[i]] <- seu_obj_data_obj_cell
  }
  seu_obj_data_B <- merge(seu_obj_data_obj_cells[[1]], seu_obj_data_obj_cells[2:length(seu_obj_data_obj_cells)])
  
  dim(seu_obj_data_B)
  table(seu_obj_data_B$celltype)
  seu_obj_data_B <- SCTransform(seu_obj_data_B)
  # seu_obj_data <- NormalizeData(seu_obj_data)
  # seu_obj_data <- FindVariableFeatures(seu_obj_data)
  # seu_obj_data <- ScaleData(seu_obj_data)
  seu_obj_data_B <- RunPCA(seu_obj_data_B, verbose = T)
  pc <- pc_num(seu_obj_data_B)
  pc.num <- 1:pc
  
  seu_obj_data_B <- RunUMAP(seu_obj_data_B, dims = pc.num) # %>% RunTSNE(reduction="harmony", dims=pc.num)
  seu_obj_data_B <- FindNeighbors(seu_obj_data_B, dims = pc.num) %>% FindClusters(resolution = 1)
  DimPlot(
    seu_obj_data_B,
    group.by = "celltype",
    label = T,
    repel = T,
    pt.size = 0.2
  ) +
    # theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
    theme_bw()
  
  DimPlot(
    seu_obj_data_B,
    group.by = "orig.ident",
    label = T,
    repel = T,
    pt.size = 0.2
  ) +
    # theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
    theme_bw()
}

# harmony -----------------------------------------------------------------
if (T) {
  stage <- c("normal", "1", "2", "3", "4")
  for (s in stage) {
    load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu.Rdata"))
    seu_obj_data <- seu_obj_filter
    rm(seu_obj_filter)
    gc()
    seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^MT-", col.name = "pMT")
    seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^HBA|^HBB", col.name = "pHB")
    seu_obj_data <- PercentageFeatureSet(seu_obj_data, pattern = "^RPS|^RPL", col.name = "pRP")
    # Normalize
    if (F) {
      seu_obj_data1 <- NormalizeData(seu_obj_data, normalization.method = "LogNormalize", scale.factor = 10000)
      seu_obj_data1 <- ScaleData(seu_obj_data1,
                                 # vars.to.regress = c("nCount_RNA", "pMT"),
                                 features = rownames(seu_obj_data)
      )
    }
    
    if (length(table(seu_obj_data$platform)) > 1) {
      # Pre-process Seurat object (sctransform)
      seu_obj_data <- SCTransform(seu_obj_data,
                                  method = "glmGamPoi",
                                  vars.to.regress = c("nCount_RNA", "pMT", "platform"),
                                  # features = rownames(seu_obj_data),
                                  conserve.memory = T
      )
    } else {
      seu_obj_data <- SCTransform(seu_obj_data,
                                  method = "glmGamPoi",
                                  vars.to.regress = c("nCount_RNA", "pMT"),
                                  # features = rownames(seu_obj_data),
                                  conserve.memory = T
      )
    }
    
    seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
    seu_obj_data <- CellCycleScoring(seu_obj_data,
                                     s.features = s.genes,
                                     g2m.features = g2m.genes,
                                     set.ident = TRUE
    )
    DimPlot(seu_obj_data)
    # ggsave2(paste0("SuppFig.2_stage-", s, "_CellCycleScoring_raw.png"),
    #   path = paste0("Results/stage-", s),
    #   width = 10, height = 8, units = "cm"
    # )
    
    seu_obj_data$CC.Difference <- seu_obj_data$S.Score - seu_obj_data$G2M.Score
    seu_obj_data <- SCTransform(seu_obj_data,
                                method = "glmGamPoi",
                                vars.to.regress = "CC.Difference",
                                # features = rownames(seu_obj_data),
                                conserve.memory = T
    )
    seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
    DimPlot(seu_obj_data)
    # ggsave2(paste0("SuppFig.2_stage-", s, "_CellCycleScoring_SCT.png"),
    #   path = paste0("Results/stage-", s),
    #   width = 10, height = 8, units = "cm"
    # )
    # 观察细胞周期基因的表达情况
    # RidgePlot(seu_obj_data, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
    # RidgePlot(seu_obj_data, features = c(s.genes[1:25], g2m.genes[1:25]), ncol = 8)
    # RidgePlot(seu_obj_data, features = c(s.genes[26:length(s.genes)], g2m.genes[26:length(g2m.genes)]), ncol = 8)
    if (F) {
      # Cell cycle scoring
      seu_obj_data <- RunPCA(seu_obj_data,
                             features = c(s.genes, g2m.genes)
      )
      DimPlot(seu_obj_data)
      seu_obj_data <- RunPCA(seu_obj_data,
                             features = VariableFeatures(seu_obj_data),
                             nfeatures.print = 10
      )
      # FeatureScatter(scRNA_harmony, "G2M.Score", "S.Score", group.by = "Phase", pt.size = .1) +
      #   coord_fixed(ratio = 1)
      DimPlot(seu_obj_data, reduction = "Phase")
      if (T) {
        seu_obj_data <- ScaleData(seu_obj_data,
                                  vars.to.regress = c("S.Score", "G2M.Score"),
                                  features = rownames(seu_obj_data)
        )
      } else {
        seu_obj_data$CC.Difference <- seu_obj_data$S.Score - seu_obj_data$G2M.Score
        seu_obj_data <- ScaleData(seu_obj_data,
                                  vars.to.regress = "CC.Difference",
                                  features = rownames(seu_obj_data)
        )
      }
    }
    # ElbowPlot(seu_obj_data)
    # save(seu_obj_data,
    #   file = paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_filter_SCT.Rdata")
    # )
    # load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_filter_SCT.Rdata")
    pc.num <- 1:pc_num(seu_obj_data)
    seu_obj_data <- seu_obj_data %>%
      RunUMAP(dims = pc.num) %>%
      FindNeighbors(dims = pc.num) %>%
      FindClusters(resolution = 0.3)
    
    if (F) {
      seu_obj_data <- annotation_celltype(seu_obj_data, method = "celltypist") # method = "celltypist" or "singleR"
      levels(seu_obj_data$seurat_clusters)
    }
    
    p1 <- DimPlot(seu_obj_data,
                  reduction = "umap",
                  group.by = "celltype"
    ) + theme_bw() + NoLegend()
    
    p2 <- DimPlot(seu_obj_data,
                  reduction = "umap",
                  group.by = "orig.ident"
    ) + theme_bw() + NoLegend()
    
    p3 <- DimPlot(seu_obj_data,
                  reduction = "umap",
                  group.by = "platform"
    ) + theme_bw() + NoLegend()
    
    p1 + p2 + p3
    ggsave2("Fig1.raw_umap.png",
            path = paste0("Results/stage-", s),
            width = 27, height = 9, units = "cm"
    )
    
    if (F) {
      Anchors <- FindIntegrationAnchors(object.list = dataset.list, dims = 1:30)
      ITG.sbj <- IntegrateData(anchorset = Anchors, dims = 1:30)
      DefaultAssay(ITG.sbj) <- "integrated"
    }
    scRNA_harmony <- RunHarmony(seu_obj_data,
                                group.by.vars = "orig.ident",
                                assay.use = "SCT",
                                # lambda = 1, # [0.5-2] The more smaller lambda value, the bigger integration efforts.
                                max.iter.harmony = 20
    )
    
    pc.num <- 1:pc_num(scRNA_harmony)
    scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num) %>%
      FindNeighbors(dims = pc.num) %>%
      FindClusters(resolution = 1)
    save(scRNA_harmony,
         file = paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_harmony.Rdata")
    )
    # load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_filter_SCT.Rdata")
    
    p1 <- DimPlot(scRNA_harmony,
                  reduction = "umap",
                  group.by = "celltype"
    ) + theme_bw() + NoLegend()
    
    p2 <- DimPlot(scRNA_harmony,
                  reduction = "umap",
                  group.by = "orig.ident"
    ) + theme_bw() + NoLegend()
    
    p3 <- DimPlot(scRNA_harmony,
                  reduction = "umap",
                  group.by = "platform"
    ) + theme_bw() + NoLegend()
    p1 + p2 + p3
    ggsave2("Fig2.harmony_umap.png",
            path = paste0("Results/stage-", s),
            width = 27, height = 9, units = "cm"
    )
    
    scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = pc.num)
    DimPlot(scRNA_harmony,
            reduction = "umap",
            group.by = "orig.ident"
    ) + theme_bw() + NoLegend()
    DimPlot(scRNA_harmony,
            reduction = "tsne",
            group.by = "celltype"
    ) + theme_bw() + NoLegend()
    
    
    if (F) {
      
      Cellratio <- prop.table(table(Idents(scRNA_harmony), scRNA_harmony$orig.ident), margin = 2)#计算各组样本不同细胞群比例
      Cellratio
      #BM1        BM2        BM3        GM1        GM2        GM3
      #  Endothelial 0.27305737 0.32663989 0.28683967 0.40820981 0.59293194 0.54664650
      #  Fibroblast  0.20733479 0.18072289 0.24096386 0.37115165 0.20418848 0.14422592
      #  Immune      0.44299201 0.19410977 0.24976830 0.15393387 0.09751309 0.18406455
      #  Epithelial  0.02505447 0.08299866 0.13253012 0.03534778 0.05366492 0.05698437
      #  Other       0.05156137 0.21552878 0.08989805 0.03135690 0.05170157 0.06807867
      Cellratio <- as.data.frame(Cellratio)
      colourCount = length(unique(Cellratio$Var1))
      library(ggplot2)
      ggplot(Cellratio) + 
        geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
        theme_classic() +
        labs(x='Sample',y = 'Ratio')+
        coord_flip()+
        theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
      
      table(scRNA_harmony$orig.ident)#查看各组细胞数
      prop.table(table(Idents(scRNA_harmony)))
      table(Idents(scRNA_harmony), scRNA_harmony$orig.ident)#各组不同细胞群细胞数 colnames(scRNA_harmony)
      Cellratio <- prop.table(table(Idents(scRNA_harmony), scRNA_harmony$orig.ident), margin = 2)#计算各组样本不同细胞群比例
      Cellratio <- data.frame(Cellratio)
      library(reshape2)
      cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
      rownames(cellper) <- cellper[,1]
      cellper <- cellper[,-1]
      
      ###添加分组信息
      sample <- scRNA_harmony$orig.ident
      group <-scRNA_harmony$stage
      samples <- data.frame(sample, group)#创建数据框
      
      rownames(samples)=samples$sample
      
      # cellper$sample <- samples[rownames(cellper),'sample']#R添加列
      cellper$sample <- rownames(cellper)#R添加列
      # cellper$group <- samples[rownames(cellper),'group']#R添加列
      cellper$group <- samples[rownames(cellper),'group']#R添加列
      
      ###作图展示
      pplist = list()
      sce_groups = c("B","Follicular B cells","Plasma","Cancer","Myeloid")
      library(ggplot2)
      library(dplyr)
      library(ggpubr)
      for(group_ in sce_groups){
        cellper_  = cellper %>% select(one_of(c('sample','group',group_)))#选择一组数据
        colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
        cellper_$percent = as.numeric(cellper_$percent)#数值型数据
        cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                            lower = quantile(percent, 0.25),
                                                            mean = mean(percent),
                                                            median = median(percent))#上下分位数
        print(group_)
        print(cellper_$median)
        
        pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
          geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
          stat_summary(fun=mean, geom="point", color="grey60") +
          theme_cowplot() +
          theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
          labs(title = group_,y='Percentage') +
          geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
        
        ###组间t检验分析
        labely = max(cellper_$percent)
        compare_means(percent ~ group,  data = cellper_)
        my_comparisons <- list( c("stage-1", "stage-2") )
        pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
        pplist[[group_]] = pp1
      }
      
      library(cowplot)
      plot_grid(pplist[['B']],
                pplist[['Follicular B cells']],
                pplist[['Plasma']],
                pplist[['Cancer']],
                pplist[['Myeloid']])
      
    }
    
    #----------------------------------------------------------------------------#
    sce.merged <- as.SingleCellExperiment(scRNA_harmony)
    # LISI index
    lisi.pca <- lisi::compute_lisi(
      reducedDim(sce.merged, "PCA"),
      colData(sce.merged), c("orig.ident", "platform")
    ) # "dataset","dataset.tech","ClusterID.pca"
    
    lisi.harmony <- lisi::compute_lisi(
      reducedDim(sce.merged, "HARMONY"),
      colData(sce.merged), c("orig.ident", "platform")
    ) # c("dataset","dataset.tech","ClusterID.harmony")
    
    lisi.pca.tb <- cbind(
      data.table(
        cellID = rownames(lisi.pca),
        rd = "Raw", method = "raw"
      ),
      lisi.pca[, c("orig.ident"), drop = F]
    )
    lisi.harmony.tb <- cbind(
      data.table(
        cellID = rownames(lisi.harmony),
        rd = "Harmony", method = "harmony"
      ),
      lisi.harmony[, c("orig.ident"), drop = F]
    )
    lisi.merge.tb <- rbind(lisi.pca.tb, lisi.harmony.tb)
    lisi.merge.tb[, .(mean(orig.ident)), by = "method"]
    lisi.merge.tb[, .(median(orig.ident)), by = "method"]
    
    p <- ggboxplot(lisi.merge.tb,
                   x = "method", y = "orig.ident",
                   fill = "method", alpha = 0.8
    ) +
      stat_compare_means(comparisons = list(c("raw", "harmony"))) +
      ylab("LISI") +
      theme(legend.position = "right") + scale_fill_tron()
    print(p)
    ggsave2(paste0("Results/stage-", s, "/Fig3.LISI.png"),
            width = 3.2, height = 4
    )
    
    # --------------------------------------------------
    if (F) {
      mainmarkers <- c(
        # Nature Medicine-Phenotype molding of stromal cells in the lung  tumor microenvironment
        # "ACAT2",
        "CLDN18", # Alveolar
        "CLDN5", # Endothelial
        "CAPS", # Epithelial
        "ALB", # Hepatocytes
        "COL1A1", # Fibroblast
        "CD79A", # B cell
        "LYZ", # Myeloid
        "CD3D", # T cell
        "EPCAM" # , # Cancer cell
        # 2022-Cancer cell-Intratumoral plasma cells predict outcomes to PD-L1 blockade in non-small cell lung cancer
        # Follicular B cells
        # "BANK1",
        # "LINC00926",
        # "MARCH1",
        # "FCER2",
        # "GAPT",
        # "HVCN1",
        # #Germinal center B cells
        # "AICDA",
        # "GCSAM",
        # "LRMP",
        # "AC023590.1",
        # "SUSD3",
        # #Plasma cell
        # "MZB1",
        # "DERL3",
        # "JSRP1",
        # "TNFRSF17",
        # "SLAMF7",
        # "IGHG2",
        # "IGHGP",
        # "IGLV3-1",
        # "IGLV6-57",
        # "IGHA2",
        # "IGKV4-1",
        # "IGKV1-12",
        # "IGLC7",
        # "IGLL5"
      )
      
      DotPlot(scRNA_harmony, features = unique(mainmarkers), group.by = "seurat_clusters") +
        RotatedAxis() +
        scale_x_discrete("") +
        scale_y_discrete("") +
        # coord_flip() +
        scale_color_viridis(discrete = F, option = "C")
      
      ggsave2(paste0("SuppFig5.FeaturePlot_mainmarkers_", s, ".pdf"),
              path = paste0("Results/stage-", s),
              width = 15, height = 20, units = "cm"
      )
      
      for (i in seq_along(mainmarkers)) {
        if (mainmarkers[i] %in% rownames(scRNA_harmony)) {
          FeaturePlot(scRNA_harmony,
                      features = mainmarkers[i],
                      coord.fixed = T,
                      order = T,
                      cols = viridis(10)
          ) +
            scale_color_viridis(discrete = F, option = "inferno")
          ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".png"),
                  path = paste0("Results/stage-", s, "/marker/"),
                  width = 10, height = 10, units = "cm"
          )
        }
      }
    }
  }
}

# annotation --------------------------------------------------------------
if (F) {
  # VlnPlot(scRNA_harmony,marker)
  # VlnPlot(scRNA_harmony,'CD79A')
  Follicular_B_cells <- c(
    "RP5-887A10.1", "AL928768.3",
    "BANK1", "LINC00926", "MARCH1", "FCER2"
  )
  
  Plasma_B_cells <- c(
    "IGHG3", "IGHG1", "IGHGP", "IGHG2",
    "IGHG4",
    "SLAMF7", "DERL3", "JSRP1"
  )
  
  Germinal_center_B_cells <- c("GCSAM", "LRMP", "AICDA", "SUSD3")
  
  MALT_B_cells <- c("IGHA1", "IGHA2")
  
  Mast_cells <- c(
    "TPSAB1", "TPSB2", "CPA3", "HPGDS", "CLU", "AREG", "MS4A2", "RGS13",
    "VWA5A", "LAPTM4A", "C1orf186", "SLC18A2", "LTC4S", "KIT", "HDC", "MAOB",
    "RGS1", "RP11-354E11.2", "SAMSN1", "RGS2", "SLC26A2", "PTGS1", "NSMCE1"
  )
  MALT_B_cells2 <- c("JCHAIN", "MZB1", "IGHM", "IGLL5", "IGLC7", "IGHD")
  Plasmacytoid_dendritic_cells <- c(
    "GZMB", "LILRA4", "CLIC3", "SMPD3", "PLD4", "LILRB4", "IRF7",
    "IL3RA", "CXCR3", "MAP1A", "PLAC8", "PTCRA", "JAML", "UGCG",
    "LAMP5", "SCT", "PPP1R14B", "IRF4", "SEC61B", "ITM2C", "C9orf142",
    "CLN8", "RNASE6", "RASD1", "IRF8", "SPIB", "TCL1A", "GPR183"
  )
  Erythroblasts <- c("HBB", "HBA2", "HBA1", "HBD", "SNCA", "SLC25A37", "ALAS2", "DCAF12", "SLC25A39")
  
  markers <- c(
    Follicular_B_cells, Plasma_B_cells, Mast_cells,
    MALT_B_cells, MALT_B_cells2, Plasmacytoid_dendritic_cells, Erythroblasts
  )
  markers <- c(
    Follicular_B_cells, Plasma_B_cells,
    Plasmacytoid_dendritic_cells
  )
  
  markers <- c(Plasma_B_cells, Follicular_B_cells, Germinal_center_B_cells)
  
  
  
  DotPlot(scRNA_harmony, features = unique(markers), group.by = "seurat_clusters") + RotatedAxis() +
    scale_x_discrete("") + scale_y_discrete("")
  
  new.cluster.ids <- c(
    "UnDef1",
    "UnDef2", "UnDef3",
    "UnDef4", "UnDef0",
    "Plasma_B_cell1", "UnDef5",
    "UnDef6", "Plasma_B_cell2",
    "UnDef7", "UnDef8",
    "UnDef9", "Plasma_B_cell3",
    "Plasma_B_cell4", "UnDef10",
    "Plasma_B_cell5", "UnDef11",
    "UnDef12", "UnDef16",
    "Plasma_B_cell6", "Plasma_B_cell7",
    "UnDef13", "Plasma_B_cell8",
    "UnDef14", "UnDef15"
  )
  
  
  scRNA_harmony@meta.data$new.cluster.ids <- NA
  
  for (i in 1:length(scRNA_harmony@meta.data$new.cluster.ids)) {
    # if (scRNA_harmony@meta.data$seurat_clusters[i]%in%c(5,8,12,13,19,20,22)) {
    #   print('1')
    #   scRNA_harmony@meta.data$new.cluster.ids[i] <- "Plasma_B_cell"
    # }else{
    #   scRNA_harmony@meta.data$new.cluster.ids[i] <- "UnDef"
    # }
    scRNA_harmony@meta.data$new.cluster.ids[i] <- new.cluster.ids[scRNA_harmony@meta.data$seurat_clusters[i]]
  }
  ###
  
  # 获取当前用的Idents
  Idents(object = scRNA_harmony)
  levels(scRNA_harmony)
  Idents(scRNA_harmony) <- "seurat_clusters"
  new.cluster.ids <- as.character(new.cluster.ids)
  names(new.cluster.ids) <- levels(scRNA_harmony)
  scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)
  
  scRNA_harmony@meta.data$new.cluster.ids
  ####
  scRNA_harmony <- RenameIdents(scRNA_harmony,
                                `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
                                `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK",
                                `7` = "T activated", `8` = "DC", `9` = "B Activated",
                                `13` = "Mono/Mk Doublets", `14` = "HSPC"
  )
  DimPlot(scRNA_harmony, label = T, cols = colP) + NoLegend()
  DimPlot(scRNA_harmony, reduction = "tsne", label = T, cols = colP, repel = T) + NoLegend()
  #
  
  DimPlot(scRNA_harmony, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5, cols = colP) + NoLegend()
  DimPlot(scRNA_harmony, reduction = "umap", group.by = "new.cluster.ids", label = TRUE, pt.size = 0.5, cols = colP) + NoLegend()
  ####
  DotPlot(scRNA_harmony, features = unique(markers), group.by = "new.cluster.ids") + RotatedAxis() +
    scale_x_discrete("") + scale_y_discrete("")
  
  print(DimPlot(scRNA_harmony, reduction = "umap", group.by = c("seurat_clusters", "new.cluster.ids"), cols = colP))
  
  p <- FeaturePlot(scRNA_harmony, features = markers, ncol = 8)
  p
}




# seurat ------------------------------------------------------------------

################################################################################
#   _____                      _
#  / ____|                    | |
# | (___   ___ _   _ _ __ __ _| |_
#  \___ \ / _ \ | | | '__/ _` | __|
#  ____) |  __/ |_| | | | (_| | |_
# |_____/ \___|\__,_|_|  \__,_|\__|
#
################################################################################

if (F) {
  ### seurat锚点整合
  rm(list = ls())
  library(future) # Seurat并行计算的一个包，不加载这个包不能进行并行计算
  options(future.globals.maxSize = 50 * 1024^3) # 将全局变量上限调至50G（锚点整合很占内存）
  
  ## 重新创建没有处理的经过降维等处理的数据
  seu_obj_filter <- readRDS("seu_obj_filter.rds")
  cellinfo <- subset(seu_obj_filter@meta.data, select = c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP"))
  
  scRNA_seurat <- CreateSeuratObject(seu_obj_filter@assays$RNA@counts, meta.data = cellinfo)
  # 做锚点整合需要把样本处理成单个的Seurat对象来两两组合
  scRNA_seuratlist <- SplitObject(scRNA_seurat, split.by = "orig.ident")
  # 也可以按别的指标（metadata中的）来进行拆分，比如可以按不同的分组来拆分样本，再进行整合。
  # SCTransform标准化(⚠️使用log标准化还是SCT标准化差别不大)
  # 如果用log标准化，后面FindIntegrationAnchors和IntegrateData函数的normalization.method参数选'LogNormalize'
  scRNA_seuratlist <- parallel::mclapply(scRNA_seuratlist, FUN = function(x) SCTransform(x), mc.cores = 10) # 10个对象最好写10个核，没有10个核少写几个也可以。top命令可以查看服务器有几个核，mc.core设置为1就每次处理一个对象。
  # mclapply是lapply的多核版本
  ### FindAnchors
  ### 每个样本的高变基因不完全一样，SelectIntegrationFeatures可以整合这些高变基因，选出3000个
  scRNA_seurat.features <- SelectIntegrationFeatures(scRNA_seuratlist, nfeatures = 3000)
  ### 将每个样本的高变基因都调整成上一步选出的3000个
  scRNA_seuratlist <- PrepSCTIntegration(scRNA_seuratlist, anchor.features = scRNA_seurat.features)
  ## 寻找锚点，运行速度非常慢，至少需要1-2小时
  plan("multisession", workers = 10)
  scRNA_seurat.anchors <- FindIntegrationAnchors(
    object.list = scRNA_seuratlist,
    normalization.method = "SCT", # 如果前面是log标准化，这里改成LogNormalize
    anchor.features = scRNA_seurat.features
  )
  ### Integrate 运行速度慢
  scRNA_seurat.sct.int <- IntegrateData(scRNA_seurat.anchors, normalization.method = "SCT") # 速度慢
  plan("sequential") # 把并行计算改为单核计算
  ### redunction
  scRNA_seurat <- RunPCA(scRNA_seurat.sct.int, npcs = 50, verbose = FALSE)
  ElbowPlot(scRNA_seurat, ndims = 50)
  pc.num <- 1:20
  scRNA_seurat <- scRNA_seurat %>%
    RunTSNE(dims = pc.num) %>%
    RunUMAP(dims = pc.num)
  
  ### Visual
  DimPlot(scRNA_seurat, group.by = "orig.ident")
  p <- DimPlot(scRNA_seurat, group.by = "orig.ident")
  ggsave("UMAP_Samples_seurat.png", p, width = 8, height = 6)
  DimPlot(scRNA_seurat, group.by = "orig.ident", split.by = "orig.ident", ncol = 4)
  p <- DimPlot(scRNA_seurat, group.by = "orig.ident", split.by = "orig.ident", ncol = 4)
  ggsave("UMAP_Samples_Split_seurat.png", p, width = 18, height = 12)
  
  
  scRNA_seurat <- FindNeighbors(scRNA_seurat, dims = 1:20) %>% FindClusters()
  ## save seurat object
  saveRDS(scRNA_seurat, "scRNA_SCT_int_seurat.rds")
  
  ### 结果评估
  # scRNA_seurat <- readRDS("scRNA_SCT_int_seurat.rds")
  
  load("ref_Hematopoietic.RData")
  DefaultAssay(scRNA_seurat) <- "RNA"
  scRNA_seurat <- cell_identify(scRNA_seurat, ref_Hematopoietic) # cell_identify是自己写的函数，ref_Hematopoietic是SingleR中的参考数据集
  p <- DimPlot(scRNA_seurat, group.by = "SingleR", label = T)
  ggsave("SingleR_Seurat.png", p, width = 8, height = 6)
  
  ###
  DefaultAssay(scRNA_seurat) <- "integrated"
  p <- FeaturePlot(scRNA_seurat, features = c(
    "CD3E", "CD3G", "CD79B", "GNLY", "NKG7",
    "CD14", "FCGR3A", "CD68", "S100A12", "CST3",
    "FCER1A", "GZMB", "IL3RA"
  ), ncol = 4)
  ggsave("Features_Seurat_int.png", p, width = 18, height = 16)
  ###
  DefaultAssay(scRNA_seurat) <- "SCT"
  p <- FeaturePlot(scRNA_seurat, features = c(
    "CD3E", "CD3G", "CD79B", "GNLY", "NKG7",
    "CD14", "FCGR3A", "CD68", "S100A12", "CST3",
    "FCER1A", "GZMB", "IL3RA"
  ), ncol = 4)
  ggsave("Features_Seurat_sct.png", p, width = 18, height = 16)
}
