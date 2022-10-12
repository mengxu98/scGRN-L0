

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
  library(reticulate)
  library(sctransform)
  library(viridis)
  library(tidyr)
  library(magrittr)
  library(reshape2)
  library(readxl)
  library(progeny)
  library(readr)
  library(stringr)
  library(SingleR)
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
  library(glmGamPoi) # 加速SCT
  source("step_function.R")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  ## 下载注释数据库
  # hpca.se <- HumanPrimaryCellAtlasData()
  # hpca.se
  # load("/data/mengxu/data/SingleR_data/HumanPrimaryCellAtlas_hpca.se_human.RData")
  # load("/data/mengxu/data/SingleR_data/BlueprintEncode_bpe.se_human.RData")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

}

# Raw
if (T) {
  stage <- c("normal", "1", "2", "3", "4")

  for (s in stage) {
    if (file.exists(paste0("/data/mengxu/data/all/lung_stage-", s, "_list_filter.Rdata")) == T) {
      load(paste0("/data/mengxu/data/all/lung_stage-", s, "_list_filter.Rdata"))
      seu_obj_filter <- merge_seu_obj(seu_obj_list_filter, samples, stage = s) # stage = "normal" or "1" or "2" or "3" or "4"
      save(seu_obj_filter, samples, file = paste0("/data/mengxu/data/all/lung_stage-", s, "_seu.Rdata"))
      rm(seu_obj_list_filter)
      rm(samples)
      gc()
    }
  }
}
# Filter
if (T) {
  stage <- c("normal", "1", "2", "3", "4")
  dataset_information <- c()
  for (s in stage) {
    if (file.exists(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu.Rdata")) == T) {
      load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu.Rdata"))

      # calculate mitochondrial, hemoglobin and ribosomal gene counts
      seu_obj_filter <- PercentageFeatureSet(seu_obj_filter, pattern = "^MT-", col.name = "pMT")
      seu_obj_filter <- PercentageFeatureSet(seu_obj_filter, pattern = "^HBA|^HBB", col.name = "pHB")
      seu_obj_filter <- PercentageFeatureSet(seu_obj_filter, pattern = "^RPS|^RPL", col.name = "pRP")

      ###
      # qcparams <- c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP")
      # for (i in seq_along(qcparams)){
      #   print(VlnPlot(object = seu_obj_filter, features = qcparams[i], group.by = "orig.ident", pt.size = 0))
      # }
      # for (i in seq_along(qcparams)){
      #   print(RidgePlot(object = seu_obj_filter, features = qcparams[i], group.by = "orig.ident"))
      # }
      #
      # VlnPlot(seu_obj_filter, features = c("nFeature_RNA", "nCount_RNA", "pMT"), pt.size = 0, group.by = "orig.ident", ncol = 1, log = T)
      # ggsave2("SuppFig1.pdf", path = "../results", width = 30, height = 25, units = "cm")
      ###

      seu_obj <- subset(seu_obj_filter,
        subset = nFeature_RNA > nFeature_lower &
          nFeature_RNA < nFeature_upper &
          nCount_RNA > nCount_lower &
          nCount_RNA < nCount_upper &
          pMT < pMT_upper &
          pHB < pHB_upper
      )
      ###
      dim(seu_obj_filter)
      dim(seu_obj)
      # qc_std_plot_nf(seu_obj_filter)
      qc_std_plot(seu_obj_filter)
      ggsave2(paste0("SuppFig.1_stage-", s, ".png"),
        path = paste0("/data/mengxu/results/figure/stage-", s),
        width = 30, height = 30, units = "cm"
      )

      # VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "pMT"), pt.size = 0, group.by = "orig.ident", ncol = 1, log = T)
      ###
      save(seu_obj_filter, samples, file = paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_filter.Rdata"))

      dataset_information_one <- c(paste0("Stage-", s), length(colnames(seu_obj_filter)), length(colnames(seu_obj)), length(rownames(seu_obj)))

      rm(seu_obj_filter)
      gc()
    }

    dataset_information <- rbind.data.frame(
      dataset_information,
      dataset_information_one
    )

    names(dataset_information) <- c("Stage", "No. of cells(raw)", "No. of cells(filter)", "No. of genes")
  }

  write.csv(dataset_information, file = paste0("/data/mengxu/data/all/dataset_information.csv"), row.names = F)
}


# harmony -----------------------------------------------------------------

if (T) {
  stage <- c("normal", "1", "2", "3", "4")

  for (s in stage) {

    # s="normal"
    # load('/data/mengxu/data/all/lung_stage-normal_seu_filter.Rdata')
    load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_filter.Rdata"))

    seu_obj_data <- seu_obj_filter
    rm(seu_obj_filter)
    gc()
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
    ggsave2(paste0("SuppFig.2_stage-", s, "_CellCycleScoring_raw.png"),
      path = paste0("/data/mengxu/results/figure/stage-", s),
      width = 10, height = 8, units = "cm"
    )

    seu_obj_data$CC.Difference <- seu_obj_data$S.Score - seu_obj_data$G2M.Score

    seu_obj_data <- SCTransform(seu_obj_data,
      method = "glmGamPoi",
      vars.to.regress = "CC.Difference",
      # features = rownames(seu_obj_data),
      conserve.memory = T
    )
    seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
    DimPlot(seu_obj_data)
    ggsave2(paste0("SuppFig.2_stage-", s, "_CellCycleScoring_SCT.png"),
      path = paste0("/data/mengxu/results/figure/stage-", s),
      width = 10, height = 8, units = "cm"
    )
    # seu_obj_data <- ScaleData(seu_obj_data,
    #                           vars.to.regress = "CC.Difference",
    #                           features = rownames(seu_obj_data))
    # 观察细胞周期基因的表达情况
    # RidgePlot(seu_obj_data, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
    # RidgePlot(seu_obj_data, features = c(s.genes[1:25], g2m.genes[1:25]), ncol = 8)
    # RidgePlot(seu_obj_data, features = c(s.genes[26:length(s.genes)], g2m.genes[26:length(g2m.genes)]), ncol = 8)
    if (F) {

      # Cell cycle scoring

      ### add cell cycle, cc.genes loaded with Seurat

      # score_cc <- function(seu_obj) {
      #   seu_obj <- CellCycleScoring(seu_obj, s.genes, g2m.genes)
      #   seu_obj@meta.data$CC.Diff <- seu_obj@meta.data$S.Score - seu_obj@meta.data$G2M.Score
      #   return(seu_obj)
      # }
      #
      # seu_obj <- score_cc(seu_obj)
      #
      # FeatureScatter(scRNA_harmony, "G2M.Score", "S.Score", group.by = "Phase", pt.size = .1) +
      #   coord_fixed(ratio = 1)


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
      DimPlot(seu_obj_data, reduction = "Phase") #+
      # scale_color_viridis(discrete = T, option = "C")
      # 计算并移除分数差异
      if (T) {
        # A
        seu_obj_data <- ScaleData(seu_obj_data,
          vars.to.regress = c("S.Score", "G2M.Score"),
          features = rownames(seu_obj_data)
        )
      } else if (F) {
        # B
        # 计算并移除分数差异
        seu_obj_data$CC.Difference <- seu_obj_data$S.Score - seu_obj_data$G2M.Score

        seu_obj_data <- ScaleData(seu_obj_data,
          vars.to.regress = "CC.Difference",
          features = rownames(seu_obj_data)
        )
      }
    }
    # ElbowPlot(seu_obj_data)
    save(seu_obj_data,
      file = paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_filter_SCT.Rdata")
    )
    # load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_filter_SCT.Rdata")

    pc <- pc_num(seu_obj_data) # 自动推断合适的'pc'值
    pc.num <- 1:pc

    seu_obj_data <- seu_obj_data %>% RunUMAP(dims = pc.num) # %>% RunTSNE(dims=pc.num) #'tsne' compute is too slow

    seu_obj_data <- FindNeighbors(seu_obj_data, dims = pc.num) %>% FindClusters(resolution = 0.3)

    save(seu_obj_data,
      file = paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_filter_SCT_PCA.Rdata")
    )
    # load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_filter_SCT_PCA.Rdata"))
    

    if (F) {
      seu_obj_data <- annotation_celltype(seu_obj_data, method = "singleR") # method = "celltypist" or "singleR"
      levels(seu_obj_data$seurat_clusters)
    }
    

    if (F) {
      ITG.suj <- RunUMAP(ITG.spj, reduction = "pca", dims = 1:10) # ???
    }


    # saveRDS(seu_obj_filter, "seu_obj_filter.rds")

    # ElbowPlot(seu_obj_data, ndims = 50)

    ### plot_tsne

    # DimPlot(seu_obj_data,
    #         reduction = "tsne",
    #         #group.by = "orig.ident",
    #         label = T,
    #         cols = colP,
    #         repel = T
    #         #pt.size = 0.2
    # )+theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid")) +NoLegend()
    # ggsave2("fig2.clusters_raw_tsne.png", path = "/data/mengxu/results/figure", width = 12, height = 12, units = "cm")
    #
    # DimPlot(seu_obj_filter,
    #         reduction = "tsne",
    #         group.by = "orig.ident",
    #         #label = T,
    #         repel = T , cols = colP
    #         #pt.size = 0.2
    # ) +theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
    # ggsave2("fig2.samples_raw_tsne.png", path = "/data/mengxu/results/figure", width = 30, height = 15, units = "cm")
    #
    # DimPlot(seu_obj_filter, reduction = "tsne",group.by = "platform", label = TRUE, pt.size = 0.5, cols = colP)+
    #   theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))#+NoLegend()
    # ggsave2("fig2.platform_raw_tsne.png", path = "/data/mengxu/results/figure", width = 14, height = 11, units = "cm")

    ### plot_umap

    DimPlot(seu_obj_data,
      reduction = "umap",
      group.by = "seurat_clusters",
      label = T,
      cols = colP,
      repel = T
      # pt.size = 0.2
    ) + theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + # NoLegend() +
      scale_colour_viridis_d(option = "inferno")
    ggsave2("SuppFig.3.clusters_raw_umap.png",
      path = paste0("/data/mengxu/results/figure/stage-", s),
      width = 13, height = 12, units = "cm"
    )

    DimPlot(seu_obj_data,
      reduction = "umap",
      group.by = "orig.ident",
      # label = T,
      repel = T,
      cols = colP
      # pt.size = 0.2
    ) + theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + # NoLegend() +
      scale_colour_viridis_d(option = "inferno")
    if (length(samples) <= 10) {
      ggsave2("SuppFig.3.samples_raw_umap.png",
        path = paste0("/data/mengxu/results/figure/stage-", s),
        width = 18, height = 12, units = "cm"
      )
    } else if (length(samples) > 10 && length(samples) <= 20) {
      ggsave2("SuppFig.3.samples_raw_umap.png",
        path = paste0("/data/mengxu/results/figure/stage-", s),
        width = 20, height = 12, units = "cm"
      )
    } else if (length(samples) > 30 && length(samples) <= 40) {
      ggsave2("SuppFig.3.samples_raw_umap.png",
        path = paste0("/data/mengxu/results/figure/stage-", s),
        width = 25, height = 12, units = "cm"
      )
    } else if (length(samples) > 40) {
      ggsave2("SuppFig.3.samples_raw_umap.png",
        path = paste0("/data/mengxu/results/figure/stage-", s),
        width = 30, height = 12, units = "cm"
      )
    }


    # DimPlot(seu_obj_data, reduction = "umap",group.by = "platform", label = TRUE, pt.size = 0.5, cols = colP)+
    #   theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))#+NoLegend()
    # ggsave2("fig2.platform_raw_umap.png",
    #         path = paste0("/data/mengxu/results/figure/stage-", s),
    #         width = 14, height = 11, units = "cm")

    ################################################################################
    #  _    _
    # | |  | |
    # | |__| | __ _ _ __ _ __ ___   ___  _ __  _   _
    # |  __  |/ _` | '__| '_ ` _ \ / _ \| '_ \| | | |
    # | |  | | (_| | |  | | | | | | (_) | | | | |_| |
    # |_|  |_|\__,_|_|  |_| |_| |_|\___/|_| |_|\__, |
    #                                           __/ |
    #                                          |___/
    ################################################################################

    ### harmony整合
    # rm(list=ls())
    # seu_obj_filter <- readRDS("seu_obj_filter.rds")
    # cellinfo <- subset(seu_obj_data@meta.data, select = c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP"))
    # scRNA_harmony <- CreateSeuratObject(seu_obj_filter@assays$RNA@counts, meta.data = cellinfo)

    ### SCT标准化数据
    # scRNA_harmony <- SCTransform(scRNA_harmony, method = "glmGamPoi") #通过'glmGamPoi'加速

    ### PCA
    # scRNA_harmony <- RunPCA(scRNA_harmony, npcs=50, verbose=FALSE)

    ### 整合方法1：单个样本间进行整合（推荐，效果更好）
    if (F) {
      Anchors <- FindIntegrationAnchors(object.list = dataset.list, dims = 1:30)

      ### Part B. Integration
      ITG.sbj <- IntegrateData(anchorset = Anchors, dims = 1:30)
      DefaultAssay(ITG.sbj) <- "integrated"
    }


    scRNA_harmony <- seu_obj_data
    scRNA_harmony <- RunHarmony(scRNA_harmony,
      group.by.vars = "orig.ident",
      assay.use = "SCT",
      max.iter.harmony = 20
    )

    # scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars="orig.ident", max.iter.harmony = 20)
    # group.by.vars参数是设置按哪个分组来整合
    # max.iter.harmony设置迭代次数，默认是10。运行RunHarmony结果会提示在迭代多少次后完成了收敛。
    # ⚠️RunHarmony函数中lambda参数，默认值是1，决定了Harmony整合的力度。
    # lambda值调小，整合力度变大，反之。（只有这个参数影响整合力度，调整范围一般在0.5-2之间）

    ElbowPlot(scRNA_harmony, ndims = 50)
    pc <- pc_num(scRNA_harmony) # 自动推断合适的'pc'值
    pc.num <- 1:pc
    scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num) # %>% RunTSNE(reduction="harmony", dims=pc.num)

    save(seu_obj_data,
      file = paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_SCT_harmony_PCA.Rdata")
    )
    # load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_SCT_harmony_PCA.Rdata"))

    scRNA_harmony <- FindNeighbors(scRNA_harmony, dims = pc.num) %>% FindClusters(resolution = 1)

    ### plot_tsne

    # DimPlot(scRNA_harmony,
    #         reduction = "tsne",
    #         #group.by = "orig.ident",
    #         label = T,
    #         cols = colP,
    #         repel = T
    #         #pt.size = 0.2
    # )+theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid")) +NoLegend()
    # ggsave2("fig3.clusters_harmony_tsne.png", path = "/data/mengxu/results/figure", width = 12, height = 12, units = "cm")
    #
    # DimPlot(scRNA_harmony,
    #         reduction = "tsne",
    #         group.by = "orig.ident",
    #         #label = T,
    #         repel = T , cols = colP
    #         #pt.size = 0.2
    # ) +theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))
    # ggsave2("fig3.samples_harmony_tsne.png", path = "/data/mengxu/results/figure", width = 30, height = 15, units = "cm")
    #
    # DimPlot(scRNA_harmony, reduction = "tsne",group.by = "platform", label = TRUE, pt.size = 0.5, cols = colP)+
    #   theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))#+NoLegend()
    # ggsave2("fig3.platform_harmony_tsne.png", path = "/data/mengxu/results/figure", width = 14, height = 11, units = "cm")

    ### plot_umap
    ###
    DimPlot(scRNA_harmony,
      reduction = "umap",
      # group.by = "orig.ident",
      label = T,
      cols = colP,
      repel = T
      # pt.size = 0.2
    ) + theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + # NoLegend() +
      scale_colour_viridis_d(option = "inferno")
    ggsave2("SuppFig3.clusters_harmony_umap.png",
      path = paste0("/data/mengxu/results/figure/stage-", s),
      width = 12, height = 12, units = "cm"
    )
    ###
    DimPlot(scRNA_harmony,
      reduction = "umap",
      group.by = "orig.ident",
      # label = T,
      repel = T,
      # cols = viridis(72)
      cols = colP
      # pt.size = 0.2
    ) + theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
      scale_colour_viridis_d(option = "inferno")
    ggsave2("SuppFig3.samples_harmony_umap.png",
      path = paste0("/data/mengxu/results/figure/stage-", s),
      width = 20, height = 15, units = "cm"
    )
    ###
    # DimPlot(scRNA_harmony, reduction = "umap",group.by = "platform", label = TRUE, pt.size = 0.5, cols = colP)+
    #   theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))#+NoLegend()
    # ggsave2("SuppFig3.platform_harmony_umap.png",
    #         path = paste0("/data/mengxu/results/figure/stage-", s),
    #         width = 14, height = 11, units = "cm")

    #----------------------------------------------------------------------------#

    sce.merged <- as.SingleCellExperiment(scRNA_harmony)

    # LISI index

    lisi.pca <- lisi::compute_lisi(
      reducedDim(sce.merged, "PCA"),
      colData(sce.merged), c("orig.ident", "platform")
    ) # "dataset","dataset.tech","ClusterID.pca" ,"RNA_snn_res.0.3"

    lisi.harmony <- lisi::compute_lisi(
      reducedDim(sce.merged, "HARMONY"),
      colData(sce.merged), c("orig.ident", "platform")
    ) # c("dataset","dataset.tech","ClusterID.harmony") ,"RNA_snn_res.0.3"

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

    # plot
    p <-
      ggboxplot(lisi.merge.tb,
        x = "method", y = "orig.ident",
        fill = "method", alpha = 0.8
      ) +
      stat_compare_means(comparisons = list(c("raw", "harmony"))) +
      ylab("LISI") +
      theme(legend.position = "right") + scale_fill_tron()

    print(p)

    ggsave2(paste0("/data/mengxu/results/figure/stage-", s, "/SuppFig4.LISI.dataset.merge.png"),
      width = 3.2, height = 4
    )

    # ggsave2("SuppFig4.LISI.dataset.merge.png",
    #         path = paste0("/data/mengxu/results/figure/stage-", s),
    #         width = 10, height = 10, units = "cm")
    #-----------------------------------------------------------------------#

    mainmarkers <- c( # Nature Medicine-Phenotype molding of stromal cells in the lung  tumor microenvironment
      # "ACAT2",
      "CLDN18", # Alveolar
      "CLDN5", # Endothelial
      "CAPS", # Epithelial
      "ALB", # Hepatocytes
      "COL1A1", # Fibroblast
      "CD79A", # B cell
      "LYZ", # Myeloid
      "CD3D", # T cell
      "EPCAM" # , #Cancer
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
      path = paste0("/data/mengxu/results/figure/stage-", s),
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
          path = paste0("/data/mengxu/results/figure/stage-", s, "/marker/"),
          width = 10, height = 10, units = "cm"
        )
      }
      ###
    }
    ###
  }
  ###
}


# annotation --------------------------------------------------------------

###
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
