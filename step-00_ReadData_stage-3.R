

# load libraries
if (T) {
  rm(list = ls())
  gc()
  library(Seurat)
  library(dplyr)
  library(reticulate)
  library(sctransform)
  library(cowplot)
  library(ggplot2)
  library(viridis)
  library(tidyr)
  library(magrittr)
  library(reshape2)
  library(readxl)
  library(progeny)
  library(readr)
  library(stringr)
  library(Matrix)
  library(scran)
  library(scuttle)
  library(harmony)
  library(tidyverse)
  library(patchwork)
  library(data.table) # fread
}
#------------------------------------------------------------------------------------------------------#
# PhilipBischoff2021_10X
if (T) {
  # Single-cell gene expression data of all patients were merged,
  # and transcriptomes were filtered for cells with 500–10,000 genes detected,
  # 1000–100,000 UMIs counted, fraction of mitochondrial reads <30%, and fraction of hemoglobin reads <5%.
  # After filtering, UMI counts were variance-stabilized using scTransform with 3000 variable features,
  # while regressing out number of UMIs and fraction of mitochondrial reads.

  ###
  nFeature_lower <- 500
  nFeature_upper <- 10000
  nCount_lower <- 1000
  nCount_upper <- 100000
  pMT_lower <- 0
  pMT_upper <- 30
  pHB_lower <- 0
  pHB_upper <- 5

  PhilipBischoff2021.P030 <- Read10X("/data/mengxu/data/PhilipBischoff2021/cellranger/p030t/filtered_feature_bc_matrix")

  PhilipBischoff2021.P030 <- CreateSeuratObject(
    counts = PhilipBischoff2021.P030,
    project = "PhilipBischoff2021.P030",
    min.features = 200,
    min.cells = 3
  )

  PhilipBischoff2021.P030$platform <- "10X"
  ### QC
  PhilipBischoff2021.P030 <- PercentageFeatureSet(PhilipBischoff2021.P030, pattern = "^MT-", col.name = "pMT")
  PhilipBischoff2021.P030 <- PercentageFeatureSet(PhilipBischoff2021.P030, pattern = "^HBA|^HBB", col.name = "pHB")
  PhilipBischoff2021.P030 <- PercentageFeatureSet(PhilipBischoff2021.P030, pattern = "^RPS|^RPL", col.name = "pRP")
  PhilipBischoff2021.P030 <- subset(PhilipBischoff2021.P030,
    subset = nFeature_RNA > nFeature_lower &
      nFeature_RNA < nFeature_upper &
      nCount_RNA > nCount_lower &
      nCount_RNA < nCount_upper &
      pMT < pMT_upper &
      pHB < pHB_upper
  )


  ###
  PhilipBischoff2021.P033 <- Read10X("/data/mengxu/data/PhilipBischoff2021/cellranger/p033t/filtered_feature_bc_matrix")

  PhilipBischoff2021.P033 <- CreateSeuratObject(
    counts = PhilipBischoff2021.P033,
    project = "PhilipBischoff2021.P033",
    min.features = 200,
    min.cells = 3
  )

  PhilipBischoff2021.P033$platform <- "10X"
  ### QC
  PhilipBischoff2021.P033 <- PercentageFeatureSet(PhilipBischoff2021.P033, pattern = "^MT-", col.name = "pMT")
  PhilipBischoff2021.P033 <- PercentageFeatureSet(PhilipBischoff2021.P033, pattern = "^HBA|^HBB", col.name = "pHB")
  PhilipBischoff2021.P033 <- PercentageFeatureSet(PhilipBischoff2021.P033, pattern = "^RPS|^RPL", col.name = "pRP")
  PhilipBischoff2021.P033 <- subset(PhilipBischoff2021.P033,
    subset = nFeature_RNA > nFeature_lower &
      nFeature_RNA < nFeature_upper &
      nCount_RNA > nCount_lower &
      nCount_RNA < nCount_upper &
      pMT < pMT_upper &
      pHB < pHB_upper
  )
}

#------------------------------------------------------------------------------------------------------#
# E-MTAB-6653_DietherLambrechts2018
if (T) {
  # From  this, all cells were removed that had either fewer than 201 UMIs,
  # over 6,000  or below 101 expressed genes, or over 10% UMIs derived from mitochondrial  genome.
  ###
  nFeature_lower <- 201
  # nFeature_upper <- 10000
  nCount_lower <- 101
  nCount_upper <- 6000
  pMT_lower <- 0
  pMT_upper <- 10
  # pHB_lower <- 0
  # pHB_upper <- 5

  DietherLambrechts2018.P06.1 <- Read10X("/data/mengxu/data/E-MTAB-6653/scrBT1429m/scrBT1429m_f/outs/filtered_feature_bc_matrix")

  DietherLambrechts2018.P06.1 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P06.1,
    project = "DietherLambrechts2018.P06.1",
    min.features = 200,
    min.cells = 3
  )

  DietherLambrechts2018.P06.1$platform <- "10X"
  ### QC
  DietherLambrechts2018.P06.1 <- PercentageFeatureSet(DietherLambrechts2018.P06.1, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P06.1 <- PercentageFeatureSet(DietherLambrechts2018.P06.1, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P06.1 <- PercentageFeatureSet(DietherLambrechts2018.P06.1, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P06.1 <- subset(DietherLambrechts2018.P06.1,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P06.2 <- Read10X("/data/mengxu/data/E-MTAB-6653/scrBT1430m/scrBT1430m_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P06.2 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P06.2,
    project = "DietherLambrechts2018.P06.2",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P06.2$platform <- "10X"
  ### QC
  DietherLambrechts2018.P06.2 <- PercentageFeatureSet(DietherLambrechts2018.P06.2, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P06.2 <- PercentageFeatureSet(DietherLambrechts2018.P06.2, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P06.2 <- PercentageFeatureSet(DietherLambrechts2018.P06.2, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P06.2 <- subset(DietherLambrechts2018.P06.2,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )

  ###
  DietherLambrechts2018.P06.3 <- Read10X("/data/mengxu/data/E-MTAB-6653/scrBT1431m/scrBT1431m_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P06.3 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P06.3,
    project = "DietherLambrechts2018.P06.3",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P06.3$platform <- "10X"
  ### QC
  DietherLambrechts2018.P06.3 <- PercentageFeatureSet(DietherLambrechts2018.P06.3, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P06.3 <- PercentageFeatureSet(DietherLambrechts2018.P06.3, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P06.3 <- PercentageFeatureSet(DietherLambrechts2018.P06.3, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P06.3 <- subset(DietherLambrechts2018.P06.3,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P06 <- merge(DietherLambrechts2018.P06.1,
    y = c(
      DietherLambrechts2018.P06.2,
      DietherLambrechts2018.P06.3
    ),
    add.cell.ids = c(
      "DietherLambrechts2018.P06",
      "DietherLambrechts2018.P06",
      "DietherLambrechts2018.P06"
    )
  )

  ### E-MTAB-6149_DietherLambrechts2018

  DietherLambrechts2018.P03.1 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT1290/BT1290_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P03.1 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P03.1,
    project = "DietherLambrechts2018.P03.1",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P03.1$platform <- "10X"
  ### QC
  DietherLambrechts2018.P03.1 <- PercentageFeatureSet(DietherLambrechts2018.P03.1, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P03.1 <- PercentageFeatureSet(DietherLambrechts2018.P03.1, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P03.1 <- PercentageFeatureSet(DietherLambrechts2018.P03.1, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P03.1 <- subset(DietherLambrechts2018.P03.1,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P03.2 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT1291/BT1291_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P03.2 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P03.2,
    project = "DietherLambrechts2018.P03.2",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P03.2$platform <- "10X"
  ### QC
  DietherLambrechts2018.P03.2 <- PercentageFeatureSet(DietherLambrechts2018.P03.2, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P03.2 <- PercentageFeatureSet(DietherLambrechts2018.P03.2, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P03.2 <- PercentageFeatureSet(DietherLambrechts2018.P03.2, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P03.2 <- subset(DietherLambrechts2018.P03.2,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P03.3 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT1292/BT1292_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P03.3 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P03.3,
    project = "DietherLambrechts2018.P03.3",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P03.3$platform <- "10X"
  ### QC
  DietherLambrechts2018.P03.3 <- PercentageFeatureSet(DietherLambrechts2018.P03.3, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P03.3 <- PercentageFeatureSet(DietherLambrechts2018.P03.3, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P03.3 <- PercentageFeatureSet(DietherLambrechts2018.P03.3, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P03.3 <- subset(DietherLambrechts2018.P03.3,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )

  ###
  DietherLambrechts2018.P03 <- merge(DietherLambrechts2018.P03.1,
    y = c(
      DietherLambrechts2018.P03.2,
      DietherLambrechts2018.P03.3
    ),
    add.cell.ids = c(
      "DietherLambrechts2018.P03",
      "DietherLambrechts2018.P03",
      "DietherLambrechts2018.P03"
    )
  )
}

#------------------------------------------------------------------------------------------------------------------#
# GSE127465_RapolasZilionis2019_IndropSeq
if (T) {
  # To retain high quality transcriptomes, total count and a mitochondrial count filters were applied.
  # For mouse, transcriptomes with more than 600 total counts and less than 15% of counts coming from mitochondrial genes were retained.
  # For human, these thresholds were 300 total counts and 20%, respectively.
  # A more permissive filtering was necessary to
  # i) avoid filtering out human neutrophils, that naturally have a lower mRNA content;
  # ii) retain non-hematopoietic cells which are more sensitive to dissociation and show a higher mitochondrial gene fraction,
  # possibly indicating premature lysis. After this initial cleanup, data was normalized by total counts as described (Klein et al., 2015).

  # nFeature_lower <- 201
  # nFeature_upper <- 10000
  nCount_lower <- 300
  # nCount_upper <- 6000
  pMT_lower <- 0
  pMT_upper <- 20
  # pHB_lower <- 0
  # pHB_upper <- 5
  ###
  RapolasZilionis2019.P2.1 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635285_human_p2t1_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )
  RapolasZilionis2019.P2.2 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635286_human_p2t2_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )
  RapolasZilionis2019.P2 <- cbind.data.frame(
    t(RapolasZilionis2019.P2.1),
    t(RapolasZilionis2019.P2.2)
  )

  RapolasZilionis2019.P2 <- CreateSeuratObject(
    counts = RapolasZilionis2019.P2,
    project = "RapolasZilionis2019.P2",
    min.features = 200,
    min.cells = 3
  )
  RapolasZilionis2019.P2$platform <- "IndropSeq"
  ### QC
  RapolasZilionis2019.P2 <- PercentageFeatureSet(RapolasZilionis2019.P2, pattern = "^MT-", col.name = "pMT")
  RapolasZilionis2019.P2 <- PercentageFeatureSet(RapolasZilionis2019.P2, pattern = "^HBA|^HBB", col.name = "pHB")
  RapolasZilionis2019.P2 <- PercentageFeatureSet(RapolasZilionis2019.P2, pattern = "^RPS|^RPL", col.name = "pRP")
  RapolasZilionis2019.P2 <- subset(RapolasZilionis2019.P2,
    subset =
    # nFeature_RNA > nFeature_lower &
    # nFeature_RNA < nFeature_upper &
      nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  RapolasZilionis2019.P7.1 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635301_human_p7t1_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )
  RapolasZilionis2019.P7.2 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635302_human_p7t2_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )
  RapolasZilionis2019.P7 <- cbind.data.frame(
    t(RapolasZilionis2019.P7.1),
    t(RapolasZilionis2019.P7.2)
  )

  RapolasZilionis2019.P7 <- CreateSeuratObject(
    counts = RapolasZilionis2019.P7,
    project = "RapolasZilionis2019.P7",
    min.features = 200,
    min.cells = 3
  )
  RapolasZilionis2019.P7$platform <- "IndropSeq"
  ### QC
  RapolasZilionis2019.P7 <- PercentageFeatureSet(RapolasZilionis2019.P7, pattern = "^MT-", col.name = "pMT")
  RapolasZilionis2019.P7 <- PercentageFeatureSet(RapolasZilionis2019.P7, pattern = "^HBA|^HBB", col.name = "pHB")
  RapolasZilionis2019.P7 <- PercentageFeatureSet(RapolasZilionis2019.P7, pattern = "^RPS|^RPL", col.name = "pRP")
  RapolasZilionis2019.P7 <- subset(RapolasZilionis2019.P7,
    subset =
    # nFeature_RNA > nFeature_lower &
    # nFeature_RNA < nFeature_upper &
      nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
}

#-------------------------------------------------------------------------------------------#
# GSE176021_JustinaXCaushi2021_10x Genomics
if (F) {
  # The quality of cells was then assessed based on (1) the number  of genes detected per cell and
  # (2) the proportion of mitochondrial gene/ ribosomal gene counts.
  # Low-quality cells were filtered if the number  of detected genes was below 250 or
  # above 3× the median absolute  deviation away from the median gene number of all cells.
  # Cells were  filtered out if the proportion of mitochondrial gene counts was higher
  # than 10% or the proportion of ribosomal genes was less than 10%.

  ### clustering.R
  # gene_count <- as.matrix(JustinaXCaushi2021.MD01.010@assays$RNA@counts)
  # gene_count_norm <- sweep(gene_count,2,colSums(gene_count),FUN="/")*1e4 #把减去均值后的矩阵在列的方向上除以极差向量
  # gene_hypervar <- hypervar(gene_count_norm,showplot=FALSE)
  # gene_hypervar_sort <- gene_hypervar$data %>% arrange(.,desc(hypervar_log2))
  # VariableFeatures(JustinaXCaushi2021.MD01.010) <- gene_hypervar_sort$feature[1:3000]
  # VariableFeatures(JustinaXCaushi2021.MD01.010)<-setdiff(VariableFeatures(JustinaXCaushi2021.MD01.010),cluster.exclude)

  nFeature_lower <- 500
  nFeature_upper <- 4000
  nCount_lower <- 300
  nCount_upper <- 10000
  pMT_lower <- 0
  pMT_upper <- 15
  # pHB_lower <- 0
  # pHB_upper <- 5

  JustinaXCaushi2021.MD01.010 <- Read10X("/data/mengxu/data/GSE176021/MD01-010_tumor_1")
  JustinaXCaushi2021.MD01.010 <- CreateSeuratObject(
    counts = JustinaXCaushi2021.MD01.010,
    project = "JustinaXCaushi2021.MD01.010",
    min.features = 200,
    min.cells = 3
  )
  ### QC
  JustinaXCaushi2021.MD01.010 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.010, pattern = "^MT-", col.name = "pMT")
  JustinaXCaushi2021.MD01.010 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.010, pattern = "^HBA|^HBB", col.name = "pHB")
  JustinaXCaushi2021.MD01.010 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.010, pattern = "^RPS|^RPL", col.name = "pRP")
  JustinaXCaushi2021.MD01.010 <- subset(JustinaXCaushi2021.MD01.010,
    subset =
      nFeature_RNA > nFeature_lower &
        nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  JustinaXCaushi2021.MD01.005_list <- list()

  for (i in 2:9) {
    JustinaXCaushi2021.MD01.005_list[[i]] <- Read10X(paste("/data/mengxu/data/GSE176021/MD01-005_tumor_", i, sep = ""))

    JustinaXCaushi2021.MD01.005_list[[i]] <- CreateSeuratObject(
      counts = JustinaXCaushi2021.MD01.005_list[[i]],
      project = "JustinaXCaushi2021.MD01.005",
      min.features = 200,
      min.cells = 3
    )

    JustinaXCaushi2021.MD01.005_list[[i]] <- RenameCells(JustinaXCaushi2021.MD01.005_list[[i]],
      add.cell.id = "JustinaXCaushi2021.MD01.005"
    )
  }

  JustinaXCaushi2021.MD01.005 <- merge(JustinaXCaushi2021.MD01.005_list[[2]],
    y = c(
      JustinaXCaushi2021.MD01.005_list[[3]],
      JustinaXCaushi2021.MD01.005_list[[4]],
      JustinaXCaushi2021.MD01.005_list[[5]],
      JustinaXCaushi2021.MD01.005_list[[6]],
      JustinaXCaushi2021.MD01.005_list[[7]],
      JustinaXCaushi2021.MD01.005_list[[8]],
      JustinaXCaushi2021.MD01.005_list[[9]]
    )
  )
  ### QC
  JustinaXCaushi2021.MD01.005 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.005, pattern = "^MT-", col.name = "pMT")
  JustinaXCaushi2021.MD01.005 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.005, pattern = "^HBA|^HBB", col.name = "pHB")
  JustinaXCaushi2021.MD01.005 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.005, pattern = "^RPS|^RPL", col.name = "pRP")
  JustinaXCaushi2021.MD01.005 <- subset(JustinaXCaushi2021.MD01.005,
    subset =
      nFeature_RNA > nFeature_lower &
        nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  JustinaXCaushi2021.NY016.021 <- Read10X("/data/mengxu/data/GSE176021/NY016-021_tumor_1")

  JustinaXCaushi2021.NY016.021 <- CreateSeuratObject(
    counts = JustinaXCaushi2021.NY016.021,
    project = "JustinaXCaushi2021.NY016.021",
    min.features = 200,
    min.cells = 3
  )
  ### QC
  JustinaXCaushi2021.NY016.021 <- PercentageFeatureSet(JustinaXCaushi2021.NY016.021, pattern = "^MT-", col.name = "pMT")
  JustinaXCaushi2021.NY016.021 <- PercentageFeatureSet(JustinaXCaushi2021.NY016.021, pattern = "^HBA|^HBB", col.name = "pHB")
  JustinaXCaushi2021.NY016.021 <- PercentageFeatureSet(JustinaXCaushi2021.NY016.021, pattern = "^RPS|^RPL", col.name = "pRP")
  JustinaXCaushi2021.NY016.021 <- subset(JustinaXCaushi2021.NY016.021,
    subset =
      nFeature_RNA > nFeature_lower &
        nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
}

#----------------------------------------------------------------------------------------------#
# GSE99254_XinyiGuo2018_Smart-Seq2
if (F) {
  # Low-quality cells were discarded if the cell library size or the number of  expressed genes (counts larger than 0) was smaller than pre-defined thresholds,
  # which were the medians of all cells minus 3×median absolute deviation.
  # Cells  were also removed if their proportions of mitochondrial gene expression were  larger than 10%.

  nFeature_lower <- 201
  nFeature_upper <- 4000
  nCount_lower <- 500
  nCount_upper <- 1000000
  pMT_lower <- 0
  pMT_upper <- 10
  # pHB_lower <- 0
  # pHB_upper <- 5

  GSE99254 <- read.table("/data/mengxu/data/GSE99254/GSE99254_NSCLC.TCell.S12346.count.txt.gz",
    row.names = 1,
    header = T
  )
  GSE99254 <- GSE99254[!duplicated(GSE99254$symbol), ] # remove 88 features
  GSE99254 <- na.omit(GSE99254) # remove 1 feature
  row.names(GSE99254) <- GSE99254$symbol
  ###
  XinyiGuo2018.P0617 <- GSE99254[c(grep("0617", colnames(GSE99254)))]
  XinyiGuo2018.P0617 <- CreateSeuratObject(
    counts = XinyiGuo2018.P0617,
    project = "XinyiGuo2018.P0617",
    min.features = 200,
    min.cells = 3
  )
  XinyiGuo2018.P0617$platform <- "SmartSeq2"
  ### QC
  XinyiGuo2018.P0617 <- PercentageFeatureSet(XinyiGuo2018.P0617, pattern = "^MT-", col.name = "pMT")
  XinyiGuo2018.P0617 <- PercentageFeatureSet(XinyiGuo2018.P0617, pattern = "^HBA|^HBB", col.name = "pHB")
  XinyiGuo2018.P0617 <- PercentageFeatureSet(XinyiGuo2018.P0617, pattern = "^RPS|^RPL", col.name = "pRP")
  XinyiGuo2018.P0617 <- subset(XinyiGuo2018.P0617,
    subset =
      nFeature_RNA > nFeature_lower &
        nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  XinyiGuo2018.P0619 <- GSE99254[c(grep("0619", colnames(GSE99254)))]
  XinyiGuo2018.P0619 <- CreateSeuratObject(
    counts = XinyiGuo2018.P0619,
    project = "XinyiGuo2018.P0619",
    min.features = 200,
    min.cells = 3
  )
  XinyiGuo2018.P0619$platform <- "SmartSeq2"
  ### QC
  XinyiGuo2018.P0619 <- PercentageFeatureSet(XinyiGuo2018.P0619, pattern = "^MT-", col.name = "pMT")
  XinyiGuo2018.P0619 <- PercentageFeatureSet(XinyiGuo2018.P0619, pattern = "^HBA|^HBB", col.name = "pHB")
  XinyiGuo2018.P0619 <- PercentageFeatureSet(XinyiGuo2018.P0619, pattern = "^RPS|^RPL", col.name = "pRP")
  XinyiGuo2018.P0619 <- subset(XinyiGuo2018.P0619,
    subset =
      nFeature_RNA > nFeature_lower &
        nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  XinyiGuo2018.P0729 <- GSE99254[c(grep("0729", colnames(GSE99254)))]
  XinyiGuo2018.P0729 <- CreateSeuratObject(
    counts = XinyiGuo2018.P0729,
    project = "XinyiGuo2018.P0729",
    min.features = 200,
    min.cells = 3
  )
  XinyiGuo2018.P0729$platform <- "SmartSeq2"
  ### QC
  XinyiGuo2018.P0729 <- PercentageFeatureSet(XinyiGuo2018.P0729, pattern = "^MT-", col.name = "pMT")
  XinyiGuo2018.P0729 <- PercentageFeatureSet(XinyiGuo2018.P0729, pattern = "^HBA|^HBB", col.name = "pHB")
  XinyiGuo2018.P0729 <- PercentageFeatureSet(XinyiGuo2018.P0729, pattern = "^RPS|^RPL", col.name = "pRP")
  XinyiGuo2018.P0729 <- subset(XinyiGuo2018.P0729,
    subset =
      nFeature_RNA > nFeature_lower &
        nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  XinyiGuo2018.P0913 <- GSE99254[c(grep("0913", colnames(GSE99254)))]
  XinyiGuo2018.P0913 <- CreateSeuratObject(
    counts = XinyiGuo2018.P0913,
    project = "XinyiGuo2018.P0913",
    min.features = 200,
    min.cells = 3
  )
  XinyiGuo2018.P0913$platform <- "SmartSeq2"
  ### QC
  XinyiGuo2018.P0913 <- PercentageFeatureSet(XinyiGuo2018.P0913, pattern = "^MT-", col.name = "pMT")
  XinyiGuo2018.P0913 <- PercentageFeatureSet(XinyiGuo2018.P0913, pattern = "^HBA|^HBB", col.name = "pHB")
  XinyiGuo2018.P0913 <- PercentageFeatureSet(XinyiGuo2018.P0913, pattern = "^RPS|^RPL", col.name = "pRP")
  XinyiGuo2018.P0913 <- subset(XinyiGuo2018.P0913,
    subset =
      nFeature_RNA > nFeature_lower &
        nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
}

#----------------------------------------------------------------------------------------------#
# GSE131907_NayoungKim2020_10x Genomics
if (T) {
  # We applied three quality measures on raw gene-cell-barcode matrix for each cell:
  #   mitochondrial genes (≤20%, unique molecular identifiers (UMIs),
  #   and gene count (ranging from 100 to 150,000 and 200 to 10,000).
  #   The UMI count for the genes in each cell was lognormalized to TPM-like values,
  #   and then used in the log2 scale transcripts per million (TPM) plus 1.
  #   For each batch, we used the filtered cells to remove genes that are expressed at
  #   low levels by counting the number of cells (min.cells) having expression of each gene i,
  #   and excluded genes with min.cells < 0.1% cells.

  nFeature_lower <- 200
  nFeature_upper <- 10000
  nCount_lower <- 100
  nCount_upper <- 150000
  pMT_lower <- 0
  pMT_upper <- 20
  # pHB_lower <- 0
  # pHB_upper <- 5

  # GSE131907 <- fread(file = "/data/mengxu/data/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
  #                sep = '\t',
  #                header = T,
  #                check.names = F)
  # GSE131907 <- as.data.frame(GSE131907) #remove 0 feature GSE131907 <- as.data.frame(na.omit(GSE131907))
  # row.names(GSE131907) <- GSE131907$Index
  # NayoungKim2020.P0028 <- GSE131907[c(grep("_LUNG_T28",colnames(GSE131907)))]
  # NayoungKim2020.P0031 <- GSE131907[c(grep("_LUNG_T31",colnames(GSE131907)))]
  #
  # fwrite(NayoungKim2020.P0028,
  #        '/data/mengxu/data/GSE131907/data_processed/NayoungKim2020.P0028_LUNG_T28.txt',
  #        sep = '\t',
  #        row.names = T) #fwrite保存行名会丢失？
  # fwrite(NayoungKim2020.P0031,
  #        '/data/mengxu/data/GSE131907/data_processed/NayoungKim2020.P0031_LUNG_T31.txt',
  #        sep = '\t',
  #        row.names = T)

  # NayoungKim2020.P0028 <- fread(file = "/data/mengxu/data/GSE131907/data_processed/NayoungKim2020.P0028_LUNG_T28.txt",
  #                sep = '\t',
  #                header = T,
  #                check.names = F)

  NayoungKim2020.P0028 <- read.table("/data/mengxu/data/GSE131907/data_processed/NayoungKim2020.P0028_LUNG_T28.txt",
    row.names = 1,
    header = T
  )

  NayoungKim2020.P0028 <- CreateSeuratObject(
    counts = NayoungKim2020.P0028,
    project = "NayoungKim2020.P0028",
    min.features = 200,
    min.cells = 3
  )
  NayoungKim2020.P0028$platform <- "10X"
  ### QC
  NayoungKim2020.P0028 <- PercentageFeatureSet(NayoungKim2020.P0028, pattern = "^MT-", col.name = "pMT")
  NayoungKim2020.P0028 <- PercentageFeatureSet(NayoungKim2020.P0028, pattern = "^HBA|^HBB", col.name = "pHB")
  NayoungKim2020.P0028 <- PercentageFeatureSet(NayoungKim2020.P0028, pattern = "^RPS|^RPL", col.name = "pRP")
  NayoungKim2020.P0028 <- subset(NayoungKim2020.P0028,
    subset =
      nFeature_RNA > nFeature_lower &
        nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  NayoungKim2020.P0031 <- read.table("/data/mengxu/data/GSE131907/data_processed/NayoungKim2020.P0031_LUNG_T31.txt",
    row.names = 1,
    header = T
  )

  NayoungKim2020.P0031 <- CreateSeuratObject(
    counts = NayoungKim2020.P0031,
    project = "NayoungKim2020.P0031",
    min.features = 200,
    min.cells = 3
  )
  NayoungKim2020.P0031$platform <- "10X"
  ### QC
  NayoungKim2020.P0031 <- PercentageFeatureSet(NayoungKim2020.P0031, pattern = "^MT-", col.name = "pMT")
  NayoungKim2020.P0031 <- PercentageFeatureSet(NayoungKim2020.P0031, pattern = "^HBA|^HBB", col.name = "pHB")
  NayoungKim2020.P0031 <- PercentageFeatureSet(NayoungKim2020.P0031, pattern = "^RPS|^RPL", col.name = "pRP")
  NayoungKim2020.P0031 <- subset(NayoungKim2020.P0031,
    subset =
      nFeature_RNA > nFeature_lower &
        nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
}

#----------------------------------------------------------------------------------------------#
# GSE154826_AndrewMLeader2021_10x Genomics
if (T) {
  nFeature_lower <- 300
  # nFeature_upper <- 3000
  nCount_lower <- 500
  # nCount_upper <- 10000
  pMT_lower <- 0
  pMT_upper <- 25
  # pHB_lower <- 0
  # pHB_upper <- 5

  patient_ids <- c("P377", "P403", "P514")
  file_ids <- c("6", "11", "29")

  AndrewMLeader2021 <- list()
  AndrewMLeader2021_samples <- c()
  for (i in 1:length(file_ids)) {
    patient <- patient_ids[i]
    file <- file_ids[i]
    AndrewMLeader2021_sample <- paste("AndrewMLeader2021", patient, file, sep = ".")
    message("[", Sys.time(), "]", " ------ Now: ", AndrewMLeader2021_sample)

    matrix <- Matrix::readMM(paste("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_", file, "/matrix.mtx", sep = ""))
    # feature <- read.table(paste("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_",file,"/features.tsv",sep = '')) #Bug
    feature <- read.delim(paste("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_", file, "/features.tsv", sep = ""), stringsAsFactors = F, header = F)
    barcode <- read.table(paste("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_", file, "/barcodes.tsv", sep = ""))
    # barcode <- read.delim(paste("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_",file,"/barcodes.tsv",sep = ''),stringsAsFactors=F,header=F)[,1]
    rownames(matrix) <- feature[, 2]
    colnames(matrix) <- barcode[, 1]

    AndrewMLeader2021_seu <- CreateSeuratObject(
      counts = matrix,
      project = AndrewMLeader2021_sample,
      min.features = 200,
      min.cells = 3
    )
    AndrewMLeader2021_seu$platform <- "10X"

    ### QC
    AndrewMLeader2021_seu <- PercentageFeatureSet(AndrewMLeader2021_seu, pattern = "^MT-", col.name = "pMT")
    AndrewMLeader2021_seu <- PercentageFeatureSet(AndrewMLeader2021_seu, pattern = "^HBA|^HBB", col.name = "pHB")
    AndrewMLeader2021_seu <- PercentageFeatureSet(AndrewMLeader2021_seu, pattern = "^RPS|^RPL", col.name = "pRP")

    AndrewMLeader2021_seu <- subset(AndrewMLeader2021_seu,
      subset =
        nFeature_RNA > nFeature_lower &
          # nFeature_RNA < nFeature_upper &
          nCount_RNA > nCount_lower &
          # nCount_RNA < nCount_upper &
          pMT < pMT_upper
    )

    AndrewMLeader2021[[i]] <- AndrewMLeader2021_seu
    AndrewMLeader2021_samples[i] <- AndrewMLeader2021_sample
  }
  # matrix <- Matrix::readMM("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_6/matrix.mtx")
  # barcode <- read.table('/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_6/barcodes.tsv')
  # feature <- read.table('/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_6/features.tsv')
  # rownames(matrix) <- feature[,2]
  # colnames(matrix) <- barcode[,1]
  # AndrewMLeader2021.P377.6 <- CreateSeuratObject(counts = matrix,
  #                                                project = 'AndrewMLeader2021.P377.6',
  #                                                min.features = 200,
  #                                                min.cells = 3)
  # AndrewMLeader2021.P377.6$platform <- '10X'
  # ###QC
  # AndrewMLeader2021.P377.6 <- PercentageFeatureSet(AndrewMLeader2021.P377.6, pattern = "^MT-", col.name = "pMT")
  # AndrewMLeader2021.P377.6 <- PercentageFeatureSet(AndrewMLeader2021.P377.6, pattern = "^HBA|^HBB", col.name = "pHB")
  # AndrewMLeader2021.P377.6 <- PercentageFeatureSet(AndrewMLeader2021.P377.6, pattern = "^RPS|^RPL", col.name = "pRP")
  # AndrewMLeader2021.P377.6 <- subset(AndrewMLeader2021.P377.6,
  #                                    subset =
  #                                      nFeature_RNA > nFeature_lower &
  #                                      #nFeature_RNA < nFeature_upper &
  #                                      nCount_RNA > nCount_lower &
  #                                      #nCount_RNA < nCount_upper &
  #                                      pMT < pMT_upper)
}

#----------------------------------------------------------------------------------------------#
# PRJNA591860_AshleyMaynard2020_Smart-Seq2
if (T) {
  # Standard procedures for filtering, variable gene selection, dimensionality reduction,
  # and clustering were performed using the Seurat v3 in RStudio using R ,
  # where cells with fewer than 500 genes and 50,000 reads were excluded.
  # We used DoubletFinder to identify potentially sorted doublet cells.
  # 218 doublets were excluded from further analysis.

  nFeature_lower <- 500
  # nFeature_upper <- 3000
  nCount_lower <- 50000
  # nCount_upper <- 50000
  # pMT_lower <- 0
  # pMT_upper <- 15
  # pHB_lower <- 0
  # pHB_upper <- 5

  if (Sys.info()[1] == "Windows") {
    # Metadata
    neo_osi_metadata <- read.csv("/\\192.168.0.109/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/neo-osi_metadata.csv")
    S01_metacells <- read.csv("/\\192.168.0.109/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/S01_metacells.csv")
    # Dataset
    # load("/data/mengxu/data/AshleyMaynard2020/Data_input/objects/S02.1_Main_Seurat_object_filtered_neo_osi.RData")
    # AshleyMaynard2020 <- read.csv('/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/S01_datafinal.csv')
    neo_osi_rawdata <- fread(
      file = "/\\192.168.0.109/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/neo-osi_rawdata.csv",
      sep = ",",
      header = T,
      check.names = F
    )
    neo_osi_rawdata <- as.data.frame(na.omit(neo_osi_rawdata))
    rownames(neo_osi_rawdata) <- neo_osi_rawdata$gene
    ###
    S01_datafinal <- fread(
      file = "/\\192.168.0.109/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/S01_datafinal.csv",
      sep = ",",
      header = T,
      check.names = F
    )
    S01_datafinal <- as.data.frame(na.omit(S01_datafinal))
    rownames(S01_datafinal) <- S01_datafinal$V1
  } else {
    # Metadata
    neo_osi_metadata <- read.csv("/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/neo-osi_metadata.csv")
    S01_metacells <- read.csv("/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/S01_metacells.csv")
    # Dataset
    neo_osi_rawdata <- fread(
      file = "/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/neo-osi_rawdata.csv",
      sep = ",",
      header = T,
      check.names = F
    )
    neo_osi_rawdata <- as.data.frame(na.omit(neo_osi_rawdata))
    rownames(neo_osi_rawdata) <- neo_osi_rawdata$gene
    ###
    S01_datafinal <- fread(
      file = "/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/S01_datafinal.csv",
      sep = ",",
      header = T,
      check.names = F
    )
    S01_datafinal <- as.data.frame(na.omit(S01_datafinal))
    rownames(S01_datafinal) <- S01_datafinal$V1
  }

  AshleyMaynard2020 <- list()
  AshleyMaynard2020_samples <- c()
  AshleyMaynard2020_id <- c("AZ_01", "AZ_03", "AZ_05", "AZ_04", "LT_S49", "LT_S69", "LT_S74")

  for (i in 1:length(AshleyMaynard2020_id)) {
    AshleyMaynard2020_seu <- c()

    if (AshleyMaynard2020_id[i] %in% neo_osi_metadata$sample_name) {
      AshleyMaynard2020_file <- neo_osi_metadata$plate[which(neo_osi_metadata$sample_name == AshleyMaynard2020_id[i], )]
      AshleyMaynard2020_file <- AshleyMaynard2020_file[!duplicated(AshleyMaynard2020_file)]

      AshleyMaynard2020_seu <- as.data.frame(neo_osi_rawdata[, 1])
      for (cell in AshleyMaynard2020_file) {
        print(cell)
        AshleyMaynard2020_dgc <- as.data.frame(neo_osi_rawdata[c(grep(cell, colnames(neo_osi_rawdata)))])
        AshleyMaynard2020_seu <- cbind.data.frame(
          AshleyMaynard2020_seu,
          AshleyMaynard2020_dgc
        )
      }
    } else {
      AshleyMaynard2020_file <- S01_metacells$plate[which(S01_metacells$sample_name == AshleyMaynard2020_id[i], )]
      AshleyMaynard2020_file <- AshleyMaynard2020_file[!duplicated(AshleyMaynard2020_file)]

      AshleyMaynard2020_seu <- as.data.frame(S01_datafinal[, 1])
      for (cell in AshleyMaynard2020_file) {
        print(cell)
        AshleyMaynard2020_dgc <- as.data.frame(S01_datafinal[c(grep(cell, colnames(S01_datafinal)))])
        AshleyMaynard2020_seu <- cbind.data.frame(
          AshleyMaynard2020_seu,
          AshleyMaynard2020_dgc
        )
      }
    }

    AshleyMaynard2020_seu <- AshleyMaynard2020_seu[, -1]
    AshleyMaynard2020_sample <- paste0("AshleyMaynard2020.", AshleyMaynard2020_id[i])

    ###
    AshleyMaynard2020_seu <- CreateSeuratObject(
      counts = AshleyMaynard2020_seu,
      project = AshleyMaynard2020_sample,
      min.features = 200,
      min.cells = 3
    )
    AshleyMaynard2020_seu$platform <- "SmartSeq2"

    ### QC
    AshleyMaynard2020_seu <- PercentageFeatureSet(AshleyMaynard2020_seu, pattern = "^MT-", col.name = "pMT")
    AshleyMaynard2020_seu <- PercentageFeatureSet(AshleyMaynard2020_seu, pattern = "^HBA|^HBB", col.name = "pHB")
    AshleyMaynard2020_seu <- PercentageFeatureSet(AshleyMaynard2020_seu, pattern = "^RPS|^RPL", col.name = "pRP")
    AshleyMaynard2020_seu <- subset(AshleyMaynard2020_seu,
      subset =
        nFeature_RNA > nFeature_lower &
          # nFeature_RNA < nFeature_upper &
          nCount_RNA > nCount_lower # &
      # nCount_RNA < nCount_upper &
      # pMT < pMT_upper
    )

    ###
    AshleyMaynard2020[[i]] <- AshleyMaynard2020_seu
    AshleyMaynard2020_samples[i] <- AshleyMaynard2020_sample
  }
}

#----------------------------------------------------------------------------------------------#

samples <- c(
  "PhilipBischoff2021.P030",
  "PhilipBischoff2021.P033",
  "DietherLambrechts2018.P06",
  "DietherLambrechts2018.P03",
  "RapolasZilionis2019.P2",
  "RapolasZilionis2019.P7",
  "NayoungKim2020.P0028",
  "NayoungKim2020.P0031",
  AndrewMLeader2021_samples,
  AshleyMaynard2020_samples
)

#----------------------------------------------------------------------------------------#
seu_obj_list <- c(
  PhilipBischoff2021.P030,
  PhilipBischoff2021.P033,
  DietherLambrechts2018.P06,
  DietherLambrechts2018.P03,
  RapolasZilionis2019.P2,
  RapolasZilionis2019.P7,
  NayoungKim2020.P0028,
  NayoungKim2020.P0031,
  AndrewMLeader2021,
  AshleyMaynard2020
)

save(seu_obj_list, samples, file = "/data/mengxu/data/all/lung_stage-3_list_raw.Rdata")

#----------------------------------------------------------------------------------------#
rm(list = ls())
gc()
