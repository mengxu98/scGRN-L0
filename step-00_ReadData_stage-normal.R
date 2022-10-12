

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
# PhilipBischoff2021_10X Genomics
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

  PhilipBischoff2021 <- list()
  PhilipBischoff2021_samples <- c()
  PhilipBischoff2021_id <- c("18", "19", "27", "28", "29", "30", "31", "32", "33", "34")
  for (i in 1:length(PhilipBischoff2021_id)) {
    PhilipBischoff2021_sample <- paste("PhilipBischoff2021.P0", PhilipBischoff2021_id[i], sep = "")

    PhilipBischoff2021_seu <- Read10X(paste("/data/mengxu/data/PhilipBischoff2021/cellranger/p0", PhilipBischoff2021_id[i], "n/filtered_feature_bc_matrix", sep = ""))

    PhilipBischoff2021_seu <- CreateSeuratObject(
      counts = PhilipBischoff2021_seu,
      project = PhilipBischoff2021_sample,
      min.features = 200,
      min.cells = 3
    )

    PhilipBischoff2021_seu$platform <- "10X"
    ### QC
    PhilipBischoff2021_seu <- PercentageFeatureSet(PhilipBischoff2021_seu, pattern = "^MT-", col.name = "pMT")
    PhilipBischoff2021_seu <- PercentageFeatureSet(PhilipBischoff2021_seu, pattern = "^HBA|^HBB", col.name = "pHB")
    PhilipBischoff2021_seu <- PercentageFeatureSet(PhilipBischoff2021_seu, pattern = "^RPS|^RPL", col.name = "pRP")
    PhilipBischoff2021_seu <- subset(PhilipBischoff2021_seu,
      subset = nFeature_RNA > nFeature_lower &
        nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper &
        pHB < pHB_upper
    )

    ###
    PhilipBischoff2021[[i]] <- PhilipBischoff2021_seu
    PhilipBischoff2021_samples[i] <- PhilipBischoff2021_sample
  }
}

#----------------------------------------------------------------------------------------------#
# CRA001963_DiHe2021_10x Genomics
if (F) {
  # Cells with low feature counts (<200) and high percent of mitochondrial genes (>10%) were removed.

  # nFeature_lower <- 201
  # nFeature_upper <- 10000
  nCount_lower <- 200
  # nCount_upper <- 6000
  pMT_lower <- 10
  # pMT_upper <- 10
  # pHB_lower <- 0
  # pHB_upper <- 5

  DiHe2021 <- list()
  DiHe2021_samples <- c()
  DiHe2021_id <- c("N1", "N2", "N3", "N4", "N5")
  DiHe2021_file <- c("CRR073027", "CRR073028", "CRR073029", "CRR073030", "CRR073031")
  for (i in 1:length(DiHe2021_id)) {
    DiHe2021_sample <- paste("DiHe2021.", DiHe2021_id[i], sep = "")
    DiHe2021_seu <- Read10X(paste("/data/mengxu/data/CRA001963/", DiHe2021_file[i], "/", DiHe2021_file[i], "_f", "/outs/filtered_feature_bc_matrix", sep = ""))

    DiHe2021_seu <- CreateSeuratObject(
      counts = DiHe2021_seu,
      project = DiHe2021_sample,
      min.features = 200,
      min.cells = 3
    )

    DiHe2021_seu$platform <- "10X"
    ### QC
    DiHe2021_seu <- PercentageFeatureSet(DiHe2021_seu, pattern = "^MT-", col.name = "pMT")
    DiHe2021_seu <- PercentageFeatureSet(DiHe2021_seu, pattern = "^HBA|^HBB", col.name = "pHB")
    DiHe2021_seu <- PercentageFeatureSet(DiHe2021_seu, pattern = "^RPS|^RPL", col.name = "pRP")
    DiHe2021_seu <- subset(DiHe2021_seu,
      subset =
      # nFeature_RNA > nFeature_lower &
      # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
          # nCount_RNA < nCount_upper &
          # pMT < pMT_upper&
          pMT > pMT_lower
    )
    ###
    DiHe2021[[i]] <- DiHe2021_seu
    DiHe2021_samples[i] <- DiHe2021_sample
  }
}

#------------------------------------------------------------------------------------------------------#
# E-MTAB-6653_DietherLambrechts2018_10X Genomics
if (F) {

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

  DietherLambrechts2018.P06 <- Read10X("/data/mengxu/data/E-MTAB-6653/scrBT1432m/scrBT1432m_f/outs/filtered_feature_bc_matrix")

  DietherLambrechts2018.P06 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P06,
    project = "DietherLambrechts2018.P06",
    min.features = 200,
    min.cells = 3
  )

  DietherLambrechts2018.P06$platform <- "10X"
  ### QC
  DietherLambrechts2018.P06 <- PercentageFeatureSet(DietherLambrechts2018.P06, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P06 <- PercentageFeatureSet(DietherLambrechts2018.P06, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P06 <- PercentageFeatureSet(DietherLambrechts2018.P06, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P06 <- subset(DietherLambrechts2018.P06,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )

  ###
  DietherLambrechts2018.P07 <- Read10X("/data/mengxu/data/E-MTAB-6653/BT1378/BT1378_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P07 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P07,
    project = "DietherLambrechts2018.P07",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P07$platform <- "10X"
  ### QC
  DietherLambrechts2018.P07 <- PercentageFeatureSet(DietherLambrechts2018.P07, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P07 <- PercentageFeatureSet(DietherLambrechts2018.P07, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P07 <- PercentageFeatureSet(DietherLambrechts2018.P07, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P07 <- subset(DietherLambrechts2018.P07,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )

  ###
  DietherLambrechts2018.P08 <- Read10X("/data/mengxu/data/E-MTAB-6653/scrBT1428m/scrBT1428m_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P08 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P08,
    project = "DietherLambrechts2018.P08",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P08$platform <- "10X"
  ### QC
  DietherLambrechts2018.P08 <- PercentageFeatureSet(DietherLambrechts2018.P08, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P08 <- PercentageFeatureSet(DietherLambrechts2018.P08, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P08 <- PercentageFeatureSet(DietherLambrechts2018.P08, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P08 <- subset(DietherLambrechts2018.P08,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
}

# E-MTAB-6149_DietherLambrechts2018_10X Genomics
if (F) {
  DietherLambrechts2018.P02 <- Read10X("/data/mengxu/data/E-MTAB-6149/1247/1247/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P02 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P02,
    project = "DietherLambrechts2018.P02",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P02$platform <- "10X"
  ### QC
  DietherLambrechts2018.P02 <- PercentageFeatureSet(DietherLambrechts2018.P02, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P02 <- PercentageFeatureSet(DietherLambrechts2018.P02, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P02 <- PercentageFeatureSet(DietherLambrechts2018.P02, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P02 <- subset(DietherLambrechts2018.P02,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P03 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT1293/BT1293_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P03 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P03,
    project = "DietherLambrechts2018.P03",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P03$platform <- "10X"
  ### QC
  DietherLambrechts2018.P03 <- PercentageFeatureSet(DietherLambrechts2018.P03, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P03 <- PercentageFeatureSet(DietherLambrechts2018.P03, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P03 <- PercentageFeatureSet(DietherLambrechts2018.P03, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P03 <- subset(DietherLambrechts2018.P03,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P04 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT1297/BT1297_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P04 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P04,
    project = "DietherLambrechts2018.P04",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P04$platform <- "10X"
  ### QC
  DietherLambrechts2018.P04 <- PercentageFeatureSet(DietherLambrechts2018.P04, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P04 <- PercentageFeatureSet(DietherLambrechts2018.P04, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P04 <- PercentageFeatureSet(DietherLambrechts2018.P04, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P04 <- subset(DietherLambrechts2018.P04,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )

  ###
  DietherLambrechts2018.P05 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT1301/BT1301_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P05 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P05,
    project = "DietherLambrechts2018.P05",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P05$platform <- "10X"
  ### QC
  DietherLambrechts2018.P05 <- PercentageFeatureSet(DietherLambrechts2018.P05, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P05 <- PercentageFeatureSet(DietherLambrechts2018.P05, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P05 <- PercentageFeatureSet(DietherLambrechts2018.P05, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P05 <- subset(DietherLambrechts2018.P05,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
}
#----------------------------------------------------------------------------------------------#
# GSE117570_QianqianSong2019_10x Genomics
if (F) {
  # Low‐quality cells were discarded if the cell number with  expressed genes was smaller than 200.
  # Cells were also removed if their proportions of mitochondrial gene expression  were larger than 40%.

  nFeature_lower <- 200
  # nFeature_upper <- 10000
  # nCount_lower <- 200
  # nCount_upper <- 6000
  pMT_lower <- 40
  # pMT_upper <- 10
  # pHB_lower <- 0
  # pHB_upper <- 5
  ###
  QianqianSong2019.P1 <- read.table("/data/mengxu/data/GSE117570/GSM3304008_P1_Normal_processed_data.txt.gz",
    row.names = 1,
    header = T
  )

  QianqianSong2019.P1 <- CreateSeuratObject(
    counts = QianqianSong2019.P1,
    project = "QianqianSong2019.P1",
    min.features = 200,
    min.cells = 3
  )
  QianqianSong2019.P1$platform <- "10X"
  ### QC
  QianqianSong2019.P1 <- PercentageFeatureSet(QianqianSong2019.P1, pattern = "^MT-", col.name = "pMT")
  QianqianSong2019.P1 <- PercentageFeatureSet(QianqianSong2019.P1, pattern = "^HBA|^HBB", col.name = "pHB")
  QianqianSong2019.P1 <- PercentageFeatureSet(QianqianSong2019.P1, pattern = "^RPS|^RPL", col.name = "pRP")
  QianqianSong2019.P1 <- subset(QianqianSong2019.P1,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        # nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        # pMT < pMT_upper&
        pMT > pMT_lower
  )

  ###
  QianqianSong2019.P2 <- read.table("/data/mengxu/data/GSE117570/GSM3304010_P2_Normal_processed_data.txt.gz",
    row.names = 1,
    header = T
  )

  QianqianSong2019.P2 <- CreateSeuratObject(
    counts = QianqianSong2019.P2,
    project = "QianqianSong2019.P2",
    min.features = 200,
    min.cells = 3
  )
  QianqianSong2019.P2$platform <- "10X"
  ### QC
  QianqianSong2019.P2 <- PercentageFeatureSet(QianqianSong2019.P2, pattern = "^MT-", col.name = "pMT")
  QianqianSong2019.P2 <- PercentageFeatureSet(QianqianSong2019.P2, pattern = "^HBA|^HBB", col.name = "pHB")
  QianqianSong2019.P2 <- PercentageFeatureSet(QianqianSong2019.P2, pattern = "^RPS|^RPL", col.name = "pRP")
  QianqianSong2019.P2 <- subset(QianqianSong2019.P2,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        # nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        # pMT < pMT_upper&
        pMT > pMT_lower
  )
  ###

  QianqianSong2019.P3 <- read.table("/data/mengxu/data/GSE117570/GSM3304012_P3_Normal_processed_data.txt.gz",
    row.names = 1,
    header = T
  )

  QianqianSong2019.P3 <- CreateSeuratObject(
    counts = QianqianSong2019.P3,
    project = "QianqianSong2019.P3",
    min.features = 200,
    min.cells = 3
  )
  QianqianSong2019.P3$platform <- "10X"
  ### QC
  QianqianSong2019.P3 <- PercentageFeatureSet(QianqianSong2019.P3, pattern = "^MT-", col.name = "pMT")
  QianqianSong2019.P3 <- PercentageFeatureSet(QianqianSong2019.P3, pattern = "^HBA|^HBB", col.name = "pHB")
  QianqianSong2019.P3 <- PercentageFeatureSet(QianqianSong2019.P3, pattern = "^RPS|^RPL", col.name = "pRP")
  QianqianSong2019.P3 <- subset(QianqianSong2019.P3,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        # nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        # pMT < pMT_upper&
        pMT > pMT_lower
  )

  ###
  QianqianSong2019.P4 <- read.table("/data/mengxu/data/GSE117570/GSM3304014_P4_Normal_processed_data.txt.gz",
    row.names = 1,
    header = T
  )

  QianqianSong2019.P4 <- CreateSeuratObject(
    counts = QianqianSong2019.P4,
    project = "QianqianSong2019.P4",
    min.features = 200,
    min.cells = 3
  )
  QianqianSong2019.P4$platform <- "10X"
  ### QC
  QianqianSong2019.P4 <- PercentageFeatureSet(QianqianSong2019.P4, pattern = "^MT-", col.name = "pMT")
  QianqianSong2019.P4 <- PercentageFeatureSet(QianqianSong2019.P4, pattern = "^HBA|^HBB", col.name = "pHB")
  QianqianSong2019.P4 <- PercentageFeatureSet(QianqianSong2019.P4, pattern = "^RPS|^RPL", col.name = "pRP")
  QianqianSong2019.P4 <- subset(QianqianSong2019.P4,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        # nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        # pMT < pMT_upper&
        pMT > pMT_lower
  )
}

#----------------------------------------------------------------------------------------------#
# GSE123902_AshleyMLaughney2020_10x Genomics
if (F) {

  # Cells with >20% of transcripts derived from mitochondria were considered  apoptotic and also excluded (Extended Data Fig. 2c).
  # Afer clustering (described  below), two minority cell clusters were additionally excluded from patient data  with characteristically
  # low library size and signatures of partial cell lysis that  escaped the automated flters described above.
  # Genes detected in fewer than ten  cells or genes with low expression levels,
  # identifed as those with count values <5  s.d. from the second mode of the log–log distribution of total transcript counts and  genes,
  # were also excluded. Tis yielded a total of 40,505 patient-derived cells with  a median library size of 4,038 transcripts per cell (Extended Data Fig. 2d,e),
  # for  downstream analysis (Fig. 1 and Extended Data Fig. 2).

  # nFeature_lower <- 200
  # nFeature_upper <- 10000
  # nCount_lower <- 200
  # nCount_upper <- 6000
  # pMT_lower <- 20
  # pMT_upper <- 10
  # pHB_lower <- 0
  # pHB_upper <- 5

  ###
  AshleyMLaughney2020 <- list()
  AshleyMLaughney2020_samples <- c()
  AshleyMLaughney2020_id <- c("LX675", "LX682", "LX684", "LX685")
  AshleyMLaughney2020_file <- c(
    "GSM3516666_MSK_LX675_NORMAL", "GSM3516673_MSK_LX682_NORMAL",
    "GSM3516675_MSK_LX684_NORMAL", "GSM3516676_MSK_LX685_NORMAL"
  )

  for (i in length(AshleyMLaughney2020_id)) {
    AshleyMLaughney2020_sample <- paste0("AshleyMLaughney2020.", AshleyMLaughney2020_id[i])
    AshleyMLaughney2020_seu <- read.table(paste0("/data/mengxu/data/GSE123902/", AshleyMLaughney2020_file[i], "_dense.csv.gz"),
      row.names = 1,
      header = T,
      sep = ","
    )

    # head(AshleyMLaughney2020_seu[1:4,1:4])
    AshleyMLaughney2020_seu <- t.data.frame(AshleyMLaughney2020_seu)

    AshleyMLaughney2020_seu <- CreateSeuratObject(
      counts = AshleyMLaughney2020_seu,
      project = AshleyMLaughney2020_sample,
      min.features = 200,
      min.cells = 3
    )
    AshleyMLaughney2020_seu$platform <- "10X"
    ### QC
    AshleyMLaughney2020_seu <- PercentageFeatureSet(AshleyMLaughney2020_seu, pattern = "^MT-", col.name = "pMT")
    AshleyMLaughney2020_seu <- PercentageFeatureSet(AshleyMLaughney2020_seu, pattern = "^HBA|^HBB", col.name = "pHB")
    AshleyMLaughney2020_seu <- PercentageFeatureSet(AshleyMLaughney2020_seu, pattern = "^RPS|^RPL", col.name = "pRP")
    # AshleyMLaughney2020_seu <- subset(AshleyMLaughney2020_seu,
    #                                   subset =
    #                                     nFeature_RNA > nFeature_lower &
    #                                     #nFeature_RNA < nFeature_upper &
    #                                     #nCount_RNA > nCount_lower &
    #                                     #nCount_RNA < nCount_upper &
    #                                     #pMT < pMT_upper&
    #                                     pMT > pMT_lower)

    AshleyMLaughney2020[[i]] <- AshleyMLaughney2020_seu
    AshleyMLaughney2020_samples[i] <- AshleyMLaughney2020_sample
  }
}

#------------------------------------------------------------------------------------------------------------------#
# GSE126030_PeterASzabo2019_10x Genomics

# PeterASzabo2019.P1	GSM3589406_Tissue donor 1
# PeterASzabo2019.P1	GSM3589407_Tissue donor 1
# PeterASzabo2019.P2	GSM3589412_Tissue donor 2
# PeterASzabo2019.P2	GSM3589413_Tissue donor 2

#------------------------------------------------------------------------------------------------------------------#

# GSE131907_NayoungKim2020_10x Genomics
# 读取数据较麻烦，且该阶段样本足够
if (F) {
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
  # NayoungKim2020.P0009 <- GSE131907[c(grep("_LUNG_T09",colnames(GSE131907)))]
  #
  # fwrite(NayoungKim2020.P0009,
  #        '/data/mengxu/data/GSE131907/data_processed/NayoungKim2020.P0009_LUNG_T09.txt',
  #        sep = '\t',
  #        row.names = T) #fwrite保存行名会丢失？

  NayoungKim2020.P0009 <- read.table("/data/mengxu/data/GSE131907/data_processed/NayoungKim2020.P0009_LUNG_T09.txt",
    row.names = 1,
    header = T
  )

  NayoungKim2020.P0009 <- CreateSeuratObject(
    counts = NayoungKim2020.P0009,
    project = "NayoungKim2020.P0009",
    min.features = 200,
    min.cells = 3
  )
  NayoungKim2020.P0009$platform <- "10X"
  ### QC
  NayoungKim2020.P0009 <- PercentageFeatureSet(NayoungKim2020.P0009, pattern = "^MT-", col.name = "pMT")
  NayoungKim2020.P0009 <- PercentageFeatureSet(NayoungKim2020.P0009, pattern = "^HBA|^HBB", col.name = "pHB")
  NayoungKim2020.P0009 <- PercentageFeatureSet(NayoungKim2020.P0009, pattern = "^RPS|^RPL", col.name = "pRP")
  NayoungKim2020.P0009 <- subset(NayoungKim2020.P0009,
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
if (F) {
  nFeature_lower <- 300
  # nFeature_upper <- 3000
  nCount_lower <- 500
  # nCount_upper <- 10000
  pMT_lower <- 0
  pMT_upper <- 25
  # pHB_lower <- 0
  # pHB_upper <- 5


  patient_ids <- c(
    "P370", "P371", "P377", "P378", "P403", "P406", "P408", "P410", "P458", "P460",
    "P464", "P514", "P522", "P532", "P564", "P569", "P570", "P572", "P578", "P581",
    "P581", "P584", "P584", "P593", "P596", "P596", "P630", "P626", "P626", "P695",
    "P695", "P695", "P695", "P695", "P522", "P714", "P729", "P706", "P706", "P706",
    "P706", "P706", "P706", "P706", "P706"
  )
  file_ids <- c(
    "1", "3", "5", "7", "10", "14", "16", "18", "22", "24",
    "26", "28", "30", "32", "37", "39", "41", "44", "46", "48",
    "49", "51", "52", "53", "54", "55", "56", "88", "89", "90",
    "91", "92", "93", "94", "244", "293", "307", "342", "343", "344",
    "345", "346", "347", "348", "350"
  )

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
}

#----------------------------------------------------------------------------------------------#
# PRJNA591860_AshleyMaynard2020_Smart-Seq2
if (F) {
  # Standard procedures for filtering, variable gene selection, dimensionality reduction,
  # and clustering were performed using the Seurat v3 in RStudio using R ,
  # where cells with fewer than 500 genes and 50,000 reads were excluded.
  # We used DoubletFinder to identify potentially sorted doublet cells.
  # 218 doublets were excluded from further analysis.

  nFeature_lower <- 500
  # nFeature_upper <- 3000
  nCount_lower <- 50000
  # nCount_upper <- 50000
  pMT_lower <- 0
  pMT_upper <- 15
  # pHB_lower <- 0
  # pHB_upper <- 5

  # load("/data/mengxu/data/AshleyMaynard2020/Data_input/objects/S02.1_Main_Seurat_object_filtered_neo_osi.RData")
  # AshleyMaynard2020 <- read.csv('/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/S01_datafinal.csv')
  AshleyMaynard2020_part1 <- fread(
    file = "/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/neo-osi_rawdata.csv",
    sep = ",",
    header = T,
    check.names = F
  )
  AshleyMaynard2020_part1 <- as.data.frame(na.omit(AshleyMaynard2020_part1))
  rownames(AshleyMaynard2020_part1) <- AshleyMaynard2020_part1$gene
  ###
  AshleyMaynard2020_part2 <- fread(
    file = "/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/S01_datafinal.csv",
    sep = ",",
    header = T,
    check.names = F
  )
  AshleyMaynard2020_part2 <- as.data.frame(na.omit(AshleyMaynard2020_part2))
  rownames(AshleyMaynard2020_part2) <- AshleyMaynard2020_part2$V1
  ######
  ### AshleyMaynard2020.AZ.01
  AshleyMaynard2020.AZ.01_1 <- AshleyMaynard2020_part1[c(grep("B001544", colnames(AshleyMaynard2020_part1)))]
  AshleyMaynard2020.AZ.01_2 <- AshleyMaynard2020_part1[c(grep("B001567", colnames(AshleyMaynard2020_part1)))] # 0 cell matched
  AshleyMaynard2020.AZ.01_3 <- AshleyMaynard2020_part1[c(grep("B001546", colnames(AshleyMaynard2020_part1)))]
  AshleyMaynard2020.AZ.01_4 <- AshleyMaynard2020_part1[c(grep("B001481", colnames(AshleyMaynard2020_part1)))] # 0 cell matched
  AshleyMaynard2020.AZ.01 <- cbind.data.frame(
    AshleyMaynard2020.AZ.01_1,
    AshleyMaynard2020.AZ.01_2,
    AshleyMaynard2020.AZ.01_3,
    AshleyMaynard2020.AZ.01_4
  )
  ###
  AshleyMaynard2020.AZ.01 <- CreateSeuratObject(
    counts = AshleyMaynard2020.AZ.01,
    project = "AshleyMaynard2020.AZ.01",
    min.features = 200,
    min.cells = 3
  )
  AshleyMaynard2020.AZ.01$platform <- "SmartSeq2"

  ### QC
  AshleyMaynard2020.AZ.01 <- PercentageFeatureSet(AshleyMaynard2020.AZ.01, pattern = "^MT-", col.name = "pMT")
  AshleyMaynard2020.AZ.01 <- PercentageFeatureSet(AshleyMaynard2020.AZ.01, pattern = "^HBA|^HBB", col.name = "pHB")
  AshleyMaynard2020.AZ.01 <- PercentageFeatureSet(AshleyMaynard2020.AZ.01, pattern = "^RPS|^RPL", col.name = "pRP")
  AshleyMaynard2020.AZ.01 <- subset(AshleyMaynard2020.AZ.01,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
}

#----------------------------------------------------------------------------------------------#
# PRJNA773987_EunYoungKim2022_10x Genomics
if (F) {
  # Filtering of low quality data was performed using the percentage of mitochondrial genes (< 20%),
  # gene counts (200 ~ 10,000), and UMI counts (100 ~ 150,000) to remove noise from empty droplets, debris, doublets and triplets [4-6].
  # Data from independent experiments were first integrated to minimize batch effects and then clustered into individual lung cell types by principal component analysis.
  # Each cell cluster was isolated and analyzed independently for additional subclustering.

  nFeature_lower <- 200
  nFeature_upper <- 10000
  nCount_lower <- 100
  nCount_upper <- 150000
  pMT_lower <- 0
  pMT_upper <- 20
  # pHB_lower <- 0
  # pHB_upper <- 5

  EunYoungKim2022 <- list()
  EunYoungKim2022_samples <- c()
  EunYoungKim2022_id <- c("05", "04", "03", "02", "06")
  EunYoungKim2022_file <- c("SRR16668029", "SRR16668031", "SRR16668033", "SRR16668035", "SRR16668037")
  for (i in 1:length(EunYoungKim2022_id)) {
    EunYoungKim2022_sample <- paste("EunYoungKim2022.", EunYoungKim2022_id[i], sep = "")
    EunYoungKim2022_seu <- Read10X(paste("/data/mengxu/data/PRJNA773987/", EunYoungKim2022_file[i], "/", EunYoungKim2022_file[i], "_f", "/outs/filtered_feature_bc_matrix", sep = ""))

    EunYoungKim2022_seu <- CreateSeuratObject(
      counts = EunYoungKim2022_seu,
      project = EunYoungKim2022_sample,
      min.features = 200,
      min.cells = 3
    )

    EunYoungKim2022_seu$platform <- "10X"
    ### QC
    EunYoungKim2022_seu <- PercentageFeatureSet(EunYoungKim2022_seu, pattern = "^MT-", col.name = "pMT")
    EunYoungKim2022_seu <- PercentageFeatureSet(EunYoungKim2022_seu, pattern = "^HBA|^HBB", col.name = "pHB")
    EunYoungKim2022_seu <- PercentageFeatureSet(EunYoungKim2022_seu, pattern = "^RPS|^RPL", col.name = "pRP")
    EunYoungKim2022_seu <- subset(EunYoungKim2022_seu,
      subset =
        nFeature_RNA > nFeature_lower &
          nFeature_RNA < nFeature_upper &
          nCount_RNA > nCount_lower &
          nCount_RNA < nCount_upper &
          pMT < pMT_upper &
          pMT > pMT_lower
    )
    # Combi.sqj <- subset(Combi.sbj, subset = percent.mt < 20 & nCount_RNA >100 & nCount_RNA < 150000 & nFeature_RNA > 200 & nFeature_RNA < 10000)
    ###
    EunYoungKim2022[[i]] <- EunYoungKim2022_seu
    EunYoungKim2022_samples[i] <- EunYoungKim2022_sample
  }
}

#----------------------------------------------------------------------------------------------#

samples <- c(
  PhilipBischoff2021_samples # ,
  # DiHe2021_samples,
  # 'DietherLambrechts2018.P06',
  # 'DietherLambrechts2018.P07',
  # 'DietherLambrechts2018.P08',
  # 'DietherLambrechts2018.P02',
  # 'DietherLambrechts2018.P03',
  # 'DietherLambrechts2018.P04',
  # 'DietherLambrechts2018.P05',
  # AndrewMLeader2021_samples,
  # EunYoungKim2022_samples
)

#----------------------------------------------------------------------------------------------#

seu_obj_list <- c(
  PhilipBischoff2021 # ,
  # DiHe2021,
  # DietherLambrechts2018.P06,
  # DietherLambrechts2018.P07,
  # DietherLambrechts2018.P08,
  # DietherLambrechts2018.P02,
  # DietherLambrechts2018.P03,
  # DietherLambrechts2018.P04,
  # DietherLambrechts2018.P05,
  # AndrewMLeader2021,
  # EunYoungKim2022
)

save(seu_obj_list, samples, file = "/data/mengxu/data/all/lung_stage-normal_list_raw.Rdata")

#----------------------------------------------------------------------------------------#
rm(list = ls())
gc()
