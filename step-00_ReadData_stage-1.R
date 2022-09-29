

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

#-------------------------------------------------------------------------------#
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
  PhilipBischoff2021_id <- c("23", "24")
  for (i in 1:length(PhilipBischoff2021_id)) {
    PhilipBischoff2021_sample <- paste("PhilipBischoff2021.P0", PhilipBischoff2021_id[i], sep = "")

    PhilipBischoff2021_seu <- Read10X(paste("/data/mengxu/data/PhilipBischoff2021/cellranger/p0", PhilipBischoff2021_id[i], "t/filtered_feature_bc_matrix", sep = ""))

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
  ###
}

#-------------------------------------------------------------------------------#
# CRA001963_DiHe2021_10x Genomics
if (T) {
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
  DiHe2021_id <- c("T3", "T4", "T1")
  DiHe2021_file <- c("CRR073024", "CRR073025", "CRR073022")
  for (i in 1:length(DiHe2021_id)) {
    DiHe2021_sample <- paste("DiHe2021.", DiHe2021_id[i], sep = "")
    DiHe2021_seu <- Read10X(paste0("/data/mengxu/data/CRA001963/", DiHe2021_file[i], "/", DiHe2021_file[i], "_f", "/outs/filtered_feature_bc_matrix"))

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

# CRA001477_DiHe2021_10x Genomics
if (T) {
  # Cells with low feature counts (<200) and high percent of mitochondrial genes (>10%) were removed.

  # nFeature_lower <- 201
  # nFeature_upper <- 10000
  nCount_lower <- 200
  # nCount_upper <- 6000
  pMT_lower <- 10
  # pMT_upper <- 10
  # pHB_lower <- 0
  # pHB_upper <- 5

  DiHe2021_CRA001477 <- list()
  DiHe2021_CRA001477_samples <- c()
  DiHe2021_CRA001477_id <- c("T6", "T7")
  DiHe2021_CRA001477_file <- c("CRR049227_28", "CRR049229_30")
  for (i in 1:length(DiHe2021_CRA001477_id)) {
    DiHe2021_CRA001477_sample <- paste("DiHe2021.", DiHe2021_CRA001477_id[i], sep = "")
    DiHe2021_CRA001477_seu <- Read10X(paste("/data/mengxu/data/CRA001477/", DiHe2021_CRA001477_file[i], "_f", "/outs/count/filtered_feature_bc_matrix", sep = ""))

    DiHe2021_CRA001477_seu <- CreateSeuratObject(
      counts = DiHe2021_CRA001477_seu,
      project = DiHe2021_CRA001477_sample,
      min.features = 200,
      min.cells = 3
    )

    DiHe2021_CRA001477_seu$platform <- "10X"
    ### QC
    DiHe2021_CRA001477_seu <- PercentageFeatureSet(DiHe2021_CRA001477_seu, pattern = "^MT-", col.name = "pMT")
    DiHe2021_CRA001477_seu <- PercentageFeatureSet(DiHe2021_CRA001477_seu, pattern = "^HBA|^HBB", col.name = "pHB")
    DiHe2021_CRA001477_seu <- PercentageFeatureSet(DiHe2021_CRA001477_seu, pattern = "^RPS|^RPL", col.name = "pRP")
    DiHe2021_CRA001477_seu <- subset(DiHe2021_CRA001477_seu,
      subset =
      # nFeature_RNA > nFeature_lower &
      # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
          # nCount_RNA < nCount_upper &
          # pMT < pMT_upper&
          pMT > pMT_lower
    )
    ###
    DiHe2021_CRA001477[[i]] <- DiHe2021_CRA001477_seu
    DiHe2021_CRA001477_samples[i] <- DiHe2021_CRA001477_sample
  }
}

#-------------------------------------------------------------------------------#
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

  DietherLambrechts2018.P07.1 <- Read10X("/data/mengxu/data/E-MTAB-6653/BT1375/BT1375_f/outs/filtered_feature_bc_matrix")

  DietherLambrechts2018.P07.1 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P07.1,
    project = "DietherLambrechts2018.P07.1",
    min.features = 200,
    min.cells = 3
  )

  DietherLambrechts2018.P07.1$platform <- "10X"
  ### QC
  DietherLambrechts2018.P07.1 <- PercentageFeatureSet(DietherLambrechts2018.P07.1, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P07.1 <- PercentageFeatureSet(DietherLambrechts2018.P07.1, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P07.1 <- PercentageFeatureSet(DietherLambrechts2018.P07.1, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P07.1 <- subset(DietherLambrechts2018.P07.1,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )

  ###
  DietherLambrechts2018.P07.2 <- Read10X("/data/mengxu/data/E-MTAB-6653/BT1376/BT1376_f/outs/filtered_feature_bc_matrix")

  DietherLambrechts2018.P07.2 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P07.2,
    project = "DietherLambrechts2018.P07.2",
    min.features = 200,
    min.cells = 3
  )

  DietherLambrechts2018.P07.2$platform <- "10X"
  ### QC
  DietherLambrechts2018.P07.2 <- PercentageFeatureSet(DietherLambrechts2018.P07.2, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P07.2 <- PercentageFeatureSet(DietherLambrechts2018.P07.2, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P07.2 <- PercentageFeatureSet(DietherLambrechts2018.P07.2, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P07.2 <- subset(DietherLambrechts2018.P07.2,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P07.3 <- Read10X("/data/mengxu/data/E-MTAB-6653/BT1377/BT1377_f/outs/filtered_feature_bc_matrix")

  DietherLambrechts2018.P07.3 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P07.3,
    project = "DietherLambrechts2018.P07.3",
    min.features = 200,
    min.cells = 3
  )

  DietherLambrechts2018.P07.3$platform <- "10X"
  ### QC
  DietherLambrechts2018.P07.3 <- PercentageFeatureSet(DietherLambrechts2018.P07.3, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P07.3 <- PercentageFeatureSet(DietherLambrechts2018.P07.3, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P07.3 <- PercentageFeatureSet(DietherLambrechts2018.P07.3, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P07.3 <- subset(DietherLambrechts2018.P07.3,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P07 <- merge(DietherLambrechts2018.P07.1,
    y = c(
      DietherLambrechts2018.P07.2,
      DietherLambrechts2018.P07.3
    ),
    add.cell.ids = c(
      "DietherLambrechts2018.P07",
      "DietherLambrechts2018.P07",
      "DietherLambrechts2018.P07"
    )
  )
}

# E-MTAB-6149_DietherLambrechts2018_10X Genomics
if (F) {
  DietherLambrechts2018.P05.1 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT1298/BT1298_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P05.1 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P05.1,
    project = "DietherLambrechts2018.P05.1",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P05.1$platform <- "10X"
  ### QC
  DietherLambrechts2018.P05.1 <- PercentageFeatureSet(DietherLambrechts2018.P05.1, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P05.1 <- PercentageFeatureSet(DietherLambrechts2018.P05.1, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P05.1 <- PercentageFeatureSet(DietherLambrechts2018.P05.1, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P05.1 <- subset(DietherLambrechts2018.P05.1,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P05.2 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT1299/BT1299_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P05.2 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P05.2,
    project = "DietherLambrechts2018.P05.2",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P05.2$platform <- "10X"
  ### QC
  DietherLambrechts2018.P05.2 <- PercentageFeatureSet(DietherLambrechts2018.P05.2, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P05.2 <- PercentageFeatureSet(DietherLambrechts2018.P05.2, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P05.2 <- PercentageFeatureSet(DietherLambrechts2018.P05.2, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P05.2 <- subset(DietherLambrechts2018.P05.2,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P05.3 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT1300/BT1300_f/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P05.3 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P05.3,
    project = "DietherLambrechts2018.P05.3",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P05.3$platform <- "10X"
  ### QC
  DietherLambrechts2018.P05.3 <- PercentageFeatureSet(DietherLambrechts2018.P05.3, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P05.3 <- PercentageFeatureSet(DietherLambrechts2018.P05.3, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P05.3 <- PercentageFeatureSet(DietherLambrechts2018.P05.3, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P05.3 <- subset(DietherLambrechts2018.P05.3,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P05 <- merge(DietherLambrechts2018.P05.1,
    y = c(
      DietherLambrechts2018.P05.2,
      DietherLambrechts2018.P05.3
    ),
    add.cell.ids = c(
      "DietherLambrechts2018.P05",
      "DietherLambrechts2018.P05",
      "DietherLambrechts2018.P05"
    )
  )

  ###
  DietherLambrechts2018.P02.1 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT2A/BT2A/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P02.1 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P02.1,
    project = "DietherLambrechts2018.P02.1",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P02.1$platform <- "10X"
  ### QC
  DietherLambrechts2018.P02.1 <- PercentageFeatureSet(DietherLambrechts2018.P02.1, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P02.1 <- PercentageFeatureSet(DietherLambrechts2018.P02.1, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P02.1 <- PercentageFeatureSet(DietherLambrechts2018.P02.1, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P02.1 <- subset(DietherLambrechts2018.P02.1,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P02.2 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT2B/BT2B/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P02.2 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P02.2,
    project = "DietherLambrechts2018.P02.2",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P02.2$platform <- "10X"
  ### QC
  DietherLambrechts2018.P02.2 <- PercentageFeatureSet(DietherLambrechts2018.P02.2, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P02.2 <- PercentageFeatureSet(DietherLambrechts2018.P02.2, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P02.2 <- PercentageFeatureSet(DietherLambrechts2018.P02.2, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P02.2 <- subset(DietherLambrechts2018.P02.2,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P02.3 <- Read10X("/data/mengxu/data/E-MTAB-6149/BT1249/BT1249/outs/filtered_feature_bc_matrix")
  DietherLambrechts2018.P02.3 <- CreateSeuratObject(
    counts = DietherLambrechts2018.P02.3,
    project = "DietherLambrechts2018.P02.3",
    min.features = 200,
    min.cells = 3
  )
  DietherLambrechts2018.P02.3$platform <- "10X"
  ### QC
  DietherLambrechts2018.P02.3 <- PercentageFeatureSet(DietherLambrechts2018.P02.3, pattern = "^MT-", col.name = "pMT")
  DietherLambrechts2018.P02.3 <- PercentageFeatureSet(DietherLambrechts2018.P02.3, pattern = "^HBA|^HBB", col.name = "pHB")
  DietherLambrechts2018.P02.3 <- PercentageFeatureSet(DietherLambrechts2018.P02.3, pattern = "^RPS|^RPL", col.name = "pRP")
  DietherLambrechts2018.P02.3 <- subset(DietherLambrechts2018.P02.3,
    subset =
      nFeature_RNA > nFeature_lower &
        # nFeature_RNA < nFeature_upper &
        nCount_RNA > nCount_lower &
        nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
  ###
  DietherLambrechts2018.P02 <- merge(DietherLambrechts2018.P02.1,
    y = c(
      DietherLambrechts2018.P02.2,
      DietherLambrechts2018.P02.3
    ),
    add.cell.ids = c(
      "DietherLambrechts2018.P02",
      "DietherLambrechts2018.P02",
      "DietherLambrechts2018.P02"
    )
  )
}

#-------------------------------------------------------------------------------#
# GSE127465_RapolasZilionis2019_IndropSeq
if (F) {

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
  RapolasZilionis2019.P4.1 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635292_human_p4t1_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )
  RapolasZilionis2019.P4.2 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635293_human_p4t2_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )
  RapolasZilionis2019.P4.3 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635294_human_p4t3_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )


  RapolasZilionis2019.P4 <- cbind.data.frame(
    t(RapolasZilionis2019.P4.1),
    t(RapolasZilionis2019.P4.2),
    t(RapolasZilionis2019.P4.3)
  )

  RapolasZilionis2019.P4 <- CreateSeuratObject(
    counts = RapolasZilionis2019.P4,
    project = "RapolasZilionis2019.P4",
    min.features = 200,
    min.cells = 3
  )
  RapolasZilionis2019.P4$platform <- "IndropSeq"
  ### QC
  RapolasZilionis2019.P4 <- PercentageFeatureSet(RapolasZilionis2019.P4, pattern = "^MT-", col.name = "pMT")
  RapolasZilionis2019.P4 <- PercentageFeatureSet(RapolasZilionis2019.P4, pattern = "^HBA|^HBB", col.name = "pHB")
  RapolasZilionis2019.P4 <- PercentageFeatureSet(RapolasZilionis2019.P4, pattern = "^RPS|^RPL", col.name = "pRP")
  RapolasZilionis2019.P4 <- subset(RapolasZilionis2019.P4,
    subset =
    # nFeature_RNA > nFeature_lower &
    # nFeature_RNA < nFeature_upper &
      nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )

  ###
  RapolasZilionis2019.P6.1 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635298_human_p6t1_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )
  RapolasZilionis2019.P6.2 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635299_human_p6t2_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )

  RapolasZilionis2019.P6 <- cbind.data.frame(
    t(RapolasZilionis2019.P6.1),
    t(RapolasZilionis2019.P6.2)
  )

  RapolasZilionis2019.P6 <- CreateSeuratObject(
    counts = RapolasZilionis2019.P6,
    project = "RapolasZilionis2019.P6",
    min.features = 200,
    min.cells = 3
  )
  RapolasZilionis2019.P6$platform <- "IndropSeq"
  ### QC
  RapolasZilionis2019.P6 <- PercentageFeatureSet(RapolasZilionis2019.P6, pattern = "^MT-", col.name = "pMT")
  RapolasZilionis2019.P6 <- PercentageFeatureSet(RapolasZilionis2019.P6, pattern = "^HBA|^HBB", col.name = "pHB")
  RapolasZilionis2019.P6 <- PercentageFeatureSet(RapolasZilionis2019.P6, pattern = "^RPS|^RPL", col.name = "pRP")
  RapolasZilionis2019.P6 <- subset(RapolasZilionis2019.P6,
    subset =
    # nFeature_RNA > nFeature_lower &
    # nFeature_RNA < nFeature_upper &
      nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
}

#-------------------------------------------------------------------------------#
# GSE176021_JustinaXCaushi2021_10x Genomics
if (F) {

  # # The quality of cells was then assessed based on (1) the number  of genes detected per cell and
  # # (2) the proportion of mitochondrial gene/ ribosomal gene counts.
  # # Low-quality cells were filtered if the number  of detected genes was below 250 or
  # # above 3× the median absolute  deviation away from the median gene number of all cells.
  # # Cells were  filtered out if the proportion of mitochondrial gene counts was higher
  # # than 10% or the proportion of ribosomal genes was less than 10%.
  #
  # ###clustering.R
  # # gene_count <- as.matrix(JustinaXCaushi2021.MD01.010@assays$RNA@counts)
  # # gene_count_norm <- sweep(gene_count,2,colSums(gene_count),FUN="/")*1e4 #把减去均值后的矩阵在列的方向上除以极差向量
  # # gene_hypervar <- hypervar(gene_count_norm,showplot=FALSE)
  # # gene_hypervar_sort <- gene_hypervar$data %>% arrange(.,desc(hypervar_log2))
  # # VariableFeatures(JustinaXCaushi2021.MD01.010) <- gene_hypervar_sort$feature[1:3000]
  # # VariableFeatures(JustinaXCaushi2021.MD01.010)<-setdiff(VariableFeatures(JustinaXCaushi2021.MD01.010),cluster.exclude)
  #
  # nFeature_lower <- 500
  # nFeature_upper <- 4000
  # nCount_lower <- 300
  # nCount_upper <- 10000
  # pMT_lower <- 0
  # pMT_upper <- 15
  # #pHB_lower <- 0
  # #pHB_upper <- 5
  #
  # JustinaXCaushi2021.MD01.010 <- Read10X("/data/mengxu/data/GSE176021/MD01-010_tumor_1")
  # JustinaXCaushi2021.MD01.010 <- CreateSeuratObject(counts = JustinaXCaushi2021.MD01.010 ,
  #                                                   project = 'JustinaXCaushi2021.MD01.010',
  #                                                   min.features = 200,
  #                                                   min.cells = 3)
  # ###QC
  # JustinaXCaushi2021.MD01.010 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.010, pattern = "^MT-", col.name = "pMT")
  # JustinaXCaushi2021.MD01.010 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.010, pattern = "^HBA|^HBB", col.name = "pHB")
  # JustinaXCaushi2021.MD01.010 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.010, pattern = "^RPS|^RPL", col.name = "pRP")
  # JustinaXCaushi2021.MD01.010 <- subset(JustinaXCaushi2021.MD01.010,
  #                                       subset =
  #                                         nFeature_RNA > nFeature_lower &
  #                                         nFeature_RNA < nFeature_upper &
  #                                         nCount_RNA > nCount_lower &
  #                                         nCount_RNA < nCount_upper &
  #                                         pMT < pMT_upper)
  # ###
  # JustinaXCaushi2021.MD01.005_list <- list()
  #
  # for (i in 2:9 ) {
  #   JustinaXCaushi2021.MD01.005_list[[i]] <- Read10X(paste("/data/mengxu/data/GSE176021/MD01-005_tumor_",i,sep = ''))
  #
  #   JustinaXCaushi2021.MD01.005_list[[i]] <- CreateSeuratObject(counts = JustinaXCaushi2021.MD01.005_list[[i]],
  #                                                       project = 'JustinaXCaushi2021.MD01.005',
  #                                                       min.features = 200,
  #                                                       min.cells = 3)
  #
  #   JustinaXCaushi2021.MD01.005_list[[i]] <- RenameCells(JustinaXCaushi2021.MD01.005_list[[i]],
  #                                                        add.cell.id = 'JustinaXCaushi2021.MD01.005')
  #
  # }
  #
  # JustinaXCaushi2021.MD01.005 <- merge(JustinaXCaushi2021.MD01.005_list[[2]],
  #                                      y=c(JustinaXCaushi2021.MD01.005_list[[3]],
  #                                          JustinaXCaushi2021.MD01.005_list[[4]],
  #                                          JustinaXCaushi2021.MD01.005_list[[5]],
  #                                          JustinaXCaushi2021.MD01.005_list[[6]],
  #                                          JustinaXCaushi2021.MD01.005_list[[7]],
  #                                          JustinaXCaushi2021.MD01.005_list[[8]],
  #                                          JustinaXCaushi2021.MD01.005_list[[9]]))
  # ###QC
  # JustinaXCaushi2021.MD01.005 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.005, pattern = "^MT-", col.name = "pMT")
  # JustinaXCaushi2021.MD01.005 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.005, pattern = "^HBA|^HBB", col.name = "pHB")
  # JustinaXCaushi2021.MD01.005 <- PercentageFeatureSet(JustinaXCaushi2021.MD01.005, pattern = "^RPS|^RPL", col.name = "pRP")
  # JustinaXCaushi2021.MD01.005 <- subset(JustinaXCaushi2021.MD01.005,
  #                                       subset =
  #                                         nFeature_RNA > nFeature_lower &
  #                                         nFeature_RNA < nFeature_upper &
  #                                         nCount_RNA > nCount_lower &
  #                                         nCount_RNA < nCount_upper &
  #                                         pMT < pMT_upper)
  # ###
  # JustinaXCaushi2021.NY016.021 <- Read10X("/data/mengxu/data/GSE176021/NY016-021_tumor_1")
  #
  # JustinaXCaushi2021.NY016.021 <- CreateSeuratObject(counts = JustinaXCaushi2021.NY016.021 ,
  #                                                   project = 'JustinaXCaushi2021.NY016.021',
  #                                                   min.features = 200,
  #                                                   min.cells = 3)
  # ###QC
  # JustinaXCaushi2021.NY016.021 <- PercentageFeatureSet(JustinaXCaushi2021.NY016.021, pattern = "^MT-", col.name = "pMT")
  # JustinaXCaushi2021.NY016.021 <- PercentageFeatureSet(JustinaXCaushi2021.NY016.021, pattern = "^HBA|^HBB", col.name = "pHB")
  # JustinaXCaushi2021.NY016.021 <- PercentageFeatureSet(JustinaXCaushi2021.NY016.021, pattern = "^RPS|^RPL", col.name = "pRP")
  # JustinaXCaushi2021.NY016.021 <- subset(JustinaXCaushi2021.NY016.021,
  #                                        subset =
  #                                          nFeature_RNA > nFeature_lower &
  #                                          nFeature_RNA < nFeature_upper &
  #                                          nCount_RNA > nCount_lower &
  #                                          nCount_RNA < nCount_upper &
  #                                          pMT < pMT_upper)
}

#-------------------------------------------------------------------------------#
# GSE99254_XinyiGuo2018_Smart-Seq2
if (F) {
  # Low-quality cells were discarded if the cell library size or the number of  expressed genes (counts larger than 0) was smaller than pre-defined thresholds,
  # which were the medians of all cells minus 3×median absolute deviation.
  # Cells  were also removed if their proportions of mitochondrial gene expression were  larger than 10%.

  # nFeature_lower <- 201
  # nFeature_upper <- 4000
  # nCount_lower <- 500
  # nCount_upper <- 1000000
  # pMT_lower <- 0
  # pMT_upper <- 10
  # #pHB_lower <- 0
  # #pHB_upper <- 5
  #
  # GSE99254 <- read.table("/data/mengxu/data/GSE99254/GSE99254_NSCLC.TCell.S12346.count.txt.gz",
  #                        row.names = 1,
  #                        header = T)
  # GSE99254 <- GSE99254[!duplicated(GSE99254$symbol), ] #remove 88 features
  # GSE99254 <- na.omit(GSE99254) #remove 1 feature
  # row.names(GSE99254) <- GSE99254$symbol
  # ###
  # XinyiGuo2018.P0617 <- GSE99254[c(grep("0617",colnames(GSE99254)))]
  # XinyiGuo2018.P0617 <- CreateSeuratObject(counts = XinyiGuo2018.P0617,
  #                                          project = 'XinyiGuo2018.P0617',
  #                                          min.features = 200,
  #                                          min.cells = 3)
  # XinyiGuo2018.P0617$platform <- 'SmartSeq2'
  # ###QC
  # XinyiGuo2018.P0617 <- PercentageFeatureSet(XinyiGuo2018.P0617, pattern = "^MT-", col.name = "pMT")
  # XinyiGuo2018.P0617 <- PercentageFeatureSet(XinyiGuo2018.P0617, pattern = "^HBA|^HBB", col.name = "pHB")
  # XinyiGuo2018.P0617 <- PercentageFeatureSet(XinyiGuo2018.P0617, pattern = "^RPS|^RPL", col.name = "pRP")
  # XinyiGuo2018.P0617 <- subset(XinyiGuo2018.P0617,
  #                              subset =
  #                                nFeature_RNA > nFeature_lower &
  #                                nFeature_RNA < nFeature_upper &
  #                                nCount_RNA > nCount_lower &
  #                                nCount_RNA < nCount_upper &
  #                                pMT < pMT_upper)
  # ###
  # XinyiGuo2018.P0619 <- GSE99254[c(grep("0619",colnames(GSE99254)))]
  # XinyiGuo2018.P0619 <- CreateSeuratObject(counts = XinyiGuo2018.P0619,
  #                                          project = 'XinyiGuo2018.P0619',
  #                                          min.features = 200,
  #                                          min.cells = 3)
  # XinyiGuo2018.P0619$platform <- 'SmartSeq2'
  # ###QC
  # XinyiGuo2018.P0619 <- PercentageFeatureSet(XinyiGuo2018.P0619, pattern = "^MT-", col.name = "pMT")
  # XinyiGuo2018.P0619 <- PercentageFeatureSet(XinyiGuo2018.P0619, pattern = "^HBA|^HBB", col.name = "pHB")
  # XinyiGuo2018.P0619 <- PercentageFeatureSet(XinyiGuo2018.P0619, pattern = "^RPS|^RPL", col.name = "pRP")
  # XinyiGuo2018.P0619 <- subset(XinyiGuo2018.P0619,
  #                              subset =
  #                                nFeature_RNA > nFeature_lower &
  #                                nFeature_RNA < nFeature_upper &
  #                                nCount_RNA > nCount_lower &
  #                                nCount_RNA < nCount_upper &
  #                                pMT < pMT_upper)
  # ###
  # XinyiGuo2018.P0729 <- GSE99254[c(grep("0729",colnames(GSE99254)))]
  # XinyiGuo2018.P0729 <- CreateSeuratObject(counts = XinyiGuo2018.P0729,
  #                                          project = 'XinyiGuo2018.P0729',
  #                                          min.features = 200,
  #                                          min.cells = 3)
  # XinyiGuo2018.P0729$platform <- 'SmartSeq2'
  # ###QC
  # XinyiGuo2018.P0729 <- PercentageFeatureSet(XinyiGuo2018.P0729, pattern = "^MT-", col.name = "pMT")
  # XinyiGuo2018.P0729 <- PercentageFeatureSet(XinyiGuo2018.P0729, pattern = "^HBA|^HBB", col.name = "pHB")
  # XinyiGuo2018.P0729 <- PercentageFeatureSet(XinyiGuo2018.P0729, pattern = "^RPS|^RPL", col.name = "pRP")
  # XinyiGuo2018.P0729 <- subset(XinyiGuo2018.P0729,
  #                              subset =
  #                                nFeature_RNA > nFeature_lower &
  #                                nFeature_RNA < nFeature_upper &
  #                                nCount_RNA > nCount_lower &
  #                                nCount_RNA < nCount_upper &
  #                                pMT < pMT_upper)
  # ###
  # XinyiGuo2018.P0913 <- GSE99254[c(grep("0913",colnames(GSE99254)))]
  # XinyiGuo2018.P0913 <- CreateSeuratObject(counts = XinyiGuo2018.P0913,
  #                                          project = 'XinyiGuo2018.P0913',
  #                                          min.features = 200,
  #                                          min.cells = 3)
  # XinyiGuo2018.P0913$platform <- 'SmartSeq2'
  # ###QC
  # XinyiGuo2018.P0913 <- PercentageFeatureSet(XinyiGuo2018.P0913, pattern = "^MT-", col.name = "pMT")
  # XinyiGuo2018.P0913 <- PercentageFeatureSet(XinyiGuo2018.P0913, pattern = "^HBA|^HBB", col.name = "pHB")
  # XinyiGuo2018.P0913 <- PercentageFeatureSet(XinyiGuo2018.P0913, pattern = "^RPS|^RPL", col.name = "pRP")
  # XinyiGuo2018.P0913 <- subset(XinyiGuo2018.P0913,
  #                              subset =
  #                                nFeature_RNA > nFeature_lower &
  #                                nFeature_RNA < nFeature_upper &
  #                                nCount_RNA > nCount_lower &
  #                                nCount_RNA < nCount_upper &
  #                                pMT < pMT_upper)
}

#-------------------------------------------------------------------------------#
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

  ###
  GSE131907 <- fread(
    file = "/data/mengxu/data/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
    sep = "\t",
    header = T,
    check.names = F
  )
  GSE131907 <- as.data.frame(GSE131907) # remove 0 feature GSE131907 <- as.data.frame(na.omit(GSE131907))
  row.names(GSE131907) <- GSE131907$Index
  ###
  patient_ids <- c("P0006", "P0018", "P0019", "P0020", "P0025", "P0030", "P0034", "P0008")
  file_ids <- c("_LUNG_T06", "_LUNG_T18", "_LUNG_T19", "_LUNG_T20", "_LUNG_T25", "_LUNG_T30", "_LUNG_T34", "_LUNG_T08")

  NayoungKim2020 <- list()
  NayoungKim2020_samples <- c()
  for (i in 1:length(file_ids)) {
    NayoungKim2020_sample <- paste0("NayoungKim2020.", patient_ids[i])

    NayoungKim2020_seu <- GSE131907[c(grep(file_ids[i], colnames(GSE131907)))]

    fwrite(NayoungKim2020_seu,
      paste0("/data/mengxu/data/GSE131907/data_processed/NayoungKim2020_seu", file_ids[i], ".txt"),
      sep = "\t",
      row.names = T
    ) # fwrite保存行名会丢失？

    NayoungKim2020_seu <- read.table(paste0("/data/mengxu/data/GSE131907/data_processed/NayoungKim2020_seu", file_ids[i], ".txt"),
      row.names = 1,
      header = T
    )

    NayoungKim2020_seu <- CreateSeuratObject(
      counts = NayoungKim2020_seu,
      project = NayoungKim2020_sample,
      min.features = 200,
      min.cells = 3
    )
    NayoungKim2020_seu$platform <- "10X"
    ### QC
    NayoungKim2020_seu <- PercentageFeatureSet(NayoungKim2020_seu, pattern = "^MT-", col.name = "pMT")
    NayoungKim2020_seu <- PercentageFeatureSet(NayoungKim2020_seu, pattern = "^HBA|^HBB", col.name = "pHB")
    NayoungKim2020_seu <- PercentageFeatureSet(NayoungKim2020_seu, pattern = "^RPS|^RPL", col.name = "pRP")
    NayoungKim2020_seu <- subset(NayoungKim2020_seu,
      subset =
        nFeature_RNA > nFeature_lower &
          nFeature_RNA < nFeature_upper &
          nCount_RNA > nCount_lower &
          nCount_RNA < nCount_upper &
          pMT < pMT_upper
    )

    NayoungKim2020[[i]] <- NayoungKim2020_seu
    NayoungKim2020_samples[i] <- NayoungKim2020_sample
  }
}

#-------------------------------------------------------------------------------#
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
    "P458", "P532", "P558", "P571", "P626", "P626", "P695", "P695", "P695", "P695",
    "P706", "P800", "P569", "P578", "P371", "P378", "P393", "P406", "P410", "P460",
    "P464", "P564", "P572", "P593", "P596", "P596", "P630", "P338", "P338", "P725"
  )
  file_ids <- c(
    "23", "33", "36", "43", "88", "89", "91", "92", "93", "95",
    "342", "668", "40", "47", "4", "8", "9", "15", "19", "25",
    "27", "38", "45", "53", "54", "55", "56", "230", "231", "309"
  )

  AndrewMLeader2021 <- list()
  AndrewMLeader2021_samples <- c()
  for (i in 1:length(file_ids)) {
    patient <- patient_ids[i]
    file <- file_ids[i]
    AndrewMLeader2021_sample <- paste("AndrewMLeader2021", patient_ids[i], file_ids[i], sep = ".")
    message("[", Sys.time(), "]", " ------ Now: ", AndrewMLeader2021_sample)

    matrix <- Matrix::readMM(paste("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_", file_ids[i], "/matrix.mtx", sep = ""))
    # feature <- read.table(paste("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_",file_ids[i],"/features.tsv",sep = '')) #Bug
    feature <- read.delim(paste("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_", file_ids[i], "/features.tsv", sep = ""), stringsAsFactors = F, header = F)
    barcode <- read.table(paste("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_", file_ids[i], "/barcodes.tsv", sep = ""))
    # barcode <- read.delim(paste("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_",file_ids[i],"/barcodes.tsv",sep = ''),stringsAsFactors=F,header=F)[,1]
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

  # AndrewMLeader2021[[i]]@project.name

  # matrix <- Matrix::readMM("/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_29/matrix.mtx")
  # barcode <- read.table('/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_29/barcodes.tsv')
  # feature <- read.table('/data/mengxu/data/GSE154826/GSE154826_amp_batch_ID_29/features.tsv')
  # rownames(matrix) <- feature[,2]
  # colnames(matrix) <- barcode[,1]
  # AndrewMLeader2021.P514.29 <- CreateSeuratObject(counts = matrix,
  #                                                 project = 'AndrewMLeader2021.P514.29',
  #                                                 min.features = 200,
  #                                                 min.cells = 3)
  # AndrewMLeader2021.P514.29$platform <- '10X'
  # ###QC
  # AndrewMLeader2021.P514.29 <- PercentageFeatureSet(AndrewMLeader2021.P514.29, pattern = "^MT-", col.name = "pMT")
  # AndrewMLeader2021.P514.29 <- PercentageFeatureSet(AndrewMLeader2021.P514.29, pattern = "^HBA|^HBB", col.name = "pHB")
  # AndrewMLeader2021.P514.29 <- PercentageFeatureSet(AndrewMLeader2021.P514.29, pattern = "^RPS|^RPL", col.name = "pRP")
  # AndrewMLeader2021.P514.29 <- subset(AndrewMLeader2021.P514.29,
  #                                     subset =
  #                                       nFeature_RNA > nFeature_lower &
  #                                       #nFeature_RNA < nFeature_upper &
  #                                       nCount_RNA > nCount_lower &
  #                                       #nCount_RNA < nCount_upper &
  #                                       pMT < pMT_upper)
}

#-------------------------------------------------------------------------------#
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
  QianqianSong2019.P1 <- read.table("/data/mengxu/data/GSE117570/GSM3304007_P1_Tumor_processed_data.txt.gz",
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

  # QianqianSong2019.P3 <- read.table("/data/mengxu/data/GSE117570/GSM3304011_P3_Tumor_processed_data.txt.gz",
  #                                   row.names = 1,
  #                                   header = T)
  #
  #
  #
  # QianqianSong2019.P3 <- CreateSeuratObject(counts = QianqianSong2019.P3,
  #                                           project = 'QianqianSong2019.P3',
  #                                           min.features = 200,
  #                                           min.cells = 3)
  # QianqianSong2019.P3$platform <- '10X'
  # ###QC
  # QianqianSong2019.P3 <- PercentageFeatureSet(QianqianSong2019.P3, pattern = "^MT-", col.name = "pMT")
  # QianqianSong2019.P3 <- PercentageFeatureSet(QianqianSong2019.P3, pattern = "^HBA|^HBB", col.name = "pHB")
  # QianqianSong2019.P3 <- PercentageFeatureSet(QianqianSong2019.P3, pattern = "^RPS|^RPL", col.name = "pRP")
  # QianqianSong2019.P3 <- subset(QianqianSong2019.P3,
  #                               subset =
  #                                 nFeature_RNA > nFeature_lower &
  #                                 #nFeature_RNA < nFeature_upper &
  #                                 #nCount_RNA > nCount_lower &
  #                                 #nCount_RNA < nCount_upper &
  #                                 #pMT < pMT_upper&
  #                                 pMT > pMT_lower)
}

#-------------------------------------------------------------------------------#
# GSE123902_AshleyMLaughney2020_10x Genomics
if (F) {

  # Cells with >20% of transcripts derived from mitochondria were considered  apoptotic and also excluded (Extended Data Fig. 2c).
  # Afer clustering (described  below), two minority cell clusters were additionally excluded from patient data  with characteristically
  # low library size and signatures of partial cell lysis that  escaped the automated flters described above.
  # Genes detected in fewer than ten  cells or genes with low expression levels,
  # identifed as those with count values <5  s.d. from the second mode of the log–log distribution of total transcript counts and  genes,
  # were also excluded. Tis yielded a total of 40,505 patient-derived cells with  a median library size of 4,038 transcripts per cell (Extended Data Fig. 2d,e),
  # for  downstream analysis (Fig. 1 and Extended Data Fig. 2).

  nFeature_lower <- 200
  # nFeature_upper <- 10000
  # nCount_lower <- 200
  # nCount_upper <- 6000
  pMT_lower <- 20
  # pMT_upper <- 10
  # pHB_lower <- 0
  # pHB_upper <- 5

  ###
  AshleyMLaughney2020.LX679 <- read.table("/data/mengxu/data/GSE123902/GSM3516669_MSK_LX679_PRIMARY_TUMOUR_dense.csv.gz",
    row.names = 1,
    header = T,
    sep = ","
  )

  head(AshleyMLaughney2020.LX679[1:4, 1:4])
  AshleyMLaughney2020.LX679 <- t.data.frame(AshleyMLaughney2020.LX679)

  AshleyMLaughney2020.LX679 <- CreateSeuratObject(
    counts = AshleyMLaughney2020.LX679,
    project = "AshleyMLaughney2020.LX679",
    min.features = 200,
    min.cells = 3
  )
  AshleyMLaughney2020.LX679$platform <- "10X"
  ### QC
  AshleyMLaughney2020.LX679 <- PercentageFeatureSet(AshleyMLaughney2020.LX679, pattern = "^MT-", col.name = "pMT")
  AshleyMLaughney2020.LX679 <- PercentageFeatureSet(AshleyMLaughney2020.LX679, pattern = "^HBA|^HBB", col.name = "pHB")
  AshleyMLaughney2020.LX679 <- PercentageFeatureSet(AshleyMLaughney2020.LX679, pattern = "^RPS|^RPL", col.name = "pRP")
  AshleyMLaughney2020.LX679 <- subset(AshleyMLaughney2020.LX679,
    subset =
    # nFeature_RNA > nFeature_lower &
    # nFeature_RNA < nFeature_upper &
      nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
}

#-------------------------------------------------------------------------------#
# GSE176021_JustinaXCaushi2021.MD01-019
if (F) {
  JustinaXCaushi2021.MD01 - 019
  JustinaXCaushi2021.NY016 - 022
  JustinaXCaushi2021.MD043 - 003
}

#-------------------------------------------------------------------------------#
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

#-------------------------------------------------------------------------------#

samples <- c(
  PhilipBischoff2021_samples,
  DiHe2021_samples,
  DiHe2021_CRA001477_samples,
  # 'DietherLambrechts2018.P07',
  # 'DietherLambrechts2018.P05',
  # 'DietherLambrechts2018.P02',
  # 'RapolasZilionis2019.P4',
  # 'RapolasZilionis2019.P6',
  # AndrewMLeader2021_samples,
  NayoungKim2020_samples
)
#-------------------------------------------------------------------------------#

seu_obj_list <- c(
  PhilipBischoff2021,
  DiHe2021,
  DiHe2021_CRA001477,
  # DietherLambrechts2018.P07,
  # DietherLambrechts2018.P05,
  # DietherLambrechts2018.P02,
  # RapolasZilionis2019.P4,
  # RapolasZilionis2019.P6,
  # AndrewMLeader2021,
  NayoungKim2020
)

save(seu_obj_list, samples, file = "/data/mengxu/data/all/lung_stage-1_list_raw.Rdata")
#-------------------------------------------------------------------------------#
rm(list = ls())
gc()
