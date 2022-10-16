

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
  RapolasZilionis2019.P3.1 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635288_human_p3t1_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )
  RapolasZilionis2019.P3.2 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635289_human_p3t2_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )
  RapolasZilionis2019.P3.3 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635290_human_p3t3_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )

  RapolasZilionis2019.P3 <- cbind.data.frame(
    t(RapolasZilionis2019.P3.1),
    t(RapolasZilionis2019.P3.2),
    t(RapolasZilionis2019.P3.3)
  )

  RapolasZilionis2019.P3 <- CreateSeuratObject(
    counts = RapolasZilionis2019.P3,
    project = "RapolasZilionis2019.P3",
    min.features = 200,
    min.cells = 3
  )
  RapolasZilionis2019.P3$platform <- "IndropSeq"
  ### QC
  RapolasZilionis2019.P3 <- PercentageFeatureSet(RapolasZilionis2019.P3, pattern = "^MT-", col.name = "pMT")
  RapolasZilionis2019.P3 <- PercentageFeatureSet(RapolasZilionis2019.P3, pattern = "^HBA|^HBB", col.name = "pHB")
  RapolasZilionis2019.P3 <- PercentageFeatureSet(RapolasZilionis2019.P3, pattern = "^RPS|^RPL", col.name = "pRP")
  RapolasZilionis2019.P3 <- subset(RapolasZilionis2019.P3,
    subset =
    # nFeature_RNA > nFeature_lower &
    # nFeature_RNA < nFeature_upper &
      nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )

  ###
  RapolasZilionis2019.P5.1 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635296_human_p5t1_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )
  RapolasZilionis2019.P5.2 <- read.table("/data/mengxu/data/GSE127465/GSE127465/GSM3635297_human_p5t2_raw_counts.tsv.gz",
    row.names = 1,
    header = T
  )

  RapolasZilionis2019.P5 <- cbind.data.frame(
    t(RapolasZilionis2019.P5.1),
    t(RapolasZilionis2019.P5.2)
  )

  RapolasZilionis2019.P5 <- CreateSeuratObject(
    counts = RapolasZilionis2019.P5,
    project = "RapolasZilionis2019.P5",
    min.features = 200,
    min.cells = 3
  )
  RapolasZilionis2019.P5$platform <- "IndropSeq"
  ### QC
  RapolasZilionis2019.P5 <- PercentageFeatureSet(RapolasZilionis2019.P5, pattern = "^MT-", col.name = "pMT")
  RapolasZilionis2019.P5 <- PercentageFeatureSet(RapolasZilionis2019.P5, pattern = "^HBA|^HBB", col.name = "pHB")
  RapolasZilionis2019.P5 <- PercentageFeatureSet(RapolasZilionis2019.P5, pattern = "^RPS|^RPL", col.name = "pRP")
  RapolasZilionis2019.P5 <- subset(RapolasZilionis2019.P5,
    subset =
    # nFeature_RNA > nFeature_lower &
    # nFeature_RNA < nFeature_upper &
      nCount_RNA > nCount_lower &
        # nCount_RNA < nCount_upper &
        pMT < pMT_upper
  )
}

#----------------------------------------------------------------------------------------------#
# GSE131907_NayoungKim2020_10x Genomics
#
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
  patient_ids <- c("P1006", "P1028", "P1049", "P1058")
  file_ids <- c("_EBUS_06", "_EBUS_28", "_EBUS_49", "_BRONCHO_58")

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

#----------------------------------------------------------------------------------------------#
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

  ######
  if (Sys.info()[1] == "Windows") {
    # Metadata
    neo_osi_metadata <- read.csv("/\\192.168.0.109/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/neo-osi_metadata.csv")
    S01_metacells <- read.csv("/\\192.168.0.109/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/S01_metacells.csv")
    # Dataset
    # load("/data/mengxu/data/AshleyMaynard2020/Data_input/objects/S02.1_Main_Seurat_object_filtered_neo_osi.RData")
    # AshleyMaynard2020 <- read.csv('/data/mengxu/data/AshleyMaynard2020/Data_input/csv_files/S01_datafinal.csv')
    neo_osi_rawdata <- fread(
      file = "/\\192.168.0.109/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/neo-osi_rawdata.csv",
      sep = ",",
      header = T,
      check.names = F
    )
    neo_osi_rawdata <- as.data.frame(na.omit(neo_osi_rawdata))
    rownames(neo_osi_rawdata) <- neo_osi_rawdata$gene
    ###
    S01_datafinal <- fread(
      file = "/\\192.168.0.109/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/S01_datafinal.csv",
      sep = ",",
      header = T,
      check.names = F
    )
    S01_datafinal <- as.data.frame(na.omit(S01_datafinal))
    rownames(S01_datafinal) <- S01_datafinal$V1
  } else {
    # Metadata
    neo_osi_metadata <- read.csv("/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/neo-osi_metadata.csv")
    S01_metacells <- read.csv("/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/S01_metacells.csv")
    # Dataset
    neo_osi_rawdata <- fread(
      file = "/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/neo-osi_rawdata.csv",
      sep = ",",
      header = T,
      check.names = F
    )
    neo_osi_rawdata <- as.data.frame(na.omit(neo_osi_rawdata))
    rownames(neo_osi_rawdata) <- neo_osi_rawdata$gene
    ###
    S01_datafinal <- fread(
      file = "/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/S01_datafinal.csv",
      sep = ",",
      header = T,
      check.names = F
    )
    S01_datafinal <- as.data.frame(na.omit(S01_datafinal))
    rownames(S01_datafinal) <- S01_datafinal$V1
  }

  AshleyMaynard2020 <- list()
  AshleyMaynard2020_samples <- c()
  AshleyMaynard2020_id <- c(
    "LT_S01", "LT_S02", "LT_S03", "LT_S05", "LT_S07", "LT_S08", "LT_S11", "LT_S13", "LT_S14", "LT_S16",
    "LT_S21", "LT_S23", "LT_S28", "LT_S34", "LT_S41", "LT_S42", "LT_S47", "LT_S48", "LT_S50", "LT_S51",
    "LT_S52", "LT_S53", "LT_S54", "LT_S56", "LT_S57", "LT_S58", "LT_S63", "LT_S67", "LT_S71", "LT_S72",
    "LT_S75", "LT_S78", "LT_S79", "LT_S80", "LT_S81", "LT_S82"
  )

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
  "RapolasZilionis2019.P3",
  "RapolasZilionis2019.P5",
  NayoungKim2020_samples,
  AshleyMaynard2020_samples
)
#----------------------------------------------------------------------------------------------#

seu_obj_list <- c(
  RapolasZilionis2019.P3,
  RapolasZilionis2019.P5,
  NayoungKim2020,
  AshleyMaynard2020
)

save(seu_obj_list, samples, file = "/data/mengxu/data/all/lung_stage-4_list_raw.Rdata")
#----------------------------------------------------------------------------------------#
rm(list = ls())
gc()
