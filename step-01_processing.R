
rm(list = ls())
gc()
### load libraries
if (T) {
  library(cowplot)
  library(data.table) # fread
  library(dplyr)
  library(ggExtra)
  library(ggplot2)
  library(gridExtra)
  library(grid, lib.loc = "/opt/R/4.1.0/lib/R/library")
  library(GSVA)
  library(harmony)
  library(knitr)
  library(magrittr)
  library(markdown)
  library(Matrix)
  library(NNLM)
  library(patchwork)
  library(pheatmap)
  library(progeny)
  library(reshape2)
  library(reticulate)
  library(readr)
  library(readxl)
  library(stringr)
  library(scds)
  library(sctransform)
  library(Seurat)
  library(SingleCellExperiment)
  library(SoupX)
  library(tidyverse)
  library(tidyr)
  library(viridis)

  source("cnvFunction.R")
  source("sc_all_function.R")
  source("scCombination.R")
  source("plot_function.R")
}

#------------------------------------------------------------------------------#
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
  #' JustinaXCaushi2021.MD01.010',
  #' JustinaXCaushi2021.MD01.005', #
  #' JustinaXCaushi2021.NY016.021',
  "XinyiGuo2018.P0617",
  "XinyiGuo2018.P0619",
  "XinyiGuo2018.P0729",
  #' XinyiGuo2018.P0913', #Error in dimnames(x) <- dn : 'dimnames'的长度[2]必需与陈列范围相等
  "NayoungKim2020.P0028",
  "NayoungKim2020.P0031",
  "AndrewMLeader2021.P377.6",
  "AndrewMLeader2021.P403.11",
  "AndrewMLeader2021.P514.29",
  "AshleyMaynard2020.AZ.01",
  "AshleyMaynard2020.AZ.03",
  "AshleyMaynard2020.AZ.04",
  "AshleyMaynard2020.AZ.05",
  #' AshleyMaynard2020.LT.S49', #Error in dimnames(x) <- dn : 'dimnames'的长度[2]必需与陈列范围相等
  "AshleyMaynard2020.LT.S69",
  "AshleyMaynard2020.LT.S74"
)

###
for (i in 1:length(samples)) {
  dataPath <- "/data/mengxu/data/all/samples" # The path of cell ranger processed data
  savePath <- paste("/data/mengxu/results/Stage-3/", samples[i], sep = "") # A path to save the results
  sampleName <- "Stage-3" # The sample name
  statPath <- paste("/data/mengxu/results/Stage-3/", samples[i], sep = "") # A path to save the results
  sample <- samples[i]
  authorName <- "mengxu"

  # Run scStatistics
  stat.results <- runScStatistics(
    dataPath = dataPath,
    sample = sample,
    savePath = savePath,
    sampleName = sampleName,
    authorName = authorName
  )

  # Run scAnnotation
  anno.resultst <- runScAnnotation(
    dataPath = dataPath,
    sample = sample,
    statPath = statPath,
    savePath = savePath,
    authorName = authorName,
    sampleName = sampleName,
    geneSet.method = "GSVA" # "average" or "GSVA"
  )
}

# The paths of all sample's "runScAnnotation" results
single.savePaths <- c()
for (i in 1:length(samples)) {
  savePath <- paste("/data/mengxu/results/Stage-3/", samples[i], sep = "")
  single.savePaths <- c(single.savePaths, savePath)
}

sampleNames <- samples # The labels for all samples
savePath <- "/data/mengxu/results/test/Combination" # A path to save the results
combName <- "Stage-3" # A label of the combined samples
authorName <- "mengxu" # The author name to mark the report
comb.method <- "Harmony" # Integration methods ("NormalMNN", "SeuratMNN", "Harmony", "Raw", "Regression", "LIGER")

# Run scCombination
comb.results <- runScCombination(
  single.savePaths = single.savePaths,
  sampleNames = sampleNames,
  savePath = savePath,
  combName = combName,
  authorName = authorName,
  comb.method = comb.method
)
