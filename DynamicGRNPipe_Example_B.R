

library(Seurat)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(patchwork)
library(reshape2)

if (F) {
  Cellratio <- prop.table(table(seu_obj_data$main_cell_type, seu_obj_data$orig.ident),
    margin = 2
  ) %>% as.data.frame()
  colourCount <- length(unique(Cellratio$Var1))
  ggplot(Cellratio) +
    geom_bar(aes(x = Var2, y = Freq, fill = Var1),
      stat = "identity",
      width = 0.8,
      size = 1,
      colour = "black"
    ) +
    labs(x = "Sample", y = "Ratio") +
    coord_flip() +
    theme_bw()
  ggsave2(paste0("Fig6.contrast1.png"),
    path = paste0("Results/"),
    width = 15, height = 10, units = "cm"
  )

  cellper <- reshape2::dcast(Cellratio, Var2 ~ Var1, value.var = "Freq")
  rownames(cellper) <- cellper[, 1]
  cellper <- cellper[, -1]
  sample <- seu_obj_data$orig.ident
  group <- seu_obj_data$tissue_type
  samples <- data.frame(sample, group)
  cellper$sample <- rownames(cellper)
  cellper$group <- ""
  for (i in 1:nrow(cellper)) {
    stage <- samples[which(samples$sample == cellper$sample[i])[1], "group"]
    cellper$group[i] <- stage
  }
  sce_groups <- c("Follicular B cells", "Germinal center B cells", "Plasma")
  pplist <- list()
  for (group_ in sce_groups) {
    cellper_ <- cellper %>% select(one_of(c("sample", "group", group_)))
    colnames(cellper_) <- c("sample", "group", "percent")
    cellper_$percent <- as.numeric(cellper_$percent)
    cellper_ <- cellper_ %>%
      group_by(group) %>%
      mutate(
        upper = quantile(percent, 0.75),
        lower = quantile(percent, 0.25),
        mean = mean(percent),
        median = median(percent)
      )
    pp1 <- ggplot(cellper_, aes(x = group, y = percent)) +
      geom_jitter(shape = 21, aes(fill = group), width = 0.25) +
      stat_summary(fun = mean, geom = "point", color = "grey60") +
      theme(
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(size = 4, face = "plain"),
        legend.position = "none"
      ) +
      labs(title = group_, y = "Percentage") +
      geom_errorbar(aes(ymin = lower, ymax = upper), col = "grey60", width = 0.5) +
      theme_bw()

    labely <- max(cellper_$percent)
    compare_means(percent ~ group, data = cellper_)
    my_comparisons <- list(c("Normal", "Tumor"))
    pp1 <- pp1 + stat_compare_means(comparisons = my_comparisons, size = 3, method = "t.test")
    pplist[[group_]] <- pp1
  }

  pplist[["Follicular B cells"]] +
    pplist[["Germinal center B cells"]] +
    pplist[["Plasma"]] +
    plot_annotation(tag_levels = "a") +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave2(paste0("Fig6.contrast.png"),
    path = paste0("Results/"),
    width = 25, height = 15, units = "cm"
  )
}

GSE131907_B <- as.matrix(seu_obj_data@assays$RNA@counts)
meta <- seu_obj_data@meta.data
CellInfor_B <- data.frame(
  UniqueCell_ID = colnames(seu_obj_data),
  Patient = meta$orig.ident,
  majorCluster = meta$main_cell_type,
  sampleType = meta$tissue_type
)
# ====0.input====
load(file = "DynamicGRNPipe_ExampleData/clusterSig.RData") # genes used to construct cell trajectories
# load(file = "../scGRN-L0_data/GSE131907_B.RData") # expression profile and cell annotation #from GSE99254

# ==== 1. Construction of cell state transformation trajectory (slingshot) ====
library(slingshot)
source("DynamicGRNPipe_1.slingshot_B.R")
source("DynamicGRNPipe_Function.R")
t.slingshot <- slingshot_run(
  scRNAseq.Exp = GSE131907_B,
  clusterLabels = CellInfor_B$majorCluster,
  ordergene = unlist(clusterSig),
  RMmethod = "pca",
  plot_output = TRUE,
  start.cluster = "Follicular B cells" # cannot determine what type of cell types as the first type
)
CellInfor.trajectory <- cbind.data.frame(CellInfor_B, t.slingshot$data)
CD8TCellExp.trajectory <- GSE131907_B

# ====2.slinding windows based on pseudotime and anotation of cells====
source("DynamicGRNPipe_2.CellWindowing.R")
cells <- rownames(CellInfor.trajectory)[grep("1", CellInfor.trajectory$Branch)] # cells in lineage1
CellInfor.trajectory$UniqueCell_ID <- rownames(CellInfor.trajectory) ########
cellInfor <- CellInfor.trajectory[cells, c("UniqueCell_ID", "majorCluster", "PseTime.Lineage1")]
table(CellInfor_B$majorCluster)

# Calculate all intersections of pseudo-time density curves of cells in different states
C1_C2_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "B_naive", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "B_mature", "PseTime.Lineage1"],
  filename = "Results/C1_C2.png"
)
C2_C3_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "B_mature", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "B_plasma", "PseTime.Lineage1"],
  filename = "Results/C2_C3.png"
)
C3_C4_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "B_plasma", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "B_plasmablast", "PseTime.Lineage1"],
  filename = "Results/C3_C4.png"
)
C4_C5_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "B_plasmablast", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "GC B cells in the DZ", "PseTime.Lineage1"],
  filename = "Results/C4_C5.png"
)

C5_C61_cross <- densityintersections(
  a = cellInfor[cellInfor$majorCluster == "Follicular B cells", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "Germinal center B cells", "PseTime.Lineage1"],
  c = cellInfor[cellInfor$majorCluster == "Plasma", "PseTime.Lineage1"],
  filename = "Results/all.png"
)
# Manually select the intersection point according to the density plot
binpoint <- c(
  0,
  max(cellInfor$PseTime.Lineage1) / 5,
  max(cellInfor$PseTime.Lineage1) * 2 / 5,
  max(cellInfor$PseTime.Lineage1) * 3 / 5,
  max(cellInfor$PseTime.Lineage1) * 4 / 5,
  max(cellInfor$PseTime.Lineage1)
)

# Two consecutive bins of cells as a window
# extract cells for each window
t1 <- cut(cellInfor$PseTime.Lineage1, binpoint)
bincells <- split(cellInfor$UniqueCell_ID, t1)
W1 <- unlist(bincells[1:2])
W2 <- unlist(bincells[2:3])
W3 <- unlist(bincells[3:4])
W4 <- unlist(bincells[4:5])
Windows1 <- list(W1 = W1, W2 = W2, W3 = W3, W4 = W4)


# ====3.constructing dynamic networks (GENIE3)======
source("DynamicGRNPipe_3.constructionNetwork_L0.R")
load("DynamicGRNPipe_ExampleData/DynamicGene1.RData")
load("DynamicGRNPipe_ExampleData/dorothea.RData")
# Construct a control network and calculate the control weight of each edge
weightofWindows <- DynNet_RF(
  Windows = Windows1,
  CD8TCellExp.trajectory = CD8TCellExp.trajectory,
  DynamicGene = DynamicGene1, # set Background genes,which used to construct the network, such as highly variable genes, dynamic genes along trajectory
  allTFs = allTFs, # set regulators
  detectNum = 10, detectPro = 0.05, meanExp = 1 # Noise filtering threshold
)

# Filter the edges for each window
lineage1dynet <- DynNet_filter(
  Windows = Windows1,
  CD8TCellExp.trajectory = CD8TCellExp.trajectory,
  weightofWindows = weightofWindows,
  weightThr = 0.02,
  nsd = 2,
  positivecor = 0 # ,
  # confidence = NULL
)
lineage1dynet1 <- lineage1dynet[1:2]
length(lineage1dynet1[[4]])
# extract active edges for each window
Dynnet_active1 <- lapply(lineage1dynet, function(x) {
  rownames(x) <- paste0(x[, 1], "_", x[, 2])
  x <- x[x[, "spearmanCor"] > 0, ]
  return(x)
})
names(Dynnet_active1) <- paste0("W", 1:length(Dynnet_active1))

source("ground-truth.R")
source("framework_main.R")
ground_truth_T(weightofWindows[[1]], dorothea_regulon_human)
evaluationObject <- prepareEval("ground_pred.txt",
  paste0("ground_truth.tsv"),
  totalPredictionsAccepted = 100000
)

GENIE3_AUROC <- calcAUROC(evaluationObject)
GENIE3_AUPR <- calcAUPR(evaluationObject)
GENIE3_AUROC

weightofWindows_L0 <- DynNet_L0(
  Windows = Windows1,
  CD8TCellExp.trajectory = CD8TCellExp.trajectory,
  DynamicGene = DynamicGene1, # set Background genes,which used to construct the network, such as highly variable genes, dynamic genes along trajectory
  allTFs = allTFs, # set regulators
  detectNum = 10, detectPro = 0.05, meanExp = 1 # Noise filtering threshold
)

# Filter the edges for each window
lineage1dynet_L0 <- DynNet_filter(
  Windows = Windows1,
  CD8TCellExp.trajectory = CD8TCellExp.trajectory,
  weightofWindows = weightofWindows_L0,
  weightThr = 0.02,
  nsd = 2,
  positivecor = 0 # ,
  # confidence = NULL
)
lineage1dynet_L0[[4]]
lineage1dynet1_L0 <- lineage1dynet_L0[1:2]
length(lineage1dynet1_L0[[4]])
# extract active edges for each window
Dynnet_active1_L0 <- lapply(lineage1dynet_L0, function(x) {
  rownames(x) <- paste0(x[, 1], "_", x[, 2])
  x <- x[x[, "spearmanCor"] > 0, ]
  return(x)
})
names(Dynnet_active1_L0) <- paste0("W", 1:length(Dynnet_active1_L0))

ground_truth_T(weightofWindows_L0[[1]], dorothea_regulon_human)
evaluationObject_L0 <- prepareEval("ground_pred.txt",
  paste0("ground_truth.tsv"),
  totalPredictionsAccepted = 100000
)

L0_AUROC <- calcAUROC(evaluationObject_L0)
L0_AUPR <- calcAUPR(evaluationObject_L0)
L0_AUROC # 0.4601513 no parallel
L0_AUPR

GENIE3_AUROC
GENIE3_AUPR

# Filter the edges for each window
lineage1dynet_L0 <- DynNet_filter(
  Windows = Windows1,
  CD8TCellExp.trajectory = CD8TCellExp.trajectory,
  weightofWindows = weightofWindows_L0,
  weightThr = 0.02,
  nsd = 2,
  positivecor = 0 # ,
  # confidence = NULL
)

ground_truth_T(weightofWindows_L0[[3]], dorothea_regulon_human)
evaluationObject_L0 <- prepareEval("ground_pred.txt",
  paste0("ground_truth.tsv"),
  totalPredictionsAccepted = 100000
)

L0_AUROC <- calcAUROC(evaluationObject_L0)
L0_AUPR <- calcAUPR(evaluationObject_L0)

GENIE3_AUROC
GENIE3_AUPR
L0_AUROC
L0_AUPR
