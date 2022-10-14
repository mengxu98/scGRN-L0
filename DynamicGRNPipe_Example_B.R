

saveRDS(seu_obj_data_B_samples, file ="../scGRN-L0_data/seu_obj_data_B_samples.rds")
data <- readRDS("../scGRN-L0_data/seu_obj_data_B_samples.rds")
meta <- data@meta.data
CellInfor_B <- data.frame(UniqueCell_ID = colnames(data),
                        Patient = meta$orig.ident,
                        majorCluster= meta$celltype,
                        sampleType="stage-1")
table(CellInfor_B$majorCluster)
GSE131907_B <- as.matrix(data@assays$SCT@counts) %>% as.data.frame()
# Example
# The example data can be downloaded at https://www.jianguoyun.com/p/DeYtp_AQ1bXiCRjXvYwE
# ====0.input====
load(file = "DynamicGRNPipe_ExampleData/clusterSig.RData") # genes used to construct cell trajectories
load(file = "../scGRN-L0_data/GSE131907_B.RData") # expression profile and cell annotation #from GSE99254

# ==== 1. Construction of cell state transformation trajectory (slingshot) ====
library(slingshot)
source("DynamicGRNPipe_1.slingshot_B.R")
source("DynamicGRNPipe_Function.R")
GSE131907_B <- as.matrix(GSE131907_B)
t.slingshot <- slingshot_run(scRNAseq.Exp = GSE131907_B,
  clusterLabels = CellInfor_B$majorCluster,
  ordergene = unlist(clusterSig),
  RMmethod = "pca",
  plot_output = TRUE
  #start.cluster = "GrB-secreting B cells" # cannot determine what type of cell types as the first type
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
  a = cellInfor[cellInfor$majorCluster == "B_mature", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "B_naive", "PseTime.Lineage1"],
  c = cellInfor[cellInfor$majorCluster == "B_plasma", "PseTime.Lineage1"],
  d = cellInfor[cellInfor$majorCluster == "B_plasmablast", "PseTime.Lineage1"],
  filename = "Results/C5_C61.png"
)
# Manually select the intersection point according to the density plot
binpoint <- c(0, 16.4, 26.29, 35.14, 43.51, max(cellInfor$PseTime.Lineage1))

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
source("DynamicGRNPipe_3.constructionNetwork.R")
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
  positivecor = 0#,
  #confidence = NULL
)
lineage1dynet1 <- lineage1dynet[1:2]
# extract active edges for each window
Dynnet_active1 <- lapply(lineage1dynet, function(x) {
  rownames(x) <- paste0(x[, 1], "_", x[, 2])
  x <- x[x[, "spearmanCor"] > 0, ]
  return(x)
})
names(Dynnet_active1) <- paste0("W", 1:length(Dynnet_active1))

source("ground-truth.R")
source("framework_main.R")
ground_truth_T(lineage1dynet[[3]], dorothea_regulon_human)
evaluationObject <- prepareEval("ground_pred.txt",
  paste0("ground_truth.tsv"),
  totalPredictionsAccepted = 100000
)

GENIE3_AUROC <- calcAUROC(evaluationObject)
GENIE3_AUPR <- calcAUPR(evaluationObject)

weightofWindows_L0 <- DynNet_L0(
  Windows = Windows1[3],
  CD8TCellExp.trajectory = CD8TCellExp.trajectory,
  DynamicGene = DynamicGene1, # set Background genes,which used to construct the network, such as highly variable genes, dynamic genes along trajectory
  allTFs = allTFs[1:200], # set regulators
  detectNum = 10, detectPro = 0.05, meanExp = 1 # Noise filtering threshold
)

ground_truth_T(weightofWindows_L0[[1]], dorothea_regulon_human)
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

# extract active edges for each window
Dynnet_active1_L0 <- lapply(lineage1dynet_L0, function(x) {
  rownames(x) <- paste0(x[, 1], "_", x[, 2])
  x <- x[x[, "spearmanCor"] > 0, ]
  return(x)
})
names(Dynnet_active1_L0) <- paste0("W", 1:length(Dynnet_active1_L0))

source("ground-truth.R")
source("framework_main.R")
ground_truth_T(weightofWindows_L0[[2]], dorothea_regulon_human)
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