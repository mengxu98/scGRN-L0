

library(Seurat)
library(destiny)
library(slingshot)
library(plotly)
library(gam)
library(RColorBrewer)
library(EBSeq)
library(Rtsne)
load("../scGRN-L0_data/data/seu_obj_B.Rdata")

# Read input files. The input file is from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75748
# Chu2017 <- read.table('../scGRN-L0_data/BEELINE-data/inputs/scRNAseq_preprocessing/data/GSE75748/GSE75748_sc_time_course_ec.csv', sep = ',', header = T, row.names = 1)
Chu2017 <- as.matrix(seu_obj_data@assays$RNA@counts)
Chu2017 <- as.matrix(seu_obj_data@assays$RNA@data)
# ChuCellTypes <- data.frame()
# for (col in colnames(Chu2017)){
#   cType <- strsplit(strsplit(col,'[.]')[[1]][2],'b4s_|_')[[1]][1]
#   ChuCellTypes[col,'Time'] <- cType
# }
ChuCellTypes <- seu_obj_data$main_cell_type
head(Chu2017)
head(ChuCellTypes)
# Compute diffusion map projection of the cells
# The scRNA-seq data were normalized by median-by-ratio normalization
Sizes <- MedianNorm(Chu2017)
if (is.na(Sizes[1])) {
  Sizes <- MedianNorm(Chu2017, alternative = TRUE)
  message("alternative normalization method is applied")
}
Chu2017MedNorm <- GetNormalizedMat(Chu2017, Sizes)
logExpression <- log2(Chu2017MedNorm + 1)
geneFiltered <- apply(Chu2017, 1, function(x) {
  sum(x >= 1) >= 1
})
logExpressionFiltered <- logExpression[which(geneFiltered), ]
logExpressionFiltered <- Chu2017
# Compute PCA to identify informative genes
pcaRes <- prcomp(t(logExpressionFiltered), scale. = FALSE)
dmapRes <- DiffusionMap(t(logExpressionFiltered), distance = "cosine", sigma = 0.25, k = 100)
rd1 <- as.data.frame(cbind(PC1 = pcaRes$x[, 1], PC2 = pcaRes$x[, 2], PC3 = pcaRes$x[, 3]))
rd2 <- as.data.frame(cbind(DC1 = dmapRes$DC1, DC2 = dmapRes$DC2))
plot_ly(as.data.frame(pcaRes$x), x = ~PC1, y = ~PC2, color = ChuCellTypes$Time, colors = brewer.pal(6, "Set1"))
plot_ly(rd2, x = ~DC1, y = ~DC2, color = ChuCellTypes$Time, colors = brewer.pal(6, "Set1"))
# Run slingshot
slingshotPT <- slingshot(rd2,
  reducedDim = rd2,
  clusterLabels = ChuCellTypes$Time, start.clus = "00h", end.clus = "96h"
)
ssPT <- slingPseudotime(slingshotPT)
ssPT <- as.data.frame(ssPT)
# plot_ly(rd2, x=~DC1, y= ~DC2,  color = ssPT$curve1)
# plot_ly(rd1, x=~PC1, y= ~PC2,  color = ssPT$curve1)
plot_ly(rd2, x = ~DC1, y = ~DC2, color = ssPT$Lineage1)
plot_ly(rd1, x = ~PC1, y = ~PC2, color = ssPT$Lineage1)
# t <- ssPT$curve1
t <- ssPT$Lineage1

# look at top variable gens
Y <- logExpressionFiltered
var1K <- names(sort(apply(Y, 1, var), decreasing = TRUE))
Y <- Y[var1K, ]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y, 1, function(z) {
  d <- data.frame(z = z, t = t)
  suppressWarnings({
    tmp <- gam(z ~ lo(t), data = d)
  })
  p <- summary(tmp)[4][[1]][1, 5]
  p
})
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:15]
heatdata <- logExpressionFiltered[topgenes, order(t, na.last = NA)]
cTypes <- as.data.frame(ChuCellTypes$Time)
rownames(cTypes) <- colnames(Chu2017)
heatclus <- as.factor(cTypes[order(t, na.last = NA), ])

heatmap(as.matrix(heatdata),
  Colv = NA,
  ColSideColors = brewer.pal(6, "Set1")[heatclus], labCol = FALSE
)
exprData <- logExpressionFiltered
colnames(exprData) <- gsub(pattern = "[.]", replacement = "_", colnames(exprData))
# ptData <- data.frame(ssPT$curve1)
ptData <- data.frame(ssPT$Lineage1)

rownames(ptData) <- colnames(exprData)
colnames(ptData) <- "PseudoTime"

geneData <- data.frame(sort(gam.pval, decreasing = FALSE))
colnames(geneData) <- "VGAMpValue"
geneData[, "Variance"] <- apply(logExpressionFiltered[rownames(geneData), ], 1, var)

write.csv(x = exprData, file = "data/GSE75748/ExpressionData.csv", quote = FALSE)
write.csv(x = ptData, file = "data/GSE75748/PseudoTime.csv", quote = FALSE)
write.csv(x = geneData, file = "data/GSE75748/GeneOrdering.csv", quote = FALSE)
write.csv(x = ChuCellTypes, file = "data/GSE75748/CellType.csv", quote = FALSE)

rdDF <- as.data.frame(rd2)

plot_ly(rdDF, x = ~DC1, y = ~DC2, color = ptData$PseudoTime)
# plot_ly(as.data.frame(rd2), x=~DC1, y= ~DC2, color = ssPT$curve1)
plot_ly(as.data.frame(rd2), x = ~DC1, y = ~DC2, color = ssPT$Lineage1)


source("Function-L0REG.R")
data <- read.csv("../scGRN-L0_data/BEELINE-data/inputs/scRNAseq_preprocessing/data/GSE75748/ExpressionData.csv", row.name = 1)
cellTypes <- read.csv("../scGRN-L0_data/BEELINE-data/inputs/scRNAseq_preprocessing/data/GSE75748/CellType.csv")
tfs <- read.csv("../scGRN-L0_data/BEELINE-Networks/human-tfs.csv")
slingshot <- read.csv("../scGRN-L0_data/BEELINE-data/inputs/scRNAseq_preprocessing/data/GSE75748/PseudoTime.csv")
CellInfor <- data.frame(
  UniqueCell_ID = colnames(data),
  Patient = cellTypes,
  majorCluster = cellTypes$Time,
  sampleType = "hESC"
)

CellInfor.trajectory <- cbind.data.frame(CellInfor, PseTime.Lineage1 = slingshot$PseudoTime)
head(CellInfor.trajectory)
# ====2.slinding windows based on pseudotime and anotation of cells====
source("DynamicGRNPipe_2.CellWindowing.R")
# cells <- rownames(CellInfor.trajectory)[grep("1", CellInfor.trajectory$Branch)] # cells in lineage1
cellInfor <- CellInfor.trajectory[, c("UniqueCell_ID", "majorCluster", "PseTime.Lineage1")]
head(cellInfor)

# Calculate all intersections of pseudo-time density curves of cells in different states
C1_C2_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "00h", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "12h", "PseTime.Lineage1"],
  filename = "C1_C2.png"
)
C2_C4_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "CD8_C2-CD28", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "CD8_C4-GZMK", "PseTime.Lineage1"],
  filename = "C2_C4.png"
)
C4_C5_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "CD8_C4-GZMK", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "CD8_C5-ZNF683", "PseTime.Lineage1"],
  filename = "C4_C5.png"
)
C5_C6_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "CD8_C5-ZNF683", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "CD8_C6-LAYN", "PseTime.Lineage1"],
  filename = "C5_C6.png"
)

# Manually select the intersection point according to the density plot
# binpoint <- c(0, 16.4, 26.29, 35.14, 43.51, max(cellInfor$PseTime.Lineage1))
C5_C61_cross <- densityintersections(
  a = cellInfor[cellInfor$majorCluster == "00h", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "12h", "PseTime.Lineage1"],
  c = cellInfor[cellInfor$majorCluster == "24h", "PseTime.Lineage1"],
  d = cellInfor[cellInfor$majorCluster == "36h", "PseTime.Lineage1"],
  e = cellInfor[cellInfor$majorCluster == "72h", "PseTime.Lineage1"],
  f = cellInfor[cellInfor$majorCluster == "96h", "PseTime.Lineage1"],
  filename = "../scGRN-L0_output/all.png"
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



cellTypesList <- unique(cellTypes$Time)
dataCellTypes <- data[tfs$TF[1:500], ] %>% na.omit()

for (i in 1:length(cellTypesList)) {
  dataCellType <- data[tfs$TF[1:500], grep(cellTypesList[i], colnames(data))] %>% na.omit()
  L0REG_L0 <- L0REG(
    matrix = dataCellType,
    # regulators = tfs$TF,
    # targets = colnames(data_GENIE3),
    penalty = "L0"
  )
  if (i == 1) {
    L0REG_L0_list <- L0REG_L0
  } else {
    L0REG_L0_list$weight <- L0REG_L0_list$weight + L0REG_L0$weight
  }
}
L0REG_L0_list <- L0REG_L0_list[order(L0REG_L0_list$weight, decreasing = TRUE), ]

write.table(L0REG_L0_list,
  paste0("../scGRN-L0_output/output_scRNA-Seq/output_L0GRN.txt"),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = F
)

ground_truth_h(
  intput = paste0("../scGRN-L0_output/output_scRNA-Seq/output_L0GRN.txt"),
  output = "../scGRN-L0_output/output_scRNA-Seq/",
  dataset_dir = hnetwork_data_dir,
  database = "All"
)

evaluationObject <- prepareEval(paste0("../scGRN-L0_output/output_scRNA-Seq/output_L0GRN.txt"),
  paste0(paste0("../scGRN-L0_output/output_scRNA-Seq/ground_truth.tsv")),
  totalPredictionsAccepted = 100000
)

AUROC_L0REG_L0_N <- calcAUROC(evaluationObject)
AUPR_L0REG_L0_N <- calcAUPR(evaluationObject)
AUROC_L0REG_L0_N



library(ppcor)
inputExpr <- dataCellTypes
geneNames <- rownames(inputExpr)
rownames(inputExpr) <- c(geneNames)
pcorResults <- pcor(x = t(as.matrix(inputExpr)), method = "spearman")
DF <- data.frame(
  Gene1 = geneNames[c(row(pcorResults$estimate))], Gene2 = geneNames[c(col(pcorResults$estimate))],
  corVal = c(pcorResults$estimate), pValue = c(pcorResults$p.value)
)
outDF <- DF[order(DF$corVal, decreasing = TRUE), ]
outDF <- outDF[-(1:length(geneNames)), ]
write.table(outDF,
  paste0("../scGRN-L0_output/output_scRNA-Seq/GRN_PPCOR.txt"),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = F
)
ground_truth_h(
  intput = paste0("../scGRN-L0_output/output_scRNA-Seq/GRN_PPCOR.txt"),
  output = "../scGRN-L0_output/output_scRNA-Seq/",
  dataset_dir = hnetwork_data_dir,
  database = "All"
)
evaluationObject <- prepareEval(paste0("../scGRN-L0_output/output_scRNA-Seq/GRN_PPCOR.txt"),
  paste0(paste0("../scGRN-L0_output/output_scRNA-Seq/ground_truth.tsv")),
  totalPredictionsAccepted = 100000
)
AUROC_PPCOR <- calcAUROC(evaluationObject)
AUPRC_PPCOR <- calcAUPR(evaluationObject)
AUROC_PPCOR

library(LEAP)
geneNames <- rownames(inputExpr)
rownames(inputExpr) <- c()
# Run LEAP's compute Max. Absolute Correlation
# MAC_cutoff is set to zero to get a score for all TFs
# max_lag_prop is set to the max. recommended value from the paper's supplementary file
# Link to paper: https://academic.oup.com/bioinformatics/article/33/5/764/2557687
MAC_results <- MAC_counter(
  data = inputExpr, # max_lag_prop=maxLag,
  MAC_cutoff = 0,
  file_name = "temp", lag_matrix = FALSE, symmetric = FALSE
)
Gene1 <- geneNames[MAC_results[, "Row gene index"]]
Gene2 <- geneNames[MAC_results[, "Column gene index"]]
Score <- MAC_results[, "Correlation"]
outDF <- data.frame(Gene1, Gene2, Score)
write.table(outDF, paste0("../scGRN-L0_output/output_scRNA-Seq/GRN_LEAP.txt"),
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = F
)
evaluationObject <- prepareEval(paste0("../scGRN-L0_output/output_scRNA-Seq/GRN_LEAP.txt"),
  paste0(paste0("../scGRN-L0_output/output_scRNA-Seq/ground_truth.tsv")),
  totalPredictionsAccepted = 100000
)
AUROC_LEAP <- calcAUROC(evaluationObject)
AUPRC_LEAP <- calcAUPR(evaluationObject)
AUROC_LEAP

library(GENIE3)
weightMat <- GENIE3(
  exprMatrix = as.matrix(dataCellTypes),
  nCores = 32
)
weightdf <- getLinkList(weightMat)
names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
write.table(weightdf, file = paste0("../scGRN-L0_output/output_scRNA-Seq/GRN_GENIE3.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
evaluationObject <- prepareEval(paste0("../scGRN-L0_output/output_scRNA-Seq/GRN_GENIE3.txt"),
  paste0(paste0("../scGRN-L0_output/output_scRNA-Seq/ground_truth.tsv")),
  totalPredictionsAccepted = 100000
)
AUROC_GENIE3 <- calcAUROC(evaluationObject)
AUPRC_GENIE3 <- calcAUPR(evaluationObject)
AUROC_GENIE3
