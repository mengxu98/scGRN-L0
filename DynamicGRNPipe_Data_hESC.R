

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

CellInfor.trajectory <- cbind.data.frame(CellInfor, PseTime.Lineage1=slingshot$PseudoTime)
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
binpoint <- c(0, 
              max(cellInfor$PseTime.Lineage1)/5, 
              max(cellInfor$PseTime.Lineage1)*2/5, 
              max(cellInfor$PseTime.Lineage1)*3/5, 
              max(cellInfor$PseTime.Lineage1)*4/5, 
              max(cellInfor$PseTime.Lineage1))
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
dataCellType <- data[tfs$TF[1:500], grep(cellTypesList[1], colnames(data))] %>% na.omit()
L0REG_L0 <- L0REG(
  matrix = dataCellType,
  # regulators = tfs$TF,
  # targets = colnames(data_GENIE3),
  penalty = "L0"
)
write.table(L0REG_L0,
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
inputExpr <- dataCellType
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
ground_truth_m(
  intput = paste0("../scGRN-L0_output/output_scRNA-Seq/GRN_PPCOR.txt"),
  output = "../scGRN-L0_output/output_scRNA-Seq/",
  dataset_dir = mnetwork_data_dir,
  database = "All"
)
evaluationObject <- prepareEval(paste0("../scGRN-L0_output/output_scRNA-Seq/GRN_PPCOR.txt"),
  paste0(paste0("../scGRN-L0_output/output_scRNA-Seq/ground_truth.tsv")),
  totalPredictionsAccepted = 100000
)
AUROC_PPCOR <- calcAUROC(evaluationObject)
AUPRC_PPCOR <- calcAUPR(evaluationObject)
AUROC_PPCOR
