

source("Function-L0REG.R")
data <- read.csv("../scGRN-L0_data/BEELINE-data/inputs/scRNAseq_preprocessing/data/GSE75748/ExpressionData.csv", row.name = 1)
cellTypes <- read.csv("../scGRN-L0_data/BEELINE-data/inputs/scRNAseq_preprocessing/data/GSE75748/CellType.csv")
tfs <- read.csv("../scGRN-L0_data/BEELINE-Networks/human-tfs.csv")
cellTypesList <- unique(cellTypes$Time)
dataCellType <- data[tfs$TF[150:200], grep(cellTypesList[1], colnames(data))] %>% na.omit()
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
