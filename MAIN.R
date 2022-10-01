

###                      SINCERITIES MAIN SCRIPT                        ###

# Please prepare the data in an Excel sheet and then save it as csv file
# using the format below
#
# Data: s-by-m+1 matrix, where s is the total number of observations/single
# cells and m is the number of genes. The first m columns contain
# the expression level of each m genes, and the last column contains the the time-stamps.
#
# Two data formats are accepted:
#
# A) with row header
# ------------------------------------
# Gene1  Gene2  Gene3  ... Genej  Time
#  27     80     56    ...  69      0
#  73     20     90    ...  45      0
#   .     .      .     ...  .       .
#   .     .      .     ...  .       .
#   .     .      .     ...  .       .
# ------------------------------------
#
# B) without row header
# ------------------------------------
#  27     80     56    ...  69      0
#  73     20     90    ...  45      0
#   .     .      .     ...  .       .
#   .     .      .     ...  .       .
#   .     .      .     ...  .       .
# ------------------------------------
#
# Please specify row header in funcion uploading. If row header is absent
# set the header argument in uploading() to false.
#
# PACKAGES required:
# kSamples
# glmnet
# ppcor

# --------------------------------------------------
library(tidyr)
library(kSamples)
library(glmnet)
library(ppcor)
library(L0Learn)
# Load all the script files
source("framework_main.R")
source("ground-truth.R")

simulation_data_dir <- "../scGRN-L0_data/BEELINE-data/inputs/Curated/mCAD/mCAD-2000-1-70/"
simulation_data_file <- "ExpressionData.csv"
simulation_PseudoTime_file <- "PseudoTime.csv"

simulation_data <- read.csv(paste0(simulation_data_dir, simulation_data_file), row.names = 1)
simulation_PseudoTime <- read.csv(paste0(simulation_data_dir, simulation_PseudoTime_file), row.names = 1)
if (ncol(simulation_PseudoTime) > 1) {
  simulation_data_news <- c()
  for (i in 1:ncol(simulation_PseudoTime)) {
    simulation_data_new <- cbind.data.frame(t(simulation_data), h = simulation_PseudoTime[, i]) %>% na.omit()
    write.csv(simulation_data_new, paste0(simulation_data_dir, "ExpressionData_", i, ".csv"), row.names = FALSE)
    simulation_data_news <- rbind.data.frame(simulation_data_news, simulation_data_new)
  }
  write.csv(simulation_data_news, paste0(simulation_data_dir, "ExpressionData_all", ".csv"), row.names = FALSE)
}

# *** Data loading ***
uploading <- dget("SINCERITIES functions/uploading.R")
# DATA <- uploading()
# DATA <- uploading(paste0(simulation_data_dir, "ExpressionData_", 1, ".csv"))
DATA <- uploading(paste0(simulation_data_dir, "ExpressionData_all", ".csv"))
# *** SINCERITIES ***

# SINCERITITES <- dget("SINCERITIES functions/SINCERITIES.R")
# result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 1)

SINCERITITES_L0 <- dget("SINCERITIES functions/SINCERITIES_L0.R")
result_L0 <- SINCERITITES_L0(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 1, penalty = "L0L2", algorithm = "CD")

# Final ranked list of regulatory edges
adj_matrix_L0 <- result_L0$adj_matrix / max(result_L0$adj_matrix)
final_ranked_predictions <- dget("SINCERITIES functions/final_ranked_predictions.R")
table_L0 <- final_ranked_predictions(adj_matrix_L0, DATA$genes, SIGN = 1, saveFile = TRUE)
table_L0 <- as.data.frame(table_L0)
write.table(table_L0[, -4], "output/test_L0.txt", row.names = F, col.names = F, sep = "\t", quote = F)

# test
ground_truth_simulation(
  intput = "output/test_L0.txt",
  output = "output/",
  dataset_dir = simulation_data_dir,
  file = "refNetwork.csv"
)
evaluationObject <- prepareEval("output/test_L0.txt",
  paste0("output/ground_truth.tsv"),
  totalPredictionsAccepted = 100000
)

L0_AUROC <- calcAUROC(evaluationObject)
L0_AUPR <- calcAUPR(evaluationObject)
L0_AUROC
L0_AUPR

a <- convertSortedRankTSVToAdjMatrix("output/ground_truth.tsv")
a
library(pracma)
auc_from_ranks_TC_sign <- dget("SINCERITIES functions/auc_from_ranks_TC_sign.R")
AUCresult_L0 <- auc_from_ranks_TC_sign(adj_matrix_L0, a, 1000)
AUROC <- AUCresult_L0$AUROC
AUPR <- AUCresult_L0$AUPR
AUROC
AUPR

# --------------------------------------------------
SINCERITITES <- dget("SINCERITIES functions/SINCERITIES.R")
result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 1)

# Final ranked list of regulatory edges
adj_matrix <- result$adj_matrix / max(result$adj_matrix)
final_ranked_predictions <- dget("SINCERITIES functions/final_ranked_predictions.R")
table <- final_ranked_predictions(adj_matrix, DATA$genes, SIGN = 1, saveFile = TRUE)
table <- as.data.frame(table)
write.table(table[, -4], "output/test.txt", row.names = F, col.names = F, sep = "\t", quote = F)

# test
ground_truth_simulation(
  intput = "output/test.txt",
  output = "output/",
  dataset_dir = simulation_data_dir,
  file = "refNetwork.csv"
)

evaluationObject <- prepareEval("output/test.txt",
  paste0("output/ground_truth.tsv"),
  totalPredictionsAccepted = 100000
)

a <- convertSortedRankTSVToAdjMatrix("output/ground_truth.tsv")
a
library(pracma)
auc_from_ranks_TC_sign <- dget("SINCERITIES functions/auc_from_ranks_TC_sign.R")
AUCresult <- auc_from_ranks_TC_sign(adj_matrix, a, 1000)
AUROC <- AUCresult$AUROC
AUPR <- AUCresult$AUPR
AUROC
AUPR

SINCERITITES_AUROC <- calcAUROC(evaluationObject)
SINCERITITES_AUPR <- calcAUPR(evaluationObject)
SINCERITITES_AUROC
SINCERITITES_AUPR

# --------------------------------------------------

# expression_dataset_test <- constructExpressionMatrixFromFile(paste0(pathway, dataway))
# expression_dataset_test <- read.table(paste0(pathway, dataway), header = T) %>% as.matrix()
expression_dataset_test <- read.csv(paste0(simulation_data_dir, "ExpressionData_all.csv"),
  header = T
) %>%
  as.matrix()

PseudoTime <- expression_dataset_test[, ncol(expression_dataset_test)]
# PseudoTime <- PseudoTime[order(PseudoTime), ] %>% as.data.frame()
# rownames(PseudoTime) <- PseudoTime$X

expression_dataset_test <- expression_dataset_test[, -ncol(expression_dataset_test)]
# expression_dataset_test <- data_filetr(expression_dataset_test,
#   dataset_dir = "../scGRN-L0_data/BEELINE-Networks/Networks/human/",
#   database = "All"
# )

# TFs <- read.csv("../scGRN-L0_data/BEELINE-Networks/human-tfs.csv")

# expression_dataset_test <- expression_dataset_test[rownames(PseudoTime)[1:50], intersect(TFs$TF, colnames(expression_dataset_test))]

# expression_dataset_test <- expression_dataset_test[, 1:100]
dim(expression_dataset_test)
min(expression_dataset_test)

# Run algorithm using default settings, writes result to output.txt
NIMEFI(expression_dataset_test,
  GENIE = T, SVM = F, EL = F, penalty = "L0",
  outputFileName = "output/output_GENIE3",
  outputFileFormat = "txt",
  SVMRankThreshold = 5, SVMEnsembleSize = 100,
  ELPredSampleMin = 20, ELPredSampleMax = 80,
  ELExpSampleMin = 20, ELExpSampleMax = 80,
  ELRankThreshold = 5, ELEnsembleSize = 200
)

evaluationObject <- prepareEval("output/output_GENIE3.txt",
  paste0("output/ground_truth.tsv"),
  totalPredictionsAccepted = 100000
)

GENIE3_AUROC <- calcAUROC(evaluationObject)
GENIE3_AUPR <- calcAUPR(evaluationObject)
GENIE3_AUROC
GENIE3_AUPR
L0_AUROC
L0_AUPR
