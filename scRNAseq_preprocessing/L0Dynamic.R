
library("tidyr")
library("kSamples")
library("glmnet")
library("ppcor")
library("L0Learn")
library("pracma")
library("ggplot2")
library("reshape2")
library("RColorBrewer")
library("patchwork")
source("G:/Research/GRN/scGRN-L0/framework_main.R")
source("G:/Research/GRN/scGRN-L0/ground-truth.R")
source("G:/Research/GRN/scGRN-L0/Function-L0REG.R")
source("G:/Research/GRN/scGRN-L0/basic_evaluation.R")

L0Dynamic <- L0REG(t(data_grn),
                   regulators = colnames(data_grn),
                   targets = colnames(data_grn),
                   maxSuppSize = 2,
                   penalty = "L0"
)
write.table(L0Dynamic,
            paste0(output, "GRN_L0.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F
)
ground_truth_simulation(
  intput = paste0(output, "GRN_GENIE3.txt"),
  output = output,
  dataset_dir = simulation_data_dir,
  file = "refNetwork.csv"
)
evaluationObject <- prepareEval(paste0(output, "GRN_L0.txt"),
                                paste0(paste0(output, "ground_truth.tsv")),
                                totalPredictionsAccepted = 100000
)
AUROC_L0 <- calcAUROC(evaluationObject)
AUPRC_L0 <- calcAUPR(evaluationObject)