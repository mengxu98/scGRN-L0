

# Load all the script files, libraries and functions
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
source("framework_main.R")
source("ground-truth.R")
source("Function-L0REG.R")

uploading <- dget("SINCERITIES functions/uploading.R")
SINCERITITES <- dget("SINCERITIES functions/SINCERITIES.R")
SINCERITITES_L0 <- dget("SINCERITIES functions/SINCERITIES_L0.R")
auc_from_ranks_TC_sign <- dget("SINCERITIES functions/auc_from_ranks_TC_sign.R")
final_ranked_predictions <- dget("SINCERITIES functions/final_ranked_predictions.R")
data_path <- c(
    "hESC"
)
output <- "../scGRN-L0_output/output_scRNA-Seq/"
hnetwork_data_dir <- "../scGRN-L0_data/BEELINE-Networks/Networks/human/"
mnetwork_data_dir <- "../scGRN-L0_data/BEELINE-Networks/Networks/mouse/"
tfs <- read.csv("../scGRN-L0_data/BEELINE-Networks/human-tfs.csv")
j=1
simulation_data_dir <- paste0("../scGRN-L0_data/BEELINE-data/inputs/scRNA-Seq/", data_path[j], "/")
simulation_data_file <- "ExpressionData.csv"
simulation_PseudoTime_file <- "PseudoTime.csv"
simulation_data <- read.csv(paste0(simulation_data_dir, simulation_data_file), row.names = 1)
simulation_data <- simulation_data[tfs$TF[100:200], ] %>% na.omit()
head(simulation_data[1:3, 1:3])
simulation_PseudoTime <- read.csv(paste0(simulation_data_dir, simulation_PseudoTime_file), row.names = 1)
if (ncol(simulation_PseudoTime) > 1) {
    simulation_data_news <- c()
    for (p in 1:ncol(simulation_PseudoTime)) {
        simulation_data_new <- cbind.data.frame(t(simulation_data), h = simulation_PseudoTime[, p]) %>% na.omit()
        write.csv(simulation_data_new, paste0(simulation_data_dir, "ExpressionData_", p, ".csv"), row.names = FALSE)
        simulation_data_news <- rbind.data.frame(simulation_data_news, simulation_data_new)
    }
    write.csv(simulation_data_news, paste0(simulation_data_dir, "ExpressionData_all", ".csv"), row.names = FALSE)
} else {
    simulation_data_new <- cbind.data.frame(t(simulation_data), h = simulation_PseudoTime$PseudoTime) %>% na.omit()
    write.csv(simulation_data_new, paste0(simulation_data_dir, "ExpressionData_all", ".csv"), row.names = FALSE)
}
data_grn_bulk <- read.csv("../scGRN-L0_data/BEELINE-data/inputs/scRNAseq_preprocessing/data/GSE75748/GSE75748_bulk_time_course_ec.csv",
    header = T, row.names = 1,
) %>% as.matrix()
TF <- intersect(tfs$TF, rownames(data_grn_bulk))
library(GENIE3)
weightMat <- GENIE3(
    exprMatrix = data_grn_bulk,
    targets = TF,
    nCores = 32
)
weightdf <- getLinkList(weightMat)
names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
write.table(weightdf[1:10000,], file = paste0(output, "GRN_GENIE3.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
if (data_path[j] %in% c("hESC", "hHep")) {
    ground_truth_h(
        intput = paste0(output, "GRN_GENIE3.txt"),
        output = output,
        dataset_dir = hnetwork_data_dir,
        database = "STRING"
    )
} else {
    ground_truth_m(
        intput = paste0(output, "GRN_GENIE3.txt"),
        output = output,
        dataset_dir = mnetwork_data_dir,
        database = "STRING"
    )
}
evaluationObject <- prepareEval(paste0(output, "GRN_GENIE3.txt"),
    paste0(paste0(output, "ground_truth.tsv")),
    totalPredictionsAccepted = 10000
)
AUROC_GENIE3 <- calcAUROC(evaluationObject)
AUPRC_GENIE3 <- calcAUPR(evaluationObject)

L0Dynamic <- L0REG(
    matrix = data_grn_bulk,
    targets = TF,
    penalty = "L0"
)
write.table(L0Dynamic[1:10000,],
    paste0(output, "GRN_L0.txt"),
    sep = "\t",
    quote = F,
    row.names = F,
    col.names = F
)
if (data_path[j] %in% c("hESC", "hHep")) {
    ground_truth_h(
        intput = paste0(output, "GRN_L0.txt"),
        output = output,
        dataset_dir = hnetwork_data_dir,
        database = "STRING"
    )
} else {
    ground_truth_m(
        intput = paste0(output, "GRN_L0.txt"),
        output = output,
        dataset_dir = mnetwork_data_dir,
        database = "STRING"
    )
}
evaluationObject <- prepareEval(paste0(output, "GRN_L0.txt"),
    paste0(paste0(output, "ground_truth.tsv")),
    totalPredictionsAccepted = 10000
)
AUROC_L0 <- calcAUROC(evaluationObject)
AUPRC_L0 <- calcAUPR(evaluationObject)
# --------------------------------------------------
L0DynamicL2 <- L0REG(t(data_grn),
    regulators = colnames(data_grn),
    targets = colnames(data_grn),
    # maxSuppSize = 5,
    penalty = "L0L2"
)
write.table(L0DynamicL2,
    paste0(output, "GRN_L0L2.txt"),
    sep = "\t",
    quote = F,
    row.names = F,
    col.names = F
)
evaluationObject <- prepareEval(paste0(output, "GRN_L0L2.txt"),
    paste0(paste0(output, "ground_truth.tsv")),
    totalPredictionsAccepted = 100000
)
AUROC_L0L2 <- calcAUROC(evaluationObject)
AUPRC_L0L2 <- calcAUPR(evaluationObject)
