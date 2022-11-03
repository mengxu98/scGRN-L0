
TF <- intersect(tfs$TF, rownames(data_grn_bulk))
library(GENIE3)
weightMat <- GENIE3(
    exprMatrix = data_grn_bulk,
    targets = TF,
    nCores = 32
)
weightdf <- getLinkList(weightMat)
names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
write.table(weightdf, file = paste0(output, "GRN_GENIE3.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
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
    totalPredictionsAccepted = 100000
)
AUROC_GENIE3 <- calcAUROC(evaluationObject)
AUPRC_GENIE3 <- calcAUPR(evaluationObject)



L0Dynamic <- L0REG(
    matrix = data_grn_bulk,
    targets = TF,
    penalty = "L0"
)
write.table(L0Dynamic,
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
    totalPredictionsAccepted = 100000
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
