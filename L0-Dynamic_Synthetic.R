

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

uploading <- dget("SINCERITIES functions/uploading.R")
SINCERITITES <- dget("SINCERITIES functions/SINCERITIES.R")
SINCERITITES_L0 <- dget("SINCERITIES functions/SINCERITIES_L0.R")
auc_from_ranks_TC_sign <- dget("SINCERITIES functions/auc_from_ranks_TC_sign.R")
final_ranked_predictions <- dget("SINCERITIES functions/final_ranked_predictions.R")

data_path <- c(
    "dyn-BF",
    "dyn-BFC",
    "dyn-CY",
    "dyn-LI",
    "dyn-LL",
    "dyn-TF"
)
cell_num <- c(
    "100",
    "200",
    "500",
    "2000",
    "5000"
)
cell_drop <- c(
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10"
)

output <- "output_Synthetic/"

evaluation_infromations_all <- c()
for (j in 1:length(data_path)) {
    evaluation_infromations2 <- c()
    for (k in 1:length(cell_num)) {
        evaluation_infromations <- c()
        for (i in 1:length(cell_drop)) {
            # if (i == 1) {
            #   message(paste0("----- Now run the GRN model with ", data_path[j], " dataset and drop-out of cell with 1", "! -----"))
            #   simulation_data_dir <- paste0("../scGRN-L0_data/BEELINE-data/inputs/Synthetic/", data_path[j], "/", data_path[j], "-2000-1", "/")
            # } else {
            #   message(paste0("----- Now run the GRN model with ", data_path[j], " dataset and drop-out of cell with ", cell_drop[i], "! -----"))
            simulation_data_dir <- paste0("../scGRN-L0_data/BEELINE-data/inputs/Synthetic/", data_path[j], "/", data_path[j], "-", cell_num[k], "-", cell_drop[i], "/")
            # }
            simulation_data_file <- "ExpressionData.csv"
            simulation_PseudoTime_file <- "PseudoTime.csv"

            simulation_data <- read.csv(paste0(simulation_data_dir, simulation_data_file), row.names = 1)
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

            DATA <- uploading(paste0(simulation_data_dir, "ExpressionData_all", ".csv"))
            # --------------------------------------------------
            if (F) {
                # result_L0 <- SINCERITITES_L0(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 0, penalty = "L0L2", algorithm = "CD")
                result_L0 <- SINCERITITES_L0(DATA, distance = 1, method = 1, noDIAG = 1, SIGN = 0, penalty = "L0", algorithm = "CD")

                # Final ranked list of regulatory edges
                adj_matrix_L0 <- result_L0$adj_matrix / max(result_L0$adj_matrix)
                table_L0 <- final_ranked_predictions(adj_matrix_L0, DATA$genes, SIGN = 1, saveFile = TRUE) %>% as.data.frame()
                write.table(table_L0[, -4], paste0(output, "test_L0.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
                ground_truth_simulation(
                    intput = paste0(output, "test_L0.txt"),
                    output = output,
                    dataset_dir = simulation_data_dir,
                    file = "refNetwork.csv"
                )
                evaluationObject <- prepareEval(paste0(output, "test_L0.txt"),
                    paste0(paste0(output, "ground_truth.tsv")),
                    totalPredictionsAccepted = 100000
                )
                evaluationObject <- prepareEval(paste0(output, "ground_pred.txt"),
                    paste0(paste0(output, "ground_truth.tsv")),
                    totalPredictionsAccepted = 100000
                )
                AUROC_L0_N <- calcAUROC(evaluationObject)
                AUPR_L0_N <- calcAUPR(evaluationObject)
                truth_network <- convertSortedRankTSVToAdjMatrix(paste0(output, "ground_truth.tsv"))
                AUCresult_L0 <- auc_from_ranks_TC_sign(adj_matrix_L0, truth_network, 1000)
                AUROC_L0_S <- AUCresult_L0$AUROC
                AUPR_L02_S <- AUCresult_L0$AUPR
                AUROC_L0_N
                AUROC_L0_S
            }

            if (F) {
                result_L0L2 <- SINCERITITES_L0(DATA, distance = 1, method = 1, noDIAG = 1, SIGN = 0, penalty = "L0L2", algorithm = "CD")
                # Final ranked list of regulatory edges
                adj_matrix_L0L2 <- result_L0L2$adj_matrix / max(result_L0L2$adj_matrix)
                table_L0L2 <- final_ranked_predictions(adj_matrix_L0, DATA$genes, SIGN = 1, saveFile = TRUE) %>% as.data.frame()
                write.table(table_L0L2[, -4], paste0(output, "test_L0L2.txt"), row.names = F, col.names = F, sep = "\t", quote = F)

                ground_truth_simulation(
                    intput = paste0(output, "test_L0L2.txt"),
                    output = output,
                    dataset_dir = simulation_data_dir,
                    file = "refNetwork.csv"
                )
                evaluationObject <- prepareEval(paste0(output, "test_L0L2.txt"),
                    paste0(paste0(output, "ground_truth.tsv")),
                    totalPredictionsAccepted = 100000
                )
                evaluationObject <- prepareEval(paste0(output, "ground_pred.txt"),
                    paste0(paste0(output, "ground_truth.tsv")),
                    totalPredictionsAccepted = 100000
                )
                AUROC_L0L2_N <- calcAUROC(evaluationObject)
                AUPR_L0L2_N <- calcAUPR(evaluationObject)
                truth_network <- convertSortedRankTSVToAdjMatrix(paste0(output, "ground_truth.tsv"))
                AUCresult_L0L2_S <- auc_from_ranks_TC_sign(adj_matrix_L0L2, truth_network, 1000)
                AUROC_L0L2_S <- AUCresult_L0L2_S$AUROC
                AUPR_L0L2_S <- AUCresult_L0L2_S$AUPR
                AUROC_L0L2_N
                AUROC_L0L2_S
            }

            # --------------------------------------------------
            if (T) {
                # result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 1)
                result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 1, SIGN = 0)
                # Final ranked list of regulatory edges
                adj_matrix <- result$adj_matrix / max(result$adj_matrix)
                if (adj_matrix[1, 1] == "NaN") {
                    SINCERITITES_AUROC_S <- NA
                    SINCERITITES_AUROC_N <- NA
                } else {
                    table <- final_ranked_predictions(adj_matrix, DATA$genes, SIGN = 1, saveFile = TRUE)
                    table <- as.data.frame(table)
                    write.table(table[, -4], paste0(output, "test_SINCERITITES.txt"), row.names = F, col.names = F, sep = "\t", quote = F)

                    ground_truth_simulation(
                        intput = paste0(output, "test_SINCERITITES.txt"),
                        output = output,
                        dataset_dir = simulation_data_dir,
                        file = "refNetwork.csv"
                    )

                    evaluationObject <- prepareEval(paste0(output, "test_SINCERITITES.txt"),
                        paste0(paste0(output, "ground_truth.tsv")),
                        totalPredictionsAccepted = 100000
                    )
                    evaluationObject <- prepareEval(paste0(output, "ground_pred.txt"),
                        paste0(paste0(output, "ground_truth.tsv")),
                        totalPredictionsAccepted = 100000
                    )
                    SINCERITITES_AUROC_N <- calcAUROC(evaluationObject)
                    SINCERITITES_AUPR_N <- calcAUPR(evaluationObject)

                    adj_matrix <- na.omit(adj_matrix)
                    truth_network <- convertSortedRankTSVToAdjMatrix(paste0(output, "ground_truth.tsv"))
                    AUCresult_SINCERITITES <- auc_from_ranks_TC_sign(adj_matrix, truth_network, 1000)
                    SINCERITITES_AUROC_S <- AUCresult_SINCERITITES$AUROC
                    SINCERITITES_AUPR_S <- AUCresult_SINCERITITES$AUPR
                }

                SINCERITITES_AUROC_N
                SINCERITITES_AUROC_S
            }
            # --------------------------------------------------
            data_GENIE3 <- read.csv(paste0(simulation_data_dir, "ExpressionData_1.csv"),
                header = T
            ) %>% as.matrix()
            PseudoTime <- data_GENIE3[, ncol(data_GENIE3)]
            data_GENIE3 <- data_GENIE3[, -ncol(data_GENIE3)]
            if (T) {
                # NIMEFI(data_GENIE3,
                #     GENIE = T, SVM = F, EL = F, penalty = "L0",
                #     outputFileName = paste0(output, "output_GENIE3"),
                #     outputFileFormat = "txt",
                #     SVMRankThreshold = 5, SVMEnsembleSize = 100,
                #     ELPredSampleMin = 20, ELPredSampleMax = 80,
                #     ELExpSampleMin = 20, ELExpSampleMax = 80,
                #     ELRankThreshold = 5, ELEnsembleSize = 500
                # )

                # evaluationObject <- prepareEval(paste0(output, "output_GENIE3.txt"),
                #     paste0(paste0(output, "ground_truth.tsv")),
                #     totalPredictionsAccepted = 100000
                # )

                # GENIE3_AUROC_NIMEFI <- calcAUROC(evaluationObject)
                # GENIE3_AUPR_NIMEFI <- calcAUPR(evaluationObject)

                library(GENIE3)
                weightMat <- GENIE3(
                    exprMatrix = t(data_GENIE3),
                    nCores = 32
                )

                AUCresult_GENIE3 <- auc_from_ranks_TC_sign(weightMat, truth_network, 1000)
                GENIE3_AUROC_S <- AUCresult_GENIE3$AUROC
                GENIE3_AUPR_S <- AUCresult_GENIE3$AUPR

                weightdf <- getLinkList(weightMat)
                # weightdf <- read.table("output_NIMEFI_L0.txt", header = F)
                names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
                write.table(weightdf, file = paste0(output, "output_GENIE32.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
                evaluationObject <- prepareEval(paste0(output, "output_GENIE32.txt"),
                    paste0(paste0(output, "ground_truth.tsv")),
                    totalPredictionsAccepted = 100000
                )
                GENIE3_AUROC <- calcAUROC(evaluationObject)
                GENIE3_AUPR <- calcAUPR(evaluationObject)
            }

            if (F) {
                NIMEFI(data_GENIE3,
                    GENIE = F, SVM = F, EL = T, penalty = "L0",
                    outputFileName = paste0(output, "output_NIMEFI_L0"),
                    outputFileFormat = "txt",
                    SVMRankThreshold = 5, SVMEnsembleSize = 100,
                    ELPredSampleMin = 20, ELPredSampleMax = 80,
                    ELExpSampleMin = 20, ELExpSampleMax = 80,
                    ELRankThreshold = 5, ELEnsembleSize = dim(data_GENIE3)[1]
                )

                evaluationObject <- prepareEval(paste0(output, "output_NIMEFI_L0.txt"),
                    paste0(output, "ground_truth.tsv"),
                    totalPredictionsAccepted = 100000
                )

                AUROC_NIMEFI_L0 <- calcAUROC(evaluationObject)
                AUPR_NIMEFI_L0 <- calcAUPR(evaluationObject)

                # --------------------------------------------------
                NIMEFI(data_GENIE3,
                    GENIE = F, SVM = F, EL = T, penalty = "L0L2",
                    outputFileName = paste0(output, "output_NIMEFI_L0"),
                    outputFileFormat = "txt",
                    SVMRankThreshold = 5, SVMEnsembleSize = 100,
                    ELPredSampleMin = 20, ELPredSampleMax = 80,
                    ELExpSampleMin = 20, ELExpSampleMax = 80,
                    ELRankThreshold = 5, ELEnsembleSize = dim(data_GENIE3)[1]
                )

                evaluationObject <- prepareEval(paste0(output, "output_NIMEFI_L0.txt"),
                    paste0(paste0(output, "ground_truth.tsv")),
                    totalPredictionsAccepted = 100000
                )

                AUROC_NIMEFI_L0L2 <- calcAUROC(evaluationObject)
                AUPR_NIMEFI_L0L2 <- calcAUPR(evaluationObject)

                AUROC_NIMEFI_L0L2
                AUROC_NIMEFI_L0
            }

            if (T) {
                source("DynamicGRNPipe_3.constructionNetwork_L0.R")
                L0REG_L0 <- L0REG(t(data_GENIE3),
                    regulators = colnames(data_GENIE3),
                    targets = colnames(data_GENIE3), penalty = "L0"
                )
                write.table(L0REG_L0,
                    paste0(output, "output_L0GRN.txt"),
                    sep = "\t",
                    quote = F,
                    row.names = F,
                    col.names = F
                )
                evaluationObject <- prepareEval(paste0(output, "output_L0GRN.txt"),
                    paste0(paste0(output, "ground_truth.tsv")),
                    totalPredictionsAccepted = 100000
                )

                AUROC_L0REG_L0_N <- calcAUROC(evaluationObject)
                AUPR_L0REG_L0_N <- calcAUPR(evaluationObject)
                AUROC_L0REG_L0_N

                # L0REG_L0 <- L0REG(t(data_GENIE3),
                #   regulators = colnames(data_GENIE3),
                #   targets = colnames(data_GENIE3), penalty = "L0"
                # )
                # AUCresult_L0REG <- auc_from_ranks_TC_sign(L0REG_L0, truth_network, 1000)
                # AUROC_L0REG_L0_S <- AUCresult_L0REG$AUROC
                # AUROC_L0REG_L0_S
                # --------------------------------------------------
                L0REG_L0L2 <- L0REG(t(data_GENIE3),
                    regulators = colnames(data_GENIE3),
                    targets = colnames(data_GENIE3), penalty = "L0L2"
                )
                write.table(L0REG_L0L2,
                    paste0(output, "output_L0GRN.txt"),
                    sep = "\t",
                    quote = F,
                    row.names = F,
                    col.names = F
                )
                evaluationObject <- prepareEval(paste0(output, "output_L0GRN.txt"),
                    paste0(paste0(output, "ground_truth.tsv")),
                    totalPredictionsAccepted = 100000
                )

                AUROC_L0REG_L0L2_N <- calcAUROC(evaluationObject)
                AUPR_L0REG_L0L2_N <- calcAUPR(evaluationObject)
                AUROC_L0REG_L0L2_N

                # L0REG_L0L2 <- L0REG(t(data_GENIE3),
                #   # regulators = colnames(data_GENIE3),
                #   targets = colnames(data_GENIE3), penalty = "L0L2"
                # )
                # AUCresult_L0REG <- auc_from_ranks_TC_sign(L0REG_L0L2, truth_network, 1000)
                # AUROC_L0REG_L0L2_S <- AUCresult_L0REG$AUROC
                # AUROC_L0REG_L0L2_S
            }

            if (T) {
                data_GENIE3 <- read.csv(paste0(simulation_data_dir, "ExpressionData_1.csv"),
                    header = T
                )
                data_GENIE3 <- data_GENIE3[order(data_GENIE3$h), ]
                PseudoTime <- data_GENIE3[, ncol(data_GENIE3)]
                # PseudoTime <- PseudoTime[order(PseudoTime), ] %>% as.data.frame()
                # rownames(PseudoTime) <- PseudoTime$X
                data_GENIE3 <- data_GENIE3[, -ncol(data_GENIE3)] %>% as.matrix()
                if (F) {
                    L0REG_L0_adjs <- matrix(0, ncol(data_GENIE3), ncol(data_GENIE3))
                    rownames(L0REG_L0_adjs) <- colnames(data_GENIE3)
                    colnames(L0REG_L0_adjs) <- colnames(data_GENIE3)
                    for (t in 1:5) {
                        s <- floor(nrow(data_GENIE3) * (t - 1) / 10) + 1
                        e <- floor(nrow(data_GENIE3) * t / 10)
                        data <- data_GENIE3[s:e, ]
                        L0REG_L0_1 <- L0REG(t(data),
                            regulators = colnames(data),
                            targets = colnames(data), penalty = "L0"
                        )
                        L0REG_L0_1$weight <- as.numeric(L0REG_L0_1$weight)
                        L0REG_L0_1 <- as.matrix(L0REG_L0_1)
                        L0REG_L0_adj <- matrix(0, ncol(data_GENIE3), ncol(data_GENIE3))
                        rownames(L0REG_L0_adj) <- colnames(data_GENIE3)
                        colnames(L0REG_L0_adj) <- colnames(data_GENIE3)
                        L0REG_L0_adj[L0REG_L0_1[, 1:2]] <- L0REG_L0_1[, 3]
                        L0REG_L0_adjs <- L0REG_L0_adjs + as.numeric(L0REG_L0_adj)
                    }
                    L0REG_L0_adjs <- L0REG_L0_adjs / max(L0REG_L0_adjs)

                    AUCresult_L0REG_L0 <- auc_from_ranks_TC_sign(L0REG_L0_adjs, truth_network, 1000)
                    L0REG_L0Dynamic_AUROC_S <- AUCresult_L0REG_L0$AUROC
                    L0REG_L0Dynamic_AUPR_S <- AUCresult_L0REG_L0$AUPR

                    weightdf <- getLinkList(L0REG_L0_adjs)
                    # weightdf <- read.table("output_NIMEFI_L0.txt", header = F)
                    names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
                    write.table(weightdf, file = paste0(output, "L0REG_L02.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
                    evaluationObject <- prepareEval(paste0(output, "L0REG_L02.txt"),
                        paste0(paste0(output, "ground_truth.tsv")),
                        totalPredictionsAccepted = 100000
                    )
                    L0REG_L0Dynamic_AUROC <- calcAUROC(evaluationObject)
                    L0REG_L0Dynamic_AUPR <- calcAUPR(evaluationObject)
                    print(
                        paste0(
                            L0REG_L0Dynamic_AUROC,
                            "---",
                            L0REG_L0Dynamic_AUROC_S
                        )
                    )
                }

                if (T) {
                    L0REG_L0_adjs <- matrix(0, ncol(data_GENIE3), ncol(data_GENIE3))
                    rownames(L0REG_L0_adjs) <- colnames(data_GENIE3)
                    colnames(L0REG_L0_adjs) <- colnames(data_GENIE3)
                    win <- 50
                    for (t in 1:(nrow(data_GENIE3) - win)) {
                        s <- t
                        e <- win + t
                        # s <- floor(nrow(data_GENIE3) * (t - 1) / 10) + 1
                        # e <- floor(nrow(data_GENIE3) * t / 10)
                        if (e < dim(data_GENIE3)[1]) {
                            data <- data_GENIE3[s:e, ]
                        } else {
                            data <- data_GENIE3
                        }
                        L0REG_L0_1 <- L0REG(t(data),
                            regulators = colnames(data),
                            targets = colnames(data), penalty = "L0"
                        )
                        L0REG_L0_1$weight <- as.numeric(L0REG_L0_1$weight)
                        L0REG_L0_1 <- as.matrix(L0REG_L0_1)
                        L0REG_L0_adj <- matrix(0, ncol(data_GENIE3), ncol(data_GENIE3))
                        rownames(L0REG_L0_adj) <- colnames(data_GENIE3)
                        colnames(L0REG_L0_adj) <- colnames(data_GENIE3)
                        L0REG_L0_adj[L0REG_L0_1[, 1:2]] <- L0REG_L0_1[, 3]
                        L0REG_L0_adjs <- L0REG_L0_adjs + as.numeric(L0REG_L0_adj)
                    }
                    L0REG_L0_adjs <- L0REG_L0_adjs / max(L0REG_L0_adjs)

                    AUCresult_L0REG_L0 <- auc_from_ranks_TC_sign(L0REG_L0_adjs, truth_network, 1000)
                    L0REG_L0Dynamic_AUROC_S <- AUCresult_L0REG_L0$AUROC
                    L0REG_L0Dynamic_AUPR_S <- AUCresult_L0REG_L0$AUPR

                    weightdf <- getLinkList(L0REG_L0_adjs)
                    # weightdf <- read.table("output_NIMEFI_L0.txt", header = F)
                    names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
                    write.table(weightdf, file = paste0(output, "L0REG_L02.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
                    evaluationObject <- prepareEval(paste0(output, "L0REG_L02.txt"),
                        paste0(paste0(output, "ground_truth.tsv")),
                        totalPredictionsAccepted = 100000
                    )
                    L0REG_L0Dynamic_AUROC <- calcAUROC(evaluationObject)
                    L0REG_L0Dynamic_AUPR <- calcAUPR(evaluationObject)
                    print(
                        paste0(
                            L0REG_L0Dynamic_AUROC,
                            "---",
                            L0REG_L0Dynamic_AUROC_S
                        )
                    )
                }
            }

            # --------------------------------------------------
            evaluation_infromation <- data.frame(
                datasets = paste0(data_path[j], "-2000-", cell_drop[i]),
                # AUROC_SINCERITITES_L0 = AUROC_L0_N,
                # AUROC_SINCERITITES_L0_S = AUROC_L0_S,
                # AUROC_SINCERITITES_L0L2 = AUROC_L0L2_N,
                # AUROC_SINCERITITES_L0L2_S = AUROC_L0L2_S,
                AUROC_L0Dynamic = L0REG_L0Dynamic_AUROC,
                # AUROC_L0Dynamic_S = L0REG_L0Dynamic_AUROC_S,
                AUROC_L0REG_L0 = AUROC_L0REG_L0_N,
                # AUROC_L0REG_L0_S = AUROC_L0REG_L0_S,
                AUROC_L0REG_L0L2 = AUROC_L0REG_L0L2_N,
                # AUROC_L0REG_L0L2_S = AUROC_L0REG_L0L2_S,
                # AUROC_NIMEFI_L0L2 = AUROC_NIMEFI_L0L2,
                # AUROC_NIMEFI_L0 = AUROC_NIMEFI_L0,
                # AUROC_GENIE3_NIMEFI = GENIE3_AUROC_NIMEFI,
                AUROC_GENIE3 = GENIE3_AUROC,
                # AUROC_GENIE3_S = GENIE3_AUROC_S,
                # AUROC_SINCERITITES_S = SINCERITITES_AUROC_S,
                AUROC_SINCERITITES = SINCERITITES_AUROC_N
            )

            message(paste0("----- ", evaluation_infromation, " -----"))
            print(evaluation_infromation)
            evaluation_infromations <- rbind.data.frame(evaluation_infromations, evaluation_infromation)
        }
        evaluation_infromations2 <- rbind.data.frame(evaluation_infromations2, evaluation_infromations)
    }
    evaluation_infromations_all <- rbind.data.frame(evaluation_infromations_all, evaluation_infromations2)
}
# evaluation_infromations_all
# evaluation_infromations_all[evaluation_infromations_all==0] <- NA
evaluation_infromations_all <- na.omit(evaluation_infromations_all)
write.csv(evaluation_infromations_all, "Results/evaluation_infromations.csv")

if (F) {
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(ggpubr)
    library(dplyr)
    library(tidyr)
    library(pheatmap)
    library(tidyverse)
    library(RColorBrewer)

    data_path <- c(
        "dyn-BF-",
        "dyn-BFC",
        "dyn-CY",
        "dyn-LI",
        "dyn-LL",
        "dyn-TF"
    )
    mean(evaluation_infromations_all$AUROC_SINCERITITES)
    evaluation_infromations_all <- read.csv("Results/evaluation_infromations.csv")
    head(evaluation_infromations_all[1:3, 1:3])
    for (i in 1:length(data_path)) {
        dataset <- data_path[i]
        evaluation_infromations_GSD <- evaluation_infromations_all[grep(dataset, evaluation_infromations_all$datasets), ]
        evaluation_infromations_GSD <- evaluation_infromations_GSD[, c("datasets", "AUROC_NIMEFI_L0", "AUROC_GENIE3", "AUROC_SINCERITITES_N")]
        names(evaluation_infromations_GSD) <- c("Dataset", "L0-Dynamic", "GENIE3", "SINCERITITES")
        methods_barplot_all <- evaluation_infromations_GSD %>%
            as.data.frame() %>%
            pivot_longer(
                cols = 2:c(ncol(evaluation_infromations_GSD)),
                names_to = "Methods",
                values_to = "AUROC"
            )

        my_comparisons <- list(
            c("L0-Dynamic", "GENIE3"),
            c("SINCERITITES", "GENIE3"),
            c("SINCERITITES", "L0-Dynamic")
        )

        p <- ggplot(
            methods_barplot_all,
            aes(
                x = Methods,
                y = AUROC
            )
        ) +
            # guides(fill = guide_legend(title = NULL)) +
            stat_compare_means(
                method = "wilcox.test",
                label = "p.signif",
                comparisons = my_comparisons,
                bracket.size = 0.6,
                sizen = 4,
                color = "#6699cc"
            ) +
            labs(
                x = "Methods",
                y = "AUROC"
            ) +
            # stat_summary(
            #   fun.data = "mean_sdl",
            #   fun.args = list(mult = 1),
            #   geom = "pointrange",
            #   color = "gray"
            # ) +
            geom_violin(aes(fill = Methods),
                trim = FALSE
            ) +
            geom_boxplot(width = 0.2) +
            theme(legend.position = "bottom") +
            # facet_wrap(~Methods) +
            # theme_gray() +
            theme_bw() +
            theme(
                axis.text.x = element_text(
                    angle = 45,
                    hjust = 1,
                    vjust = 1,
                    size = 10
                )
            )
        p
        ggsave(paste0("Results/Methods-contrast-", dataset, "-1.png"), width = 3, height = 4, dpi = 600)

        # methods_barplot_all %>% ggplot(., aes(x = Methods, y = AUROC, colour = Methods)) +
        #   geom_boxplot() +
        #   theme_bw() +
        #   theme(legend.position = "none")

        P1 <- ggplot(methods_barplot_all, aes(x = Methods, y = AUROC, fill = Methods)) + # ”fill=“设置填充颜色
            stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) + # 由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
            geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white") + # size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属性
            geom_jitter(aes(fill = Methods), width = 0.2, shape = 21, size = 2.5) + # 设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
            scale_fill_manual(values = c("black", "gray", "white", "#c6524a", "#eabf00", "#696969")) + # 设置填充的颜色
            scale_color_manual(values = c("black", "#2874c5", "#008a00", "#c6524a", "#eabf00", "#696969")) + # 设置散点图的圆圈的颜色为黑色
            # scale_x_discrete(labels = c("L0-Dynamic", "GENIE3","SINCERITITES"))+
            ggtitle(" ") + # 设置总的标题
            theme_bw() +
            theme(legend.position = "bottom") +
            # theme(legend.position="none",
            #       axis.text.x=element_text(colour="black",family="Times",size=14), #设置x轴刻度标签的字体属性
            #       axis.text.y=element_text(family="Times",size=14,face="plain"), #设置x轴刻度标签的字体属性
            #       axis.title.y=element_text(family="Times",size = 14,face="plain"), #设置y轴的标题的字体属性
            #       axis.title.x=element_text(family="Times",size = 14,face="plain"), #设置x轴的标题的字体属性
            #       plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
            #       # panel.grid.major = element_blank(), #不显示网格线
            #       panel.grid.minor = element_blank())+
            ylab("AUROC") +
            xlab("Method")
        # P1
        # ggsave("Results/Methods-contrast-1.png", width = 4, height = 4, dpi = 600)

        mycol <- c("black", "gray", "white")
        results_10nets <- evaluation_infromations_GSD
        results_10nets$Dataset <- c(
            "100-1", "100-2", "100-3", "100-4", "100-5", "100-6", "100-7", "100-8", "100-9", "100-10",
            "200-1", "200-2", "200-3", "200-4", "200-5", "200-6", "200-7", "200-8", "200-9", "200-10",
            "300-1", "300-2", "300-3", "300-4", "300-5", "300-6", "300-7", "300-8", "300-9", "300-10",
            "400-1", "400-2", "400-3", "400-4", "400-5", "400-6", "400-7", "400-8", "400-9", "400-10",
            "500-1", "500-2", "500-3", "500-4", "500-5", "500-6", "500-7", "500-8", "500-9", "500-10"
        )
        df_res10 <- melt(results_10nets, id = "Dataset", variable.name = "Method", value.name = "AUROC")
        df_res10$Method <- factor(df_res10$Method,
            levels = c("L0-Dynamic", "GENIE3", "SINCERITITES"),
            labels = c("L0-Dynamic", "GENIE3", "SINCERITITES")
        )

        p2 <- ggplot(df_res10, aes(x = Dataset, y = AUROC, fill = Method)) +
            geom_bar(stat = "identity", position = position_dodge(), color = "black", width = .6) +
            scale_fill_manual(values = mycol) +
            # geom_errorbar(aes(ymin=AUROC - Sd, ymax=AUROC + Sd), position = position_dodge(.6), width=.2)
            scale_x_discrete(labels = results_10nets$Dataset) +
            theme_bw() +
            theme(
                axis.text.x = element_text(
                    angle = 45,
                    hjust = 1,
                    vjust = 1,
                    size = 10
                )
            )
        p2
        ggsave(paste0("Results/Methods-contrast-", dataset, "-2.png"), width = 8, height = 4, dpi = 600)
    }
}
