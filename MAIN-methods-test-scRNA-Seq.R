

# 将单细胞时间数据进行分片，比如0.2，0.5归为0-1之间
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
  "GSD",
  # "HSC",
  # "mCAD",
  "VSC"
)
cell_drop <- c(
  "1", "1-50", "1-70",
  # "2", "2-50", "2-70",
  # "3", "3-50", "3-70",
  # "4", "4-50", "4-70",
  # "5", "5-50", "5-70",
  # "6", "6-50", "6-70",
  # "7", "7-50", "7-70",
  # "8", "8-50", "8-70",
  # "9", "9-50", "9-70",
  "10", "10-50", "10-70"
)
evaluation_infromations_all <- c()
for (j in 1:length(data_path)) {
  evaluation_infromations <- c()
  for (i in 1:length(cell_drop)) {
    if (i == 1) {
      message(paste0("----- Now run the GRN model with ", data_path[j], " dataset and drop-out of cell with 1", "! -----"))
      simulation_data_dir <- paste0("../scGRN-L0_data/BEELINE-data/inputs/Curated/", data_path[j], "/", data_path[j], "-2000-1", "/")
    } else {
      message(paste0("----- Now run the GRN model with ", data_path[j], " dataset and drop-out of cell with ", cell_drop[i], "! -----"))
      simulation_data_dir <- paste0("../scGRN-L0_data/BEELINE-data/inputs/Curated/", data_path[j], "/", data_path[j], "-2000-", cell_drop[i], "/")
    }
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

    DATA <- uploading(paste0(simulation_data_dir, "ExpressionData_all", ".csv"))
    # result_L0 <- SINCERITITES_L0(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 0, penalty = "L0L2", algorithm = "CD")
    result_L0 <- SINCERITITES_L0(DATA, distance = 1, method = 1, noDIAG = 1, SIGN = 0, penalty = "L0L2", algorithm = "CD")

    # Final ranked list of regulatory edges
    adj_matrix_L0 <- result_L0$adj_matrix / max(result_L0$adj_matrix)
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
    evaluationObject <- prepareEval("output/ground_pred.txt",
      paste0("output/ground_truth.tsv"),
      totalPredictionsAccepted = 100000
    )
    L0_AUROC <- calcAUROC(evaluationObject)
    L0_AUPR <- calcAUPR(evaluationObject)

    truth_network <- convertSortedRankTSVToAdjMatrix("output/ground_truth.tsv")

    AUCresult_L0 <- auc_from_ranks_TC_sign(adj_matrix_L0, truth_network, 1000)
    L0_AUROC2 <- AUCresult_L0$AUROC
    L0_AUPR2 <- AUCresult_L0$AUPR
    L0_AUROC
    L0_AUROC2

    # --------------------------------------------------
    # result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 1)
    result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 1, SIGN = 0)
    # Final ranked list of regulatory edges
    adj_matrix <- result$adj_matrix / max(result$adj_matrix)
    table <- final_ranked_predictions(adj_matrix, DATA$genes, SIGN = 1, saveFile = TRUE)
    table <- as.data.frame(table)
    write.table(table[, -4], "output/test.txt", row.names = F, col.names = F, sep = "\t", quote = F)

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

    SINCERITITES_AUROC <- calcAUROC(evaluationObject)
    SINCERITITES_AUPR <- calcAUPR(evaluationObject)

    AUCresult_SINCERITITES <- auc_from_ranks_TC_sign(adj_matrix, truth_network, 1000)
    SINCERITITES_AUROC2 <- AUCresult_SINCERITITES$AUROC
    SINCERITITES2_AUPR <- AUCresult_SINCERITITES$AUPR
    SINCERITITES_AUROC
    SINCERITITES_AUROC2
    # --------------------------------------------------
    if (T) {
      data_GENIE3 <- read.csv(paste0(simulation_data_dir, "ExpressionData_all.csv"),
        header = T
      ) %>% as.matrix()

      PseudoTime <- data_GENIE3[, ncol(data_GENIE3)]
      # PseudoTime <- PseudoTime[order(PseudoTime), ] %>% as.data.frame()
      # rownames(PseudoTime) <- PseudoTime$X

      data_GENIE3 <- data_GENIE3[, -ncol(data_GENIE3)]

      NIMEFI(data_GENIE3,
        GENIE = T, SVM = F, EL = F, penalty = "L0",
        outputFileName = "output/output_GENIE3",
        outputFileFormat = "txt",
        SVMRankThreshold = 5, SVMEnsembleSize = 100,
        ELPredSampleMin = 20, ELPredSampleMax = 80,
        ELExpSampleMin = 20, ELExpSampleMax = 80,
        ELRankThreshold = 5, ELEnsembleSize = 500
      )

      evaluationObject <- prepareEval("output/output_GENIE3.txt",
        paste0("output/ground_truth.tsv"),
        totalPredictionsAccepted = 100000
      )

      GENIE3_AUROC <- calcAUROC(evaluationObject)
      GENIE3_AUPR <- calcAUPR(evaluationObject)
      GENIE3_AUROC
      GENIE3_AUPR
    }
    if (T) {
      data_GENIE3 <- read.csv(paste0(simulation_data_dir, "ExpressionData_all.csv"),
        header = T
      ) %>% as.matrix()

      PseudoTime <- data_GENIE3[, ncol(data_GENIE3)]
      # PseudoTime <- PseudoTime[order(PseudoTime), ] %>% as.data.frame()
      # rownames(PseudoTime) <- PseudoTime$X

      data_GENIE3 <- data_GENIE3[, -ncol(data_GENIE3)]

      NIMEFI(data_GENIE3,
        GENIE = F, SVM = F, EL = T, penalty = "L0",
        outputFileName = "output/output_GENIE3",
        outputFileFormat = "txt",
        SVMRankThreshold = 5, SVMEnsembleSize = 100,
        ELPredSampleMin = 20, ELPredSampleMax = 80,
        ELExpSampleMin = 20, ELExpSampleMax = 80,
        ELRankThreshold = 5, ELEnsembleSize = 500
      )

      evaluationObject <- prepareEval("output/output_GENIE3.txt",
        paste0("output/ground_truth.tsv"),
        totalPredictionsAccepted = 100000
      )

      L0_AUROC3_L0 <- calcAUROC(evaluationObject)
      L0_AUPR3_L0 <- calcAUPR(evaluationObject)
      L0_AUROC3_L0
      L0_AUPR3_L0


      NIMEFI(data_GENIE3,
        GENIE = F, SVM = F, EL = T, penalty = "L0L2",
        outputFileName = "output/output_GENIE3",
        outputFileFormat = "txt",
        SVMRankThreshold = 5, SVMEnsembleSize = 100,
        ELPredSampleMin = 20, ELPredSampleMax = 80,
        ELExpSampleMin = 20, ELExpSampleMax = 80,
        ELRankThreshold = 5, ELEnsembleSize = 500
      )

      evaluationObject <- prepareEval("output/output_GENIE3.txt",
        paste0("output/ground_truth.tsv"),
        totalPredictionsAccepted = 100000
      )

      L0_AUROC3_L0L2 <- calcAUROC(evaluationObject)
      L0_AUPR3_L0L2 <- calcAUPR(evaluationObject)
      L0_AUROC3_L0L2
      L0_AUPR3_L0L2
    }
    # --------------------------------------------------
    evaluation_infromation <- data.frame(
      datasets = paste0(data_path[j], "_", cell_drop[i]),
      L0_AUROC = L0_AUROC2, # L0_AUROC2 = L0_AUROC2,
      L0_AUROC3 = L0_AUROC3_L0,
      L0_AUROC3 = L0_AUROC3_L0L2,
      SINCERITITES_AUROC = SINCERITITES_AUROC2, # SINCERITITES_AUROC2 = SINCERITITES_AUROC2,
      GENIE3_AUROC = GENIE3_AUROC
    )
    message(paste0("----- ", evaluation_infromation, " -----"))
    evaluation_infromations <- rbind.data.frame(evaluation_infromations, evaluation_infromation)
  }
  evaluation_infromations_all <- rbind.data.frame(evaluation_infromations_all, evaluation_infromations)
}

if (F) {
  mycol <- c("black", "gray", "white")
  data_path <- c(
    "GSD",
    # "HSC",
    # "mCAD",
    "VSC"
  )
  cell_drop <- c(
    "1", "1-50", "1-70",
    # "2", "2-50", "2-70",
    # "3", "3-50", "3-70",
    # "4", "4-50", "4-70",
    # "5", "5-50", "5-70",
    # "6", "6-50", "6-70",
    # "7", "7-50", "7-70",
    # "8", "8-50", "8-70",
    # "9", "9-50", "9-70",
    "10", "10-50", "10-70"
  )
  results_10nets <- evaluation_infromations_all[, c("datasets", "L0_AUROC3", "SINCERITITES_AUROC", "GENIE3_AUROC")]
  results_10nets$datasets <- c("GSD_1", "GSD_1-50", "GSD_1-70", "GSD_10", "GSD_10-50", "GSD_10-70", "VSC_1", "VSC_1-50", "VSC_1-70", "VSC_10", "VSC_10-50", "VSC_10-70")
  names(results_10nets) <- c("Dataset", "L0-Dynamic", "SINCERITITES", "GENIE3-Dynamic")
  df_res10 <- melt(results_10nets, id = "Dataset", variable.name = "Method", value.name = "AUROC")
  df_res10$Method <- factor(df_res10$Method,
    levels = c("L0-Dynamic", "SINCERITITES", "GENIE3-Dynamic"),
    labels = c("L0-Dynamic", "SINCERITITES", "GENIE3-Dynamic")
  )

  p2 <- ggplot(df_res10, aes(x = Dataset, y = AUROC, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black", width = .6) +
    scale_fill_manual(values = mycol) +
    # geom_errorbar(aes(ymin=AUROC - Sd, ymax=AUROC + Sd), position = position_dodge(.6), width=.2)
    scale_x_discrete(labels = c("GSD_1", "GSD_1-50", "GSD_1-70", "GSD_10", "GSD_10-50", "GSD_10-70", "VSC_1", "VSC_1-50", "VSC_1-70", "VSC_10", "VSC_10-50", "VSC_10-70")) +
    theme_bw()
  p2
  ggsave("Results/Methods-contrast.png", width = 12, height = 5, dpi = 600)
  # results_5nets <- read.csv(paste0("evaluation_gnw_5_", 2, ".csv"))
  # results_5nets <- results_5nets[, -1]

  # df_res5 <- melt(results_5nets, id = "Dataset", variable.name = "Method", value.name = "AUROC")
  # df_res5$Method <- factor(df_res5$Method,
  #   levels = c("L0Reg.framework", "GENIE3"),
  #   labels = c("L0Reg framework", "GENIE3")
  # )

  # p1 <- ggplot(df_res5, aes(x = Dataset, y = AUROC, fill = Method)) +
  #   geom_bar(stat = "identity", position = position_dodge(), color = "black", width = .6) +
  #   scale_fill_manual(values = mycol) +
  #   # geom_errorbar(aes(ymin=AUROC - Sd, ymax=AUROC + Sd), position = position_dodge(.6), width=.2)
  #   theme_bw() +
  #   ylim(0, 1)

  # p1 + p2 + plot_layout(widths = c(1, 2)) +
  #   plot_annotation(tag_levels = "a") +
  #   # plot_layout(ncol = 2) +
  #   plot_annotation(tag_levels = "a") +
  #   plot_layout(guides = "collect") &
  #   theme(legend.position = "bottom") +
  #     theme(text = element_text(size = 12)) +
  #     # theme(legend.title = element_text(color="134", size=16, face="bold"))+
  #     theme(text = element_text(family = "Times New Roman")) # ,face = "bold"

  # ggsave("Results/L0-time_VS_SINCERITITES.png",width = 8, height = 4, dpi =600)
}
