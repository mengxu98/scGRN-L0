

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
  "GSD",
  "HSC",
  "mCAD",
  "VSC"
)
cell_drop <- c(
  "1", "1-50", "1-70",
  "2", "2-50", "2-70",
  "3", "3-50", "3-70",
  "4", "4-50", "4-70",
  "5", "5-50", "5-70",
  "6", "6-50", "6-70",
  "7", "7-50", "7-70",
  "8", "8-50", "8-70",
  "9", "9-50", "9-70",
  "10", "10-50", "10-70"
)
output <- "../scGRN-L0_output/output_Curated/"
evaluation_AUROC_all <- c()
evaluation_AUPRC_all <- c()
for (j in 1:length(data_path)) {
  evaluation_AUROC_one <- c()
  evaluation_AUPRC_one <- c()
  for (i in 1:length(cell_drop)) {
    simulation_data_dir <- paste0("../scGRN-L0_data/BEELINE-data/inputs/Curated/", data_path[j], "/", data_path[j], "-2000-", cell_drop[i], "/")
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
    # --------------------------------------------------
    if (T) {
      DATA <- uploading(paste0(simulation_data_dir, "ExpressionData_all", ".csv"))
      # result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 1)
      result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 1, SIGN = 0)
      adj_matrix <- result$adj_matrix / max(result$adj_matrix)
      if (adj_matrix[1, 1] == "NaN") {
        SINCERITITES_AUROC_S <- 0
        AUROC_SINCERITITES <- 0
      } else {
        table <- final_ranked_predictions(adj_matrix, DATA$genes, SIGN = 1, saveFile = FALSE)
        table <- as.data.frame(table)
        write.table(table[, -4], paste0(output, "GRN_SINCERITITES.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
        ground_truth_simulation(
          intput = paste0(output, "GRN_SINCERITITES.txt"),
          output = output,
          dataset_dir = simulation_data_dir,
          file = "refNetwork.csv"
        )
        evaluationObject <- prepareEval(paste0(output, "GRN_SINCERITITES.txt"),
          paste0(paste0(output, "ground_truth.tsv")),
          totalPredictionsAccepted = 100000
        )
        # evaluationObject <- prepareEval(paste0(output, "ground_pred.txt"),
        #   paste0(paste0(output, "ground_truth.tsv")),
        #   totalPredictionsAccepted = 100000
        # )
        AUROC_SINCERITITES <- calcAUROC(evaluationObject)
        AUPRC_SINCERITITES <- calcAUPR(evaluationObject)
        # adj_matrix <- na.omit(adj_matrix)
        # truth_network <- convertSortedRankTSVToAdjMatrix(paste0(output, "ground_truth.tsv"))
        # AUCresult_SINCERITITES <- auc_from_ranks_TC_sign(adj_matrix, truth_network, 1000)
        # SINCERITITES_AUROC_S <- AUCresult_SINCERITITES$AUROC
        # SINCERITITES_AUPR_S <- AUCresult_SINCERITITES$AUPR
      }
    }
    # --------------------------------------------------
    data_grn <- read.csv(paste0(simulation_data_dir, "ExpressionData_all.csv"),
      header = T
    ) %>% as.matrix()
    PseudoTime <- data_grn[, ncol(data_grn)]
    data_grn <- data_grn[, -ncol(data_grn)]
    if (F) {
      NIMEFI(data_grn,
        GENIE = T, SVM = F, EL = F,
        outputFileName = paste0(output, "GRN_GENIE3_N"),
        outputFileFormat = "txt"
      )
      evaluationObject <- prepareEval(paste0(output, "GRN_GENIE3_N.txt"),
        paste0(paste0(output, "ground_truth.tsv")),
        totalPredictionsAccepted = 100000
      )
      AUROC_GENIE3_N <- calcAUROC(evaluationObject)
      AUPRC_GENIE3_N <- calcAUPR(evaluationObject)
    }
    if (T) {
      library(GENIE3)
      weightMat <- GENIE3(
        exprMatrix = t(data_grn),
        nCores = 32
      )
      # AUCresult_GENIE3 <- auc_from_ranks_TC_sign(weightMat, truth_network, 1000)
      # AUROC_GENIE3_S <- AUCresult_GENIE3$AUROC
      # AUPRC_GENIE3_S <- AUCresult_GENIE3$AUPR
      weightdf <- getLinkList(weightMat)
      names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
      write.table(weightdf, file = paste0(output, "GRN_GENIE3.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
      evaluationObject <- prepareEval(paste0(output, "GRN_GENIE3.txt"),
        paste0(paste0(output, "ground_truth.tsv")),
        totalPredictionsAccepted = 100000
      )
      AUROC_GENIE3 <- calcAUROC(evaluationObject)
      AUPRC_GENIE3 <- calcAUPR(evaluationObject)
    }

    # --------------------------------------------------
    if (F) {
      NIMEFI(data_grn,
        GENIE = F, SVM = F, EL = T, penalty = "L0",
        outputFileName = paste0(output, "output_NIMEFI_L0"),
        outputFileFormat = "txt",
        ELPredSampleMin = 20, ELPredSampleMax = 80,
        ELExpSampleMin = 20, ELExpSampleMax = 80,
        ELRankThreshold = 5, ELEnsembleSize = dim(data_grn)[1]
      )
      evaluationObject <- prepareEval(paste0(output, "output_NIMEFI_L0.txt"),
        paste0(output, "ground_truth.tsv"),
        totalPredictionsAccepted = 100000
      )
      AUROC_L0_N <- calcAUROC(evaluationObject)
      AUPRC_L0_N <- calcAUPR(evaluationObject)
      NIMEFI(data_grn,
        GENIE = F, SVM = F, EL = T, penalty = "L0L2",
        outputFileName = paste0(output, "output_NIMEFI_L0"),
        outputFileFormat = "txt",
        SVMRankThreshold = 5, SVMEnsembleSize = 100,
        ELPredSampleMin = 20, ELPredSampleMax = 80,
        ELExpSampleMin = 20, ELExpSampleMax = 80,
        ELRankThreshold = 5, ELEnsembleSize = dim(data_grn)[1]
      )
      evaluationObject <- prepareEval(paste0(output, "output_NIMEFI_L0.txt"),
        paste0(paste0(output, "ground_truth.tsv")),
        totalPredictionsAccepted = 100000
      )
      AUROC_L0L2_N <- calcAUROC(evaluationObject)
      AUPRC_L0L2_N <- calcAUPR(evaluationObject)
    }

    if (T) {
      L0REG_L0 <- L0REG(
        matrix = t(data_grn),
        regulators = colnames(data_grn),
        targets = colnames(data_grn),
        # maxSuppSize = 5,
        penalty = "L0"
      )
      write.table(L0REG_L0,
        paste0(output, "GRN_L0.txt"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = F
      )
      evaluationObject <- prepareEval(paste0(output, "GRN_L0.txt"),
        paste0(paste0(output, "ground_truth.tsv")),
        totalPredictionsAccepted = 100000
      )
      AUROC_L0 <- calcAUROC(evaluationObject)
      AUPRC_L0 <- calcAUPR(evaluationObject)
      # L0REG_L0$weight <- as.numeric(L0REG_L0$weight)
      # L0REG_L0$regulatoryGene <- as.factor(L0REG_L0$regulatoryGene)
      # L0REG_L0$targetGene <- as.factor(L0REG_L0$targetGene)
      # L0REG_L0 <- as.matrix(L0REG_L0)
      # L0REG_L0_adj <- matrix(0, ncol(data_grn), ncol(data_grn))
      # rownames(L0REG_L0_adj) <- colnames(data_grn)
      # colnames(L0REG_L0_adj) <- colnames(data_grn)
      # L0REG_L0_adj[L0REG_L0[, 1:2]] <- L0REG_L0[, 3]
      # L0REG_L0_adj <- as.numeric(L0REG_L0_adj)
      # truth_network <- convertSortedRankTSVToAdjMatrix(paste0(output, "ground_truth.tsv"))
      # AUCresult_L0REG_L0 <- auc_from_ranks_TC_sign(L0REG_L0_adj, truth_network, 1000)
      # L0REG_L0Dynamic_AUROC_S <- AUCresult_L0REG_L0$AUROC
      # L0REG_L0Dynamic_AUPR_S <- AUCresult_L0REG_L0$AUPR
      # for (g in 1:nrow(L0REG_L0)) {
      #   L0REG_L0_adj[L0REG_L0$regulatoryGene[g], L0REG_L0$targetGene[g]] <- L0REG_L0$weight[g]
      # }
      # L0REG_L0 <- L0REG(t(data_grn),
      #   regulators = colnames(data_grn),
      #   targets = colnames(data_grn), penalty = "L0"
      # )
      # AUCresult_L0REG <- auc_from_ranks_TC_sign(L0REG_L0, truth_network, 1000)
      # AUROC_L0REG_L0_S <- AUCresult_L0REG$AUROC
      # AUROC_L0REG_L0_S
      # --------------------------------------------------
      L0REG_L0L2 <- L0REG(t(data_grn),
        regulators = colnames(data_grn),
        targets = colnames(data_grn),
        # maxSuppSize = 5,
        penalty = "L0L2"
      )
      write.table(L0REG_L0L2,
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
      # L0REG_L0L2 <- L0REG(t(data_grn),
      #   # regulators = colnames(data_grn),
      #   targets = colnames(data_grn), penalty = "L0L2"
      # )
      # AUCresult_L0REG <- auc_from_ranks_TC_sign(L0REG_L0L2, truth_network, 1000)
      # AUROC_L0REG_L0L2_S <- AUCresult_L0REG$AUROC
    }
    # --------------------------------------------------
    if (T) {
      library(ppcor)
      inputExpr <- t(data_grn)
      geneNames <- rownames(inputExpr)
      rownames(inputExpr) <- c(geneNames)
      pcorResults <- pcor(x = t(as.matrix(inputExpr)), method = "spearman")
      DF <- data.frame(
        Gene1 = geneNames[c(row(pcorResults$estimate))], Gene2 = geneNames[c(col(pcorResults$estimate))],
        corVal = c(pcorResults$estimate), pValue = c(pcorResults$p.value)
      )
      outDF <- DF[order(DF$corVal, decreasing = TRUE), ]
      write.table(outDF[-(1:8), ],
        paste0(output, "GRN_PPCOR.txt"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = F
      )
      evaluationObject <- prepareEval(paste0(output, "GRN_PPCOR.txt"),
        paste0(paste0(output, "ground_truth.tsv")),
        totalPredictionsAccepted = 100000
      )
      AUROC_PPCOR <- calcAUROC(evaluationObject)
      AUPRC_PPCOR <- calcAUPR(evaluationObject)
    }
    # --------------------------------------------------
    if (T) {
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
      write.table(outDF,
        paste0(output, "GRN_LEAP.txt"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = F
      )
      evaluationObject <- prepareEval(paste0(output, "GRN_LEAP.txt"),
        paste0(paste0(output, "ground_truth.tsv")),
        totalPredictionsAccepted = 100000
      )
      AUROC_LEAP <- calcAUROC(evaluationObject)
      AUPRC_LEAP <- calcAUPR(evaluationObject)
    }
    # --------------------------------------------------
    evaluation_AUROC <- data.frame(
      Dataset = paste0(data_path[j], "-2000-", cell_drop[i]),
      L0REG_L0 = AUROC_L0,
      L0REG_L0L2 = AUROC_L0L2,
      GENIE3 = AUROC_GENIE3,
      SINCERITITES = AUROC_SINCERITITES,
      PPCOR = AUROC_PPCOR,
      LEAP = AUROC_LEAP
    )
    evaluation_AUPRC <- data.frame(
      Dataset = paste0(data_path[j], "-2000-", cell_drop[i]),
      L0REG_L0 = AUPRC_L0,
      L0REG_L0L2 = AUPRC_L0L2,
      GENIE3 = AUPRC_GENIE3,
      SINCERITITES = AUPRC_SINCERITITES,
      PPCOR = AUPRC_PPCOR,
      LEAP = AUPRC_LEAP
    )
    evaluation_AUROC_one <- rbind.data.frame(evaluation_AUROC_one, evaluation_AUROC)
    evaluation_AUPRC_one <- rbind.data.frame(evaluation_AUPRC_one, evaluation_AUPRC)
    print(evaluation_AUROC_one)
  }
  evaluation_AUROC_all <- rbind.data.frame(evaluation_AUROC_all, evaluation_AUROC_one)
  evaluation_AUPRC_all <- rbind.data.frame(evaluation_AUPRC_all, evaluation_AUPRC_one)
}

write.csv(evaluation_AUROC_all[, -1], paste0(output, "evaluation_AUROC.csv"), row.names = F)
write.csv(evaluation_AUPRC_all[, -1], paste0(output, "evaluation_AUPRC.csv"), row.names = F)

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
    "GSD",
    "HSC",
    "mCAD",
    "VSC"
  )
  evaluation_AUROC_all <- read.csv(paste0(output, "evaluation_AUROC.csv"))
  head(evaluation_AUROC_all[1:3, 1:3])

  for (i in 1:length(data_path)) {
    dataset <- data_path[i]
    evaluation_AUROC_dataset <- evaluation_AUROC_all[grep(dataset, evaluation_AUROC_all$Dataset), ]
    names(evaluation_AUROC_dataset) <- c("Dataset", "L0-Dynamic", "L0-Dynamic2", "GENIE3", "SINCERITITES")
    methods_barplot_all <- evaluation_AUROC_dataset %>%
      as.data.frame() %>%
      pivot_longer(
        cols = 2:c(ncol(evaluation_AUROC_dataset)),
        names_to = "Methods",
        values_to = "AUROC"
      )

    my_comparisons <- list(
      c("L0REG_L0", "GENIE3"),
      c("L0REG_L0", "PPCOR"),
      c("L0REG_L0", "SINCERITITES"),
      c("L0REG_L0", "LEAP")
    )

    mycol <- c("#3399cc", "#3366cc", "#ff00cc", "#cc0033")
    mycol <- c("black", "black", "gray", "white")
    p <- ggplot(
      methods_barplot_all,
      aes(x = Methods, y = AUROC)
    ) +
      stat_compare_means(
        method = "wilcox.test",
        label = "p.signif",
        comparisons = my_comparisons,
        bracket.size = 0.6,
        sizen = 4,
        color = "#6699cc"
      ) +
      labs(x = "Methods", y = "AUROC") +
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

    P1 <- ggplot(methods_barplot_all, aes(x = Methods, y = AUROC, fill = Methods)) +
      stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) + # 由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
      geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white") +
      geom_jitter(aes(fill = Methods), width = 0.2, shape = 21, size = 2.5) +
      scale_fill_manual(values = mycol) +
      scale_color_manual(values = mycol) +
      # scale_x_discrete(labels = c("L0-Dynamic", "GENIE3","SINCERITITES"))+
      ggtitle(" ") +
      theme_bw() +
      theme(legend.position = "bottom") +
      ylab("AUROC") +
      xlab("Methods")
    # P1

    results_10nets <- evaluation_AUROC_dataset
    results_10nets$Dataset <- c(
      "1", "1-50", "1-70",
      "2", "2-50", "2-70",
      "3", "3-50", "3-70",
      "4", "4-50", "4-70",
      "5", "5-50", "5-70",
      "6", "6-50", "6-70",
      "7", "7-50", "7-70",
      "8", "8-50", "8-70",
      "9", "9-50", "9-70",
      "10", "10-50", "10-70"
    )
    df_res10 <- melt(results_10nets, id = "Dataset", variable.name = "Method", value.name = "AUROC")
    df_res10$Method <- factor(df_res10$Method,
      levels = c("L0-Dynamic", "L0-Dynamic2", "GENIE3", "SINCERITITES")
    )

    p2 <- ggplot(df_res10, aes(x = Dataset, y = AUROC, fill = Method)) +
      geom_bar(stat = "identity", position = position_dodge(), color = "black", width = 0.8) +
      scale_fill_manual(values = mycol) +
      # geom_errorbar(aes(ymin=AUROC - 0.1, ymax=AUROC + 0.1), position = position_dodge(.6), width=.2)+
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

# p+p2
# ggsave("Results/Methods-contrast-3.png", width = 12, height = 4, dpi = 600)
