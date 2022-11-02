

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
  "hESC",
  "hHep",
  "mDC",
  "mESC",
  "mHSC-E",
  "mHSC-GM",
  "mHSC-L"
)
output <- "../scGRN-L0_output/output_scRNA-Seq/"
hnetwork_data_dir <- "../scGRN-L0_data/BEELINE-Networks/Networks/human/"
mnetwork_data_dir <- "../scGRN-L0_data/BEELINE-Networks/Networks/mouse/"
tfs <- read.csv("../scGRN-L0_data/BEELINE-Networks/human-tfs.csv")
evaluation_AUROC_all <- c()
evaluation_AUPRC_all <- c()
for (j in 1:length(data_path)) {
  simulation_data_dir <- paste0("../scGRN-L0_data/BEELINE-data/inputs/scRNA-Seq/", data_path[j], "/")
  simulation_data_file <- "ExpressionData.csv"
  simulation_PseudoTime_file <- "PseudoTime.csv"
  simulation_data <- read.csv(paste0(simulation_data_dir, simulation_data_file), row.names = 1)
  simulation_data <- simulation_data[tfs$TF[1:500], ] %>% na.omit()
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
      if (data_path[j] %in% c("hESC", "hHep")) {
        ground_truth_h(
          intput = paste0(output, "GRN_SINCERITITES.txt"),
          output = output,
          dataset_dir = hnetwork_data_dir,
          database = "STRING"
        )
      } else {
        ground_truth_m(
          intput = paste0(output, "GRN_SINCERITITES.txt"),
          output = output,
          dataset_dir = mnetwork_data_dir,
          database = "STRING"
        )
      }
      evaluationObject <- prepareEval(paste0(output, "GRN_SINCERITITES.txt"),
        paste0(paste0(output, "ground_truth.tsv")),
        totalPredictionsAccepted = 100000
      )
      AUROC_SINCERITITES <- calcAUROC(evaluationObject)
      AUPRC_SINCERITITES <- calcAUPR(evaluationObject)
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
  }
  if (T) {
    L0REG_L0_adjs <- matrix(0, ncol(data_grn), ncol(data_grn))
    rownames(L0REG_L0_adjs) <- colnames(data_grn)
    colnames(L0REG_L0_adjs) <- colnames(data_grn)
    n <- 10
    for (t in 1:n) {
      s <- floor(nrow(data_grn) * (t - 1) / 10) + 1
      e <- floor(nrow(data_grn) * t / 10)

      # s <- which(data_GENIE3$h <= max(PseudoTime) / n * (t - 1))
      # if (length(s) == 0) {
      #     s <- 1
      # } else {
      #     s <- s[length(s)]
      # }
      # e <- which(data_GENIE3$h <= max(PseudoTime) / n * t)
      # e <- e[length(e)]

      data <- data_grn[s:e, ]

      L0REG_L0_1 <- L0REG(
        matrix = t(data),
        regulators = colnames(data),
        targets = colnames(data),
        penalty = "L0"
      )
      # L0REG_L0_1$weight <- as.numeric(L0REG_L0_1$weight)
      # L0REG_L0_1 <- as.matrix(L0REG_L0_1)
      L0REG_L0_adj <- matrix(0, ncol(data_grn), ncol(data_grn))
      rownames(L0REG_L0_adj) <- colnames(data_grn)
      colnames(L0REG_L0_adj) <- colnames(data_grn)
      # L0REG_L0_adj[L0REG_L0_1[, 1:2]] <- L0REG_L0_1[, 3]
      for (g in 1:nrow(L0REG_L0_1)) {
        L0REG_L0_adj[L0REG_L0_1$regulatoryGene[g], L0REG_L0_1$targetGene[g]] <- L0REG_L0_1$weight[g]
      }
      L0REG_L0_adjs <- L0REG_L0_adjs + as.numeric(L0REG_L0_adj)
    }
    L0REG_L0_adjs <- L0REG_L0_adjs / max(L0REG_L0_adjs)
    weightdf <- GENIE3::getLinkList(L0REG_L0_adjs)
    # weightdf <- read.table("output_NIMEFI_L0.txt", header = F)
    names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
    weightdf <- weightdf[order(weightdf$weight), ]
    write.table(weightdf, file = paste0(output, "L0REG_L02.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
    evaluationObject <- prepareEval(paste0(output, "L0REG_L02.txt"),
      paste0(paste0(output, "ground_truth.tsv")),
      totalPredictionsAccepted = 100000
    )
    truth_network <- convertSortedRankTSVToAdjMatrix(paste0(output, "ground_truth.tsv"))
    AUCresult_L0REG_L0 <- auc_from_ranks_TC_sign(L0REG_L0_adjs, truth_network, 1000)
    L0REG_L0Dynamic_AUROC_S <- AUCresult_L0REG_L0$AUROC
    L0REG_L0Dynamic_AUPR_S <- AUCresult_L0REG_L0$AUPR
    evaluationObject <- prepareEval(paste0(output, "L0REG_L02.txt"),
      paste0(paste0(output, "ground_truth.tsv")),
      totalPredictionsAccepted = 100000
    )
    L0REG_L0Dynamic_AUROC <- calcAUROC(evaluationObject)
    L0REG_L0Dynamic_AUPRC <- calcAUPR(evaluationObject)
  }
  # --------------------------------------------------
  if (F) {
    NIMEFI(data_grn,
      GENIE = F, SVM = F, EL = T, penalty = "L0",
      outputFileName = paste0(output, "GRN_L0_N"),
      outputFileFormat = "txt",
      ELPredSampleMin = 20, ELPredSampleMax = 80,
      ELExpSampleMin = 20, ELExpSampleMax = 80,
      ELRankThreshold = 5, ELEnsembleSize = dim(data_grn)[2]
    )
    evaluationObject <- prepareEval(paste0(output, "GRN_L0_N.txt"),
      paste0(output, "ground_truth.tsv"),
      totalPredictionsAccepted = 100000
    )
    AUROC_L0_N <- calcAUROC(evaluationObject)
    AUPRC_L0_N <- calcAUPR(evaluationObject)
    NIMEFI(data_grn,
      GENIE = F, SVM = F, EL = T, penalty = "L0L2",
      outputFileName = paste0(output, "GRN_L0L2_N"),
      outputFileFormat = "txt",
      SVMRankThreshold = 5, SVMEnsembleSize = 100,
      ELPredSampleMin = 20, ELPredSampleMax = 80,
      ELExpSampleMin = 20, ELExpSampleMax = 80,
      ELRankThreshold = 5, ELEnsembleSize = dim(data_grn)[2]
    )
    evaluationObject <- prepareEval(paste0(output, "GRN_L0L2_N.txt"),
      paste0(paste0(output, "ground_truth.tsv")),
      totalPredictionsAccepted = 100000
    )
    AUROC_L0L2_N <- calcAUROC(evaluationObject)
    AUPRC_L0L2_N <- calcAUPR(evaluationObject)
  }
  if (T) {
    L0Dynamic <- L0REG(
      matrix = t(data_grn),
      regulators = colnames(data_grn),
      targets = colnames(data_grn),
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
  }
  # --------------------------------------------------
  library(dyngen)
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
    Dataset = paste0(data_path[j]),
    L0Dynamic = AUROC_L0,
    L0DynamicL2 = AUROC_L0L2,
    L0Dynamic2 = L0REG_L0Dynamic_AUROC,
    GENIE3 = AUROC_GENIE3,
    SINCERITITES = AUROC_SINCERITITES,
    PPCOR = AUROC_PPCOR,
    LEAP = AUROC_LEAP
  )
  evaluation_AUPRC <- data.frame(
    Dataset = paste0(data_path[j]),
    L0Dynamic = AUPRC_L0,
    L0DynamicL2 = AUPRC_L0L2,
    L0Dynamic2 = L0REG_L0Dynamic_AUPRC,
    GENIE3 = AUPRC_GENIE3,
    SINCERITITES = AUPRC_SINCERITITES,
    PPCOR = AUPRC_PPCOR,
    LEAP = AUPRC_LEAP
  )
  print(evaluation_AUROC)
  evaluation_AUROC_all <- rbind.data.frame(evaluation_AUROC_all, evaluation_AUROC)
  evaluation_AUPRC_all <- rbind.data.frame(evaluation_AUPRC_all, evaluation_AUPRC)
}
write.csv(evaluation_AUROC_all, paste0(output, "evaluation_AUROC_scRNA-Seq.csv"))
write.csv(evaluation_AUPRC_all, paste0(output, "evaluation_AUPRC_scRNA-Seq.csv"))
