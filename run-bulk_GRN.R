

library("patchwork")
library("ggplot2")
library("reshape2")
library("ggpubr")
library("dplyr")
library("tidyr")
library("pheatmap")
library("tidyverse")
library("RColorBrewer")
library("cowplot")
library("bslib")
library("ggthemes")
theme_set(theme_pubclean())

# Load all the script files
source("framework_main.R")

# Dataset: GNW
evaluation_gnw_5_list <- list()
for (n in 1:2) {
  evaluation_gnw_5 <- c()
  for (i in 1:5) {
    if (T) {
      pathway <- paste0("gnw_5_compendia_5_gold_networks/", i) # i=1
      dataway <- "/measure/measure.tsv" #
      goldway <- "/gold/gold.tsv"
    } else {
      pathway <- paste0("gnw_10_compendia_1_gold_network/", i) # i=1
      dataway <- "/Ecoli-20_nonoise_multifactorial.tsv" #
      goldway <- "/Ecoli-20_goldstandard.tsv"
    }

    expression_dataset_test <- constructExpressionMatrixFromFile(paste0(pathway, dataway))
    # L0 0.7357716
    expression_dataset_test <- read.table(paste0(pathway, dataway), header = T) %>% as.matrix()
    # L0 0.7315969 GENIE3 0.6787165
    min(expression_dataset_test)

    # Run algorithm using default settings, writes result to output.txt
    NIMEFI(expression_dataset_test,
      GENIE = F, SVM = F, EL = T, # penalty = "L0L2",
      outputFileName = "output/output_net_raw_L0.txt",
      ELPredSampleMin = 20, ELPredSampleMax = 80,
      ELExpSampleMin = 20, ELExpSampleMax = 80,
      ELRankThreshold = 5, ELEnsembleSize = 200
    )

    evaluationObject <- prepareEval("output/output_net_raw_L0.txt", paste0(pathway, goldway))
    L0_AUROC <- calcAUROC(evaluationObject)
    L0_AUPR <- calcAUPR(evaluationObject)
    L0_AUROC # (abs) 0.7412103 # 0.5796524 # 0.6807636 0.7605543 0.8545814
    L0_AUPR
    # Run algorithm using default settings, writes result to output.txt
    NIMEFI(expression_dataset_test,
      GENIE = T, SVM = F, EL = F,
      outputFileName = "output/output_net_raw_GENIE3.txt"
    )

    evaluationObject <- prepareEval("output/output_net_raw_GENIE3.txt", paste0(pathway, goldway))
    GENIE3_AUROC <- calcAUROC(evaluationObject)
    GENIE3_AUPR <- calcAUPR(evaluationObject)
    GENIE3_AUROC

    evaluation_gnw_1 <- data.frame(
      Dataset = paste0("Net-", i),
      `L0Reg framework` = L0_AUROC,
      GENIE3 = GENIE3_AUROC
    )
    evaluation_gnw_5 <- rbind.data.frame(evaluation_gnw_5, evaluation_gnw_1)

    L0_AUROC
    GENIE3_AUROC

  }
  write.csv(evaluation_gnw_5, paste0("evaluation_gnw_5_", n, ".csv"))
  evaluation_gnw_5_list[[n]] <- evaluation_gnw_5
}

# Dataset: GNW
evaluation_gnw_10_list <- list()
for (n in 1:3) {
  evaluation_gnw_10 <- c()
  for (i in 1:10) {
    if (F) {
      pathway <- paste0("gnw_5_compendia_5_gold_networks/", i) # i=1
      dataway <- "/measure/measure.tsv" #
      goldway <- "/gold/gold.tsv"
    } else {
      pathway <- paste0("gnw_10_compendia_1_gold_network/", i) # i=1
      dataway <- "/Ecoli-20_nonoise_multifactorial.tsv" #
      goldway <- "/Ecoli-20_goldstandard.tsv"
    }

    expression_dataset_test <- constructExpressionMatrixFromFile(paste0(pathway, dataway))
    # L0: 0.7357716; GENIE3: 0.6944884
    expression_dataset_test <- read.table(paste0(pathway, dataway), header = T) %>% as.matrix()
    # L0: 0.7315969; GENIE3: 0.6787165
    min(expression_dataset_test)

    # Run algorithm using default settings, writes result to output.txt
    NIMEFI(expression_dataset_test,
      GENIE = F, SVM = F, EL = T,
      outputFileName = "output/output_net_raw_L0.txt",
      ELPredSampleMin = 20, ELPredSampleMax = 80,
      ELExpSampleMin = 20, ELExpSampleMax = 80,
      ELRankThreshold = 5, ELEnsembleSize = 100
    )

    evaluationObject <- prepareEval("output/output_net_raw_L0.txt", paste0(pathway, goldway))
    L0_AUROC <- calcAUROC(evaluationObject)
    L0_AUPR <- calcAUPR(evaluationObject)
    L0_AUROC

    # Run algorithm using default settings, writes result to output.txt
    NIMEFI(expression_dataset_test,
      GENIE = T, SVM = F, EL = F,
      outputFileName = "output/output_net_raw_GENIE3.txt"
    )

    evaluationObject <- prepareEval("output/output_net_raw_GENIE3.txt", paste0(pathway, goldway))
    GENIE3_AUROC <- calcAUROC(evaluationObject)
    GENIE3_AUPR <- calcAUPR(evaluationObject)
    GENIE3_AUROC

    evaluation_gnw_1 <- data.frame(
      Dataset = paste0("Net-", i),
      `L0Reg framework` = L0_AUROC,
      GENIE3 = GENIE3_AUROC
    )
    evaluation_gnw_10 <- rbind.data.frame(evaluation_gnw_10, evaluation_gnw_1)

    L0_AUROC
    GENIE3_AUROC
  }
  write.csv(evaluation_gnw_10, paste0("evaluation_gnw_10_", n, "_L0L2.csv"))
  evaluation_gnw_10_list[[n]] <- evaluation_gnw_10
}

# Dataset: GNW
if (F) {
  i <- 1
  pathway <- paste0("gnw_5_compendia_5_gold_networks/", i)
  dataway <- "/measure/measure.tsv"
  goldway <- "/gold/gold.tsv"
} else {
  i <- 1
  pathway <- paste0("gnw_10_compendia_1_gold_network/", i) # i=1
  dataway <- "/Ecoli-20_nonoise_multifactorial.tsv"
  goldway <- "/Ecoli-20_goldstandard.tsv"
}
expression_dataset_test <- constructExpressionMatrixFromFile(paste0(pathway, dataway))
expression_dataset_test <- read.table(paste0(pathway, dataway), header = T) %>% as.matrix()

# Run algorithm using default settings, writes result to output.txt
NIMEFI(expression_dataset_test,
  GENIE = FALSE,
  SVM = F,
  EL = T, penalty = "L0",
  outputFileName = "output/output_net1_raw_L0.txt",
  ELPredSampleMin = 20, ELPredSampleMax = 80,
  ELExpSampleMin = 20, ELExpSampleMax = 80,
  ELRankThreshold = 5, ELEnsembleSize = 100
)

evaluationObject <- prepareEval("output/output_net1_raw_L0.txt", paste0(pathway, goldway))
L0_AUROC <- calcAUROC(evaluationObject)
L0_AUPR <- calcAUPR(evaluationObject)
L0_AUROC
# L0: 0.85; L0L1: 0.6741838; L0L2: 0.8047406

# Run algorithm using default settings, writes result to output.txt
NIMEFI(expression_dataset_test,
  GENIE = T, SVM = F, EL = F,
  outputFileName = "output/output_net_raw_GENIE3.txt"
)

evaluationObject <- prepareEval("output/output_net_raw_GENIE3.txt", paste0(pathway, goldway))
GENIE3_AUROC <- calcAUROC(evaluationObject)
GENIE3_AUPR <- calcAUPR(evaluationObject)
GENIE3_AUROC

# test --------------------------------------------------------------------
for (d in 1:5) {
  # Read the expression dataset
  expression_dataset_test <- constructExpressionMatrixFromFile("gnw_10_compendia_1_gold_network/1/Ecoli-20_nonoise_multifactorial.tsv")
  expression_dataset_test <- read.table("gnw_10_compendia_1_gold_network/1/Ecoli-20_nonoise_multifactorial.tsv", header = T) %>% as.matrix()
  gold <- read.table("gnw_10_compendia_1_gold_network/9/Ecoli-20_goldstandard.tsv", header = F)

  dim(expression_dataset_test)
  expression_dataset_test <- abs(expression_dataset_test)
  min(expression_dataset_test)
  set.seed(2022)

  L0_AUROC_all <- c()
  GENIE3_AUROC_all <- c()
  for (i in 1:10) {
    # train_idx <- sample(ncol(expression_dataset_test), 50/dim(expression_dataset_test)[2] * ncol(expression_dataset_test))
    train_idx <- sample(ncol(expression_dataset_test), 0.5 * ncol(expression_dataset_test))
    train_data <- expression_dataset_test[, train_idx]

    gold_gene <- colnames(train_data)

    gold_new_genes <- c()
    for (g in 1:length(gold_gene)) {
      gold_gene1 <- gold[which(gold$V1 == gold_gene[g]), ]
      # gold_gene2 <- gold[which(gold$V2==gold_gene[g]),]
      gold_new_genes <- rbind.data.frame(gold_new_genes, gold_gene1)
    }

    gold_new <- c()
    for (g in 1:length(gold_gene)) {
      gold_gene2 <- gold_new_genes[which(gold_new_genes$V2 == gold_gene[g]), ]
      # gold_gene2 <- gold[which(gold$V2==gold_gene[g]),]
      gold_new <- rbind.data.frame(gold_new, gold_gene2)
    }

    write.table(gold_new, "exampleData/DREAM5Network1_new.tsv", quote = F, sep = "\t", row.names = F, col.names = F)
    gold_gene_new <- c(gold_new$V1, gold_new$V2)
    gold_gene_new <- unique(gold_gene_new)
    test_data <- train_data[, gold_gene_new]
    # Optional: Read a list of TF
    # tf <- getIndicesOfGenesInMatrix(expression_dataset,inputFile = "")

    # Run algorithm using default settings, writes result to output.txt
    NIMEFI(test_data,
      GENIE = F, SVM = F, EL = TRUE, outputFileName = "output/output_net_raw_test.txt",
      ELPredSampleMin = 20, ELPredSampleMax = 80,
      ELExpSampleMin = 20, ELExpSampleMax = 80,
      ELRankThreshold = 5, ELEnsembleSize = 200
    )

    evaluationObject <- prepareEval("output/output_net_raw_test.txt", "exampleData/DREAM5Network1_new.tsv")
    L0_AUROC_filter <- calcAUROC(evaluationObject)
    L0_AUPR <- calcAUPR(evaluationObject)
    L0_AUROC_filter

    # Run algorithm using default settings, writes result to output.txt
    NIMEFI(expression_dataset_test,
      GENIE = F, SVM = F, EL = TRUE, outputFileName = "output/output_net_raw.txt",
      ELPredSampleMin = 20, ELPredSampleMax = 80,
      ELExpSampleMin = 20, ELExpSampleMax = 80,
      ELRankThreshold = 5, ELEnsembleSize = 100
    )

    evaluationObject <- prepareEval("output/output_net_raw.txt", "gnw_10_compendia_1_gold_network/9/Ecoli-20_goldstandard.tsv")
    L0_AUROC <- calcAUROC(evaluationObject)
    L0_AUPR <- calcAUPR(evaluationObject)
    L0_AUROC

    L0_AUROC_all <- rbind.data.frame(L0_AUROC_all, L0_AUROC) %>% na.omit()

    NIMEFI(test_data, GENIE = T, SVM = F, EL = F, outputFileName = "output/output_net_raw.txt")
    # Calculate AUROC/AUPR
    evaluationObject <- prepareEval("output/output_net_raw.txt", "exampleData/DREAM5Network1_new.tsv")
    # Print AUROC/AUPR
    GENIE3_AUROC <- calcAUROC(evaluationObject)
    GENIE3_AUPR <- calcAUPR(evaluationObject)
    GENIE3_AUROC

    GENIE3_AUROC_all <- rbind.data.frame(GENIE3_AUROC_all, GENIE3_AUROC) %>% na.omit()
  }

  evaluation_AUROC <- cbind.data.frame(L0Reg = L0_AUROC_all, GENIE3 = GENIE3_AUROC_all)
  names(evaluation_AUROC) <- c("L0Reg framework", "GENIE3")
  mean(evaluation_AUROC$`L0Reg framework`)
  mean(evaluation_AUROC$GENIE3)


  evaluation_AUROC_box <- evaluation_AUROC %>%
    as.data.frame() %>%
    # rownames_to_column("sample") %>% #bug
    pivot_longer(
      cols = 1:2,
      names_to = "Method",
      values_to = "AUROC"
    )
  my_comparisons <- list(c("L0Reg framework", "GENIE3"))

  p1 <- ggplot(data = evaluation_AUROC_box, aes(x = Method, y = AUROC)) +
    geom_boxplot(aes(fill = Method)) +
    # ylab(paste0(target_gene," Expression"))+
    stat_compare_means() + # Calculate P-value
    theme_bw() +
    # geom_jitter(color="gray")+
    geom_jitter()

  Contrast_all <- read.table("G:/Research/ImmuCycReg/L0/results/Contrast_all.csv", header = T, row.names = 1, sep = ",")
  Contrast_all_box_data <- Contrast_all[, -c(1, 2)]
  names(Contrast_all_box_data) <- c("L0Reg framework", "GENIE3")
  mean(Contrast_all_box_data$`L0Reg framework`)
  mean(Contrast_all_box_data$GENIE3)

  Contrast_all_box <- Contrast_all_box_data %>%
    as.data.frame() %>%
    pivot_longer(
      cols = 1:2,
      names_to = "Method",
      values_to = "Correlation"
    )
  my_comparisons <- list(c("L0Reg framework", "GENIE3"))

  p2 <- ggplot(data = Contrast_all_box, aes(x = Method, y = Correlation)) +
    geom_boxplot(aes(fill = Method)) +
    # ylab(paste0(target_gene," Expression"))+
    stat_compare_means() + # Calculate P-value
    theme_bw() +
    # geom_jitter(color="gray")+
    geom_jitter()

  p1 + p2 +
    plot_annotation(tag_levels = "a") +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}
# example -----------------------------------------------------------------

NIMEFI(expression_dataset, GENIE = T, SVM = F, EL = F, outputFileName = "output/output_net_raw.txt")
# Calculate AUROC/AUPR
evaluationObject <- prepareEval("output/output_net_raw.txt", "exampleData_raw/gold.tsv")
# Print AUROC/AUPR
GENIE3_AUROC <- calcAUROC(evaluationObject)
GENIE3_AUPR <- calcAUPR(evaluationObject)
GENIE3_AUROC
# raw: 0.7330426 0.6924941
# abs: 0.6196811

evaluation_all <- c()
RankThreshold <- c(5, 10)
EnsembleSize <- c(100, 200, 500, 1000, 1500)
for (r in 1:length(RankThreshold)) {
  rank <- RankThreshold[r]
  evaluation_rank <- c()
  for (i in 1:length(EnsembleSize)) {
    ensemble <- EnsembleSize[i]
    evaluation_dataset <- c()
    for (d in 1:10) {
      expression_dataset <- constructExpressionMatrixFromFile(paste0(
        "gnw_10_compendia_1_gold_network/",
        d,
        "/Ecoli-20_multifactorial.tsv"
      ))

      # Run algorithm using default settings, writes result to output.txt
      NIMEFI(expression_dataset,
        GENIE = F, SVM = F, EL = T,
        outputFileName = "output/output_L0.txt",
        ELPredSampleMin = 20, ELPredSampleMax = 80,
        ELExpSampleMin = 20, ELExpSampleMax = 80,
        ELRankThreshold = rank,
        ELEnsembleSize = ensemble
      )

      evaluationObject_L0 <- prepareEval(
        "output/output_L0.txt",
        paste0(
          "gnw_10_compendia_1_gold_network/",
          d,
          "/Ecoli-20_goldstandard.tsv"
        )
      )

      AUROC_L0 <- calcAUROC(evaluationObject_L0)
      AUPR_L0 <- calcAUPR(evaluationObject_L0)

      message("AUROC_L0: ", AUROC_L0, "; ", "AUPR_L0: ", AUPR_L0, "; ")
      evaluation <- data.frame(
        RankThreshold = rank,
        EnsembleSize = ensemble,
        Dataset = paste0("gnw_10_compendia_1_gold_network/", d),
        AUROC_L0 = AUROC_L0, AUPR_L0 = AUPR_L0
      )
      evaluation_dataset <- rbind.data.frame(evaluation_dataset, evaluation)
    }
    evaluation_rank <- rbind.data.frame(evaluation_rank, evaluation_dataset)
  }

  evaluation_all <- rbind.data.frame(evaluation_all, evaluation_rank)
}

evaluation_all_gnw_10 <- c()
for (i in 1:10) {
  expression_dataset <- constructExpressionMatrixFromFile(paste0(
    "gnw_10_compendia_1_gold_network/",
    i,
    "/Ecoli-20_multifactorial.tsv"
  ))
  # expression_dataset <- expression_dataset - min(expression_dataset)
  # Run algorithm using default settings, writes result to output.txt
  NIMEFI(expression_dataset,
    GENIE = F, SVM = F, EL = T,
    outputFileName = "output/output_L0.txt",
    ELPredSampleMin = 20, ELPredSampleMax = 80,
    ELExpSampleMin = 20, ELExpSampleMax = 80,
    ELRankThreshold = 5, ELEnsembleSize = 100
  )

  NIMEFI(expression_dataset,
    GENIE = T, SVM = F, EL = F,
    outputFileName = "output/output_GENIE3.txt",
    ELPredSampleMin = 20, ELPredSampleMax = 80,
    ELExpSampleMin = 20, ELExpSampleMax = 80,
    ELRankThreshold = 5, ELEnsembleSize = 100
  )

  evaluationObject_L0 <- prepareEval(
    "output/output_L0.txt",
    paste0(
      "gnw_10_compendia_1_gold_network/",
      i,
      "/Ecoli-20_goldstandard.tsv"
    )
  )

  evaluationObject_GENIE3 <- prepareEval(
    "output/output_GENIE3.txt",
    paste0(
      "gnw_10_compendia_1_gold_network/",
      i,
      "/Ecoli-20_goldstandard.tsv"
    )
  )

  AUROC_L0 <- calcAUROC(evaluationObject_L0)
  AUPR_L0 <- calcAUPR(evaluationObject_L0)

  AUROC_GENIE3 <- calcAUROC(evaluationObject_GENIE3)
  AUPR_GENIE3 <- calcAUPR(evaluationObject_GENIE3)

  message(
    "AUROC_L0: ", AUROC_L0, "; ", "AUPR_L0: ", AUPR_L0, "; ",
    "AUROC_GENIE3: ", AUROC_GENIE3, "; ", "AUPR_GENIE3: ", AUPR_GENIE3
  )
  evaluation <- data.frame(
    Dataset = paste0("Ecoli-20_net", i),
    AUROC_L0 = AUROC_L0, AUPR_L0 = AUPR_L0,
    AUROC_GENIE3 = AUROC_GENIE3, AUPR_GENIE3 = AUPR_GENIE3
  )
  evaluation_all_gnw_10 <- rbind.data.frame(evaluation_all_gnw_10, evaluation)
}

evaluation_all_gnw_5 <- c()
for (i in 1:5) {
  expression_dataset <- constructExpressionMatrixFromFile(paste0(
    "gnw_5_compendia_5_gold_networks/",
    i,
    "/measure/measure.tsv"
  ))
  expression_dataset <- expression_dataset - min(expression_dataset)
  # Run algorithm using default settings, writes result to output.txt
  NIMEFI(expression_dataset,
    GENIE = F, SVM = F, EL = T, outputFileName = "output/output_L0.txt",
    ELPredSampleMin = 20, ELPredSampleMax = 80,
    ELExpSampleMin = 20, ELExpSampleMax = 100,
    ELRankThreshold = 5, ELEnsembleSize = 100
  )

  NIMEFI(expression_dataset,
    GENIE = T, SVM = F, EL = F, outputFileName = "output/output_GENIE3.txt",
    ELPredSampleMin = 20, ELPredSampleMax = 80,
    ELExpSampleMin = 20, ELExpSampleMax = 80,
    ELRankThreshold = 5, ELEnsembleSize = 100
  )

  evaluationObject_L0 <- prepareEval(
    "output/output_L0.txt",
    paste0(
      "gnw_5_compendia_5_gold_networks/",
      i,
      "/gold/gold.tsv"
    )
  )

  evaluationObject_GENIE3 <- prepareEval(
    "output/output_GENIE3.txt",
    paste0(
      "gnw_5_compendia_5_gold_networks/",
      i,
      "/gold/gold.tsv"
    )
  )

  AUROC_L0 <- calcAUROC(evaluationObject_L0)
  AUPR_L0 <- calcAUPR(evaluationObject_L0)

  AUROC_GENIE3 <- calcAUROC(evaluationObject_GENIE3)
  AUPR_GENIE3 <- calcAUPR(evaluationObject_GENIE3)

  message(
    "AUROC_L0: ", AUROC_L0, "; ", "AUPR_L0: ", AUPR_L0, "; ",
    "AUROC_GENIE3: ", AUROC_GENIE3, "; ", "AUPR_GENIE3: ", AUPR_GENIE3
  )
  evaluation <- data.frame(
    Dataset = paste0("Ecoli-20_net", i),
    AUROC_L0 = AUROC_L0, AUPR_L0 = AUPR_L0,
    AUROC_GENIE3 = AUROC_GENIE3, AUPR_GENIE3 = AUPR_GENIE3
  )
  evaluation_all_gnw_5 <- rbind.data.frame(evaluation_all_gnw_5, evaluation)
}
