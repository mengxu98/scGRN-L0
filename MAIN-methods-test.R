

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

output <- "output_scRNA-Seq/"

evaluation_infromations_all <- c()
simulation_data_dir <- paste0("../scGRN-L0_data/mESC/")

simulation_data_file <- "counts.csv"
simulation_PseudoTime_file <- "dpt_time.txt"
simulation_data_celltype <- "anno.csv"

simulation_data <- read.csv(paste0(simulation_data_dir, simulation_data_file), row.names = 1)
head(simulation_data[1:3, 1:3])
simulation_PseudoTime <- read.table(paste0(simulation_data_dir, simulation_PseudoTime_file), row.names = 1)
simulation_celltype <- read.csv(paste0(simulation_data_dir, simulation_data_celltype), row.names = 1)
truth_network <- read.csv("../scGRN-L0_data/mESC/GRN.csv", row.names = 1)
truth_network_list <- GENIE3::getLinkList(as.matrix(truth_network))
for (i in 1:nrow(truth_network_list)) {
  if (truth_network_list$regulatoryGene[i] == truth_network_list$targetGene[i]) {
    print(1)
  }
}
apply(truth_network_list,1,function(x){
  if (truth_network_list$regulatoryGene[i] == truth_network_list$targetGene[i]) {
    print(1)
  }
})
write.table(truth_network_list, "../scGRN-L0_data/mESC/ground_truth.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
simulation_data <- cbind(simulation_data, celltype = simulation_celltype$label)
simulation_data <- cbind(simulation_data, h = simulation_PseudoTime$V2)
celltypes <- unique(simulation_celltype$label)
for (c in 1:length(celltypes)) {
  celltype <- celltypes[c]
  simulation_data_new <- simulation_data[which(simulation_data$celltype == celltype), ]
  simulation_data_new <- simulation_data_new[, -(ncol(simulation_data_new) - 1)]
  write.csv(simulation_data_new, paste0(simulation_data_dir, "ExpressionData_", celltype, ".csv"), row.names = FALSE)
}

DATA <- uploading(paste0(simulation_data_dir, "ExpressionData_", celltype, ".csv"))
# --------------------------------------------------
if (T) {
  # result_L0 <- SINCERITITES_L0(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 0, penalty = "L0L2", algorithm = "CD")
  result_L0 <- SINCERITITES_L0(DATA, distance = 1, method = 1, noDIAG = 1, SIGN = 0, penalty = "L0", algorithm = "CD")
  result_L0 <- SINCERITITES_L0(DATA, distance = 1, method = 1, noDIAG = 1, SIGN = 0, penalty = "L0L2", algorithm = "CD")
  adj_matrix_L0 <- result_L0$adj_matrix / max(result_L0$adj_matrix)
  table_L0 <- final_ranked_predictions(adj_matrix_L0, DATA$genes, SIGN = 1, saveFile = TRUE) %>% as.data.frame()
  # write.table(table_L0[, -4], paste0(output, "test_L0.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
  AUCresult_L0 <- auc_from_ranks_TC_sign(adj_matrix_L0, truth_network, 10000)
  AUROC_L0_S <- AUCresult_L0$AUROC
  AUPR_L0_S <- AUCresult_L0$AUPR
  AUROC_L0_S
  AUPR_L0_S
}
# --------------------------------------------------
if (T) {
  # result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 1)
  result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 1, SIGN = 0)
  adj_matrix <- result$adj_matrix / max(result$adj_matrix)
  table <- final_ranked_predictions(adj_matrix, DATA$genes, SIGN = 1, saveFile = TRUE)
  AUCresult_SINCERITITES <- auc_from_ranks_TC_sign(adj_matrix, truth_network, 1000)
  SINCERITITES_AUROC_S <- AUCresult_SINCERITITES$AUROC
  SINCERITITES_AUPR_S <- AUCresult_SINCERITITES$AUPR
  SINCERITITES_AUROC_S
}
# --------------------------------------------------
if (T) {
  data_GENIE3 <- read.csv(paste0(simulation_data_dir, "ExpressionData_", celltype, ".csv"),
    header = T
  ) %>% as.matrix()
  PseudoTime <- data_GENIE3[, ncol(data_GENIE3)]
  data_GENIE3 <- data_GENIE3[, -ncol(data_GENIE3)]
  library(GENIE3)
  weightMat <- GENIE3(
    exprMatrix = t(data_GENIE3),
    regulators = rownames(t(data_GENIE3)),
    nCores = 32
  )
  weightdf <- getLinkList(weightMat)
  names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
  write.table(weightdf, file = paste0("output_GENIE32.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  evaluationObject <- prepareEval(paste0("output_GENIE32.txt"), # paste0(output, "output_GENIE32.txt")
    paste0(paste0("../scGRN-L0_data/mESC/ground_truth.tsv")),
    totalPredictionsAccepted = 100000
  )
  GENIE3_AUROC2 <- calcAUROC(evaluationObject)
  GENIE3_AUPR2 <- calcAUPR(evaluationObject)
  GENIE3_AUROC2
  GENIE3_AUPR2
  AUCresult_GENIE3 <- auc_from_ranks_TC_sign(weightMat, truth_network, 1000)
  GENIE3_AUROC_S <- AUCresult_GENIE3$AUROC
  GENIE3_AUPR_S <- AUCresult_GENIE3$AUPR
  GENIE3_AUROC_S
}

if (T) {
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

  L0REG_L0 <- L0REG(t(data_GENIE3),
    # regulators = colnames(data_GENIE3),
    targets = colnames(data_GENIE3), penalty = "L0"
  )
  AUCresult_L0REG <- auc_from_ranks_TC_sign(L0REG_L0, truth_network, 1000)
  AUROC_L0REG_L0_S <- AUCresult_L0REG$AUROC
  AUROC_L0REG_L0_S
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

  L0REG_L0L2 <- L0REG(t(data_GENIE3),
    # regulators = colnames(data_GENIE3),
    targets = colnames(data_GENIE3), penalty = "L0L2"
  )
  AUCresult_L0REG <- auc_from_ranks_TC_sign(L0REG_L0L2, truth_network, 1000)
  AUROC_L0REG_L0L2_S <- AUCresult_L0REG$AUROC
  AUROC_L0REG_L0L2_S
}

# --------------------------------------------------
evaluation_infromation <- data.frame(
  datasets = paste0(data_path[j]),
  AUROC_SINCERITITES_L0 = AUROC_L0_N,
  # AUROC_SINCERITITES_L0_S = AUROC_L0_S,
  AUROC_SINCERITITES_L0L2 = AUROC_L0L2_N,
  # AUROC_SINCERITITES_L0L2_S = AUROC_L0L2_S,
  # AUROC_L0REG_L0_N = AUROC_L0REG_L0_N,
  # AUROC_L0REG_L0_S = AUROC_L0REG_L0_S,
  # AUROC_L0REG_L0L2_N = AUROC_L0REG_L0L2_N,
  # AUROC_L0REG_L0L2_S = AUROC_L0REG_L0L2_S,
  AUROC_NIMEFI_L0L2 = AUROC_NIMEFI_L0L2,
  AUROC_NIMEFI_L0 = AUROC_NIMEFI_L0,
  AUROC_GENIE3 = GENIE3_AUROC,
  AUROC_SINCERITITES = SINCERITITES_AUROC_N # ,
  # AUROC_SINCERITITES_S = SINCERITITES_AUROC_S
)

message(paste0("----- ", evaluation_infromation, " -----"))
# evaluation_infromations <- rbind.data.frame(evaluation_infromations, evaluation_infromation)

evaluation_infromations_all <- rbind.data.frame(evaluation_infromations_all, evaluation_infromation)
# evaluation_infromations_all
# evaluation_infromations_all[evaluation_infromations_all==0] <- NA
na.omit(evaluation_infromations_all)
write.csv(evaluation_infromations_all, "Results/evaluation_infromations_scRNA-Seq.csv")

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

  evaluation_infromations_all <- read.csv("Results/evaluation_infromations.csv")
  head(evaluation_infromations_all[1:3, 1:3])
  methods_barplot_all <- evaluation_infromations_all %>%
    as.data.frame() %>%
    pivot_longer(
      cols = 2:c(ncol(evaluation_infromations_all)),
      names_to = "Methods",
      values_to = "AUROC"
    )


  my_comparisons <- list(
    c("AUROC_NIMEFI_L0", "AUROC_GENIE3")
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
    # ylim(0, max(methods_barplot_all$AUROC) + 0.07) +
    labs(
      x = "Methods",
      y = "Percentage"
    ) +
    # stat_summary(
    #   fun.data = "mean_sdl",
    #   fun.args = list(mult = 1),
    #   geom = "pointrange",
    #   color = "gray"
    # ) +
    # geom_violin(aes(fill = datasets),
    #   trim = FALSE
    # ) +
    geom_boxplot(width = 0.5) +
    scale_fill_manual(values = c(
      "#2874c5",
      "#008a00",
      "#c6524a",
      "#eabf00",
      "#696969",
      "#008a00",
      "#c6524a",
      "#eabf00",
      "#696969",
      "#008a00",
      "#c6524a",
      "#eabf00",
      "#696969"
    )) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c(
      "#2874c5",
      "#008a00",
      "#c6524a",
      "#eabf00",
      "#696969",
      "#008a00",
      "#c6524a",
      "#eabf00",
      "#696969",
      "#008a00",
      "#c6524a",
      "#eabf00",
      "#696969"
    )) +
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
  ggsave(paste0("Results/", "Fig. 7.png"), width = 11, height = 8)
  ggsave(p, filename = paste0("Results/", "Fig. 7.png"))
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
  results_10nets <- evaluation_infromations_all[, c("datasets", "AUROC_L03", "SINCERITITES_AUROC", "GENIE3_AUROC")]
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
