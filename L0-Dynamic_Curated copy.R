

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
source("Function-L0REG.R", chdir = TRUE)


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
output <- "../scGRN-L0_output/test"
evaluation_infromations_all <- c()
for (j in 1:length(data_path)) {
    evaluation_infromations <- c()
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
        data_GENIE3 <- read.csv(paste0(simulation_data_dir, "ExpressionData_all.csv"),
            header = T
        )
        data_GENIE3 <- data_GENIE3[order(data_GENIE3$h), ]
        PseudoTime <- data_GENIE3[, ncol(data_GENIE3)]
        data_GENIE31 <- data_GENIE3[, -ncol(data_GENIE3)] %>% as.matrix()
        if (T) {
            if (T) {
                L0REG_L0_adjs <- matrix(0, ncol(data_GENIE31), ncol(data_GENIE31))
                rownames(L0REG_L0_adjs) <- colnames(data_GENIE31)
                colnames(L0REG_L0_adjs) <- colnames(data_GENIE31)
                n <- 5
                for (t in 1:n) {
                    # s <- floor(nrow(data_GENIE31) * (t - 1) / 10) + 1
                    # e <- floor(nrow(data_GENIE31) * t / 10)
                    # data <- data_GENIE31[s:e, ]
                    # data <- data_GENIE31[which(data_GENIE3$h <= max(PseudoTime) / nrow(data_GENIE31) * 10 * t), ]
                    s <- which(data_GENIE3$h <= max(PseudoTime) / n * (t - 1))
                    if (length(s) == 0) {
                        s <- 1
                    } else {
                        s <- s[length(s)]
                    }

                    e <- which(data_GENIE3$h <= max(PseudoTime) / n * t)
                    e <- e[length(e)]
                    data <- data_GENIE31[s:e, ]

                    L0REG_L0_1 <- L0REG(t(data),
                        regulators = colnames(data),
                        targets = colnames(data),
                        penalty = "L0"
                    )
                    # L0REG_L0_1$weight <- as.numeric(L0REG_L0_1$weight)
                    # L0REG_L0_1 <- as.matrix(L0REG_L0_1)
                    L0REG_L0_adj <- matrix(0, ncol(data_GENIE31), ncol(data_GENIE31))
                    rownames(L0REG_L0_adj) <- colnames(data_GENIE31)
                    colnames(L0REG_L0_adj) <- colnames(data_GENIE31)
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
                write.table(weightdf, file = paste0(output, "L0REG_L02.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
                ground_truth_simulation(
                    intput = paste0(output, "L0REG_L02.txt"),
                    output = output,
                    dataset_dir = simulation_data_dir,
                    file = "refNetwork.csv"
                )
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
                L0REG_L0Dynamic_AUPR <- calcAUPR(evaluationObject)
            }
            L0Dynamic <- L0REG(
                matrix = t(data_GENIE31),
                regulators = colnames(data_GENIE31),
                targets = colnames(data_GENIE31),
                # maxSuppSize = 5,
                penalty = "L0"
            )
            write.table(L0Dynamic,
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
            # data_GENIE31 <- FQnorm(data_GENIE31)
            if (F) {
                L0REG_L0_adjs <- matrix(0, ncol(data_GENIE31), ncol(data_GENIE31))
                rownames(L0REG_L0_adjs) <- colnames(data_GENIE31)
                colnames(L0REG_L0_adjs) <- colnames(data_GENIE31)

                library(mclust, quietly = TRUE)
                cl1 <- Mclust(data_GENIE31)$classification
                library(RColorBrewer)
                # plot(data_GENIE31, col = brewer.pal(9, "Set1")[cl1], pch = 16, asp = 1)
                n <- length(table(cl1))

                # library(NMF)
                # seed <- "20210101"
                # data_nmf <- log(t(data_GENIE31) + 1, 10)
                # gene_no <- ncol(data_GENIE31)
                # mads <- apply(data_nmf, 1, mad)
                # data_nmf <- data_nmf[rev(order(mads)), ]
                # dataset <- data_nmf[1:gene_no, ]
                # res <- nmf(dataset, 2:6, nrun = 50, seed = seed)
                # plot(res)
                for (t in 1:n) {
                    data <- data_GENIE31[which(cl1 == t), ]
                    L0REG_L0_1 <- L0REG(t(data),
                        regulators = colnames(data),
                        targets = colnames(data),
                        penalty = "L0"
                    )
                    # L0REG_L0_1$weight <- as.numeric(L0REG_L0_1$weight)
                    # L0REG_L0_1 <- as.matrix(L0REG_L0_1)
                    L0REG_L0_adj <- matrix(0, ncol(data_GENIE31), ncol(data_GENIE31))
                    rownames(L0REG_L0_adj) <- colnames(data_GENIE31)
                    colnames(L0REG_L0_adj) <- colnames(data_GENIE31)
                    # L0REG_L0_adj[L0REG_L0_1[, 1:2]] <- L0REG_L0_1[, 3]
                    for (g in 1:nrow(L0REG_L0_1)) {
                        L0REG_L0_adj[L0REG_L0_1$regulatoryGene[g], L0REG_L0_1$targetGene[g]] <- L0REG_L0_1$weight[g]
                    }
                    L0REG_L0_adjs <- L0REG_L0_adjs + as.numeric(L0REG_L0_adj)
                }
                # L0REG_L0_adjs <- L0REG_L0_adjs / max(L0REG_L0_adjs)

                weightdf <- GENIE3::getLinkList(L0REG_L0_adjs)
                # weightdf <- read.table("output_NIMEFI_L0.txt", header = F)
                names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
                write.table(weightdf, file = paste0(output, "L0REG_L02.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
                ground_truth_simulation(
                    intput = paste0(output, "L0REG_L02.txt"),
                    output = output,
                    dataset_dir = simulation_data_dir,
                    file = "refNetwork.csv"
                )
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
                L0REG_L0Dynamic_AUPR <- calcAUPR(evaluationObject)
                print(
                    paste0(
                        L0REG_L0Dynamic_AUROC,
                        "---",
                        L0REG_L0Dynamic_AUROC_S
                    )
                )
            }

            if (F) {
                L0REG_L0_adjs <- matrix(0, ncol(data_GENIE31), ncol(data_GENIE31))
                rownames(L0REG_L0_adjs) <- colnames(data_GENIE31)
                colnames(L0REG_L0_adjs) <- colnames(data_GENIE31)
                win <- 100
                for (t in 1:(nrow(data_GENIE31) - win)) {
                    s <- t
                    e <- win + t
                    if (e < dim(data_GENIE31)[1]) {
                        data <- data_GENIE31[s:e, ]
                    } else {
                        data <- data_GENIE31
                    }
                    L0REG_L0_1 <- L0REG(t(data),
                        regulators = colnames(data),
                        targets = colnames(data), penalty = "L0"
                    )
                    L0REG_L0_1$weight <- as.numeric(L0REG_L0_1$weight)
                    L0REG_L0_1 <- as.matrix(L0REG_L0_1)
                    L0REG_L0_adj <- matrix(0, ncol(data_GENIE31), ncol(data_GENIE31))
                    rownames(L0REG_L0_adj) <- colnames(data_GENIE31)
                    colnames(L0REG_L0_adj) <- colnames(data_GENIE31)
                    L0REG_L0_adj[L0REG_L0_1[, 1:2]] <- L0REG_L0_1[, 3]
                    L0REG_L0_adjs <- L0REG_L0_adjs + as.numeric(L0REG_L0_adj)
                }
                # L0REG_L0_adjs <- L0REG_L0_adjs / max(L0REG_L0_adjs)
                AUCresult_L0REG_L0 <- auc_from_ranks_TC_sign(L0REG_L0_adjs, truth_network, 1000)
                L0REG_L0Dynamic_AUROC_S <- AUCresult_L0REG_L0$AUROC
                L0REG_L0Dynamic_AUPR_S <- AUCresult_L0REG_L0$AUPR

                weightdf <- GENIE3::getLinkList(L0REG_L0_adjs)
                # weightdf <- read.table("output_NIMEFI_L0.txt", header = F)
                names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
                write.table(weightdf, file = paste0(output, "L0REG_L02.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
                evaluationObject <- prepareEval(paste0(output, "L0REG_L02.txt"),
                    paste0(paste0(output, "ground_truth.tsv")),
                    totalPredictionsAccepted = 100000
                )
                L0REG_L0Dynamic_AUROC <- calcAUROC(evaluationObject)
                L0REG_L0Dynamic_AUPR <- calcAUPR(evaluationObject)

            }
        }

        # --------------------------------------------------
        evaluation_infromation <- data.frame(
            datasets = paste0(data_path[j], "-2000-", cell_drop[i]),
            AUROC_L0Dynamic = L0REG_L0Dynamic_AUROC,
            AUROC_L0Dynamic_S = L0REG_L0Dynamic_AUROC_S,
            AUROC_L0 = AUROC_L0
        )
        message(paste0("----- ", evaluation_infromation, " -----"))
        print(evaluation_infromation)
        evaluation_infromations <- rbind.data.frame(evaluation_infromations, evaluation_infromation)
    }
    evaluation_infromations_all <- rbind.data.frame(evaluation_infromations_all, evaluation_infromations)
    mean(evaluation_infromations_all$AUROC_L0Dynamic)
    mean(evaluation_infromations_all$AUROC_L0Dynamic_S)
}
evaluation_infromations_all[evaluation_infromations_all == 0] <- NA
na.omit(evaluation_infromations_all)
write.csv(evaluation_infromations_all, paste0(output, "evaluation_infromations_Curated_L0REG.csv"))

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
    evaluation_infromations_all <- read.csv("Results/evaluation_infromations_Curated_L0REG.csv")
    head(evaluation_infromations_all[1:3, 1:3])

    for (i in 1:length(data_path)) {
        dataset <- data_path[i]
        evaluation_infromations_GSD <- evaluation_infromations_all[grep(dataset, evaluation_infromations_all$datasets), ]
        # evaluation_infromations_GSD <- evaluation_infromations_GSD[, c("datasets", "AUROC_L0REG_L0_N", "AUROC_L0REG_L0L2_N", "AUROC_GENIE3", "AUROC_SINCERITITES")]
        evaluation_infromations_GSD <- evaluation_infromations_GSD[, -1]
        # names(evaluation_infromations_GSD) <- c("Dataset", "L0-Dynamic", "L0-Dynamic2", "GENIE3", "SINCERITITES")
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
            aes(x = Methods, y = AUROC)
        ) +
            # guides(fill = guide_legend(title = NULL)) +
            # stat_compare_means(
            #     method = "wilcox.test",
            #     label = "p.signif",
            #     comparisons = my_comparisons,
            #     bracket.size = 0.6,
            #     sizen = 4,
            #     color = "#6699cc"
            # ) +
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
        ggsave(paste0("Results/Methods-contrast-", dataset, "-1_L0REG.png"), width = 9, height = 4, dpi = 600)

        # methods_barplot_all %>% ggplot(., aes(x = Methods, y = AUROC, colour = Methods)) +
        #   geom_boxplot() +
        #   theme_bw() +
        #   theme(legend.position = "none")

        P1 <- ggplot(methods_barplot_all, aes(x = Methods, y = AUROC, fill = Methods)) + # ???fill=?????????????????????
            stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) + # ??????????????????????????????????????????????????????????????????????????????????????????
            geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white") + # size???????????????????????????????????????????????????fill?????????????????????outlier.fill???outlier.color????????????????????????
            geom_jitter(aes(fill = Methods), width = 0.2, shape = 21, size = 2.5) + # ?????????????????????????????????????????????width??????????????????????????????????????????????????????
            scale_fill_manual(values = c("black", "gray", "white", "#c6524a", "#eabf00", "#696969")) + # ?????????????????????
            scale_color_manual(values = c("black", "#2874c5", "#008a00", "#c6524a", "#eabf00", "#696969")) + # ??????????????????????????????????????????
            # scale_x_discrete(labels = c("L0-Dynamic", "GENIE3","SINCERITITES"))+
            ggtitle(" ") + # ??????????????????
            theme_bw() +
            theme(legend.position = "bottom") +
            # theme(legend.position="none", #???????????????
            #       axis.text.x=element_text(colour="black",family="Times",size=14), #??????x??????????????????????????????
            #       axis.text.y=element_text(family="Times",size=14,face="plain"), #??????x??????????????????????????????
            #       axis.title.y=element_text(family="Times",size = 14,face="plain"), #??????y???????????????????????????
            #       axis.title.x=element_text(family="Times",size = 14,face="plain"), #??????x???????????????????????????
            #       plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #??????????????????????????????
            #       # panel.grid.major = element_blank(), #??????????????????
            #       panel.grid.minor = element_blank())+
            ylab("AUROC") +
            xlab("Methods") # ??????x??????y????????????
        # P1
        # ggsave("Results/Methods-contrast-1.png", width = 4, height = 4, dpi = 600)
        # ggsave(paste0("Results/", "Methods-contrast-cur.png"), width = 11, height = 8, dpi = 600)
        # ggsave(p, filename = paste0("Results/", "Methods-contrast-cur.png"))

        mycol <- c("black", "gray", "white")
        results_10nets <- evaluation_infromations_GSD
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
            levels = c("L0-Dynamic", "L0-Dynamic2", "GENIE3", "SINCERITITES"),
            labels = c("L0-Dynamic", "L0-Dynamic2", "GENIE3", "SINCERITITES")
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
        ggsave(paste0("Results/Methods-contrast-", dataset, "-2_L0REG.png"), width = 8, height = 4, dpi = 600)
    }
}

# p+p2
# ggsave("Results/Methods-contrast-3.png", width = 12, height = 4, dpi = 600)
