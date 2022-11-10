
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
    output <- "../scGRN-L0_output/output_Curated/"
    data_path <- c(
        "GSD",
        "HSC",
        "mCAD",
        "VSC"
    )
    cell_num <- c(
        "50",
        "70"
    )
    evaluation_AUROC_all <- read.csv(paste0(output, "evaluation_AUROC.csv"))
    names(evaluation_AUROC_all) <- c("Dataset", "L0DWGRN", "GENIE3", "SINCERITITES", "PPCOR", "LEAP")
    my_comparisons <- list(
        c("L0DWGRN", "GENIE3"),
        c("L0DWGRN", "SINCERITITES"),
        c("L0DWGRN", "PPCOR"),
        c("L0DWGRN", "LEAP")
    )
    mycol <- c( "#3366cc","#008B00", "#008B8B", "#7171C6","#79CDCD")
    head(evaluation_AUROC_all[1:3, 1:3])
    for (d in 1:length(data_path)) {
        dataset <- data_path[d]
        evaluation_AUROC_dataset <- evaluation_AUROC_all[grep(dataset, evaluation_AUROC_all$Dataset), ]
        methods_barplot_all <- evaluation_AUROC_dataset %>%
            as.data.frame() %>%
            pivot_longer(
                cols = 2:c(ncol(evaluation_AUROC_dataset)),
                names_to = "Method",
                values_to = "AUROC"
            )
        methods <- names(evaluation_AUROC_dataset[, -1])
        methods_barplot_all$Method <- factor(methods_barplot_all$Method,
            levels = methods
        )
        p <- ggplot(
            methods_barplot_all,
            aes(x = Method, y = AUROC)
        ) +
            # geom_violin(aes(fill = Method),
            #     trim = FALSE
            # ) +
            geom_boxplot(aes(fill = Method),
                width = 0.8
            ) +
            stat_compare_means(
                method = "wilcox.test",
                label = "p.signif",
                comparisons = my_comparisons,
                bracket.size = 0.6,
                sizen = 4,
                color = "#6699cc"
            ) +
            scale_fill_manual(values = mycol) +
            # scale_color_manual(values = mycol) +
            scale_x_discrete(labels = methods) +
            labs(x = "Methods", y = "AUROC") +
            theme(legend.position = "bottom") +
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

        ggsave(paste0("../scGRN-L0_output/output_Curated/dataset/Methods-contrast-AUROC-", dataset, "-1.png"), width = 4, height = 4, dpi = 600)

        # methods_barplot_all %>% ggplot(., aes(x = Methods, y = AUROC, colour = Methods)) +
        #   geom_boxplot() +
        #   theme_bw() +
        #   # scale_x_discrete(labels = methods)+
        #   theme(legend.position = "none")

        P1 <- ggplot(methods_barplot_all, aes(x = Methods, y = AUROC, fill = Methods)) +
            stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
            geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white") +
            geom_jitter(aes(fill = Methods), width = 0.2, shape = 21, size = 2.5) +
            scale_fill_manual(values = mycol) +
            scale_color_manual(values = mycol) +
            scale_x_discrete(labels = methods) +
            ggtitle(" ") +
            theme_bw() +
            theme(legend.position = "bottom") +
            ylab("AUROC") +
            xlab("Methods")
        # P1
        df_res10 <- melt(evaluation_AUROC_dataset, id = "Dataset", variable.name = "Method", value.name = "AUROC")
        df_res10$Method <- factor(df_res10$Method, levels = methods)

        p2 <- ggplot(df_res10, aes(x = Dataset, y = AUROC, fill = Method)) +
            geom_bar(stat = "identity", position = position_dodge(), color = "black", width = 0.8) +
            scale_fill_manual(values = mycol) +
            # geom_errorbar(aes(ymin=AUROC - 0.1, ymax=AUROC + 0.1), position = position_dodge(.6), width=.2)+
            scale_x_discrete(labels = evaluation_AUROC_dataset$Dataset) +
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
        ggsave(paste0("../scGRN-L0_output/output_Curated/dataset/Methods-contrast-AUROC-", dataset, "-2.png"),
            width = 15, height = 4, dpi = 600
        )
    }

    for (c in 1:length(data_path)) {
        dataset <- cell_num[c]
        evaluation_AUROC <- evaluation_AUROC_all[grep(dataset, evaluation_AUROC_all$Dataset), ]
        methods_barplot_all <- evaluation_AUROC %>%
            as.data.frame() %>%
            pivot_longer(
                cols = 2:c(ncol(evaluation_AUROC)),
                names_to = "Method",
                values_to = "AUROC"
            )
        methods <- names(evaluation_AUROC[, -1])
        methods_barplot_all$Method <- factor(methods_barplot_all$Method,
            levels = methods
        )
        p <- ggplot(
            methods_barplot_all,
            aes(x = Method, y = AUROC)
        ) +
            # geom_violin(aes(fill = Method),
            #     trim = FALSE
            # ) +
            geom_boxplot(aes(fill = Method),
                width = 0.8
            ) +
            stat_compare_means(
                method = "wilcox.test",
                label = "p.signif",
                comparisons = my_comparisons,
                bracket.size = 0.6,
                sizen = 4,
                color = "#6699cc"
            ) +
            scale_fill_manual(values = mycol) +
            # scale_color_manual(values = mycol) +
            scale_x_discrete(labels = methods) +
            labs(x = "Methods", y = "AUROC") +
            theme(legend.position = "bottom") +
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
        ggsave(paste0("../scGRN-L0_output/output_Curated/cell/Methods-contrast-AUROC-", dataset, "-1.png"), width = 4, height = 4, dpi = 600)

        # methods_barplot_all %>% ggplot(., aes(x = Methods, y = AUROC, colour = Methods)) +
        #   geom_boxplot() +
        #   theme_bw() +
        #   # scale_x_discrete(labels = methods)+
        #   theme(legend.position = "none")

        P1 <- ggplot(methods_barplot_all, aes(x = Methods, y = AUROC, fill = Methods)) +
            stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
            geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white") +
            geom_jitter(aes(fill = Methods), width = 0.2, shape = 21, size = 2.5) +
            scale_fill_manual(values = mycol) +
            scale_color_manual(values = mycol) +
            scale_x_discrete(labels = methods) +
            ggtitle(" ") +
            theme_bw() +
            theme(legend.position = "bottom") +
            ylab("AUROC") +
            xlab("Methods")
        # P1
        df_res10 <- melt(evaluation_AUROC_dataset, id = "Dataset", variable.name = "Method", value.name = "AUROC")
        df_res10$Method <- factor(df_res10$Method,
            levels = methods
        )

        p2 <- ggplot(df_res10, aes(x = Dataset, y = AUROC, fill = Method)) +
            geom_bar(stat = "identity", position = position_dodge(), color = "black", width = 0.8) +
            scale_fill_manual(values = mycol) +
            # geom_errorbar(aes(ymin=AUROC - 0.1, ymax=AUROC + 0.1), position = position_dodge(.6), width=.2)+
            scale_x_discrete(labels = evaluation_AUROC_dataset$Dataset) +
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
        ggsave(paste0("../scGRN-L0_output/output_Curated/cell/Methods-contrast-AUROC-", dataset, "-2.png"), width = 15, height = 4, dpi = 600)
    }
}

# --------------------------------------------------
if (T) {
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(ggpubr)
    library(dplyr)
    library(tidyr)
    library(pheatmap)
    library(tidyverse)
    library(RColorBrewer)
    output <- "../scGRN-L0_output/output_Synthetic/"
    data_path <- c(
        "dyn-BF",
        # "dyn-BFC",
        "dyn-CY",
        "dyn-LI",
        # "dyn-TF",
        "dyn-LL"
    )
    cell_num <- c(
        "100",
        "200",
        "500",
        "2000",
        "5000"
    )
    evaluation_AUROC_all <- read.csv(paste0(output, "evaluation_AUROC.csv"))
    names(evaluation_AUROC_all) <- c("Dataset", "L0DWGRN", "GENIE3", "SINCERITITES", "PPCOR", "LEAP")
    head(evaluation_AUROC_all[1:3, 1:3])
    my_comparisons <- list(
        c("L0DWGRN", "GENIE3"),
        c("L0DWGRN", "SINCERITITES"),
        c("L0DWGRN", "PPCOR"),
        c("L0DWGRN", "LEAP")
    )
    mycol <- c( "#3366cc","#008B00", "#008B8B", "#7171C6","#79CDCD")
    for (c in 1:length(cell_num)) {
        cell <- cell_num[c]
        evaluation_AUROC_cell <- evaluation_AUROC_all[grep(cell, evaluation_AUROC_all$Dataset), ]
        for (d in 1:length(data_path)) {
            dataset <- data_path[d]
            evaluation_AUROC <- evaluation_AUROC_cell[grep(dataset, evaluation_AUROC_cell$Dataset), ]
            methods_barplot_all <- evaluation_AUROC %>%
                as.data.frame() %>%
                pivot_longer(
                    cols = 2:c(ncol(evaluation_AUROC)),
                    names_to = "Method",
                    values_to = "AUROC"
                )
            methods <- names(evaluation_AUROC[, -1])
            methods_barplot_all$Method <- factor(methods_barplot_all$Method,
                levels = methods
            )
            p <- ggplot(
                methods_barplot_all,
                aes(x = Method, y = AUROC)
            ) +
                # geom_violin(aes(fill = Method),
                #     trim = FALSE
                # ) +
                geom_boxplot(aes(fill = Method),
                    width = 0.8
                ) +
                stat_compare_means(
                    method = "wilcox.test",
                    label = "p.signif",
                    comparisons = my_comparisons,
                    bracket.size = 0.5,
                    tip.length = 0,
                    # color = "#6699cc",
                    sizen = 4
                ) +
                scale_fill_manual(values = mycol) +
                # scale_color_manual(values = mycol) +
                scale_x_discrete(labels = methods) +
                labs(x = "Methods", y = "AUROC") +
                theme(legend.position = "bottom") +
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
            ggsave(paste0("../scGRN-L0_output/output_Synthetic/cell-dataset/Methods-contrast-AUROC-", cell, "-", dataset, "-1.png"), width = 4, height = 4, dpi = 600)

            df_res10 <- melt(evaluation_AUROC_cell, id = "Dataset", variable.name = "Method", value.name = "AUROC")
            df_res10$Method <- factor(df_res10$Method,
                levels = methods
            )

            p2 <- ggplot(df_res10, aes(x = Dataset, y = AUROC, fill = Method)) +
                geom_bar(stat = "identity", position = position_dodge(), color = "black", width = 0.8) +
                scale_fill_manual(values = mycol) +
                # geom_errorbar(aes(ymin=AUROC - 0.1, ymax=AUROC + 0.1), position = position_dodge(.6), width=.2)+
                scale_x_discrete(labels = evaluation_AUROC_cell$Dataset) +
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
            # ggsave(paste0("../scGRN-L0_output/output_Synthetic/cell-dataset/Methods-contrast-AUROC-", cell, "-", dataset, "-2.png"), width = 15, height = 4, dpi = 600)
        }
    }

    for (d in 1:length(data_path)) {
        dataset <- data_path[d]
        evaluation_AUROC_dataset <- evaluation_AUROC_all[grep(dataset, evaluation_AUROC_all$Dataset), ]
        for (c in 1:length(cell_num)) {
            cell <- cell_num[c]
            evaluation_AUROC <- evaluation_AUROC_dataset[grep(cell, evaluation_AUROC_dataset$Dataset), ]

            methods_barplot_all <- evaluation_AUROC %>%
                as.data.frame() %>%
                pivot_longer(
                    cols = 2:c(ncol(evaluation_AUROC)),
                    names_to = "Method",
                    values_to = "AUROC"
                )

            methods <- names(evaluation_AUROC[, -1])
            methods_barplot_all$Method <- factor(methods_barplot_all$Method,
                levels = methods
            )
            p <- ggplot(
                methods_barplot_all,
                aes(x = Method, y = AUROC)
            ) +
                # geom_violin(aes(fill = Method),
                #     trim = FALSE
                # ) +
                geom_boxplot(aes(fill = Method),
                    width = 0.8
                ) +
                stat_compare_means(
                    method = "wilcox.test",
                    label = "p.signif",
                    comparisons = my_comparisons,
                    bracket.size = 0.6,
                    sizen = 4,
                    color = "#6699cc"
                ) +
                scale_fill_manual(values = mycol) +
                # scale_color_manual(values = mycol) +
                scale_x_discrete(labels = methods) +
                labs(x = "Methods", y = "AUROC") +
                theme(legend.position = "bottom") +
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
            ggsave(paste0("../scGRN-L0_output/output_Synthetic/dataset-cell/Methods-contrast-AUROC-", cell, "-", dataset, "-1.png"), width = 4, height = 4, dpi = 600)

            # methods_barplot_all %>% ggplot(., aes(x = Methods, y = AUROC, colour = Methods)) +
            #   geom_boxplot() +
            #   theme_bw() +
            #   # scale_x_discrete(labels = methods)+
            #   theme(legend.position = "none")

            P1 <- ggplot(methods_barplot_all, aes(x = Methods, y = AUROC, fill = Methods)) +
                stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
                geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white") +
                geom_jitter(aes(fill = Methods), width = 0.2, shape = 21, size = 2.5) +
                scale_fill_manual(values = mycol) +
                scale_color_manual(values = mycol) +
                scale_x_discrete(labels = methods) +
                theme_bw() +
                theme(legend.position = "bottom") +
                ylab("AUROC") +
                xlab("Methods")
            # P1
            df_res10 <- melt(evaluation_AUROC_dataset, id = "Dataset", variable.name = "Method", value.name = "AUROC")
            df_res10$Method <- factor(df_res10$Method,
                levels = methods
            )

            p2 <- ggplot(df_res10, aes(x = Dataset, y = AUROC, fill = Method)) +
                geom_bar(stat = "identity", position = position_dodge(), color = "black", width = 0.8) +
                scale_fill_manual(values = mycol) +
                # geom_errorbar(aes(ymin=AUROC - 0.1, ymax=AUROC + 0.1), position = position_dodge(.6), width=.2)+
                scale_x_discrete(labels = evaluation_AUROC$Dataset) +
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
            # ggsave(paste0("../scGRN-L0_output/output_Synthetic/dataset-cell/Methods-contrast-AUROC-", cell, "-", dataset, "-2.png"), width = 15, height = 4, dpi = 600)
        }
    }

    datasetPlot_AUROC_list <- list()
    for (d in 1:length(data_path)) {
        dataset <- data_path[d]
        if (dataset == "dyn-BF") title <- "BF"
        if (dataset == "dyn-CY") title <- "CY"
        if (dataset == "dyn-LI") title <- "LI"
        if (dataset == "dyn-LL") title <- "LL"
        evaluation_AUROC <- evaluation_AUROC_all[grep(dataset, evaluation_AUROC_all$Dataset), ]
        methods_barplot_all <- evaluation_AUROC %>%
            as.data.frame() %>%
            pivot_longer(
                cols = 2:c(ncol(evaluation_AUROC)),
                names_to = "Method",
                values_to = "AUROC"
            )

        methods <- names(evaluation_AUROC[, -1])
        methods_barplot_all$Method <- factor(methods_barplot_all$Method,
            levels = methods
        )

        p <- ggplot(
            methods_barplot_all,
            aes(x = Method, y = AUROC)
        ) +
            geom_boxplot(aes(fill = Method),
                width = 0.8
            ) +
            stat_compare_means(
                method = "wilcox.test",
                label = "p.signif",
                comparisons = my_comparisons,
                bracket.size = 0.6,
                tip.length = 0,
                vjust = 1.5,
                # color = "#6699cc",
                sizen = 4
            ) +
            scale_fill_manual(values = mycol) +
            # scale_color_manual(values = mycol) +
            scale_x_discrete(labels = methods) +
            labs(x = "Methods", y = "AUROC") +
            theme(legend.position = "bottom") +
            theme_bw() +
            ggtitle(title) +
            theme(
                axis.text.x = element_text(
                    angle = 45,
                    hjust = 1,
                    vjust = 1,
                    size = 10
                )
            )
        p
        datasetPlot_AUROC_list[[d]] <- p
        ggsave(paste0("../scGRN-L0_output/output_Synthetic/dataset/Methods-contrast-AUROC-", dataset, "-1.png"), width = 4, height = 5, dpi = 600)

        df_res10 <- melt(evaluation_AUROC, id = "Dataset", variable.name = "Method", value.name = "AUROC")
        df_res10$Method <- factor(df_res10$Method,
            levels = methods
        )

        p2 <- ggplot(df_res10, aes(x = Dataset, y = AUROC, fill = Method)) +
            geom_bar(stat = "identity", position = position_dodge(), color = "black", width = 0.8) +
            scale_fill_manual(values = mycol) +
            # geom_errorbar(aes(ymin=AUROC - 0.1, ymax=AUROC + 0.1), position = position_dodge(.6), width=.2)+
            scale_x_discrete(labels = evaluation_AUROC$Dataset) +
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
        # ggsave(paste0("../scGRN-L0_output/output_Synthetic/dataset/Methods-contrast-AUROC-", dataset, "-2.png"), width = 15, height = 4, dpi = 600)
    }
    datasetPlot_AUROC_list[[1]] + datasetPlot_AUROC_list[[2]] + datasetPlot_AUROC_list[[3]] + datasetPlot_AUROC_list[[4]] +
        plot_annotation(tag_levels = "A") +
        plot_layout(ncol = 4) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom") +
            theme(text = element_text(size = 12)) +
            theme(text = element_text(family = "Times New Roman"))
    ggsave(paste0("../scGRN-L0_output/output_Synthetic/Methods-contrast-AUROC-dataset.png"), width = 10, height = 4.5, dpi = 600)

    evaluation_AUROC_mean <- matrix(0, nrow = 5, ncol = 5)
    colnames(evaluation_AUROC_mean) <- cell_num
    rownames(evaluation_AUROC_mean) <- methods
    for (c in 1:length(cell_num)) {
        cell <- cell_num[c]
        evaluation_AUROC <- evaluation_AUROC_all[grep(cell, evaluation_AUROC_all$Dataset), ]
        for (m in 1:5) {
            evaluation_AUROC_mean[m, c] <- mean(evaluation_AUROC[, (m + 1)])
        }
        methods_barplot_all <- evaluation_AUROC %>%
            as.data.frame() %>%
            pivot_longer(
                cols = 2:c(ncol(evaluation_AUROC)),
                names_to = "Method",
                values_to = "AUROC"
            )

        methods <- names(evaluation_AUROC[, -1])
        methods_barplot_all$Method <- factor(methods_barplot_all$Method,
            levels = methods
        )
        p <- ggplot(
            methods_barplot_all,
            aes(x = Method, y = AUROC)
        ) +
            # geom_violin(aes(fill = Method),
            #     trim = FALSE
            # ) +
            geom_boxplot(aes(fill = Method),
                width = 0.8
            ) +
            stat_compare_means(
                method = "wilcox.test",
                label = "p.signif",
                comparisons = my_comparisons,
                bracket.size = 0.6,
                sizen = 4,
                color = "#6699cc"
            ) +
            scale_fill_manual(values = mycol) +
            # scale_color_manual(values = mycol) +
            scale_x_discrete(labels = methods) +
            labs(x = "Methods", y = "AUROC") +
            theme(legend.position = "bottom") +
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
        ggsave(paste0("../scGRN-L0_output/output_Synthetic/cell/Methods-contrast-AUROC-", cell, "-1.png"), width = 4, height = 4, dpi = 600)

        df_res10 <- melt(evaluation_AUROC, id = "Dataset", variable.name = "Method", value.name = "AUROC")
        df_res10$Method <- factor(df_res10$Method,
            levels = methods
        )

        p2 <- ggplot(df_res10, aes(x = Dataset, y = AUROC, fill = Method)) +
            geom_bar(stat = "identity", position = position_dodge(), color = "black", width = 0.8) +
            scale_fill_manual(values = mycol) +
            # geom_errorbar(aes(ymin=AUROC - 0.1, ymax=AUROC + 0.1), position = position_dodge(.6), width=.2)+
            scale_x_discrete(labels = evaluation_AUROC$Dataset) +
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
        # ggsave(paste0("../scGRN-L0_output/output_Synthetic/cell/Methods-contrast-AUROC-", cell, "-2.png"), width = 15, height = 4, dpi = 600)
    }

    mydata <- melt(evaluation_AUROC_mean, id = rownames(evaluation_AUROC_mean))
    colnames(mydata) <- c("Method", "CellNum", "AUROC")
    cellnum_con <- ggplot(data = mydata, aes(x = factor(CellNum), y = AUROC, group = Method, color = Method, shape = Method)) +
        geom_point(
            size = 2.5
        ) +
        geom_line() +
        scale_color_manual(values = mycol) +
        scale_fill_manual(values = mycol) +
        xlab("Cell Number") +
        ylab("AUROC") +
        theme_bw() +
        theme(
            text = element_text(family = "Times New Roman"),
            # legend.position = c(.75, .25),
            # legend.box.background = element_rect(color = "black"),
            legend.position = "right"
        )

    ggsave(paste0("../scGRN-L0_output/output_Synthetic/Methods-contrast-AUROC-cell.png"), width = 5, height = 2.5, dpi = 600)
    cellnum_con_list <- list()
    cellnum_con_list[[1]] <- cellnum_con
}

# # --------------------------------------------------
# output <- "../scGRN-L0_output/output_Curated/"
# evaluation_AUROC_all1 <- read.csv(paste0(output, "evaluation_AUROC.csv"))
output <- "../scGRN-L0_output/output_Synthetic/"
evaluation_AUROC_all2 <- read.csv(paste0(output, "evaluation_AUROC.csv"))
# evaluation_AUROC_all <- rbind.data.frame(evaluation_AUROC_all1, evaluation_AUROC_all2)
# head(evaluation_AUROC_all[1:3, 1:3])
# methods_barplot_all <- evaluation_AUROC %>%
#     as.data.frame() %>%
#     pivot_longer(
#         cols = 2:c(ncol(evaluation_AUROC)),
#         names_to = "Method",
#         values_to = "AUROC"
#     )

# methods_barplot_all$Method <- factor(methods_barplot_all$Method,
#     levels = methods
# )
# p <- ggplot(
#     methods_barplot_all,
#     aes(x = Method, y = AUROC)
# ) +
#     # geom_violin(aes(fill = Method),
#     #     trim = FALSE
#     # ) +
#     geom_boxplot(aes(fill = Method),
#         width = 0.8
#     ) +
#     stat_compare_means(
#         method = "wilcox.test",
#         label = "p.signif",
#         comparisons = my_comparisons,
#         bracket.size = 0.6,
#         sizen = 4,
#         color = "#6699cc"
#     ) +
#     scale_fill_manual(values = mycol) +
#     # scale_color_manual(values = mycol) +
#     scale_x_discrete(labels = methods) +
#     labs(x = "Methods", y = "AUROC") +
#     theme(legend.position = "bottom") +
#     theme_bw() +
#     theme(
#         axis.text.x = element_text(
#             angle = 45,
#             hjust = 1,
#             vjust = 1,
#             size = 10
#         )
#     )
# p
# ggsave(paste0("../scGRN-L0_output/Methods-contrast-AUROC-all-1.png"), width = 4, height = 4, dpi = 600)

# # methods_barplot_all %>% ggplot(., aes(x = Methods, y = AUROC, colour = Methods)) +
# #   geom_boxplot() +
# #   theme_bw() +
# #   # scale_x_discrete(labels = methods)+
# #   theme(legend.position = "none")

# P1 <- ggplot(methods_barplot_all, aes(x = Method, y = AUROC, fill = Method)) +
#     stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
#     geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white") +
#     geom_jitter(aes(fill = Method), width = 0.2, shape = 21, size = 2.5) +
#     scale_fill_manual(values = mycol) +
#     scale_color_manual(values = mycol) +
#     scale_x_discrete(labels = methods) +
#     theme_bw() +
#     theme(legend.position = "none") +
#     ylab("AUROC") +
#     xlab("Methods") +
#     theme(
#         axis.text.x = element_text(
#             angle = 45,
#             hjust = 1,
#             vjust = 1,
#             size = 10
#         )
#     )
# P1
# ggsave(paste0("../scGRN-L0_output/Methods-contrast-AUROC-all-2.png"), width = 4, height = 4, dpi = 600)

# # --------------------------------------------------
# output <- "../scGRN-L0_output/output_Curated/"
# evaluation_AUROC_all1 <- read.csv(paste0(output, "evaluation_AUROC.csv"))
# output <- "../scGRN-L0_output/output_Synthetic/"
# evaluation_AUROC_all2 <- read.csv(paste0(output, "evaluation_AUROC.csv"))
# evaluation_AUROC_all <- rbind.data.frame(evaluation_AUROC_all1, evaluation_AUROC_all2)

# head(evaluation_AUROC_all1[1:3, 1:3])
# methods_barplot_all <- evaluation_AUROC %>%
#     as.data.frame() %>%
#     pivot_longer(
#         cols = 2:c(ncol(evaluation_AUROC)),
#         names_to = "Method",
#         values_to = "AUROC"
#     )

# methods_barplot_all$Method <- factor(methods_barplot_all$Method,
#     levels = methods
# )
# p <- ggplot(
#     methods_barplot_all,
#     aes(x = Method, y = AUROC)
# ) +
#     # geom_violin(aes(fill = Method),
#     #     trim = FALSE
#     # ) +
#     geom_boxplot(aes(fill = Method),
#         width = 0.8
#     ) +
#     stat_compare_means(
#         method = "wilcox.test",
#         label = "p.signif",
#         comparisons = my_comparisons,
#         bracket.size = 0.6,
#         sizen = 4,
#         color = "#6699cc"
#     ) +
#     scale_fill_manual(values = mycol) +
#     # scale_color_manual(values = mycol) +
#     scale_x_discrete(labels = methods) +
#     labs(x = "Methods", y = "AUROC") +
#     theme(legend.position = "bottom") +
#     theme_bw() +
#     theme(
#         axis.text.x = element_text(
#             angle = 45,
#             hjust = 1,
#             vjust = 1,
#             size = 10
#         )
#     )
# p
# ggsave(paste0("../scGRN-L0_output/Methods-contrast-AUROC-Curated-1.png"), width = 4, height = 4, dpi = 600)

# # methods_barplot_all %>% ggplot(., aes(x = Methods, y = AUROC, colour = Methods)) +
# #   geom_boxplot() +
# #   theme_bw() +
# #   # scale_x_discrete(labels = methods)+
# #   theme(legend.position = "none")

# P1 <- ggplot(methods_barplot_all, aes(x = Method, y = AUROC, fill = Method)) +
#     stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
#     geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white") +
#     geom_jitter(aes(fill = Method), width = 0.2, shape = 21, size = 2.5) +
#     scale_fill_manual(values = mycol) +
#     scale_color_manual(values = mycol) +
#     scale_x_discrete(labels = methods) +
#     theme_bw() +
#     theme(legend.position = "none") +
#     ylab("AUROC") +
#     xlab("Methods") +
#     theme(
#         axis.text.x = element_text(
#             angle = 45,
#             hjust = 1,
#             vjust = 1,
#             size = 10
#         )
#     )
# P1
# ggsave(paste0("../scGRN-L0_output/Methods-contrast-AUROC-Curated-2.png"), width = 4, height = 4, dpi = 600)

head(evaluation_AUROC_all2[1:3, 1:3])
methods_barplot_all <- evaluation_AUROC %>%
    as.data.frame() %>%
    pivot_longer(
        cols = 2:c(ncol(evaluation_AUROC)),
        names_to = "Method",
        values_to = "AUROC"
    )

methods_barplot_all$Method <- factor(methods_barplot_all$Method,
    levels = methods
)
p <- ggplot(
    methods_barplot_all,
    aes(x = Method, y = AUROC)
) +
    # geom_violin(aes(fill = Method),
    #     trim = FALSE
    # ) +
    geom_boxplot(aes(fill = Method),
        width = 0.8
    ) +
    stat_compare_means(
        method = "wilcox.test",
        label = "p.signif",
        comparisons = my_comparisons,
        bracket.size = 0.6,
        sizen = 4,
        color = "#6699cc"
    ) +
    scale_fill_manual(values = mycol) +
    # scale_color_manual(values = mycol) +
    scale_x_discrete(labels = methods) +
    labs(x = "Methods", y = "AUROC") +
    theme(legend.position = "bottom") +
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
ggsave(paste0("../scGRN-L0_output/Methods-contrast-AUROC-Synthetic-1.png"), width = 4, height = 4, dpi = 600)

# methods_barplot_all %>% ggplot(., aes(x = Methods, y = AUROC, colour = Methods)) +
#   geom_boxplot() +
#   theme_bw() +
#   # scale_x_discrete(labels = methods)+
#   theme(legend.position = "none")

P1 <- ggplot(methods_barplot_all, aes(x = Method, y = AUROC, fill = Method)) +
    stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
    geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white") +
    geom_jitter(aes(fill = Method), width = 0.2, shape = 21, size = 2.5) +
    scale_fill_manual(values = mycol) +
    scale_color_manual(values = mycol) +
    scale_x_discrete(labels = methods) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab("AUROC") +
    xlab("Methods") +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1,
            size = 10
        )
    )
P1
ggsave(paste0("../scGRN-L0_output/Methods-contrast-AUROC-Synthetic-2.png"), width = 4, height = 4, dpi = 600)
