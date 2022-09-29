<<<<<<< HEAD
# TODO: Add comment
#
# Author: xiaoyun
###############################################################################
# This algorithm is a trajectory construction algorithm based on cell clusters.
# First, the cells are dimensionality reduction (PCA), and then the center of each cluster is identified, and the center points are connected to determine the preliminary trajectory
# then smooth the trajectory; finally project all points on the branch, and the distance between each point and the starting point is the pseudo time of the cell state.
###
# Using slingshot to infer cell trajectories and calculate pseudo-time
###
# input
#' @param scRNAseq.Exp：matrix，normalized count expression profile
#' @param clusterLabels：cell annotation
#' @param ordergene：character vector. feature genes, which used to construct the trajectory
#' @param RMmethod: character,Methods of dimensionality reduction, pca/tsne
#' @param start.cluster: character
#' @param end.cluster:   character
#' @param  stretch:default 2, numeric factor by which curves can be extrapolated beyond endpoints.
#'                  Must be a numeric value between 0 and 2.
#' @return
#' 1) PC：top 3 PCs of PCA；Branch：Which branch the cell belongs to；W.curve：weighs of cells in lineage；PseTime.curve：pseudotime of cells
#' 2) Curves：the coordinates of different lineages
#' 3) Centers：the coordinates of centers for all clusters
#' 4) Lineages：different lineages
slingshot_run <- function(scRNAseq.Exp,
                          clusterLabels,
                          ordergene,
                          RMmethod = "pca",
                          start.cluster = NULL,
                          end.cluster = NULL,
                          stretch = 2,
                          plot_output = FALSE) {
    library(Rtsne)
    library(slingshot)
    library(SingleCellExperiment)
    cat(format(Sys.time(), "[%b-%d,%H:%M:%S] "), "start to run......", "\n")
    # 2)Feature selection
    ordergene <- intersect(ordergene, rownames(scRNAseq.Exp))
    sim <- SingleCellExperiment(assays = List(counts = scRNAseq.Exp[ordergene, ]), colData = CellInfor_B) ##### revise
    sim <- sc_obj_filter(sim)
    # 3)Dimensionality reduction and cell trajectory construction
    if (RMmethod == "pca") {
        # pca <- prcomp(t(scRNAseq.Exp[ordergene, ]), scale. = T)
        #pca <- prcomp(t(scRNAseq.Exp[ordergene, ]), scale. = FALSE) ##### revise #for a count
        # pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE) ##### revise #for a seurat objective
        pca <- prcomp(t(log1p(assays(sim)$counts)), scale. = FALSE) ##### revise #for a seurat objective
        reducedDims(sim) <- SimpleList(PCA = pca$x[, 1:30]) ##### revise
        # top 30 PCs
        sds <- slingshot(sim, clusterLabels = clusterLabels, start.clus = start.cluster, end.clus = end.cluster, stretch = stretch, reducedDim = "PCA") ##### revise
        sds <- slingshot(sim, clusterLabels = "majorCluster", reducedDim = "UMAP")
        # sds <- slingshot(pca$x[,1:30], clusterLabels = clusterLabels, start.clus = start.cluster, end.clus=end.cluster, stretch=stretch)
    }

    if (RMmethod == "tsne") {
        tsne <- Rtsne(t(scRNAseq.Exp[ordergene, ]))
        tsne <- tsne$Y
        rownames(tsne) <- colnames(scRNAseq.Exp)
        sds <- slingshot(tsne, clusterLabels = clusterLabels, start.clus = start.cluster, end.clus = end.cluster, stretch = stretch)
    }

    if (plot_output == TRUE) {
        library(scater)
        library(scRNAseq)
        library(RColorBrewer)
        colors <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)
        plotcol <- colors[cut(sds$slingPseudotime_1, breaks = 100)]
        plotPCA(sds,
                # ncomponents = 4,
                colour_by = "majorCluster"
        )
        plotUMAP(sds,
                # ncomponents = 4,
                colour_by = "majorCluster"
        )
        png("Results/Slingshot1.png")
        plot(reducedDims(sds)$UMAP, col = plotcol, pch = 16, asp = 1,xlab="UMAP_1",ylab="UMAP_2")
        lines(SlingshotDataSet(sds), lwd = 2, col = "black")
        dev.off()
        png("Results/Slingshot2.png")
        plot(reducedDims(sds)$UMAP, col = brewer.pal(9, "Set1"), pch = 16, asp = 1 ,xlab="UMAP_1",ylab="UMAP_2")
        # plot(reducedDims(sds)$UMAP, col = brewer.pal(9, "Set1")[sds$majorCluster], pch = 16, asp = 1)
        lines(SlingshotDataSet(sds), lwd = 2, type = "lineages", col = "black")
        dev.off()
        # plotColData(sds, x = "sum", y="detected", colour_by="Patient")
        # plotColData(sds, x = "sum", y="detected", colour_by="majorCluster")
        # plotColData(sds, x = "sum", y="detected", colour_by="sampleType")
        # plotColData(sds, x = "sum", y="subsets_Mito_percent", 
        #             other_fields="majorCluster") + facet_wrap(~majorCluster)
        
        # plotHighestExprs(sds, exprs_values = "counts")
        p <- plotUMAP(sds,
                      # ncomponents = 4,
                      colour_by = "Patient"
        )+theme_bw()
        ggsave(p, filename = "Results/Patient.png")
        p <- plotUMAP(sds,
                      # ncomponents = 4,
                      colour_by = "majorCluster"
        )+theme_bw()
        ggsave(p, filename = "Results/Cluster.png")
        p <- plotUMAP(sds,
                # ncomponents = 4,
                colour_by = "sampleType"
        )+theme_bw()
        ggsave(p, filename = "Results/PatientType.png")
        # plotTSNE(sds,
        #          # ncomponents = 4,
        #          colour_by = "Patient"
        # )
        # plotTSNE(sds,
        #          # ncomponents = 4,
        #          colour_by = "majorCluster"
        # )
        # plotTSNE(sds,
        #          # ncomponents = 4,
        #          colour_by = "sampleType"
        # )
        
    }

    # top 3 PCs
    t.PC <- reducedDim(sds)[, 1:3]

    # Extract the branch to which the cell belongs
    t.Branch <- slingBranchID(sds)

    # Extract the weight of cells
    t.Weight <- slingCurveWeights(sds)
    colnames(t.Weight) <- paste("W", colnames(t.Weight), sep = ".")

    # extract the pseudotime
    traj <- slingPseudotime(sds)
    colnames(traj) <- paste("PseTime", colnames(traj), sep = ".")

    t.lineages <- slingLineages(sds)

    t1 <- slingCurves(sds)

    t.curve <- lapply(t1, function(x) {
        x$s[x$ord, 1:3]
    })

    # Determine the center  of each cluster to plot the main trajectory
    centers <- apply(t.PC, 2, function(pc) {
        tapply(pc, clusterLabels, mean)
    })
    t.centers <- lapply(t.lineages, function(x) centers[x, ])

    # 5)output
    df <- list(data = data.frame(t.PC, Branch = t.Branch, t.Weight, traj), Curves = t.curve, Centers = t.centers, Lineages = t.lineages)
    cat(format(Sys.time(), "[%b-%d,%H:%M:%S] "), "complete.", "\n")
    return(df)
}
=======
# TODO: Add comment
#
# Author: xiaoyun
###############################################################################
# This algorithm is a trajectory construction algorithm based on cell clusters.
# First, the cells are dimensionality reduction (PCA), and then the center of each cluster is identified, and the center points are connected to determine the preliminary trajectory
# then smooth the trajectory; finally project all points on the branch, and the distance between each point and the starting point is the pseudo time of the cell state.
###
# Using slingshot to infer cell trajectories and calculate pseudo-time
###
# input
#' @param scRNAseq.Exp：matrix，normalized count expression profile
#' @param clusterLabels：cell annotation
#' @param ordergene：character vector. feature genes, which used to construct the trajectory
#' @param RMmethod: character,Methods of dimensionality reduction, pca/tsne
#' @param start.cluster: character
#' @param end.cluster:   character
#' @param  stretch:default 2, numeric factor by which curves can be extrapolated beyond endpoints.
#'                  Must be a numeric value between 0 and 2.
#' @return
#' 1) PC：top 3 PCs of PCA；Branch：Which branch the cell belongs to；W.curve：weighs of cells in lineage；PseTime.curve：pseudotime of cells
#' 2) Curves：the coordinates of different lineages
#' 3) Centers：the coordinates of centers for all clusters
#' 4) Lineages：different lineages
slingshot_run <- function(scRNAseq.Exp,
                          clusterLabels,
                          ordergene,
                          RMmethod = "pca",
                          start.cluster = NULL,
                          end.cluster = NULL,
                          stretch = 2,
                          plot_output = FALSE) {
    library(Rtsne)
    library(slingshot)
    cat(format(Sys.time(), "[%b-%d,%H:%M:%S] "), "start to run......", "\n")
    # 2)Feature selection
    ordergene <- intersect(ordergene, rownames(scRNAseq.Exp))
    sim <- SingleCellExperiment(assays = List(counts = scRNAseq.Exp[ordergene, ]), colData = CellInfor_B) ##### revise

    # 3)Dimensionality reduction and cell trajectory construction
    if (RMmethod == "pca") {
        # pca <- prcomp(t(scRNAseq.Exp[ordergene, ]), scale. = T)
        pca <- prcomp(t(scRNAseq.Exp[ordergene, ]), scale. = FALSE) ##### revise #for a count
        # pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE) ##### revise #for a seurat objective
        reducedDims(sim) <- SimpleList(PCA = pca$x[, 1:30]) ##### revise
        # top 30 PCs
        sds <- slingshot(sim, clusterLabels = clusterLabels, start.clus = start.cluster, end.clus = end.cluster, stretch = stretch, reducedDim = "PCA") ##### revise
        # sds <- slingshot(pca$x[,1:30], clusterLabels = clusterLabels, start.clus = start.cluster, end.clus=end.cluster, stretch=stretch)
    }

    if (RMmethod == "tsne") {
        tsne <- Rtsne(t(scRNAseq.Exp[ordergene, ]))
        tsne <- tsne$Y
        rownames(tsne) <- colnames(scRNAseq.Exp)
        sds <- slingshot(tsne, clusterLabels = clusterLabels, start.clus = start.cluster, end.clus = end.cluster, stretch = stretch)
    }

    if (plot_output == TRUE) {
        library(RColorBrewer)
        sds@assays
        example_sce <- ZeiselBrainData()
        class(example_sce)
        library(scater)
        library(scRNAseq)
        example_sce <- runPCA(example_sce)
        example_sce <- runTSNE(example_sce, perplexity = 10)
        # uses the top 500 genes with the highest variances to compute the first PCs
        str(reducedDim(example_sce, "PCA"))

        plotPCA(example_sce,
            ncomponents = 4,
            colour_by = "Patient"
        )

        colors <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)
        plotcol <- colors[cut(sds$slingPseudotime_1, breaks = 100)]
        plot(reducedDims(sds)$PCA, col = plotcol, pch = 16, asp = 1)
        lines(SlingshotDataSet(sds), lwd = 2, col = "black")
        plot(reducedDims(sds)$PCA, col = brewer.pal(9, "Set1")[sds$GMM], pch = 16, asp = 1)
        lines(SlingshotDataSet(sds), lwd = 2, type = "lineages", col = "black")
    }

    # top 3 PCs
    t.PC <- reducedDim(sds)[, 1:3]

    # Extract the branch to which the cell belongs
    t.Branch <- slingBranchID(sds)

    # Extract the weight of cells
    t.Weight <- slingCurveWeights(sds)
    colnames(t.Weight) <- paste("W", colnames(t.Weight), sep = ".")

    # extract the pseudotime
    traj <- slingPseudotime(sds)
    colnames(traj) <- paste("PseTime", colnames(traj), sep = ".")

    t.lineages <- slingLineages(sds)

    t1 <- slingCurves(sds)

    t.curve <- lapply(t1, function(x) {
        x$s[x$ord, 1:3]
    })

    # Determine the center  of each cluster to plot the main trajectory
    centers <- apply(t.PC, 2, function(pc) {
        tapply(pc, clusterLabels, mean)
    })
    t.centers <- lapply(t.lineages, function(x) centers[x, ])

    # 5)output
    df <- list(data = data.frame(t.PC, Branch = t.Branch, t.Weight, traj), Curves = t.curve, Centers = t.centers, Lineages = t.lineages)
    cat(format(Sys.time(), "[%b-%d,%H:%M:%S] "), "complete.", "\n")
    return(df)
}
>>>>>>> origin/master
