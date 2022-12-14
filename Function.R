

# Color--------------------------------------------------
colP <- c(
  "#3366cc", "#66ff66", "#003366", "#66ffcc", "#ffffcc", "#66ffff",
  "#ffcc00", "#66cc66", "#ffcc66", "#66cccc", "#ffcccc", "#66ccff",
  "#ff9900", "#669966", "#ff9966", "#6699cc", "#ff99cc", "#6699ff",
  "#ff6600", "#666666", "#ff6666", "#6666cc", "#ff66cc", "#6666ff",
  "#ff3300", "#663366", "#ff3366", "#6633cc", "#ff33cc", "#6633ff",
  "#ff0000", "#660066", "#ff0066", "#6600cc", "#ff00cc", "#6600ff",
  "#000099", "#33ff66", "#ccff66", "#66ff99", "#ccffcc", "#33ffff",
  "#ffccff", "#33cc66", "#cccc66", "#66cc99", "#cccccc", "#33ccff",
  "#ff99ff", "#339966", "#cc9966", "#669999", "#cc99cc", "#3399ff",
  "#ff66ff", "#336666", "#cc6666", "#666699", "#cc66cc", "#3366ff",
  "#ff33ff", "#333366", "#cc3366", "#663399", "#cc33cc", "#3333ff",
  "#ff00ff", "#330066", "#cc0066", "#660099", "#cc00cc", "#3300ff",
  "#99ff00", "#00ff66", "#99ff66", "#00ffcc", "#99ffcc", "#00ffff",
  "#99cc00", "#00cc66", "#99cc66", "#00cccc", "#99cccc", "#00ccff",
  "#999900", "#009966", "#999966", "#0099cc", "#9999cc", "#0099ff",
  "#996600", "#006666", "#996666", "#0066cc", "#9966cc", "#0066ff",
  "#993300", "#ffff00", "#993366", "#0033cc", "#9933cc", "#0033ff",
  "#990000", "#000066", "#990066", "#0000cc", "#9900cc", "#0000ff",
  "#66ff00", "#ffff99", "#33ffcc", "#ffff33", "#ccff00", "#66ff33",
  "#66cc00", "#ffcc99", "#33cccc", "#ffcc33", "#cccc00", "#66cc33",
  "#669900", "#ff9999", "#3399cc", "#ff9933", "#cc9900", "#669933",
  "#666600", "#ff6699", "#ffff66", "#ff6633", "#cc6600", "#666633",
  "#663300", "#ff3399", "#3333cc", "#ff3333", "#cc3300", "#663333",
  "#660000", "#ff0099", "#3300cc", "#ff0033", "#cc0000", "#660033",
  "#33ff00", "#ccff99", "#33ff99", "#ccff33", "#ccffff", "#33ff33",
  "#33cc00", "#cccc99", "#33cc99", "#cccc33", "#ccccff", "#33cc33",
  "#339900", "#cc9999", "#339999", "#cc9933", "#cc99ff", "#339933",
  "#336600", "#cc6699", "#336699", "#cc6633", "#cc66ff", "#336633",
  "#333300", "#cc3399", "#333399", "#cc3333", "#cc33ff", "#333333",
  "#330000", "#cc0099", "#330099", "#cc0033", "#cc00ff", "#330033",
  "#00ff00", "#99ff99", "#00ff99", "#99ff33", "#99ffff", "#00ff33",
  "#00cc00", "#99cc99", "#00cc99", "#99cc33", "#99ccff", "#00cc33",
  "#009900", "#999999", "#009999", "#999933", "#9999ff", "#009933",
  "#006600", "#996699", "#006699", "#996633", "#9966ff", "#006633",
  "#003300", "#993399", "#003399", "#993333", "#9933ff", "#003333",
  "#000000", "#990099", "#ffffff", "#990033", "#9900ff", "#000033"
)

# Packages check, download and library --------------------------------------------------
#' Title
#'
#' @param packages 
#'
#' @return
#' @export
#'
#' @examples
package.check <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      if (!requireNamespace("dplyr", quietly = TRUE)) {
        install.packages("dplyr")
      }
      library(dplyr)
      if (!requireNamespace("rvest", quietly = TRUE)) {
        install.packages("rvest")
      }
      library(rvest)
      message("[", Sys.time(), "] -----: No package: ", package, " in R environment!")
      CRANpackages <- available.packages() %>%
        as.data.frame() %>%
        select(Package) %>%
        mutate(source = "CRAN")
      url <- "https://www.bioconductor.org/packages/release/bioc/"
      biocPackages <- url %>%
        read_html() %>%
        html_table() %>%
        .[[1]] %>%
        select(Package) %>%
        mutate(source = "BioConductor")
      if (package %in% CRANpackages$Package) {
        message("[", Sys.time(), "] -----: Now install package: ", package, " from CRAN!")
        install.packages(package)
        library(package, character.only = T)
      } else if (package %in% biocPackages$Package) {
        message("[", Sys.time(), "] -----: Now install package: ", package, " from BioConductor!")
        BiocManager::install(package)
        library(package, character.only = T)
      } else { # Bug
        if (!requireNamespace("githubinstall", quietly = TRUE)) {
          install.packages("githubinstall")
        }
        library(githubinstall)
        # githubinstall(package)
        gh_suggest(package)
      }
    } else {
      library(package, character.only = T)
    }
  }
}

# Auto compute PC value --------------------------------------------------
#' Title
#'
#' @param sce 
#'
#' @return
#' @export
#'
#' @examples
pc_num <- function(sce) {
  # Judgment criteria:
  # 1. The cumulative contribution of principal components is greater than 90%
  # 2. The contribution of PC itself to the other party's difference is less than 5%
  # 3. The difference between two consecutive PCs is less than 0.1%
  message("[", Sys.time(), "]", " --- Now compute 'pc' number!")
  # Determine percent of variation associated with each PC
  pct <- sce[["pca"]]@stdev / sum(sce[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC, last point where change of variation is more than 0.1%.
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  message("[", Sys.time(), "]", " --- PC num should set as 1:", pcs)
  # Create a dataframe with values
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
  # Elbow plot to visualize
  p <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
    geom_text() +
    geom_vline(xintercept = cumu[pcs], color = "blue") + # 90
    geom_hline(yintercept = min(pct[pct > pct[pcs]]), color = "blue") + # pct > 5
    theme_bw()
  print(p)
  return(pcs)
}

# ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#' Title
#'
#' @param counts 
#'
#' @return
#' @export
#'
#' @examples
FQnorm <- function(counts) {
  rk <- apply(counts, 2, rank, ties.method = "min")
  counts.sort <- apply(counts, 2, sort)
  refdist <- apply(counts.sort, 1, median)
  norm <- apply(rk, 2, function(r) {
    refdist[r]
  })
  rownames(norm) <- rownames(counts)
  return(norm)
}

# Normalization for sample --------------------------------------------------
#' Title
#'
#' @param seu_obj 
#' @param platform 
#'
#' @return
#' @export
#'
#' @examples
# This function only support 3 platforms (10X, IndropSeq and SmartSeq2) of scRNA-seq dataset, please pay attention to select one of these.
normalize_data <- function(seu_obj, platform = NULL) {
  if (length(table(seu_obj$platform)) == 1) {
    if (is.null(platform)) {
      tryCatch(
        {
          is.null(seu_obj$platform)
        },
        # warning = function(w) {
        # },
        error = function(e) {
          platform <- "10X"
          seu_obj$platform <- platform
          message("[", Sys.time(), "] -----: There is no platform information, the default setting is \"10X\"!")
        },
        finally = {
          platform <- seu_obj$platform[1]
          project.name <- seu_obj@project.name
        }
      )
    }
    if (platform == "10X") {
      message("[", Sys.time(), "] -----: Platform 10x Genomics")
      seu_obj <- NormalizeData(object = seu_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
      seu_obj <- CreateSeuratObject(
        counts = seu_obj[["RNA"]]@data,
        project = project.name,
        min.features = 200,
        min.cells = 3
      )
    } else if (platform == "IndropSeq") {
      message("[", Sys.time(), "] -----: Platform IndropSeq")
      seu_obj <- NormalizeData(object = seu_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
      seu_obj <- CreateSeuratObject(
        counts = seu_obj[["RNA"]]@data,
        project = project.name,
        min.features = 200,
        min.cells = 3
      )
    } else if (platform == "SmartSeq2") {
      message("[", Sys.time(), "] -----: Platform Smart-seq2")
      sce_obj <- as.SingleCellExperiment(seu_obj)
      # QC
      qcstats <- perCellQCMetrics(sce_obj)
      qcfilter <- quickPerCellQC(qcstats)
      sce_obj <- sce_obj[, !qcfilter$discard]
      clusters <- quickCluster(sce_obj)
      # Normalization
      sce_obj <- computeSumFactors(sce_obj, clusters = clusters)
      summary(sizeFactors(sce_obj))
      sce_obj <- logNormCounts(sce_obj)
      seu_obj <- CreateSeuratObject(
        counts = sce_obj@assays@data$logcounts,
        project = project.name,
        min.features = 200,
        min.cells = 3
      )
    }
    seu_obj$platform <- platform
    return(seu_obj)
  } else {
    message("[", Sys.time(), "] -----: There are multiple platform information, please consider removing batch effect!")
  }
}

# Seurat object to 10x files --------------------------------------------------
# seu_obj <- ScaleData(seu_obj)  #seu_obj@assays$RNA@data
# seu_obj <- SCTransform(seu_obj)
# # switch assay
# DefaultAssay(object = seu_obj) <- "SCT"
# # Seurat object as 10x files, and the Seurat must be construst from 10X files, or will error: path exist
# package.check("DropletUtils")
# write10xCounts(x = seu_obj@assays$RNA@counts, path = '10x', version="3")

# Doublets --------------------------------------------------
#' Title
#'
#' @param seu_obj 
#' @param doublet_rate 
#' @param plot 
#' @param filename 
#'
#' @return
#' @export
#'
#' @examples
# Doublets filter --------------------------------------------------
doublets_filter <- function(seu_obj, doublet_rate = 0.039, plot = FALSE, filename = NULL) {
  package.check(c("scDblFinder", "DoubletFinder")) # devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
  pc.num <- 1:pc_num(seu_obj)
  seu_obj_raw <- seu_obj
  # Seek optimal pK value
  sweep.res.list <- paramSweep_v3(seu_obj, PCs = pc.num, sct = T)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>%
    as.character() %>%
    as.numeric()
  # doublet_rate = ncol(seu_obj)*8*1e-6 #The doublets rate is 3.9% of 5000 cells. For every 1000 additional cells, the double cell ratio increases by 8%
  homotypic.prop <- modelHomotypic(seu_obj$celltype)
  nExp_poi <- round(doublet_rate * ncol(seu_obj))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  seu_obj <- doubletFinder_v3(seu_obj, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
  seu_obj$doubFind.class <- seu_obj@meta.data %>% select(contains("DF.classifications"))
  seu_obj$doubFind.score <- seu_obj@meta.data %>% select(contains("pANN"))
  # table(seu_obj$doubFind.class)
  if (plot) {
    p <- DimPlot(seu_obj,
      reduction = "umap",
      group.by = "doubFind.class"
    ) +
      # theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
      theme_bw()
    # print(p)
    p
    if (is.null(filename)) {
      filename <- "Doublets.pdf"
    }
    ggsave(p, filename = filename)
  }

  sce_obj <- as.SingleCellExperiment(seu_obj)
  sce_obj <- scDblFinder(sce_obj, dbr = 0.1)
  # plotDoubletMap(sce_obj)
  # table(truth = sce_obj$scDblFinder.class, call = sce_obj$scDblFinder.class)
  # table(sce_obj$scDblFinder.class)
  sce_obj <- as.Seurat(sce_obj)
  intersection_doublet <- intersect(
    colnames(seu_obj[, which(seu_obj@meta.data$doubFind.class == "Doublet")]), # Singlet
    colnames(sce_obj[, which(sce_obj@meta.data$scDblFinder.class == "doublet")])
  )
  if (length(intersection_doublet) == 0) {
    seu_obj <- seu_obj_raw[, colnames(seu_obj_raw[, which(seu_obj@meta.data$doubFind.class == "Singlet")])]
  } else {
    seu_obj <- seu_obj_raw[, setdiff(colnames(seu_obj_raw), intersection_doublet)]
  }
  raw_cell_num <- length(colnames(seu_obj_raw))
  filter_cell_num <- length(colnames(seu_obj))
  cell_num_filtered <- raw_cell_num - filter_cell_num
  message("[", Sys.time(), "] -----: ", cell_num_filtered, " cells had filtered as doublets!")
  return(seu_obj)
}

# Annotation --------------------------------------------------
# https://github.com/Teichlab/celltypist
# celltypist$models$download_models(force_update = T) # First run needs to download the trained model data.
# source_python('celltypist_model_download.py')
# rownames(seu_obj@assays$RNA@counts)
#' Title
#'
#' @param seu_obj 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
annotation_celltype <- function(seu_obj, method = "celltypist") {
  if (method == "celltypist") {
    package.check("reticulate") # For Python packages
    pandas <- import("pandas")
    numpy <- import("numpy")
    scanpy <- import("scanpy")
    celltypist <- import("celltypist")
    tryCatch(
      {
        matrix <- t(as.matrix(seu_obj[["RNA"]]@counts))
        # matrix <- t(as.matrix(seu_obj@assays$SCT@counts))
      },
      # warning = function(w) {
      #   source("Function-as_matrix_cpp.R")
      #   matrix <- t(as_matrix(seu_obj[["RNA"]]@counts))
      # },
      error = function(e) {
        source("Function-as_matrix_cpp.R")
        matrix <- t(as_matrix(seu_obj[["RNA"]]@counts))
      }
    )
    message("[", Sys.time(), "] -----: Run 'celltypist'!")
    adata <- scanpy$AnnData(
      X = numpy$array(matrix),
      obs = pandas$DataFrame(seu_obj@meta.data),
      var = pandas$DataFrame(data.frame(
        gene = rownames(seu_obj[["RNA"]]@counts),
        row.names = rownames(seu_obj[["RNA"]]@counts)
      ))
    )
    model <- celltypist$models$Model$load(model = "Cells_Lung_Airway.pkl")
    model$cell_types
    scanpy$pp$normalize_total(adata, target_sum = 1e4)
    scanpy$pp$log1p(adata)
    predictions <- celltypist$annotate(adata, model = "Cells_Lung_Airway.pkl", majority_voting = FALSE)
    seu_obj <- AddMetaData(seu_obj, predictions$predicted_labels)
    seu_obj$celltype <- seu_obj$predicted_labels
    return(seu_obj)
  } else if (method == "SingleR") {
    package.check("SingleR")
    # hpca.se <- HumanPrimaryCellAtlasData() # Bug
    load("/data/mengxu/data/SingleR_data/HumanPrimaryCellAtlas_hpca.se_human.RData")
    # load("/data/mengxu/data/SingleR_data/BlueprintEncode_bpe.se_human.RData")
    message("[", Sys.time(), "] -----: Run 'singleR'!")
    matrix <- seu_obj@assays$RNA@counts
    singler_obj <- SingleR(test = matrix, ref = hpca.se, labels = hpca.se$label.main)
    seu_obj$celltype <- singler_obj$labels
    return(seu_obj)
  } else if (!(method %in% c("celltypist", "SingleR"))) {
    message("[", Sys.time(), "] -----: Please choose one of methods: 'celltypist' or 'SingleR!'")
  }
}

# Merge multiple Seurat objects --------------------------------------------------
#' Title
#'
#' @param seu_obj_list 
#' @param samples 
#' @param stage 
#'
#' @return
#' @export
#'
#' @examples
merge_seu_obj <- function(seu_obj_list, samples, stage) {
  package.check("Seurat")
  seu_obj <- merge(seu_obj_list[[1]],
    y = c(
      seu_obj_list[2:length(seu_obj_list)]
    ),
    add.cell.ids = samples,
    project = paste0("NSCLC-stage-", stage)
  )
  seu_obj$stage <- paste0("stage-", stage)
  return(seu_obj)
}

# Plot function --------------------------------------------------
nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 25
pHB_lower <- 0
pHB_upper <- 5

#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
qc_std_plot_helper <- function(x) {
  x +
    scale_color_viridis() +
    geom_point(size = 0.01, alpha = 0.3)
}

#' Title
#'
#' @param seu_obj 
#'
#' @return
#' @export
#'
#' @examples
qc_std_plot_nf <- function(seu_obj) {
  qc_data <- as_tibble(FetchData(seu_obj, c("nCount_RNA", "nFeature_RNA", "pMT", "pHB", "pRP")))
  plot_grid(
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pMT))),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pHB))),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pRP))),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pMT, color = nFeature_RNA))),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pHB, color = nFeature_RNA))),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pRP, color = nFeature_RNA))),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pMT, color = nCount_RNA))),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pHB, color = nCount_RNA))),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pRP, color = nCount_RNA))),
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nCount_RNA))),
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nFeature_RNA))),
    ggplot(gather(qc_data, key, value), aes(key, value)) +
      geom_violin() +
      facet_wrap(~key, scales = "free", ncol = 5),
    ncol = 3, align = "hv"
  )
}

#' Title
#'
#' @param seu_obj 
#'
#' @return
#' @export
#'
#' @examples
qc_std_plot <- function(seu_obj) {
  qc_data <- as_tibble(FetchData(seu_obj, c("nCount_RNA", "nFeature_RNA", "pMT", "pHB", "pRP")))
  plot_grid(
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pMT))) +
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pHB))) +
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pRP))) +
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pMT, color = nFeature_RNA))) +
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pHB, color = nFeature_RNA))) +
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pRP, color = nFeature_RNA))) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pMT, color = nCount_RNA))) +
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pHB, color = nCount_RNA))) +
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pRP, color = nCount_RNA))) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nCount_RNA))) +
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nFeature_RNA))) +
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),
    ggplot(gather(qc_data, key, value), aes(key, value)) +
      geom_violin() +
      facet_wrap(~key, scales = "free", ncol = 5),
    ncol = 3, align = "hv"
  )
}
