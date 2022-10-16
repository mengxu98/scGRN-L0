

if (T) {
  rm(list = ls())
  gc()
  library(tidyverse)
  library(patchwork)
  library(scuttle)
  library(Seurat)
  library(scran)
  library(celldex)
  source("Function.R")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
}

stage <- c("normal", "1", "2", "3", "4")
for (s in stage) {
  seu_obj_list_filter <- list()
  cell_num_information <- c()
  samples_delete <- c()
  if (file.exists(paste0("/data/mengxu/data/all/lung_stage-", s, "_list_raw.Rdata")) == T) {
    load(paste0("/data/mengxu/data/all/lung_stage-", s, "_list_raw.Rdata"))
    for (i in 1:length(seu_obj_list)) {
      seu_obj_data <- seu_obj_list[[i]]
      seu_obj_data$orig.ident <- as.factor(seu_obj_data$orig.ident)
      sample <- levels(seu_obj_data$orig.ident)[1]
      samples[i] <- sample
      message("[", Sys.time(), "] -----: Now sample ", samples[i], " data normalization!")
      if (length(samples) == length(seu_obj_list_filter)) {
        print("END")
      } else if (length(colnames(seu_obj_data)) > 100) {
        seu_obj_data <- normalize_data(seu_obj_data)
        seu_obj_data <- CellCycleScoring(seu_obj_data,
          s.features = s.genes,
          g2m.features = g2m.genes,
          set.ident = TRUE
        )
        seu_obj_data$CC.Difference <- seu_obj_data$S.Score - seu_obj_data$G2M.Score
        # 观察细胞周期基因的表达情况？是否会影响轨迹推断和拟时间分析呢？如果有影响，是否应该去除呢？
        # RidgePlot(seu_obj_data, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
        # RidgePlot(seu_obj_data, features = c(s.genes[1:25], g2m.genes[1:25]), ncol = 8)
        # RidgePlot(seu_obj_data, features = c(s.genes[26:length(s.genes)], g2m.genes[26:length(g2m.genes)]), ncol = 8)

        if (F) {
          seu_obj_data <- RunPCA(seu_obj_data,
            features = c(s.genes, g2m.genes)
          )
          DimPlot(seu_obj_data)
          seu_obj_data <- RunPCA(seu_obj_data,
            features = VariableFeatures(seu_obj_data),
            nfeatures.print = 10
          )
          DimPlot(seu_obj_data)
          if (T) {
            seu_obj_data <- ScaleData(seu_obj_data,
              vars.to.regress = c("S.Score", "G2M.Score"),
              features = rownames(seu_obj_data)
            )
          } else if (F) {
            seu_obj_data$CC.Difference <- seu_obj_data$S.Score - seu_obj_data$G2M.Score
            seu_obj_data <- ScaleData(seu_obj_data,
              vars.to.regress = "CC.Difference",
              features = rownames(seu_obj_data)
            )
          }
          seu_obj_data <- RunPCA(seu_obj_data,
            features = VariableFeatures(seu_obj_data),
            nfeatures.print = 10
          )
          DimPlot(seu_obj_data)
          seu_obj_data <- RunPCA(seu_obj_data,
            features = c(s.genes, g2m.genes)
          )
          DimPlot(seu_obj_data)
        }

        seu_obj_data <- ScaleData(seu_obj_data,
          # vars.to.regress = c("nCount_RNA", "pMT"), #针对每个样本处理是是否还需要消除测序深度和线粒体的影响呢？
          # vars.to.regress = "CC.Difference", #too slow
          features = rownames(seu_obj_data)
        )

        seu_obj_data <- FindVariableFeatures(seu_obj_data,
          # nfeatures=2000,
          selection.method = "vst"
        )

        seu_obj_data <- seu_obj_data %>% RunPCA(
          verbose = T,
          features = VariableFeatures(seu_obj_data),
          # ndims.print = 6:10,
          nfeatures.print = 10
        )

        DimHeatmap(seu_obj_data, dims = c(1:10))
        DimHeatmap(seu_obj_data, dims = 1:9, nfeatures = 10, cells = 500, balanced = TRUE)
        pc.num <- 1:pc_num(seu_obj_data)
        seu_obj_data <- seu_obj_data %>%
          RunUMAP(dims = pc.num)

        seu_obj_data <- FindNeighbors(seu_obj_data, dims = pc.num) %>% FindClusters(resolution = 0.5)
        seu_obj_data <- annotation_celltype(seu_obj_data, method = "celltypist") # method = "celltypist" or "SingleR"

        raw_cell_num <- length(colnames(seu_obj_data))
        seu_obj_data <- doublets_filter(seu_obj_data, doublet_rate = 0.05)
        remain_cell_num <- length(colnames(seu_obj_data))

        filtered_cell_num <- raw_cell_num - remain_cell_num
        cell_doublets_filtered <- c(
          samples[i],
          raw_cell_num,
          remain_cell_num,
          filtered_cell_num
        )
        cell_num_information <- rbind(
          cell_num_information,
          cell_doublets_filtered
        )
        # GetAssayData(seu_obj_data,slot="counts",assay="RNA") 原始表达矩阵
        # GetAssayData(seu_obj_data,slot="data",assay="RNA") 标准化之后的表达矩阵
        # GetAssayData(seu_obj_data,slot="scale.data",assay="RNA") 中心化之后的表达矩阵
        platform <- seu_obj_data$platform
        seu_obj_data <- CreateSeuratObject(
          counts = seu_obj_data@assays$RNA@counts,
          project = samples[i],
          min.features = 200,
          min.cells = 3
        )
        seu_obj_data$platform <- platform
        seu_obj_data <- RenameCells(seu_obj_data, new.names = gsub("-1", "_1", colnames(seu_obj_data)))
        seu_obj_data <- RenameCells(object = seu_obj_data, add.cell.id = samples[i])
        
        seu_obj_data <- NormalizeData(seu_obj_data)
        seu_obj_data <- FindVariableFeatures(seu_obj_data)
        seu_obj_data <- ScaleData(seu_obj_data)
        seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
        pc.num <- 1:pc_num(seu_obj_data)
        
        seu_obj_data <- RunUMAP(seu_obj_data, dims = pc.num)
        seu_obj_data <- FindNeighbors(seu_obj_data, dims = pc.num) %>% FindClusters(resolution = 0.5)
        seu_obj_data <- annotation_celltype(seu_obj_data, method = "celltypist") # method = "celltypist" or "SingleR"
        
        seu_obj_list_filter[[i]] <- seu_obj_data
        # saveRDS(seu_obj_data, paste0("/data/mengxu/data/all/samples/",samples[i], ".rds")) #Save every sample
      } else {
        message("[", Sys.time(), "] -----: The cells of sample:", samples[i], " that in 'lung_stage-", s, " _list_raw.Rdata' is fewer 100!")
        # samples <- samples[-i]
        # seu_obj_list <- seu_obj_list[-i]
        samples_delete_filtered <- c(
          i,
          samples[i],
          length(colnames(seu_obj_data))
        )
        samples_delete <- rbind(
          samples_delete,
          samples_delete_filtered
        )
      }
    }
  } else {
    message("[", Sys.time(), "] -----: Not found 'lung_stage-", s, "_list_raw.Rdata' file, please run step00 ", s, " code!")
    # source(paste0('step-00_ReadData_stage-',s,'.R'))
  }

  message("[", Sys.time(), "] -----: Save samples information were selected!")
  cell_num_information <- as.data.frame(cell_num_information)
  names(cell_num_information) <- c("Sample", "Num. of cells-raw", "Num. of cells-remian", "Num. of cells-filtered")
  # write.csv(cell_num_information, file = paste0("/data/mengxu/data/all/lung_stage-", s, "_cell_num_information.csv"), row.names = F)

  if (length(samples_delete) != 0) {
    message("[", Sys.time(), "] -----: Save samples' information with less than 100 cells!")
    samples_delete <- as.data.frame(samples_delete)
    names(samples_delete) <- c("NO.", "Sample", "Num. of cells-raw")
    # write.csv(samples_delete, file = paste0("/data/mengxu/data/all/lung_stage-", s, "_samples_delete.csv"), row.names = F)
    message("[", Sys.time(), "] -----: Delete samples with less than 100 cells!")
    for (j in 1:length(samples_delete$NO.)) {
      print(j)
      n <- as.numeric(samples_delete$NO.[j])
      print(n)
      n <- n + 1 - j
      samples <- samples[-n]
      seu_obj_list_filter <- seu_obj_list_filter[-n]
    }
  }
  save(seu_obj_list_filter, samples, file = paste0("/data/mengxu/data/all/lung_stage-", s, "_list_filter.Rdata"))
}
