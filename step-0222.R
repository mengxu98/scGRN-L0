

### load libraries
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)

nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 30
pHB_lower <- 0
pHB_upper <- 5

theme_set(theme_cowplot())

#color scheme
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  Plasma = "#0000CD",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E")

# Data loading and QC

### sample list
samples <- read_excel("/data/mengxu/data/PhilipBischoff2021/metadata/patients_metadata_raw.xlsx", range = cell_cols("A:A")) %>% .$sample_id

### import cellranger files from different data sets
for (i in seq_along(samples)){
  assign(paste0("scs_data", i), Read10X(data.dir = paste0("/data/mengxu/data/PhilipBischoff2021/cellranger/", samples[i], "/filtered_feature_bc_matrix")))
}

### create seurat objects from cellranger files
for (i in seq_along(samples)){
  assign(paste0("seu_obj", i), CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), project = samples[i], min.cells = 3))
}

### merge data sets
seu_obj <- merge(seu_obj1, y = c(seu_obj2, seu_obj3, seu_obj4, seu_obj5, seu_obj6, seu_obj7, seu_obj8, seu_obj9, seu_obj10, 
                                 seu_obj11, seu_obj12, seu_obj13, seu_obj14, seu_obj15, seu_obj16, seu_obj17, seu_obj18, 
                                 seu_obj19, seu_obj20), add.cell.ids = samples, project = "lung")

### calculate mitochondrial, hemoglobin and ribosomal gene counts
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^MT-", col.name = "pMT")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^HBA|^HBB", col.name = "pHB")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^RPS|^RPL", col.name = "pRP")

qcparams <- c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP")
for (i in seq_along(qcparams)){
  print(VlnPlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident", pt.size = 0))
}
for (i in seq_along(qcparams)){
  print(RidgePlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident"))
}

VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "pMT"), pt.size = 0, group.by = "orig.ident", ncol = 1, log = T)
### clear environment
remove(seu_obj1)
remove(seu_obj2)
remove(seu_obj3)
remove(seu_obj4)
remove(seu_obj5)
remove(seu_obj6)
remove(seu_obj7)
remove(seu_obj8)
remove(seu_obj9)
remove(seu_obj10)
remove(seu_obj11)
remove(seu_obj12)
remove(seu_obj13)
remove(seu_obj14)
remove(seu_obj15)
remove(seu_obj16)
remove(seu_obj17)
remove(seu_obj18)
remove(seu_obj19)
remove(seu_obj20)

remove(scs_data1)
remove(scs_data2)
remove(scs_data3)
remove(scs_data4)
remove(scs_data5)
remove(scs_data6)
remove(scs_data7)
remove(scs_data8)
remove(scs_data9)
remove(scs_data10)
remove(scs_data11)
remove(scs_data12)
remove(scs_data13)
remove(scs_data14)
remove(scs_data15)
remove(scs_data16)
remove(scs_data17)
remove(scs_data18)
remove(scs_data19)
remove(scs_data20)
gc()

seu_obj <- subset(seu_obj, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & pMT < pMT_upper & pHB < pHB_upper)
seu_obj <- SCTransform(seu_obj, verbose = T, vars.to.regress = c("nCount_RNA", "pMT"), conserve.memory = T)
seu_obj <- RunPCA(seu_obj)
seu_obj <- RunUMAP(seu_obj, dims = 1:15, verbose = T)
seu_obj <- FindNeighbors(seu_obj, dims = 1:15)
seu_obj <- FindClusters(seu_obj, resolution = 0.2)

# for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
#   seu_obj <- FindClusters(seu_obj, resolution = i)
#   print(DimPlot(seu_obj, reduction = "umap") + labs(title = paste0("resolution: ", i)))
# }
# 
# for (i in c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP")) {
#   print(FeaturePlot(seu_obj, features = i, coord.fixed = T, sort.cell = T))
# }

# Add metadata
metatable <- read_excel("/data/mengxu/data/PhilipBischoff2021/metadata/patients_metadata_raw.xlsx")

metadata <- FetchData(seu_obj, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata$sample_id <- metadata$orig.ident
metadata <- left_join(x = metadata, y = metatable, by = "sample_id")
rownames(metadata) <- metadata$cell_id

seu_obj <- AddMetaData(seu_obj, metadata = metadata)

p1 <- DimPlot(seu_obj,
              reduction = "umap",
              group.by = "orig.ident"
) +
  theme_bw()
p2 <- DimPlot(seu_obj,
              reduction = "umap",
              group.by = "patient_id"
) +
  theme_bw()
p3 <- DimPlot(seu_obj,
              reduction = "umap",
              group.by = "tissue_type"
) +
  theme_bw()
p4 <- DimPlot(seu_obj,
              reduction = "umap",
              group.by = "SCT_snn_res.0.2"
) +
  theme_bw()+
  NoLegend()
p1 + p2 + p3+ p4

library(harmony)
seu_obj_harmony <- RunHarmony(seu_obj,
                            group.by.vars = "orig.ident",
                            assay.use = "SCT",
                            lambda = 1, # [0.5-2] The more smaller lambda value, the bigger integration efforts.
                            max.iter.harmony = 20)

seu_obj_harmony <- RunPCA(seu_obj_harmony)
seu_obj_harmony <- RunUMAP(seu_obj_harmony, reduction = "harmony",dims = 1:15, verbose = T)
seu_obj_harmony <- FindNeighbors(seu_obj_harmony, dims = 1:15)
seu_obj_harmony <- FindClusters(seu_obj_harmony, resolution = 0.2)

p1 <- DimPlot(seu_obj_harmony,
              reduction = "umap",
              group.by = "orig.ident"
) +
  theme_bw()
p2 <- DimPlot(seu_obj_harmony,
              reduction = "umap",
              group.by = "patient_id"
) +
  theme_bw()
p3 <- DimPlot(seu_obj_harmony,
              reduction = "umap",
              group.by = "tissue_type"
) +
  theme_bw()
p4 <- DimPlot(seu_obj_harmony,
              reduction = "umap",
              group.by = "SCT_snn_res.0.2"
) +
  theme_bw()+
  NoLegend()
p1+p2+p3+p4

save(seu_obj, file = "../scGRN-L0_data/data/seu_obj.Rdata")
save(seu_obj_harmony, file = "../scGRN-L0_data/data/seu_obj_harmony.Rdata")
load("../scGRN-L0_data/data/seu_obj_harmony.Rdata")

# Main cell type annotation

# seu_obj_harmony <- annotation_celltype(seu_obj_harmony, method = "celltypist")
# 
# table(seu_obj_harmony$celltype)
# 
# DotPlot(seu_obj_harmony, features = mainmarkers, group.by = "SCT_snn_res.0.2") + 
#   coord_flip() + 
#   scale_color_viridis()

mainmarkers <- c("PECAM1", "VWF", "ACTA2", "JCHAIN", "MS4A1", "PTPRC", "CD68", "KIT", "EPCAM", "CDH1", "KRT7", "KRT19")

mainmarkers <- c(
  "CD79A", # B cell
  "CD19",
  # 2022-Cancer cell-Intratumoral plasma cells predict outcomes to PD-L1 blockade in non-small cell lung cancer
  # Follicular B cells
  "BANK1",
  "CD83",
  "CD69",
  "SELL",
  # "LINC00926",
  # "MARCH1",
  # "FCER2",
  # "GAPT",
  # "HVCN1",
  #Germinal center B cells
  "AICDA",
  "HMGA1",
  "RGS13",
  # "GCSAM",
  # "LRMP",
  # "AC023590.1",
  # "SUSD3",
  #Plasma cell
  "MZB1",
  "DERL3",
  "XBP1",
  "IGHG2",
  "IGHGP",
  "IGHA2",
  # "SDC1",
  # "DERL3",
  # "JSRP1",
  # "TNFRSF17",
  # "SLAMF7",
  # "IGLV3-1",
  # "IGLV6-57",
  # "IGKV4-1",
  # "IGKV1-12",
  # "IGLC7",
  # "IGLL5"
  # Other papers
  # "MSA41", # B(MS4A1+)
  "IGHG1", # naive B cells (MS4A1+IGHG1-)
  #plasma cells (MZB1+IGHG1+)
  #cycling plasma cells(MZB1+IGHG1+MKI67+TOP2A+)
  "MKI67",
  "TOP2A",
  #memory B cells(CD27+MS4A1+IGHG1+)
  "CD27",
  #Bn (TCL1A, IGHD 和 IL4R)
  "TCL1A",
  "IGHD",
  "IL4R",
  "TNFRSF13B",
  "JCHAIN",
  "IGHG3",
  #Bm (AIM2, TNFRSF13B 和 CD27)
  "AIM2"
  #浆细胞 (MZB1, IGHG3 和 JCHAIN)
)
DotPlot(seu_obj_harmony, features = mainmarkers, group.by = "SCT_snn_res.0.2") +theme_bw()+ 
  coord_flip() + 
  scale_color_viridis()

ggsave2(paste0("Fig5.FeaturePlot_mainmarkers.png"),
        path = paste0("Results/"),
        width = 20, height = 12, units = "cm"
)

for (i in seq_along(mainmarkers)) {
  if (mainmarkers[i] %in% rownames(seu_obj_harmony)) {
      FeaturePlot(seu_obj_harmony, features = mainmarkers[i], coord.fixed = T, order = T, cols = viridis(10))
      ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".png"),
              path = paste0("Results/marker/"),
              width = 10, height = 10, units = "cm"
      )
    
  }
}

DotPlot(seu_obj, features = mainmarkers, group.by = "SCT_snn_res.0.2") + 
  coord_flip() + 
  scale_color_viridis()

DimPlot(seu_obj_harmony, group.by = "SCT_snn_res.0.2", label = T, label.size = 5)
#ggsave2("DimPlot_all_clusters.pdf", path = "../results", width = 20, height = 20, units = "cm")

Idents(seu_obj_harmony) <- seu_obj_harmony$SCT_snn_res.0.2
annotation_curated_main <- read_excel("/data/mengxu/data/PhilipBischoff2021/curated_annotation/curated_annotation_main2.xlsx")
new_ids_main <- annotation_curated_main$main_cell_type
names(new_ids_main) <- levels(seu_obj_harmony)
seu_obj_harmony <- RenameIdents(seu_obj_harmony, new_ids_main)
seu_obj_harmony@meta.data$main_cell_type <- Idents(seu_obj_harmony)
DimPlot(seu_obj_harmony, label = T, label.size = 5)



obj_cells <- c("Follicular B cells", "Plasma", "Germinal center B cells")
# subset(x = pbmc, idents = c("CD4 T cells", "CD8 T cells"), invert = TRUE)
seu_obj_data_obj_cells <- list()
for (i in 1:length(obj_cells)) {
  obj_cell <- obj_cells[i]
  seu_obj_data_obj_cell <- seu_obj_harmony[, (seu_obj_harmony$main_cell_type == obj_cell)]
  seu_obj_data_obj_cells[[i]] <- seu_obj_data_obj_cell
  rm(seu_obj_data_obj_cell)
}
seu_obj_data <- merge(seu_obj_data_obj_cells[[1]], seu_obj_data_obj_cells[2:length(seu_obj_data_obj_cells)])
rm(seu_obj_data_obj_cells)
dim(seu_obj_data)
table(seu_obj_data$main_cell_type)













