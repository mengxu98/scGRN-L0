

source("step_function.R")
library(data.table)
library(tidyr)
library(Seurat)
# Python
# library(reticulate)
# pandas <- import("pandas")
# numpy <- import("numpy")
# scanpy <- import("scanpy")
# celltypist <- import("celltypist")



nFeature_lower <- 200
nFeature_upper <- 10000
nCount_lower <- 100
nCount_upper <- 150000
pMT_lower <- 0
pMT_upper <- 20
GSE131907 <- fread(
    file = "/data/mengxu/data/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
    sep = "\t",
    header = T,
    check.names = F
) %>% as.data.frame()
row.names(GSE131907) <- GSE131907$Index
head(GSE131907[1:3, 1:3])

GSE131907_annotation <- fread(
    file = "/data/mengxu/data/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    sep = "\t",
    header = T,
    check.names = F
) %>% as.data.frame()
row.names(GSE131907_annotation) <- GSE131907_annotation$Index
GSE131907_annotation_B <- GSE131907_annotation[which(GSE131907_annotation$Cell_type == "B lymphocytes"), ]

filter_samples <- c("mBrain", "mLN", "nLN", "PE")
for (i in 1:length(filter_samples)) {
    GSE131907_annotation_B <- GSE131907_annotation_B[-which(GSE131907_annotation_B$Sample_Origin == filter_samples[i]), ]
}
GSE131907_annotation_B <- GSE131907_annotation_B[-which(GSE131907_annotation_B$Cell_subtype == "Undetermined"), ]
table(GSE131907_annotation_B$Sample_Origin)
table(GSE131907_annotation_B$Cell_subtype)
dim(GSE131907_annotation_B)
GSE131907_B <- GSE131907[, rownames(GSE131907_annotation_B)]
dim(GSE131907_B)
CellInfor_B <- data.frame(
    UniqueCell_ID = GSE131907_annotation_B$Barcode,
    Patient = GSE131907_annotation_B$Sample,
    majorCluster = GSE131907_annotation_B$Cell_subtype,
    sampleType = GSE131907_annotation_B$Sample_Origin
)
rownames(CellInfor_B) <- CellInfor_B$UniqueCell_ID
<<<<<<< HEAD

save(CellInfor_B, GSE131907_B,file = "../scGRN-L0_data/DynamicGRNPipe_ExampleData/GSE131907_B.RData")
=======
scRNAseq.Exp <- as.matrix(scRNAseq.Exp)
>>>>>>> origin/master
# GSE131907_B_seu <- CreateSeuratObject(
#     counts = GSE131907_B,
#     project = "GSE131907",
#     min.features = 200,
#     min.cells = 3
# )
# GSE131907_B_seu$Platform <- "10X"
# ### QC
# GSE131907_B_seu <- PercentageFeatureSet(GSE131907_B_seu, pattern = "^MT-", col.name = "pMT")
# GSE131907_B_seu <- PercentageFeatureSet(GSE131907_B_seu, pattern = "^HBA|^HBB", col.name = "pHB")
# GSE131907_B_seu <- PercentageFeatureSet(GSE131907_B_seu, pattern = "^RPS|^RPL", col.name = "pRP")
# GSE131907_B_seu <- subset(GSE131907_B_seu,
#     subset =
#         nFeature_RNA > nFeature_lower &
#             nFeature_RNA < nFeature_upper &
#             nCount_RNA > nCount_lower &
#             nCount_RNA < nCount_upper &
#             pMT < pMT_upper
# )

# GSE131907_B_seu <- annotation_celltype(GSE131907_B_seu, method = "celltypist") # method = "celltypist" or "singleR"
# levels(GSE131907_B_seu$predicted_labels)
# table(GSE131907_B_seu$predicted_labels)

# ###
# neo_osi_metadata <- read.csv("/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/neo-osi_metadata.csv")
# S01_metacells <- read.csv("/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/S01_metacells.csv")
# # Dataset
# neo_osi_rawdata <- fread(
#     file = "/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/neo-osi_rawdata.csv",
#     sep = ",",
#     header = T,
#     check.names = F
# ) %>% as.data.frame(na.omit())
# # neo_osi_rawdata <- as.data.frame(na.omit(neo_osi_rawdata))
# rownames(neo_osi_rawdata) <- neo_osi_rawdata$gene
# ###
# S01_datafinal <- fread(
#     file = "/data/mengxu/data/AshleyMaynard2020(PRJNA591860)/Data_input/csv_files/S01_datafinal.csv",
#     sep = ",",
#     header = T,
#     check.names = F
# ) %>% as.data.frame(na.omit())
# # S01_datafinal <- as.data.frame(na.omit(S01_datafinal))
# rownames(S01_datafinal) <- S01_datafinal$V1
