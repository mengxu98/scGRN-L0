

library(Seurat)
library(tidyr)
# data <- readRDS("../scGRN-L0_data/seu_obj_data_B_samples.rds")
# meta <- data@meta.data

load("../scGRN-L0_data/seu_obj_data_B_samples.Rdata")

load(paste0("/data/mengxu/data/all/lung_seu.Rdata"))
table(seu_obj_data$celltype)
table(seu_obj_data$orig.ident)

obj_cells <- c("B_mature", "B_naive", "B_plasma", "B_plasmablast")
# subset(x = pbmc, idents = c("CD4 T cells", "CD8 T cells"), invert = TRUE)
seu_obj_data_obj_cells <- list()
for (i in 1:length(obj_cells)) {
  obj_cell <- obj_cells[i]
  seu_obj_data_obj_cell <- seu_obj_data[, (seu_obj_data$celltype == obj_cell)]
  seu_obj_data_obj_cells[[i]] <- seu_obj_data_obj_cell
  rm(seu_obj_data_obj_cell)
}
seu_obj_data <- merge(seu_obj_data_obj_cells[[1]], seu_obj_data_obj_cells[2:length(seu_obj_data_obj_cells)])
rm(seu_obj_data_obj_cells)
rm(samples_list)
dim(seu_obj_data)

if (F) {
  
  Cellratio <- prop.table(table(seu_obj_data$celltype, seu_obj_data$orig.ident), margin = 2)#计算各组样本不同细胞群比例
  Cellratio

  Cellratio <- as.data.frame(Cellratio)
  colourCount = length(unique(Cellratio$Var1))
  library(ggplot2)
  ggplot(Cellratio) + 
    geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
    theme_classic() +
    labs(x='Sample',y = 'Ratio')+
    coord_flip()+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
  
  table(seu_obj_data$orig.ident)#查看各组细胞数
  prop.table(table(Idents(seu_obj_data)))
  table(Idents(seu_obj_data), seu_obj_data$orig.ident)#各组不同细胞群细胞数 colnames(scRNA_harmony)
  Cellratio <- prop.table(table(seu_obj_data$celltype, seu_obj_data$orig.ident), margin = 2)#计算各组样本不同细胞群比例
  Cellratio <- data.frame(Cellratio)
  library(reshape2)
  cellper <- reshape2::dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
  rownames(cellper) <- cellper[,1]
  cellper <- cellper[,-1]
  
  ###添加分组信息
  sample <- seu_obj_data$orig.ident
  group <-seu_obj_data$stage
  samples <- data.frame(sample, group)#创建数据框
  
  # rownames(samples)=samples$sample
  
  # cellper$sample <- samples[rownames(cellper),'sample']#R添加列
  cellper$sample <- rownames(cellper)#R添加列
  # cellper$group <- samples[rownames(cellper),'group']#R添加列
  cellper$group <- ""
  for (i in 1:nrow(cellper)) {
    stage <- samples[which(samples$sample== cellper$sample[i])[1] , "group"]
    cellper$group[i] <- stage
  }
  table(seu_obj_data$celltype)
  ###作图展示
  
  sce_groups = c("B_mature","B_naive","B_plasma", "B_plasmablast")
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  pplist = list()
  for(group_ in sce_groups){
    cellper_  = cellper %>% select(one_of(c('sample','group',group_)))#选择一组数据
    colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
    cellper_$percent = as.numeric(cellper_$percent)#数值型数据
    cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                        lower = quantile(percent, 0.25),
                                                        mean = mean(percent),
                                                        median = median(percent))#上下分位数
    print(group_)
    print(cellper_$median)
    
    pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
      geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
      stat_summary(fun=mean, geom="point", color="grey60") +
      theme_cowplot() +
      theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
      labs(title = group_, y='Percentage') +
      geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  0.5)
    
    ###组间t检验分析
    labely = max(cellper_$percent)
    compare_means(percent ~ group,  data = cellper_)
    my_comparisons <- list( c("stage-1", "stage-normal"),
                            c("stage-2", "stage-normal"),
                            c("stage-3", "stage-normal"),
                            c("stage-4", "stage-normal") )
    pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
    pplist[[group_]] = pp1
  }
  
  library(cowplot)
  plot_grid(pplist[["B_mature"]],
            pplist[["B_naive"]],
            pplist[["B_plasma"]],
            pplist[["B_plasmablast"]])

}

meta <- seu_obj_data@meta.data

CellInfor_B <- data.frame(UniqueCell_ID = colnames(seu_obj_data),
                        Patient = meta$orig.ident,
                        majorCluster= meta$celltype,
                        sampleType=meta$stage)
table(CellInfor_B$majorCluster)
GSE131907_B <- as.matrix(seu_obj_data@assays$RNA@counts)
# Example
# The example data can be downloaded at https://www.jianguoyun.com/p/DeYtp_AQ1bXiCRjXvYwE
# ====0.input====
load(file = "DynamicGRNPipe_ExampleData/clusterSig.RData") # genes used to construct cell trajectories
# load(file = "../scGRN-L0_data/GSE131907_B.RData") # expression profile and cell annotation #from GSE99254

# ==== 1. Construction of cell state transformation trajectory (slingshot) ====
library(slingshot)
source("DynamicGRNPipe_1.slingshot_B.R")
source("DynamicGRNPipe_Function.R")
t.slingshot <- slingshot_run(scRNAseq.Exp = GSE131907_B,
  clusterLabels = CellInfor_B$majorCluster,
  ordergene = unlist(clusterSig),
  RMmethod = "pca",
  plot_output = TRUE,
  start.cluster = "B_naive" # cannot determine what type of cell types as the first type
)
CellInfor.trajectory <- cbind.data.frame(CellInfor_B, t.slingshot$data)
CD8TCellExp.trajectory <- GSE131907_B

# ====2.slinding windows based on pseudotime and anotation of cells====
source("DynamicGRNPipe_2.CellWindowing.R")
cells <- rownames(CellInfor.trajectory)[grep("1", CellInfor.trajectory$Branch)] # cells in lineage1
CellInfor.trajectory$UniqueCell_ID <- rownames(CellInfor.trajectory) ########
cellInfor <- CellInfor.trajectory[cells, c("UniqueCell_ID", "majorCluster", "PseTime.Lineage1")]
table(CellInfor_B$majorCluster)

# Calculate all intersections of pseudo-time density curves of cells in different states
C1_C2_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "B_naive", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "B_mature", "PseTime.Lineage1"],
  filename = "Results/C1_C2.png"
)
C2_C3_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "B_mature", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "B_plasma", "PseTime.Lineage1"],
  filename = "Results/C2_C3.png"
)
C3_C4_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "B_plasma", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "B_plasmablast", "PseTime.Lineage1"],
  filename = "Results/C3_C4.png"
)
C4_C5_cross <- densityintersection(
  a = cellInfor[cellInfor$majorCluster == "B_plasmablast", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "GC B cells in the DZ", "PseTime.Lineage1"],
  filename = "Results/C4_C5.png"
)

C5_C61_cross <- densityintersections(
  a = cellInfor[cellInfor$majorCluster == "B_mature", "PseTime.Lineage1"],
  b = cellInfor[cellInfor$majorCluster == "B_naive", "PseTime.Lineage1"],
  c = cellInfor[cellInfor$majorCluster == "B_plasma", "PseTime.Lineage1"],
  d = cellInfor[cellInfor$majorCluster == "B_plasmablast", "PseTime.Lineage1"],
  filename = "Results/all.png"
)
# Manually select the intersection point according to the density plot
binpoint <- c(0, 
              max(cellInfor$PseTime.Lineage1)/5, 
              max(cellInfor$PseTime.Lineage1)*2/5, 
              max(cellInfor$PseTime.Lineage1)*3/5, 
              max(cellInfor$PseTime.Lineage1)*4/5, 
              max(cellInfor$PseTime.Lineage1))

# Two consecutive bins of cells as a window
# extract cells for each window
t1 <- cut(cellInfor$PseTime.Lineage1, binpoint)
bincells <- split(cellInfor$UniqueCell_ID, t1)
W1 <- unlist(bincells[1:2])
W2 <- unlist(bincells[2:3])
W3 <- unlist(bincells[3:4])
W4 <- unlist(bincells[4:5])
Windows1 <- list(W1 = W1, W2 = W2, W3 = W3, W4 = W4)


# ====3.constructing dynamic networks (GENIE3)======
source("DynamicGRNPipe_3.constructionNetwork_L0.R")
load("DynamicGRNPipe_ExampleData/DynamicGene1.RData")
load("DynamicGRNPipe_ExampleData/dorothea.RData")
# Construct a control network and calculate the control weight of each edge
weightofWindows <- DynNet_RF(
  Windows = Windows1[1],
  CD8TCellExp.trajectory = CD8TCellExp.trajectory,
  DynamicGene = DynamicGene1, # set Background genes,which used to construct the network, such as highly variable genes, dynamic genes along trajectory
  allTFs = allTFs[100:200], # set regulators
  detectNum = 10, detectPro = 0.05, meanExp = 1 # Noise filtering threshold
)

# Filter the edges for each window
lineage1dynet <- DynNet_filter(
  Windows = Windows1,
  CD8TCellExp.trajectory = CD8TCellExp.trajectory,
  weightofWindows = weightofWindows,
  weightThr = 0.02,
  nsd = 2,
  positivecor = 0#,
  #confidence = NULL
)
lineage1dynet1 <- lineage1dynet[1:2]
length(lineage1dynet1[[4]])
# extract active edges for each window
Dynnet_active1 <- lapply(lineage1dynet, function(x) {
  rownames(x) <- paste0(x[, 1], "_", x[, 2])
  x <- x[x[, "spearmanCor"] > 0, ]
  return(x)
})
names(Dynnet_active1) <- paste0("W", 1:length(Dynnet_active1))

source("ground-truth.R")
source("framework_main.R")
ground_truth_T(weightofWindows[[1]], dorothea_regulon_human)
evaluationObject <- prepareEval("ground_pred.txt",
  paste0("ground_truth.tsv"),
  totalPredictionsAccepted = 100000
)

GENIE3_AUROC <- calcAUROC(evaluationObject)
GENIE3_AUPR <- calcAUPR(evaluationObject)
GENIE3_AUROC

weightofWindows_L0 <- DynNet_L0(
  Windows = Windows1[1],
  CD8TCellExp.trajectory = CD8TCellExp.trajectory,
  DynamicGene = DynamicGene1, # set Background genes,which used to construct the network, such as highly variable genes, dynamic genes along trajectory
  allTFs = allTFs, # set regulators
  detectNum = 10, detectPro = 0.05, meanExp = 1 # Noise filtering threshold
)

# Filter the edges for each window
lineage1dynet_L0 <- DynNet_filter(
  Windows = Windows1,
  CD8TCellExp.trajectory = CD8TCellExp.trajectory,
  weightofWindows = weightofWindows_L0,
  weightThr = 0.02,
  nsd = 2,
  positivecor = 0#,
  #confidence = NULL
)
lineage1dynet_L0[[4]]
lineage1dynet1_L0 <- lineage1dynet_L0[1:2]
length(lineage1dynet1_L0[[4]])
# extract active edges for each window
Dynnet_active1_L0 <- lapply(lineage1dynet_L0, function(x) {
  rownames(x) <- paste0(x[, 1], "_", x[, 2])
  x <- x[x[, "spearmanCor"] > 0, ]
  return(x)
})
names(Dynnet_active1_L0) <- paste0("W", 1:length(Dynnet_active1_L0))

ground_truth_T(weightofWindows_L0[[1]], dorothea_regulon_human)
evaluationObject_L0 <- prepareEval("ground_pred.txt",
                                   paste0("ground_truth.tsv"),
                                   totalPredictionsAccepted = 100000
)

L0_AUROC <- calcAUROC(evaluationObject_L0)
L0_AUPR <- calcAUPR(evaluationObject_L0)
L0_AUROC # 0.4601513 no parallel
L0_AUPR

GENIE3_AUROC
GENIE3_AUPR

# Filter the edges for each window
lineage1dynet_L0 <- DynNet_filter(
  Windows = Windows1,
  CD8TCellExp.trajectory = CD8TCellExp.trajectory,
  weightofWindows = weightofWindows_L0,
  weightThr = 0.02,
  nsd = 2,
  positivecor = 0 # ,
  # confidence = NULL
)

ground_truth_T(weightofWindows_L0[[3]], dorothea_regulon_human)
evaluationObject_L0 <- prepareEval("ground_pred.txt",
  paste0("ground_truth.tsv"),
  totalPredictionsAccepted = 100000
)

L0_AUROC <- calcAUROC(evaluationObject_L0)
L0_AUPR <- calcAUPR(evaluationObject_L0)

GENIE3_AUROC
GENIE3_AUPR
L0_AUROC
L0_AUPR