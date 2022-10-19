

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

samples_inf <- as.data.frame(table(seu_obj_data$orig.ident))
samples <- c()
for (i in 1:nrow(samples_inf)) {
  if (samples_inf$Freq[i]>100 & samples_inf$Freq[i]<2000) {
    sample <- as.character(samples_inf$Var1[i])
  }
  samples[i] <- sample
}
samples <- as.character(samples)
seu_obj_data_obj_cells <- list()
for (i in 1:length(samples)) {
  obj_cell <- samples[i]
  seu_obj_data_obj_cell <- seu_obj_data[, (seu_obj_data$orig.ident == obj_cell)]
  seu_obj_data_obj_cells[[i]] <- seu_obj_data_obj_cell
  rm(seu_obj_data_obj_cell)
}
seu_obj_data <- merge(seu_obj_data_obj_cells[[1]], seu_obj_data_obj_cells[2:length(seu_obj_data_obj_cells)])
rm(seu_obj_data_obj_cells)
dim(seu_obj_data)

seu_obj_data <- SCTransform(seu_obj_data,
                            method = "glmGamPoi")
seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
pc.num <- 1:pc_num(seu_obj_data)
seu_obj_data <- RunUMAP(seu_obj_data, dims = pc.num)
seu_obj_data <- FindNeighbors(seu_obj_data, dims = pc.num) %>% FindClusters(resolution = 1)
p1 <- DimPlot(seu_obj_data,
              reduction = "umap",
              group.by = "orig.ident"
) +
  theme_bw() +
  NoLegend()
p2 <- DimPlot(seu_obj_data,
              reduction = "umap",
              group.by = "stage"
) +
  theme_bw()
p3 <- DimPlot(seu_obj_data,
              reduction = "umap",
              group.by = "platform"
) +
  theme_bw()
p1 + p2 + p3
ggsave2("Fig1.raw_umap.png",
        path = paste0("Results/"),
        width = 30, height = 9, units = "cm"
)

scRNA_harmony <- RunHarmony(seu_obj_data,
                            group.by.vars = "stage",
                            lambda = 1, # [0.5-2] The more smaller lambda value, the bigger integration efforts.
                            max.iter.harmony = 20)
scRNA_harmony <- RunHarmony(scRNA_harmony,
                            group.by.vars = "orig.ident",
                            lambda = 1, # [0.5-2] The more smaller lambda value, the bigger integration efforts.
                            max.iter.harmony = 20)
pc.num <- 1:pc_num(scRNA_harmony)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num) %>%
  FindNeighbors(dims = pc.num) %>%
  FindClusters(resolution = 1)

p1 <- DimPlot(scRNA_harmony,
              reduction = "umap",
              group.by = "orig.ident"
) + theme_bw() + NoLegend()

p2 <- DimPlot(scRNA_harmony,
              reduction = "umap",
              group.by = "stage"
) + theme_bw()

p3 <- DimPlot(scRNA_harmony,
              reduction = "umap",
              group.by = "platform"
) + theme_bw() + NoLegend()
p1 + p2 + p3
ggsave2("Fig3.harmony_umap.png",
        path = paste0("Results/"),
        width = 30, height = 9, units = "cm"
)
DimPlot(scRNA_harmony,
        reduction = "umap",
        group.by = "celltype"
) + theme_bw()
scRNA_harmony <- annotation_celltype(scRNA_harmony, method = "celltypist")
DimPlot(scRNA_harmony,
        reduction = "umap",
        group.by = "celltype"
) + theme_bw()
table(scRNA_harmony$celltype)

obj_cells <- c("B_mature", "B_naive", "B_plasma", "B_plasmablast")
# subset(x = pbmc, idents = c("CD4 T cells", "CD8 T cells"), invert = TRUE)
seu_obj_data_obj_cells <- list()
for (i in 1:length(obj_cells)) {
  obj_cell <- obj_cells[i]
  seu_obj_data_obj_cell <- scRNA_harmony[, (scRNA_harmony$celltype == obj_cell)]
  seu_obj_data_obj_cells[[i]] <- seu_obj_data_obj_cell
  rm(seu_obj_data_obj_cell)
}
seu_obj_data <- merge(seu_obj_data_obj_cells[[1]], seu_obj_data_obj_cells[2:length(seu_obj_data_obj_cells)])
rm(seu_obj_data_obj_cells)
dim(seu_obj_data)

seu_obj_data <- SCTransform(seu_obj_data,
                            method = "glmGamPoi")
seu_obj_data <- RunPCA(seu_obj_data, verbose = T)
pc.num <- 1:pc_num(seu_obj_data)
seu_obj_data <- RunUMAP(seu_obj_data, dims = pc.num)
seu_obj_data <- FindNeighbors(seu_obj_data, dims = pc.num) %>% FindClusters(resolution = 1)
p1 <- DimPlot(seu_obj_data,
              reduction = "umap",
              group.by = "orig.ident"
) +
  theme_bw() +
  NoLegend()
p2 <- DimPlot(seu_obj_data,
              reduction = "umap",
              group.by = "stage"
) +
  theme_bw()
p3 <- DimPlot(seu_obj_data,
              reduction = "umap",
              group.by = "platform"
) +
  theme_bw()
p1 + p2 + p3

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
  "MSA41", # B(MS4A1+)
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
library(ggplot2)
FeaturePlot(seu_obj_data,
            features = "CD19",
            reduction = "umap",
            coord.fixed = T,
            order = T,
            cols = viridis(10)
) +
  scale_color_viridis(discrete = F, option = "inferno")

DotPlot(seu_obj_data, features = unique(mainmarkers), group.by = "seurat_clusters") +theme_bw()+
  RotatedAxis() +
  scale_x_discrete("") +
  scale_y_discrete("") +
  # coord_flip() +
  scale_color_viridis(discrete = F, option = "C")

ggsave2(paste0("Fig5.FeaturePlot_mainmarkers.png"),
        path = paste0("Results/"),
        width = 20, height = 12, units = "cm"
)

for (i in seq_along(mainmarkers)) {
  if (mainmarkers[i] %in% rownames(scRNA_harmony)) {
    FeaturePlot(scRNA_harmony,
                features = mainmarkers[i],
                coord.fixed = T,
                order = T,
                cols = viridis(10)
    ) +
      scale_color_viridis(discrete = F, option = "inferno")
    ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".png"),
            path = paste0("Results/marker/"),
            width = 12, height = 12, units = "cm"
    )
  }
}


celltypes <- c(
  "Follicular B cells", 
  "Follicular B cells",
  "Follicular B cells",
  "Follicular B cells",  
  "Plasma", 
  "UND", # 05
  "Plasma",
  "Plasma",
  "Follicular B cells", 
  "Plasma", 
  "Plasma", # 10
  "UND", 
  "Follicular B cells",
  "Follicular B cells", 
  "Follicular B cells", 
  "Plasma", # 15
  "UND", 
  "Plasma", 
  "UND",
  "Plasma", 
  "UND",# 20
  "Plasma", 
  "UND",
  "UND",
  "Memory B"
)


scRNA_harmony@meta.data$celltypes <- NA
for (i in 1:length(scRNA_harmony@meta.data$celltypes)) {
  scRNA_harmony@meta.data$celltypes[i] <- celltypes[scRNA_harmony@meta.data$seurat_clusters[i]]
}
table(scRNA_harmony$celltypes)

p1 <- DimPlot(scRNA_harmony,
        reduction = "umap",
        group.by = "celltypes"
) + theme_bw()
p2 <- DimPlot(scRNA_harmony,
        reduction = "umap",
        group.by = "orig.ident"
) + theme_bw()
p3 <- DimPlot(scRNA_harmony,
        reduction = "umap",
        group.by = "stage"
) + theme_bw()
p1+p2+p3

if (F) {
  
  Cellratio <- prop.table(table(scRNA_harmony$celltypes, scRNA_harmony$orig.ident), margin = 2)#计算各组样本不同细胞群比例
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
  Cellratio <- prop.table(table(Idents(seu_obj_data), seu_obj_data$orig.ident), margin = 2)#计算各组样本不同细胞群比例
  Cellratio <- data.frame(Cellratio)
  library(reshape2)
  cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
  rownames(cellper) <- cellper[,1]
  cellper <- cellper[,-1]
  
  ###添加分组信息
  sample <- seu_obj_data$orig.ident
  group <-seu_obj_data$stage
  samples <- data.frame(sample, group)#创建数据框
  
  rownames(samples)=samples$sample
  
  # cellper$sample <- samples[rownames(cellper),'sample']#R添加列
  cellper$sample <- rownames(cellper)#R添加列
  # cellper$group <- samples[rownames(cellper),'group']#R添加列
  cellper$group <- samples[rownames(cellper),'group']#R添加列
  cellper$group <- ""
  for (i in 1:nrow(cellper)) {
    stage <- samples[which(samples$sample== cellper$sample[i])[1] , "group"]
    cellper$group[i] <- stage
  }
  
  ###作图展示
  pplist = list()
  sce_groups = c("B","Follicular B cells","Plasma")
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
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
      labs(title = group_,y='Percentage') +
      geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
    
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
  plot_grid(pplist[['B']],
            pplist[['Follicular B cells']],
            pplist[['Plasma']])
  
}

# 获取当前用的Idents
Idents(object = scRNA_harmony)
levels(scRNA_harmony)
Idents(scRNA_harmony) <- "seurat_clusters"
new.cluster.ids <- as.character(new.cluster.ids)
names(new.cluster.ids) <- levels(scRNA_harmony)
scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)

scRNA_harmony@meta.data$new.cluster.ids

scRNA_harmony <- RenameIdents(scRNA_harmony,
                              `0` = "B", `1` = "B", `2` = "Follicular B cells",
                              `3` = "Follicular B cells", `4` = "B", `5` = "B", `6` = "Plasma",
                              `7` = "Plasma", `8` = "Plasma", `9` = "Plasma",`10` = "B", `11` = "B",`12` = "B",
                              `13` = "Cancer", `14` = "Plasma", `15` = "Myeloid",
                              `16` = "Plasma", `17` = "Myeloid", `18` = "B",
                              `19` = "B", `20` = "Plasma", `21` = "Plasma",`22` = "B",`23` = "Follicular B cells"
)

DimPlot(scRNA_harmony, label = T) + NoLegend()+theme_bw()
DimPlot(scRNA_harmony, group.by = "orig.ident") +theme_bw()+ NoLegend()
DotPlot(scRNA_harmony, features = unique(mainmarkers), group.by = "new.cluster.ids") + RotatedAxis() +
  scale_x_discrete("") + scale_y_discrete("")

print(DimPlot(scRNA_harmony, reduction = "umap", group.by = c("new.cluster.ids")))

p <- FeaturePlot(scRNA_harmony, features = markers, ncol = 8)
p

obj_cells <- c("B", "Plasma", "Follicular B cells")
seu_obj_data_obj_cells <- list()
for (i in 1:length(obj_cells)) {
  obj_cell <- obj_cells[i]
  seu_obj_data_obj_cell <- scRNA_harmony[, (scRNA_harmony$new.cluster.ids == obj_cell)]
  seu_obj_data_obj_cells[[i]] <- seu_obj_data_obj_cell
}
seu_obj_data <- merge(seu_obj_data_obj_cells[[1]], seu_obj_data_obj_cells[2:length(seu_obj_data_obj_cells)])
dim(seu_obj_data)
table(seu_obj_data$new.cluster.ids)
save(seu_obj_data, file="../scGRN-L0_data/seu_obj_data_B_samples.Rdata")



