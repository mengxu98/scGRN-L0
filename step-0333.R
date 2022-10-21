

cellInfor1 <- cellInfor[order(cellInfor$PseTime.Lineage1),]
samples_labels <- c()
for (i in 1:10) {
  s <- floor(length(cellInfor1$PseTime.Lineage1)*(i-1)/10)+1
  e <- floor(length(cellInfor1$PseTime.Lineage1)*i/10)
  sample <- cellInfor1$UniqueCell_ID[s:e]
  samples_label <- data.frame(cell_id=sample,time_label=paste0("t",i))
  samples_labels <- rbind.data.frame(samples_labels,samples_label)
}

metadata <- seu_obj_data@meta.data
metadata$cell_id <- rownames(metadata)
metadata$sample_id <- metadata$orig.ident
metadata <- left_join(x = metadata, y = samples_labels, by = "cell_id")
rownames(metadata) <- metadata$cell_id
seu_obj_data <- AddMetaData(seu_obj_data, metadata = metadata)


Cellratio1 <- prop.table(table(seu_obj_data$main_cell_type, seu_obj_data$time_label), margin = 2)#计算各组样本不同细胞群比例
Cellratio1
Cellratio1 <- as.data.frame(table(seu_obj_data$main_cell_type, seu_obj_data$time_label))
Cellratio1 <- as.data.frame(Cellratio1)
colourCount1 = length(unique(Cellratio1$Var1))
library(ggplot2)
ggplot(Cellratio1) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),position="dodge",stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_x_discrete(labels = c("t1","t2","t3","t4","t5","t6","t7","t8","t9","t10"))+
  # coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

Cellratio2 <- prop.table(table(seu_obj_data$tissue_type, seu_obj_data$time_label), margin = 2)#计算各组样本不同细胞群比例
Cellratio2

Cellratio2 <- as.data.frame(table(seu_obj_data$tissue_type, seu_obj_data$time_label))
Cellratio2 <- as.data.frame(Cellratio2)
colourCount2 = length(unique(Cellratio2$Var1))
library(ggplot2)
ggplot(Cellratio2) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),position="dodge",stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_x_discrete(labels = c("t1","t2","t3","t4","t5","t6","t7","t8","t9","t10"))+
  # coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))







cell_types <- FetchData(seu_obj_data, vars = c("sample_id", "main_cell_type", "time_label")) %>%
  mutate(main_cell_type = factor(main_cell_type, levels = c("Follicular B cells", "Germinal center B cells","Plasma"))) %>%
  mutate(sample_id = factor(sample_id, levels = rev(c("p018t", "p019t", "p023t", "p024t", "p027t", "p028t",
                                                      "p030t", "p031t", "p032t", "p033t", "p034t", "p018n",
                                                      "p019n", "p027n", "p028n", "p029n", "p030n", "p031n",
                                                      "p032n", "p033n", "p034n"))))

ggplot(data = cell_types) + 
  geom_bar(mapping = aes(x = time_label, fill = main_cell_type, ), position = "fill", width = 0.75) +
  # scale_fill_manual(values = use_colors) +
  scale_x_discrete(labels = c("t1","t2","t3","t4","t5","t6","t7","t8","t9","t10"))+
  coord_flip()

seu_obj_data$tissue_type
cell_types <- FetchData(seu_obj_data, vars = c("tissue_type", "main_cell_type", "time_label")) %>%
  mutate(main_cell_type = factor(main_cell_type, levels = c("Follicular B cells", "Germinal center B cells","Plasma"))) %>%
  mutate(tissue_type = factor(tissue_type, levels = rev(c("Normal", "Tumor"))))

ggplot(data = cell_types) + 
  geom_bar(mapping = aes(x = time_label, fill = tissue_type, ), position = "fill", width = 0.75) +
  # scale_fill_manual(values = use_colors) +
  scale_x_discrete(labels = c("t1","t2","t3","t4","t5","t6","t7","t8","t9","t10"))+
  coord_flip()

p <- ggplot(cell_types,aes(x=time_label,y=tissue_type,fill=tissue_type))+geom_bar(position="dodge",stat="identity")
p+xlab("time_label") + ylab("tissue_type") + labs(fill="count")




