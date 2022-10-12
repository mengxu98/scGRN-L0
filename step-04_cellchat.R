

save(scRNA_harmony_data, file = paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_SCT_harmony_PCA.Rdata"))
s=1
load(paste0("/data/mengxu/data/all/lung_stage-", s, "_seu_SCT_harmony_PCA.Rdata"))

scRNA_harmony <- FindNeighbors(scRNA_harmony, dims = pc.num) %>% FindClusters(resolution = 1)


Idents(scRNA_harmony) <- scRNA_harmony$SCT_snn_res.1
annotation_curated_main <- read_excel(paste0("/data/mengxu/results/figure/stage-", s, "/curated_annotation_main.xlsx"))
new_ids_main <- annotation_curated_main$main_cell_type
names(new_ids_main) <- levels(scRNA_harmony)
scRNA_harmony <- RenameIdents(scRNA_harmony, new_ids_main)
scRNA_harmony@meta.data$main_cell_type <- Idents(scRNA_harmony)

DotPlot(scRNA_harmony, features = unique(mainmarkers), group.by = "main_cell_type") +
  RotatedAxis() +
  scale_x_discrete("") +
  scale_y_discrete("") +
  # coord_flip() +
  scale_color_viridis(discrete = F, option = "C")

### plot_umap
###
DimPlot(scRNA_harmony,
  reduction = "umap",
  group.by = "main_cell_type",
  label = T,
  label.box = T,
  cols = colP,
  repel = T
  # pt.size = 0.2
) + theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + NoLegend() +
  scale_color_viridis(discrete = T, option = "C")
ggsave2("SuppFig6.clusters_harmony_umap.png",
  path = paste0("/data/mengxu/results/figure/stage-", s),
  width = 15, height = 15, units = "cm"
)
###
seu_harmony <- subset(scRNA_harmony, idents = c("B", "T", "Myeloid", "Cancer", "Fibroblast", "Endothelial", "Alveolar", "Hepatocytes"))

seu_harmony@meta.data$main_cell_type <- Idents(seu_harmony)

cellchat <- createCellChat(object = seu_harmony, meta = seu_harmony@meta.data, group.by = "main_cell_type")
levels(cellchat@idents)


###
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # Cell-Cell Contact, ECM-Receptor,  Secreted Signaling
cellchat@DB <- CellChatDB.use
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 8)

# future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)


cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2, 3), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("CXCL")
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells
vertex.receiver <- seq(1, 4) # a numeric vector.
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Heatmap
par(mfrow = c(1, 1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("CMP", 6), rep("DC", 6), rep("GMP", 6)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1, ] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver <- seq(1, 4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, vertex.receiver = vertex.receiver)


# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver

levels(cellchat@idents)
# vertex.receiver = seq(1,4)
# for (i in 1:length(pathways.show.all)) {
#   # Visualize communication network associated with both signaling pathway and individual L-R pairs
#   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
#   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#   ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
# }

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL", "CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL", "CXCL", "FGF"))
netVisual_bubble(cellchat, sources.use = c(3, 4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5, legend.pos.y = 30)

# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1, 2, 3, 4), targets.use = 8, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1, 2, 3, 4), targets.use = c(5:11), signaling = c("CCL", "CXCL"), legend.pos.x = 8)
#> Note: The second link end is drawn out of sector 'CXCR4 '.
#> Note: The first link end is drawn out of sector 'CXCL12 '.

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1, 2, 3, 4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)
#> Note: The first link end is drawn out of sector 'CXCL '.

plotGeneExpression(cellchat, signaling = "CXCL")
#> Registered S3 method overwritten by 'spatstat':
#>   method     from
#>   print.boxx cli
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.

plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)



# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))




library(NMF)
selectK(cellchat, pattern = "outgoing")


nPatterns <- 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)


library(ggalluvial)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)


cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

# saveRDS(cellchat, file = "cellchat_humanSkin_LS.rds")


remotes::install_github("hafen/rminiconda")
py <- rminiconda::find_miniconda_python("my_python")
reticulate::use_python(py, required = TRUE)
reticulate::py_install(packages = "umap-learn")
