

output <- "../scGRN-L0_output/"
sc_T_cell <- sc[which(sc$CellType == "T_cell"), ]
table(sc$CellType)
sc_T_cell= sc[, Idents(sc) %in% c( "T_cell")]
table(MC.seurat$CellType)
mc_T_cell= MC.seurat[, Idents(MC.seurat) %in% c( "T_cell")]

matrix <- as.data.frame(mc_T_cell@assays$RNA@scale.data)
matrix <- as.data.frame(sc_T_cell@assays$RNA@scale.data)
head(matrix[1:5,1:5])
max(matrix)
L0Dynamic <- L0REG(
  matrix = matrix,
  regulators = hvg,
  targets = hvg,
  penalty = "L0"
)
max(L0Dynamic$weight)
write.table(L0Dynamic,
            paste0(output, "GRN_L0.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F
)

ground_truth_h(
  intput = paste0(output, "GRN_L0.txt"),
  output = output,
  dataset_dir = hnetwork_data_dir,
  database = "STRING"
)

evaluationObject <- prepareEval(paste0(output, "GRN_L0.txt"),
                                paste0(paste0(output, "ground_truth.tsv")),
                                totalPredictionsAccepted = 100000)
AUROC_L0 <- calcAUROC(evaluationObject)
AUPRC_L0 <- calcAUPR(evaluationObject)
AUROC_L0
AUPRC_L0
