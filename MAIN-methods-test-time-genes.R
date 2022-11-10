

# Load all the script files, libraries and functions
library("tidyr")
library("kSamples")
library("glmnet")
library("ppcor")
library("L0Learn")
library("pracma")
library("ggplot2")
library("reshape2")
library("RColorBrewer")
library("patchwork")
source("framework_main.R")
source("ground-truth.R")
source("Function-L0REG.R")
# source("Function-L0REG_delay.R")

uploading <- dget("SINCERITIES functions/uploading.R")
SINCERITITES <- dget("SINCERITIES functions/SINCERITIES.R")
SINCERITITES_L0 <- dget("SINCERITIES functions/SINCERITIES_L0.R")
auc_from_ranks_TC_sign <- dget("SINCERITIES functions/auc_from_ranks_TC_sign.R")
final_ranked_predictions <- dget("SINCERITIES functions/final_ranked_predictions.R")

gene_num <- c(
  "20",
  "40",
  "60",
  "80",
  "100"
)
runtime <- matrix(0, nrow = 5, ncol = 5)
colnames(runtime) <- gene_num
rownames(runtime) <- c("L0DWGRN", "GENIE3", "SINCERITITES", "PPCOR", "LEAP")
output <- "../scGRN-L0_output/output_Synthetic/"
for (i in 1:length(gene_num)) {
  simulation_data_dir <- paste0("../scGRN-L0_data/BEELINE-data/inputs/scRNA-Seq/hESC/")
  simulation_data_file <- "ExpressionData.csv"
  simulation_PseudoTime_file <- "PseudoTime.csv"
  simulation_data <- read.csv(paste0(simulation_data_dir, simulation_data_file), row.names = 1)
  head(simulation_data[1:3, 1:3])

  simulation_PseudoTime <- read.csv(paste0(simulation_data_dir, simulation_PseudoTime_file), row.names = 1)
  simulation_data_a <- cbind.data.frame(t(simulation_data), h = simulation_PseudoTime$PseudoTime) %>% na.omit()
  simulation_data_new <- simulation_data_a[1:500, (ncol(simulation_data_a) - as.numeric(gene_num[i]) - 1):ncol(simulation_data_a)]
  write.csv(simulation_data_new, paste0(simulation_data_dir, "ExpressionData_all_genes", ".csv"), row.names = FALSE)
  data_grn <- read.csv(paste0(simulation_data_dir, "ExpressionData_all_genes.csv"),
    header = T
  ) %>% as.matrix()
  PseudoTime <- data_grn[, ncol(data_grn)]
  data_grn <- data_grn[, -ncol(data_grn)]
  # --------------------------------------------------
  ptm <- proc.time()
  DATA <- uploading(paste0(simulation_data_dir, "ExpressionData_all_genes", ".csv"))
  # result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 1)
  result <- SINCERITITES(DATA, distance = 1, method = 1, noDIAG = 1, SIGN = 0)
  time <- proc.time() - ptm
  runningtime_SINCERITITES <- time[3]
  # --------------------------------------------------
  ptm <- proc.time()
  library(GENIE3)
  weightMat <- GENIE3(
    exprMatrix = t(data_grn),
    nCores = 8,
    nTrees = 1000,
    verbose = TRUE
  )
  time <- proc.time() - ptm
  runningtime_GENIE3 <- time[3]
  # --------------------------------------------------
  ptm <- proc.time()
  L0Dynamic <- L0REG(
    matrix = t(data_grn),
    regulators = colnames(data_grn),
    targets = colnames(data_grn),
    # maxSuppSize = (nrow(data_grn) * 0.5),
    penalty = "L0"
  )
  time <- proc.time() - ptm
  runningtime_L0 <- time[3]
  # --------------------------------------------------
  library(ppcor)
  ptm <- proc.time()
  inputExpr <- t(data_grn)
  geneNames <- rownames(inputExpr)
  rownames(inputExpr) <- c(geneNames)
  pcorResults <- pcor(x = t(as.matrix(inputExpr)), method = "spearman")
  DF <- data.frame(
    Gene1 = geneNames[c(row(pcorResults$estimate))], Gene2 = geneNames[c(col(pcorResults$estimate))],
    corVal = c(pcorResults$estimate), pValue = c(pcorResults$p.value)
  )
  outDF <- DF[order(DF$corVal, decreasing = TRUE), ]
  time <- proc.time() - ptm
  runningtime_ppcor <- time[3]
  # --------------------------------------------------
  ptm <- proc.time()
  library(LEAP)
  geneNames <- rownames(inputExpr)
  rownames(inputExpr) <- c()
  MAC_results <- MAC_counter(
    data = inputExpr, # max_lag_prop=maxLag,
    MAC_cutoff = 0,
    file_name = "temp", lag_matrix = FALSE, symmetric = FALSE
  )
  Gene1 <- geneNames[MAC_results[, "Row gene index"]]
  Gene2 <- geneNames[MAC_results[, "Column gene index"]]
  Score <- MAC_results[, "Correlation"]
  outDF <- data.frame(Gene1, Gene2, Score)
  time <- proc.time() - ptm
  runningtime_LEAP <- time[3]

  runtime[1, i] <- runningtime_L0
  runtime[2, i] <- runningtime_GENIE3
  runtime[3, i] <- runningtime_SINCERITITES
  runtime[4, i] <- runningtime_ppcor
  runtime[5, i] <- runningtime_LEAP
}

runtime
mydata <- melt(runtime, id = rownames(runtime))
colnames(mydata) <- c("Method", "Cellnum", "Runtime")
runtime_gene <- ggplot(data = mydata, aes(x = factor(Cellnum), y = Runtime, group = Method, color = Method, shape = Method)) +
  geom_point(
    size = 2.5
  ) +
  geom_line() +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  xlab("Gene Number") +
  ylab("Run Time (Second)") +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman"),
    # legend.position = c(.75, .25),
    # legend.box.background = element_rect(color = "black"),
    legend.position = "right"
  )
runtime_gene
ggsave(paste0("../scGRN-L0_output/output_Synthetic/Methods-contrast-Runtime-gene.png"), width = 5, height = 2, dpi = 600)

runtime_cell + runtime_gene +
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") +
    theme(text = element_text(size = 10)) +
    theme(text = element_text(family = "Times New Roman"))
ggsave(paste0("../scGRN-L0_output/output_Synthetic/Methods-contrast-runningtime.png"), width = 8, height = 3, dpi = 600)
