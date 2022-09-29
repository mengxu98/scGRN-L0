

# Load all ibraries and script files
library("tidyr")
source("framework_main.R")
source("ground-truth.R")

pathway <- paste0("../scGRN-L0_data/BEELINE-data/inputs/scRNA-Seq/", "mDC", "/") # i=1
dataway <- "ExpressionData.csv" #
dataway_PseudoTime <- "PseudoTime.csv"
# goldway <- "/gold/gold.tsv"

# expression_dataset_test <- constructExpressionMatrixFromFile(paste0(pathway, dataway))
# expression_dataset_test <- read.table(paste0(pathway, dataway), header = T) %>% as.matrix()
expression_dataset_test <- read.csv(paste0(pathway, dataway),
  header = T,
  row.names = 1
) %>%
  t() %>%
  as.matrix()
PseudoTime <- read.csv(paste0(pathway, dataway_PseudoTime),
  header = T
)
PseudoTime <- PseudoTime[order(PseudoTime$PseudoTime), ] %>% as.data.frame()
rownames(PseudoTime) <- PseudoTime$X

expression_dataset_test <- data_filetr(expression_dataset_test,
  dataset_dir = "../scGRN-L0_data/BEELINE-Networks/Networks/human/",
  database = "All"
)

TFs <- read.csv("../scGRN-L0_data/BEELINE-Networks/human-tfs.csv")

expression_dataset_test <- expression_dataset_test[rownames(PseudoTime)[1:50], intersect(TFs$TF, colnames(expression_dataset_test))]

expression_dataset_test <- expression_dataset_test[, 1:100]
dim(expression_dataset_test)
min(expression_dataset_test)

# Run algorithm using default settings, writes result to output.txt
NIMEFI(expression_dataset_test,
  GENIE = F, SVM = F, EL = T, penalty = "L0L2",
  outputFileName = "output/output_net_hESC_L0",
  outputFileFormat = "txt",
  SVMRankThreshold = 5, SVMEnsembleSize = 100,
  ELPredSampleMin = 20, ELPredSampleMax = 80,
  ELExpSampleMin = 20, ELExpSampleMax = 80,
  ELRankThreshold = 5, ELEnsembleSize = 1
)

ground_truth("output/output_net_hESC_L0.txt",
  "output/",
  dataset_dir = "../scGRN-L0_data/BEELINE-Networks/Networks/human/",
  database = "All"
)

evaluationObject <- prepareEval("output/output_net_hESC_L0.txt",
  paste0("output/ground_truth.tsv"),
  totalPredictionsAccepted = 100000
)

# test
ground_truth_simulation(intput = "output/test1.txt",
  output = "output/",
  dataset_dir = "../scGRN-L0_data/BEELINE-data/inputs/Synthetic/dyn-BF/dyn-BF-100-1/",
  file = "refNetwork.csv"
)
evaluationObject <- prepareEval("output/test1.txt",
  paste0("output/ground_truth.tsv"),
  totalPredictionsAccepted = 100000
)

L0_AUROC <- calcAUROC(evaluationObject)
L0_AUPR <- calcAUPR(evaluationObject)
L0_AUROC
L0_AUPR

