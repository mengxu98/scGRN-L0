

# Load packages--------------------------------------------------
library("kSamples")
library("glmnet")
library("ppcor")
library("R.matlab")
library("L0Learn")

# *** Data loading ***
mat <- readMat("In silico single cell data/20_nets_10genes_6UNEVENtime_sigma01B_no_initial_points2.mat")
time <- as.vector(mat$time.points)
numGENES <- as.vector(mat$n)
AUROC <- vector()
AUPR <- vector()
SIGN <- 1
for (numEXAMPLES in 1:dim(mat$networks)[3]) {
  data_time_series <- mat$data.tot.array[[numEXAMPLES]][[1]]
  singleCELLdata <- list()
  for (i in 1:mat$num.time.points) {
    singleCELLdata[[i]] <- data_time_series[, i, ]
  }
  genes <- vector(length = numGENES)
  for (i in 1:numGENES) {
    genes[i] <- sprintf("Gene %d", i)
  }
  totDATA <- matrix(nrow = 0, ncol = dim(data_time_series)[3])
  for (i in 1:mat$num.time.points) {
    totDATA <- rbind(totDATA, data_time_series[, i, ])
  }

  DATA <- list(time = time, numGENES = numGENES, singleCELLdata = singleCELLdata, genes = genes, totDATA = totDATA)

  # *** SINCERITIES ***
  SINCERITITES <- dget("SINCERITIES functions/SINCERITIES_L0.R")
  result <- SINCERITITES(DATA, distance = 3, method = 1, noDIAG = 1, SIGN = 1)
  # result <- SINCERITITES(DATA,distance=1,method = "L0",noDIAG = 1,SIGN = 1)
  adj_matrix <- result$adj_matrix

  # *** Performance Evaluation ***

  # Gold standard GRN
  a <- mat$networks[, , numEXAMPLES]
  a[row(a) == col(a)] <- 0
  if (SIGN == 0) {
    a[which(a != 0)] <- 1
  }

  # Final ranked list, AUROC and AUPR
  adj_matrix <- adj_matrix / max(adj_matrix)
  library(pracma)
  auc_from_ranks_TC_sign <- dget("SINCERITIES functions/auc_from_ranks_TC_sign.R")
  AUCresult <- auc_from_ranks_TC_sign(adj_matrix, a, 1000)
  AUROC[numEXAMPLES] <- AUCresult$AUROC
  AUPR[numEXAMPLES] <- AUCresult$AUPR
  final_ranked_predictions <- dget("SINCERITIES functions/final_ranked_predictions.R")
  table <- final_ranked_predictions(adj_matrix, DATA$genes, SIGN = SIGN, fileNAME = sprintf("results4insilicoNETWORK%d", numEXAMPLES), saveFile = TRUE)
}
AUC <- cbind(AUROC, AUPR)
m <- apply(AUC, 2, mean)
s <- apply(AUC, 2, std)
AUC <- rbind(AUC, rbind(m, s))
print(AUC)

<<<<<<< HEAD
#       AUROC      AUPR
#   0.7394444 0.1795571
#   0.7378627 0.3232821
#   0.7744659 0.4422365
#   0.7610599 0.5214426
#   0.7837154 0.2523979
#   0.6440977 0.1323813
#   0.5536117 0.2325753
#   0.8583028 0.2509836
#   0.7922222 0.2964554
#   0.9000000 0.3686655
#   0.7860061 0.2319159
#   0.9870924 0.6996460
#   0.8475379 0.6187193
#   0.8811888 0.6843920
#   0.8532197 0.5825992
#   0.9264556 0.5444362
#   0.9372222 0.5181557
#   0.8941176 0.5647393
#   0.7973485 0.5517647
#   0.9474510 0.8410191
# m 0.8201211 0.4418682
# s 0.1057922 0.1987979

# distance = 3, method = 1, noDIAG = 1,SIGN = 1
#       AUROC      AUPR
=======
# AUROC      AUPR
# 0.6490390 0.1356164
# 0.7136752 0.2588762
# 0.7732199 0.1538346
# 0.8529900 0.6390968
# 0.8019780 0.2354546
# 0.5876404 0.1221352
# 0.5375999 0.1740624
# 0.8252015 0.2707332
# 0.7917393 0.2904240
# 0.8756883 0.3542469
# 0.6879469 0.2623950
# 0.6358696 0.1267483
# 0.8371212 0.6787861
# 0.8360784 0.5724397
# 0.8655303 0.4411362
# 0.9039837 0.4834967
# 0.7833333 0.3859044
# 0.9192157 0.7238272
# 0.7267992 0.2296734
# 0.9360784 0.7289519
# m 0.7770364 0.3633920
# s 0.1122153 0.2076085

# distance = 3, method = 1, noDIAG = 1,SIGN = 1
#       AUROC      AUPR
#   0.6426231 0.1215555
#   0.7545788 0.2778658 0
#   0.7487050 0.2307468
#   0.8343023 0.6116839
#   0.8418152 0.2515986 0
#   0.5629934 0.1189508
#   0.3924566 0.1378816
#   0.7902869 0.1969678
#   0.8336806 0.4842825
#   0.8977778 0.3299649
#   0.6547497 0.3641774
#   0.7038043 0.2434772
#   0.8096591 0.6492946
#   0.8321569 0.5459008
#   0.8797348 0.5228053
#   0.8917263 0.4945470
#   0.7311111 0.3587481
#   0.9333333 0.7788478
#   0.7727273 0.4789326
#   0.9317647 0.7158874
# m 0.7719994 0.3957058
# s 0.1331083 0.2016986

#       AUROC      AUPR
>>>>>>> origin/master
#   0.7105556 0.1712579
#   0.6962298 0.3079329
#   0.7244855 0.4258128
#   0.7481670 0.5147832
#   0.7672366 0.2490852
#   0.6036659 0.1236446
#   0.3804910 0.1416284
#   0.8565263 0.2500017
#   0.7783333 0.2428947
#   0.9009471 0.3751932
#   0.7609806 0.2207239
#   0.9891304 0.7158234
#   0.8271780 0.6112080
#   0.8622621 0.6651726
#   0.8257576 0.5607282
#   0.9203269 0.5308015
#   0.9383333 0.5231933
#   0.8803922 0.5343286
#   0.7698864 0.5437078
#   0.9419608 0.8395081
# m 0.7941423 0.4273715
# s 0.1375813 0.2064932
