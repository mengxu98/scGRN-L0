####################
# PACKAGES required:
# kSamples
# glmnet
# ppcor
# pracma
# R.matlab
####################
library(kSamples)
library(glmnet)
library(ppcor)
library(R.matlab)

## *** Data loading ***

mat <- readMat('In silico single cell data/20_nets_10genes_8UNEVENtime_sigma01B_no_initial_points2.mat')
time <- as.vector(mat$time.points)
numGENES <- as.vector(mat$n)
AUROC <- vector()
AUPR <- vector()
for (numEXAMPLES in 1:dim(mat$networks)[3]) {
  data_time_series <- mat$data.tot.array[[numEXAMPLES]][[1]]
  singleCELLdata <- list()
  for (i in 1:mat$num.time.points) {
    singleCELLdata[[i]] <- data_time_series[,i,]
  }
  genes <- vector(length=numGENES)
  for (i in 1:numGENES) {
    genes[i] <- sprintf('Gene %d',i)
  }
  totDATA <- matrix(nrow = 0, ncol = dim(data_time_series)[3])
  for (i in 1:mat$num.time.points) {
    totDATA <- rbind(totDATA,data_time_series[,i,])
  }
  
  DATA <- list(time=time, numGENES=numGENES, singleCELLdata=singleCELLdata, genes=genes, totDATA=totDATA)
  
  ## *** SINCERITIES ***
  
  library(kSamples)
  library(glmnet)
  library(ppcor)
  SINCERITITES <- dget("SINCERITIES functions/SINCERITIES_L0.R")
  result <- SINCERITITES(DATA,distance=1, method = 1,noDIAG = 1) #,SIGN = 0
  # result <- SINCERITITES(DATA,distance=1,method = "L0",noDIAG = 1,SIGN = 1)
  adj_matrix <- result$adj_matrix
  
  ## *** Performance Evaluation ***
  
  #Gold standard GRN
  a <- mat$networks[,,numEXAMPLES]
  a[row(a)==col(a)] <- 0
  if(SIGN==0){
    a[which(a!=0)] <- 1
  }
  
  #Final ranked list, AUROC and AUPR
  adj_matrix <- adj_matrix/max(adj_matrix)
  library(pracma)
  auc_from_ranks_TC_sign <- dget("SINCERITIES functions/auc_from_ranks_TC_sign.R")
  AUCresult <- auc_from_ranks_TC_sign(adj_matrix,a,1000)
  AUROC[numEXAMPLES] <- AUCresult$AUROC
  AUPR[numEXAMPLES] <- AUCresult$AUPR
  final_ranked_predictions <- dget("SINCERITIES functions/final_ranked_predictions.R")
  table <- final_ranked_predictions(adj_matrix,DATA$genes,SIGN=SIGN,fileNAME=sprintf('results4insilicoNETWORK%d',numEXAMPLES),saveFile = TRUE)
}
AUC <- cbind(AUROC,AUPR)
m <- apply(AUC,2,mean)
s <- apply(AUC,2,std)
AUC <- rbind(AUC,rbind(m,s))
AUC

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

# AUROC      AUPR
# 0.6567235 0.1653079
# 0.6898657 0.1450217
# 0.7343750 0.1665757
# 0.8417774 0.3530668
# 0.7944444 0.1909405
# 0.5633299 0.1189544
# 0.5486597 0.2268660
# 0.8233333 0.2301572
# 0.7957099 0.4422236
# 0.8733401 0.3804112
# 0.6593463 0.2608745
# 0.6019022 0.1178857
# 0.8153409 0.6646166
# 0.8137255 0.5594679
# 0.8484848 0.4288358
# 0.8687436 0.4329460
# 0.7322222 0.2658719
# 0.8960784 0.6878457
# 0.7078598 0.2329139
# 0.9266667 0.6946476
# m 0.7595965 0.3382715
# s 0.1114742 0.1914145

# !Spearman SIGN==1 !abs
# m 0.6884634 0.25856028

# !Spearman SIGN==0 !abs
# m 0.6884634 0.25856028

# !Spearman SIGN==0 abs
# m 0.6884634 0.25856028

# !Spearman SIGN==0 abs
# m 0.6884634 0.25856028

