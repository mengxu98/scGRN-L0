

rm(list=ls())
library(igraph)
library(RColorBrewer)
library(broom)
library(olsrr)
library(car)
library(HH)
library(Metrics)
library(plotmo)
library(L0Learn)
library(tidyverse)
library(caret)
library(glmnet)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(ggthemes)
library(psych)
library("ggpubr")
library("GENIE3")

TRN_L0_all_list <- list()
links_list <- list()
for (d in c(1,3,4)) {
  
  dataset <- read.table(paste0("DREAM5_network_inference_challenge/Network",
                               d,
                               "/input data/net",
                               d,
                               "_expression_data.tsv"), 
                        header = T)
  # dataset <- dataset[,1:500]
  #GENIE3
  if (T) {
    genie3_data <- as.matrix(dataset) %>% t()
    row.names(genie3_data) <-  row.names(genie3_data)
    pre_time <- proc.time()
    weightMat <- GENIE3(genie3_data, 
                        nCores = 4, 
                        verbose=TRUE)
    time_GRN <- proc.time() - pre_time
    
    links <- getLinkList(weightMat)
      
      write.table(links, 
                  paste0("DREAM5_network_inference_challenge/Evaluation scripts/INPUT/predictions/myteam/GENIE3/DREAM5_NetworkInference_myteam_Network",
                         d,
                         ".txt"),
                  col.names = F,
                  row.names = F,
                  sep = "\t",
                  quote = F)
      
      links_list[[d]] <- links
      
  }
  
  if (T) {
    maxSNVSize <- ncol(dataset)
    TRN_L0_all <- c()
    evaluate <- c()
    for(i in 1:ncol(dataset)){
      
      target_gene <- colnames(dataset)[i]
      
      # cat(paste0("Computing gene ",i ,"/", maxSNVSize, ": ",target_gene))
      message(paste0("Computing gene ",i ,"/", maxSNVSize, ": ",target_gene))
      
      Y <- dataset[, i]
      X <- dataset[, -i] %>% as.matrix()
      row.names(X) = row.names(X)
      
      # set.seed(2022)
      # train_idx <- sample(nrow(X), 0.7 * nrow(X))
      # train_data_X <- X[train_idx, ]
      # train_data_Y <- Y[train_idx]
      # test_data_X <- X[-train_idx, ]
      # test_data_Y <- Y[-train_idx]
      # 
      # cvfit_L0_validation <- L0Learn.cvfit(train_data_X, train_data_Y, nFolds = 10)
      # # plot(cvfit_L0_validation)
      # fit_L0_validation <- L0Learn.fit(train_data_X, train_data_Y,
      #                                  penalty="L0",
      #                                  maxSuppSize = maxSNVSize)
      # # plot(fit_L0_validation)
      # 
      # # print(fit_L0_validation)
      # 
      # # Extract coefficient at middle lambda
      # fit_L0_information <- as.data.frame(print(fit_L0_validation))
      # fit_L0_information <- fit_L0_information[order(fit_L0_information$suppSize, decreasing = TRUE), ]
      # lambda_L0 <- fit_L0_information$lambda[1]
      # gamma_L0 <- fit_L0_information$gamma[1]
      # 
      # test_data_y_cat = predict(fit_L0_validation,
      #                           newx=test_data_X,
      #                           lambda=lambda_L0,
      #                           gamma=gamma_L0)
      # test_data_y_hat=as.vector(test_data_y_cat)
      # ###
      # 
      # re_status <- cbind(test_data_Y, test_data_y_hat)
      # # head(re_status)
      # re_status <- as.data.frame(re_status)
      # colnames(re_status) <- c('raw', 'pre')
      # 
      # peak_corr <- corr.test(re_status$raw,
      #                        re_status$pre,
      #                        # method = "spearman",
      #                        adjust = "none")
      # 
      # if (peak_corr$p < 0.05 && peak_corr$r > 0.8) {
        
      nfolds <- round(sqrt(dim(X)[1]))
      
      # Cross validate
      cross <-glmnet::cv.glmnet(X,Y,nfolds=nfolds,family="gaussian",standardize = FALSE)
      # Retrain optimal 
      bestModel<-glmnet::glmnet(X,Y,family="gaussian",lambda=cross$lambda.min,standardize =FALSE)	
      wghts<-abs(as.vector(bestModel$beta))
      
        ####L0####
        cvfit_L0 <- L0Learn.cvfit(X, Y, nFolds = 10)
        fit_L0 <- L0Learn.fit(X, Y,
                              penalty="L0",
                              maxSuppSize = maxSNVSize)
        # print(fit_L0)
        
        # Extract coefficient at middle lambda
        fit_L0_information <- as.data.frame(print(fit_L0))
        fit_L0_information <- fit_L0_information[order(fit_L0_information$suppSize, decreasing = TRUE), ]
        lambda_L0 <- fit_L0_information$lambda[1]
        gamma_L0 <- fit_L0_information$gamma[1]
        
        y_cat = predict(fit_L0,
                        newx=X,
                        lambda=lambda_L0,
                        gamma=gamma_L0)
        
        y_hat=as.vector(y_cat)
        # plot(fit_L0, gamma = gamma_L0, showLines=TRUE)
        
        temp = coef(fit_L0,
                    lambda=lambda_L0,
                    gamma=gamma_L0)
        
        temp = as.vector(temp)
        temp = temp[-1] # 排除第一个位置上的intercept
        temp = which(temp!=0) # 排除系数为0的冗余特征
        temp = colnames(X)[temp]
        temp <- na.omit(temp)
        # temp
        
        if (length(temp) == 1) {
          X_Y=cbind(X[, temp], Y)
          colnames(X_Y)[1] <- temp
          X_Y_frame = as.data.frame(X_Y)
        }else{
          X_Y=cbind(X[, temp],Y)
          X_Y_frame = as.data.frame(X_Y)
        }
        
        lmfit = lm(Y~., data=X_Y_frame)
        fit_temp=summary(lmfit)
        res_data <- as.matrix(fit_temp$coefficients)
        # res_data
        
        # L0_Rsquare = 1-mse(Y,y_hat)/var(Y)
        # ####评价指标####
        # # evaluate_L0 <- data.frame(L0_Rsquare = 1-mse(Y,y_hat)/var(Y), L0_RMSE = RMSE(Y,y_hat))
        # evaluate_L0 <- data.frame(Gene = target_gene,
        #                           L0_Rsquare = 1-mse(Y,y_hat)/var(Y), 
        #                           L0_RMSE = RMSE(Y,y_hat))
        # evaluate <- rbind(evaluate, evaluate_L0)
        
        ####只保存P小于0.05的结果####
        res_data_f <- as.data.frame(res_data)
        res_data_f <- res_data_f[which(res_data_f$`Pr(>|t|)` <= 0.05),]
        
        if (nrow(res_data_f)>0) {
          
          if (rownames(res_data_f)[1]=="(Intercept)") {
            res_data_f <- res_data_f[-1,]
          }
          
          res_data_f$strength <- abs(res_data_f$Estimate)
          
          for (k in 1:nrow(res_data_f)) {
            
            # res_data_f$strength[k] <- abs(res_data_f$Estimate[k])/sum(abs(res_data_f$Estimate))
            if (res_data_f$Estimate[k] < 0) {
              res_data_f$reg[k] <- '2'
            }else{
              res_data_f$reg[k] <- '1'
            }
          }
          
        }
        
        TRN_L0 <- cbind("TF"=row.names(res_data_f),
                        "Gene"=target_gene,
                        "Strength"=res_data_f$strength) %>% as.data.frame()
        
        if (i>2) {
          if (ncol(TRN_L0_all)==ncol(TRN_L0)) {
            TRN_L0_all <- rbind.data.frame(TRN_L0_all, TRN_L0)
          }
        }else{
          TRN_L0_all <- rbind.data.frame(TRN_L0_all, TRN_L0)
        }
        
      }
      
    }
    TRN_L0_all$Strength <- as.numeric(TRN_L0_all$Strength)
    TRN_L0_all <- TRN_L0_all[order(-TRN_L0_all$Strength),]
    TRN_L0_all_list[[d]] <- TRN_L0_all
    write.table(TRN_L0_all, 
                paste0("DREAM5_network_inference_challenge/Evaluation scripts/INPUT/predictions/myteam/L0/DREAM5_NetworkInference_myteam_Network",
                       d,
                       ".txt"),
                col.names = F,
                row.names = F,
                sep = "\t",
                quote = F)
  # }
  
  
}

