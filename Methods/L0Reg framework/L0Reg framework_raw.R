

rm(list=ls())

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
library(ggpubr)

load('data/luad-rsem-count-tcga-t.Rdata')
load('data/luad-rsem-fpkm-tcga-t_normlized.Rdata')
#load('../data/luad-rsem-fpkm-tcga-t_unnormlized.Rdata')
#load('../data/luad-rsem-TPM-tcga-t_normlized.Rdata')
#load('../data/luad-rsem-TPM-tcga-t_unnormlized.Rdata')
raw_tcga <- scale(raw_tcga)
# raw_tcga <- t(scale(t(raw_tcga)))
# raw_tcga[raw_tcga >= 2] = 2
# raw_tcga[raw_tcga <= 0] = 0

sample_label <- read.table("../results/2191/NMF/cluster-nk-ligands-k=8/sample_cluster_4.csv",
                           header = FALSE,
                           sep=',',
                           check.names=FALSE)
table(sample_label$V2)
deg.data_clusters <- c()

evaluate_all <- c()
res_TRN_all <- c()

for (j in 1:4) {
  
  survival_tfs <- c()
  
  res_TRN_cluster <- c()
  
  deg.data_cluster <- c()
  
  if (dir.exists(paste0('results/cluster', j)) == F) {
    dir.create(file.path(paste0('results/cluster', j)))
  }

  print(paste("------------cluster", j, "start------------", sep = ' '))
  
  if (T) {
    genes_list <- read.table(paste("data/genes_list", ".txt", sep = ""),
                             header = TRUE,
                             row.names = 1)
  }else{
    genes_list <- read.table(paste("data/genes_list_cluster", j, ".txt", sep = ""),
                             header = TRUE,
                             row.names = 1)
  }
  
  cluster <- sample_label[which(sample_label$V2 == j), 1]
  
  Y_raw <- raw_tcga[row.names(genes_list), cluster]
  
  maxSNVSize <- 100
  evaluate_cluster <- c()
  for(i in 1:nrow(genes_list)){
    
    target_gene <- row.names(genes_list)[i]
    if (file.exists(paste('data/', target_gene, '_TFs_list.txt', sep = '')) == T) {
      
      print(paste('---------- the cluster', j, ',', 'number', i, 'gene:', target_gene, 'computation is start ! ----------', sep = ' '))
      
      TFs_list  <- read.table(paste('data/', target_gene, '_TFs_list.txt', sep = ''),
                              header = T,
                              row.names = 1)
      
      X <- raw_tcga[row.names(TFs_list), cluster]
      X <- t(X) %>% as.matrix()
      row.names(X) = row.names(X)
      
      Y <- raw_tcga[target_gene, cluster]
      Y <- t(Y) %>% as.numeric()
      
      # Validation --------------------------------------------------------------
      set.seed(2022)
      train_idx <- sample(nrow(X), 0.7 * nrow(X))
      train_data_X <- X[train_idx, ]
      train_data_Y <- Y[train_idx]
      test_data_X <- X[-train_idx, ]
      test_data_Y <- Y[-train_idx]
      
      cvfit_L0_validation <- L0Learn.cvfit(train_data_X, train_data_Y, nFolds = 10)
      plot(cvfit_L0_validation)
      fit_L0_validation <- L0Learn.fit(train_data_X, train_data_Y,
                            penalty="L0",
                            maxSuppSize = maxSNVSize)
      plot(fit_L0_validation)
      
      print(fit_L0_validation)
      
      # Extract coefficient at middle lambda
      fit_L0_information <- as.data.frame(print(fit_L0_validation))
      fit_L0_information <- fit_L0_information[order(fit_L0_information$suppSize, decreasing = TRUE), ]
      lambda_L0 <- fit_L0_information$lambda[1]
      gamma_L0 <- fit_L0_information$gamma[1]
      
      test_data_y_cat = predict(fit_L0_validation,
                      newx=test_data_X,
                      lambda=lambda_L0,
                      gamma=gamma_L0)
      test_data_y_hat=as.vector(test_data_y_cat)
      ###
      
      re_status <- cbind(test_data_Y, test_data_y_hat)
      head(re_status)
      re_status <- as.data.frame(re_status)
      colnames(re_status) <- c('raw', 'pre')
      
      peak_corr <- corr.test(re_status$raw,
                             re_status$pre,
                             # method = "spearman",
                             adjust = "none")
      
      deg.data <- data.frame(Cluster = paste('Cluster',j),
                             Gene = target_gene,
                             corr = peak_corr$r,
                             pval = peak_corr$p)
      deg.data_cluster <- rbind.data.frame(deg.data_cluster, deg.data)
      
      
      cor(re_status)
      ggplot(data = re_status, mapping = aes(x = raw, y = pre)) + 
        geom_point()+
        # stat_bin_2d(binwidth = c(.012,.012))+
        # guides(fill = guide_colorbar(title = "Counts", title.position = "top",title.hjust = .5,ticks = T))+
        # scale_fill_gradientn(colours = rev(palette),limits=c(0,4000), breaks=c(0,1000,2000,3000,4000),
        #                      labels=c("0","1000",'2000','3000','>4000'))+ 
        geom_smooth(method = 'lm', se = F, color = 'red')+
        theme_bw()+
        stat_cor(data=re_status)+ #, method = "spearman"
        #对角线
        # geom_abline(slope = 1, intercept = 0, color='#B03060', linetype = "dashed", size=1)+
        #拟合线
        geom_smooth(method = 'lm',se = F,color='#006699',size=1)+
        labs(x ='True Expression',y="L0Reg framework",
             title = paste0("Gene: ",target_gene))
      
      ggsave(paste("results/cluster", j,"/",target_gene,'_corrgram.png',sep = ""),
             width = 3,
             height = 3)
      
      # ggpairs(re_status,
      #         title = target_gene ,
      #         upper = list(continuous = wrap("cor",method="spearman")))

      if (peak_corr$p < 1 && peak_corr$r > 0) {
        
        ####L0####
        fit_L0 <- L0Learn.fit(X, Y,
                              penalty="L0",
                              maxSuppSize = maxSNVSize)
        print(fit_L0)
        # Extract the coefficients at middle lambda
        # lambda_value <- as.data.frame(fit_L0$lambda)
        # names(lambda_value) <- 'lambda'
        # lambda_value[nrow(lambda_value)/2,]
        
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
        plot(fit_L0, gamma = gamma_L0, showLines=TRUE)
        
        temp = coef(fit_L0,
                    lambda=lambda_L0,
                    gamma=gamma_L0)
        temp = as.vector(temp)
        temp = temp[-1] # 排除第一个位置上的intercept
        temp = which(temp!=0) # 排除系数为0的冗余特征
        temp = colnames(X)[temp]
        temp
        
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
        
        ####打印模型信息####
        # print(tidy(lmfit))
        
        ####ElasticNet####
        # # Build the model using the training set
        # set.seed(2022)
        # ElasticNet_model <- train(X, Y, method = "glmnet")
        # # Best tuning parameter
        # ElasticNet_model$bestTune
        # # Coefficient of the final model. You need to specify the best lambda
        # coef(ElasticNet_model$finalModel, ElasticNet_model$bestTune$lambda)
        # print(ElasticNet_model)
        # # Extract the coefficients at middle lambda
        # ElasticNet_cat = predict(ElasticNet_model,
        #                          newx=X,
        #                          lambda=ElasticNet_model$bestTune$lambda,
        #                          gamma=0)
        # ElasticNet_hat=as.vector(ElasticNet_cat)
        # plot(ElasticNet_model)
        # ElasticNet_temp = coef(ElasticNet_model$finalModel, ElasticNet_model$bestTune$lambda)
        # ElasticNet_temp = as.vector(ElasticNet_temp)
        # ElasticNet_temp = ElasticNet_temp[-1]
        # ElasticNet_temp = which(ElasticNet_temp!=0)
        # ElasticNet_temp = colnames(X)[ElasticNet_temp]
        # ElasticNet_temp
        # if (length(ElasticNet_temp) == 1) {
        #   ElasticNet_X_Y=cbind(X[, ElasticNet_temp], Y)
        #   colnames(ElasticNet_X_Y)[1] <- ElasticNet_temp
        #   ElasticNet_X_Y_frame = as.data.frame(ElasticNet_X_Y)
        # }else{
        #   ElasticNet_X_Y=cbind(X[,ElasticNet_temp],Y)
        #   ElasticNet_X_Y_frame = as.data.frame(ElasticNet_X_Y)
        # }
        # ElasticNet_lmfit = lm(Y~.,data=ElasticNet_X_Y_frame)
        # fit_ElasticNet_temp=summary(ElasticNet_lmfit)
        # ElasticNet_res <- as.matrix(fit_ElasticNet_temp$coefficients)
        # ElasticNet_res
        # 
        # ####打印模型信息####
        # print(tidy(ElasticNet_lmfit))
        # 
        L0_Rsquare = 1-mse(Y,y_hat)/var(Y)
        ####评价指标####
        # evaluate_L0 <- data.frame(L0_Rsquare = 1-mse(Y,y_hat)/var(Y), L0_RMSE = RMSE(Y,y_hat))
        evaluate_L0 <- data.frame(Cluster = paste0('Cluster',j),
                                  Gene = target_gene,
                                  L0_Rsquare = 1-mse(Y,y_hat)/var(Y), 
                                  L0_RMSE = RMSE(Y,y_hat))
        # evaluate_ElasticNet <- data.frame(ElasticNet_Rsquare = 1-mse(Y,ElasticNet_hat)/var(Y), ElasticNet_RMSE = RMSE(Y,ElasticNet_hat))
        # evaluate_L0 <- data.frame(L0_Rsquare = R2(Y,y_hat),L0_RMSE = RMSE(Y,y_hat))
        # evaluate_ElasticNet <- data.frame(ElasticNet_Rsquare = R2(Y,ElasticNet_hat),ElasticNet_RMSE = RMSE(Y,ElasticNet_hat))
        ####保存模型的评价指标####
        # evaluate <- data.frame(R_square,res_rmse)
        # evaluate <- cbind(evaluate_L0,evaluate_ElasticNet)
        evaluate <- evaluate_L0
        # row.names(evaluate) <- paste('cluster',j, target_gene,sep = '_')
        evaluate_cluster <- rbind(evaluate_cluster, evaluate)
        
        ####只保存P小于0.05的结果####
        res_data_f <- as.data.frame(res_data)
        res_data_f <- res_data_f[which(res_data_f$`Pr(>|t|)` <= 0.05),]
        
        # ElasticNet_res_f <- as.data.frame(ElasticNet_res)
        # ElasticNet_res_f <- ElasticNet_res[which(ElasticNet_res_f$`Pr(>|t|)` <= 0.05),]
        
        if (nrow(res_data_f)>0) {
          for (k in 1:nrow(res_data_f)) {
            if (res_data_f$Estimate[k] < 0) {
              res_data_f$reg[k] <- '2'
            }else
              res_data_f$reg[k] <- '1'
          }
          
          write.csv(res_data_f,sprintf(paste0("results/cluster",j,'/%s_TFs_Selected_L0.csv'),target_gene))
          # write.csv(ElasticNet_res_f,sprintf('%s_TFs_Selected_Selected_ElasticNet.csv',target_gene))
          
          res_TRN <- cbind(row.names(res_data_f),res_data_f$reg, target_gene, paste('cluster', j, sep = ''))
          #res_TRN <- res_TRN[-1, ]
          res_TRN_cluster <- rbind(res_TRN_cluster, res_TRN)
          survival_tfs_gene <- rownames(res_data_f) %>% as.data.frame()
          survival_tfs <- rbind.data.frame(survival_tfs, survival_tfs_gene, target_gene)
          
        }
      }
      
      print(paste('---------- The cluster',j,',','number',i,'gene:', target_gene,'computation is completed ! ----------',sep = ' '))
    }else{
      print(paste('---------- No', target_gene,'TFs file ! ----------',sep = ' '))
    }
  }
  
  deg.data_clusters <- rbind.data.frame(deg.data_clusters, deg.data_cluster)
  
  write.csv(evaluate_cluster,sprintf('results/cluster%s_evaluate.csv',j))
  
  evaluate_all <- rbind(evaluate_all, evaluate_cluster)
  res_TRN_cluster <- as.data.frame(res_TRN_cluster)
  names(res_TRN_cluster) <- c('TFs','reg', 'gene', 'cluster')
  write.csv(res_TRN_cluster,sprintf('results/res_TRN_cluster%s.csv',j))
  res_TRN_all <- rbind(res_TRN_all, res_TRN_cluster)
  survival_tfs <- survival_tfs[!duplicated(survival_tfs),] %>% as.data.frame()
  write.csv(survival_tfs,paste0("../SurvivalAnalysis/genes_list_cluster", j, ".csv"))
  print(paste('---------- cluster',j,'done ! ----------', sep = ' '))
}

res_TRN_all <- res_TRN_all[-which(res_TRN_all$TFs=='(Intercept)'),]
write.csv(res_TRN_all, 'results/res_TRN_all.csv',row.names = T)
write.csv(evaluate_all, 'results/evaluate_all.csv',row.names = F)

# write.csv(deg.data_clusters, "results/deg.data_clusters.csv")

