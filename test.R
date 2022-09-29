


# Initialization
single_cell_data <- t(expression_dataset_test)
time <- PseudoTime[rownames(single_cell_data),]
time <- time$PseudoTime %>% as.numeric()
numGENES <- dim(single_cell_data)[1]
num_time_points <- length(time)

# Distribution Distance
DISTANCE_matrix <- matrix(data = 0, nrow = num_time_points - 1, ncol = numGENES)
# h <- matrix(data=0,nrow=num_time_points-1,ncol = numGENES)
totalDATA <- single_cell_data

cmtest2 <- dget("SINCERITIES functions/cmtest2.R")

distance=1
for (ti in 1:(num_time_points - 1)) {
    totalDATA <- rbind(totalDATA, single_cell_data[[ti + 1]])
    data_ti <- t(single_cell_data[[ti]])
    data_ti_plus1 <- t(single_cell_data[[ti + 1]])

    for (gi in 1:numGENES) {
        p1 <- data_ti[gi, ]
        p2 <- data_ti_plus1[gi, ]
        if (distance == 1) {
            test.stat <- ks.test(p1, p2)
            DISTANCE_matrix[ti, gi] <- test.stat$statistic
        } else if (distance == 2) {
            DISTANCE_matrix[ti, gi] <- cmtest2(p1, p2)$CM_limiting_stat
        } else if (distance == 3) {
            test.stat <- ad.test(p1, p2)
            DISTANCE_matrix[ti, gi] <- test.stat$ad[2, 1]
        } else if (distance == 4) {
            DISTANCE_matrix[ti, gi] <- abs(mean(p1) - mean(p2))
        }
    }
}

# normalization
deltaT <- replicate(dim(DISTANCE_matrix)[2], time[2:length(time)] - time[1:(length(time) - 1)])
DISTANCE_matrix_normed <- DISTANCE_matrix / deltaT

# Generate Y and X_matrix for glmnet
if (method == 1) {
    alphas <- 0
} else if (method == 2) {
    alphas <- seq(0, 1, 0.1)
} else if (method == 3) {
    alphas <- 1
} else if (method == "L0") {
    alphas <- 0
    message("Choose L0 model!")
} else {
    input <- readline(" *** Please input manually the alpha values (between 0 and 1) separated by comma: ")
    alphas <- as.numeric(unlist(strsplit(input, ",")))
}
DISTANCE_matrix <- DISTANCE_matrix_normed
X_matrix <- DISTANCE_matrix[1:(num_time_points - 2), ]

# LOOCV settings
nfold <- dim(X_matrix)[1]
foldid <- 1:nfold
keep <- TRUE
pred_lambda_min <- matrix(0, nrow = numGENES, ncol = numGENES)

lambda_res <- vector()
alpha_res <- vector()

for (gi in 1:numGENES) {
    lambda <- vector()
    cvERROR <- vector()
    beta <- matrix(data = 0, nrow = dim(X_matrix)[2], ncol = length(alphas))

    # for (test in 1:length(alphas)) {
    #   Y_vector <- DISTANCE_matrix[2:(num_time_points-1),gi]
    #   if(noDIAG==1){
    #     CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],exclude=gi,nfolds = nfold, foldid = foldid,
    #                             keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
    #   }else{
    #     CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],nfolds = nfold, foldid = foldid,
    #                             keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
    #   }
    #   lambda[test] <- CV_results$lambda.min
    #   cvERROR[test] <- CV_results$cvm[CV_results$lambda==CV_results$lambda.min]
    #   # coef.CV_results <- coef.cv.glmnet(CV_results,s='lambda.min') # glmnet::
    #   coef.CV_results <- coef(CV_results,s='lambda.min') # glmnet::
    #   beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
    # }
    #
    # minIdx <- max(which(cvERROR==min(cvERROR)))
    # lambda_res[gi] <- lambda[minIdx]
    # alpha_res[gi] <- alphas[minIdx]
    # pred_lambda_min[,gi] <- beta[,minIdx]

    for (test in 1:length(alphas)) {
        Y_vector <- DISTANCE_matrix[2:(num_time_points - 1), gi]
        L0_Model <- L0Learn.fit(X_matrix, Y_vector,
            penalty = "L0L2",
            algorithm = "CD", # "CD" and "CDPSI"
            maxSuppSize = ncol(X_matrix)
        )

        L0_Model_Information <- as.data.frame(print(L0_Model))
        L0_Model_Information <- L0_Model_Information[order(L0_Model_Information$suppSize,
            decreasing = TRUE
        ), ]
        lambda_L0 <- L0_Model_Information$lambda[1]
        gamma_L0 <- L0_Model_Information$gamma[1]
        temp <- coef(L0_Model,
            lambda = lambda_L0,
            gamma = gamma_L0
        )
        temp <- as.vector(temp)
        wghts <- temp[-1]
        # wghts <- abs(wghts)
    }
    pred_lambda_min[, gi] <- t(wghts)
}

if (SIGN == 1) {
    parcorr_matrix <- pcor(DATA$totDATA, method = "spearman")$estimate # pearson
    pred_lambda_min <- pred_lambda_min * sign(parcorr_matrix)
}

result <- list(DISTANCE_matrix = DISTANCE_matrix, adj_matrix = pred_lambda_min)
return(result)
