

L0REG <- function(matrix,
                  targets = rownames(matrix),
                  penalty = penalty,
                  regulators = NULL,
                  suppSize=NULL) {
  library(L0Learn)
  matrix <- as.data.frame(t(matrix))
  if (!is.null(regulators)) {
    weightdf <- c()
    for (i in 1:length(regulators)) {
      X <- as.matrix(matrix[, -which(colnames(matrix) == regulators[i])])
      Y <- matrix[, regulators[i]]
      if (T) {
        fit <- L0Learn.cvfit(X, Y,
          penalty = penalty,
          nFolds = 10, seed = 1, 
          maxSuppSize = 20, 
          nGamma = 5,
          gammaMin = 0.0001, gammaMax = 10
        )
        fit_inf <- print(fit)
        optimalGammaIndex <- which(unlist(lapply(fit$cvMeans, min)) == min(unlist(lapply(fit$cvMeans, min))))
        gamma <- fit$fit$gamma[optimalGammaIndex]
        lambda_list <- fit_inf[which(fit_inf$gamma == gamma), ]

        if (is.null(suppSize)) {
          lambda <- min(lambda_list$lambda)
        } else {
          if (suppSize %in% lambda_list$suppSize) {
            lambda <- lambda_list$suppSize[which(lambda_list$suppSize == suppSize)]
          } else {
            lambda <- min(lambda_list$lambda)
          }
        }
        temp <- coef(fit, lambda = lambda, gamma = gamma)
      } else {
        L0_Model <- L0Learn.fit(X, Y,
          penalty = penalty,
          maxSuppSize = dim(matrix)[2]
        )

        L0_Model_Information <- as.data.frame(print(L0_Model))
        L0_Model_Information <- L0_Model_Information[order(L0_Model_Information$suppSize, decreasing = TRUE), ]
        lambda_L0 <- L0_Model_Information$lambda[1]
        gamma_L0 <- L0_Model_Information$gamma[1]
        # lambda_L0 <- L0_Model_Information$lambda[ceiling(nrow(L0_Model_Information))]
        # gamma_L0 <- L0_Model_Information$gamma[ceiling(nrow(L0_Model_Information))]
        temp <- coef(L0_Model,
          lambda = lambda_L0,
          gamma = gamma_L0
        )
      }

      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      # wghts <- wghts / max(wghts)

      if (F) {
        wghts <- wghts / max(wghts)

        # Now sort the wghts
        indices <- sort.list(wghts, decreasing = TRUE)
        # Check for zero entries
        zeros <- which(wghts == 0)
        # Now replace by ones that are in the top and are non-zero
        wghts[1:length(wghts)] <- 0
        wghts[indices[1:(0.25 * length(wghts))]] <- 1
        # wghts[indices[1:5]] <- 1
        # Set the ones that were zero to zero anyway
        wghts[zeros] <- 0
      }
      if (sum(wghts) == 0 & length(wghts) != nrow(matrix)) {
        # weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = 0)
        weightd <- data.frame(regulatoryGene = regulators[i], targetGene = colnames(X), weight = 0)
      } else {
        # weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = wghts)
        weightd <- data.frame(regulatoryGene = regulators[i], targetGene = colnames(X), weight = wghts)
      }
      # weightd$weight <- weightd$weight / max(weightd$weight)
      weightdf <- rbind.data.frame(weightdf, weightd)
      if (i == length(regulators)) {
        weightdf <- weightdf[order(weightdf$weight, decreasing = TRUE), ]
      }
    }
  } else {
    regulators <- rownames(matrix)
    weightdf <- c()
    for (i in 1:length(regulators)) {
      Y <- matrix[, regulators[i]]
      X <- as.matrix(matrix[, -which(colnames(matrix) == regulators[i])])
      L0_Model <- L0Learn.fit(X, Y,
        penalty = penalty,
        maxSuppSize = dim(matrix)[2]
      )

      L0_Model_Information <- as.data.frame(print(L0_Model))
      L0_Model_Information <- L0_Model_Information[order(L0_Model_Information$suppSize, decreasing = TRUE), ]
      lambda_L0 <- L0_Model_Information$lambda[1]
      gamma_L0 <- L0_Model_Information$gamma[1]
      # lambda_L0 <- L0_Model_Information$lambda[ceiling(nrow(L0_Model_Information))]
      # gamma_L0 <- L0_Model_Information$gamma[ceiling(nrow(L0_Model_Information))]
      temp <- coef(L0_Model,
        lambda = lambda_L0,
        gamma = gamma_L0
      )
      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      # wghts <- wghts / max(wghts)

      if (F) {
        wghts <- wghts / max(wghts)

        # Now sort the wghts
        indices <- sort.list(wghts, decreasing = TRUE)
        # Check for zero entries
        zeros <- which(wghts == 0)
        # Now replace by ones that are in the top and are non-zero
        wghts[1:length(wghts)] <- 0
        wghts[indices[1:(0.25 * length(wghts))]] <- 1
        # wghts[indices[1:5]] <- 1
        # Set the ones that were zero to zero anyway
        wghts[zeros] <- 0
      }

      # weightd <- data.frame(regulatoryGene = colnames(X), regulators[i], weight = wghts)
      weightd <- data.frame(regulatoryGene = regulators[i], colnames(X), weight = wghts)
      weightd$weight <- weightd$weight / max(weightd$weight)
      weightdf <- rbind.data.frame(weightdf, weightd)
      if (i == length(regulators)) {
        # weightdf$weight <- weightdf$weight / max(weightdf$weight)
        names(weightdf) <- c("regulatoryGene", "targetGene", "weight")
        weightdf <- weightdf[order(weightdf$weight, decreasing = TRUE), ]
      }
    }

    if (F) {
      weightdf <- matrix(1, dim(matrix)[2], dim(matrix)[2])
      rownames(weightdf) <- colnames(matrix)
      colnames(weightdf) <- colnames(matrix)
      for (i in 1:dim(matrix)[2]) {
        Y <- matrix[, i]
        X <- as.matrix(matrix[, -i])
        L0_Model <- L0Learn.fit(X, Y,
          penalty = penalty,
          maxSuppSize = dim(matrix)[2]
        )
        L0_Model_Information <- as.data.frame(print(L0_Model))
        L0_Model_Information <- L0_Model_Information[order(L0_Model_Information$suppSize, decreasing = TRUE), ]
        lambda_L0 <- L0_Model_Information$lambda[1]
        gamma_L0 <- L0_Model_Information$gamma[1]
        temp <- coef(L0_Model,
          lambda = lambda_L0,
          gamma = gamma_L0
        )
        temp <- as.vector(temp)
        wghts <- temp[-1]
        wghts <- wghts / max(abs(wghts))
        # if (F) {
        #   wghts <- abs(wghts)
        # } else {
        #   wghts <- abs(wghts)
        #   # Now sort the wghts
        #   indices <- sort.list(wghts, decreasing = TRUE)
        #   # Check for zero entries
        #   zeros <- which(wghts == 0)
        #
        #   # Now replace by ones that are in the top and are non-zero
        #   wghts[1:length(wghts)] <- 0
        #   wghts[indices[1:rankThreshold]] <- 1
        #
        #   # Set the ones that were zero to zero anyway
        #   wghts[zeros] <- 0
        # }
        # write result to matrix
        weightdf[colnames(matrix)[-i], colnames(matrix)[i]] <- wghts
        rownames(weightdf) <- colnames(matrix)
        colnames(weightdf) <- colnames(matrix)
      }
    }
  }
  # max(resultMatrix)
  return(weightdf)
}