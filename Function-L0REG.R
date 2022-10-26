

LO_fit <- function(X, Y,
                   penalty = penalty,
                   nFolds = 10, seed = 1,
                   maxSuppSize = 20,
                   nGamma = 5,
                   gammaMin = 0.0001, gammaMax = 10) {
  fit <- L0Learn.cvfit(X, Y,
    penalty = penalty,
    nFolds = 10, seed = 1,
    maxSuppSize = dim(matrix)[2],
    nGamma = 5,
    gammaMin = 0.0001, gammaMax = 10
  )
  fit_inf <- print(fit)
  optimalGammaIndex <- which(unlist(lapply(fit$cvMeans, min)) == min(unlist(lapply(fit$cvMeans, min))))
  gamma <- fit$fit$gamma[optimalGammaIndex]
  lambda_list <- fit_inf[which(fit_inf$gamma == gamma), ]
  if (is.null(maxSuppSize)) {
    lambda <- min(lambda_list$lambda)
  } else {
    if (maxSuppSize %in% lambda_list$maxSuppSize) {
      lambda <- lambda_list$maxSuppSize[which(lambda_list$maxSuppSize == maxSuppSize)]
    } else {
      lambda <- min(lambda_list$lambda)
    }
  }
  temp <- coef(fit, lambda = lambda, gamma = gamma)
}

L0REG <- function(matrix,
                  targets = rownames(matrix),
                  penalty = penalty,
                  regulators = NULL,
                  maxSuppSize = dim(matrix)[2]) {
  library(L0Learn)
  matrix <- as.data.frame(t(matrix))
  if (!is.null(regulators)) {
    weightdf <- c()
    for (i in 1:length(regulators)) {
      X <- as.matrix(matrix[, -which(colnames(matrix) == regulators[i])])
      Y <- matrix[, regulators[i]]
      temp <- LO_fit(X, Y,
        penalty = penalty,
        nFolds = 10, seed = 1,
        maxSuppSize = dim(matrix)[2],
        nGamma = 5,
        gammaMin = 0.0001, gammaMax = 10
      )
      temp <- as.vector(temp)

      if (F) {
        fit <- L0Learn.fit(X, Y,
          penalty = penalty,
          maxSuppSize = dim(matrix)[2]
        )
        fit_inf <- as.data.frame(print(fit))
        fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
        lambda <- fit_inf$lambda[1]
        gamma <- fit_inf$gamma[1]
        # lambda <- fit_inf$lambda[ceiling(nrow(fit_inf))]
        # gamma <- fit_inf$gamma[ceiling(nrow(fit_inf))]
        temp <- coef(fit,
          lambda = lambda,
          gamma = gamma
        )
        temp <- as.vector(temp)
      }

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
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = 0)
        # weightd <- data.frame(regulatoryGene = regulators[i], targetGene = colnames(X), weight = 0)
      } else {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = wghts)
        # weightd <- data.frame(regulatoryGene = regulators[i], targetGene = colnames(X), weight = wghts)
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
      fit <- L0Learn.fit(X, Y,
        penalty = penalty,
        maxSuppSize = dim(matrix)[2]
      )

      fit_inf <- as.data.frame(print(fit))
      fit_inf <- fit_inf[order(fit_inf$maxSuppSize, decreasing = TRUE), ]
      lambda <- fit_inf$lambda[1]
      gamma <- fit_inf$gamma[1]
      # lambda <- fit_inf$lambda[ceiling(nrow(fit_inf))]
      # gamma <- fit_inf$gamma[ceiling(nrow(fit_inf))]
      temp <- coef(fit,
        lambda = lambda,
        gamma = gamma
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
        fit <- L0Learn.fit(X, Y,
          penalty = penalty,
          maxSuppSize = dim(matrix)[2]
        )
        fit_inf <- as.data.frame(print(fit))
        fit_inf <- fit_inf[order(fit_inf$maxSuppSize, decreasing = TRUE), ]
        lambda <- fit_inf$lambda[1]
        gamma <- fit_inf$gamma[1]
        temp <- coef(fit,
          lambda = lambda,
          gamma = gamma
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
