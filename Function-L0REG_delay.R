

LO_fit <- function(X, Y,
                   penalty = penalty,
                   nFolds = 10,
                   seed = 1,
                   maxSuppSize = maxSuppSize,
                   nGamma = 5,
                   gammaMin = 0.0001, gammaMax = 10) {
  tryCatch(
    {
      fit <- L0Learn.cvfit(X, Y,
        penalty = penalty,
        maxSuppSize = maxSuppSize,
        nFolds = 10,
        seed = 1,
        nGamma = 5,
        gammaMin = 0.0001,
        gammaMax = 10
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
    },
    error = function(e) {
      fit <- L0Learn.fit(X, Y,
        penalty = penalty,
        maxSuppSize = maxSuppSize,
        nGamma = 5,
        gammaMin = 0.0001,
        gammaMax = 10
      )
      fit_inf <- print(fit)
      fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
      lambda <- fit_inf$lambda[1]
      gamma <- fit_inf$gamma[1]
      temp <- coef(fit,
        lambda = lambda,
        gamma = gamma
      )
    }
  )
}

L0REG <- function(matrix,
                  penalty = NULL,
                  regulators = NULL,
                  targets = NULL,
                  maxSuppSize = NULL) {
  library(L0Learn)
  matrix <- as.data.frame(t(matrix))
  weightdf <- c()
  if (is.null(penalty)) {
    penalty <- "L0"
  }
  if (is.null(maxSuppSize)) {
    maxSuppSize <- dim(matrix)[2]
  }
  if (is.null(targets)) {
    targets <- colnames(matrix)
  }
  if (!is.null(regulators)) {
    matrix_reg <- matrix[, regulators]
    for (i in 1:length(targets)) {
      if (targets[i] %in% regulators) {
        X <- as.matrix(matrix_reg[, -which(colnames(matrix_reg) == targets[i])])
      } else {
        X <- as.matrix(matrix_reg)
      }
      Y <- matrix[, targets[i]]
      temp <- LO_fit(X, Y,
        penalty = penalty,
        nFolds = 10,
        seed = 1,
        maxSuppSize = maxSuppSize,
        nGamma = 5,
        gammaMin = 0.0001,
        gammaMax = 10
      )
      # --------------------------------------------------
      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      # wghts <- wghts / max(wghts)
      # mean(wghts)
      # median(wghts)
      index <- which(wghts >= mean(wghts))
      X1 <- X[, index]
      if (ncol(X1) > 1) {
        temp1 <- LO_fit(X1, Y,
          penalty = penalty,
          nFolds = 10,
          seed = 1,
          maxSuppSize = maxSuppSize,
          nGamma = 5,
          gammaMin = 0.0001,
          gammaMax = 10
        )
        # --------------------------------------------------
        temp1 <- as.vector(temp1)
        wghts1 <- temp1[-1]
        wghts1 <- abs(wghts1)
        # wghts1 <- wghts1 / max(wghts)
        wghts[1:length(wghts)] <- 0
        for (ii in 1:length(index)) {
          wghts[index[ii]] <- wghts1[ii]
        }
      }
      if (F) {
        wghts <- wghts / max(wghts)
        indices <- sort.list(wghts, decreasing = TRUE)
        zeros <- which(wghts <= 0.8)
        # wghts[1:length(wghts)] <- 1
        wghts[zeros] <- 0
      }
      if (length(wghts) != dim(X)[2]) {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = targets[i], weight = 0)
        # weightd <- data.frame(regulatoryGene = targets[i], targetGene = colnames(X), weight = 0)
      } else {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = targets[i], weight = wghts)
        # weightd <- data.frame(regulatoryGene = targets[i], targetGene = colnames(X), weight = wghts)
      }
      # weightd$weight <- weightd$weight / max(weightd$weight)
      weightdf <- rbind.data.frame(weightdf, weightd)
      if (i == length(regulators)) {
        weightdf <- weightdf[order(weightdf$weight, decreasing = TRUE), ]
      }
    }
  } else {
    regulators <- colnames(matrix)
    for (i in 1:length(regulators)) {
      X <- as.matrix(matrix[, -which(colnames(matrix) == regulators[i])])
      Y <- matrix[, regulators[i]]
      temp <- LO_fit(X, Y,
        penalty = penalty,
        nFolds = 10, seed = 1,
        maxSuppSize = maxSuppSize,
        nGamma = 5,
        gammaMin = 0.0001, gammaMax = 10
      )
      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      wghts <- wghts / max(wghts)
      if (F) {
        wghts <- wghts / max(wghts)
        indices <- sort.list(wghts, decreasing = TRUE)
        zeros <- which(wghts <= 0.8)
        # wghts[1:length(wghts)] <- 1
        wghts[zeros] <- 0
      }
      if (length(wghts) != dim(X)[2]) {
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
  }
  return(weightdf)
}

# --------------------------------------------------
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

  weightdf[colnames(matrix)[-i], colnames(matrix)[i]] <- wghts
  rownames(weightdf) <- colnames(matrix)
  colnames(weightdf) <- colnames(matrix)
}
