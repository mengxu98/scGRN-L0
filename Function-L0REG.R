#' LO_fit
#'
#' @param X The rows are samples and the columns are genes of the matrix
#' @param Y 
#' @param penalty 
#' @param nFolds 
#' @param seed 
#' @param maxSuppSize 
#' @param nGamma 
#' @param gammaMin 
#' @param gammaMax 
#'
#' @return
#' @export
#'
#' @examples
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

#' L0REG
#'
#' @param matrix The rows are samples and the columns are genes of the matrix
#' @param penalty 
#' @param regulators 
#' @param targets 
#' @param maxSuppSize 
#'
#' @return
#' @export
#'
#' @examples
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
      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      # wghts <- wghts / max(wghts)
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
      # wghts <- wghts / max(wghts)
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

#' compute.gene.rank
#' 
#' Function to compute page rank of TF+target networks
#'
#' @param weightdf Result of GRN reconstruction
#' @param directedGraph If GRN is directed or not
#'
#' @return
#' @export
#'
#' @examples
compute.gene.rank <- function(weightdf, directedGraph = FALSE) {
  library("igraph")
  colnames(weightdf) <- c("regulatoryGene", "targetGene", "weight")
  tfnet <- graph_from_data_frame(weightdf, directed = directedGraph)
  pageRank <- data.frame(page_rank(tfnet, directed = directedGraph)$vector)
  colnames(pageRank) <- c("pageRank")
  pageRank$gene <- rownames(pageRank)
  pageRank <- pageRank[, c("gene", "pageRank")]
  pageRank <- pageRank[order(pageRank$pageRank, decreasing = TRUE), ]
  pageRank$is_regulator <- FALSE
  pageRank$is_regulator[pageRank$gene %in% unique(weightdf$regulatoryGene)] <- TRUE
  return(pageRank)
}

#' Title
#'
#' @param grn
#'
#' @return
#' @export
#'
#' @examples
net.format <- function(grn){
  colnames(grn) <- c("TF","TG","weight")
  grn$weight <- as.numeric(grn$weight)
  grn$interaction <- "activation"
  grn$interaction[grn$weight < 0] <- "repression"
  return(grn)
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


# Plot function --------------------------------------------------
#' Plot of dynamic networks
#'
#'
#' @param grn
#' @param tfs tfs
#' @param only_TFs plot only TF network
#' 
#' @return 
#' 
#' @export
#'
plot.dynamic.network <- function(weightdf, tfs = NULL, onlyTFs = TRUE, order = NULL, thresh = NULL){
  library("ggnetwork")
  df <- net.format(weightdf)
  net <- graph_from_data_frame(df[, c("TF", "TG", "interaction")], directed=FALSE)
  #tfnet<-ggnetwork(net,layout="fruchtermanreingold",cell.jitter=0)
  layout <- layout_with_fr(net)
  rownames(layout) <- V(net)$name
  layout_ordered <- layout[V(net)$name,]
  tfnet <- ggnetwork(net, layout = layout_ordered, cell.jitter=0)
  tfnet$is_regulator <- as.character(tfnet$name %in% tfs)
  cols<-c("activation" = "blue", "repression" = "red")
  g<-ggplot()+
    geom_edges(data=tfnet,
               aes(x=x, y=y, xend=xend, yend=yend, color=interaction),
               size=0.75,
               curvature=0.1,
               alpha=.6)+
    geom_nodes(data=tfnet,
               aes(x=x, y=y, xend=xend, yend=yend),
               color="darkgray",
               size=6,
               alpha=.5)+
    geom_nodes(data=tfnet[tfnet$is_regulator=="TRUE",],
               aes(x=x, y=y, xend=xend, yend=yend),
               color="#8C4985",
               size=6,
               alpha=.8)+
    #geom_nodelabel_repel(data=tfnet,aes(x=x, y=y, label=vertex.names),size=6, color="#8856a7")+
    scale_color_manual(values=cols)+
    geom_nodelabel_repel(data=tfnet,
                         aes(x=x, y=y, label=name),
                         size=2.5,
                         color="#5A8BAD")+
    theme_blank() #+
  #ggtitle(names(grn)[i])
  g <- g + theme(legend.position="none")
}

L0Plot <- function(data, plotType = NULL) {
  if (is.null(plotType)) {
    plotType <- boxplot
  }
  if (plotType == boxplot) {
    p <- ggplot(
      data,
      aes(x = Method, y = AUPRC)
    ) +
      # geom_violin(aes(fill = Method),
      #     trim = FALSE
      # ) +
      geom_boxplot(aes(fill = Method),
        width = 0.8
      ) +
      stat_compare_means(
        method = "wilcox.test",
        label = "p.signif",
        comparisons = my_comparisons,
        bracket.size = 0.6,
        sizen = 4,
        color = "#6699cc"
      ) +
      scale_fill_manual(values = mycol) +
      # scale_color_manual(values = mycol) +
      scale_x_discrete(labels = methods) +
      labs(x = "Methods", y = "AUPRC") +
      theme(legend.position = "bottom") +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 10
        )
      )
  }

  p
}
