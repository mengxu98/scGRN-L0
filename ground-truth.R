

# intput = "output/output_net_hESC_L0.txt"
# dataset = "STRING"
# output = "output/"
# dataset_dir = "../scGRN-L0_data/BEELINE-Networks/Networks/human/"

data_filetr <- function(data, dataset_dir, database = "ALL") {
  if (database %in% c("STRING", "ChIP-seq", "HepG2-ChIP-se", "hESC-ChIP-seq", "All")) {
    if (database == "STRING") {
      dataset <- read.csv(paste0(dataset_dir, "STRING-network.csv"))
    } else if (database == "ChIP-seq") {
      dataset <- read.csv(paste0(dataset_dir, "Non-specific-ChIP-seq-network.csv"))
    } else if (database == "HepG2-ChIP-seq") {
      dataset <- read.csv(paste0(dataset_dir, "HepG2-ChIP-seq-network.csv"))
    } else if (database == "hESC-ChIP-seq") {
      dataset <- read.csv(paste0(dataset_dir, "hESC-ChIP-seq-network.csv"))
    } else if (database == "All") {
      dataset1 <- read.csv(paste0(dataset_dir, "STRING-network.csv"))
      dataset2 <- read.csv(paste0(dataset_dir, "Non-specific-ChIP-seq-network.csv"))
      dataset3 <- read.csv(paste0(dataset_dir, "HepG2-ChIP-seq-network.csv"))
      dataset4 <- read.csv(paste0(dataset_dir, "hESC-ChIP-seq-network.csv"))
      dataset <- rbind.data.frame(dataset1, dataset2, dataset3, dataset4)
    }
    genes_tfs <- c(dataset$Gene1, dataset$Gene2)

    if (!exists("data")) {
      message("----- Expression matrix does not exist! Please read dataset! -----")
      message("----- The rows of the expression matrix are samples listed as genes! -----")
    } else {
      data <- data[, intersect(colnames(expression_dataset_test), genes_tfs)]
    }
  } else {
    message("----- Regulatory network dataset does not exist! -----")
    message("----- Please choose one of STRING, ChIP-seq, HepG2-ChIP-seq and hESC-ChIP-seq, or All! -----")
  }
}


ground_truth_h <- function(intput, output, dataset_dir, database = "STRING") {
  if (database %in% c("STRING", "ChIP-seq", "HepG2-ChIP-se", "hESC-ChIP-seq", "All")) {
    if (database == "STRING") {
      dataset <- read.csv(paste0(dataset_dir, "STRING-network.csv"))
    } else if (database == "ChIP-seq") {
      dataset <- read.csv(paste0(dataset_dir, "Non-specific-ChIP-seq-network.csv"))
    } else if (database == "HepG2-ChIP-seq") {
      dataset <- read.csv(paste0(dataset_dir, "HepG2-ChIP-seq-network.csv"))
    } else if (database == "hESC-ChIP-seq") {
      dataset <- read.csv(paste0(dataset_dir, "hESC-ChIP-seq-network.csv"))
    } else if (database == "All") {
      dataset1 <- read.csv(paste0(dataset_dir, "STRING-network.csv"))
      dataset2 <- read.csv(paste0(dataset_dir, "Non-specific-ChIP-seq-network.csv"))
      dataset3 <- read.csv(paste0(dataset_dir, "HepG2-ChIP-seq-network.csv"))
      dataset4 <- read.csv(paste0(dataset_dir, "hESC-ChIP-seq-network.csv"))
      dataset <- rbind.data.frame(dataset1, dataset2, dataset3, dataset4)
    }

    if (!exists("output")) {
      message("----- Output folder does not exist! Please reset! -----")
    }

    if (!exists("intput")) {
      message("----- Regulatory network does not exist! Please run NIMEFI! -----")
    } else {
      grn <- read.table(paste0(intput))
      tf_genes <- c()
      pred_tf_genes <- c()
      for (i in 1:nrow(grn)) {
        if (grn$V1[i] != grn$V2[i]) {
          dataset_list <- dataset$Gene2[which(dataset$Gene1 == grn$V1[i])]
          if (grn$V2[i] %in% dataset_list) {
            tf_gene <- data.frame(Gene1 = grn$V1[i], Gene2 = grn$V2[i], label = "1")
          } else {
            tf_gene <- data.frame(Gene1 = grn$V1[i], Gene2 = grn$V2[i], label = "0")
          }
          tf_genes <- rbind.data.frame(tf_genes, tf_gene)
          pred_tf_gene <- data.frame(Gene1 = grn$V1[i], Gene2 = grn$V2[i], rank = i)
          pred_tf_genes <- rbind.data.frame(pred_tf_genes, pred_tf_gene)
        }
      }
      write.table(pred_tf_genes,
        paste0(output, "ground_pred.txt"),
        quote = F,
        sep = "\t",
        row.names = F,
        col.names = F
      )
      write.table(tf_genes,
        paste0(output, "ground_truth.tsv"),
        quote = F,
        sep = "\t",
        row.names = F,
        col.names = F
      )
    }
  } else {
    message("----- Regulatory network dataset does not exist! -----")
    message("----- Please choose one of STRING, ChIP-seq, HepG2-ChIP-seq and hESC-ChIP-seq, or All! -----")
  }
}
# --------------------------------------------------
ground_truth_m <- function(intput, output, dataset_dir, database = "STRING") {
  if (database %in% c("STRING", "mDC-ChIP-seq", "mESC-ChIP-seq", "mESC-lofgof", "mHSC-ChIP-seq", "ChIP-seq", "All")) {
    if (database == "STRING") {
      dataset <- read.csv(paste0(dataset_dir, "STRING-network.csv"))
    } else if (database == "ChIP-seq") {
      dataset <- read.csv(paste0(dataset_dir, "Non-specific-ChIP-seq-network.csv"))
    } else if (database == "HepG2-ChIP-seq") {
      dataset <- read.csv(paste0(dataset_dir, "mDC-ChIP-seq-network.csv"))
    } else if (database == "ChIP-seq") {
      dataset <- read.csv(paste0(dataset_dir, "mESC-ChIP-seq-network.csv"))
    } else if (database == "HepG2-ChIP-seq") {
      dataset <- read.csv(paste0(dataset_dir, "mESC-lofgof-network.csv"))
    } else if (database == "hESC-ChIP-seq") {
      dataset <- read.csv(paste0(dataset_dir, "mHSC-ChIP-seq-network.csv"))
    } else if (database == "All") {
      dataset1 <- read.csv(paste0(dataset_dir, "STRING-network.csv"))
      dataset2 <- read.csv(paste0(dataset_dir, "Non-specific-ChIP-seq-network.csv"))
      dataset3 <- read.csv(paste0(dataset_dir, "mDC-ChIP-seq-network.csv"))
      dataset4 <- read.csv(paste0(dataset_dir, "mESC-ChIP-seq-network.csv"))
      dataset5 <- read.csv(paste0(dataset_dir, "mESC-lofgof-network.csv"))
      dataset6 <- read.csv(paste0(dataset_dir, "mHSC-ChIP-seq-network.csv"))
      dataset <- rbind.data.frame(dataset1, dataset2, dataset3, dataset4, dataset5, dataset6)
    }

    if (!exists("output")) {
      message("----- Output folder does not exist! Please reset! -----")
    }

    if (!exists("intput")) {
      message("----- Regulatory network does not exist! Please run NIMEFI! -----")
    } else {
      grn <- read.table(paste0(intput))
      tf_genes <- c()
      pred_tf_genes <- c()
      for (i in 1:nrow(grn)) {
        if (grn$V1[i] != grn$V2[i]) {
          dataset_list <- dataset$Gene2[which(dataset$Gene1 == grn$V1[i])]
          if (grn$V2[i] %in% dataset_list) {
            tf_gene <- data.frame(Gene1 = grn$V1[i], Gene2 = grn$V2[i], label = "1")
          } else {
            tf_gene <- data.frame(Gene1 = grn$V1[i], Gene2 = grn$V2[i], label = "0")
          }
          tf_genes <- rbind.data.frame(tf_genes, tf_gene)
          pred_tf_gene <- data.frame(Gene1 = grn$V1[i], Gene2 = grn$V2[i], rank = i)
          pred_tf_genes <- rbind.data.frame(pred_tf_genes, pred_tf_gene)
        }
      }
      write.table(pred_tf_genes,
        paste0(output, "ground_pred.txt"),
        quote = F,
        sep = "\t",
        row.names = F,
        col.names = F
      )
      write.table(tf_genes,
        paste0(output, "ground_truth.tsv"),
        quote = F,
        sep = "\t",
        row.names = F,
        col.names = F
      )
    }
  } else {
    message("----- Regulatory network dataset does not exist! -----")
    message("----- Please choose one of STRING, ChIP-seq, HepG2-ChIP-seq and hESC-ChIP-seq, or All! -----")
  }
}
# --------------------------------------------------
ground_truth_simulation <- function(intput, output, dataset_dir, file) {
  if (!is.null(dataset_dir) && !is.null(file)) {
    dataset <- read.csv(paste0(dataset_dir, file))

    if (!exists("output")) {
      message("----- Output folder does not exist! Please reset! -----")
    }

    if (!exists("intput")) {
      message("----- Regulatory network does not exist! Please run NIMEFI! -----")
    } else {
      grn <- read.table(paste0(intput))
      tf_genes <- c()
      pred_tf_genes <- c()
      for (i in 1:nrow(grn)) {
        if (grn$V1[i] != grn$V2[i]) {
          dataset_list <- dataset$Gene2[which(dataset$Gene1 == grn$V1[i])]
          if (grn$V2[i] %in% dataset_list) {
            tf_gene <- data.frame(Gene1 = grn$V1[i], Gene2 = grn$V2[i], label = "1")
          } else {
            tf_gene <- data.frame(Gene1 = grn$V1[i], Gene2 = grn$V2[i], label = "0")
          }
          tf_genes <- rbind.data.frame(tf_genes, tf_gene)
          pred_tf_gene <- data.frame(Gene1 = grn$V1[i], Gene2 = grn$V2[i], rank = i)
          pred_tf_genes <- rbind.data.frame(pred_tf_genes, pred_tf_gene)
        }
      }
      write.table(pred_tf_genes,
        paste0(output, "ground_pred.txt"),
        quote = F,
        sep = "\t",
        row.names = F,
        col.names = F
      )
      write.table(tf_genes,
        paste0(output, "ground_truth.tsv"),
        quote = F,
        sep = "\t",
        row.names = F,
        col.names = F
      )
    }
  } else {
    message("----- Regulatory network dataset does not exist! -----")
    message("----- Please choose one of STRING, ChIP-seq, HepG2-ChIP-seq and hESC-ChIP-seq, or All! -----")
  }
}


# For T cell dynamic regulation

ground_truth_T <- function(input, dataset,output=NULL) {
    if (!exists("input")) {
      message("----- Regulatory network does not exist! Please check input or dataset! -----")
    } else {
      grn <- input
      tf_genes <- c()
      pred_tf_genes <- c()
      for (i in 1:nrow(grn)) {
        if (grn$regulatoryGene[i] != grn$targetGene[i]) {
          dataset_list <- dataset$target[which(dataset$tf == grn$regulatoryGene[i])]
          if (grn$targetGene[i] %in% dataset_list) {
            tf_gene <- data.frame(Gene1 = grn$regulatoryGene[i], Gene2 = grn$targetGene[i], label = "1")
          } else {
            tf_gene <- data.frame(Gene1 = grn$regulatoryGene[i], Gene2 = grn$targetGene[i], label = "0")
          }
          tf_genes <- rbind.data.frame(tf_genes, tf_gene)
          pred_tf_gene <- data.frame(Gene1 = grn$regulatoryGene[i], Gene2 = grn$targetGene[i], rank = i)
          pred_tf_genes <- rbind.data.frame(pred_tf_genes, pred_tf_gene)
        }
      }
      write.table(pred_tf_genes,
                  paste0(output, "ground_pred.txt"),
                  quote = F,
                  sep = "\t",
                  row.names = F,
                  col.names = F
      )
      write.table(tf_genes,
                  paste0(output, "ground_truth.tsv"),
                  quote = F,
                  sep = "\t",
                  row.names = F,
                  col.names = F
      )
    }
}
