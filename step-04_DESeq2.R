

# Load libraries
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

counts <- scRNA_harmony@assays$RNA@counts 

metadata <- scRNA_harmony@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(scRNA_harmony@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "orig.ident")]

# Explore the raw counts for the dataset

## Check the assays present
assays(sce)

## Explore the raw counts for the dataset
dim(counts(sce))

counts(sce)[1:6, 1:6]

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(as.factor(sce$orig.ident)))

# Total number of samples 
ns <- length(sids)
ns


# Generate sample level metadata

## Determine the number of cells per sample
table(sce$orig.ident)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$orig.ident))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$orig.ident)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")
ei

# Perform QC if not already performed
dim(sce)

# Calculate quality control (QC) metrics
sce <- calculateQCMetrics(sce)

# Get cells w/ few/many detected genes
sce$is_outlier <- isOutlier(
  metric = sce$total_features_by_counts,
  nmads = 2, type = "both", log = TRUE)

# Remove outlier cells
sce <- sce[, !sce$is_outlier]
dim(sce)

## Remove lowly expressed genes which have less than 10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

dim(sce)



# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "orig.ident")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)

dim(pb)

pb[1:6, 1:6]

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb)

# Explore the different components of list
str(pb)

options(width = 100)
table(sce$cluster_id, sce$orig.ident)


# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("orig.ident", "main_cell_type")]) 


metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id) 

metadata 

# Generate vector of cluster IDs
clusters <- levels(metadata$cluster_id)
clusters

# Subset the metadata to only the B cells
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$orig.ident
head(cluster_metadata)

# Subset the counts to only the B cells
counts <- pb[[clusters[1]]]

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))  


dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group_id)



