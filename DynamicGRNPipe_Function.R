

sc_obj_filter <- function(sim){
        library(scater)
        library(scRNAseq)
        example_sce <- addPerCellQC(sim,
                                    subsets=list(Mito=grep("mt-", rownames(sim))))
        
        # rownames(example_sce)[grep("^ERCC-", rownames(example_sce))]
        # keep.total <- example_sce$sum > 5e3 # 选择总counts数大于5,000，表达的基因数大于500的细胞
        # keep.n <- example_sce$`total mRNA mol` > 500
        # example_sce <- example_sce[,keep.total & keep.n]
        # keep_feature <- nexprs(example_sce, byrow=TRUE) >= 3 # 选择至少在三个细胞中表达的基因
        # example_sce <- example_sce[keep_feature,]
        
        example_sce <- logNormCounts(example_sce)
        
        example_sce <- runPCA(example_sce)
        example_sce <- runUMAP(example_sce)
        example_sce <- runTSNE(example_sce, perplexity = 10)
        # uses the top 500 genes with the highest variances to compute the first PCs
        str(reducedDim(example_sce, "UMAP"))
        return(example_sce)
}
