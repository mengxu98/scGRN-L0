

<<<<<<< HEAD
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
=======
## Diagnostic plots for quality control
GSE131907_B <- matrix(GSE131907_B)

example_sce <- addPerCellQC(sim, 
                            subsets=list(Mito=grep("mt-", rownames(sim))))

plotColData(example_sce, x = "sum", y="detected", colour_by="Patient")
plotColData(example_sce, x = "sum", y="detected", colour_by="majorCluster")
plotColData(example_sce, x = "sum", y="detected", colour_by="sampleType")
plotColData(example_sce, x = "sum", y="subsets_Mito_percent", 
            other_fields="majorCluster") + facet_wrap(~majorCluster)
 
plotHighestExprs(example_sce, exprs_values = "counts")
 
## remove damaged cells and poorly sequenced libraries.自己制定筛选条件
 
# 查看是否有spike-ins
rownames(example_sce)[grep("^ERCC-", rownames(example_sce))]
 
## 选择总counts数大于5,000，表达的基因数大于500的细胞。
keep.total <- example_sce$sum > 5e3
#table(keep.total)
 
keep.n <- example_sce$`total mRNA mol` > 500
# table(keep.n)
 
# 根据设定的条件进行过滤
example_sce <- example_sce[,keep.total & keep.n]
#dim(example_sce)
 
## 选择至少在三个细胞中表达的基因
keep_feature <- nexprs(example_sce, byrow=TRUE) >= 3
# table(keep_feature)
example_sce <- example_sce[keep_feature,]
 
# log转化并归一化
example_sce <- logNormCounts(example_sce)
 
# plotExplanatoryVariables(example_sce)
 
vars <- getVarianceExplained(example_sce, 
                             variables=c("majorCluster", "Patient", "sampleType"))
head(vars)
class(vars)
plotExplanatoryVariables(vars)
>>>>>>> origin/master
