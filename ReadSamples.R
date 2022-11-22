

library(Seurat)
# Set QC values according Paper or your need
nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 30
pHB_lower <- 0
pHB_upper <- 5

pathSamples <- getwd() # Add the pathway of sample files
samplesList <- list.files(pathSamples, pattern = "*.txt.gz")
seuratObjList <- list()
for (i in 1:length(samplesList)) {
    seuratObj <- read.table(paste0(samplesList[i]),
        row.names = 1,
        header = T,
        sep = "\t"
    )
    seuratObj <- CreateSeuratObject(
        counts = seuratObj,
        project = samplesList[i],
        min.features = 200,
        min.cells = 3
    )
    ### QC
    seuratObj <- PercentageFeatureSet(seuratObj, pattern = "^MT-", col.name = "pMT")
    seuratObj <- PercentageFeatureSet(seuratObj, pattern = "^HBA|^HBB", col.name = "pHB")
    seuratObj <- PercentageFeatureSet(seuratObj, pattern = "^RPS|^RPL", col.name = "pRP")
    seuratObj <- subset(seuratObj,
        subset = nFeature_RNA > nFeature_lower &
            nFeature_RNA < nFeature_upper &
            nCount_RNA > nCount_lower &
            nCount_RNA < nCount_upper &
            pMT < pMT_upper &
            pHB < pHB_upper
    )
    seuratObjList[[i]] <- seuratObj
}
seuratObjMer <- merge(seuratObjList[[1]],
    y = c(
        seuratObjList[2:length(seuratObjList)]
    )
)
