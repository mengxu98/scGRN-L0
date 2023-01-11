

expDat <- as.data.frame(t(data_grn))
expDat <- data_grn
t1C <- PseudoTime

genes <- rownames(expDat)
expMat <- expDat
celltime <- PseudoTime
ans <- apply(expMat[genes2,],1,function(z){
  d <- data.frame(z=z, t=celltime)
  tmp <- gam(z ~ lo(celltime), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
  print(tmp)
  print(p)
})
ans
ansd <- as.data.frame(ans)
# for slingshot
gamFit<-function(expMat,
                 genes, # genes to test
                 celltime){
  library("gam")
  genes2 <- intersect(genes, rownames(expMat))
  # could print out if any missing
  ans <- apply(expMat[genes2,names(celltime)],1,function(z){
    d <- data.frame(z=z, t=celltime)
    tmp <- gam(z ~ lo(celltime), data=d)
    
    p <- summary(tmp)[4][[1]][1,5]
    p
  })
  ans
}

#' Find genes expressed dynamically
#'
#' uses slingshot approach
#'
#'
#' @param expDat properly normalized expression matrix
#' @param sampTab sample table that includes pseudotime, rownames = cell_name, and a group column
#' @param path vector of group names to include
#' @param group_column column name in sampTab annotating groups in the path
#' @param pseudotime_column column name in sampTab annotating pseudotime or latent time
#' @param method method to find dynamic genes. Defaults to "gam", if "tradeseq", must provide expRaw
#' @param expRaw raw expression matrix, required if method="tradeseq"

#'
#' @return pvals and cell info
#' 
#' @export
#'
findDynGenes<-function(expDat, 
                       sampTab, 
                       path=NULL, 
                       group_column="dpt_groups", 
                       pseudotime_column="pseudotime",
                       method = "gam",
                       expRaw=NULL){
  
  sampTab$dpt_groups<-sampTab[,group_column]
  sampTab$pseudotime<-sampTab[,pseudotime_column]
  sampTab$cell_name<-rownames(sampTab)
  if (is.null(path)){
    path<-unique(sampTab[,group_column])
  }
  ids = vector()
  for(grp in path){
    ids = c(ids, as.vector(sampTab[sampTab$dpt_groups==grp,]$cell_name))
  }
  #assumes cell_name is rownames of sampTab and colnames of expDat
  sampTab = sampTab[ids,]
  expDat = expDat[,ids]
  if (method!="tradeseq"){
    t1 = sampTab$pseudotime
    names(t1) = as.vector(sampTab$cell_name)
    
    t1C = t1[ids]
    cat("starting gammma...\n")
    gpChr <-gamFit(expDat[,names(t1C)], rownames(expDat), t1C)
    gpChr <-gamFit(expDat, rownames(expDat), t1C)
    cells = data.frame(cell_name = names(t1), pseudotime = t1, group = as.vector(sampTab$dpt_groups))
    rownames(cells)= names(t1)
    cells = cells[order(cells$pseudotime),]
    
    ans <- list( genes = gpChr, cells = cells)
  }else{
    if (is.null(expRaw)){
      stop("Must provide expRaw for TradeSeq.")
    }
    require(tradeSeq)
    # subset raw data based on normalized data
    expRaw<-expRaw[,rownames(sampTab)]
    expRaw<-expRaw[rownames(expDat),]
    pt<-as.data.frame(sampTab[,pseudotime_column])
    rownames(pt)<-rownames(sampTab)
    colnames(pt)<-"pseudotime"
    cw<-as.matrix(rep(1,nrow(pt)))
    rownames(cw)<-rownames(sampTab)
    ts<-tradeSeq::fitGAM(as.matrix(expRaw),pseudotime=as.matrix(pt),cellWeights=cw)
    ATres<-associationTest(ts)
    genes<-ATres$pvalue
    names(genes)<-rownames(ATres)
    cells<-data.frame(cell_name = sampTab$cell_name, pseudotime = sampTab$pseudotime, group = as.vector(sampTab$dpt_groups))
    rownames(cells)<-sampTab$cell_name
    cells<-cells[order(cells$pseudotime),]
    ans<-list(genes=genes, cells=cells)
  }
  ans
}

