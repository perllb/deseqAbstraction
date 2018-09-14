#' @name GO_geneFC
#' @description Do GO term enrichement and analysis for a list of genes
#' @param dds: deseq object
#' @param species: only human (hsa) annotation currently supported, mm10 coming soon..
#' @param BP: test for biological process terms
#' @param MF: test for molecular function terms
#' @param CC: test for cellular compartment terms
#' @param sameDir: test one-directional, terms enriched with genes changed in only same direction 
#' @title GOanalysis in R - GO terms
#' @export GO_geneFC
#' @examples
#' dabs <- deseqAbs$new(name="drugTest",colData=colDat)
#' GO_pathview(dabs = dabs,species = "hsa",BP=T,MF=F,CC=F,sameDir=T)


GO_geneFC <- function(dabs=NULL,species="hsa",BP=T,MF=F,CC=F,sameDir=T) {
  
  ### Dependencies:
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("pathview")
  #biocLite("gage")
  #biocLite("gageData")
  #if(species=="hg38"){ biocLite("org.Hs.eg.db") }
  
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  library(pathview)
  library(gage)
  library(gageData)
  
  # make diffex if not done
  if(is.null(dabs$test$Default)) {
    dabs$makeDiffex()
  }
  res <- dabs$test$Default
  
  # convert symbol to entrez and name
  res$entrez = mapIds(org.Hs.eg.db,
                      keys=row.names(res), 
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
  res$name =   mapIds(org.Hs.eg.db,
                      keys=row.names(res), 
                      column="GENENAME",
                      keytype="SYMBOL",
                      multiVals="first")
  
  
  #### GO terms
  fc <- res$log2FoldChange
  names(fc) <- res$entrez
  
  data(go.sets.hs)
  data(go.subs.hs)
  
  if(!dir.exists(paths = "GO/GageGO")){
    dir.create("GO/GageGO")
  }
  
  if(BP){
    print("Testing for enrichment of Biological Process (BP).. \n")
    goSets = go.sets.hs[go.subs.hs$BP]
    goRes = gage(fc, gsets=goSets, same.dir = F)
    
    if(sameDir) {
      #write less
      write.table(x = goRes$less,file = "GO/GageGo/goBP_less.SameDir.txt",quote = F,sep = "\t",row.names = T)
    }
    
    outfile <- ifelse(sameDir,yes =  "GO/GageGo/goBP_greater.SameDir.txt",no="GO/GageGo/goBP_greater.nSame.txt")
    # write greater
    write.table(x = goRes$greater,file =outfile,quote = F,sep = "\t",row.names = T)
    
    
    
  }
  if(MF){
    print("Testing for enrichment of Molecular Functions (MF) ..\n")
    goMFsets = go.sets.hs[go.subs.hs$MF]
    goMFres = gage(fc, gsets=goMFsets, same.dir=TRUE)
    if(sameDir) {
      #write less
      write.table(x = goRes$less,file = "GO/GageGo/goMF_less.SameDir.txt",quote = F,sep = "\t",row.names = T)
    }
    
    outfile <- ifelse(sameDir,yes =  "GO/GageGo/goMF_greater.SameDir.txt",no="GO/GageGo/goMF_greater.nSame.txt")
    # write greater
    write.table(x = goRes$greater,file =outfile,quote = F,sep = "\t",row.names = T)
    
  }
  if(CC){
    print("Testing for enrichment of Cellular Compartment (CC) ..\n")
    goCCsets = go.sets.hs[go.subs.hs$CC]
    goCCres = gage(fc, gsets=goCCsets, same.dir=TRUE)
    if(sameDir) {
      #write less
      write.table(x = goRes$less,file = "GO/GageGo/goCC_less.SameDir.txt",quote = F,sep = "\t",row.names = T)
    }
    
    outfile <- ifelse(sameDir,yes =  "GO/GageGo/goCC_greater.SameDir.txt",no="GO/GageGo/goCC_greater.nSame.txt")
    # write greater
    write.table(x = goRes$greater,file =outfile,quote = F,sep = "\t",row.names = T)
    
  }
  
  
}