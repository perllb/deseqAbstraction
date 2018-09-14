#' @name GO_topGO_rank
#' @description Do GO term enrichement and analysis for a list of genes
#' @param dabs: deseq object
#' @param org: only hg38 annotation currently supported, mm10 coming soon..
#' @param term: which go term ("BP" "MF" or "CC")
#' @param nodeSize: Set smallest included node size (number of genes in term) in enrichment test
#' @param rank: How to rank genes (either "log2fc" or "padj")
#' @param sigCut: where to cut for sign genes. when rank="log2fc", it is the absolute log2FC. when rank="padj", it is p-adj value..
#' @import topGO
#' @import PANTHER.db
#' @import org.Hs.eg.db
#' @import AnnotationDbi
#' @import pathview
#' @import gage
#' @import gageData
#' @import genefilter
#' @title GOanalysis in R - topGO enrichment test terms RANK based
#' @export GO_topGO_rank
#' @examples
#' 
#' geneSet <- getSignName(x = dabs$test$Default,p=0.01)$up # get upregulated genes
#' GO_topGO_rank(dabs = dabs,org = "hsa",BP=T,MF=F,CC=F,geneSet=geneSet)



GO_topGO_rank <- function(dabs=NULL,org="hsa",term="BP",nodeSize=5,rank="log2fc",sigCut=1) {
  
  if(rank=="padj" && sigCut==1){
    print(">> WARNING: using p-adj as ranking and sigCut is set to 1 (default)! This will select all genes as significant.. better change rank to 'log2fc' or change sigCut to e.g. 0.01")
  }else{
    print(">> RANK GO term enrichment!")
    print(paste(" - Rank used   :   ",rank,sep = ""))
    print(paste(" - cutoff used :   ",sigCut,sep = ""))
  }
  
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
  
  
  BPterms <- ls(GOBPTerm)
  head(BPterms)
  
  if(is.null(dabs$normCounts)){
    stop("> ERROR: dabs$normCounts is NULL.. counts should be normalized when creating new deseqAbs object.. Check your data")
  }
  
  ## Filter low abundancy genes
  selProbes <- genefilter(dabs$normCounts, filterfun(pOverA(0.20, log2(40)), function(x) (IQR(x) > 0.25)))
  eset <- dabs$normCounts[selProbes, ]
  print(paste("> Filtering low abundance reads: ",nrow(eset)," out of ",nrow(dabs$normCounts)," genes remain..",sep = ""))
  
  ## panther annotation

  if(org=="hsa") {
    pthOrganisms(PANTHER.db) <- "HUMAN"
  }else {
    print("> ERROR: Currently, only human is supported")
  }
  
  keytypes(PANTHER.db)
  # Get mapping of all entrez ids to GO terms
  allEntrez <- as.vector(res$entrez)
  library(tidyverse)
  selection <- select(PANTHER.db,keytype = "ENTREZ",columns = c("GOSLIM_ID","GOSLIM_TERM"),keys=allEntrez)
  
  #BP
  ## Select only BP and collapse on entrez ID
  goSelection <- selection %>% 
    filter(GOSLIM_TERM==term) %>%
    group_by(ENTREZ) %>%
    summarise_each(funs(toString))
  
  #make data frame of data
  geneToGO <- data.frame(ENTREZ=goSelection$ENTREZ,GOSLIM_ID=goSelection$GOSLIM_ID)
  
  if(!dir.exists("GO/topGO")){
    dir.create("GO/topGO")
  }
  # write mapping 
  gene2gofile<-paste("GO/topGO/geneToGO_entrez-panther.",term,".txt",sep = "")
  write.table(x = geneToGO,file = gene2gofile,quote = F,sep="\t",row.names = F,col.names = F)
  
  # read mapping to correct format
  geneID2GO <- readMappings(file = gene2gofile)
  # inverse mappings
  GO2geneID <- inverseList(geneID2GO)
  
  ## get gene set to be tested
  geneNames <- names(geneID2GO)
  if(rank == "padj") {
    geneList <- res$padj
    names(geneList) <- res$entrez
    geneList[is.na(geneList)] <- 1
  } else if( rank == "log2fc") {
    geneList <- res$log2FoldChange
    names(geneList) <- res$entrez
    geneList[is.na(geneList)] <- 0
  }
  
  topDiffGenes <- function(allScore) {
    return(abs(allScore) < sigCut)
  }
  
  # create GO object
  GOdata <- new("topGOdata", ontology = term,geneSel=topDiffGenes, allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=nodeSize)  
  
  print(paste("> ",numGenes(GOdata)," out of ",nrow(res)," (from deseqAbs object) genes are included in GOdata object.. for analysis ",sep = ""))  
  
  print(graph(GOdata))
  
  sel.terms <- sample(usedGO(GOdata), 10)
  termStat(GOdata, sel.terms)  
  
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(GOdata, test.stat)
  resultWeight <- getSigGroups(GOdata, test.stat)
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test.stat)
  
  hist(resultKS@score,50,xlab = "p-values - test")
  
  
  allRes <- GenTable(GOdata, classic = resultFisher, KS = resultKS, weight = resultWeight,
                     orderBy = "weight", ranksOf = "classic", topNodes = length(resultFisher@score))
  allRes[1:5,]
  
  goID <- allRes[1, "GO.ID"]
  print(showGroupDensity(GOdata, goID, ranks = TRUE))
  return(allRes)
}





