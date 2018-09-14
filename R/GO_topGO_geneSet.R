#' @name GO_topGO_geneSet
#' @description Do GO term enrichement and analysis for a list of genes
#' @param dabs: deseq object
#' @param org: only hg38 annotation currently supported, mm10 coming soon..
#' @param term: which go term ("BP" "MF" or "CC")
#' @param geneSet: list of genes for test
#' @param nodeSize: Set smallest included node size (number of genes in term) in enrichment test
#' @param outdir: name of output directory
#' @import topGO
#' @import PANTHER.db
#' @import org.Hs.eg.db
#' @import AnnotationDbi
#' @import pathview
#' @import gage
#' @import gageData
#' @import genefilter
#' @import RCurl
#' @title GOanalysis in R - topGO enrichment test terms
#' @export GO_topGO_geneSet
#' @examples
#' 
#' dabs$makeDiffex
#' geneSet <- getSignName(x = dabs$test$Default,p=0.01)$up # get upregulated genes
#' GO_topGO_geneSet(org = "hsa",BP=T,MF=F,CC=F,geneSet=geneSet,outdir=currdir,nodeSize=3)

GO_topGO_geneSet <- function(dabs=NULL,geneSet=NULL,org="hsa",term="BP",nodeSize=5,outdir=".") {
  
  if(is.null(geneSet)){
    stop("> ERROR: you need to enter genes to test..")
  }
  geneSet <- data.frame(symbol = geneSet,stringsAsFactors = F)
  
  # Get annotation mapping
  print("> Get annotation mapping symbol to entrezID")
  mapping <- read.delim(text = getURL("https://raw.githubusercontent.com/perllb/deseqabstraction/master/annotation/genenames.org_entrez.genesymbol.ensembl.txt"))
  
  print("> Mapping file downloaded.. ")
  mergeGenes <- base::merge.data.frame(x=geneSet$symbol,y=data.frame(mapping),by.x=1,by.y=2)
  # update geneSet with new mapping
  geneSet <- data.frame(symbol=mergeGenes$x,entrezid=mergeGenes$Entrez.Gene.ID,ensemblid=mergeGenes$Ensembl.ID.supplied.by.Ensembl.,stringsAsFactors = F)

  if(is.null(dabs$normCounts)){
    dabs$makeDESeq()
    if(is.null(dabs$normCounts)){
      stop("> ERROR: dabs$normCounts is NULL.. counts should be normalized when creating new deseqAbs object.. Check your data")
    }  
  }
  
  ## Filter low abundancy genes
  selProbes <- genefilter(dabs$normCounts, filterfun(pOverA(0.20, log2(40)), function(x) (IQR(x) > 0.25)))
  eset <- dabs$normCounts[selProbes, ]
  meset <- base::merge.data.frame(x=rownames(eset),y=mapping,by.x=1,by.y=2)
  eset_entrez <- meset$Entrez.Gene.ID
  
  print(paste("> Filtering low abundance reads: ",length(eset_entrez)," out of ",nrow(dabs$normCounts)," genes remain..",sep = ""))
  
  
  if(org=="hsa") {
    pthOrganisms(PANTHER.db) <- "HUMAN"
  }else {
    print("> ERROR: Currently, only human is supported")
  }
  
  # Get mapping of all entrez ids to GO terms
  allEntrez <- as.vector(eset_entrez)
  unloadNamespace("tidyverse")
  unloadNamespace("modelr")
  unloadNamespace("broom")
  unloadNamespace("dplyr")
  selection <- select(PANTHER.db,keytype = "ENTREZ",columns = c("GOSLIM_ID","GOSLIM_TERM"),keys=allEntrez)
  
  #BP
  ## Select only BP and collapse on entrez ID
  goSelection <- selection %>% 
    filter(GOSLIM_TERM==term) %>%
    group_by(ENTREZ) %>%
    summarise_all(funs(toString))
  
  #make data frame of data
  geneToGO <- data.frame(ENTREZ=goSelection$ENTREZ,GOSLIM_ID=goSelection$GOSLIM_ID)
  
  if(!dir.exists(paste(outdir,"/GO/topGO",sep = ""))){
   if(!dir.exists(paste(outdir,"/GO/",sep = ""))){
     dir.create(paste(outdir,"/GO/",sep = ""))
   }
    dir.create(paste(outdir,"/GO/topGO",sep = ""))
  }
  
  # write mapping 
  print("> Gene2GO mapping .. ")
  gene2gofile<-paste(outdir,"/GO/topGO/geneToGO_entrez-panther.",term,".txt",sep = "")
  write.table(x = geneToGO,file = gene2gofile,quote = F,sep="\t",row.names = F,col.names = F)
  
  # read mapping to correct format
  if(!file.exists(gene2gofile)){
    stop("> ERROR: No gene2go mapping file written.. ")
  }
  geneID2GO <- readMappings(file = gene2gofile)
  # inverse mappings
  GO2geneID <- inverseList(geneID2GO)
  
  print("> Defining gene lists.. ")
  ## get gene set to be tested
  geneNames <- names(geneID2GO)
  geneList <- factor(as.integer(geneNames %in% geneSet$entrez))
  names(geneList) <- geneNames
  
  print("> Create GO object")
  # create GO object
  GOdata <- new("topGOdata", ontology = term, allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=nodeSize)  
  
  print(paste("> ",numGenes(GOdata)," out of ",nrow(res)," (from deseqAbs object) genes are included in GOdata object.. for analysis ",sep = ""))  
  
  #significant genes
  print(paste("> ",numSigGenes(GOdata)," out of ",numGenes(GOdata)," (all with GO terms) genes are in significant terms ",sep = ""))  
  
  print(graph(GOdata))
  
  print("> GO Fisher test ..")
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test.stat)
  
  hist(resultFisher@score,50,xlab = "p-values")
  
  allRes <- GenTable(object = GOdata, classic = resultFisher, ranksOf = "classic", topNodes = length(resultFisher@score)) 
  print(head(allRes))
  
  return(allRes)
}


