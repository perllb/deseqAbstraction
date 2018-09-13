#' @name GO_topGO_geneSet
#' @description Do GO term enrichement and analysis for a list of genes
#' @param dds: deseq object
#' @param org: only hg38 annotation currently supported, mm10 coming soon..
#' @param term: which go term ("BP" "MF" or "CC")
#' @param geneSet: list of genes for test
#' @param nodeSize: Set smallest included node size (number of genes in term) in enrichment test
#' @title GOanalysis in R - topGO enrichment test terms
#' @export GO_topGO_geneSet
#' @examples
#' dabs <- deseqAbs$new(name="drugTest",colData=colDat,file=pathToFeatureCountsOutput)
#' dabs$makeDiffex
#' geneSet <- getSignName(x = dabs$test$Default,p=0.01)$up # get upregulated genes
#' GO_topGO_geneSet(dabs = dabs,org = "hsa",BP=T,MF=F,CC=F,geneSet=geneSet)

GO_topGO_geneSet <- function(dabs=NULL,geneSet=NULL,org="hsa",term="BP",nodeSize=5,outdir=".") {
  
  if(is.null(geneSet)){
    stop("> ERROR: you need to enter genes to test..")
  }
  
  geneSet <- data.frame(symbol = geneSet,stringsAsFactors = F)
  
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  # convert symbol to entrez and name
  geneSet$entrez = mapIds(org.Hs.eg.db,
                      keys=geneSet$symbol, 
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
  
  
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("topGO")
  library(topGO)
  
  BPterms <- ls(GOBPTerm)
  head(BPterms)
  
  if(is.null(dabs$normCounts)){
    dabs$makeDESeq()
    if(is.null(dabs$normCounts)){
      stop("> ERROR: dabs$normCounts is NULL.. counts should be normalized when creating new deseqAbs object.. Check your data")
    }  
  }
  
  ## Filter low abundancy genes
  library(genefilter)
  selProbes <- genefilter(dabs$normCounts, filterfun(pOverA(0.20, log2(40)), function(x) (IQR(x) > 0.25)))
  eset <- dabs$normCounts[selProbes, ]
  print(paste("> Filtering low abundance reads: ",nrow(eset)," out of ",nrow(dabs$normCounts)," genes remain..",sep = ""))
  
  ## panther annotation
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("PANTHER.db")
  library(PANTHER.db,quietly=T)
  
  if(org=="hsa") {
    pthOrganisms(PANTHER.db) <- "HUMAN"
  }else {
    print("> ERROR: Currently, only human is supported")
  }
  
  # Get mapping of all entrez ids to GO terms
  allEntrez <- as.vector(res$entrez)
  detach("package:dplyr",unload = T)
  library(tidyverse)
  selection <- select(PANTHER.db,keytype = "ENTREZ",columns = c("GOSLIM_ID","GOSLIM_TERM"),keys=allEntrez)
  
  #BP
  ## Select only BP and collapse on entrez ID
  library(dplyr)
  goSelection <- selection %>% 
    filter(GOSLIM_TERM==term) %>%
    group_by(ENTREZ) %>%
    summarise_all(funs(toString))
  
  #make data frame of data
  geneToGO <- data.frame(ENTREZ=goSelection$ENTREZ,GOSLIM_ID=goSelection$GOSLIM_ID)
  
  if(!dir.exists(paste(dir,"/GO/topGO",sep = ""))){
   if(!dir.exists(paste(dir,"/GO/",sep = ""))){
     dir.create(paste(dir,"/GO/",sep = ""))
   }
    dir.create(paste(dir,"/GO/topGO",sep = ""))
  }
  # write mapping 
  gene2gofile<-paste(dir,"/GO/topGO/geneToGO_entrez-panther.",term,".txt",sep = "")
  write.table(x = geneToGO,file = gene2gofile,quote = F,sep="\t",row.names = F,col.names = F)
  
  # read mapping to correct format
  geneID2GO <- readMappings(file = gene2gofile)
  # inverse mappings
  GO2geneID <- inverseList(geneID2GO)
  
  ## get gene set to be tested
  geneNames <- names(geneID2GO)
  geneList <- factor(as.integer(geneNames %in% geneSet$entrez))
  names(geneList) <- geneNames
  str(geneList)
  
  # create GO object
  GOdata <- new("topGOdata", ontology = term, allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=nodeSize)  
  
  genes(GOdata)
  numGenes(GOdata)  
  print(paste("> ",numGenes(GOdata)," out of ",nrow(res)," (from deseqAbs object) genes are included in GOdata object.. for analysis ",sep = ""))  
  
  #significant genes
  sg <- sigGenes(GOdata)
  numSigGenes(GOdata)
  print(paste("> ",numSigGenes(GOdata)," out of ",numGenes(GOdata)," (all with GO terms) genes are in given gene set ",sep = ""))  
  
  print(graph(GOdata))
  
  sel.terms <- sample(usedGO(GOdata), 10)
  termStat(GOdata, sel.terms)  
  
  goID <- "GO:0044255"
  gene.universe <- genes(GOdata)
  go.genes <- genesInTerm(GOdata, goID)[[1]]
  sig.genes <- sigGenes(GOdata)
  
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test.stat)
  
  hist(resultFisher@score,50,xlab = "p-values")
  
  allRes <- GenTable(object = GOdata, classic = resultFisher, ranksOf = "classic", topNodes = length(resultFisher@score)) 
  head(allRes)
  
  return(allRes)
}

library(DESeq2) 
library(devtools)
#install_github(repo = "perllb/deseqAbstraction",username = "perllb")
library(deseqAbstraction)

## FeatureCount raw file
path <- "~/Dropbox (MN)/Per/PhD/Projects/KRAB-ZNF/Quant_20180423/hg38.s2.multi.Gencode27.Exon.sjdb.txt"
## Get the sampleFile name header
header <- read.delim(path,nrows = 1,skip=1,header=T)
## Define sample names, by removing pre- and suffix 
samples <- gsub(pattern ="X.projects.fs1.medpvb.backup.projects.ChimpHuman.RNAseq.Aligned_hg38_STAR_mMap.sjdb.genc.v27.hg38.mMap.genc.v27.",replacement = "",gsub(pattern = "Aligned.out.bam",replacement = "",x = names(header)))

#### 1.3 Define metadata
##### Note that 'condition' and 'samples' are needed here. Conditions define the design of DESeq object, defines what can be tested for and how the data is normalized.

## fix colData
line <- c(rep(c("pt_iPS_Sandra","hs_iPS_6","h9"),8),c("hs_hIPS_48","hs_hIPS_48","pt_ciPS_PT5","pt_ciPS_PT5"))
species <- c(rep(c("chimp","human","human"),8),c("human","human","chimp","chimp"))
time <- c(rep(c(rep("d13",3),rep("d14",3),rep("d15",3),rep("d16",3)),2),rep("d14",4))
cond <- paste(line,time,sep = "-")
batch <- c(rep(c("batch_1","batch_2"),each=12),rep("batch_3",4))
colDat <- data.frame(species=species,time=time,condition=species,line=line,samples=samples[-c(1:6)],batch=batch)


##install_github(repo = "perllb/deseqAbstraction",username = "perllb")
#library(deseqAbstraction)

dabs <- deseqAbs$new("humanChimp",colData=colDat,filename=path)
dabs$makeDiffex()
res <- dabs$test$Default
genes.up <- getSignName(x = res,p = .01,l = .2)$up

GO_topGO_geneSet(dabs = dabs,geneSet = genes.up,outdir="~/Dropbox (MN)/Per/PhD/Projects/Chimp/Proteome_RNAseq.comparison/plots/")

