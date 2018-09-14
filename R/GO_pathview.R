#' @name GO_pathview
#' @description Do GO term enrichement and analysis for a list of genes (pathway analysis and plots)
#' @param dds: deseq object
#' @param species: only human (hsa) annotation currently supported, mm10 coming soon
#' @param Npathways: How many pathways to plot
#' @title GOanalysis in R  - pathways
#' @export GO_pathview
#' @examples
#' dabs <- deseqAbs$new(name="drugTest",colData=colDat)
#' GO_pathview(dabs = dabs,species = "hsa",Npathways = 4)

GO_pathview <- function(dabs=NULL,species="hsa",Npathways=5) {
  
  ### Dependencies:
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("pathview")
  #biocLite("gage")
  #biocLite("gageData")
  #if(species=="hg38"){ biocLite("org.Hs.eg.db") }
  
  getwd()
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
  
  
  #### KEGG
  data("kegg.sets.hs")
  data("sigmet.idx.hs")  
  
  fc <- res$log2FoldChange
  names(fc) <- res$entrez
  
  # Get the results
  kegg.sets.hs.2 = kegg.sets.hs[sigmet.idx.hs]
  
  keggres2 = gage(fc, gsets=kegg.sets.hs.2, same.dir=TRUE)
  
  ### Greater
  # Get the pathways
  library(dplyr)
  
  ## plot greater
  keggrespathways = data.frame(id=rownames(keggres2$greater), keggres2$greater) %>% 
    tbl_df() %>% 
    filter(row_number()<=Npathways) %>% 
    .$id %>% 
    as.character()
  keggrespathways
  
  # Get the IDs.
  keggresids = substr(keggrespathways, start=1, stop=8)
  
  # Define plotting function for applying later
  
  plot_pathway = function(pid) pathview(gene.data=fc, pathway.id=pid, species=species, new.signature=FALSE)
  
  # plot multiple pathways (plots saved to disk and returns a throwaway list object)
  detach("package:dplyr",unload = T)
  tmp = sapply(keggresids, function(pid) pathview(gene.data=fc, pathway.id=pid, species=species))
  
  #### Less
  # Get the pathways
  
  library(dplyr)
  ## plot lesser
  keggrespathways = data.frame(id=rownames(keggres2$less), keggres2$less) %>% 
    tbl_df() %>% 
    filter(row_number()<=Npathways) %>% 
    .$id %>% 
    as.character()
  keggrespathways
  
  # Get the IDs.
  keggresids = substr(keggrespathways, start=1, stop=8)
  keggresids
  
  # Define plotting function for applying later
  
  plot_pathway = function(pid) pathview(gene.data=fc, pathway.id=pid, species=species, new.signature=FALSE)
  
  detach("package:dplyr",unload = T)
  # plot multiple pathways (plots saved to disk and returns a throwaway list object)
  tmp = sapply(keggresids, function(pid) pathview(gene.data=fc, pathway.id=pid, species=species))
  
}
