#' @name mostVariableHeat
#' @description Plots heatmap showing the expression of most variable genes
#' @param data: normalized deseq matrix data (with varianceStablizing, rlog etc - eg. assay(vsd)
#' @param ntop: how many genes to display (the <ntop> most variable genes)
#' @param a1: annotation of the samples
#' @param a2: annotation of the samples
#' @param n1: name of annotation in a1
#' @param n2: name of annotation in a2
#' @title mostVariableHeat - heatmap of variable genes
#' @export mostVariableHeat
#' @examples
#' vst <- varianceStabilizingTransformation(dds)
#' mostVariableHeat(data = assay(vst),ntop = 100, a1 = colData$cellLine, a2 = colData$treatment,n1="Cell Line",n2="Treatment")

mostVariableHeat <- function(data,ntop=50,a1=NULL,a2=NULL,n1=NULL,n2=NULL) {

  library(graphics)
  library(pheatmap)
  library(RColorBrewer)

  if(!is.matrix(data) & !is.data.frame(data)) {
    data <- assay(data)
  }

  if(!is.matrix(data) & !is.data.frame(data)) {
    cat("ERROR: Data is not in correct format. Must be matrix or DESeq object")
  } else {

    sd.ordered <- data[order(-apply(data,1,sd)),]
    plotData <- sd.ordered[1:ntop,]

    heatGenes(data=plotData,genes = rownames(plotData),a1 = a1,a2 = a2,n1 = n1,n2 = n2)
   
    return(rownames(plotData))
  }
}
