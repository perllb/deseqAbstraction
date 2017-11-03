#' @name maPlot
#' @description Plots scatter plot with log2FC on log2(baseMean), and marks significant genes
#' @param test: output object from results() of DESeq2
#' @param c1: string describing condition 1 in results deseq test
#' @param c2: string describing condition 2 in results deseq test
#' @title maPlot - my awesome function #2
#' @export maPlot
#' @example
#' test <- results(dds,contrast = c("condition","genex-KO","WT"))
#' maPlot(test = test, c1 = "KO", c2 = "WT" )

maPlot <- function(test,c1,c2) {

  colVec <- ifelse(test = test$padj<0.001, yes = ifelse(test = (is.na(test$padj)),yes = "black",no = "red"), no= "black")
  colVec[is.na(colVec)] <- "black"
  cexVec <- ifelse(test = test$padj<0.001, yes = ifelse(test = (is.na(test$padj)),yes = 0.15,no = 0.4), no= 0.15)
  cexVec[is.na(cexVec)] <- 0.15
  head(test)
  plot(log2(test[,1]),test[,2],
       col=colVec,
       cex=cexVec,
       pch=16,
       ylab=paste("log2(FC: [ ",c1," / ",c2," ])",sep=""),
       xlab="log2(mean expression)",
       main=paste(c1," / ",c2,sep=""))
  legend("topleft",legend = c("padj < 1e-3","padj >= 1e-3"),pch=16,col=c("red","black"),bty='n')

}

