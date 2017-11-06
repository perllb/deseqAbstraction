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

maPlot <- function(test,c1,c2,p=.5,l=0) {

  #color up and down sign..
  colVec <- ifelse(test = resSai$padj<p,
                   yes=ifelse(test = resSai$log2FoldChange>l,
                              yes="firebrick3",
                              no=ifelse(resSai$log2FoldChange< -l,
                                        yes = "steelblue4",
                                        no = "black")),
                   no = "black")
  cexVec <- ifelse(test = resSai$padj<p,yes = 0.4,no = 0.15)
  colVec[is.na(colVec)] <- "black" ## if NA make sure it's not counted as <p
  #size of points
  cexVec <- ifelse(test = test$padj<p, yes = ifelse(test = (is.na(test$padj)),yes = 0.15,no = 0.4), no= 0.15)

  sign <- getSign(x = test,p = p,l = l)
  u <- nrow(sign$up)
  d <- nrow(sign$down)

  plot(log2(test[,1]),test[,2],
       col=colVec,
       cex=cexVec,
       pch=16,
       ylab=paste("log2(FC: [ ",c1," / ",c2," ])",sep=""),
       xlab="log2(mean expression)",
       main=paste(c1," / ",c2,sep=""),
       sup=paste("p-adj < ",p,", log2(fc) < ",l,sep("")))

  legend("topleft",legend = c(paste("up (",u,")",sep=""),paste("down (",d,")",sep = ""),"not significant"),pch=16,col=c("firebrick3","steelblue4","black"),bty='n')

}

