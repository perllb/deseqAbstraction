#' @name PCAplotter
#' @description creates PCA plot
#' @param dat: normalized deseq data (with varianceStablizing, rlog etc)
#' @param ntop: how many genes to include
#' @param color: how to group colors
#' @param shape: how to group shapes
#' @param title: title of plot
#' @title PCAplotter my awesome function #1
#' @export PCAplotter
#' @example
#' vst <- varianceStabilizingTransformation(dds)
#' PCAplotter(dat = vst,ntop = 1000,color = colData$cellType, shape = colData$treatment,title = "PCA plot top1000 genes")



PCAplotter <- function(dat,ntop,color,shape,title) {

  data <- plotPCA(dat,returnData = T,ntop = ntop)
  #store the percentage variance for each PC
  percentVar <- round(100*attr(data,"percentVar"))
  #load ggplot2
  library("ggplot2")
  #plot
  ## (Cell and Vector are two datavectors in the colData data.frame of dds)
  ggplot(data,aes(PC1,PC2,color=color,shape=shape)) +
    geom_point(size=5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text(aes(label=colData$cell,hjust=rep(c(1,2,3),8),vjust=rep(c(3,1,2),8))) +
    ggtitle(title)

}
