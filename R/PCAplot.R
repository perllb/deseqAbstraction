#' @name PCAplotter
#' @description creates PCA plot
#' @param dat: normalized deseq data (with varianceStablizing, rlog etc)
#' @param ntop: how many genes to include
#' @param color: how to group colors
#' @param shape: how to group shapes. (by default NULL if only color by groups wanted)
#' @param title: title of plot
#' @param label: vector of labels for each point (optional)
#' @title PCAplotter my awesome function #1
#' @export PCAplotter
#' @examples
#' vst <- varianceStabilizingTransformation(dds);
#' PCAplotter(dat = vst,ntop = 1000,color = colData$cellType, shape = colData$treatment,title = "PCA plot top1000 genes", label = colData$sample)



PCAplotter <- function(dat,ntop=5000,color=NULL,shape=NULL,title="",label=NULL) {

  data <- plotPCA(dat,returnData = T,ntop = ntop)
  #store the percentage variance for each PC
  percentVar <- round(100*attr(data,"percentVar"))
  #load ggplot2
  library("ggplot2")
  #plot
  ## (Cell and Vector are two datavectors in the colData data.frame of dds)

  if(!is.null(shape)) { p <- ggplot(data,aes(PC1,PC2,color=color,shape=shape)) }
  if(is.null(shape)) { p <- ggplot(data,aes(PC1,PC2,color=color)) }

  if(!is.null(label)) {

    p + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      geom_text(aes(label=label),col="black",hjust=2) +
      scale_x_continuous(expand = c(.2,.2)) +
      scale_y_continuous(expand = c(.2,.2)) +
      geom_point(size=3) +
      ggtitle(title)

  }  else {

    p + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      ggtitle(title) +
      geom_point(size=5)
  }
}
