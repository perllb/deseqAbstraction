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

    #show rownames if 40 or less genes are plotted
    rowShow <- T
    if(ntop>100) { rowShow <- F }

    if (!is.null(a1) & is.null(a2)) {

      df <- data.frame(Var1 = factor(a1))
      rownames(df) <- colnames(data)
      colnames(df) <- n1

      cols <- colorRampPalette(brewer.pal(8, "Set1"))
      mycolors <- cols(length(unique(a1)))
      names(mycolors) <- unique(a1)
      mycolors <- list(a = mycolors)
      names(mycolors) <- n1

      pheatmap(plotData, annotation_col = df, annotation_colors = mycolors,fontsize_row = 4,border_color = NA,
               cluster_rows = T, show_rownames = rowShow, cluster_cols = T,color = rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)))

    } else if (!is.null(a1) & !is.null(a2)) {

      df <- data.frame(Var1 = factor(a1),Var2 = factor(a2))
      rownames(df) <- colnames(data)
      colnames(df) <- c(n1,n2)

      cols <- colorRampPalette(brewer.pal(9, "Set1"))
      mycolors <- cols(length(unique(a1)))
      names(mycolors) <- unique(a1)

      cols <- colorRampPalette(brewer.pal(9, "Dark2"))
      mycolors2 <- cols(length(unique(a2)))
      names(mycolors2) <- unique(a2)

      mycolors <- list(a = mycolors,b = mycolors2)
      names(mycolors) <- c(n1,n2)

      pheatmap(plotData, annotation_col = df, annotation_colors = mycolors,fontsize_row = 4, cluster_rows = T,border_color = NA,
               show_rownames = rowShow, cluster_cols = T,color = rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)))


    } else if (!is.null(a2) & is.null(a1)) {

      a1<-a2
      df <- data.frame(Var1 = factor(a1))
      rownames(df) <- colnames(data)
      colnames(df) <- n1

      cols <- colorRampPalette(brewer.pal(9, "Set1"))
      mycolors <- cols(length(unique(a1)))
      names(mycolors) <- unique(a1)
      mycolors <- list(a = mycolors)
      names(mycolors) <- n1

      pheatmap(plotData, annotation_col = df, annotation_coloddsrs = mycolors,fontsize_row = 4, cluster_rows = T,border_color = NA,
               show_rownames = rowShow, cluster_cols = T,color = rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)))

    } else {
      pheatmap(plotData, cluster_rows = T, fontsize_row = 4,  show_rownames = rowShow,border_color = NA, cluster_cols = T,color = rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)))
    }
    return(rownames(plotData))
  }
}
