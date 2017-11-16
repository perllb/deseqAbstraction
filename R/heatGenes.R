#' @name heatGenes
#' @description Plots heatmap showing the expression of the genes of given gene-list and having sd above given value
#' @param data: deseq matrix data (with varianceStablizing, rlog etc - eg. assay(vsd)
#' @param genes: a vector of gene IDs to be extracted from rownames of data
#' @param a1: annotation of the samples
#' @param a2: annotation of the samples
#' @param n1: name of annotation in a1
#' @param n2: name of annotation in a2
#' @param sd: standard deviation cutoff for heatmap
#' @param z: logical operator. if TRUE, normalize by row (z-score). if FALSE (default), then plot normalized expression values
#' @param cluster_col: logical operator: if TRUE (default), cluster columns, if FALSE do not cluster columns
#' @title HeatGenes: Heatmap of your genes of interest!
#' @export heatGenes
#' @example
#' vst <- varianceStabilizingTransformation(dds)
#' genes <- c("TRIM28","DNMT1","ZNF52")
#' heatGenes(data = assay(vst),genes = genes, sd = 1, a1 = colData$cellLine, a2 = colData$treatment,n1 = "Cell Line",n2 = "Treatment",z = T)


heatGenes <- function(data,genes,a1=NULL,a2=NULL,n1=NULL,n2=NULL,sd=1,z=FALSE,cluster_col=T) {

  library(pheatmap)
  library(graphics)
  library(RColorBrewer)

  if(!is.matrix(data) & !is.data.frame(data)) {
    data <- assay(data)
  }

  if(!is.matrix(data) & !is.data.frame(data)) {
    cat("ERROR: Data is not in correct format. Must be matrix or DESeq object")
  } else {

    match <- paste(genes,collapse = "$|^")
    genes.exp <- data[grep(match,rownames(data)),]

    ## get standard deviation of each gene
    sd.exp <- apply(genes.exp,1,sd)

    # data to plot
    plotData <- genes.exp[sd.exp > sd,]

    ## Show rownames? only if less that 41 genes are plotted
    rowShow <- T; if(nrow(plotData)>80) { rowShow <- F }
    print(paste("There are ",nrow(plotData)," genes in you gene-set with sd > ",sd,".",sep=""))

    ## scale by row?
    scale <- 'none' ; if(z) { scale <- 'row' }

    if (!is.null(a1) & is.null(a2)) {

      df <- data.frame(Var1 = factor(a1))
      rownames(df) <- colnames(data)
      colnames(df) <- n1

      cols <- colorRampPalette(brewer.pal(8, "Set2"))
      mycolors <- cols(length(unique(a1)))
      names(mycolors) <- unique(a1)
      mycolors <- list(a = mycolors)
      names(mycolors) <- n1

      pheatmap(plotData, annotation_col = df, annotation_colors = mycolors,fontsize_row = 4, cluster_rows = T, show_rownames = rowShow, cluster_cols = cluster_col,scale = scale)

    } else if (!is.null(a1) & !is.null(a2)) {

      df <- data.frame(Var1 = factor(a1),Var2 = factor(a2))
      rownames(df) <- colnames(data)
      colnames(df) <- c(n1,n2)

      cols <- colorRampPalette(brewer.pal(9, "Set1"))
      mycolors <- cols(length(unique(a1)))
      names(mycolors) <- unique(a1)

      cols <- colorRampPalette(brewer.pal(7, "Set3"))
      mycolors2 <- cols(length(unique(a2)))
      names(mycolors2) <- unique(a2)

      mycolors <- list(a = mycolors,b = mycolors2)
      names(mycolors) <- c(n1,n2)

      pheatmap(plotData, annotation_col = df, annotation_colors = mycolors,fontsize_row = 4, cluster_rows = T, show_rownames = rowShow, cluster_cols = cluster_col,scale = scale)


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

      pheatmap(plotData, annotation_col = df, annotation_colors = mycolors,fontsize_row = 4, cluster_rows = T, show_rownames = rowShow, cluster_cols = cluster_col,scale = scale)

    } else {
      pheatmap(plotData, cluster_rows = T, fontsize_row = 4,  show_rownames = rowShow, cluster_cols = cluster_col,scale = scale)
    }

    return(rownames(plotData))
  }
}
