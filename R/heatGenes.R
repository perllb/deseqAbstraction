#' @name heatGenes
#' @description Plots heatmap showing the expression of the genes of given gene-list and having sd above given value
#' @param data: deseq matrix data (with varianceStablizing, rlog etc - eg. assay(vsd)
#' @param genes: a vector of gene IDs that can be grepped from rownames of data
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
#' mostVariableHeat(data = assay(vst),genes = genes, sd = 1, a1 = colData$cellLine, a2 = colData$treatment,n1 = "Cell Line",n2 = "Treatment",z = T)


heatGenes <- function(data,genes,a1=NULL,a2=NULL,n1=NULL,n2=NULL,sd,z=FALSE,cluster_col=T) {

  match <- paste(genes,collapse = "$|^")
  genes.exp <- data[grep(match,rownames(data)),]

  ## get standard deviation of each gene
  sd.exp <- apply(genes.exp,1,sd)

  ## Show rownames? only if less that 41 genes are plotted
  rowShow <- T; if(length(sd.exp[sd.exp>sd])>40) { rowShow <- F }
  print(paste("There are ",length(sd.exp[sd.exp>sd])," genes in you gene-set with sd > ",sd,".",sep=""))

  ## scale by row?
  scale <- 'none' ; if(z) { scale <- 'row' }


  if (!is.null(a1) & is.null(a2)) {
    df <- data.frame(Var1 = factor(a1))
    rownames(df) <- colnames(data)
    colnames(df) <- n1
    pheatmap(genes.exp[sd.exp>sd, ], annotation_col = df, cluster_rows = T, show_rownames = rowShow, cluster_cols = cluster_col,scale = scale)
  }
  else if (!is.null(a1) & !is.null(a2)) {
    df <- data.frame(Var1 = factor(a1), Var2 = factor(a2))
    rownames(df) <- colnames(data)
    colnames(df) <- c(n1, n2)
    pheatmap(genes.exp[sd.exp>sd, ], annotation_col = df, cluster_rows = T, show_rownames = rowShow, cluster_cols = cluster_col,scale = scale)
  }
  else if (!is.null(a2) & is.null(a1)) {
    df <- data.frame(Var1 = factor(a2))
    rownames(df) <- colnames(data)
    colnames(df) <- n2
    pheatmap(genes.exp[sd.exp>sd, ], annotation_col = df, cluster_rows = T, show_rownames = rowShow, cluster_cols = cluster_col,scale = scale)
  }
  else {
    pheatmap(genes.exp[sd.exp>sd, ], cluster_rows = T, show_rownames = rowShow, cluster_cols = cluster_col,scale = scale)
  }
}
