#' @name mostVariableHeat
#' @description Plots heatmap showing the expression of most variable genes
#' @param data: normalized deseq data (with varianceStablizing, rlog etc)
#' @param ntop: how many genes to display (the <ntop> most variable genes)
#' @param a1: annotation of the samples
#' @param a2: annotation of the samples
#' @title mostVariableHeat - heatmap of variable genes
#' @export mostVariableHeat
#' @example
#' vst <- varianceStabilizingTransformation(dds)
#' mostVariableHeat(data = vst,ntop = 100, a1 = colData$cellLine, a2 = colData$treatment)

mostVariableHeat <- function(data,ntop=50,a1=NULL,a2=NULL) {

  ## Get most variable genes
  sd <- assay(data)[order(-apply(assay(data),1,sd)),]
  head(sd)

  if(!is.na(a1) & is.na(a2)) {

    df <- as.data.frame(a1=a1)
    pheatmap(sd[1:ntop,],annotation_col = df,cluster_rows = T,show_rownames = T,
             cluster_cols = T)
  } else if(!is.na(a1) & !is.na(a2)) {
    df <- as.data.frame(a1=a1,a2=a2)
    pheatmap(sd[1:ntop,],annotation_col = df,cluster_rows = T,show_rownames = T,
             cluster_cols = T)
  } else if(!is.na(a2) & is.na(a1)) {
    df <- as.data.frame(a2=a2)
    pheatmap(sd[1:ntop,],annotation_col = df,cluster_rows = T,show_rownames = T,
             cluster_cols = T)
  }else{

    pheatmap(sd[1:ntop,],cluster_rows = T,show_rownames = T,
             cluster_cols = T)
  }


}
