#' @name getSign
#' @description Get the significantly changed genes, returning deseq-table with up and down-regulated genes separately in a list.
#' @param x: object of results(dds), or a data.frame with column named "padj" with adjusted p-values, and one named "log2FoldChange" with log2(fold changes)
#' @param p: p-adj cutoff to use
#' @param l: log2(fc) cutoff to use
#' @title getSign genes as table with test data
#' @export getSign
#' @example
#' test <- results(dds,contrast = c("condition","genex-KO","WT"))
#' sign <- getSign(x = test, p = 0.01, l = 0.5)
#' up <- sign$up
#' down <- sign$down

getSign <- function(x,p = .01,l = .2) {

  up <- x[!is.na(x$padj) & x$padj < p & x$log2FoldChange > l,]
  down <- x[!is.na(x$padj) & x$padj < p & x$log2FoldChange < -l,]
  return(list(up=up,down=down))

}

#' @name getSignName
#' @description Get the significantly changed genes, returning IDs of up and down-regulated genes separately in a list.
#' @param x: object of results(dds)
#' @param p: p-adj cutoff to use
#' @param l: log2(fc) cutoff to use
#' @title getSign IDs of genes that are sign. up or down
#' @export getSignName
#' @example
#' test <- results(dds,contrast = c("condition","genex-KO","WT"))
#' sign <- getSign(x = test, p = 0.01, l = 0.5)
#' up <- sign$up
#' down <- sign$down
#'
getSignName <- function(x,p,l) {

  up <- x[!is.na(x$padj) & x$padj < p & x$log2FoldChange > l,]
  down <- x[!is.na(x$padj) & x$padj < p & x$log2FoldChange < -l,]
  return(list(up=rownames(up),down=rownames(down)))

}
