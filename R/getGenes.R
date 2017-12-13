#' @name getGenes
#' @description gets data for given genes. data can be of any format, as long as the genes IDs are rownames
#' @param data: matrix data (can be any format, but MUST have geneID to match 'genes' as rownames)
#' @param genes: a vector of gene IDs that can be matched and merged with rownames of data
#' @title getGenes: Get data for your genes of interest!
#' @export getGenes
#' @examples
#' getGenes(data = assay(vst),genes = genes)


getGenes <- function(data,genes) {
  return(merge(data.frame(genes=as.character(genes)),data.frame(data),by.x=1,by.y=0))
}
