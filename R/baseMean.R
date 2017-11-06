#' @name getAverage
#' @description Calculates average expression of each group, based on condition of dds object. Returns list with Mean and SD values for each gene in each condition
#' @param dds: dds object. MUST have condition as colData grouping!
#' @title Get average read numbers within each condition
#' @export getAverage
#' @examples
#' dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, desing =~ condition)
#' baseMeans <- getAverage(dds)
#' means <- baseMeans$Mean
#' SD <- baseMeans$SD


getAverage <- function(dds) {

  baseMeanPerLvl <- sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )
  baseSDPerLvl <- sapply( levels(dds$condition), function(lvl) apply( counts(dds,normalized=TRUE)[,dds$condition == lvl],1,sd ) )
  colnames(baseSDPerLvl) <- paste("st.err:",colnames(baseSDPerLvl),sep="")
  return(list(Mean=baseMeanPerLvl,SD=baseSDPerLvl))

}
