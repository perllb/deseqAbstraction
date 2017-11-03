#' @name sampleToSample
#' @description Plots heatmap showing sample-to-sample distances
#' @param data: normalized deseq data (with varianceStablizing, rlog etc)
#' @title sampleToSample my awesome function #2
#' @export sampleToSample
#' @example
#' vst <- varianceStabilizingTransformation(dds)
#' sampleToSample(data = vst)

sampleToSample <- function(data) {

  #sample-to-sample distance
  sampleDist <- dist(t(assay(data)))
  library("RColorBrewer")
  sampleDistMatrix <- as.matrix(sampleDist)
  rownames(sampleDistMatrix) <- data$condition
  colnames(sampleDistMatrix) <- NULL

  colors <- colorRampPalette(colors = c("red","white","Blue"))(25)
  par(mar=c(17,5,5,7))
  pheatmap(sampleDistMatrix,
           clustring_distance_rows=sampleDist,
           clustring_distance_cols=sampleDist,
           col=colors,
           main = "Sample to Sample distances")


}
