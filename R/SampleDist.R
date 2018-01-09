#' @name sampleToSample
#' @description Plots heatmap showing sample-to-sample distances
#' @param data: normalized deseq data (with varianceStablizing, rlog etc)
#' @param samples: vector of sample names to show in map
#' @title sampleToSample my awesome function #2
#' @export sampleToSample
#' @examples
#' vst <- varianceStabilizingTransformation(dds)
#' sampleToSample(data = vst)

sampleToSample <- function(data,samples=NULL) {

  library(graphics)
  library(pheatmap)
  library(RColorBrewer)

  #sample-to-sample distance
  sampleDist <- dist(t(assay(data)))
  sampleDistMatrix <- as.matrix(sampleDist)
  if(is.null(samples)){
    rownames(sampleDistMatrix) <- colnames(assay(vst))
  }else{
    rownames(sampleDistMatrix) <- samples
  }
  colnames(sampleDistMatrix) <- NULL

  colors <- colorRampPalette(brewer.pal(11,"RdBu"))(100)
  par(mar=c(17,5,5,7))
  pheatmap(sampleDistMatrix,
           clustring_distance_rows=sampleDist,
           clustring_distance_cols=sampleDist,
           col=colors,
           main = "Sample to Sample distances")
}
