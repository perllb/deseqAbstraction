#' @name getVariableGenes
#' @description returns the most variable genes in data
#' @param data: normalized deseq matrix data (with varianceStablizing, rlog etc - eg. assay(vsd)
#' @param ntop: how many genes to display (the <ntop> most variable genes)
#' @param sdcut: return all genes having higher sd than this cutoff
#' @title getVariableGenes - get genes with high variation
#' @export getVariableGenes
#' @example
#' vst <- varianceStabilizingTransformation(dds)
#' getVariableGenes(data = assay(vst),ntop = 100)

## get most variable genes
getVariableGenes = function(data,ntop=100,sdcut = 0) {

  if(!is.matrix(data) | !is.data.frame(data)) {
    data <- assay(data)
  }

  if(!is.matrix(data) | !is.data.frame(data)) {
    cat("ERROR: Data is not in correct format. Must be matrix or DESeq object")
  } else {
    sd <- apply(data,1,sd)

    # if user specifices cutoff of sd for genes to return
    if(sdcut > 0) {

      return(rownames(self$VST[sd>sdcut,]))

    } else { ## if no sdcut defined, return top 50 (or user specified ntop genes)

      tmp <- assay(self$VST)
      sd.ordered <- tmp[order(-apply(tmp,1,sd)),]
      return(rownames(sd.ordered[1:ntop,]))

    }
  }
}
