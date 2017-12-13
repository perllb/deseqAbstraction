#' @name genesClose
#' @description Get genes (or other features) close to some other features (.bed). Must have $Chr $Start $End $Strand $ID as column names. In case of collaped positions, such as from featureCounts output, the first TSS will be used, and alternative transcripts are ignored. For including all transcript variants, a bed file with one line for each transcript of each gene is optimal.
#' @param genPos: position of all genes of interest in Bed format. Must have $Chr, $Start, $End $ID $Strand columns. If the annotation is collapsed (such as after featureCounts), only on tss (the most extreme) is used for each gene.
#' @param featPos: a .bed file with position that genes should be close to (bed: chr - start - end - ID - . - strand). Must have $Chr, $Start, $End $ID $Strand column names.
#' @param dist: how far from features should a gene be to be included?
#' @title genesClose : Get genes close to features
#' @export genesClose
#' @examples
#' closeGenes <- enesClose(pos = gencode.genes, feat=L1.up, distance = 10000 )

genesClose <- function(genPos,featPos,dist=10000) {

  if('Chr' %in% colnames(genPos) & 'Strand' %in% colnames(genPos) & 'Start' %in% colnames(genPos) & 'End' %in% colnames(genPos) & 'ID' %in% colnames(genPos)){
    if('Chr' %in% colnames(featPos) & 'Strand' %in% colnames(featPos) & 'Start' %in% colnames(featPos) & 'End' %in% colnames(featPos) & 'ID' %in% colnames(featPos)){

      ## data frame to build
      allclose <- data.frame()
      # iterate through each feature -> get genes with TSS <50kb away)
      cat(">> Getting genes with TSS within",dist,"bps from each given feature. This might take a couple of minutes. Patience..\n")

      for (i in 1:nrow(featPos)){

        # get position of current L1
        feat.chr <- as.character(featPos$Chr[i])
        feat.tss <- ifelse(as.character(featPos$Strand[i])=="+",as.numeric(as.character(featPos$Start[i])),as.numeric(as.character(featPos$End[i])))
        feat.end <- ifelse(as.character(featPos$Strand[i])=="+",as.numeric(as.character(featPos$End[i])),as.numeric(as.character(featPos$Start[i])))
        feat.str <- as.character(featPos$Strand[i])
        feat.name <- as.character(featPos$ID[i])
        # get genes on same chromosome
        genes.chr <- genPos[genPos$Chr == feat.chr,]
        # get genes with TSS < 50kb from L1 start
        closeGenes <-genes.chr[abs(as.numeric(as.character(genes.chr$tss))-feat.tss)<dist,]
        ## add L1 data
        closeGenes[,5] <- rep(as.character(featPos$V4[i]),nrow(closeGenes))
        closeGenes[,6] <- rep(feat.str,nrow(closeGenes))
        closeGenes[,7] <- rep(feat.tss,nrow(closeGenes))
        allclose <- rbind(allclose,closeGenes)

      }

      colnames(allclose) <- c('geneClose','gene_chr','gene_strand','gene_tss','feature','feature_strand','feature_tss')
      allclose[,'distance'] <- as.numeric(as.character(allclose$gene_tss))-as.numeric(as.character(allclose$feature_tss))
      cat("- ..complete! Genes with TSS within",dist,"bps from each features TSS is computed.\n")
      return(allclose)
    }else{
      print("> ERROR: Colnames of your featPos table/data.frame must be with Chr, Start, End, Strand, ID colnames")
      print("> ERROR: Please change your column names.")
    }
  }else{
    print("> ERROR: Colnames of your genPos table/data.frame must be with Chr, Start, End, Strand, ID colnames")
    print("> ERROR: Please change your column names.")
  }
}
