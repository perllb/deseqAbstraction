#' @name genesClose
#' @description Plots scatter plot with log2FC on log2(baseMean), and marks significant genes
#' @param pos: position of all genes of interest (e.g. deseqAbsObj$pos). Must have $Chr, $Start, $End columns and gene names as rownames. If the annotation is collapsed (such as after featureCounts), only on tss (the most extreme) is used for each gene.
#' @param feat: a .bed file with position that genes should be close to (bed: chr - start - end - ID - . - strand)
#' @param distance: how far from features should a gene be to be included?
#' @title genesClose : Get genes close to features
#' @export genesClose
#' @examples
#' closeGenes <- enesClose(pos = desqAobj$pos, feat=L1.up.bed, distance = 10000 )

genesClose <- function(genPos,featPos,dist) {

  # for each gene, get chromsome, remove redundant terms
  genes.chr <- sapply(sapply(X = genPos$Chr,FUN = function(x) strsplit(as.character(x),';')),FUN = function(l) l[[1]])

  # for each gene, get strand, remove redundant terms
  genes.strand <- sapply(sapply(X = genPos$Strand,FUN = function(x) strsplit(as.character(x),';')),FUN = function(l) l[[1]])

  # for each gene, get TSS, remove redundant terms
  # depending on strand, get start of first or last
  genes.tss <- ifelse(test = (genes.strand == '+'),
                      # if + strand, get first exon start
                      yes = sapply(sapply(X = genPos$Start,FUN = function(x) strsplit(as.character(x),';')),FUN = function(l) l[[1]]),
                      # if - strand, get last exon end
                      no = sapply(sapply(X = genPos$End,FUN = function(x) strsplit(as.character(x),';')),FUN = function(l) l[[length(x = l)]]))

  genes <- data.frame(gene=rownames(genPos),chr=genes.chr,strand=genes.strand,tss=genes.tss)

  ## data frame to build
  allclose <- data.frame()
  # iterate through each feature -> get genes with TSS <50kb away)
  cat(">> Getting genes with TSS within",dist,"bps from each given feature. This might take a couple of minutes. Patience..")

  for (i in 1:nrow(featPos)){

    # get position of current L1
    feat.chr <- as.character(featPos[i,1])
    feat.tss <- ifelse(as.character(featPos[i,6])=="+",as.numeric(as.character(featPos[i,2])),as.numeric(as.character(featPos[i,3])))
    feat.end <- ifelse(as.character(featPos[i,6])=="+",as.numeric(as.character(featPos[i,3])),as.numeric(as.character(featPos[i,2])))
    feat.str <- as.character(featPos[i,6])
    feat.name <- as.character(featPos[i,4])
    # get genes on same chromosome
    genes.chr <- genes[genes$chr == feat.chr,]
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
  cat("- ..complete! Genes with TSS within",dist,"bps from each features TSS is computed.")
  return(allclose)
}
