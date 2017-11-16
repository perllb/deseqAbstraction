require(R6)
require(DESeq2)
#' Class a simple interface to deseq with transposons - inherits from deseqAbstraction
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom grDevices dev.off png rainbow x11
#' @importFrom graphics legend axis barplot image par pie abline hist layout lines mtext plot plot.new rect text title
#' @importFrom stats quantile as.dendrogram density dist hclust median order.dendrogram reorder sd
#' @importFrom DESeq2 DESeq results DESeqDataSetFromMatrix varianceStabilizingTransformation plotPCA
#' @export
#' @keywords Deseq
#' @return Object of \code{\link{R6Class}} to store Deseq data.
#' @format \code{\link{R6Class}} object.
#' @examples
#' dnmt <- deseqAbs$new("DNMT1KO",fc.file)
#' dnmt$filename
#' head(dnmt$rawfile)
#' head(dnmt$pos)
#' head(dnmt$length)
#' head(dnmt$rawCounts)
#' ## sampleNames
#' dnmt$sampleNames <- colnames(dnmt$rawCounts)
#' ## define coldata conditions
#' dnmt$colData <- data.frame(condition=c(rep("CTR",3),rep("KO",3)))
#' dnmt$makeDESeq()
#' dnmt$deseq
#' dnmt$makeDiffex()
#' dnmt$test
#' dnmt$makeVST()
#' dnmt$makeRPKM()
#' @field name the name of the experiment
#' @field filename the name of the raw featurecount output
#' @field rawfile the raw featurecount file
#' @field rawCounts the count matrix (removing position and length info column 1-6) with ID rownames
#' @field baseMean list of data.frames: Mean and SD: mean and SD of normalized count data
#' @field rpkmMean list of data.frames: Mean and SD: mean and SD of RPKM
#' @field geneID geneIDs
#' @field colData a data.frame with condition information
#' @field sampleNames a vector given by user to provide suitable sampleNames to replace filenames from featureCounts
#' @field VST the output from varianceStabilizingTransformation(dds)
#' @field deseq the output from DESeq(dds)
#' @field rpkm matrix with rpkm values for all genes
#' @field test output from diffex analysis results(dds)
#' @field pos position data for each gene
#' @field length of each gene
#' @export
#'

deseqTE <- R6Class("deseqTE",
                   inherit = deseqAbs,
                   public = list(
                     familyDF = NULL,
                     filteredRawfile = NULL,

                     read_file = function(filename) {
                       if(!is.na(filename)) {
                         cat("- reading featureCount file\n")
                         tmp <- read.csv(filename,header=T,sep = "\t",skip=1)
                         cat("- ..featureCount file reading done. \n")
                         cat("- Filter low read elements.. ")
                         self$rawfile <- tmp[rowMeans(tmp[,7:ncol(tmp)])>3,]
                         cat("- ..featureCount file filtering done. Access rawdata with $rawfile\n")

                       } else {
                         cat("- You must add name of raw featurecount file.\n")
                         cat("> dnmt$filename <- \"<path>/<featureCountOutput>\"\n")
                       }
                     },


                     percGenome = function() {



                     }





                   )
)
