require(R6)

#' Class a simple interface to deseq data
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom grDevices dev.off png rainbow x11
#' @importFrom graphics legend axis barplot image par pie abline hist layout lines mtext plot plot.new rect text title
#' @importFrom stats quantile as.dendrogram density dist hclust median order.dendrogram reorder sd
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
#' @field colData a data.frame with condition information
#' @field sampleNames a vector given by user to provide suitable sampleNames to replace filenames from featureCounts
#' @field VST the output from varianceStabilizingTransformation(dds)
#' @field deseq the output from DESeq(dds)
#' @field rpkm matrix with rpkm values for all genes
#' @field test output from diffex analysis results(dds)
#' @field pos position data for each gene
#' @field length of each gene
#' @export



deseqAbs <- R6Class("deseqAbs",
                    public = list(
                      name = "character",
                      filename = NULL,
                      rawfile = "data.frame",
                      rawCounts = "matrix",
                      colData = "data.frame",
                      sampleNames = "vector",
                      VST = "matrix",
                      deseq = "DESeq",
                      rpkm = "matrix",
                      test = "DESeqResults",
                      pos = "matrix",
                      length = "vector",
                      initialize = function(name = NA,filename = NA) {
                        self$name <- name
                        self$filename <- filename
                        self$greet()
                        if(!is.na(filename)){
                          self$read_file(filename)
                          self$getPos()
                          self$getRawCounts()
                        } else {
                          cat("-Please provide name of raw featureCounts output file upon creating this object\n\n")
                        }
                      },

                      read_file = function(filename) {
                        if(!is.na(filename)) {
                          cat("- reading featureCount file\n")
                          self$rawfile <- read.csv(filename,header=T,sep = "\t",skip=1)
                        } else {
                          cat("- You must add name of raw featurecount file.\n")
                          cat("> dnmt$filename <- \"<path>/<featureCountOutput>\"\n")
                        }
                      },

                      getPos = function() {
                        cat("- Fetching Positional info from file\n")
                        self$length <- self$rawfile$Length
                        names(self$length) <- self$rawfile$Geneid
                        self$pos <- self$rawfile[,2:5]
                        rownames(self$pos) <- self$rawfile$Geneid
                      },

                      getRawCounts = function() {

                        cat("- Getting countData matrix\n")
                        self$rawCounts <- self$rawfile[,-c(1:6)]
                        rownames(self$rawCounts) <- self$rawfile$Geneid

                      },

                      greet = function() {
                        cat("- DESeq object created..\n")
                      },

                      makeDESeq = function() {

                        dds <- DESeqDataSetFromMatrix(countData = self$rawCounts,
                                                      colData = self$colData,
                                                      design =~ condition)
                        self$deseq <- DESeq(dds)
                      },

                      makeVST = function() {

                        self$VST <- varianceStabilizingTransformation(self$deseq)

                      },

                      makeDiffex = function() {

                        self$test <- results(self$deseq)

                      },

                      makeRPKM = function() {

                        cat("computing RPKM\n")

                      }
                    )
)
