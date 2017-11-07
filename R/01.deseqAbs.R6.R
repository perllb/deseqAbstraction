require(R6)
require(DESeq2)
#' Class a simple interface to deseq data
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



deseqAbs <- R6Class("deseqAbs",
                    public = list(
                      name = "character",
                      filename = NULL,
                      rawfile = NULL,
                      rawCounts = NULL,
                      geneID = NULL,
                      colData = NULL,
                      sampleNames = NULL,
                      VST = NULL,
                      deseq = NULL,
                      rpkm = NULL,
                      test = NULL,
                      pos = NULL,
                      length = NULL,
                      initialize = function(name = NA,filename = NA) {
                        self$name <- name
                        self$filename <- filename
                        self$greet()
                        if(!is.na(filename)){
                          self$read_file(filename)
                          self$geneID <- as.character(self$rawfile[,1])
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
                          cat("- ..featureCount file reading done. Access rawdata with $rawfile\n")
                        } else {
                          cat("- You must add name of raw featurecount file.\n")
                          cat("> dnmt$filename <- \"<path>/<featureCountOutput>\"\n")
                        }
                      },

                      getPos = function() {
                        cat("- Fetching Positional info from file\n")
                        self$length <- self$rawfile$Length
                        names(self$length) <- self$geneID
                        self$pos <- self$rawfile[,2:5]
                        rownames(self$pos) <- self$geneID
                        cat("- ..done. Get position of genes with $pos\n")
                      },

                      getRawCounts = function() {

                        cat("- Getting countData matrix\n")
                        self$rawCounts <- self$rawfile[,-c(1:6)]
                        rownames(self$rawCounts) <- self$geneID
                        if(!is.null(self$sampleNames)) {
                          colnames(self$rawCounts) <- self$sampleNames
                        }
                        cat("- ..done. access raw countData with $rawCounts\n")

                      },

                      greet = function() {
                        cat("- ..deseqAbs object created..\n")
                      },

                      makeDESeq = function() {

                        if(!is.data.frame(self$colData)) {
                          cat("= Add colData data.frame! E.g: >deseqAbs$colData <- data.frame(condition = x)\n")
                        } else {
                          cat("- Running DESeq")
                          dds <- DESeqDataSetFromMatrix(countData = self$rawCounts,
                                                        colData = self$colData,
                                                        design =~ condition)
                          self$deseq <- DESeq(dds)
                          cat("- ..DESeq done. Access object with $deseq \n")
                        }
                      },

                      makeVST = function() {

                        self$VST <- varianceStabilizingTransformation(self$deseq)

                      },

                      makeDiffex = function() {

                        cat("- Testing for differential expression..\n")
                        self$test <- results(self$deseq)
                        cat("- ..Diffex done. Access with $test.\n")

                      },

                      makeRPKM = function() {

                        cat("- Computing RPKM..\n")
                        self$rpkm <- self$rawCounts/ (self$length/1000) / (colSums(self$rawCounts/1000000))
                        rownames(self$rpkm) <- self$geneID
                        if(!is.null(self$sampleNames)) {
                          colnames(self$rpkm) <- self$sampleNames
                        }
                        cat("- ..RPKM computed. Access with $rpkm.\n")

                      },

                      fullAuto = function() {

                        if(!is.null(self$colData) & !is.null(self$rawCounts)) {

                          if(!is.null(self$sampleNames)) {
                            colnames(self$rawCounts) <- make.names(self$sampleNames,unique = T)
                          }
                          self$makeDESeq()
                          self$makeDiffex()
                          self$makeVST()
                          self$makeRPKM()

                        }else {

                          cat("= Cannot do this yet.. \n")

                          }
                        if(is.null(self$colData)){

                          cat("= make colData data.frame!\n")

                        }
                        if(is.null(self$rawCounts)) {
                          cat("= add rawCounts matrix!\n")
                        }
                    }
))
