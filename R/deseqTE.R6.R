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

                     filteredRawfile = NULL,
                     TE.features = NULL,
                     genome = NULL,

                     initialize = function(name = NA,filename = NA,genome=NULL) {

                       if(is.null(genome)) {
                         cat("== you did not enter genome. please insert genome = hg38 or mm10!")
                       } else {

                         super$initialize(name = name,filename = filename)

                         cat("- Reading genomic RepeatMasker feature for ",genome)
                         self$genome = genome
                         self$TE.features <- self$getFeatures(genome)
                         cat("- genomic RepeatMasker feature for ",genome,"reading completed.. stored in $TE.features")
                       }
                      },

                     percentTE = function() {

                       te.summary <- paste(self$filename,".summary",sep="")
                       sum <- read.delim(te.summary)
                       tot.map <- colSums(sum[,-1])
                       map.rm <- sum[1,-1]

                       map.LINE <- colSums(getClass(self$rawCounts,"LINE"))
                       map.SINE <- colSums(getClass(self$rawCounts,"SINE"))
                       map.LTR <- colSums(getClass(self$rawCounts,"LTR"))
                       map.SVA <- colSums(getClass(self$rawCounts,"Retroposon"))


                       df <- t(data.frame(LINE=map.LINE,SINE=map.SINE,LTR=map.LTR,SVA=map.SVA))
                       colnames(df) <- make.names(names = te$colData$condition,unique = T)
                       plotPerc <- df*100/tot.map
                       col <- c("darkolivegreen3","indianred4","steelblue","tan4")
                       x <- barplot(plotPerc,ylim = c(0,max(colSums(plotPerc))*1.5),col=col,ylab="% reads mapping TE / mapping to genome")
                       legend("topleft",legend = rev(c("LINE","SINE","LTR","SVA")),fill=rev(col))
                       title("Percentage of reads mapping to EREs / genome")

                     },

                     getClass <- function(data,class) {

                       class.features <- as.character(self$TE.features[grep(class,self$TE.features$V2),1])
                       return(data[grep(paste(class.features,collapse = "|^"),rownames(data)),])

                     },

                     getFamily <- function(data,family) {

                       family.features <- as.character(self$TE.features[grep(class,self$TE.features$V3),1])
                       return(data[grep(paste(family.features,collapse = "|^"),rownames(data)),])

                     },

                     getSubFamily <- function(data,subfamily) {

                       sfamily.features <- as.character(self$TE.features[grep(class,self$TE.features$V3),1])
                       return(data[grep(paste(sfamily.features,collapse = "|^"),rownames(data)),])

                     },

                     getFeatures = function(genome) {

                       library(RCurl)
                       return(read.delim(text = getURL(paste("https://raw.githubusercontent.com/perllb/deseqAbstraction/master/data/",genome,".repeats.features.txt",sep = "")),header=F))

                     },

                     read_file = function(filename,filter=3) {
                       if(!is.na(filename)) {
                         cat("- reading featureCount file\n")
                         tmp <- read.csv(filename,header=T,sep = "\t",skip=1)
                         len <- nrow(tmp)
                         cat("- ..featureCount file reading done. \n")
                         cat("- Filter low read elements.. ")
                         self$rawfile <- tmp[rowMeans(tmp[,7:ncol(tmp)])>filter,]
                         lenf <- nrow(self$rawfile)
                         rem <- len-lenf
                         cat("- ..featureCount file filtering done.\n --Original rawfile had",len,"elements \n --After filtering",lenf,"elements remain.\n --",rem,"elements removed due to < ",filter," reads on average.\n --- If you want another cutoff for filtering, enter [filter = x] in call to method. \n --- Access filtered file with $rawfile\n")

                       } else {
                         cat("- You must add name of raw featurecount file.\n")
                         cat("> dnmt$filename <- \"<path>/<featureCountOutput>\"\n")
                       }
                     }
                   )
)
