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



deseqAbs <- R6Class("deseqAbs",
                    public = list(
                      name = "character",
                      filename = NULL,
                      rawfile = NULL,
                      rawCounts = NULL,
                      baseMean = NULL,
                      rpkmMean = NULL,
                      geneID = NULL,
                      colData = NULL,
                      sampleNames = NULL,
                      VST = NULL,
                      deseq = NULL,
                      rpkm = NULL,
                      test = NULL,
                      pos = NULL,
                      length = NULL,

                      initialize = function(name = NA,filename = NA,colData = NA) {

                        ### Check if all required parameters are set!
                        if(is.null(filename)) {
                          cat("-- ERROR: No filename (with full path) given! \n")
                          cat("--- Please provide name (and full path) of raw featureCounts output file upon creating this object\n\n")
                        }
                        if(is.null(colData)) {
                          cat("-- ERROR: No colData given! \n")
                        } else {
                          if(is.null(colData$condition)) {
                            cat("-- ERROR: colData has no 'condition' column. \n")
                          }
                          if(is.null(colData$samples)) {
                            cat("-- ERROR: colData has no 'samples' column. \n")
                          }
                        }

                        ## If all parameters are set, initialize
                        if(!is.null(filename) & !is.null(colData) & !is.null(colData$condition) & !is.null(colData$samples)) {

                          # message
                          self$greet()

                          self$name <- name
                          self$filename <- filename
                          self$colData <- colData

                          self$test <- list()

                          self$read_file(filename)
                          self$geneID <- as.character(self$rawfile[,1])
                          self$getPos()
                          self$getRawCounts()

                        }
                        else {
                          cat("Could not initialize object..")
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

                      readsAssigned = function(summaryFile=NULL) {

                        if(is.null(summaryFile)) {
                          # read summary file
                          sum <- read.delim(paste(self$filename,".summary",sep=""))
                        }else {
                          sum <- read.delim(summaryFile)
                        }

                        if(!is.null(self$colData$samples)) {
                          colnames(sum) <- c('a',as.character(self$colData$samples))
                        }

                        ## get total reads mapping to the genome
                        tot.map <- colSums(sum[,-1])
                        ## get number of reads assigned and not
                        assigned <- sum[1,-1]
                        notassigned <- tot.map-assigned

                        # plot
                        plot <- as.matrix(rbind(assigned = assigned,not.assigned = notassigned))

                        x <- barplot(plot,col=c("blue","grey80"),ylim=c(0,max(tot.map)*1.2),ylab="total read number",las=2)
                        legend("topleft",legend = c("not assigned","assigned"),fill=c("grey80","blue"))
                        title(main = "Reads assigned to annotation out of all mapped reads")

                      },

                      getPos = function() {
                        cat("- Fetching Positional info from file\n")
                        self$length <- self$rawfile[,6]
                        names(self$length) <- self$rawfile[,1]
                        self$pos <- self$rawfile[,2:5]
                        rownames(self$pos) <- self$rawfile[,1]
                        cat("- ..done. Get position of genes with $pos\n")
                      },

                      getAverage = function() {

                        self$getAverageReads()
                        self$getAverageRPKM()

                      },

                      getAverageReads = function() {

                        if(is.null(self$deseq)) {
                          cat("You must run deseq to get normalized read counts..\n")
                          self$makeDESeq()
                        }

                        if(length(levels(self$deseq$condition))>sum(duplicated(self$deseq$condition))) {
                          cat("Some of your levels do not have replicates.. \n")
                        } else {
                          ## normalized counts
                          cat("- Computing mean normalized counts of each condition\n")
                          baseMeanPerLvl <- sapply( levels(self$deseq$condition), function(lvl) rowMeans( counts(self$deseq,normalized=TRUE)[,self$deseq$condition == lvl] ) )
                          baseSDPerLvl <- sapply( levels(self$deseq$condition), function(lvl) apply( counts(self$deseq,normalized=TRUE)[,self$deseq$condition == lvl],1,sd ) )
                          colnames(baseSDPerLvl) <- paste("st.dev:",colnames(baseSDPerLvl),sep="")
                          self$baseMean <- list(Mean=baseMeanPerLvl,SD=baseSDPerLvl)
                          cat("- ..mean normalized counts computed for each condition. access mean with $baseMean$Mean, and st.dev with $baseMean$SD \n")
                        }
                      },

                      getAverageRPKM = function() {

                        ## RPKM
                        if(is.null(self$rpkm)) {self$makeRPKM() }

                        if(length(levels(self$deseq$condition))>sum(duplicated(self$deseq$condition))) {
                          cat("Some of your levels do not have replicates.. ")
                        } else {

                          if(!is.null(self$rpkm)) {
                          cat("- Computing mean RPKM of each condition\n")
                          baseMeanPerLvl <- sapply( levels(self$colData$condition), function(lvl) rowMeans( self$rpkm[,self$colData$condition == lvl] ) )
                          baseSDPerLvl <- sapply( levels(self$colData$condition), function(lvl) apply( self$rpkm[,self$colData$condition == lvl],1,sd ) )
                          colnames(baseSDPerLvl) <- paste("st.dev:",colnames(baseSDPerLvl),sep="")
                          self$rpkmMean <- list(Mean=baseMeanPerLvl,SD=baseSDPerLvl)
                          cat("- ..mean normalized expression computed for each condition. access mean with $baseMean$Mean, and st.dev with $baseMean$SD \n")
                          }
                        }
                      },

                      getRawCounts = function() {

                        cat("- Getting countData matrix\n")
                        self$rawCounts <- self$rawfile[,-c(1:6)]
                        rownames(self$rawCounts) <- self$geneID

                        if(!is.null(self$sampleNames)) {
                          colnames(self$rawCounts) <- self$sampleNames
                        } else if(!is.null(self$colData)) {
                          colnames(self$rawCounts) <- make.names(names = self$colData$condition,unique = T)
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

                      makeVST = function(blind=NULL) {


                        if(is.null(blind)) {

                          cat("== You need to define if you want to do blind dispersion estimates!\n")
                          cat("== set blind=FALSE if you expect a large fraction of genes to have large differences in counts explainable by experimental design!\n")
                          cat("== set blind=TRUE otherwise.\n")
                          cat("== If you are not sure, try both, and compare clustring results\n")

                        } else {

                          bl <- "blind"
                          if(!blind) { bl <- "not-blind" }

                          cat(paste(" - performing ",bl," variance stabilizing transformation \n",sep = ""))
                          self$VST <- varianceStabilizingTransformation(self$deseq,blind = blind)

                          if(!is.null(self$sampleNames)) {
                            colnames(assay(self$VST)) <- make.names(self$sampleNames,unique = T)
                          }else if(!is.null(self$colData)) {
                            colnames(assay(self$VST)) <- make.names(names = as.character(self$colData$condition),unique = T)
                          }

                          cat(paste(" - ..completed ",bl," variance stabilizing transformation \n",sep = ""))

                        }



                      },

                      # name: name of test
                      # c1: condition 1, has to be of the form of output of resultsNames(dds)
                      # c2: condition 2
                      # n1: condition 1, as index (integer) of the condition in vector from resultsNames(x$deseq)
                      # n1: condition 2, as index (integer) of the condition in vector from resultsNames(x$deseq)
                      makeDiffex = function(name=NULL,c1=NULL,c2=NULL,n1=NULL,n2=NULL) {

                        # if no specific conditions input, do default conditions
                        if(is.null(n1) & is.null(c1)) {

                          cat("- Testing for differential expression..\n")
                          self$test[["Default"]] <- results(self$deseq)
                          cat("- ..Diffex completed with default values. Access with $test.\n")

                        } else if(!is.null(n1) & !is.null(n2)) {

                          c1 <- gsub(pattern = "condition",replacement = "",resultsNames(self$deseq)[n1])
                          c2 <- gsub(pattern = "condition",replacement = "",resultsNames(self$deseq)[n2])

                          cat("- Testing for differential expression..\n")
                          cat(paste("-- Testing ",c1," vs. ",c2,"..\n",sep = ""))

                          if(!is.null(name)) {

                            self$test[[name]] <- results(self$deseq,contrasts = c("condition",c1,c2))

                          } else {
                            self$test[[paste("Test:",c1,"_vs._",c2,"",sep = "")]] <- results(self$deseq,contrasts = c("condition",c1,c2))
                          }

                          cat(paste("-- Test completed for ",c1," vs. ",c2,"..\n",sep = ""))

                        } else if(!is.null(c1) & !is.null(c2)) {

                          c1 <- gsub(pattern = "condition",replacement = "",c1)
                          c2 <- gsub(pattern = "condition",replacement = "",c2)

                          cat("- Testing for differential expression..\n")
                          cat(paste("-- Testing ",c1," vs. ",c2,"..\n",sep = ""))

                          if(!is.null(name)) {

                            self$test[[name]] <- results(self$deseq,contrasts = c("condition",c1,c2))

                          } else {
                            self$test[[paste("Test:",c1,"_vs._",c2,"",sep = "")]] <- results(self$deseq,contrasts = c("condition",c1,c2))
                          }
                          cat(paste("-- Test completed for ",c1," vs. ",c2,"..\n",sep = ""))
                        }

                      },

                      makeRPKM = function() {

                        cat("- Computing RPKM..\n")
                        self$rpkm <- self$rawCounts/ (self$length/1000) / (colSums(self$rawCounts/1000000))
                        rownames(self$rpkm) <- self$geneID
                        if(!is.null(self$sampleNames)) {
                          colnames(self$rpkm) <- make.names(names = self$sampleNames,unique = T)
                        }else if(!is.null(self$colData)) {
                          colnames(self$rpkm) <- make.names(names = self$colData$condition,unique = T)
                        }
                        cat("- ..RPKM computed. Access with $rpkm.\n")

                      },

                      pca = function(ntop=1000,title=NULL,label=NULL) {

                        if(is.null(title)) {title = "PCA"}
                        if(is.null(label)) {label = rep("",length(self$colData$condition))}

                        PCAplotter(dat = self$VST,
                                   color = self$colData$condition,
                                   title = title,
                                   ntop = ntop,
                                   label = label)

                      },

                      fullAuto = function() {

                        if(!is.null(self$colData) & !is.null(self$rawCounts)) {

                          if(!is.null(self$sampleNames)) {
                            colnames(self$rawCounts) <- make.names(self$sampleNames,unique = T)
                          }else if(!is.null(self$colData)) {
                            colnames(self$rawCounts) <- make.names(names = self$colData$condition,unique = T)
                          }
                          self$makeDESeq()
                          self$makeDiffex()
                          self$makeVST(blind = F)
                          self$makeRPKM()
                          self$getAverage()

                        }else {

                          cat("= Cannot do this yet.. \n")

                          }
                        if(is.null(self$colData)){

                          cat("= make colData data.frame!\n")

                        }
                        if(is.null(self$rawCounts)) {
                          cat("= add rawCounts matrix!\n")
                        }
                      },

                      ## get most variable genes
                      getVariable = function(ntop=50,sdcut = 0) {

                        sd <- apply(assay(self$VST),1,sd)

                        # if user specifices cutoff of sd for genes to return
                        if(sdcut > 0) {

                          return(rownames(self$VST[sd>sdcut,]))

                        } else { ## if no sdcut defined, return top 50 (or user specified ntop genes)

                          tmp <- assay(self$VST)
                          sd.ordered <- tmp[order(-apply(tmp,1,sd)),]
                          return(rownames(sd.ordered[1:ntop,]))

                        }
                      }
                    )
)

