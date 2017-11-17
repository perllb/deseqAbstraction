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

                       map.LINE <- colSums(self$getTEClass(self$rawCounts,"LINE"))
                       map.SINE <- colSums(self$getTEClass(self$rawCounts,"SINE"))
                       map.LTR <- colSums(self$getTEClass(self$rawCounts,"LTR"))

                       if ( self$genome == "hg38") {
                         map.SVA <- colSums(self$getTEClass(self$rawCounts,"Retroposon"))
                         df <- t(data.frame(LINE=map.LINE,SINE=map.SINE,LTR=map.LTR,SVA=map.SVA))
                         colnames(df) <- make.names(names = self$colData$condition,unique = T)
                         plotPerc <- df*100/tot.map
                         col <- c("darkolivegreen3","indianred4","steelblue","tan4")
                         x <- barplot(plotPerc,ylim = c(0,max(colSums(plotPerc))*1.5),col=col,ylab="% reads mapping TE / mapping to genome",las=2)
                         legend("topleft",legend = rev(c("LINE","SINE","LTR","SVA")),fill=rev(col))
                         title("Percentage of reads mapping to EREs / genome")

                       } else {

                         df <- t(data.frame(LINE=map.LINE,SINE=map.SINE,LTR=map.LTR))
                         colnames(df) <- make.names(names = self$colData$condition,unique = T)
                         plotPerc <- df*100/tot.map
                         col <- c("darkolivegreen3","indianred4","steelblue","tan4")
                         x <- barplot(plotPerc,ylim = c(0,max(colSums(plotPerc))*1.5),col=col,ylab="% reads mapping TE / mapping to genome",las=2)
                         legend("topleft",legend = rev(c("LINE","SINE","LTR")),fill=rev(col))
                         title("Percentage of reads mapping to EREs / genome")

                       }
                     },

                     upSubFamily = function(p=.05,l=0,n=10) {

                       # get upregulated elements
                       up <- getSignName(x = self$test$Default,p = p,l = l)$up

                       # if human genome, include SVAs..
                       if(self$genome=="hg38") {

                         up.sva <- up[grep("SVA",up)]
                         sva.up.fam <- substr(x = up.sva,start = 0,stop=5)
                         up.noSVA <- up[grep("SVA",up,invert = T)]
                         names <- gsub("\\_.*","",up.noSVA)
                         names <- c(names,sva.up.fam)

                         RE.features <- self$TE.features[grep("LINE|SINE|LTR|Retroposon",self$TE.features$V2),]

                         tab <- table(names)
                         tab <- tab[order(-tab)][1:n]
                         tab.df <- data.frame(names=names(tab),count=tab)

                         # fix color to color family by class
                         tab.fam <- merge(tab.df,unique(RE.features[,c(1,2)]),by.x=1,by.y=1,all.x=T,sort=F)
                         col <- ifelse(tab.fam$V2=="LINE",yes="darkolivegreen",
                                       no=ifelse(tab.fam$V2=="LTR",yes="darkblue",
                                                 no = ifelse(tab.fam$V2=="Retroposon",yes="red",
                                                             no=ifelse(tab.fam$V2=="SINE",yes="orange",no="black"))))

                         par(mar=c(8,4,4,4))
                         x <- barplot(tab,ylim = c(0,max(tab)*1.4),ylab = "Number of upregulated elements",col = col,las = 2)
                         title(paste("Number of upregulated elements in each family\np-adj<",p,", log2(fc)<",l,sep = ""))

                       } else {

                         RE.features <- self$TE.features[grep("LINE|SINE|LTR",self$TE.features$V2),]

                         tab <- table(names)
                         tab <- tab[order(-tab)][1:n]
                         tab.df <- data.frame(names=names(tab),count=tab)

                         # fix color to color family by class
                         tab.fam <- merge(tab.df,unique(RE.features[,1:2]),by.x=1,by.y=1,all.x=T,sort=F)
                         col <- ifelse(tab.fam$V2=="LINE",yes="darkolivegreen",
                                       no=ifelse(tab.fam$V2=="LTR",yes="darkblue",
                                                 no=ifelse(tab.fam$V2=="SINE",yes="orange",no="black")))

                         par(mar=c(8,4,4,4))
                         x <- barplot(tab,ylim = c(0,max(tab)*1.4),ylab = "Number of upregulated elements",col = col,las = 2)
                         title(paste("Number of upregulated elements in each family\np-adj<",p,", log2(fc)<",l,sep = ""))

                       }
                     },

                     upFamilyBar = function(p=.05,l=0,n=10) {

                       # get upregulated elements
                       up <- getSignName(x = self$test$Default,p = p,l = l)$up

                       # if human genome, include SVAs..
                       if(self$genome=="hg38") {

                         up.sva <- up[grep("SVA",up)]
                         sva.up.fam <- substr(x = up.sva,start = 0,stop=5)
                         up.noSVA <- up[grep("SVA",up,invert = T)]
                         names <- gsub("\\_.*","",up.noSVA)
                         names <- c(names,sva.up.fam)

                         RE.features <- self$TE.features[grep("LINE|SINE|LTR|Retroposon",self$TE.features$V2),]
                         name.family <- as.character(merge(names,RE.features,by=1)$V3)
                         tab <- table(name.family)
                         tab <- tab[order(-tab)][1:n]
                         tab.df <- data.frame(names=names(tab),count=tab)

                         # fix color to color family by class
                         tab.fam <- merge(tab.df,unique(RE.features[,2:3]),by.x=1,by.y=2,all.x=T,sort=F)
                         col <- ifelse(tab.fam$V2=="LINE",yes="darkolivegreen",
                                       no=ifelse(tab.fam$V2=="LTR",yes="darkblue",
                                                 no = ifelse(tab.fam$V2=="Retroposon",yes="red",
                                                             no=ifelse(tab.fam$V2=="SINE",yes="orange",no="black"))))

                         par(mar=c(8,4,4,4))
                         x <- barplot(tab,ylim = c(0,max(tab)*1.4),ylab = "Number of upregulated elements",col = col,las = 2)
                         title(paste("Number of upregulated elements in each family\np-adj<",p,", log2(fc)<",l,sep = ""))

                       } else {

                         RE.features <- self$TE.features[grep("LINE|SINE|LTR",self$TE.features$V2),]
                         name.family <- as.character(merge(names,RE.features,by=1)$V3)
                         tab <- table(name.family)
                         tab <- tab[order(-tab)][1:n]
                         tab.df <- data.frame(names=names(tab),count=tab)

                         # fix color to color family by class
                         tab.fam <- merge(tab.df,unique(RE.features[,2:3]),by.x=1,by.y=2,all.x=T,sort=F)
                         col <- ifelse(tab.fam$V2=="LINE",yes="darkolivegreen",
                                       no=ifelse(tab.fam$V2=="LTR",yes="darkblue",
                                                 no=ifelse(tab.fam$V2=="SINE",yes="orange",no="black")))

                         par(mar=c(8,4,4,4))
                         x <- barplot(tab,ylim = c(0,max(tab)*1.4),ylab = "Number of upregulated elements",col = col,las = 2)
                         title(paste("Number of upregulated elements in each family\np-adj<",p,", log2(fc)<",l,sep = ""))

                       }
                     },

                     upClassBar = function(p=.05,l=0) {

                       # get upregulated elements
                       up <- getSignName(x = self$test$Default,p = p,l = l)$up

                       # if human genome, include SVAs..
                       if(self$genome=="hg38") {

                         up.sva <- up[grep("SVA",up)]
                         sva.up.fam <- substr(x = up.sva,start = 0,stop=5)
                         up.noSVA <- up[grep("SVA",up,invert = T)]
                         names <- gsub("\\_.*","",up.noSVA)
                         names <- c(names,sva.up.fam)

                         name.class <- as.character(merge(names,self$TE.features,by=1)$V2)
                         tab <- table(name.class)
                         tab.RE <- tab[grep("LINE$|SINE$|LTR$|Retroposon$",names(tab))]
                         names(tab.RE)[names(tab.RE)=="Retroposon"] <- "SVA"
                         col <- c("darkolivegreen3","indianred4","steelblue","tan4")
                         x <- barplot(tab.RE,ylim = c(0,max(tab)*1.4),ylab = "Number of upregulated elements",col = col,las = 2)
                         title(paste("Number of upregulated elements in each class\np-adj<",p,", log2(fc)<",l,sep = ""))

                       } else { ## mouse (no SVA)

                         names <- gsub("\\_.*","",up)

                         name.class <- as.character(merge(names,self$TE.features,by=1)$V2)
                         tab <- table(name.class)
                         tab.RE <- tab[grep("LINE$|SINE$|LTR$",names(tab))]
                         col <- c("darkolivegreen3","indianred4","steelblue")
                         x <- barplot(tab.RE,ylim = c(0,max(tab)*1.4),ylab = "Number of upregulated elements",col = col,las = 2)
                         title(paste("Number of upregulated elements in each class\np-adj<",p,", log2(fc)<",l,sep = ""))

                       }
                     },

                     maPlotTE = function(testData,c1="cond1",c2="cond2",p=.05,l=0) {

                       ## plot standard maPlot
                       maPlot(testData,c1,c2)

                       ## define colors for TEs
                       col <- c("darkolivegreen3","indianred4","steelblue","tan4")
                       idx <- 1 # index to increment pr TE class

                       for (ere in c("LINE","SINE","LTR","SVA")) {

                         # for SVA, grep on Retroposon
                         if(ere == "SVA") { ere <- "Retroposon" }

                         curr <- self$getTEClass(data = testData,ere)
                         curr.up <- getSign(x = curr,p = p,l=l)$up
                         curr.down <- getSign(x = curr,p = p,l=l)$down
                         curr.sig <- rbind(curr.up,curr.down)

                          if(nrow(curr.sig)>0){
                           points(log2(curr.sig$baseMean),curr.sig$log2FoldChange,col=col[idx],pch=16)
                         }

                         idx <- idx+1
                       }

                       legend("topright",legend = rev(c("LINE","SINE","LTR","SVA")),col=rev(col),pch=16)
                     },

                     getTEClass = function(data,TEclass) {

                       TEclass.features <- as.character(self$TE.features[grep(TEclass,self$TE.features$V2),1])
                       return(data[grep(paste(TEclass.features,collapse = "|^"),rownames(data)),])

                     },

                     getFamily = function(data,family) {

                       family.features <- as.character(self$TE.features[grep(class,self$TE.features$V3),1])
                       return(data[grep(paste(family.features,collapse = "|^"),rownames(data)),])

                     },

                     getSubFamily = function(data,subfamily) {

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
