#' @name meanBar
#' @description barplot of gene expression
#' @param deseqAbs: deseqAbs object with rpkmMean and test data
#' @param genes: a vector of gene IDs to be extracted from rownames of data
#' @param cond: a vector of conditions to plot
#' @param rpkm: Set to TRUE if RPKM should be plotted instead of baseMean
#' @title meanBar: Barplot of your genes of interest!
#' @export meanBar
#' @examples
#' genes <- c("DNMT1","TRIM28","PAX6","DCX","SOX2","AGO2")
#' meanBar(dnmt,genes)

meanBar <- function(deseqAbs,genes,cond=NULL,rpkm=FALSE) {
  
  cat(">>> meanBar: plot your genes:\n")
  cat(">>",genes)
  
  ## set graphical area
  row <- ifelse(test = sqrt(length(genes))%%1 > .5,yes = floor(sqrt(length(genes)))+1,no = floor(sqrt(length(genes)))) 
  par(mfrow=c(row,ceiling(sqrt(length(genes)))),mar=c(7,5,4,4))
  library(RColorBrewer)
  cols <- colorRampPalette(brewer.pal(8, "Greys"))
  
  # if no condition defined, run default
  if ( is.null(cond) ) {
    data <- ifelse(rpkm,deseqAbs$rpkmMean$Mean,deseqAbs$baseMean$Mean)
    mycolors <- cols(length(unique(deseqAbs$colData$condition)))
    padj.a <- deseqAbs$test$Default$padj
    names(padj.a) <- rownames(deseqAbs$test$Default)
  } else {
    data <- ifelse(rpkm,deseqAbs$rpkmMean$Mean[,cond],deseqAbs$baseMean$Mean[,cond])
    mycolors <- cols(length(cond))
    str <- paste()
    deseqAbs$makeDiffex(name='tmptest',c1=cond[1],c2=cond[2])
    padj.a <- deseqAbs$test$tmptest$padj
    names(padj.a) <- rownames(deseqAbs$test$tmptest)
  }
  
  # plot one for each gene 
  for(gene in genes) {
    
    plot <- data[gene,]
    sd <- data[gene,]
    x <- barplot(plot,ylim=c(0,max(plot+sd)*1.25),ylab="",col = mycolors,las=2)
    ylab <- ifelse(rpkm,"RPKM mean","Mean normalized read counts")
    mtext(ylab,side = 2,line = 4,cex = .6)
    arrows(x0 = x,y0 = plot,x1 = x,y1 = plot+sd,length = .1,angle = 90)
    title(main = gene)
    
    padj <- padj.a[gene]
    ## if more than two conditions, then skip plotting errorbars  
    if(ncol(data)<3) {
      
      lab <- ifelse(test = is.na(padj) | padj > .05,yes = "NA",
                    no = ifelse(test = padj<.0001,yes = "***",
                                no = ifelse(test = padj<.01,yes = "**",
                                            no = ifelse(test = padj<.05,yes = "*",no = "NA"))))
      arrows(x0 = x[1],y0 = max(plot+sd)*1.1,x1 = x[2],y1 = max(plot+sd)*1.1,code=0)
      text(x = x[1]+((x[2]-x[1])/2),y = max(plot+sd)*1.2,labels = lab,cex = 1)
    }
  }
  par(mar=c(4,4,4,4))
}
