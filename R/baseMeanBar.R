#' @name baseMeanBar
#' @description barplot of gene expression
#' @param deseqAbs: deseqAbs object with baseMean and test data
#' @param genes: a vector of gene IDs to be extracted from rownames of data
#' @param cond: a vector with condition names to use (must match column names of baseMean df)
#' @title baseMeanBar: Barplot of your genes of interest!
#' @export baseMeanBar
#' @examples
#' deseqAbs <- dnmt
#' genes <- c("DNMT1","TRIM28","PAX6","DCX","SOX2","AGO2")
#' baseMeanBar(dnmt,genes)

baseMeanBar <- function(deseqAbs,genes,cond=NULL) {
  
  ## set graphical area
  row <- ifelse(test = sqrt(length(genes))%%1 > .5,yes = floor(sqrt(length(genes)))+1,no = floor(sqrt(length(genes)))) 
  par(mfrow=c(row,ceiling(sqrt(length(genes)))))
  library(RColorBrewer)
  cols <- colorRampPalette(brewer.pal(8, "Greys"))

  # if no condition defined, run default
  if ( is.null(cond) ) {
    data <- deseqAbs$baseMean$Mean
    mycolors <- cols(length(unique(deseqAbs$colData$condition)))
    padj <- deseqAbs$test$Default[gene,]$padj 
  } else {
    data <- deseqAbs$baseMean$Mean[,cond]
    mycolors <- cols(length(cond))
    str <- paste()
    deseqAbs$makeDiffex(name='tmptest',c1=cond[1],c2=cond[2])
    padj <- deseqAbs$test$tmptest$padj
  }

  # plot one for each gene 
  for(gene in genes) {

    plot <- data[gene,]
    sd <- data[gene,]
    x <- barplot(plot,ylim=c(0,max(plot+sd)*1.25),ylab="mean normalized read counts",col = mycolors,las=2)
    arrows(x0 = x,y0 = plot,x1 = x,y1 = plot+sd,length = .1,angle = 90)
    title(main = gene)

    ## if more than two conditions, then skip plotting errorbars  
    if(ncol(data)<3) {
      
      lab <- ifelse(test = is.na(padj) | padj > .05,yes = "NA",
                    no = ifelse(test = padj<.0001,yes = "***",
                                no = ifelse(test = padj<.01,yes = "**",
                                            no = ifelse(test = padj<.05,yes = "*",no = "NA"))))
      arrows(x0 = x[1],y0 = max(plot+sd)*1.1,x1 = x[2],y1 = max(plot+sd)*1.1,code=0)
      text(x = x[1]+((x[2]-x[1])/2),y = max(plot)*1.2,labels = lab,cex = 1)
    }
  }
}
