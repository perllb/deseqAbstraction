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
  par(mfrow=c(floor(sqrt(length(genes))),ceiling(sqrt(length(genes)))))
  library(RColorBrewer)
  cols <- colorRampPalette(brewer.pal(8, "Dark2"))
  mycolors <- cols(length(unique(deseqAbs$colData$condition)))
  
  # if no condition defined, run default
  if ( is.null(cond) ) {
    # plot one for each gene 
    for(gene in genes) {
  
      plot <- deseqAbs$baseMean$Mean[gene,]
      sd <- deseqAbs$baseMean$SD[gene,]
      x <- barplot(plot,ylim=c(0,max(plot+sd)*1.5),ylab="mean normalized read counts",col = cols)
      arrows(x0 = x,y0 = plot,x1 = x,y1 = plot+sd,length = .1,angle = 90)
      title(main = gene)
  
      ## if more than two conditions, then skip plotting errorbars  
      if(ncol(deseqAbs$baseMean$Mean)<1) {
        padj <- deseqAbs$test$Default[gene,]$padj
    
        lab <- ifelse(test = is.na(padj) | padj > .05,yes = "NA",
                      no = ifelse(test = padj<.0001,yes = "***",
                                  no = ifelse(test = padj<.01,yes = "**",
                                              no = ifelse(test = padj<.05,yes = "*",no = "NA"))))
        arrows(x0 = x[1],y0 = max(plot)*1.4,x1 = x[2],y1 = max(plot)*1.4,code=0)
        text(x = x[1]+((x[2]-x[1])/2),y = max(plot)*1.5,labels = lab,cex = 1)
      }
    }
  }
  else {
    data <- deseqAbs$baseMean$Mean[,cond]
    # plot one for each gene 
    for(gene in genes) {
      
      plot <- data[gene,]
      sd <- data[gene,]
      x <- barplot(plot,ylim=c(0,max(plot+sd)*1.5),ylab="mean normalized read counts",col = cols)
      arrows(x0 = x,y0 = plot,x1 = x,y1 = plot+sd,length = .1,angle = 90)
      title(main = gene)
      
      ## if more than two conditions, then skip plotting errorbars  
      if(ncol(data)<1) {
        padj <- deseqAbs$test$Default[gene,]$padj
        
        lab <- ifelse(test = is.na(padj) | padj > .05,yes = "NA",
                      no = ifelse(test = padj<.0001,yes = "***",
                                  no = ifelse(test = padj<.01,yes = "**",
                                              no = ifelse(test = padj<.05,yes = "*",no = "NA"))))
        arrows(x0 = x[1],y0 = max(plot)*1.4,x1 = x[2],y1 = max(plot)*1.4,code=0)
        text(x = x[1]+((x[2]-x[1])/2),y = max(plot)*1.5,labels = lab,cex = 1)
      }
    }
  }
}
