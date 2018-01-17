#' @name meanBar
#' @description barplot of gene expression
#' @param deseqAbs: deseqAbs object with FPKMMean and test data
#' @param genes: a vector of gene IDs to be extracted from rownames of data
#' @param cond: a vector of conditions to plot
#' @param FPKM: Set to TRUE if FPKM should be plotted instead of baseMean
#' @param points: Set to TRUE if points should be added
#' @title meanBar: Barplot of your genes of interest!
#' @export meanBar
#' @examples
#' genes <- c("DNMT1","TRIM28","PAX6","DCX","SOX2","AGO2")
#' meanBar(dnmt,genes)

meanBar <- function(deseqAbs,genes,cond=NULL,FPKM=FALSE,points=FALSE) {
  
  cat(">>> meanBar: plot your genes:\n")
  cat(">>",genes,"\n")
  
  ## set graphical area
  row <- ifelse(test = sqrt(length(genes))%%1 > .5,yes = floor(sqrt(length(genes)))+1,no = floor(sqrt(length(genes)))) 
  par(mfrow=c(row,ceiling(sqrt(length(genes)))),mar=c(7,5,4,4))
  library(RColorBrewer)
  cols <- colorRampPalette(brewer.pal(8, "Greys"))
  
  # if no condition defined, run default
  if ( is.null(cond) ) {
    if(FPKM) {
      data <- deseqAbs$FPKMMean$Mean
      se.a <- deseqAbs$FPKMMean$SE
    } else {
      data <- deseqAbs$baseMean$Mean 
      se.a <- deseqAbs$baseMean$SE
    }
    mycolors <- cols(length(unique(deseqAbs$colData$condition)))
    if(is.null(deseqAbs$test$Default)){ 
      cat(">> To plot significance stars between the two conditions, Diffex must be run..\n")
      deseqAbs$makeDiffex() 
    }
    padj.a <- deseqAbs$test$Default$padj
    names(padj.a) <- rownames(deseqAbs$test$Default)
  } else {
    if(FPKM) {
      data <- deseqAbs$FPKMMean$Mean[,cond]
      se.a <- deseqAbs$FPKMMean$SE[,cond]
    } else {
      data <- deseqAbs$baseMean$Mean[,cond]
      se.a <- deseqAbs$baseMean$SE[,cond]
    }
    mycolors <- cols(length(cond))
    str <- paste()
    deseqAbs$makeDiffex(name='tmptest',c1=cond[1],c2=cond[2])
    padj.a <- deseqAbs$test$tmptest$padj
    names(padj.a) <- rownames(deseqAbs$test$tmptest)
  }
  
  # plot one for each gene 
  for(gene in genes) {
    
    mean <- data[gene,]
    se <- se.a[gene,]
    my_data <- data.frame(cond=colnames(data),mean=mean,se=se)
    p <- ggplot() +
      geom_bar(data=my_data, aes(y=mean,x=cond,ymin=mean-se,ymax=mean+se), stat="identity", width = 0.1) + 
      geom_errorbar(data=my_data, aes(y=mean,x=cond,ymin=mean-se,ymax=mean+se), width = 0.05) 
    
    if(FPKM){
      data.p <- melt(deseqAbs$FPKM[gene,])
      p2 <- p +
        geom_point(data=data.p,aes(y=value, x=X2))  
      ylab <- "FPKM mean"
    }else{
      data.p <- melt(deseqAbs$normCounts[gene,])
      p2 <- p +
        geom_point(data=data.p,aes(y=value, x=X2))  
      ylab <- "Mean normalized read counts"
    }
    
    p2 + theme_classic() + ylim(c(0,max(mean+se))) + labs(y=ylab,title=gene) 
    
    padj <- padj.a[gene]
    ## if more than two conditions, then skip plotting errorbars  
    if(ncol(data)<3) {
      
      lab <- ifelse(test = is.na(padj) | padj > .05,yes = "NA",
                    no = ifelse(test = padj<.0001,yes = "***",
                                no = ifelse(test = padj<.01,yes = "**",
                                            no = ifelse(test = padj<.05,yes = "*",no = "NA"))))
      arrows(x0 = x[1],y0 = max(plot+se)*1.1,x1 = x[2],y1 = max(plot+se)*1.1,code=0)
      text(x = x[1]+((x[2]-x[1])/2),y = max(plot+se)*1.2,labels = lab,cex = 1)
    }
  }
  par(mar=c(4,4,4,4))
}
