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
  par(mfrow=c(row,ceiling(sqrt(length(genes)))),mar=c(10,5,4,4))
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
    
    ymax <- ifelse(FPKM,yes=max(deseqAbs$FPKM[gene,]),no = max(deseqAbs$normCounts[gene,]))
    plot <- data[gene,]
    se <- se.a[gene,]
    x <- barplot(plot,ylim=c(0,ymax*1.3),ylab="",col = mycolors,las=2,space = 0)
    # If points set to T, then add points of each sample to plot
    if(points) {
      if(FPKM) {
        x.s <- jitter(rep(x,times=table(deseqAbs$colData$condition)),factor=0.4)
        points(x=x.s,y=deseqAbs$FPKM[gene,],col="black",pch=16,cex=1.3)
        points(x=x.s,y=deseqAbs$FPKM[gene,],col="white",pch=16)
      }else{
        x.s <- jitter(rep(x,times=table(deseqAbs$colData$condition)),factor=0.4)
        points(x=x.s,y=deseqAbs$normCounts[gene,],col="black",pch=16,cex=1.5)
        points(x=x.s,y=deseqAbs$normCounts[gene,],col="white",pch=16)
      }
    }
    ylab <- ifelse(FPKM,"FPKM mean","Mean normalized read counts")
    mtext(ylab,side = 2,line = 4,cex = .6)
    arrows(x0 = x,y0 = plot,x1 = x,y1 = plot+se,length = .1,angle = 90)
    title(main = gene)
    
    padj <- padj.a[gene]
    ## if more than two conditions, then skip plotting errorbars  
    if(ncol(data)<3) {
      
      lab <- ifelse(test = is.na(padj) | padj > .05,yes = "NA",
                    no = ifelse(test = padj<.0001,yes = "***",
                                no = ifelse(test = padj<.01,yes = "**",
                                            no = ifelse(test = padj<.05,yes = "*",no = "NA"))))
      arrows(x0 = x[1],y0 = ymax*1.1,x1 = x[2],y1 = ymax*1.1,code=0)
      text(x = x[1]+((x[2]-x[1])/2),y = ymax*1.2,labels = lab,cex = 1)
    }
  }
}

