#' @name meanBarGeneric
#' @description barplot of gene expression
#' @param tab: data.frame to plot!
#' @param cond: a vector of conditions to group data
#' @param points: Set to TRUE if points should be added
#' @param jitter: set to TRUE if jitter should be applied to points 
#' @param ylab: enter ylab for plot!
#' @title meanBarGeneric: Barplot of your genes of interest!
#' @export meanBarGeneric
#' @example meanBarGeneric(table,points=deeq$colData$cond,points=T)

meanBarGeneric <- function(tab,cond=NULL,points=FALSE,jitter=F,ylab='') {
  
  cond <- factor(cond)
  par(mar=c(6,6,6,6))
  
  if(is.null(tab)){ 
    stop(">>> ERROR: Insert tab!")
  } 
  if(!is.data.frame(tab)){ 
    stop('>>> ERROR: tab must be data frame! Each row is a feature, Columns are samples')
  }
  
  cat(">>> meanBarGeneric: plot your data frame!\n")
  
  if(is.vector(tab)){n.genes<-1} else { n.genes <- nrow(tab) }
  
  ## set graphical area and colors
  row <- ifelse(test = sqrt(n.genes)%%1 > .5,yes = floor(sqrt(n.genes))+1,no = floor(sqrt(n.genes))) 
  par(mfrow=c(row,ceiling(sqrt(n.genes))),mar=c(7,5,4,4))
  library(RColorBrewer)
  cols <- colorRampPalette(brewer.pal(8, "Greys"))
  mycolors <- cols(length(cond))
  
  data <-sapply( levels(cond), function(lvl) rowMeans( tab[,cond == lvl] ) )
  se.a <-sapply( levels(cond), 
                  function(lvl) apply( tab[,cond == lvl],1,
                                       function(dat) sd(dat)/sqrt(length(dat))))
  if(!is.matrix(data)) {
    data <- data.frame(t(data))
    se.a <- data.frame(t(se.a))
  }
  
  # plot one for each gene 
  for(gene in rownames(tab)) {
    
    ymax <- max(tab[gene,])
    plot <- data[gene,]
    se <- se.a[gene,]
    x <- barplot(plot,ylim=c(0,ymax*1.3),ylab="",col = mycolors,las=2,space = 0,main=gene)
    # If points set to T, then add points of each sample to plot
    if(points) {
      if(jitter){
        x.s <- jitter(as.numeric(cond)-x[1],factor=0.4)
      }else {
        x.s <- as.numeric(cond)-x[1]
      }
        
        points(x=x.s,y=tab[gene,],col="black",pch=16,cex=1.5)
        points(x=x.s,y=tab[gene,],col="white",pch=16)
    }
    
    mtext(ylab,side = 2,line = 4,cex = .6)
    arrows(x0 = x,y0 = plot,x1 = x,y1 = plot+se,length = .1,angle = 90)
    title(main = gene)
    
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
  par(mfrow=c(1,1))
}

