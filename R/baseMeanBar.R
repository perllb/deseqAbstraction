#' @name baseMeanBar
#' @description barplot of gene expression
#' @param deseqAbs: deseqAbs object with baseMean and test data
#' @param genes: a vector of gene IDs to be extracted from rownames of data
#' @title baseMeanBar: Barplot of your genes of interest!
#' @export baseMeanBar
#' @examples
#' deseqAbs <- dnmt
#' genes <- c("DNMT1","TRIM28","PAX6","DCX","SOX2","AGO2")
#' baseMeanBar(dnmt,genes)

baseMeanBar <- function(deseqAbs,genes) {

  ## set graphical area
  par(mfrow=c(floor(sqrt(length(genes))),ceiling(sqrt(length(genes)))))

  cols <- c("white","black","grey","blue","green")

  for(gene in genes) {

    plot <- deseqAbs$baseMean$Mean[gene,]
    sd <- deseqAbs$baseMean$SD[gene,]
    x <- barplot(plot,ylim=c(0,max(plot+sd)*1.5),ylab="mean normalized read counts",col = cols)
    arrows(x0 = x,y0 = plot,x1 = x,y1 = plot+sd,length = .1,angle = 90)
    title(main = gene)

    padj <- deseqAbs$test$Default[gene,]$padj

    lab <- ifelse(test = is.na(padj) | padj > .05,yes = "NA",
                  no = ifelse(test = padj<.0001,yes = "***",
                              no = ifelse(test = padj<.01,yes = "**",
                                          no = ifelse(test = padj<.05,yes = "*",no = "NA"))))
    arrows(x0 = x[1],y0 = max(plot)*1.4,x1 = x[2],y1 = max(plot)*1.4,code=0)
    text(x = x[1]+((x[2]-x[1])/2),y = max(plot)*1.5,labels = lab,cex = 2)


  }
}
