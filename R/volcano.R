#' @name volcanoPlot
#' @description Generates volcano plot and marks significant genes. Scales the y axis if promted
#' @param test: output object from results() of DESeq2
#' @param p: cutoff of significance (p-adj)
#' @param title: string describing plot
#' @title volcano plot
#' @export volcanoPlot
#' @example
#' test <- results(dds,contrast = c("condition","genex-KO","WT"))
#' volcanoPlot(test = test, max=9,p=.001,title="transcriptome changes upon KO")

volcanoPlot <- function(test,max=NULL,p=.5,title="transcriptome changes") {

  # get values of y axis
  yaxis <- -log10(test$padj)
  # size of dots according to significane
  cex <- ifelse(test = yaxis<p,yes = .9,no = .4)
  #color up and down sign..
  col <- ifelse(test = test$padj<p,
                   yes=ifelse(test = test$log2FoldChange>0,
                              yes="firebrick3",
                              no=ifelse(test$log2FoldChange<0,
                                        yes = "steelblue4",
                                        no = "black")),
                   no = "black")
  col[is.na(col)] <- "black" ## if NA make sure it's not counted as <p

  # if user puts a max on -log10(padj), adjust plot
  if(!is.null(max)) {
    # change y value of points that are above max to y = max + 0.1 instead)
    yaxis[yaxis>max] <- max+0.1
    # change point type of those points
    pchY <- ifelse(test = (yaxis==max+0.1),yes = 18,no = 16)
  }

  log2FC.max <- max(abs(test$log2FoldChange))

  plot(x = test$log2FoldChange,
        y = yaxis,
       main=title,
       col=col,
       pch=pchY,
       cex=cex,
       ylab = "-log10(p-adj)",
       xlab = "log2FC",
       xlim = c(-log2FC.max*1.1,log2FC.max*1.1),
       ylim = c(0,max(yaxis)*1.2))

  sign <- getSign(x = test,p = p,l = 0)
  u <- nrow(sign$up)
  d <- nrow(sign$down)

  legend("topleft",legend = c(paste("up (",u,")",sep=""),paste("down (",d,")",sep = ""),"not significant"),pch=16,col=c("firebrick3","steelblue4","black"),bty='n')
  abline(v = 0,lty=2)

}
