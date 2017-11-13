#' @name volcanoPlot
#' @description Generates volcano plot and marks significant genes. Scales the y axis if promted
#' @param test: output object from results() of DESeq2
#' @param p: cutoff of significance (p-adj)
#' @param title: string describing plot
#' @param id: If TRUE, you can identify points and label their names. FALSE by default
#' @title volcano plot
#' @export volcanoPlot
#' @example
#' test <- results(dds,contrast = c("condition","genex-KO","WT"))
#' volcanoPlot(test = test, max=9,p=.001,title="transcriptome changes upon KO")

volcanoPlot <- function(test,max=NULL,p=.05,title="transcriptome changes",id=F) {

  # get genes with p-adj 0. (these are so significant that there are not a low enough number for the p-adj in R)
  t <- which(is.infinite(-log10(test$padj)))
  # give these genes p-adj = the lowest non-zero number in R
  if(length(t)>0) { test[t,]$padj <- min(test[!is.na(test$padj) & test$padj>0,]$padj) }

  # create -log10 ( padj ) yaxis
  yaxis <- -log10(test$padj)

  cex <- ifelse(test = yaxis > p, yes = 0.6, no = 0.3)
  col <- ifelse(test = test$padj < p, yes = ifelse(test = test$log2FoldChange >
                                                     0, yes = "firebrick3", no = ifelse(test$log2FoldChange <
                                                                                          0, yes = "steelblue4", no = "black")), no = "black");col[is.na(col)] <- "black"

  pch <- 16
  # if user specifies max ylim
  if (!is.null(max)) {
    yaxis[yaxis > max] <- max + 0.1
    pch <- ifelse(test = (yaxis == max + 0.1), yes = 18,
                   no = 16)
  }

  # xaxis scaling
  log2FC.max <- max(abs(test[is.finite(test$log2FoldChange),]$log2FoldChange))

  plot(x = test$log2FoldChange, y = yaxis, main = title, col = col,
       pch = pch, cex = cex, ylab = "-log10(p-adj)", xlab = "log2FC",
       xlim = c(-log2FC.max, log2FC.max), ylim = c(0, max(yaxis[is.finite(yaxis)]) * 1.05))

  # subtitle of plot
  mtext(text = paste("p-adj < ", p, sep = ""), side = 3)

  # get sign changed genes
  sign <- getSignName(x = test, p = p, l = 0)
  u <- length(sign$up)
  d <- length(sign$down)
  legend("topleft", legend = c(paste("up (", u, ")", sep = ""),
                               paste("down (", d, ")", sep = ""), "not significant"),
         pch = 16, col = c("firebrick3", "steelblue4", "black"),
         bty = "n", cex = 0.7)
  abline(v = 0, lty = 2)

  if(id==T) {

    identify(x = test$log2FoldChange, y = yaxis,labels = rownames(test))

  }

}

