#' @name meanPlot
#' @description Plots scatter plot with log2(cond1) vs. log2(cond2), and marks significant genes
#' @param exp: data.frame with average expression of genes (rows) in conditions (cols). Two columns only! if you experiment has more conditions, then select columns first.
#' @param c1: string describing condition 1 in results deseq test
#' @param c2: string describing condition 2 in results deseq test
#' @title meanPlot compare gene expression
#' @export meanPlot
#' @example
#' exp <- getAverage(dds)
#' maPlot(exp = exp, c1 = "KO", c2 = "CTR" )

meanPlot <- function(exp,c1 = "condition 1",c2 = "condition 2",p=.5,l=0) {

  #color up and down sign..
  colVec <- ifelse(test = exp$padj<p,
                   yes=ifelse(test = exp$log2FoldChange>l,
                              yes="firebrick3",
                              no=ifelse(exp$log2FoldChange< -l,
                                        yes = "steelblue4",
                                        no = "black")),
                   no = "black")
  colVec[is.na(colVec)] <- "black" ## if NA make sure it's not counted as <p
  #size of points
  cexVec <- ifelse(test = exp$padj<p, yes = ifelse(test = (is.na(exp$padj)),yes = 0.15,no = 0.4), no= 0.15)

  sign <- getSign(x = exp,p = p,l = l)
  u <- nrow(sign$up)
  d <- nrow(sign$down)
  n <- nrow(exp) - u - d

  plot(log2(exp[,1]),log2(exp[,2]),
       col=colVec,
       cex=cexVec,
       pch=16,
       xlab=paste("log2(mean ",c1,")",sep=""),
       ylab=paste("log2(mean ",c2,")",sep=""))
  title(main=paste(c1," vs. ",c2,sep=""))
  mtext(text = paste("p-adj < ",p,", log2(fc) > ",l,sep=""),side = 3)
  legend("topleft",legend = c(paste("up (",u,")",sep=""),paste("down (",d,")",sep = ""),paste("not significant (",n,")",sep = "")),pch=16,col=c("firebrick3","steelblue4","black"),bty='n')

}


