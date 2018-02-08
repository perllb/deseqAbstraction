#' @name meanPlot
#' @description Plots scatter plot with log2(cond1) vs. log2(cond2), and marks significant genes
#' @param deseqAbs: deseqAbs object
#' @param c1: string describing condition 1 in results deseq test
#' @param c2: string describing condition 2 in results deseq test
#' @param id: If TRUE, you can identify points and label their names. FALSE by default#'
#' @title meanPlot compare gene expression
#' @export meanPlot
#' @examples
#' exp <- getAverage(dds)
#' maPlot(exp = exp, c1 = "KO", c2 = "CTR" )

meanPlot <- function(deseqAbs,c1 = NULL,c2 = NULL,p=.05,l=0,id=F) {

  if(is.null(c1) || is.null(c2)){
    stop("Set conditions!")
  }
  
  name <- paste(c1,c2,sep = ".vs.")
  deseqAbs$makeDiffex(name = name,c1 = c1,c2 = c2)
  test <- deseqAbs$test$name
  
  sign <- getSignName(x = test,p = p,l = l)
  u <- sign$up
  d <- sign$down
  n <- nrow(exp) - length(u) - length(d)

  exp <- deseqAbs$baseMean$Mean[,c(c2,c1)]
  
  #color up and down sign..
  colVec <- ifelse(test = (rownames(exp) %in% u),
                   yes = "firebrick3",
                   no = ifelse(test = (rownames(exp) %in% d),
                               yes = "steelblue4", no ="black"))
  colVec[is.na(colVec)] <- "black" ## if NA make sure it's not counted as <p
  #size of points
  cexVec <- ifelse(test = (rownames(exp) %in% u),
                   yes = .5,
                   no = ifelse(test = (rownames(exp) %in% d),
                               yes = .5, no = .3))

  plot(log2(exp[,1]),log2(exp[,2]),
       col=colVec,
       cex=cexVec,
       pch=16,
       xlab=paste("log2(mean ",c1,")",sep=""),
       ylab=paste("log2(mean ",c2,")",sep=""))
  title(main=paste(c2," vs. ",c1,sep=""))
  mtext(text = paste("p-adj < ",p,", log2(fc) > ",l,sep=""),side = 3)
  legend("bottomright",legend = c(paste("up (",length(u),")",sep=""),paste("down (",length(d),")",sep = ""),paste("not significant (",n,")",sep = "")),pch=16,col=c("firebrick3","steelblue4","black"),cex=.5)

  if(id==T) {

    identify(log2(exp[,1]),log2(exp[,2]),labels = rownames(exp))

  }

}
