#' @name heatGenes
#' @description Plots heatmap showing the expression of the genes of given gene-list and having sd above given value
#' @param data: deseq matrix data (with varianceStablizing, rlog etc - eg. assay(vsd)
#' @param genes: a vector of gene IDs to be extracted from rownames of data
#' @param a1: annotation of the samples
#' @param a2: annotation of the samples
#' @param n1: name of annotation in a1
#' @param n2: name of annotation in a2
#' @param sd: standard deviation cutoff for heatmap
#' @param z: logical operator. if TRUE, normalize by row (z-score). if FALSE (default), then plot normalized expression values
#' @param cluster_col: logical operator: if TRUE (default), cluster columns, if FALSE do not cluster columns
#' @param k: number of clusters in k-means
#' @param cutreeR: number of cuts in rowtree
#' @param cutreeC: number of cuts in Coltree
#' @title HeatGenes: Heatmap of your genes of interest!
#' @export heatGenes
#' @examples
#' vst <- varianceStabilizingTransformation(dds)
#' genes <- c("TRIM28","DNMT1","ZNF52")
#' heatGenes(data = assay(vst),genes = genes, sd = 1, a1 = colData$cellLine, a2 = colData$treatment,n1 = "Cell Line",n2 = "Treatment",z = T)


heatGenes <- function(data,genes,a1=NULL,a2=NULL,n1=NULL,n2=NULL,sd=.001,z=FALSE,cluster_col=T,k=NA,cutreeR=1,cutreeC=1,redBlue=T,breaks=NA) {

  library(pheatmap)
  library(graphics)
  library(RColorBrewer)
  
  if(!is.matrix(data) & !is.data.frame(data)) {
    data <- assay(data)
  }

  if(!is.matrix(data) & !is.data.frame(data)) {
    cat("ERROR: Data is not in correct format. Must be matrix or DESeq object")
  } else {

    # Set annotation colors (9 colors)
    a1col <- c( "#808000", 	"#FFD700",	"#20B2AA",	"#D2691E","#BC8F8F",	"#FFE4B5", "#BD1212",	"#00008B")
            #  olive green, gold,ligth sea, forest green, rosy brown, moccasin  , red ,   dark blue, 
    a2col <- c("#FAEBD7",      "#8B4513" ,    "#B0C4DE",         "#B0C4DE",  "#000080" ,"	#6495ED",   "#008080","#00FF00", 	"#F0E68C")
            # antique white, saddle brown , ligth steel blue,  slate brue,  navy , corn flower blue, teal  ,  lime,      khaki
    ## Change to RdYlBu (if RedBlue == F)
    heatCol <- ifelse(redBlue,yes = "RdBu",no = "RdYlBu")
    
    if(!is.na(breaks[1])) {
      heatScaleCol <- rev(colorRampPalette(brewer.pal(10,heatCol))(length(breaks)))
    }else {
      heatScaleCol <- rev(colorRampPalette(brewer.pal(10,heatCol))(200))
    }
    
    genes.exp <- getGenes(data = data,genes = genes)
    rownames(genes.exp) <- make.names(genes.exp[,1],unique = T)
    genes.exp <- genes.exp[,-1]

    ## get standard deviation of each gene
    sd.exp <- apply(genes.exp,1,sd)

    # data to plot
    plotData <- genes.exp[sd.exp > sd,]

    ## Show rownames? only if less that 41 genes are plotted
    rowShow <- T; if(nrow(plotData)>100) { rowShow <- F }
    print(paste("There are ",nrow(plotData)," genes in you gene-set with sd > ",sd,".",sep=""))

    ## scale by row?
    scale <- 'none' ; if(z) { scale <- 'row' }

    if (!is.null(a1) & is.null(a2)) {

      df <- data.frame(Var1 = factor(a1))
      rownames(df) <- colnames(data)
      colnames(df) <- n1
      names(plotData) <- rownames(df)
      
      mycolors <- a1col[1:length(unique(a1))]
      names(mycolors) <- unique(a1)
      mycolors <- list(a = mycolors)
      names(mycolors) <- n1
  
      pheatmap(plotData, annotation_col = df, breaks = breaks, annotation_colors = mycolors,border_color = NA,kmeans_k=k,cutree_rows=cutreeR,cutree_cols=cutreeC, cluster_rows = T, show_rownames = rowShow, cluster_cols = cluster_col,scale = scale,color = heatScaleCol)

    } else if (!is.null(a1) & !is.null(a2)) {

      df <- data.frame(Var1 = factor(a1),Var2 = factor(a2))
      rownames(df) <- colnames(data)
      colnames(df) <- c(n1,n2)
      names(plotData) <- rownames(df)
      
      mycolors <- a1col[1:length(unique(a1))]
      names(mycolors) <- unique(a1)

      mycolors2 <- a2col[1:length(unique(a2))]
      names(mycolors2) <- unique(a2)

      mycolors <- list(a = mycolors,b = mycolors2)
      names(mycolors) <- c(n1,n2)

     pheatmap(plotData, annotation_col = df, annotation_colors = mycolors,border_color = NA,kmeans_k=k,cutree_rows=cutreeR,cutree_cols=cutreeC,cluster_rows = T, show_rownames = rowShow, cluster_cols = cluster_col,scale = scale,color = heatScaleCol,breaks = breaks)


    } else if (!is.null(a2) & is.null(a1)) {

      a1<-a2
      df <- data.frame(Var1 = factor(a1))
      rownames(df) <- colnames(data)
      colnames(df) <- n1
      names(plotData) <- rownames(df)
      
      mycolors <- a1col[1:length(unique(a1))]
      names(mycolors) <- unique(a1)
      mycolors <- list(a = mycolors)
      names(mycolors) <- n1

      pheatmap(plotData, annotation_col = df, annotation_colors = mycolors,border_color = NA,kmeans_k=k,cutree_rows=cutreeR,cutree_cols=cutreeC,cluster_rows = T, show_rownames = rowShow, cluster_cols = cluster_col,scale = scale,color = heatScaleCol,breaks = breaks)

    } else {
      pheatmap(plotData, cluster_rows = T,  show_rownames = rowShow,border_color = NA,kmeans_k=k,cutree_rows=cutreeR,cutree_cols=cutreeC, cluster_cols = cluster_col,scale = scale,color = heatScaleCol,breaks = breaks)
    }

    return(rownames(plotData))
  }
}
