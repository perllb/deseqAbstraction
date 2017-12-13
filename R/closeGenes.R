#' @name closeGenes
#' @description Get genes (or other features) close to some other features (.bed). TSS used.
#' @param a: by default the gencode proteincoding transcripts will be used. Has to be a .bed position of all genes of interest. Must have $Chr, $Start, $End $ID $Strand column names! TSS of the transcrips will be used.
#' @param b: a .bed position of all features of interest. Must have $Chr, $Start, $End $ID $Strand column names! TSS of the features will be used
#' @param d: how many basepairs from features B should a feature A be to be included? Default = 10000bps
#' @title closeGenes : Get close features (genes )
#' @export closeGenes
#' @examples
#' closeGenes <- enesClose(a = gencode.transcripts, b = L1s.up, d = 50000 )

closeGenes <- function(a=NULL,b,d=10000) {

  # If a is NULL, then get gencode proteincoding transcript bed file
  if(is.null(a)){
    cat("> Downloading the gencode v25 annotation of protein-coding transcripts..\n")
    library(RCurl)
    a <- read.delim(text = getURL(paste("https://raw.githubusercontent.com/perllb/deseqAbstraction/master/data/gencode.v25.annotation.proteinCoding.Transcript.bed",sep = "")),header=F)
    colnames(a) <- c("Chr","Start","End","ID",".","Strand")
    cat("> Gencode v25 protein-coding transcripts .bed downloaded!\n")
    cat("> (By default, this bed file is used, as it provides coordinates for each protein-coding transcript of gencode.v25 genes.. If you wish to use an alternative gene file, it can be provided as an a = <annotation file> argument to this function. Beware of format..\n")
  }
  if('Chr' %in% colnames(a) & 'Strand' %in% colnames(a) & 'Start' %in% colnames(a) & 'End' %in% colnames(a) & 'ID' %in% colnames(a)){
    if('Chr' %in% colnames(b) & 'Strand' %in% colnames(b) & 'Start' %in% colnames(b) & 'End' %in% colnames(b) & 'ID' %in% colnames(b)){

      ## data frame to build
      allclose <- data.frame()
      # iterate through each feature -> get genes with TSS <50kb away)
      cat(">> Getting genes with TSS within",d,"bps from each given feature. This might take a couple of minutes. Patience..\n")

      for (i in 1:nrow(b)){

        # get position of current L1
        feat.tss <- ifelse(as.character(b$Strand[i])=="+",as.numeric(as.character(b$Start[i])),as.numeric(as.character(b$End[i])))
        # get genes on same chromosome
        genes.chr <- a[a$Chr == as.character(b$Chr[i]),]
        genes.chr.tss <- ifelse(genes.chr$Strand=="+",yes = genes.chr$Start,no = genes.chr$End)
        # get genes with TSS < d from L1 features
        closeGenes <-genes.chr[abs(as.numeric(as.character(genes.chr.tss))-feat.tss)<d,]

        atss <- ifelse(closeGenes$Strand=="+",yes = closeGenes$Start,no = closeGenes$End)
        closedf <- data.frame(A_ID=trimws(closeGenes$ID),A_chr=closeGenes$Chr,A_Start=closeGenes$Start,
                              A_End=closeGenes$End,A_TSS=atss,A_Strand=closeGenes$Strand,
                              B_ID=rep(as.character(b$ID[i]),nrow(closeGenes)),
                              B_TSS=rep(feat.tss,nrow(closeGenes)),B_Strand=rep(b$Strand[i],nrow(closeGenes)),Distance=atss-feat.tss)
        allclose <- rbind(allclose,closedf)

      }

      cat("- ..complete! Genes A with TSS within",d,"bps from each features B TSS is computed.\n")
      return(allclose)
    }else{
      print("> ERROR: Colnames of your 'b' table/data.frame must be with Chr, Start, End, Strand, ID colnames")
      print("> ERROR: Please change your column names.")
    }
  }else{
    print("> ERROR: Colnames of your 'a' table/data.frame must be with Chr, Start, End, Strand, ID colnames")
    print("> ERROR: Please change your column names.")
  }
}
