# deseqAbstraction

Easy analysis and visualization of featureCount data in R. 

- Object-oriented analysis of RNA-seq data in R
- Wrapper for DESeq2
- deseqAbs object to ease and automate every analysis and visualization of RNA-seq data
- <b>Main usage tutorial: Usage.pdf or Usage.html in deseqAbstraction (current) directory. </b>
- More usage tutorials found in usageTutorials/deseqAbs repository!


## Setup of workspace

```{r setup workspace, message = FALSE}
# clean environment
rm(list=ls())
gc()

# install and load packages

# load deseqAbstraction
library(devtools)
install_github("perllb/deseqAbstraction")
library(deseqAbstraction)
library(DESeq2)
```


## Create deseqAbs object 
 Define path to count file
```{r, create object}
path <- "/home/pbrattaas//Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/Quant/hg38.gencode.exon.primary.txt"
```

## Define colData (describe samples)
`condition` and `samples` needed.
```{r }
colData <- data.frame(condition=c(rep("CTR",3),rep("KO",3)),
                      samples=c(107,108,109,110,111,112))
```

## Initialize deseqAbs object
```{r }
dnmt <- deseqAbs$new(name="DNMT1KO",filename=path,colData = colData)
```

## Run DESeq pipeline
 - You can try to do this in one function $fullAuto(). 
 - Full auto has four main steps:
 1. create DEseq object (access with $deseq)
 2. Do diffex analysis (access with $test)
 3. Do varianceStablizingTransformation (access with $VST)
 4. Do RPKM normalization (access with $rpkm)

However, for this automation, it is critical that your deseqAbs object has
- the raw featureCount file in the correct format
- the colData data.frame with at least one column called "condition", which contain the grouping of your samples, that will later be used for diffex anaysis

```{r, automation }
# run the automated pipe from creating deseq object to normalization (rpkm, vst) and diffex analysis!
dnmt$fullAuto()
```

## Basic QC 
```{r, qc}
dnmt$sampleQC()
```


## Custom diffex analysis
```{r}
## By default, DESeq2 will make diffex test between two conditions in your colData. To see which conditions are tested, type
dnmt$test$Default
## here you can see the two conditions, in this case: KO vs CTR
```

If you want to make a custom diffex test
```{r}
## If you have more conditions, do the following to test between two other conditions
# get the names of the conditions 
conds <- gsub(pattern = "condition",replacement = "",resultsNames(dnmt$deseq))
conds
# list condition names with index
cbind(conds,1:length(conds))

# Test KO vs. CTR
c1 <- conds[3] # since KO is the 3rd entry of the conds vector
c2 <- conds[2] # since CTR is the 2nd entry of the conds vector
# make diffex analysis , and give it a name (without space is the best) that makes it easy to access later
dnmt$makeDiffex(name = "KOvsWT",c1 = c1,c2 = c2)
## access the test data with 
dnmt$test$KOvsWT
## dnmt$test is a list of diffex objects, which can be retrieved one by one
## type dnmt$test to get a list of the diffex objects created
```

## Visualize data
### PCA 

```{r}
# Default PCA plot 
dnmt$pca()

# plot PCA without labels on points
PCAplotter(dat = dnmt$VST,title = "NCBI top 5000",ntop = 5000,color = dnmt$colData$condition,shape=dnmt$colData$condition)

# plot PCA with labels
PCAplotter(dat = dnmt$VST,title = "NCBI top 5000",ntop = 5000,color = dnmt$colData$condition,shape=dnmt$colData$condition,label = dnmt$colData$condition)
```

### sample-to-sample distance
```{r }
dnmt$sampleToSample()
```

### maPlot
```{r}
# make maPlot: color genes if they have p-adj and log2fc below/above given values (p and l) 
maPlot(test = dnmt$test$KOvsWT,"DNMT1 KO","CTR",l = .5,p = 1e-5)

# set id = T to label points
maPlot(test = dnmt$test$KOvsWT,"DNMT1 KO","CTR",l = .5,p = 1e-5,id=T)
# now click on the points you want, and press ESC when done
```

### Volcano plot

```{r}
# generate volcanoplot. 
volcanoPlot(test = dnmt$test$KOvsWT)
# cut y-axis 
volcanoPlot(test = dnmt$test$KOvsWT,max = 100)
# set id = T to label points
volcanoPlot(test = dnmt$test$KOvsWT,max = 100,id = T)
# now click on the points you want, and press ESC when done
```

### Heatmaps
#### Most significant genes for a given test
```{r}
# plot top significant genes
mostSignificantHeat(data = assay(dnmt$VST),test = dnmt$test$KOvsWT)

# plot top significant genes, top 20
mostSignificantHeat(data = assay(dnmt$VST),test = dnmt$test$KOvsWT,ntop = 20)

# plot top significant genes, top 20, with annotation of columns
mostSignificantHeat(data = assay(dnmt$VST),test = dnmt$test$KOvsWT,ntop = 20,a1 = dnmt$colData$condition,n1 = "Treatment")

## plot top significant genes, top 70, with two annotation of columns
mostSignificantHeat(data = assay(dnmt$VST),test = dnmt$test$KOvsWT,ntop = 70,a1 = dnmt$colData$condition,n1 = "Treatment",a2 = dnmt$colData$samples,n2 = "sampleNames")
```

#### Most variable genes (over all conditions)
```{r }
# plot most variable genes
mostVariableHeat(data = assay(dnmt$VST))

# plot most variable genes
mostVariableHeat(data = assay(dnmt$VST),ntop = 100,a1 = dnmt$colData$condition,n1 = "Treatment")
```

#### Heat a set of genes
```{r}
genes <- c("LPPR4", "NEFH", "BOC", "TRPC6", "GPRIN1", "ST8SIA2", "PTK2B", "UCHL1", "CTNNA1", "PLCG1", "SEMA6C", "DSCAML1")

# set sd cutoff > 0 to be sure not to get error
heatGenes(data = assay(dnmt$VST),genes = genes,sd = .01,a1 =  dnmt$colData$condition,n1 = "condition")

```

#### Mean Plot (mean of condition A vs B)
```{r}
## exp: mean normalized reads for each condition
## test: diffex test between the conditions
## c1 and c2: specify the name of each condition
meanPlot(exp = dnmt$baseMean$Mean,test = dnmt$test$Default,c1="DNMT1 KO",c2="CTR",p = 1e-5)
## the same, but plot in FPKM instead
meanPlot(exp = dnmt$FPKMMean$Mean,test = dnmt$test$Default,c1="DNMT1 KO",c2="CTR",p = 1e-5)

```

#### B barplots
```{r}
## give a vector of any length, with the genes you want to plot
genes <- c("DNMT1","DNMT3A","TRIM28")
dnmt$meanBars(genes = genes )
# to plot in FPKM instead of mean normalized reads:
dnmt$meanBars(genes = genes,FPKM = T )

```

## Get data 

### Get normalized expression counts of all genes 
```{r }
head(dnmt$normCounts)

```
### Get average normalized expression in conditions 
```{r}
# This is done automatically in $fullAuto(), but can also be done manually like this
dnmt$getAverage()

# Retrieve :
head(dnmt$baseMean$Mean)
head(dnmt$baseMean$SD)
head(dnmt$baseMean$SE)

```
### Get FPKM
```{r }
# done automatically by fullAuto(), but can be done manually by:
dnmt$makeFPKM()

# retrieve:
head(dnmt$FPKM)
```

### Get average FPKM in conditions
```{r }
# done automatically by fullAuto(), but can be done manually by:
dnmt$getAverageFPKM()

# retrieve:
head(dnmt$FPKMMean$Mean)
head(dnmt$FPKMMean$SE)
head(dnmt$FPKMMean$SD)

```

### Get significant genes
#### Get diffex data
```{r}
# Get diffex-data of genes with p-adj < 0.001 and log2FC more/less than 0.2
sign <- getSign(x = dnmt$test$KOvsWT,p = .001,l = .2)
up <- sign$up
down <- sign$down

head(up)
head(down)
```

#### Get names only
```{r }
# Get names of genes with p-adj < 0.01 and log2FC more/less than .1
signID <- getSignName(x = dnmt$test$KOvsWT,p = .01,l = .2)
upID <- signID$up
downID <- signID$down

head(upID)
head(downID)
```

#### Get most variable genes
```{r }
var.50 <- dnmt$getVariable(ntop = 50)
var.50
```

## Get genes close to features (in bed file)
```{r}
## takes as input a bed table
## E.g. a SVA-F elements
## NB: Must have col.names as defined here col.names = c("Chr","Start","End","ID",".","Strand")
svaf <- read.delim(file = "~/genomicData/Repeats/hg38/SVA/SVA_F.hg38.bed",col.names = c("Chr","Start","End","ID",".","Strand"),header=F)
head(svaf)

# calculate genes within 15000bp from SVA-F! 
# By default, this function search in Gencode v27 (hg38) annotation. This can be changed with setting a = <bedfile>.
close.svaf <- closeGenes(b=svaf,d = 15000)

head(close.svaf)
```


