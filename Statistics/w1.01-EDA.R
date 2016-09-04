## ----global_palette, results = 'asis'------------------------------------
rm(list=ls())
tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

## ----global_options,warning=FALSE,message=FALSE--------------------------
## see ch. 10 Hooks of Xie's knitr book
library(knitr)
knit_hooks$set(setPch = function(before, options, envir) {
    if(before) par(pch = 19)
})
opts_chunk$set(setPch = TRUE)

## ----global_plot,warning=FALSE, message=FALSE----------------------------
knitr::opts_chunk$set(fig.width=5, fig.height=5, size="footnotesize",
                      warning=FALSE, message=FALSE)
knitr::knit_hooks$set(small.mar = function(before, options, envir) {
    if (before) graphics::par(mar = c(5,5,1.5,1))
})

## ----load_hidden, echo=FALSE, results="hide", warning=FALSE--------------
suppressPackageStartupMessages({
    library(devtools)
    library(Biobase)
    library(RSkittleBrewer)
    library(gplots)
    library(dplyr)
    library(AnnotationDbi)
})

## ----load----------------------------------------------------------------
library(gplots)
library(devtools)
library(Biobase)
library(RSkittleBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)

## ----install_packages, eval=FALSE----------------------------------------
#  install.packages(c("devtools","gplots"))
#  source("http://www.bioconductor.org/biocLite.R")
#  biocLite(c("Biobase","org.Hs.eg.db","AnnotationDbi"))
#  biocLite("alyssafrazee/RSkittleBrewer")

## ----pretty, eval=FALSE--------------------------------------------------
#  library(RSkittleBrewer)
#  # Make the colors pretty
#  trop = RSkittleBrewer("tropical")
#  palette(trop)
#  par(pch=19)

## ----load_data-----------------------------------------------------------
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
ls()
bm = bodymap.eset # simplify the name
pdata=pData(bm) # extract phenotipe data
edata=exprs(bm) # extract expression data
fdata = fData(bm) #extract feature data
ls()              # see the new variables

## ----tables--------------------------------------------------------------
table(pdata$gender) # tables for phenotype data
table(pdata$race)   # tables for phenotype data
table(pdata$gender,pdata$race) # cross-tables

## ----summary-------------------------------------------------------------
summary(edata) #summary of the distribution for every column

## ----missing-------------------------------------------------------------
# Use option useNA to include NA's in table
table(pdata$age)   # table looks fine: No NA
table(pdata$age,useNA="ifany") # NA values show up

# is.na checks for NA values
table(is.na(pdata$age))

# Check for other common missing names
sum(pdata$age==" ") #empty space; common missing value
sum(pdata$age==" ", na.rm = TRUE) #if there is no empty space values = 0

# Check genomic data for NAs
sum(is.na(edata))

# Make the distribution of NA's by genes
gene_na = rowSums(is.na(edata))
table(gene_na) #there are 52580 genes with 0 NA values

# Make the distribution of NA's by samples
sample_na = rowSums(is.na(edata))
table(sample_na)


## ----dimensions----------------------------------------------------------
dim(fdata) #check that all the dimensions match over the 3 tables 
dim(pdata)
dim(edata)

## ----boxplot-------------------------------------------------------------
boxplot(edata[,1]) #First set of values. most of values are down at 0
boxplot(log2(edata[,1]+1)) #transform: log2+1
boxplot(log2(edata+1),col=2,range=0) #plot more than 1 value at once

## ----histograms----------------------------------------------------------
par(mfrow=c(1,2))
hist(log2(edata[,1]+1),col=2) #allmost all the values = 0 
hist(log2(edata[,2]+1),col=2) #same plot for the 2nd sample

## ----densities-----------------------------------------------------------
plot(density(log2(edata[,1]+1)),col=2)
lines(density(log2(edata[,2]+1)),col=3)

## ------------------------------------------------------------------------
qqplot(log2(edata[,1]+1), log2(edata[,2]+1),col=3, pch=19)
abline(c(0,1),lw=2, col=2) #45º line

## ------------------------------------------------------------------------
mm = log2(edata[,1]+1) - log2(edata[,2]+1)
aa = log2(edata[,1]+1) + log2(edata[,2]+1)
plot(aa,mm,col=2, pch=19) #handy for tech replicates
abline(h= 0, col=4, lw=2)

## ------------------------------------------------------------------------
# '''To make figures better, we sometimes have to remove the low features
# '''to see the real distribution of the data
df_edata = as.data.frame(edata)
filt_edata = filter(df_edata,rowMeans(df_edata) > 1) #just keep features with mean >1
dim(filt_edata)
boxplot(as.matrix(log2(filt_edata+1)),col=2) #same boxplot than before
boxplot(as.matrix(log2(filt_edata+1)),col=2, range=0)

## ------------------------------------------------------------------------
aeid = as.character(fdata[,1]) # get Id for the features(genes)
chr = AnnotationDbi::select(org.Hs.eg.db,keys=aeid,keytype="ENSEMBL",columns="CHR")
head(chr) # chr information for those Id

## ------------------------------------------------------------------------
dim(edata) 
dim(chr) #it does not have the same dimensions that the data set
# '''Problem: some of the samples are duplicated label.

# Take non-duplicated chromsomes
chr = chr[!duplicated(chr[,1]),]
dim(chr) #the chr dataset has the same dimension as the expression dataset

# Confirm that the annotation still is in the right order
all(chr[,1] == rownames(edata))

# Select the chromosome Y samples
edatay = dplyr::filter(edata,chr$CHR=="Y") 
edatay = dplyr::filter(df_edata, chr$CHR=="Y")
dim(edatay)#small data frame with genes on chrY

# Males have Y chromsome expression as expected
boxplot(colSums(edatay) ~ pdata$gender, col = 2)
boxplot(colSums(edatay) ~ pdata$gender)
points(colSums(edatay) ~ jitter(as.numeric(pdata$gender)),
       col=as.numeric(pdata$gender),
       pch=19) #overlay points


## ------------------------------------------------------------------------
ematrix = as.matrix(edata)[rowMeans(edata) > 10000,] #filter for genes with high rowMeans
heatmap(ematrix) #multi-variate plots of the data

## ------------------------------------------------------------------------
colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(ematrix,col=colramp)

## ------------------------------------------------------------------------
heatmap(ematrix,col=colramp,Rowv=NA,Colv=NA) #remove the clustering

## ------------------------------------------------------------------------
heatmap.2(ematrix,col=colramp,Rowv=NA,Colv=NA,
          dendrogram="none", scale="row",trace="none")

## ----session_info--------------------------------------------------------
devtools::session_info()
