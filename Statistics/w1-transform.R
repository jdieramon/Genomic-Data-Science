# Data Transforms
## transforms and filtering low-counts observations  

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
})

## ----load----------------------------------------------------------------
library(devtools)
library(Biobase)

## ----install_packages, eval=FALSE----------------------------------------
#  install.packages(c("devtools"))
#  source("http://www.bioconductor.org/biocLite.R")
#  biocLite(c("Biobase"))

## ------------------------------------------------------------------------
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
fdata = fData(bm)
ls()

## ------------------------------------------------------------------------
hist(rnorm(1000),col=2) # typical plot of normal data (simulation)

## ------------------------------------------------------------------------
# '''histogram for the 1st sample in the expression dataset: almost all values are = 0 
hist(edata[,1],col=2,breaks=100) 
hist(edata[,1],col=2,breaks=1000) #increse the bins 

## ------------------------------------------------------------------------
# '''transformation I: log scale(put things on a multiplicative scale)
hist(log(edata[,1]),col=2,breaks=100)

## ------------------------------------------------------------------------
min(log(edata)) #log0 is undefined, so the plot looks better
quantile(log(edata[,1])) #it can't calculate it (becasue mostly we have 0 in the dataset)

## ------------------------------------------------------------------------
# '''You can add a small number to avoid that problem
min(log(edata[,1] + 1))
hist(log(edata[,1] + 1),breaks=100,col=2)
hist(log(edata[,1] + 1),breaks=50,col=2)

## ------------------------------------------------------------------------
# '''transformation II: log2 scale(put things on a multiplicative scale)
hist(log2(edata[,1] + 1),breaks=100,col=2)

## ------------------------------------------------------------------------
# '''zoom in (=ignore the 0s)
hist(log2(edata[,1] + 1),breaks=100,col=2,xlim=c(1,15),ylim=c(0,400))

## ------------------------------------------------------------------------
# '''Calculate for each row(gene) how many values = 0
hist(rowSums(edata==0),col=2) #~40000 genes have basically 0(almost every value is 0)

## ------------------------------------------------------------------------
# '''remove the low value genes : taking the mean for each row
low_genes = rowMeans(edata) < 5 #you can pick a different number
table(low_genes)
filt_edata = filter(as.data.frame(edata),!low_genes) #keep the ones that are not low genes 
dim(filt_edata)

# ''' The average isn't always necessarily great for count data because you often have these really high values.
summary(edata)
# '''for a lot of genes the mean is high just because the max. value is high. 
# '''But the median is not affected for this, so the median is a bit more robust

# '''remove the low value genes : taking the median for each row
low_genes2 = rowMedians(as.matrix(edata)) < 5
table(low_genes2,low_genes)
filt_edata2 = filter(as.data.frame(edata),!low_genes2)
dim(filt_edata2)

## ------------------------------------------------------------------------
hist(log2(filt_edata[,1] + 1),col=2) #plot after removing the low value genes

## ----session_info--------------------------------------------------------
devtools::session_info()
