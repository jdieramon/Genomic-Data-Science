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

## ---load data -------------------------------------------------------
#We will use this expression set that combines two studies [Transcriptome genetics using second generation sequencing in a 
#Caucasian population.](http://www.ncbi.nlm.nih.gov/pubmed?term=20220756%5Buid%5D) and [Understanding mechanisms underlying 
#human gene expression variation with RNA sequencing.](http://www.ncbi.nlm.nih.gov/pubmed?term=20220758). 
#These studies are different populations but we counted the same genes for both. Then we'll explore the differences. 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
rm(con,montpick.eset, mp)
ls()

## ------------------------------------------------------------------------
edata = edata[rowMeans(edata) > 100, ]   #reduce number of rows
edata = log2(edata + 1)                  #log2 transform : scale , easier to work with
edata_centered = edata - rowMeans(edata) #center the data
svd1 = svd(edata_centered)
names(svd1)

## ------------------------------------------------------------------------
#The percent of variance explained is given by $\frac{d_{ii}}{\sum_{j}d_{jj}^2}$
plot(svd1$d,ylab="Singular value",col=2)
plot(svd1$d^2/sum(svd1$d^2),ylab="Percent Variance Explained",col=2)
#The first singular value explains > 50% of the variance. So, it's a highly explanatory variable.  
x = sum(svd1$d[1:2]^2/sum(svd1$d^2))
#The two first dimensions explain `r x` of the total variance. 

## ------------------------------------------------------------------------
par(mfrow=c(1,2), mar=c(4.5,4.1,1,1))
plot(svd1$v[,1],col=2,ylab="1st PC")
plot(svd1$v[,2],col=2,ylab="2nd PC")

## ---------PC1 vs. PC2 ----------------------------------------------------
# A very common plot is to plot PC1 versus PC2 to see if you can see any "clusters" or "groups"
par(mfrow=c(1,1))
plot(svd1$v[,1],svd1$v[,2],col=2,ylab="2nd PC",xlab="1st PC", pch = 19)

## ------------------------------------------------------------------------
# One thing you can do is color them by different variables to see if clusters stand out. 
plot(svd1$v[,1],svd1$v[,2],ylab="2nd PC",
     xlab="1st PC",col=as.numeric(pdata$study), pch=19)


# '''If you look in the PC1 axis, the two studies have very different 
# ''' values of the PC. So it seems that one of the big sources of signals in
# ''' the data set is which study the two data sets come from

## ------------------------------------------------------------------------
# Another common plot is to make boxplots comparing the PC for different levels 
# of known covariates (don't forget to show the actual data!).
boxplot(svd1$v[,1] ~ pdata$study,border=c(1,2))
points(svd1$v[,1] ~ jitter(as.numeric(pdata$study)),col=as.numeric(pdata$study))

# '''there's a big difference in the PC1 between the Montgomery and the 
# ''' Pickrell studies.

## --------PCs versus SVs -------------------------------------------------
### What we have been plotting is not exactly the PC. 
### They are not exactly the same (they are not scaled in the same way)
pc1 = prcomp(edata)
plot(pc1$rotation[,1],svd1$v[,1])


## ------------------------------------------------------------------------
# To get the actual PCs you have to subtract the column means rather than the 
# row means when normalizing. 
edata_centered2 = t(t(edata) - colMeans(edata)) #data centered by columns
svd2 = svd(edata_centered2)
plot(pc1$rotation[,1],svd2$v[,1],col=2)

# Despite this, it is most common for people to perform row-centering and then 
# plot the singular vectors (sometimes labeling them PCs like I have done in this document)

## ------------------------------------------------------------------------
# Outliers can drive these decompositions
# What happens if we introduce a single outlying gene
edata_outlier = edata_centered
edata_outlier[6,] = edata_centered[6,] * 10000
svd3 = svd(edata_outlier)
plot(svd1$v[,1],svd3$v[,1],xlab="Without outlier",ylab="With outlier")

# It turns out the new top singular vector is perfectly correlated with the outlying gene

## ------------------------------------------------------------------------
plot(svd3$v[,1],edata_outlier[6,],col=4)

## ----Further resources --------------------------------------------------
#There are a large number of resources available about PCA and SVD but the 
#lecture notes from [Advanced Statistics for the Life Sciences](http://genomicsclass.github.io/book/)
#are the best set of lecture notes focused on genomics currently available. 
                                                                                                                                                                                           
                                                                                              
## ----session_info--------------------------------------------------------
devtools::session_info()
