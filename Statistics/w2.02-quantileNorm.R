# quantile normalization

Normalization : common distribution across samples

We want to normalize the samples, so they have common distributions. 
We´ll do this when we think the distributions are driven by technical variables. 

## ----global_palette, results = 'asis'------------------------------------
rm(list=ls())
tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

## ----global_options,warning=FALSE,message=FALSE--------------------------
## see ch. 10 Hooks of Xie's knitr book
# library(knitr)
# knit_hooks$set(setPch = function(before, options, envir) {
#   if(before) par(pch = 19)
# })
# opts_chunk$set(setPch = TRUE)

## ----global_plot,warning=FALSE, message=FALSE----------------------------
# knitr::opts_chunk$set(fig.width=5, fig.height=5, size="footnotesize",
#                       warning=FALSE, message=FALSE)
# knitr::knit_hooks$set(small.mar = function(before, options, envir) {
#   if (before) graphics::par(mar = c(5,5,1.5,1))
# })

## ----load_hidden, echo=FALSE, results="hide", warning=FALSE--------------
# suppressPackageStartupMessages({
#   library(devtools)
#   library(Biobase)
#   library(preprocessCore)
# })

## ----load----------------------------------------------------------------
library(devtools)
library(Biobase)
library(preprocessCore)

## ----install_packages, Bioc v3.8 ----------------------------------------
# install.packages(c("BiocManager"))
# BiocManager::install(c("preprocessCore"))



## ------------------------------------------------------------------------
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
rm(con,montpick.eset, mp)
ls()

## ----transform data and remove outlying values --------------------------
edata = log2(edata + 1)
dim(edata)
edata = edata[rowMeans(edata) > 3, ]
dim(edata)
# Show the distributions of each of these samples
colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:20){lines(density(edata[,i]),lwd=3,col=colramp[i])}
# these differences are probably due to technology and not biology

## --quantilenormalization-------------------------------------------------
norm_edata = normalize.quantiles(as.matrix(edata))
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.20))
for(i in 2:20){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}

#There is only a  bit of variability on the left low end because The
#quantiles for the very low and high values are difficultu to match up. 

## ------------------------------------------------------------------------
# We have removes differences in variability but not the gene by gene 
# variability.  So, we can study it. 
plot(norm_edata[1,],col=as.numeric(pdata$study), pch=19)
# We still have differences between the two studies. 

## ------------------------------------------------------------------------
#To see this more clearly we can do decomposition
svd1 = svd(norm_edata - rowMeans(norm_edata))
plot(svd1$v[,1],svd1$v[,2],xlab="PC1",ylab="PC2",
     col=as.numeric(pdata$study), pch=19)

# important : even though we've normalized out the total distribution, 
# we can still have artifacts like batch effects or other types of artifacts 
# in the data set. 

## ----session_info--------------------------------------------------------
devtools::session_info()
