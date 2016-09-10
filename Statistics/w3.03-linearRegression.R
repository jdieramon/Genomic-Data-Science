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
    library(broom)
})

## ----load----------------------------------------------------------------
library(devtools)
library(Biobase)
library(broom)

## ----install_packages, eval=FALSE----------------------------------------
#  install.packages(c("devtools","broom"))
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
edata = as.matrix(edata)
# Fit the linear model : relationship between expression and age
lm1 = lm(edata[1,] ~ pdata$age) # edata[1,] is the first gene
tidy(lm1)

## ------------------------------------------------------------------------
plot(pdata$age,edata[1,], col=1, pch=19)
abline(lm1$coeff[1],lm1$coeff[2], col=2,lwd=3) # add fitted line
#abline(lm1, col=2) #this is the same 

## ------------------------------------------------------------------------
# What if you have a factor covariate? 
pdata$gender
table(pdata$gender)

## ------------------------------------------------------------------------
# How can we quantify that variability?
boxplot(edata[1,] ~ pdata$gender)
points(edata[1,] ~ jitter(as.numeric(pdata$gender)),
       col=as.numeric(pdata$gender), pch=19)

## ------------------------------------------------------------------------
dummy_m = pdata$gender=="M"
dummy_m *1

dummy_f = pdata$gender=="F"
dummy_f *1

## ------------------------------------------------------------------------
#lm does a dummy variable for you
lm2 = lm(edata[1,] ~ pdata$gender)
tidy(lm2)

## ------------------------------------------------------------------------
mod2 = model.matrix(~pdata$gender)
mod2

## ------------------------------------------------------------------------
# You can do the same even if you have more categories
table(pdata$tissue.type)
pdata$tissue.type == "adipose"
pdata$tissue.type == "adrenal"

## ------------------------------------------------------------------------
tidy(lm(edata[1,] ~ pdata$tissue.type ))

## ------------------------------------------------------------------------
# You can adjust for variables (expanding the model formula)
lm3 = lm(edata[1,] ~ pdata$age + pdata$gender)
tidy(lm3)

## ------------------------------------------------------------------------
# Interaction models (* command in the formula)
lm4 = lm(edata[1,] ~ pdata$age*pdata$gender)
tidy(lm4)

## ------------------------------------------------------------------------
lm4 = lm(edata[6,] ~ pdata$age)
plot(pdata$age,edata[6,],col=2, pch=19)
abline(lm4,col=1,lwd=3)
#In this example the outlier does not affect the line
## ------------------------------------------------------------------------
index = 1:19
lm5 = lm(edata[6,] ~ index)
plot(index,edata[6,],col=2, pch=19)
abline(lm5,col=1,lwd=3)

lm6 = lm(edata[6,-19] ~ index[-19])
abline(lm6,col=3,lwd=3)

legend(5,1000,c("With outlier","Without outlier"),col=c(1,3),lwd=3)
#You have to be careful when you have an outlier

## ------------------------------------------------------------------------
par(mfrow=c(1,2))
hist(lm6$residuals,col=2)
hist(lm5$residuals,col=3)

## ------------------------------------------------------------------------
gene1 = log2(edata[1,]+1)
lm7 = lm(gene1 ~ index)
hist(lm7$residuals,col=4)

## ------------------------------------------------------------------------
lm8 = lm(gene1 ~ pdata$tissue.type + pdata$age)
tidy(lm8)

## ------------------------------------------------------------------------
colramp = colorRampPalette(1:4)(17)
lm9 = lm(edata[2,] ~ pdata$age)
plot(lm9$residuals,col=colramp[as.numeric(pdata$tissue.type)])

## ----session_info--------------------------------------------------------
devtools::session_info()
