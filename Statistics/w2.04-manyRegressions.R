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
    library(limma)
    library(edge)
})

## ----load----------------------------------------------------------------
library(devtools)
library(Biobase)
library(limma)
library(edge)

## ----install_packages, eval=FALSE----------------------------------------
#  install.packages(c("devtools"))
#  source("http://www.bioconductor.org/biocLite.R")
#  biocLite(c("Biobase","limma","jdstorey/edge"))

## ------------------------------------------------------------------------
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)
ls()

## ------------------------------------------------------------------------
edata = log2(as.matrix(edata) + 1)    # log transform
edata = edata[rowMeans(edata) > 10, ] # remove low expessed genes
dim(edata)

## ------------------------------------------------------------------------
# You fit the exact same model for every single gene:
mod = model.matrix(~ pdata$strain) # build a model matrix
#fit all those models with lm() function
# you pass the model you've just created and the expression data set that you have
fit = lm.fit(mod,t(edata))  
names(fit)

## ------------------------------------------------------------------------
# Look at the coefficients from the 1st model fit
fit$coefficients[,1]
# You compare that to fitting the linar model to each data set individually
tidy(lm(as.numeric(edata[1, ]) ~ pdata$strain))

# So, "lm" fits one model at a time. "lm.fit" fits many models at a time.
# So the nice thing about lmfit is it's much, much faster. If you wanted to 
# fit every single model with lm, it would take a long time to do that model 
# fitting. 


## ------------------------------------------------------------------------

# Look at the distributions. Since I've hit many regression models I can 
# look at the distribution of coefficients, to see if I see anything interesting or 
# funny about them. 
par(mfrow=c(1,2))
# So first I can look at the distribution of the coefficients for the intercept. 
# So you can see that they tend to be positive here, and that's because we're modeling 
# these sort of account data. 
hist(fit$coefficients[1,],breaks=100,col=2,xlab="Intercept")
# Then, you can model the string coefficients. And so you can see those tend to be close
# to that centered around zero. And so the further from zero they are, the more likely there is to be an association with string. 
hist(fit$coefficients[2,],breaks=100,col=2,xlab="Strain")
abline(v=0,lwd=3,col=1)

## ------------------------------------------------------------------------
# You can plot the residuals from the 1st model fit or from the 2nd model fit
par(mfrow=c(1,2))
plot(fit$residuals[,1],col=2, pch=19)
plot(fit$residuals[,2],col=2, pch=19)

## ------------------------------------------------------------------------
mod_adj = model.matrix(~ pdata$strain + as.factor(pdata$lane.number))
fit_adj = lm.fit(mod_adj,t(edata))
fit_adj$coefficients[,1]

## ------------------------------------------------------------------------
fit_limma = lmFit(edata,mod_adj)
names(fit_limma)
fit_limma$coefficients[1,]
fit_adj$coefficients[,1]

## ------------------------------------------------------------------------
edge_study = build_study(data=edata,grp=pdata$strain,adj.var=as.factor(pdata$lane.number))
fit_edge = fit_models(edge_study)
summary(fit_edge)
fit_edge@beta.coef[1,]
fit_limma$coefficients[1,]

## ----session_info--------------------------------------------------------
devtools::session_info()

