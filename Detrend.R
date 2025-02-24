
library(rgdal)
library(MCMCglmm)
require(lme4)
require(MuMIn)
library(nlme)

set.seed(100)

setwd("~/Documents/Collab/CUE_revision/Nature_communications")

d<-read.csv('data_original_CUE_env_submission.csv')

#Fit the model using MCMCglmm
CUE_glmm <-MCMCglmm(CUE ~ fb_value + temperature, random = ~substrate, data = d, pr = T)
CUE_glmm_res<-summary.MCMCglmm(CUE_glmm)

#check residuals and overall fit - they look good.----
fitted <- predict(CUE_glmm)
resid <- (d$CUE) - fitted
plot(resid ~ fitted, bty = 'l')

#detrend influence of random effects of substrate.----
pars <- colMeans(CUE_glmm$Sol) #grab parameter means
plot.dat <- fastDummies::dummy_cols(d, select_columns = c('substrate')) #get dummy columns for random effect variable.
plot.dat <- as.data.frame(plot.dat)
plot.dat <- plot.dat[,grep('substrate_',colnames(plot.dat))]
colnames(plot.dat) <- gsub('substrate_','substrate.',colnames(plot.dat))

pars <- pars[names(pars) %in% colnames(plot.dat)]
pars <- pars[order(match(names(pars), colnames(plot.dat)))]

#multiply matrix of random effects by the vector of RE parameters.
adjust <- as.matrix(plot.dat) %*% pars

#detrend CUE for effect of random effects.
d$CUE.detrend <- (d$CUE) - adjust

head(d)

write.csv(d,'data_CUE_original_detrended_env_submission.csv')





