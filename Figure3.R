rm(list=ls())
library(nlme)
library(lme4)
library(ggplot2)
library('randomForest')
library(pdp)           # for partial dependence plots
library(vip)           # for variable importance plots
library(segmented)
library(chngpt)
library(maptools)
library(maps)

r2.mixed<-function(mF){
  mFX<-model.matrix(mF)
  VarF <- var(as.vector(fixef(mF) %*% t(mFX)))
  VarR<-sum(as.numeric(VarCorr(mF))) 
  VarResid<-attr(VarCorr(mF), "sc")^2
  fR2<-VarF/(VarF + VarR + VarResid)
  rfR2<-(VarF + VarR)/(VarF + VarR + VarResid)
  list(fR2=fR2,rfR2=rfR2)
}

#vif function for evaluating models:

vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

# for figure 3a
setwd("~/Documents/Collab/CUE_revision/Nature_communications")
data_Dacal<-read.csv('DacaletalRespData.csv')

# view data distribution
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
points(data_Dacal$Long_decimal,data_Dacal$Lat_decimal, col="red", pch=6)

# use CFE
data_Dacal0<-data_Dacal[,c('NAM','Temperature',
                                 'MAT','CFE',
                                 'Resp','CLC','pH','SOC')]

data_Dacal0$Resp<-log(data_Dacal0$Resp)

Rmodel1min<-lmer(Resp~ Temperature+MAT+CFE+Resp+CLC+pH+SOC+(1|NAM),data=data_Dacal0)
summary(Rmodel1min)
plot(Rmodel1min)
r2.mixed(Rmodel1min)
sqrt(vif.mer(Rmodel1min))

# use partial dependence plot in random forest and plot figure 3a
data_Dacal1<-data_Dacal[,c('Temperature',
                                 'MAT','CFE',
                                 'Resp','CLC','pH','SOC')]

data_Dacal1$Resp<-log(data_Dacal1$Resp)


set.seed(100)
mean_ran_for<-randomForest(Resp ~., 
                           data = data_Dacal1, ntree=300, importance=TRUE,na.action=na.roughfix, mtry=3)

# Default PDP
pdf(file = paste("Figure3a.pdf",sep=""), width = 4.4, height = 3.6)

mean_ran_for %>%  # the %>% operator is read as "and then"
  partial(pred.var = "MAT") %>%
  autoplot(smooth = FALSE, lwd = 1,xlab = 'MAT (°C)',ylab ='Log (Potential respiration rate)',color = "red") + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 0, hjust = 1,size=14),
        axis.text.y = element_text(size=14),axis.title=element_text(size=16)) + xlim(0,25)

dev.off()


df <- partial(mean_ran_for, pred.var = "MAT")
colnames(df) <- c("MAT","Respiration")
plot(df$MAT,df$Respiration,type="l",lwd=2,col="red")


#=================================================================
#Then the "changpt" package in R was used to fit segmented regression and estimate threshold.
lm <- lm(Respiration~MAT,df)
os <-selgmented(lm, Kmax=2, type="score")

summary(os)

intercept(os)
slope(os)

pdf(file = paste("Figure S10a.pdf",sep=""), width = 4.4, height = 3.6)

plot(os, conf.level=0.95, shade=TRUE,xlab='MAT (°C)',ylab='Log (Potential respiration rate)',
     ylim=c(1.5,2.3),xlim=c(-3,25),col = c("red"))

points(df$MAT,df$Respiration,pch=16,cex=0.5)  ## add points

dev.off()

# use data to plot figure 3b
dataBradford<-read.csv('BradfordRespData.csv')

# use original analysis by Bradford 
#Note also that I ran MAT vs. MAT2 as measures of mean annual temperature. Both had the same qualitative effects, 
#but MAT2 was more informative and was collected in the same manner for all locations 
# (see paper; MAT was taken from Crowther et al. 2014 Global Change Biol but mixes broader and site-specific climate data; see below).
#Running the same full model with the outliers identified later in this script removed.
Rmodel1min<-lmer(LogResp~AssTemp+Gluc+Oxal+MAT2+TotPLFA+Clay+Hion+(1|SITE:COVER:YEAR),data=dataBradford[-c(110,370,494),])
summary(Rmodel1min)
plot(Rmodel1min)
r2.mixed(Rmodel1min)
sqrt(vif.mer(Rmodel1min))

# Following the original study by Bradford et al, delete the extreme data points 
data=dataBradford[-c(110,370,494),]

# view data distribution
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
points(data$LongE,data$LatN, col="red", pch=6)

# use partial dependence plot in random forest and plot figure 3b
data_ran<-data[,c('LogResp','AssTemp','Gluc','Oxal',
                  'MAT2','TotPLFA','Clay','Hion')]

data_ran$Gluc<-factor(data_ran$Gluc)
data_ran$Oxal<-factor(data_ran$Oxal)

# use random fores

set.seed(100)

mean_ran_for<-randomForest(LogResp ~., 
                           data = data_ran, ntree=300, importance=TRUE,na.action=na.roughfix, mtry=3)

mean_ran_for


# Default PDP
pdf(file = paste("Figure3b.pdf",sep=""), width = 4.4, height = 3.6)

mean_ran_for %>%  # the %>% operator is read as "and then"
  partial(pred.var = "MAT2") %>%
  autoplot(smooth = FALSE, lwd = 1,xlab = 'MAT (°C)',ylab ='Log (Potential respiration rate)',color = "red") + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 0, hjust = 1,size=14),
        axis.text.y = element_text(size=14),axis.title=element_text(size=16)) + xlim(0,25)

dev.off()


df <- partial(mean_ran_for, pred.var = "MAT2")
colnames(df) <- c("MAT","Respiration")
plot(df$MAT,df$Respiration,type="l",lwd=2,col="red")

lm <- lm(Respiration~MAT,df)
os <-selgmented(lm, Kmax=2, type="score")

summary(os)

intercept(os)
slope(os)

pdf(file = paste("Figure S10b.pdf",sep=""), width = 4.4, height = 3.6)

plot(os, conf.level=0.95, shade=TRUE,xlab='MAT (°C)',ylab='Log (Potential respiration rate)',
     ylim=c(-0.7,0),xlim=c(-3,25),col = c("red"))

points(df$MAT,df$Respiration,pch=16,cex=0.5)  ## add points

dev.off()




