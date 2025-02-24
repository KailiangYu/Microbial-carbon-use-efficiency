rm(list=ls())
library(ggplot2)
library(ggpmisc)
library(randomForest)
library(rfPermute)
library(segmented)
library(pdp)           # for partial dependence plots
library(vip)           # for variable importance plots
library(chngpt)

# figure 2a
setwd("~/Documents/Collab/CUE_revision/Nature_communications")
CUE<-read.csv('data_CUE_original_detrended_env_submission.csv')

# see results without 18O-H2O method 
CUE<-CUE[-c(which(CUE$substrate=='H2O')),]

# for figure 2a
set.seed(1)
CUE1<-na.omit(subset(CUE,select=c(CUE,AI_value,MAT_value,temperature,fb_value,CN_value,SOC_value,substrate)))

CUE1$substrate<-factor(CUE1$substrate)
# scale environmental variabels 
CUE1[,c(2:7)]<-scale(CUE1[,c(2:7)])

mean_ran_for<-randomForest(CUE ~., 
                           data = CUE1, ntree=300, importance=TRUE,na.action=na.roughfix, mtry=3)

# use random forest
impToPlot<-list()
impToPlot_purity<-list()
for (i in 1:100) {
  set.seed(i)
  mean_ran_for<-randomForest(CUE ~., 
                             data = CUE1, ntree=300, importance=TRUE,na.action=na.roughfix, mtry=3)
  impToPlot[[i]]<-unname(randomForest::importance(mean_ran_for)[,1])
}
str(impToPlot)
mat_mean_MSE<-matrix(unlist(impToPlot), ncol = length(CUE1)-1, byrow = TRUE)
mean_clo_MSE<-colMeans(mat_mean_MSE)

mf <- function(x){sd(mat_mean_MSE[,x])}
mean_st_clo_MSE<-sapply(1:length(mat_mean_MSE[1,]),mf) 

Type_gg<-c('Aridity index','MAT','Incubation temperature','F:B ratio','Soil C:N','SOC','Substrate')

# get data.frame
Type<-rep('%IncMSE',length(Type_gg))
ran_gg1<-data.frame(mean_clo_MSE,mean_st_clo_MSE,Type_gg,Type)
colnames(ran_gg1)<-c('mean_clo','mean_st_clo','Type_gg','Type')

Bomperatest<-factor(ran_gg1$Type_gg)
ran_gg1$Type_gg<-factor(Bomperatest,levels(Bomperatest)[rev(c(7,2,3,4,1,5,6))])

pdf(file = paste("Figure 2a.pdf",sep=""), width = 4.4, height = 3.6)

p <- ggplot(ran_gg1, aes(x=factor(Type_gg), y=mean_clo)) + 
  geom_point(size=3, shape=21,position=position_dodge(.6),color=(c('#D95F02','#D95F02','#D95F02','#D95F02','#D95F02','#D95F02','#D95F02'))) +
  geom_errorbar(aes(ymin=(mean_clo)-(mean_st_clo), ymax=(mean_clo)+(mean_st_clo)), width=.2,
                position=position_dodge(.6),color=(c('#D95F02','#D95F02','#D95F02','#D95F02','#D95F02','#D95F02','#D95F02')))


p + scale_x_discrete(labels=)+ labs(x=c(),y="%IncMSE") + coord_flip() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=14)) + theme(legend.title=element_blank(),legend.key=element_rect(fill='white')) + theme(legend.position=c(0.13,0.88))    

dev.off()

names(CUE1)

# use rfPermute to estimate significance of varibales
P_cue<-rfPermute(CUE ~ ., data = CUE1, na.action = na.omit, ntree = 100, num.rep = 50)

plotImportance(P_cue, scale = TRUE) # Significant (p <= 0.05) predictors are in red


# plot for figure 2b
CUE.value<-c(CUE$CUE,CUE$CUE,CUE$CUE,CUE$CUE,CUE$CUE)
value<-c(scale(CUE$fb_value),scale(CUE$MAT_value),scale(CUE$AI_value),scale(CUE$temperature),scale(CUE$CN_value))
Vari<-c(rep('F:B ratio',length(CUE$fb_value)),rep('MAT',length(CUE$MAT_value)),rep('Aridity index',length(CUE$AI_value)),
        rep('Incubation temperature',length(CUE$AI_value)),rep('Soil C:N',length(CUE$AI_value)))

data_all<-as.data.frame(cbind(CUE.value,value,Vari))

data_all$CUE.value<-as.numeric(data_all$CUE.value)
data_all$value<-as.numeric(data_all$value)

pdf(file = paste("Figure 2b.pdf",sep=""), width = 4.4, height = 3.6)

my.formula <- y ~ x
ggplot(data_all, aes(x=value, y=CUE.value, color=Vari, shape=Vari))+ labs(x=c('Scale (Variable)'),y="CUE") +
  geom_point(aes(size=Vari)) + 
  geom_smooth(method=lm, aes(fill=Vari), formula = my.formula)+
#  stat_poly_eq(formula = my.formula, 
#              aes(label = paste(..rr.label..,sep = "~~~")),
#              parse = TRUE,size = 4) +   
  scale_shape_manual(values=c(16,3,7,10,8))+ 
  scale_color_manual(values=c('darkred','blue','black','red','purple'))+
  scale_fill_manual(values = c('darkred','blue','black','red','purple'))+
  scale_size_manual(values=c(0.5,0.5,0.5,0.5,0.5))+
  theme_classic() + theme(legend.position='none') + theme(legend.background = element_rect(color = NA)) +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) + theme(legend.title=element_blank(),legend.key=element_rect(fill='white')) + xlim(-3,3) + ylim(0.3,1.1)

dev.off()

# note that R2 and legend for figure 2b was manually added


#####################################################################################################
# plot for figure 2c
set.seed(1)
CUE$substrate<-factor(CUE$substrate)
# use partial independence plot and see non-linear effects 
CUE1_pdp<-na.omit(subset(CUE,select=c(CUE,AI_value,MAT_value,temperature,CN_value,fb_value,SOC_value,substrate)))

mean_ran_for<-randomForest(CUE ~., 
                           data = CUE1_pdp, ntree=200, importance=TRUE,na.action=na.roughfix, mtry=3)

# Default PDP
pdf(file = paste("Figure 2c.pdf",sep=""), width = 4.4, height = 3.6)

mean_ran_for %>%  # the %>% operator is read as "and then"
  partial(pred.var = "MAT_value") %>%
  autoplot(smooth = FALSE, lwd = 1,xlab = 'MAT (°C)',ylab ='CUE',color = "red") + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 0, hjust = 1,size=14),
        axis.text.y = element_text(size=14),axis.title=element_text(size=14)) +
  ylim(0.45, 0.6) 

dev.off()

############################################################
############################################################
############################################################
# use segmented methods to quantitatively determine the threshold points
df <- partial(mean_ran_for, pred.var = "MAT_value")
colnames(df) <- c("MAT","CUE")
plot(df$MAT,df$CUE,type="l",lwd=2,col="red")


lm <- lm(CUE~MAT,df)
os <-selgmented(lm, Kmax=2, type="score")

summary(os)

pdf(file = paste("Figure S4a.pdf",sep=""), width = 4.4, height = 3.6)

plot(os, conf.level=0.95, shade=TRUE,xlab='MAT (°C)',ylab='CUE',
     ylim=c(0.4,0.65),xlim=c(-3,20),col = c("red"))

points(df$MAT,df$CUE,pch=16,cex=0.5)  ## add points

dev.off()

# get equations for figure 4
intercept(os)
slope(os)


#  CUE <-  0.012107 * MAT + 0.50065       for MAT < 1.190
#  CUE <-  0.001456 * MAT + 0.51333       for 1.190 < MAT < 16.264
#  CUE <-  0.012842 * MAT + 0.32814       for MAT > 16.264

data_CUE_MAT <- seq(min(CUE$MAT_value), max(CUE$MAT_value), by = 0.1)
data_CUE_CUE <- rep(1,length(data_CUE_MAT))
data_CUE <- data.frame(data_CUE_MAT, data_CUE_CUE); colnames(data_CUE) <- c('MAT', 'CUE')

for (i in 1:length(data_CUE_MAT)){
  if (data_CUE$MAT[i] < 1.19){
    data_CUE$CUE[i] <- 0.012107 * data_CUE$MAT[i] + 0.50065
  } else if (data_CUE$MAT[i]>=1.190 & data_CUE$MAT[i]<=16.264){
    data_CUE$CUE[i] <- 0.001456 * data_CUE$MAT[i] + 0.51333 
  } else {
    data_CUE$CUE[i] <- 0.012842 * data_CUE$MAT[i] + 0.32814
  }
}

plot(data_CUE$MAT, data_CUE$CUE,ylim = c(0.4,0.65))





