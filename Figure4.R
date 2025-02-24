rm(list = ls())
library(dplyr)
library(zoo)
library(ggplot2)
library(egg)
library(randomForest)

#####################################################################################################
data_CUE_MAT <- seq(-5, 25, by = 0.05)
data_CUE_CUE <- rep(1,length(data_CUE_MAT))
CUE_MAT <- data.frame(data_CUE_MAT, data_CUE_CUE); colnames(CUE_MAT) <- c('MAT', 'CUE')

for (i in 1:length(data_CUE_CUE)){
  if (CUE_MAT$MAT[i] < 1.19){
    CUE_MAT$CUE[i] <- 0.012107 * CUE_MAT$MAT[i] + 0.50065
  } else if (CUE_MAT$MAT[i]>=1.190 & CUE_MAT$MAT[i]<=16.264){
    CUE_MAT$CUE[i] <- 0.001456 * CUE_MAT$MAT[i] + 0.51333 
  } else {
    CUE_MAT$CUE[i] <- 0.012842 * CUE_MAT$MAT[i] + 0.32814
  }
}


#plot(CUE_MAT$MAT,CUE_MAT$CUE)
#============================================================
#microbial model=============================================
#============================================================
#model
R_T   <- c()
CR_T  <- c()
MBC_T <- c()
SOC_T <- c()
DOC_T <- c()
for (tassay in c(12,20,28)){
  
  CR_output  <- c()
  MBC_output <- c()
  SOC_output <- c()
  DOC_output <- c()
  MAT <- CUE_MAT$MAT
  for (j in MAT) {
    
    Is <- 0 # 0.00015                            #SOC input rate  mg C g-1 soil h-1
    ID <- 0 # 0.00001                            #DOC input rate  mg C g-1 soil h-1
    
    rB <- 0.00028                            #MBC turnover rate (same as kB,ref in CON)
    rE <- 0.0000056                          #Enzyme production rate
    rL <- 0.001                              #Enzyme loss rate
    aBS <- 0.5                               #Fraction of dead MBC transferred to SOC 
    
    B <- 0.459217                            #D, S, or B representing DOC, SOC, and MBC pools
    S <- 18.76216
    D <- 0.01318385
    E <- 0.002571615                         #Extracellular enzyme
    
    Ec <- CUE_MAT$CUE[CUE_MAT$MAT == j]                 #CUE
    
    Tref <- 20 + 273.15                     #Reference temperature, degree C.  + 273.15 convert to K
    T <- tassay + 273.15                        ##need to be changed. Temperature -> assay temperature 
    R <- 8.314/1000                          #J mol-1 K-1, the ideal gas constant  /1000 J to kJ
    Vref <- 1                                #SOC reference Vmax
    VUref <- 0.01                            #DOC uptake reference Vmax (similar to Vref in GER)
    Kref <- 250                              #SOC reference Km
    KUref <- 0.26                            #DOC uptake reference Km 
    EaV <- 47                                #SOC Vmax activation energy
    EaVU <- 47                               #Uptake Vmax activation energy
    EaK <- 30                                #Uptake Km activation energy
    EaKU <- 30                               #DOC activation energy
    
    V <- Vref * exp(-EaV/R*(1/T-1/Tref))     #extracellular enzyme Vmax and Km
    K <- Kref * exp(-EaK/R*(1/T-1/Tref))
    VU <- VUref * exp(-EaVU/R*(1/T-1/Tref))  #reference extracellular enzyme Vmax and Km
    KU <- KUref * exp(-EaKU/R*(1/T-1/Tref))
    
    CR0 <- 0
    B0 <- 0
    simutime <- 1:(24*365) # 1 year
    for (i in simutime) {
      FS <- V * E * S / (K + S)            #Decomposition of SOC pool is catalyzed according to Michaelis-Menten kinetics by the enzyme pool
      FB <- rB * B^2                       #Microbial biomass death is modeled as a first-order process with a rate constant rB
      FE <- rE * B                         #Enzyme production is modeled as a constant fraction (rE) of microbial biomass
      FU <- VU * B * D / (KU + D)          #DOC uptake (FU) by microbes
      FL <- rL * E                         #Enzyme death is modeled as a first-order process with a rate constant rL
      
      deltaS <- Is + FB * aBS - FS
      deltaD <- ID + FB * (1-aBS) + FS + FL - FU
      deltaE <- FE - FL
      deltaB <- FU * Ec - FB - FE         #Microbial biomass increases with DOC uptake (FU) times C use efficiency
      #and declines with death (FB) and enzyme production (FE)
      S <- S + deltaS 
      D <- D + deltaD 
      E <- E + deltaE 
      B <- B + deltaB 
      CR  <- FU * (1-Ec)                   #CO2 respiration is the fraction of DOC uptake that is not assimilated into MBC
      CR0 <- CR                    
      B0  <- B
      S0  <- S
      D0  <- D
    }
    CR_output  <- append(CR_output,CR0)
    MBC_output <- append(MBC_output,B0)
    SOC_output <- append(SOC_output,S0)
    DOC_output <- append(DOC_output,D0)
  }
  CR  <- CR_output
  MBC <- MBC_output
  SOC <- SOC_output
  DOC <- DOC_output
  R  <- CR
  R  <- R*1000*44/12
  CR <- CR/MBC           #mass specific
  CR <- CR*1000*44/12    #*44/12 CO2weight,
  R_T   <- cbind(R_T,R)
  CR_T  <- cbind(CR_T,CR)
  MBC_T <- cbind(MBC_T,MBC)
  SOC_T <- cbind(SOC_T,SOC)
  DOC_T <- cbind(DOC_T,DOC)
}

dat1 <- cbind(R_T,CR_T,MBC_T,SOC_T,DOC_T)
dat1 <- cbind(CUE_MAT,dat1)
colnames(dat1)  <- c("MAT","CUE",paste0(rep(c("R","CR","MBC","SOC","DOC"),each=3),rep(c(12,20,28),times=4)))

dat1$k_eff12 <- -log(dat1$SOC12/18.76216) 
dat1$k_eff20 <- -log(dat1$SOC20/18.76216) 
dat1$k_eff28 <- -log(dat1$SOC28/18.76216) 
#============================================================
#first order model===========================================
#============================================================
#=========================================================
R_T  <- c()
CR_T <- c()
MBC_T <- c()
SOC_T <- c()
DOC_T <- c()
for (tassay in c(12,20,28)){
  
  CR_output  <- c()
  MBC_output <- c()
  SOC_output <- c()
  DOC_output <- c()
  MAT <- CUE_MAT$MAT
  for (j in MAT) {
    
    Is <- 0 # 0.00015                            #SOC input rate  mg C g-1 soil h-1
    ID <- 0 # 0.00001                            #DOC input rate  mg C g-1 soil h-1
    
    B <- 0.07142857                         #D, S, or B representing DOC, SOC, and MBC pools mg C g-1 soil
    D <- 0.03747493
    S <- 32.68743
    
    Ec <- CUE_MAT$CUE[CUE_MAT$MAT == j]                 #CUE
    
    Tref <- 20 + 273.15                     #Reference temperature, degree C.  + 273.15 convert to K
    T <- tassay + 273.15                        #Temperature -> assay temperature 
    
    aDS <- Ec                               #Transfer coefficient from the DOC to the SOC pool
    aSD <- Ec                               #Transfer coefficient from the SOC to the DOC pool   
    aB  <- Ec                               #Transfer coefficient from the MBC to the DOC and SOC pools  
    aBS <- 0.5                               #Partition coefficient for dead MBC between the SOC and DOC pools
    #CON does not include an explicit CUE, but the coefficients that specify partitioning of fluxes
    #into CO2 versus C pools are analogous.
    #lower values of transfer coefficients indicating a larger fraction of C respired as CO2.
    u <- 0.0005                              #DOC uptake rate by microbes
    
    R <- 8.314/1000                          #J mol-1 K-1, the ideal gas constant  /1000 J to kJ
    KSref <- 0.000005                        #SOC decay rate
    KDref <- 0.001                           #DOC decay rate
    KBref <- 0.00028                         #MBC turnover rate
    EaS <- 47                                #SOC activation energy  kJ mol-1
    EaD <- 47                                #DOC activation energy
    EaB <- 20                                #MBC activation energy
    KS <- KSref * exp(-EaS/R*(1/T-1/Tref))   #The decay constant
    KD <- KDref * exp(-EaD/R*(1/T-1/Tref))
    KB <- KBref * exp(-EaB/R*(1/T-1/Tref))
    
    CR0 <- 0
    B0 <- 0
    simutime <- 1:(365*24)
    for (i in simutime) {
      FS <- KS * S                        #Decomposition of each pool
      FD <- KD * D
      FB <- KB * B
      deltaS <- Is + aDS * FD + aB * aBS * FB - FS
      deltaD <- ID + aSD * FS + aB * (1-aBS) * FB - u * D - FD
      deltaB <- u * D - FB
      CR <- FS * (1-aSD) + FD * (1-aDS) + FB * (1-aB)
      S <- S + deltaS 
      D <- D + deltaD 
      B <- B + deltaB 
      CR0 <- CR                    
      B0 <- B
      S0 <- S
      D0 <- D
    }
    CR_output <- append(CR_output,CR0)
    MBC_output <- append(MBC_output,B0)
    SOC_output <- append(SOC_output,S0)
    DOC_output <- append(DOC_output,D0)
  }
  
  CR <- CR_output
  MBC <- MBC_output
  SOC <- SOC_output
  DOC <- DOC_output
  
  R<- CR*1000*44/12    #*44/12 CO2weight,
  CR <- CR/MBC           #mass specific
  CR <- CR*1000*44/12    #*44/12 CO2weight,
  R_T   <- cbind(R_T,R)
  CR_T  <- cbind(CR_T,CR)
  MBC_T <- cbind(MBC_T,MBC)
  SOC_T <- cbind(SOC_T,SOC)
  DOC_T <- cbind(DOC_T,DOC)
}

dat2 <- cbind(R_T,CR_T,MBC_T,SOC_T,DOC_T)
dat2 <- cbind(CUE_MAT,dat2)
colnames(dat2)  <- c("MAT","CUE",paste0(rep(c("R","CR","MBC","SOC","DOC"),each=3),rep(c(12,20,28),times=4)))


dat2$k_eff12 <- -log(dat2$SOC12/32.68743) 
dat2$k_eff20 <- -log(dat2$SOC20/32.68743) 
dat2$k_eff28 <- -log(dat2$SOC28/32.68743) 

#====================================================================
#k_eff-R
k_eff20_m <- dat1[,c("MAT","k_eff20")]
k_eff20_m$model <- "Microbial"
colnames(k_eff20_m) <- c("MAT","k_eff","model")
k_eff20_m$k_eff <- (k_eff20_m$k_eff - k_eff20_m$k_eff[k_eff20_m$MAT == 20])/k_eff20_m$k_eff[k_eff20_m$MAT == 20] * 100


k_eff20_f <- dat2[,c("MAT","k_eff20")]
k_eff20_f$model <- "First-order"
colnames(k_eff20_f) <- c("MAT","k_eff","model")
k_eff20_f$k_eff <- (k_eff20_f$k_eff - k_eff20_f$k_eff[k_eff20_f$MAT == 20])/k_eff20_f$k_eff[k_eff20_f$MAT == 20] * 100

#relative change
k_eff <- rbind(k_eff20_m,k_eff20_f)


p <- ggplot(k_eff,aes(x=MAT,y=k_eff)) + 
  geom_point(aes(colour = factor(model)),size = 1.5, shape=1) + #display group points
  geom_line(aes(colour = factor(model)), size=0.8)+
  coord_cartesian(xlim = c(-5,25), ylim = c(-20, 20))+           #axix range 
  scale_y_continuous(breaks = seq(from = -20, to = 20, by = 10)) + #axix interval
  scale_x_continuous(breaks = seq(from = -5, to = 25, by = 5),position = "bottom")
p <- p + scale_color_manual(values=c("#5e3c99","#e66101"),breaks=c("Microbial","First-order"),labels=c("Microbial","First-order"))
p <- p + theme_bw() #delete gray backound
p <- p + theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank())
p <- p + labs(y = "Changes in effective K (%)") + labs(x = "MAT (°C)")
p <- p + annotate(geom = 'text', x=-Inf, y=Inf,hjust=-0.1,vjust=1.5,
                  label=" ", size = 5.5)   #annotation
p <- p + theme(axis.title.x = element_text(size = 18)) + #x axix title, font size and type 
  theme(axis.title.y = element_text(size = 18)) +
  theme(axis.text = element_text(size = 18,colour="black"))
p <- p + theme(legend.title = element_blank())+ #delete legend title
  theme(legend.text = element_text(size = 18)) 
p <- p + ggtitle("a") + theme(plot.title = element_text(hjust = -0.17)) +
  theme(title = element_text(size = 18))
p <- p + guides(colour = guide_legend(ncol = 2)) # two columns, legend
p1 <- p + theme(legend.position = c(0.4,0.92))
# p1 <- p + theme(aspect.ratio=3/4,
#                plot.margin = margin(2, 2, 2, 10))  #  top, right, bottom,left ))
p1
#====================================================================
#SOC-R
SOC20_m <- dat1[,c("MAT","SOC20")]
SOC20_m$model <- "Microbial"
colnames(SOC20_m) <- c("MAT","SOC","model")
SOC20_m$SOC <- (SOC20_m$SOC - SOC20_m$SOC[SOC20_m$MAT == 20])/SOC20_m$SOC[SOC20_m$MAT == 20] * 100

SOC20_f <- dat2[,c("MAT","SOC20")]
SOC20_f$model <- "First-order"
colnames(SOC20_f) <- c("MAT","SOC","model")
SOC20_f$SOC <- (SOC20_f$SOC - SOC20_f$SOC[SOC20_f$MAT == 20])/SOC20_f$SOC[SOC20_f$MAT == 20] * 100

SOC <- rbind(SOC20_m,SOC20_f)

p <- ggplot(SOC,aes(x=MAT,y=SOC)) + 
  geom_point(aes(colour = factor(model)),size = 1.5, shape=1) + #display group points
  geom_line(aes(colour = factor(model)), size=0.8)+
  coord_cartesian(xlim = c(-5, 25), ylim = c(-0.6, 0.8))+           #axix range 
  scale_y_continuous(breaks = seq(from = -0.6, to = 0.8, by = 0.3)) + #axix interval
  scale_x_continuous(breaks = seq(from = -5, to = 25, by = 5),position = "bottom")
p <- p + scale_color_manual(values=c("#5e3c99","#e66101"),breaks=c("Microbial","First-order"),labels=c("Microbial","First-order"))
p <- p + theme_bw() #delete gray backound
p <- p + theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank())
p <- p + labs(y = "Change in SOC (%)") + labs(x = "MAT (°C)")
p <- p + annotate(geom = 'text', x=-Inf, y=Inf,hjust=-0.1,vjust=1.5,
                  label=" ", size = 5.5)   #annotation
p <- p + theme(axis.title.x = element_text(size = 18)) + #x axix title, font size and type 
  theme(axis.title.y = element_text(size = 18)) +
  theme(axis.text = element_text(size = 18,colour="black"))
p <- p + theme(legend.title = element_blank())+ #delete legend title
  theme(legend.text = element_text(size = 18)) 
p <- p + ggtitle("b") + theme(plot.title = element_text(hjust = -0.17)) +
  theme(title = element_text(size = 18))
p2 <- p + theme(legend.position = "none")       #legend position
p2
# p2 <- p + theme(aspect.ratio=3/4,
#                plot.margin = margin(2, 2, 2, 2))  #  top, right, bottom,left ))
# p2
#============================================================

setwd("~/Documents/Collab/CUE_revision/Nature_communications")

pdf("Figure4.pdf",width=10,height=5,onefile=F)

ggarrange(p1, p2, ncol  = 2)

dev.off()

