rm(list=ls())
library(dplyr)
library(rgdal)
library(maptools)
library(maps)
library(ggplot2)
library(car)
library(RColorBrewer)
library(rapport)
library(MuMIn)
library(egg)
library(showtext)

#====================================================
#plot Figure 1a, distribution of CUE sites
setwd("~/Documents/Collab/CUE_revision/Nature_communications")

cue <- read.csv("data_CUE_original_detrended_env_submission.csv")

table(cue$substrate)
cue <- cue[,c("latitude","longitude","CUE","substrate")]
cue$Source <- NA
cue$Source[cue$substrate == "H2O"] <- "Water-based"
cue$Source[cue$substrate != "H2O"] <- "C-based"
cue <- cue[cue$Source != "Water-based",]
cue <- cue[,-4]
colnames(cue) <- c("lat","lon","CUE","Source")

map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
points(cue$lon,cue$lat, col="red", pch=6)

pal <- brewer.pal(8,"YlGn")
scaleFUN <- function(x) sprintf("%.1f", x)
#pdf(file = paste("FigureS1.pdf",sep=""), width = 5, height = 3.6)
# Plot 1a - map
map_cue <- ggplot() + 
  ylim(-55, 90) + xlim(-200,180) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "#bababa",
               color = NA,
               size = 0.1) +
  coord_fixed(1.1) +
  geom_point(data = cue,
             aes(x = lon, y = lat, shape=Source, color = CUE),#change 'fill' to 'color' to remove black outline
             size = 1.3) + guides(shape = 'none') +
  scale_color_gradientn(colors = pal, #change 'fill' to 'color' to remove black outline
                        limits = c(0, 1),
                        breaks = seq(from = 0, to = 1, by = 0.5),
                        name = "CUE ",labels = scaleFUN) + 
  theme_minimal() +
  theme(legend.title = element_text(size = 12, hjust = 1.5, vjust = 1),
        legend.text = element_text(size = 12),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.6,"cm"),
        legend.position = c(0.5,0.05),
        legend.direction ="horizontal",
        legend.box = "horizontal",
        panel.grid = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  ggtitle("a") + theme(plot.title = element_text(hjust = 0.05,vjust = -5)) +
  theme(title = element_text(size = 14))

map_cue

map_cue2 <- map_cue +
  annotation_custom(
    ggplotGrob(
      ggplot(cue, aes(x = lat)) + xlab("Latitude") + ylab("Density") +
        geom_density(fill = "#5d8549", color = "white") + scale_x_continuous(labels = function(x) {
          ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N"))
        }, breaks = c(-40, 0, 40, 80),limits = c(-40, 80)) +coord_flip() +theme_classic() + theme(panel.grid.major=element_line(colour=NA),
                                                                                                  panel.background = element_rect(fill = "transparent",colour = NA),
                                                                                                  plot.background = element_rect(fill = "transparent",colour = NA),
                                                                                                  panel.grid.minor = element_blank()) +
        theme(axis.text.y = element_text(size = 8),
              axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
              axis.title = element_text(size = 10))
    ),
    xmin = -200, xmax = -110, ymin = -40, ymax = 50
  )

map_cue2
#dev.off()
#====================================================
#plot Figure 1b, distribution of respiration sites
#====================================================
#Dacal
data_Dacal <- read.csv('DacaletalRespdata.csv')
# data_Dacal <- data_Dacal[data_Dacal$Temperature == 20,]
# data_Dacal <- data_Dacal[data_Dacal$Microsite == 0,]
data_Dacal <- data_Dacal[,c("Lat_decimal","Long_decimal","Resp")]
colnames(data_Dacal) <- c("lat","lon","Resp")
data_Dacal  <- data_Dacal %>% group_by(lat,lon) %>% 
  summarise(Resp=mean(Resp),.groups = 'drop') %>% as.data.frame()
data_Dacal$source <- "Dacal"
hist(data_Dacal$Resp)
# 
# Dacal <- data_Dacal
# coordinates(Dacal) <- ~lon + lat
# proj4string(Dacal) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# plot(Dacal)
# writeOGR(Dacal, ".","geo.Dacal", driver="ESRI Shapefile",overwrite_layer = TRUE)

#Mark
data_Mark  <- read.csv("BradfordRespData.csv")
data_Mark  <- data_Mark[data_Mark$AssSub == "Glucose",]
data_Mark  <- data_Mark[data_Mark$AssTemp == 20,]
#data_Mark  <- data_Mark[data_Mark$Cover == "forest",]
data_Mark  <- data_Mark[,c("Site","LatN","LongE","Resp")]
data_Mark  <- data_Mark %>% group_by(Site,LatN,LongE) %>% 
  summarise(Resp=mean(Resp),.groups = 'drop') %>% as.data.frame()
colnames(data_Mark) <- c("Site","lat","lon","Resp")
data_Mark <- data_Mark[,-1]
data_Mark$source <- "Bradford"
hist(data_Mark$Resp)
# 
# Mark <- data_Mark
# coordinates(Mark) <- ~lon + lat
# proj4string(Mark) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# plot(Mark)
# writeOGR(Mark, ".","geo.Mark", driver="ESRI Shapefile",overwrite_layer = TRUE)

#====================================================
#plot Figure 1b, distribution of respiration sites
resp <- rbind(data_Dacal,data_Mark)
colnames(resp) <- c("lat","lon","Resp","Source")

map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
points(resp$lon,resp$lat, col="red", pch=6)

pal <- brewer.pal(8,"YlGnBu")
#pdf(file = paste("Figure1_b.pdf",sep=""), width = 5, height = 3.6)

map_resp <- ggplot() + 
  ylim(-55, 90) + xlim(-200,180) +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "#bababa",
               color = NA,
               size = 0.1) +
  coord_fixed(1.1) +
  scale_color_gradientn(colors = pal, #change 'fill' to 'color' to remove black outline
                        limits = c(0, 18),
                        name = expression('Respiration ('*mu*g~'C-C'*O[2]~'dry-weight'~soil^-1~h^-1*')')) + 
  geom_point(data = resp,
             aes(x = lon, y = lat, shape=Source, color = Resp),#change 'fill' to 'color' to remove black outline
             size = 1.3) + theme_minimal() + 
  guides(shape = guide_legend(title = NULL,order = 1),color = guide_colorbar(title.position = "top",order = 2,title.theme = element_text(face = 'bold'),barwidth = 15,barheight = 1)) +
  theme(legend.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 12),
        #legend.key.height = unit(0.3, "cm"),
        #legend.key.width = unit(0.6,"cm"),
        legend.position = c(0.4,0.01),
        legend.direction ="horizontal",
        legend.box = "horizontal",
        legend.margin = margin(t = 20, r = 20, b = 20, l = 50),
        panel.grid = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  ggtitle("b") + theme(plot.title = element_text(hjust = 0.05,vjust = -5)) +
  theme(title = element_text(size = 14))
map_resp
map_resp2 <- map_resp +
  annotation_custom(
    ggplotGrob(
      ggplot(resp, aes(x = lat)) + xlab("Latitude") + ylab("Density") +
        geom_density(fill = "#5d8549", color = "white") + scale_x_continuous(labels = function(x) {
          ifelse(x < 0, paste0(abs(x), "°S"), paste0(x, "°N"))
        }, breaks = c(-40, 0, 40, 80),limits = c(-40, 80))+coord_flip() + theme_classic() + theme(panel.grid.major=element_line(colour=NA),
                                                                                                  panel.background = element_rect(fill = "transparent",colour = NA),
                                                                                                  plot.background = element_rect(fill = "transparent",colour = NA),
                                                                                                  panel.grid.minor = element_blank()) +
        theme(axis.text.y = element_text(size = 8),
              axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
              axis.title = element_text(size = 10))
    ),
    xmin = -200, xmax = -110, ymin = -40, ymax = 50
  )



map_resp2
#dev.off()
#=================================
showtext_auto()
pdf(file = paste("Figure1.pdf"), width = 7.8, height = 8, onefile=FALSE)
ggarrange(map_cue2,map_resp2,ncol = 1)
dev.off()

