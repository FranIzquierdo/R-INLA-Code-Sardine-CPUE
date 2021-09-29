#*********************************
# Paper Sardine CPUE             #
# INLA best model and prediction #
#*********************************
# Francisco Izquierdo        #
# francisco.izqtar@gmail.com #
#*****************************

## https://becarioprecario.bitbucket.io/spde-gitbook/
## To see the index press: Ctrl + shift + O

# *Best model* ----------------------------------------------------------------

rm(list=ls()) # Clean workspace
library(INLA)
data <- read.table("./data/working_data/5_final_dataobs_esc_GPS.txt" , 
                   header=TRUE, dec=".")
data$year<-as.factor(data$year)

## ---- mesh
coords4<-as.matrix(data[,1:2])
length(coords4[,1])
bound=inla.nonconvex.hull(as.matrix(data[,1:2]), convex=-0.07, eps=0.05, 
                          resolution=40)
mesh.s <- inla.mesh.2d(loc=coords4, boundary=bound, max.edge=c(0.05, 0.6), 
                       offset=c(0.1,0.6), cutoff=0.08, min.angle = 0.05)
plot(mesh.s,asp=1)
points(coords4, pch=19, cex=0.7)
mesh.s$n

## Prediction data -------------------------------------------------------------

# Month and year
month_pred<-c(rep(1:12, times=1, each=mesh.s$n))
length(month_pred)
table(month_pred)
year_pred<-c(rep(1, times=mesh.s$n*12))
length(year_pred)
# Length (mean length fleet)
unique(inla.group(data$length, n=6, method="cut"))
mean(data$length)
round(mean(data$length))
# so mean vessel is on group of 22 m length
length_pred<-rep(round(mean(data$length)), times=mesh.s$n*12)
length(length_pred)
# Coordinates(mesh coords repeated each month and year)
lat_pred<-rep(mesh.s$loc[,1], times=12)
lon_pred<-rep(mesh.s$loc[,2], times=12)
# Createe datapred
datapred<-cbind(year_pred,month_pred, length_pred,lat_pred, lon_pred)
datapred<-as.data.frame(datapred)
head(datapred)

#datapred<-read.table(file="./data/working_data/5_final_datapred_esc_GPS.txt", dec=".", header=TRUE)

## ---- spde
spde <- inla.spde2.pcmatern(mesh.s, prior.range=c(0.5, 0.05), prior.sigma=c(0.6, 0.05))
(m <- spde$n.spde)

# ----- AMat
Ast <- inla.spde.make.A(mesh = mesh.s, loc = coords4,
                        group = data$month_id,
                        n.group = length(unique(data$month_id)))

Astp <- inla.spde.make.A(mesh = mesh.s,
                         index=rep(1:mesh.s$n, times=12),
                         group=datapred$month_pred,
                         n.group = length(unique(datapred$month_pred)))# mesh es mesh.s, localizaciones son los puntos

dim(Ast)
dim(Astp)

## ---- idx
idx<- inla.spde.make.index("s", n.spde=spde$n.spde, n.group=max(data$month_id))

# ----- Priors
hyper.prec <- list(theta = list(prior="pc.prec", 
                                param = c(1, 0.05))) # allow smaller values prec iid
hyper.prec.bath <- list(theta = list(prior="pc.prec", 
                                     param = c(1, 0.05))) # Good prior for bath in Gamma 
pcrho <- list(prior = 'pccor1', param = c(0.7, 0.7))

# Formula ----------------------------------------------------------------------

## Model eff + ves + st + year + month rw2

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          month=data$month_id)), tag="est")

stk.p <- inla.stack(  data = list(y = NA, link=1),
                      A=list(Astp, 1),
                      effects=list(idx,list(a0=rep(1, length=length(datapred$lat)), 
                                            length=datapred$length_pred, 
                                            month=datapred$month_pred)), 
                      tag="pred")# ID vars not necessary

stk.full <- inla.stack(stk.e, stk.p)

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, 
    control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  f(month, model="rw2", hyper=hyper.prec, cyclic=TRUE)

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.full),
                control.predictor = list(A = inla.stack.A(stk.full),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, 
                num.threads=2)

saveRDS(abu.res, "./5_Best_model_prediction/rds/Best_model_len_ves_st_yeariid_monthrw2.rds")

# Plot RW2-iid --------------------------------------------------------------

shared<-readRDS("./5_Best_model_prediction/rds/Best_model_len_ves_st_yeariid_monthrw2.rds")# M1

# Year iid
library(ggplot2)
suabinm <- shared$summary.random$year$mean
suabin2 <- shared$summary.random$year$`0.025quant`
suabin9 <-shared$summary.random$year$`0.975quant`
suabinID<-shared$summary.random$year$ID
suabin<-data.frame(suabinm, suabin2,suabin9,suabinID)

m1<-ggplot(data = suabin, aes(x = suabinID, y = suabinm))+
  geom_line(aes(x = suabinID, y = suabinm), color="#29AF7FFF", size=0.9)+ 
  geom_ribbon(aes(x = suabinID, ymin = (suabin2), ymax = (suabin9)), 
              alpha = 0.25, fill="gray70", linetype=1)+
  ggtitle(" ")+
  xlab("Year ID")+
  ylab("Year effect ")+
  theme_light() + scale_x_continuous(breaks=c(1,2,3),labels=c("2011","2012","2013"),
                                     limits=c(1, 3)) + 
  theme(plot.title = element_text(hjust=0.5))

m1

# Vessel iid
library(ggplot2)
suabinm <- shared$summary.random$vessel$mean
suabin2 <- shared$summary.random$vessel$`0.025quant`
suabin9 <-shared$summary.random$vessel$`0.975quant`
suabinID<-shared$summary.random$vessel$ID
suabin<-data.frame(suabinm, suabin2,suabin9,suabinID)

m2<-ggplot(data = suabin, aes(x = suabinID, y = suabinm))+
  
  geom_line(aes(x = suabinID, y = suabinm), color="#29AF7FFF", size=0.9)+
  geom_ribbon(aes(x = suabinID, ymin = suabin2, ymax = suabin9), alpha = 0.25,
              fill="gray70")+
  ggtitle(" ")+
  xlab("Vessel ID")+
  ylab("Vessel effect ")+
  theme_light()+ scale_x_continuous(breaks=seq(1,66,4), limits=c(1, 66)) + 
  theme(plot.title = element_text(hjust=0.5))

m2

# Month RW2
library(ggplot2)
suabinm <- shared$summary.random$month$mean
suabin2 <- shared$summary.random$month$`0.025quant`
suabin9 <-shared$summary.random$month$`0.975quant`
suabinID<-shared$summary.random$month$ID
suabin<-data.frame(suabinm, suabin2,suabin9,suabinID)

m3<-ggplot(data = suabin, aes(x = suabinID, y = suabinm))+
  
  geom_line(aes(x = suabinID, y = suabinm), color="#33638DFF", size=0.9)+
  geom_line(aes(suabinID, (suabin2)), color = "grey50", size = 0.1, linetype="dashed") + 
  geom_line(aes(suabinID, (suabin9)), color = "grey50", size = 0.1, linetype="dashed") +
  ggtitle(" ")+
  xlab("Month")+
  ylab("Month effect ")+
  theme_light() + ylim(-0.6,0.6) + theme(plot.title = element_text(hjust=0.5))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),
                     labels=c("1","2","3","4","5","6","7","8","9","10","11","12"), 
                     limits=c(1, 12))+ theme(plot.title = element_text(hjust=0.5))

m3

# Length RW2
library(ggplot2)
suabinm <- shared$summary.random$length$mean
suabin2 <- shared$summary.random$length$`0.025quant`
suabin9 <-shared$summary.random$length$`0.975quant`
suabinID<-shared$summary.random$length$ID
suabin<-data.frame(suabinm, suabin2,suabin9,suabinID)

m4<-ggplot(data = suabin, aes(x = suabinID, y = suabinm))+
  
  geom_line(aes(x = suabinID, y = suabinm), color="#33638DFF", size=0.9)+
  geom_line(aes(suabinID, (suabin2) ), color = "grey50", size = 0.1, linetype="dashed") + 
  geom_line(aes(suabinID, (suabin9)), color = "grey50", size = 0.1, linetype="dashed") +
  ggtitle(" ")+
  xlab("Length (m)")+
  ylab("Length effect ")+
  theme_light() + theme(plot.title = element_text(hjust=0.5))+
  scale_x_continuous(breaks=c(16,18,20,22,24,26),
                     labels=c("16","18","20","22","24","26"), limits=c(16, 26))+ 
  theme(plot.title = element_text(hjust=0.5))

m4

library(patchwork)
(m4+m2)/(m3+m1)
ggsave("5_Best_model_prediction/Final_effects_pr.jpg")

## Plot spatial--------------------------------------------------------------

library(sf)
library(sp)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)

shape_path <- "./data/working_data/2_oceandatafiles/50m_physical_separated_zip/"
ocean_shapefile <- paste(shape_path, "ne_50m_ocean/ne_50m_ocean.shp", sep="")
layer <- ogrListLayers(ocean_shapefile)
ogrInfo(ocean_shapefile, layer=layer)
ocean_poly <- readOGR(ocean_shapefile, layer=layer)# read the shape file
ocean <- ocean_poly
bbx <- readWKT("POLYGON((-9.899 41.86, -6.136 41.86, -6.136 36.9, -9.899 36.9, -9.899 41.86))") 
proj4string(bbx)<-proj4string(ocean)#Le damos CRS de pen_rec (WGS84)
ocean.cut <- gIntersection(ocean, bbx)
pen.cut<-gDifference(bbx, ocean.cut)

library(ggplot2)
library(dplyr)
require(maps)
require(viridis)
library(marmap)
theme_set(
  theme_light()
)

world_map <- map_data("world")
some.eu.countries <- c(
  "Portugal", "Spain"
)
# Retrievethe map data
reg <- map_data("world", region = some.eu.countries)

# create the breaks- and label vectors
ewbrks <- seq(-10.5,-5.5,1)
nsbrks <- seq(37,42,1)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(-x, "°W"), 
                                                   ifelse(x > 0, paste(x, "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "°S"), 
                                                   ifelse(x > 0, paste(x, "°N"),x))))

# get bathymetry data
b = getNOAA.bathy(lon1 = -9.84, lon2 = -8.5, lat1 = 41.82, lat2 = 37.1, 
                  resolution = 1)

# set map limits
lons = c(-10.5, -7)
lats = c(37, 42 )
bf = fortify.bathy(b)

library(splancs)
box<-bbox(ocean.cut)
Xrange<-c(box[1,1], box[1,2])
Yrange<-c(box[2,1], box[2,2])
nxy=c(120,120)

prj <- inla.mesh.projector(mesh.s, xlim=range(Xrange),
                           ylim=range(Yrange), dims=nxy) 

# mean
k<-max(unique(data$month_id))
m<-mesh.s$n
m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         shared$summary.ran$s$mean[1:m + (j - 1) * m])
  return(r) 
})

library(fields)
library(RColorBrewer)
#par(mfrow = c(3, 4), mar = c(0, 1.3, 1.2, 0.8))
(zlm <- range(unlist(m.prj), na.rm = TRUE)) # para definir breaks

library(RColorBrewer)
palette_blues <- colorRampPalette(colors = c("lightgrey", "#009FFD"))(3)
scales::show_col(palette_blues)

for (j in 1:12) {
  
  
  library(raster)
  sp.mean.raster1<-raster(list(x=prj$x,
                               y = prj$y,
                               z = m.prj[[j]]))
  bath<-raster("./data/working_data/raster_env/GPS_envars/Bath.tif")
  res<-resample(bath,sp.mean.raster1)
  aa<-mask(sp.mean.raster1,res)
  aa <- rasterToPoints(aa)
  aa<-data.frame(aa)
  colnames(aa)<-c("Longitude","Latitude","MAP")
  #Now make the map
  pr<- ggplot() +
    geom_tile(data=aa, aes(y=Latitude, x=Longitude, fill=MAP)) +
    scale_fill_gradient2(low="#440154FF",high="white",mid="#2D708EfF",
                         midpoint = -0.5, na.value = "transparent",
                         breaks=c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1),
                         labels=c("-2.5","-2","-1.5","-1","-0.5","0","0.5","1"),
                         limits=c(-2.5,1), name="Mean")+
    geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
                 fill= "grey", color = "darkgrey") +
    coord_map(xlim = lons, ylim = lats)+
    # formatting
    scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
    ylab(" ")+xlab(" ")+
    theme_light() + ggtitle(paste0( j))+ labs(fill = "CPUE") +
    theme(axis.text.x= element_text(size=8), axis.text.y= element_text(size=8),
          plot.title = element_text(hjust = 0.5))
  
ggsave(paste("./5_Best_model_prediction/ggplot_maps/spatial/mean/",j,"pr.jpg"))
  
}

# sd
k<-max(unique(data$month_id))
m<-mesh.s$n
m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         shared$summary.ran$s$sd[1:m + (j - 1) * m])
  return(r) 
})


library(fields)
library(RColorBrewer)
#par(mfrow = c(3, 4), mar = c(0, 1.3, 1.2, 0.8))
(zlm <- range(unlist(m.prj), na.rm = TRUE)) # para definir breaks

for (j in 1:12) {
  
  
  library(raster)
  sp.mean.raster1<-raster(list(x=prj$x,
                               y = prj$y,
                               z = m.prj[[j]]))
  bath<-raster("./data/working_data/raster_env/GPS_envars/Bath.tif")
  res<-resample(bath,sp.mean.raster1)
  aa<-mask(sp.mean.raster1,res)
  aa <- rasterToPoints(aa)
  aa<-data.frame(aa)
  colnames(aa)<-c("Longitude","Latitude","MAP")
  #Now make the map
  pr<- ggplot() +
    geom_tile(data=aa, aes(y=Latitude, x=Longitude, fill=MAP)) +
    scale_fill_gradient2(low="#440154FF",high="white",mid="#2D708EfF",
                         midpoint = 0.5, na.value = "transparent",
                         breaks=c(0.5,0.75,1,1.25,1.5),
                         labels=c("0.5","0.75","1","1.25","1.5"),
                         limits=c(0.4,1.5), name="SD")+
    # add coastline
    geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
                 fill= "grey", color = "darkgrey") +
    coord_map(xlim = lons, ylim = lats)+
    # formatting
    scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
    ylab(" ")+xlab(" ")+
    theme_light() + ggtitle(paste0( j))+ labs(fill = "CPUE") +
    theme(axis.text.x= element_text(size=8), axis.text.y= element_text(size=8),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(paste("./5_Best_model_prediction/ggplot_maps/spatial/sd/",j,"pr.jpg"))
  
}

# Prediction -----------------------------------------------------------------

#shared<-readRDS("./5_Best_model_prediction/rds/Best_model_len_ves_st_yeariid_bathrw2.rds")# M1

index.p <-
  inla.stack.index(stk.full, tag = "pred")$data # indexes of the predicted values
output_mean<-exp(shared$summary.linear.predictor[index.p,"0.5quant"])- 0.0001# 
output_25<-exp(shared$summary.linear.predictor[index.p,"0.025quant"])- 0.0001# 
output_975<-exp(shared$summary.linear.predictor[index.p,"0.975quant"])- 0.0001# 

# Plot prediction --------------------------------------------------------------

library(ggplot2)
library(dplyr)
require(maps)
require(viridis)
library(marmap)

world_map <- map_data("world")
some.eu.countries <- c("Portugal", "Spain")

# Retrievethe map data
reg <- map_data("world", region = some.eu.countries)

# create the breaks- and label vectors
ewbrks <- seq(-10.5,-5.5,1)
nsbrks <- seq(37,42,1)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(-x, "°W"), 
                                                   ifelse(x > 0, paste(x, "°E"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "°S"), 
                                                   ifelse(x > 0, paste(x, "°N"),x))))

# get bathymetry data
b = getNOAA.bathy(lon1 = -9.84, lon2 = -8.5, lat1 = 41.82, lat2 = 37.1, 
                  resolution = 1)

# set map limits
# set map limits
lons = c(-10.5, -7)
lats = c(37, 42 )
bf = fortify.bathy(b)

# GAMMA MEDIAN: Change txt fil for pred_abun_median, pred_abun_q25, pred_abun_q75

abun_median<-output_mean
summary(abun_median)
m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         abun_median[1:m + (j - 1) * m])
  return(r) 
})
library(fields)
library(RColorBrewer)
#par(mfrow = c(3, 4), mar = c(0, 1.3, 1.2, 0.8))
(zlm <- range(unlist(m.prj), na.rm = TRUE)) # para definir breaks

for (j in 1:12) {
  
  
  library(raster)
  sp.mean.raster1<-raster(list(x=prj$x,
                               y = prj$y,
                               z = m.prj[[j]]))
  bath<-raster("./data/working_data/raster_env/GPS_envars/Bath.tif")
  res<-resample(bath,sp.mean.raster1)
  aa<-mask(sp.mean.raster1,res)
  aa <- rasterToPoints(aa)
  aa<-data.frame(aa)
  colnames(aa)<-c("Longitude","Latitude","MAP")
  #Now make the map
  pr<- ggplot() +
    geom_tile(data=aa, aes(y=Latitude, x=Longitude, fill=MAP)) +
    #scale_fill_distiller(palette = "YlOrRd", direction=1,
    scale_fill_viridis(alpha = 0.9, begin = 0, end = 1, direction = 1,
                      discrete = FALSE, option = "D",                    
     breaks=c(0,1000,2000,3000,4000,5000,6000,7000),
     labels=c("0","1000","2000","3000","4000","5000","6000","7000"),
                         limits=c(0,7050)) +
     geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
                 fill= "grey", color = "darkgrey") +
    coord_map(xlim = lons, ylim = lats)+
    # formatting
    scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
    ylab(" ")+xlab(" ")+
    theme_light() + ggtitle(paste0( j))+ labs(fill = "CPUE") +
    theme(axis.text.x= element_text(size=8), axis.text.y= element_text(size=8),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(paste("./5_Best_model_prediction/ggplot_maps/prediction/",j,"pr.jpg"))
  
  
}
