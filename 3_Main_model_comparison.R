rm(list=ls()) # Clean workspace
library(INLA)

# 1) Base model selection

# 2) Environmental variable selection

# 3) Model comparison

# Model structure ---------------------------------------------------------

## ----data-----------

data <- read.table("./data/working_data/2_data_obs_joined_5_covs_scaled.txt", header=TRUE, dec=".")

## ----mesh-----------

coords4<-as.matrix(data[,1:2])
RangeGuess<- 100
MaxEdge<- RangeGuess/5
ConvHull <- inla.nonconvex.hull(coords4,convex=-0.06)
mesh.s <- inla.mesh.2d(coords4,
                       boundary = ConvHull, max.edge = c(1,4) * MaxEdge, cutoff=MaxEdge/5)
plot(mesh.s, asp=1) 
points(coords4[,1:2], pch=19, cex=0.7)
mesh.s$n

## ----spde-----------

Range0<-90# we force inla to come up with a range larger than 90 km
AlphaRange0<- 0.05
Sigma0<- 0.6 
AlphaSigma0 <-0.6
spde<-inla.spde2.pcmatern(mesh.s, prior.range = c(Range0, AlphaRange0), 
                          prior.sigma = c(Sigma0,AlphaSigma0))
m <- spde$n.spde
m

# A MAT -------------

Ast <- inla.spde.make.A(mesh = mesh.s, loc = coords4,
                        group = data$month_id,
                        n.group = length(unique(data$month_id)))
dim(Ast)

## ----idx----------

idx<- inla.spde.make.index("s", n.spde=spde$n.spde, n.group=max(data$month_id))

# PRIORS  ----------

hyper.prec <- list(theta = list(prior="pc.prec", 
                                param = c(1, 0.05)))#allow smaller values prec iid
hyper.prec.bath <- list(theta = list(prior="pc.prec", 
                                     param = c(1, 0.05)))# good for bathy
hyper.prec.chl <- list(theta = list(prior="pc.prec", 
                                    param = c(0.5, 0.05)))# Good results with U = 1 o 0.5
hyper.prec.sst <- list(theta = list(prior="pc.prec", 
                                    param = c(1, 0.05))) # U = 1
pcrho <- list(prior = 'pccor1', param = c(0.7, 0.7)) # rho


#1) BASE MODEL SELECTION --------------------------------------------------------------

# All variables were scaled before including them into the modeling process

# Length is the effort variable as fixed

# Vessel ID is the IID effect as random

# We try the spatio-temporal structured AR(1) progressive term (Paradinas et al., 2015)

# We try the inclusion of year as iid term

# We try the inclusion of month as RW2 term

# Year as RW2 and some other possibilities were also tested (not included)

# M 0 -------------------------------------------------------------------

# Length esc

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=data$length_esc, 
                                          vessel=data$vessel_id)), tag="est")


form <- y ~ 0 + a0  + length  

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./5_model_vessel_standarization/5_rds_nuevo/M0_0_length.rds")


# M 1 -----------------------------------------------------------------------

# Length + vessel 

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=data$length_esc, 
                                          vessel=data$vessel_id)), tag="est")


form <- y ~ 0 + a0  + length  + f(vessel, model="iid", hyper=hyper.prec) 

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./5_model_vessel_standarization/5_rds_nuevo/M0_length_vessel.rds")

# M 2  -----------------------------------------------------------------------

#length + vessel + spatio-temporal

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=data$length_esc, 
                                          vessel=data$vessel_id)), tag="est")

form <- y ~ 0 + a0  + 
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho)))+
  f(vessel, model="iid", hyper=hyper.prec)

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./5_model_vessel_standarization/5_rds_nuevo/M2_length_spat_vessel.rds")




# M 3 -----------------------------------------------------------------------

# Length + vessel + spatio-temporal + year

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=data$length_esc, 
                                          vessel=data$vessel_id,year=data$year_id)), tag="est")


form <- y ~ 0 + a0  + length  + f(vessel, model="iid", hyper=hyper.prec) +  
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +  f(year, model="iid", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./5_model_vessel_standarization/5_rds_nuevo/M3_length_spat_vessel.rds")


# M 4 (Best Model) -----------------------------------------------------------------------

# Length + vessel + spatio-temporal + year + month

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=data$length_esc, 
                                          vessel=data$vessel_id,year=data$year_id, month=data$month_id)), tag="est")


form <- y ~ 0 + a0  + length  + f(vessel, model="iid", hyper=hyper.prec) +  
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +  f(year, model="iid", hyper=hyper.prec) +
  f(month, model="rw2", hyper=hyper.prec)  


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./5_model_vessel_standarization/5_rds_nuevo/M4_length_spat_vessel_month.rds")




# 2) ENVIRONMENTAL VARIABLE SELECTION --------------------------------------------------------

# NOW, THE BEST SELECTED MODEL IS M4

# M4 includes: LENGTH + VESSEL + SPATIO-TEMPORAL + YEAR + MONTH

# We try the 5 environmental variables (Chl-a, SST, Bathymetry, Intensity and Angle)

# We tried them extracting the sampling points over the MEAN and also the MEDIAN raster maps

# M 5 -----------------------------------------------------------------------

# CHL

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          chl=inla.group(data$chl_mean_esc,n=9), length=data$length_esc, 
                                          vessel=data$vessel_id, year=data$year_id, month=data$month_id)), tag="est")


form <- y ~ 0 + a0  + length  + f(vessel, model="iid", hyper=hyper.prec) +  
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(chl, model="rw2", hyper=hyper.prec.chl) +
  f(year, model="iid", hyper=hyper.prec) +
  f(month, model="rw2", hyper=hyper.prec)  

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./5_model_vessel_standarization/5_rds_nuevo/M4_chl.rds")

# M 6 -----------------------------------------------------------------------

# SST

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          sst=inla.group(data$sst_mean_esc,n=9), length=data$length_esc, 
                                          vessel=data$vessel_id, year=data$year_id, month=data$month_id)), tag="est")


form <- y ~ 0 + a0  + length  + f(vessel, model="iid", hyper=hyper.prec) +  
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(sst, model="rw2", hyper=hyper.prec.sst) +
  f(year, model="iid", hyper=hyper.prec) +
  f(month, model="rw2", hyper=hyper.prec)  

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./5_model_vessel_standarization/5_rds_nuevo/M5_sst.rds")

# M 7 -----------------------------------------------------------------------

# BATHYMETRY

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          bath=inla.group(data$bath_esc,n=6), length=data$length_esc, 
                                          vessel=data$vessel_id, year=data$year_id, month=data$month_id)), tag="est")


form <- y ~ 0 + a0  + length  + f(vessel, model="iid", hyper=hyper.prec) +  
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(bath, model="rw2", hyper=hyper.prec.bath) +
  f(year, model="iid", hyper=hyper.prec) +
  f(month, model="rw2", hyper=hyper.prec)  

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./5_model_vessel_standarization/5_rds_nuevo/M6_bath.rds")


# M 8  -----------------------------------------------------------------------

# INTENSITY

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          int=inla.group(data$int_mean_esc,n=10), length=data$length, 
                                          vessel=data$vessel_id, year=data$year_id, month=data$month_id)), tag="est")


form <- y ~ 0 + a0  + length  + f(vessel, model="iid", hyper=hyper.prec) +  
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(int, model="rw2", hyper=hyper.prec.chl)  +
  f(year, model="iid", hyper=hyper.prec) +
  f(month, model="rw2", hyper=hyper.prec)  

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./5_model_vessel_standarization/5_rds_nuevo/M7_int.rds")




# M 9 -----------------------------------------------------------------------

# ANGLE

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          ang=inla.group(data$ang_mean_esc,n=10), length=data$length, 
                                          vessel=data$vessel_id, year=data$year_id, month=data$month_id)), tag="est")


form <- y ~ 0 + a0  + length  + f(vessel, model="iid", hyper=hyper.prec) +  
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(ang, model="rw2", hyper=hyper.prec.chl)  +
  f(year, model="iid", hyper=hyper.prec) +
  f(month, model="rw2", hyper=hyper.prec)  

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./5_model_vessel_standarization/5_rds_nuevo/M8_angle.rds")




# M 10 -----------------------------------------------------------------------

# int + angle

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          int=inla.group(data$int_mean_esc,n=10), ang=inla.group(data$ang_mean_esc, n=10),
                                          length=data$length_esc, 
                                          vessel=data$vessel_id, year=data$year_id, month=data$month_id)), tag="est")


form <- y ~ 0 + a0  + length  + f(vessel, model="iid", hyper=hyper.prec) +  
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(int, model="rw2", hyper=hyper.prec.chl) +   f(ang, model="rw2", hyper=hyper.prec.chl) +
  f(year, model="iid", hyper=hyper.prec) +
  f(month, model="rw2", hyper=hyper.prec)  


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./5_model_vessel_standarization/5_rds_nuevo/M9_int_angle.rds")



# 3) MODEL COMPARISON -----------------------------------------------------------------------

library(INLA)

mod1<- readRDS("./5_model_vessel_standarization/5_rds_nuevo/M0_0_length.rds")
mod2<- readRDS("./5_model_vessel_standarization/5_rds_nuevo/M2_length_spat_vessel.rds")
mod3<- readRDS("./5_model_vessel_standarization/5_rds_nuevo/M3_length_spat_vessel.rds")
mod4<- readRDS("./5_model_vessel_standarization/5_rds_nuevo/M4_length_spat_vessel_month.rds")
mod5<- readRDS("./5_model_vessel_standarization/5_rds_nuevo/M5_sst.rds")
mod6<- readRDS("./5_model_vessel_standarization/5_rds_nuevo/M6_bath.rds")
mod7<- readRDS("./5_model_vessel_standarization/5_rds_nuevo/M7_int.rds")
mod8<- readRDS("./5_model_vessel_standarization/5_rds_nuevo/M8_angle.rds")
mod9<- readRDS("./5_model_vessel_standarization/5_rds_nuevo/M9_int_angle.rds")


COMP_iid<-data.frame(
  model=c("m1_length_vessel","m2_length_spat_vessel","m3_length_spat_vessel_year","m4_month","m5_sst","m6_bath","m7_int","m8_ang","m9_int_ang"),
  dic=c(mod1$dic$dic,mod2$dic$dic,mod3$dic$dic,mod4$dic$dic,mod5$dic$dic,mod6$dic$dic,mod7$dic$dic,mod8$dic$dic,mod9$dic$dic),
  waic=c(mod1$waic$waic,mod2$waic$waic,mod3$waic$waic,mod4$waic$waic,mod5$waic$waic,mod6$waic$waic,mod7$waic$waic,mod8$waic$waic,mod9$waic$waic),
  lcpo=c(-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T),-mean(log(mod4$cpo$cpo),na.rm=T),-mean(log(mod5$cpo$cpo),na.rm=T),-mean(log(mod6$cpo$cpo),na.rm=T),-mean(log(mod7$cpo$cpo),na.rm=T),-mean(log(mod8$cpo$cpo),na.rm=T),-mean(log(mod9$cpo$cpo),na.rm=T)),
  failure=c(sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T),sum((mod4$cpo$failure>0)*1,na.rm=T),sum((mod5$cpo$failure>0)*1,na.rm=T),sum((mod6$cpo$failure>0)*1,na.rm=T),sum((mod7$cpo$failure>0)*1,na.rm=T),sum((mod8$cpo$failure>0)*1,na.rm=T),sum((mod9$cpo$failure>0)*1,na.rm=T)),
  time=c(mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4],mod4$cpu.used[4],mod5$cpu.used[4],mod6$cpu.used[4],mod7$cpu.used[4],mod8$cpu.used[4],mod9$cpu.used[4])
)
COMP_iid

write.table(COMP_iid ,"./5_model_vessel_standarization/5_rds_nuevo/comp_mods.txt" )

# THE FINAL SELECTED MODEL OF THIS STUDY WAS M4 (Length + vessel + spatio-temporal + year + month)
