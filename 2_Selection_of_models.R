#************************
# Paper Sardine CPUE    #
# INLA model comparison #
#*****************************
# Francisco Izquierdo        #
#*****************************

## https://becarioprecario.bitbucket.io/spde-gitbook/
## To see the index press: Ctrl + shift + O

## All variables were scaled before including them into the modelling process
## Each of the next sections are additive components into the model
## 0) Select effort variable and also if linear or rw2 effect
## 1) Select the best ID variable (iid random effect)
## 234) Select different spatio-temporal variables
## We try the spatio-temporal structured AR(1) progressive term (Paradinas et al., 2015)
## We try the inclusion of year as iid term
## We try the inclusion of month as RW2 cyclic term
## 5) Select environmental variables
## 6) Select environmental variables but removing the ST component to avoid confounding
## Finally, *Best model* is selected

rm(list=ls()) # Clean workspace
library(INLA)
inla.setOption(scale.model.default = TRUE) # set scale.model=TRUE, see scale tutorial

# 1) Base model selection

# 2) Environmental variable selection

# 3) Model comparison

# Model structure ---------------------------------------------------------

data <- read.table("./data/working_data/2_data_obs_covs_GPS_escaled.txt" , 
                   header=TRUE, dec=".")
data$year<-as.factor(data$year)

## ---- mesh
coords4<-as.matrix(data[,12:13])
length(coords4[,1])
bound=inla.nonconvex.hull(as.matrix(data[,12:13]), convex=-0.07, eps=0.05, resolution=40)
mesh.s <- inla.mesh.2d(loc=coords4, boundary=bound, max.edge=c(0.05, 0.6), offset=c(0.1,0.6), cutoff=0.08, min.angle = 0.05)
plot(mesh.s,asp=1)
points(coords4, pch=19, cex=0.7)
mesh.s$n

## ---- spde
spde <- inla.spde2.pcmatern(mesh.s, prior.range=c(0.5, 0.05), prior.sigma=c(0.6, 0.05))
(m <- spde$n.spde)

# ----- AMat
Ast <- inla.spde.make.A(mesh = mesh.s, loc = coords4,
                        group = data$month_id,
                        n.group = length(unique(data$month_id)))
dim(Ast)

## ---- idx
idx<- inla.spde.make.index("s", n.spde=spde$n.spde, n.group=max(data$month_id))

# ----- Priors
hyper.prec <- list(theta = list(prior="pc.prec", 
                                param = c(1, 0.05)))#allow smaller values prec iid
hyper.prec.bath <- list(theta = list(prior="pc.prec", 
                                     param = c(1, 0.05)))# good for bathy
hyper.prec.chl <- list(theta = list(prior="pc.prec", 
                                    param = c(0.5, 0.05)))# good results with U = 0.5
hyper.prec.sst <- list(theta = list(prior="pc.prec", 
                                    param = c(1, 0.05))) # U = 1
pcrho <- list(prior = 'pccor1', param = c(0.7, 0.7)) # rho


# 0) Effort vars --------------------------------------------------------------

## M 0 length fixed -----------------------------------------------------------------------

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=data$length)), tag="est")

form <- y ~ 0 + a0  + length  

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M0_length_fix.rds")



## M 0 length rw2 -----------------------------------------------------------------------

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6))), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) 

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M0_length_rw2.rds")

## M 0 power fixed -----------------------------------------------------------------------

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          power=data$power)), tag="est")

form <- y ~ 0 + a0  + power  

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M0_power_fix.rds")



## M 0 power rw2 -----------------------------------------------------------------------

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          power=inla.group(data$power, n=6))), tag="est")

form <- y ~ 0 + a0 + f(power, model="rw2", hyper=hyper.prec) 

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M0_power_rw2.rds")


## M 0 ton fixed -----------------------------------------------------------------------

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          ton=data$ton)), tag="est")

form <- y ~ 0 + a0  + ton  

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M0_ton_fix.rds")


## M 0 ton rw2 -----------------------------------------------------------------------

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          ton=inla.group(data$ton, n=6))), tag="est")

form <- y ~ 0 + a0 + f(ton, model="rw2", hyper=hyper.prec) 

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M0_ton_rw2.rds")

# *0 Comp ---------------------------------------------------------------------------------

library(INLA)

mod1<- readRDS("./4_base_model_standardization/models_rds/M0_length_fix.rds")
mod2<- readRDS("./4_base_model_standardization/models_rds/M0_length_rw2.rds")
mod3<- readRDS("./4_base_model_standardization/models_rds/M0_power_fix.rds")
mod4<- readRDS("./4_base_model_standardization/models_rds/M0_power_rw2.rds")
mod5<- readRDS("./4_base_model_standardization/models_rds/M0_ton_fix.rds")
mod6<- readRDS("./4_base_model_standardization/models_rds/M0_ton_rw2.rds")

COMP_iid<-data.frame(
  model=c("len_fix","len_rw2","pow_fix","pow_rw2","ton_fix","ton_rw2"),
  dic=c(mod1$dic$dic,mod2$dic$dic,mod3$dic$dic,mod4$dic$dic,mod5$dic$dic,mod6$dic$dic),
  waic=c(mod1$waic$waic,mod2$waic$waic,mod3$waic$waic,mod4$waic$waic,mod5$waic$waic,mod6$waic$waic),
  lcpo=c(-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T),-mean(log(mod4$cpo$cpo),na.rm=T),-mean(log(mod5$cpo$cpo),na.rm=T),-mean(log(mod6$cpo$cpo),na.rm=T)),
  failure=c(sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T),sum((mod4$cpo$failure>0)*1,na.rm=T),sum((mod5$cpo$failure>0)*1,na.rm=T),sum((mod6$cpo$failure>0)*1,na.rm=T)),
  time=c(mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4],mod4$cpu.used[4],mod5$cpu.used[4],mod6$cpu.used[4])
)
COMP_iid

write.table(COMP_iid ,"./4_base_model_standardization/models_rds/Comp_0.txt" )

# 1) ID vars -----------------------------------------------------------------------

## M 1 eff + ves ID -----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec)

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M1_eff_vessel.rds")

## M 1 eff + har ID -----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          harbour=data$harbour_id)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(harbour, model="iid", hyper=hyper.prec)

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M1_eff_harbour.rds")


## M 1 eff + reg ID -----------------------------------------------------------------------

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),region=data$region_id)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(region, model="iid", hyper=hyper.prec)

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M1_eff_region.rds")

## M 1 eff + vessel + vessinhar-----------------------------------------------------------------------

unique(data$vessel_id)
length(unique(data$vessel_id))
unique(data$harbour_id)
length(unique(data$harbour_id))
(data$vessinhar<-factor(paste0(data$harbour_id,data$vessel_id)))
length(unique(data$vessinhar))

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id,
                                          vessinhar=data$vessel_id)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(vessinhar, model="iid", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M1_eff_vess_vessinhar.rds")


## M 1 eff + vessel + vessinreg-----------------------------------------------------------------------

# NESTED: harbour associated to the correspondent vessels
unique(data$vessel_id)
length(unique(data$vessel_id))
unique(data$harbour_id)
length(unique(data$harbour_id))
(data$vessinreg<-factor(paste0(data$region_id,data$vessel_id)))
length(unique(data$vessinhar))

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          vessinreg=data$vessinreg)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(vessinreg, model="iid", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M1_eff_vess_vessinreg.rds")

## M 1 eff + harbour + harinreg-----------------------------------------------------------------------

# NESTED: each region associated to the correspondent vessels

data$harinreg<-factor(paste0(data$region_id,data$harbour_id))

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          harinreg=data$harinreg)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(harinreg, model="iid", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M1_eff_har_harinreg.rds")

# *1 Comp ------------------------------------------------------------------------------------

mod1<- readRDS("./4_base_model_standardization/models_rds/M1_eff_vessel.rds")
mod2<- readRDS("./4_base_model_standardization/models_rds/M1_eff_vess_vessinreg.rds")
mod3<- readRDS("./4_base_model_standardization/models_rds/M1_eff_vess_vessinhar.rds")
mod4<- readRDS("./4_base_model_standardization/models_rds/M1_eff_region.rds")
mod5<- readRDS("./4_base_model_standardization/models_rds/M1_eff_harbour.rds")
mod6<- readRDS("./4_base_model_standardization/models_rds/M1_eff_har_harinreg.rds")

COMP_iid<-data.frame(
  model=c("vess","vessinreg","vessinhar","region","harbour","harinreg"),
  dic=c(mod1$dic$dic,mod2$dic$dic,mod3$dic$dic,mod4$dic$dic,mod5$dic$dic,mod6$dic$dic),
  waic=c(mod1$waic$waic,mod2$waic$waic,mod3$waic$waic,mod4$waic$waic,mod5$waic$waic,mod6$waic$waic),
  lcpo=c(-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T),-mean(log(mod4$cpo$cpo),na.rm=T),-mean(log(mod5$cpo$cpo),na.rm=T),-mean(log(mod6$cpo$cpo),na.rm=T)),
  failure=c(sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T),sum((mod4$cpo$failure>0)*1,na.rm=T),sum((mod5$cpo$failure>0)*1,na.rm=T),sum((mod6$cpo$failure>0)*1,na.rm=T)),
  time=c(mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4],mod4$cpu.used[4],mod5$cpu.used[4],mod6$cpu.used[4])
)
COMP_iid

write.table(COMP_iid ,"./4_base_model_standardization/models_rds/Comp_1.txt" )

# 234) STemp vars -----------------------------------------------------------------------

## M 2 eff + ves + st -----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) 


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M2_eff_ves_st.rds")

## M 3 eff + ves + st + year iid-----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, year=data$year_id)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M3_eff_ves_st_yeariid.rds")


## M 3 eff + ves + st + year rw2-----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, year=data$year_id)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="rw2", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M3_eff_ves_st_yearrw2.rds")


## M 4 eff + ves + st + year + month rw2-----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          month=data$month_id)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  f(month, model="rw2", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M4_eff_ves_st_yeariid_monthrw2.rds")


## M 4 eff + ves + st + year + month rw2 cyclic-----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          month=data$month_id)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  f(month, model="rw2", hyper=hyper.prec, cyclic=TRUE)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M4_eff_ves_st_yeariid_monthrw2_cyclic.rds")


# *234 Comp -------------------------------------------------------------------------------------------


mod1<- readRDS("./4_base_model_standardization/models_rds/M2_eff_ves_st.rds")
mod2<- readRDS("./4_base_model_standardization/models_rds/M3_eff_ves_st_yearrw2.rds")
mod3<- readRDS("./4_base_model_standardization/models_rds/M3_eff_ves_st_yeariid.rds")
mod4<- readRDS("./4_base_model_standardization/models_rds/M4_eff_ves_st_yeariid_monthrw2.rds")
mod5<- readRDS("./4_base_model_standardization/models_rds/M4_eff_ves_st_yeariid_monthrw2_cyclic.rds")

COMP_iid<-data.frame(
  model=c("st","st yearrw2","st yeariid","st yeariid month","st yeariid month cycl"),
  dic=c(mod1$dic$dic,mod2$dic$dic,mod3$dic$dic,mod4$dic$dic,mod5$dic$dic),
  waic=c(mod1$waic$waic,mod2$waic$waic,mod3$waic$waic,mod4$waic$waic,mod5$waic$waic),
  lcpo=c(-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T),-mean(log(mod4$cpo$cpo),na.rm=T),-mean(log(mod5$cpo$cpo),na.rm=T)),
  failure=c(sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T),sum((mod4$cpo$failure>0)*1,na.rm=T),sum((mod5$cpo$failure>0)*1,na.rm=T)),
  time=c(mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4],mod4$cpu.used[4],mod5$cpu.used[4])
)
COMP_iid

write.table(COMP_iid ,"./4_base_model_standardization/models_rds/Comp_1.txt" )

# 5) Envars ----------------------------------------------------------------

## M 5 eff + ves + st + year + bath rw2----------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          bath=inla.group(data$bath,n=10, 
                                                          method="quantile"))), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  f(bath, model="rw2", hyper=hyper.prec.bath)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M5_eff_ves_st_yeariid_bathrw2.rds")


## M 5 eff + ves + st + year + bath L----------------------------------------------------------------

table(inla.group(data$bath,n=12, idx=FALSE))

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          bath=data$bath)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  bath

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M5_eff_ves_st_yeariid_bathL.rds")



## M 5 eff + ves + st + year + chl rw2----------------------------------------------------------------

table(inla.group(data$chl_esc,n=10, idx=TRUE, method="quantile"))

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          chl=inla.group(data$chl_esc,n=10, 
                                                         method="quantile"))), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  f(chl, model="rw2", hyper=hyper.prec.bath)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M5_eff_ves_st_yeariid_chlrw2.rds")

## M 5 eff + ves + st + year + chl L-----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          chl=data$chl_esc)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  chl


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M5_eff_ves_st_yeariid_chlL.rds")


## M 5 eff + ves + st + year + sst rw2----------------------------------------------------------------

table(inla.group(data$sst_esc,n=10, idx=TRUE, method="quantile"))


effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          sst=inla.group(data$sst_esc,n=10, 
                                                         method="quantile"))), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  f(sst, model="rw2", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M5_eff_ves_st_yeariid_sstrw2.rds")


## M 5 eff + ves + st + year + sst L-----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          sst=data$sst_esc)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  sst


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M5_eff_ves_st_yeariid_sstL.rds")


## M 5 eff + ves + st + year + U rw2----------------------------------------------------------------



stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          currU=inla.group(data$curr_u,n=10, 
                                                           method="quantile"))), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  f(currU, model="rw2", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M5_eff_ves_st_yeariid_currUrw2.rds")


## M 5 eff + ves + st + year + currU L-----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          currU=data$curr_u)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  currU


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M5_eff_ves_st_yeariid_currUL.rds")


## M 5 eff + ves + st + year + V rw2----------------------------------------------------------------



stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          currV=inla.group(data$curr_v,n=10, 
                                                           method="quantile"))), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  f(currV, model="rw2", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M5_eff_ves_st_yeariid_currVrw2.rds")


## M 5 eff + ves + st + year + currV L-----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          year=data$year_id,
                                          currV=data$curr_v)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(year, model="iid", hyper=hyper.prec) +
  currV


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M5_eff_ves_st_yeariid_currVL.rds")



# *5 Comp ---------------------------------------------------------------------------------------


mod1<- readRDS("./4_base_model_standardization/models_rds/M6_eff_ves_st_bathrw2.rds")
mod2<- readRDS("./4_base_model_standardization/models_rds/M6_eff_ves_st_bathL.rds")
mod3<- readRDS("./4_base_model_standardization/models_rds/M6_eff_ves_st_chlrw2.rds")
mod4<- readRDS("./4_base_model_standardization/models_rds/M6_eff_ves_st_chlL.rds")
mod5<- readRDS("./4_base_model_standardization/models_rds/M6_eff_ves_st_sstrw2.rds")
mod6<- readRDS("./4_base_model_standardization/models_rds/M6_eff_ves_st_sstL.rds")
mod7<- readRDS("./4_base_model_standardization/models_rds/M6_eff_ves_st_currUrw2.rds")
mod8<- readRDS("./4_base_model_standardization/models_rds/M6_eff_ves_st_currUL.rds")
mod9<- readRDS("./4_base_model_standardization/models_rds/M6_eff_ves_st_currVrw2.rds")
mod10<- readRDS("./4_base_model_standardization/models_rds/M6_eff_ves_st_currVL.rds")


COMP_iid<-data.frame(
  model=c("b","br","c","cr","s","sr","u","ur","v","vr"),
  dic=c(mod1$dic$dic,mod2$dic$dic,mod3$dic$dic,mod4$dic$dic,mod5$dic$dic,mod6$dic$dic,mod7$dic$dic,mod8$dic$dic,mod9$dic$dic,mod10$dic$dic),
  waic=c(mod1$waic$waic,mod2$waic$waic,mod3$waic$waic,mod4$waic$waic,mod5$waic$waic,mod6$waic$waic,mod7$waic$waic,mod8$waic$waic,mod9$waic$waic,mod10$waic$waic),
  lcpo=c(-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T),-mean(log(mod4$cpo$cpo),na.rm=T),-mean(log(mod5$cpo$cpo),na.rm=T),-mean(log(mod6$cpo$cpo),na.rm=T),-mean(log(mod7$cpo$cpo),na.rm=T),-mean(log(mod8$cpo$cpo),na.rm=T),-mean(log(mod9$cpo$cpo),na.rm=T),-mean(log(mod10$cpo$cpo),na.rm=T)),
  failure=c(sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T),sum((mod4$cpo$failure>0)*1,na.rm=T),sum((mod5$cpo$failure>0)*1,na.rm=T),sum((mod6$cpo$failure>0)*1,na.rm=T),sum((mod7$cpo$failure>0)*1,na.rm=T),sum((mod8$cpo$failure>0)*1,na.rm=T),sum((mod9$cpo$failure>0)*1,na.rm=T),sum((mod10$cpo$failure>0)*1,na.rm=T)),
  time=c(mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4],mod4$cpu.used[4],mod5$cpu.used[4],mod6$cpu.used[4],mod7$cpu.used[4],mod8$cpu.used[4],mod9$cpu.used[4],mod10$cpu.used[4])
)
COMP_iid

write.table(COMP_iid ,"./5_model_vessel_standarization/5_rds_nuevo/comp_mods.txt" )

# 6) Envars no ST -------------------------------------------------------------------------

## We try environmental variables without Spatio-Temporal term
## In order to check that there is no confounding among ST and envars
## Environmental variables were not relevant in any of the cases

## M 6 eff + ves  + year + bath rw2----------------------------------------------------------------

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          bath=inla.group(data$bath,n=10, 
                                                          method="quantile"))), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(bath, model="rw2", hyper=hyper.prec.bath)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M6_eff_ves_st__bathrw2.rds")


## M 6 eff + ves  + year + bath L----------------------------------------------------------------

table(inla.group(data$bath,n=12, idx=FALSE))

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          bath=data$bath)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  bath

abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M6_eff_ves_st_bathL.rds")



## M 6 eff + ves  + year + chl rw2----------------------------------------------------------------

table(inla.group(data$chl_esc,n=10, idx=TRUE, method="quantile"))

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          chl=inla.group(data$chl_esc,n=10, 
                                                         method="quantile"))), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(chl, model="rw2", hyper=hyper.prec.bath)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M6_eff_ves_st_chlrw2.rds")

## M 6 eff + ves  + year + chl L-----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          chl=data$chl_esc)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  chl


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M6_eff_ves_st_chlL.rds")


## M 6 eff + ves  + year + sst rw2----------------------------------------------------------------

table(inla.group(data$sst_esc,n=10, idx=TRUE, method="quantile"))


effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          sst=inla.group(data$sst_esc,n=10, 
                                                         method="quantile"))), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(sst, model="rw2", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M6_eff_ves_st_sstrw2.rds")


## M 6 eff + ves  + year + sst L-----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          sst=data$sst_esc)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  sst


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M6_eff_ves_st_sstL.rds")


## M 6 eff + ves  + year + U rw2----------------------------------------------------------------



stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          currU=inla.group(data$curr_u,n=10, 
                                                           method="quantile"))), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(currU, model="rw2", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M6_eff_ves_st_currUrw2.rds")


## M 6 eff + ves  + year + currU L-----------------------------------------------------------------------

effvar<-data$len_esc # chose selected one, this time scaled

stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          currU=data$curr_u)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  currU


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M6_eff_ves_st_currUL.rds")


## M 6 eff + ves  + year + V rw2----------------------------------------------------------------



stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          currV=inla.group(data$curr_v,n=10, 
                                                           method="quantile"))), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  f(currV, model="rw2", hyper=hyper.prec)


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M6_eff_ves_st_currVrw2.rds")


## M 6 eff + ves  + year + currV L-----------------------------------------------------------------------


stk.e <- inla.stack(data=list(y=data$sardine, link=1), 
                    A=list(Ast, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$sardine)), 
                                          length=inla.group(data$length, n=6),
                                          vessel=data$vessel_id, 
                                          currV=data$curr_v)), tag="est")

form <- y ~ 0 + a0 + f(length, model="rw2", hyper=hyper.prec) +
  f(vessel, model="iid", hyper=hyper.prec) +
  currV


abu.res <- inla(form, family = 'gamma', 
                data = inla.stack.data(stk.e),
                control.predictor = list(A = inla.stack.A(stk.e),link=1),
                control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
                control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./4_base_model_standardization/models_rds/M6_eff_ves_st_currVL.rds")


# *Best model*------------------------------------------------------------------

## More models were tried (not included here)
## M4 (Length (rw2) + vessel(iid) + spatio-temporal(AR1) + year(iid) + month(rw2 cyclic)
