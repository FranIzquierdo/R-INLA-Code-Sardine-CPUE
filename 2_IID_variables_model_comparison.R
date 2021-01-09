rm(list=ls()) # Clean workspace
library(INLA)


#1) IID MODELS  -----------------------------------------------

## PRIOR for iid effects
## This prior allow smaller values prec iid (from SPDE book)

hyper.prec <- list(theta = list(prior="pc.prec", 
                                param = c(1, 0.05)))

### Region ID

I1 <- inla(sardine ~ 1 +
             length + 
             f(region_id, model = "iid", hyper = hyper.prec), family="Gamma",
           control.predictor=list(compute=TRUE),
           control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE),
           num.threads=2,
           control.inla=list(strategy="gaussian"),
           verbose=TRUE,
           data=data)

saveRDS(I1, "./5_model_vessel_standarization/random_effects/length_region_iid.rds")

### Vessel ID

I2 <- inla(sardine ~ 1 +
             length + 
             f(vessel_id, model = "iid", hyper = hyper.prec), family="Gamma",
           control.predictor=list(compute=TRUE),
           control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE),
           num.threads=2,
           control.inla=list(strategy="gaussian"),
           verbose=TRUE,
           data=data)

saveRDS(I2, "./5_model_vessel_standarization/random_effects/length_vessel_iid.rds")

# Harbour ID

I3 <- inla(sardine ~ 1 +
             length + 
             f(harbour_id, model = "iid", hyper = hyper.prec), family="Gamma",
           control.predictor=list(compute=TRUE),
           control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE),
           num.threads=2,
           control.inla=list(strategy="gaussian"),
           verbose=TRUE,
           data=data)

saveRDS(I3, "./5_model_vessel_standarization/random_effects/length_harbour_iid.rds")


# 2) IID NESTED MODELS -----------------------------------------------------

# Just the most important model combinations are presented here (more models were tried)

length(unique(data$vessel_id))# 66
length(unique(data$harbour_id))# 13
length(unique(data$region_id))# 3

### NESTED: harbour associated to the correspondent vessels

data$vessinhar<-factor(paste0(data$harbour_id,data$vessel_id))

I41 <- inla(sardine ~ 1 +
              length + 
              f(vessinhar, model="iid"), family="Gamma",
            control.predictor=list(compute=TRUE),
            control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE),
            num.threads=2,
            control.inla=list(strategy="gaussian"),
            verbose=TRUE,
            data=data)

saveRDS(I41, "./5_model_vessel_standarization/random_effects/length_vessinhar_iid.rds")

# NESTED harbour vessel + harbour

I42 <- inla(sardine ~ 1 +
              length + 
              f(harbour_id, model = "iid", hyper = hyper.prec) +
              f(vessinhar, model="iid"), family="Gamma",
            control.predictor=list(compute=TRUE),
            control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE),
            num.threads=2,
            control.inla=list(strategy="gaussian"),
            verbose=TRUE,
            data=data)

saveRDS(I42, "./5_model_vessel_standarization/random_effects/ength_harbour_vessinhar_iid.rds")

# NESTED harbour vessel + vessel

I43 <- inla(sardine ~ 1 +
              length + 
              f(vessel_id, model = "iid", hyper = hyper.prec) +
              f(vessinhar, model="iid"), family="Gamma",
            control.predictor=list(compute=TRUE),
            control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE),
            num.threads=2,
            control.inla=list(strategy="gaussian"),
            verbose=TRUE,
            data=data)

saveRDS(I43, "./5_model_vessel_standarization/random_effects/length_vessel_vessinhar_iid.rds")

### NESTED: each region associated to the correspondent vessels

data$vessinreg<-factor(paste0(data$region_id,data$vessel_id))

I51 <- inla(sardine ~ 1 +
              length + 
              f(vessinreg, model="iid"), family="Gamma",
            control.predictor=list(compute=TRUE),
            control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE),
            num.threads=2,
            control.inla=list(strategy="gaussian"),
            verbose=TRUE,
            data=data)

saveRDS(I51, "./5_model_vessel_standarization/random_effects/length_vessinreg_iid.rds")

# NESTED vessel in region + region

I52 <- inla(sardine ~ 1 +
              length + 
              f(region_id, model = "iid", hyper = hyper.prec) +
              f(vessinreg, model="iid"), family="Gamma",
            control.predictor=list(compute=TRUE),
            control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE),
            num.threads=2,
            control.inla=list(strategy="gaussian"),
            verbose=TRUE,
            data=data)

saveRDS(I52, "./5_model_vessel_standarization/random_effects/length_region_vessinreg_iid.rds")


# NESTED vessel in region + vessel

I53 <- inla(sardine ~ 1 +
              length + 
              f(vessel_id, model = "iid", hyper = hyper.prec) +
              f(vessinreg, model="iid"), family="Gamma",
            control.predictor=list(compute=TRUE),
            control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE),
            num.threads=2,
            control.inla=list(strategy="gaussian"),
            verbose=TRUE,
            data=data)

saveRDS(I53, "./5_model_vessel_standarization/random_effects/length_vessel_vessinreg_iid.rds")

# 3) MODEL COMPARISON --------------------------------------------------------------

mod1<- readRDS("./5_model_vessel_standarization/random_effects/length_region_iid.rds")
mod2<- readRDS("./5_model_vessel_standarization/random_effects/length_vessel_iid.rds")
mod3<- readRDS("./5_model_vessel_standarization/random_effects/length_harbour_iid.rds")
mod4<- readRDS("./5_model_vessel_standarization/random_effects/length_vessel_vessinhar_iid.rds")
mod5<- readRDS("./5_model_vessel_standarization/random_effects/length_harbour_vessinhar_iid.rds")
mod6<- readRDS("./5_model_vessel_standarization/random_effects/length_vessinhar_iid.rds")
mod7<- readRDS("./5_model_vessel_standarization/random_effects/length_vessinreg_iid.rds")
mod8<- readRDS("./5_model_vessel_standarization/random_effects/length_vessel_vessinreg_iid.rds")
mod9<- readRDS("./5_model_vessel_standarization/random_effects/length_region_vessinreg_iid.rds")

COMP_iid<-data.frame(
  model=c("reg","vess","har","vess_vesshar","har_vesshar","vesshar","vessreg","vess_vessreg","reg_vessreg"),
  dic=c(mod1$dic$dic,mod2$dic$dic,mod3$dic$dic,mod4$dic$dic,mod5$dic$dic,mod6$dic$dic,mod7$dic$dic,mod8$dic$dic,mod9$dic$dic),
  waic=c(mod1$waic$waic,mod2$waic$waic,mod3$waic$waic,mod4$waic$waic,mod5$waic$waic,mod6$waic$waic,mod7$waic$waic,mod8$waic$waic,mod9$waic$waic),
  lcpo=c(-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T),-mean(log(mod4$cpo$cpo),na.rm=T),-mean(log(mod5$cpo$cpo),na.rm=T),-mean(log(mod6$cpo$cpo),na.rm=T),-mean(log(mod7$cpo$cpo),na.rm=T),-mean(log(mod8$cpo$cpo),na.rm=T),-mean(log(mod9$cpo$cpo),na.rm=T))
)
COMP_iid

write.table(COMP_iid ,"./5_model_vessel_standarization/random_effects/comp_iid.txt" )


