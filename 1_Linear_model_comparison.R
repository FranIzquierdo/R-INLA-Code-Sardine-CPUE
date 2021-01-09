rm(list=ls()) # Clean workspace
library(INLA)

# Read dataset ------------------------------------------------------------
rm(list=ls()) # Clean workspace
data <- read.table("./data/working_data/2_data_obs_joined_5_covs_scaled.txt", header=TRUE, dec=".")

# Loop function ------------------------------------------------------------

## This loop function allows to run and compare multiple bayesian linear models
## Note all possible linear combination between variables will be tried
## Previous exploratory analyses to avoid correlation and collinearity issues must be done

Bdiclcpomodel_stack<-function(resp, variables, datos, n, family="Gamma",...)
{  
  #Terms
  sel.terms <- switch('terms',terms=variables)
  
  # All possible combinations from the m elements of v
  comb.terms <- function(m, v=sel.terms) {
    if(m==0) return('resp ~ -1')
    else {
      combis <- apply(combn(v, m), 2, paste, collapse=' + ')
      return(paste('resp ~ -1', combis, sep=' + '))
    }
  }
  
  #Creates a list with all possible models
  f.list <- unlist(sapply(0:length(sel.terms), comb.terms))
  
  # Run all models and save them in 'f.list'. Then save DIC,LCPO and WAIC
  dic<-numeric()
  LCPO<-numeric()
  waic<-numeric()
  for(i in 1:length(f.list)){
    res =inla(formula = eval(parse(text=f.list[i])), family=family, data=datos, ...)
    dic[i] <- res$dic$dic
    LCPO[i] = -mean(log(res$cpo$cpo))
    waic[i]<-res$waic$waic
    print(c(i, dic[i], waic[i], LCPO[i]))
  }
  
  # Show all models ordered by DIC
  modelos_dic<-data.frame(f.list[order(dic)[1:n]], dic[order(dic)[1:n]], waic[order(dic)[1:n]], LCPO[order(dic)[1:n]])
  colnames(modelos_dic)<-c("Modelos", "Dic", "Waic", "LCPO")
  
  #Show all models ordered by WAIC
  modelos_waic<-data.frame(f.list[order(waic)[1:n]], dic[order(waic)[1:n]], waic[order(waic)[1:n]], LCPO[order(waic)[1:n]])
  colnames(modelos_waic)<-c("Modelos", "Dic", "Waic", "LCPO")
  
  
  #Show all models ordered by LCPO
  modelos_lcpo<-data.frame(f.list[order(LCPO)[1:n]], dic[order(LCPO)[1:n]], waic[order(LCPO)[1:n]], LCPO[order(LCPO)[1:n]])
  colnames(modelos_lcpo)<-c("Modelos", "Dic", "Waic", "LCPO")
  
  
  modelos<-list(modelos_dic, modelos_waic, modelos_lcpo)
  names(modelos)<-c("Modelos dic", "Modelos waic", "Modelos lcpo")
  modelos
  
}

#1) EFFORT variable ---------------------------------------------------

## We try among the available effort variables on the sardine's dataset
## Distance to harbour and total catch were discarded
## Note that the 3 variables were correlated, so just models with 1 variable will be selected
## Note also that all variables were scaled before the modeling process

variables <- c("power_esc","length_esc", "tonn_esc")#f indica el efecto espacial, model=spde y si fuese aleatorio seria model="iid"

### --- Response variable --- ###
resp=data$sardine

### --- Call the function --- ###
models_gam<-Bdiclcpomodel_stack(resp=resp, variables=variables, datos=data, n=15,
                                family="gamma",
                                control.predictor=list(compute=TRUE),
                                control.compute = list(config=TRUE, dic=TRUE, cpo=TRUE, waic=TRUE),
                                num.threads=3,
                                control.inla=list(strategy="gaussian"),
                                verbose=FALSE)

models_gam #Select the best model which has "Length of the vessel" as the best variable
saveRDS(models_gam, "./5_model_vessel_standarization/models_fixed_var_selection.rds")


#2) TRY environmental variable  ---------------------------------------------------

# Just have an idea, once we have selected length, we try the environmental covs
# Note this models contains just linear effects
# Note all variables were previously scaled to avoid numerical confounding

variables <- c("length_esc","chl_mean_esc","bath_esc","sst_mean_esc", "int_mean_esc","ang_mean_esc")#f indica el efecto espacial, model=spde y si fuese aleatorio seria model="iid"

### --- Response variable --- ###
resp=data$sardine

### --- Call the function --- ###
models_gam<-Bdiclcpomodel_stack(resp=resp, variables=variables, datos=data, n=25,
                                family="gamma",
                                control.predictor=list(compute=TRUE),
                                control.compute = list(config=TRUE, dic=TRUE, cpo=TRUE, waic=TRUE),
                                num.threads=3,
                                control.inla=list(strategy="gaussian"),
                                verbose=FALSE)

models_gam
saveRDS(models_gam, "./5_model_vessel_standarization/models_fixed_env_var_selection.rds")


