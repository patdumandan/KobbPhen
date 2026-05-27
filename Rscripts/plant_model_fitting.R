library(cmdstanr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(bayesplot)
library(posterior)
library(ggridges)

plant_mod=cmdstan_model("C:/pdumandanSLU/PatD-SLU/SLU/phenology-project/KobbPhen/RScripts/plant_model.stan")
plant_mod2=cmdstan_model("C:/pdumandanSLU/PatD-SLU/SLU/phenology-project/KobbPhen/RScripts/plant_mod2.stan")
plant_mod3=cmdstan_model("C:/pdumandanSLU/PatD-SLU/SLU/phenology-project/KobbPhen/RScripts/plant_mod3.stan")

plant_fits= lapply(plant_list, fit_plant_model,data= kobb_dat2, model= plant_mod)
plant_fits2= lapply(plant_list, fit_plant_model,data= kobb_dat2, model= plant_mod2)
plant_fits3= lapply(plant_list, fit_plant_model,data= kobb_dat2, model= plant_mod3)

#predictions
preds_plant1=lapply(plant_list, plot_plant_preds1, data=kobb_dat2)

#diagnistics
peak_diag=lapply(plant_list, diagnose_peak_parameters, kobb_plant_dat)
