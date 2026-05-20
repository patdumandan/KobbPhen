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

plant_fits= lapply(plant_list, fit_plant_model,data= kobb_plant_dat, model= plant_mod)
plant_fits2= lapply(plant_list, fit_plant_model,data= kobb_plant_dat, model= plant_mod2)

#predictions
preds_plant=lapply(plant_list, plot_plant_preds, data=kobb_plant_dat)

#diagnistics
peak_diag <- diagnose_peak_parameters(species_name = "Angelica_archangelica",
                                      data = kobb_plant_dat)

peak_diag=lapply(plant_list, diagnose_peak_parameters, kobb_plant_dat)
