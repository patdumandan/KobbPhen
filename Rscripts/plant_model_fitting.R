library(cmdstanr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(bayesplot)
library(posterior)
library(ggridges)

plant_mod=cmdstan_model("C:/pdumandanSLU/PatD-SLU/SLU/phenology-project/KobbPhen/plant_model.stan")

plant_fits= lapply(plant_list, fit_plant_model,data= kobb_plant_dat, model= plant_mod)

preds_plant=lapply(plant_list, plot_plant_preds, data=kobb_plant_dat)

