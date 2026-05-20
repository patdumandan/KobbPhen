require(tidyr)
require(dplyr)
require(lubridate)

#data####
file_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\KobbPhen"
plant_file=paste(file_path, '\\Kobb_plant_raw','.csv', sep='')
kobb_plant=read.csv(plant_file, header=T,  stringsAsFactors = F)

#restructure

kobb_plant$date=as.Date.character(kobb_plant$Date)
kobb_plant_dat=kobb_plant%>%
  mutate(year=year(date), month=month(date), dia=day(date),
         DOY=yday(date))%>%
  pivot_longer(cols=13:27, names_to = "taxon")%>%
  rename(abundance=value, PlotID=new.overall.plot.name)%>%
  mutate(plot_id = dense_rank(PlotID),
         taxon=ifelse(as.character(taxon) == "Taraxacum_sp.", "Taraxacum_sp", as.character(taxon)))%>% #to create species-specific plot ID
  select(-Transect_and_corresponding_malaise_trap_name, -Biobasis_Plot_name_that_is_best_representative)

kobb_plant_dat$abundance[is.na(kobb_plant_dat$abundance)] <- 0
kobb_plant_dat$DOYs = scale(kobb_plant_dat$DOY, center = TRUE, scale = TRUE)[,1]

plant_list=unique(kobb_plant_dat$taxon)

