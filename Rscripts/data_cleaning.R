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

#data filtering####

#remove spp. with <1000 individuals total

kobb_plant_dat_tot=kobb_plant_dat%>%
  group_by(taxon)%>%summarise(tot=sum(abundance))%>%
  filter(tot< 1000)

kobb_list1=c(unique(kobb_plant_dat_tot$taxon))

kobb_dat1=kobb_plant_dat%>%filter(!taxon%in%kobb_list1)

#remove spp with <50+ indivs per year

kobb_dat1_sum=kobb_dat1%>%group_by(taxon, year)%>%
  summarise(tot=sum(abundance))%>%
  mutate(include=if_else(tot<50, "no", "yes"))

kobb_dat2=kobb_plant_dat%>%left_join(kobb_dat1_sum)%>%
  filter(!include=="no")

plant_list=unique(kobb_dat2$taxon)
