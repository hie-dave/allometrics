# Exampine_wide_species

rm(list = ls())

library(dplyr)
library(base)
library(stats)
library(ggplot2)


setwd("/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom")

# what is the DBH limit that I am looking at
DBH_lim <-200

Austallo <- read.csv("/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom/Aus-Tallodf.csv")
# limit to below 2m DBH
Austallo_lim200 <- subset(Austallo,stem_diameter_cm <= DBH_lim) 
# create a vector with species names that their DBH can be higher than 2m
Wide_Species_list<-unique(subset(Austallo,stem_diameter_cm > DBH_lim)$species)
# a df that contains ONLY the 
Wide_Species<-Austallo[Austallo$species %in% Wide_Species_list,]
# create a vector with the lits of PFTs that the above species belong to
Wide_Species_pft_list<-unique(Wide_Species$pft)


plot(Wide_Species$stem_diameter_cm,Wide_Species$height_m)

# FACTES : PFTs, In each panel, species per PFT
ggplot(Wide_Species ,aes(x=stem_diameter_cm, y=height_m, col=species)) +
  geom_point()+
  facet_wrap(pft ~ .) 
# FACTES : species, In each panel, species by itself
ggplot(Wide_Species ,aes(x=stem_diameter_cm, y=height_m, col=species)) +
  geom_point()+
  facet_wrap(species ~ .) 

