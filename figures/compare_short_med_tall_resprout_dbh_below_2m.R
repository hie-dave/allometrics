# explore and plot dbh-height relationships from the shortlisted Aus-Tallo database

rm(list = ls())

library(dplyr)
library(base)
library(stats)
library(ggplot2)

setwd("/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom")

Austallo <- read.csv("/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom/Aus-Tallodf.csv")

# split dataframe according to PFT and set the name of the dataframe according to the PFT
PFTdata<-split(Austallo, Austallo$pft)
pftnames <-unique(Austallo$pft)
numspecies <- rep(0,length(pftnames))
# find number of shortlisted species per PFT
for (i in 1:length(pftnames)) {
  temppftname<-pftnames[i]

  numspecies[i]<-length(unique(Austallo$species[Austallo$pft==temppftname]))
}
# create pft df with details
PFTs <- data.frame(name=pftnames, numspecies=numspecies)
# a
# exponent

# plot for all DBH<200cm
pft_all_col_200 <- ggplot(subset(Austallo,stem_diameter_cm<200), aes(x=stem_diameter_cm, y=height_m, col=pft)) +
  geom_point() +
  facet_wrap(pft ~ .) 
pft_all_nocol_200 <- ggplot(subset(Austallo,stem_diameter_cm<200), aes(x=stem_diameter_cm, y=height_m)) +
  geom_point() +
  facet_wrap(pft ~ .) 
# plot for all data
pft_all_col <- ggplot(Austallo, aes(x=stem_diameter_cm, y=height_m, col=pft)) +
  geom_point() +
  facet_wrap(pft ~ .) 
pft_all_nocol <- ggplot(Austallo, aes(x=stem_diameter_cm, y=height_m)) +
  geom_point() +
  facet_wrap(pft ~ .) 

####

# plot all species for tall-sclerophyl-resprout in one grqph
tsr_all_col <- ggplot(subset(Austallo,pft=="tall_sclerophyll_resprout"),
                      aes(x=stem_diameter_cm, y=height_m, col=species)) +
  geom_point() 

a<-subset(Austallo,(pft=="tall_sclerophyll_resprout"))
b<-subset(Austallo,(pft=="med_sclerophyll_resprout"))
c<-subset(Austallo,(pft=="short_sclerophyll_resprout"))
d<-rbind(a,b,c)
          
ggplot(subset(d,stem_diameter_cm<200) ,aes(x=stem_diameter_cm, y=height_m, col=pft)) +
  geom_point()



