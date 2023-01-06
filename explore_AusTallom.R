# explore and plot dbh-height relationships from the shortlisted Aus-Tallo database

rm(list = ls())

library(dplyr)
library(base)
library(stats)
library(ggplot2)

setwd("/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom")


### DEFINE FUNCTIONS
df_power_coefficients <- function(input_df, minimum_datapoints) {
  output_df <- data.frame(pft=pftnames,
                          numspecies=numspecies,
                          maxdbh=NA,
                          maxheight=NA,
                          base=NA,
                          power=NA,
                          r_sqr=NA)
  
  numpfts<-length(unique(input_df$pft))
  max_stem_diameter_cm<-0
  max_height_m<-0
  
  for (ii in 1:numpfts){
    # save the maxmun DBH and height
    max_stem_diameter_cm[ii]<-max(subset(input_df,pft==unique(input_df$pft)[ii])$stem_diameter_cm)
    max_height_m[ii]<-max(subset(input_df,pft==unique(input_df$pft)[ii])$height_m)
    # X vector (DBH,cm)
    X<-subset(input_df,pft==unique(input_df$pft)[ii])$stem_diameter_cm
    # Y Vvector (hright,m)
    Y<-subset(input_df,pft==unique(input_df$pft)[ii])$height_m
    powerpft = lm(log(Y) ~ log(X))
    print(paste0('maxdbh= ',max(subset(input_df,pft==unique(input_df$pft)[ii])$stem_diameter_cm))) # debug
    print(paste0('maxheight= ',max(subset(input_df,pft==unique(input_df$pft)[ii])$height_m))) # debug
    output_df$base[ii]<-exp(powerpft$coefficients[1]) # exp(Intercept)
    output_df$power[ii]<-powerpft$coefficients[2]
    output_df$r_sqr[ii]<-summary(powerpft)$r.squared
  }
  # save maximum DBH and height vectors into appropriate place in the output_df dataframe
  output_df$maxdbh<-max_stem_diameter_cm
  output_df$maxheight<-max_height_m
  
  return(output_df)
}


Austallo <- read.csv("/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom/Aus-Tallodf.csv")
# limit to below 2m DBH
Austallo_lim200 <- subset(Austallo,stem_diameter_cm <= 200) 

# split dataframe according to PFT and set the name of the dataframe according to the PFT
PFTdata<-split(Austallo, Austallo$pft)
pftnames <-unique(Austallo$pft)
specieslist <- unique(Austallo$species) # create species list
numspeciesVec <- rep(0,length(pftnames))
numspecies <- sum(numspeciesVec)

# find number of shortlisted species per PFT
for (i in 1:length(pftnames)) {
  temppftname<-pftnames[i]
  numspecies[i]<-length(unique(Austallo$species[Austallo$pft==temppftname]))
}

# CREATE df WITH COEFFICIENTS OF THE POWER LAW PER SPECIES AND PFT
species_coeff_df <- data.frame(species=specieslist,
                               pft=NA,
                               numberobserv=NA,
                               base=NA,
                               power=NA,
                               r_sqr=NA)
# minimum number of datapoints belowhich species is ignored
minimum_datapoints <- 20
for (j in 1:sum(numspecies)){
  # fill the $pft column
  speciestemp<-subset(Austallo, species==species_coeff_df$species[j])
  species_coeff_df$numberobserv[j]<-dim(speciestemp)[1]
  # fill details only if there is a minimum umber of data points to form a curve
  if (dim(speciestemp)[1]>minimum_datapoints){
      species_coeff_df$pft[j]<-speciestemp$pft[1]
      ## Fill coefficients of a power function
      y<-speciestemp$height_m
      x<-speciestemp$stem_diameter_cm
      powerspecies = lm(log(y) ~ log(x))
      species_coeff_df$base[j]<-exp(powerspecies$coefficients[1]) # exp(Intercept)
      species_coeff_df$power[j]<-powerspecies$coefficients[2]
      species_coeff_df$r_sqr[j]<-summary(powerspecies)$r.squared
  }
}
# get rid of lines with NA
species_coeff_noNA<-na.omit(species_coeff_df)
num_species_removed<-dim(species_coeff_df)[1]-dim(species_coeff_noNA)[1]
species_coeff_df<-species_coeff_noNA

###################################### K_allom2 and K_allom3
# Find paremters (coefficients) for the power function per PFT
pft_coeff_df <- df_power_coefficients(Austallo, 20)
# 
# pft_coeff_df <- data.frame(pft=pftnames, 
#                            numspecies=numspecies,
#                            maxdbh=NA,
#                            maxheight=NA,
#                            base=NA,
#                            power=NA,
#                            r_sqr=NA)
# max_stem_diameter_cm<-0
# max_height_m<-0
# for (ii in 1:length(pftnames)){
#   # save the maxmun DBH and height
#   max_stem_diameter_cm[ii]<-max(subset(Austallo,pft==unique(Austallo$pft)[ii])$stem_diameter_cm)
#   max_height_m[ii]<-max(subset(Austallo,pft==unique(Austallo$pft)[ii])$height_m)
#   # X vector (DBH,cm)
#   X<-subset(Austallo,pft==unique(Austallo$pft)[ii])$stem_diameter_cm
#   # Y Vvector (hright,m)
#   Y<-subset(Austallo,pft==unique(Austallo$pft)[ii])$height_m
#   powerpft = lm(log(Y) ~ log(X))
#   # print(paste0('maxdbh= ',max(subset(Austallo,pft==unique(Austallo$pft)[ii])$stem_diameter_cm))) % debug
#   # print(paste0('maxheight= ',max(subset(Austallo,pft==unique(Austallo$pft)[ii])$height_m))) % debug
#   pft_coeff_df$base[ii]<-exp(powerpft$coefficients[1]) # exp(Intercept)
#   pft_coeff_df$power[ii]<-powerpft$coefficients[2]
#   pft_coeff_df$r_sqr[ii]<-summary(powerpft)$r.squared
# }
# # save maximum DBH and height vectors into appropriate place in the pft_coeff_df dataframe
# pft_coeff_df$maxdbh<-max_stem_diameter_cm
# pft_coeff_df$maxheight<-max_height_m

###################################### Plot power function models
DBH<-seq(0, 200)
height_predicted <- data.frame(pred_stem_diamter_cm=NA,
                               pred_height_m=NA,
                               pft=NA)
for (iii in 1:length(pftnames)){
  temp_height_pred <- data.frame(pred_stem_diamter_cm=DBH,
                                 pred_height_m=NA,
                                 pft=NA)
  Base<-pft_coeff_df$base[iii]
  Power<-pft_coeff_df$power[iii]
  # insert PFT name
  temp_height_pred$pft<-pftnames[iii]
  # insert Y data (predicted height) per PFT
  temp_height_pred$pred_height_m<-Base*((temp_height_pred$pred_stem_diamter_cm)^Power)
  # print('here') # debug line
  height_predicted<-rbind(height_predicted,temp_height_pred)
}
# remove first line from df {is an NA}
height_predicted = height_predicted[-1,]

floodplain_sclerophyll<-subset(Austallo,(pft==pftnames[1]))
cool_rainforest_tree<-subset(Austallo,(pft==pftnames[2]))
med_sclerophyll_resprout<-subset(Austallo,(pft==pftnames[3]))
med_sclerophyll_seed<-subset(Austallo,(pft==pftnames[4]))
mesic_shrub<-subset(Austallo,(pft==pftnames[5]))
Nfix_mesic_midstory<-subset(Austallo,(pft==pftnames[6]))
Nfix_mod_rain<-subset(Austallo,(pft==pftnames[7]))
Nfix_xeric<-subset(Austallo,(pft==pftnames[8]))
savanna_tree<-subset(Austallo,(pft==pftnames[9]))
short_sclerophyll_resprout<-subset(Austallo,(pft==pftnames[10]))
short_sclerophyll_seed<-subset(Austallo,(pft==pftnames[11]))
subalpine_tree<-subset(Austallo,(pft==pftnames[12]))
tall_sclerophyll_resprout<-subset(Austallo,(pft==pftnames[13]))
tall_sclerophyll_seed<-subset(Austallo,(pft==pftnames[14]))
warm_rainforest_tree<-subset(Austallo,(pft==pftnames[15]))
xeric_shrub<-subset(Austallo,(pft==pftnames[16]))
xeric_shrub_Nfix<-subset(Austallo,(pft==pftnames[17]))

h_floodplain_sclerophyll<-subset(height_predicted,(pft==pftnames[1]))
h_cool_rainforest_tree<-subset(height_predicted,(pft==pftnames[2]))
h_med_sclerophyll_resprout<-subset(height_predicted,(pft==pftnames[3]))
h_med_sclerophyll_seed<-subset(height_predicted,(pft==pftnames[4]))
h_mesic_shrub<-subset(height_predicted,(pft==pftnames[5]))
h_Nfix_mesic_midstory<-subset(height_predicted,(pft==pftnames[6]))
h_Nfix_mod_rain<-subset(height_predicted,(pft==pftnames[7]))
h_Nfix_xeric<-subset(height_predicted,(pft==pftnames[8]))
h_savanna_tree<-subset(height_predicted,(pft==pftnames[9]))
h_short_sclerophyll_resprout<-subset(height_predicted,(pft==pftnames[10]))
h_short_sclerophyll_seed<-subset(height_predicted,(pft==pftnames[11]))
h_subalpine_tree<-subset(height_predicted,(pft==pftnames[12]))
h_tall_sclerophyll_resprout<-subset(height_predicted,(pft==pftnames[13]))
h_tall_sclerophyll_seed<-subset(height_predicted,(pft==pftnames[14]))
h_warm_rainforest_tree<-subset(height_predicted,(pft==pftnames[15]))
h_xeric_shrub<-subset(height_predicted,(pft==pftnames[16]))
h_xeric_shrub_Nfix<-subset(height_predicted,(pft==pftnames[17]))

# combine the sclerophylls together
sclerophyll<-rbind(h_tall_sclerophyll_resprout,h_med_sclerophyll_resprout,h_short_sclerophyll_resprout,
                   h_tall_sclerophyll_seed,h_med_sclerophyll_seed,h_short_sclerophyll_seed)

# create a df that splits all sclerophyll pfts to "short", "med"  and "tall"
sclerophyll_data<-rbind(tall_sclerophyll_resprout,med_sclerophyll_resprout,short_sclerophyll_resprout,
                   tall_sclerophyll_seed,med_sclerophyll_seed,short_sclerophyll_seed)
sclerophyll_data[sclerophyll_data == "tall_sclerophyll_resprout"] <- "tall_sclerophyll"
sclerophyll_data[sclerophyll_data == "tall_sclerophyll_seed"] <- "tall_sclerophyll"
sclerophyll_data[sclerophyll_data == "med_sclerophyll_resprout"] <- "med_sclerophyll"
sclerophyll_data[sclerophyll_data == "med_sclerophyll_seed"] <- "med_sclerophyll"
sclerophyll_data[sclerophyll_data == "short_sclerophyll_resprout"] <- "short_sclerophyll"
sclerophyll_data[sclerophyll_data == "short_sclerophyll_seed"] <- "short_sclerophyll"

### PLOT
# pots all the data faceted by PFT
ggplot(Austallo ,aes(x=log(stem_diameter_cm), y=log(height_m))) +
  geom_point()+
  facet_wrap(pft ~ .)+
  stat_smooth(method = "lm", formula = y ~ x , size = 1)

ggplot(Austallo ,aes(x=(stem_diameter_cm), y=(height_m))) +
  geom_point()+
  facet_wrap(pft ~ .)

ggplot(height_predicted ,aes(x=pred_stem_diamter_cm, y=pred_height_m, col=pft)) +
  geom_line()+
  facet_wrap(pft ~ .)+
  ylim(0, 100)

ggplot(sclerophyll,aes(x=pred_stem_diamter_cm, y=pred_height_m, col=pft)) +
  geom_line()+
  ylim(0, 60)+
  xlim(0,100)

ggplot(sclerophyll_data,aes(x=(stem_diameter_cm), y=(height_m), col=pft)) +
  geom_point()

ggplot(sclerophyll_data,aes(x=log(stem_diameter_cm), y=log(height_m), col=pft)) +
  geom_point()+
  stat_smooth(method = "lm", formula = y ~ x , size = 1)


ggplot(data = species_coeff_df) + geom_point(aes(x = base, y = power, color = pft, size = numberobserv))
ggplot(data = species_coeff_df) + geom_point(aes(x = base, y = power, color = pft, size = r_sqr))


par(mfrow=c(1,5)) 
hist(species_coeff_df$power, xlab="exponent")
hist(species_coeff_df$base, xlab="base")
plot(species_coeff_df$power, species_coeff_df$r_sqr, ylab="R^2", xlab="exponent")
plot(species_coeff_df$base, species_coeff_df$r_sqr, ylab="R^2", xlab="base")
plot(species_coeff_df$base, species_coeff_df$power, ylab="power", xlab="base")


par(mfrow=c(1,3)) 
plot(species_coeff_df$numberobserv, species_coeff_df$r_sqr, ylab="R^2", xlab="number of observations")
plot(species_coeff_df$numberobserv, species_coeff_df$base, ylab="base", xlab="number of observations")
plot(species_coeff_df$numberobserv, species_coeff_df$power, ylab="power", xlab="number of observations")
########

par(mfrow=c(,2)) # similar to subplot in MATLAB
plot(xeric_shrub_Nfix$stem_diameter_cm,xeric_shrub_Nfix$height_m, pch=0, cex=0.8)
xa<-0:100
xb<-pft_coeff_df$base[17]*xa^pft_coeff_df$power[17]
lines(xa,xb,col='red')




# create pft df with details
PFTs <- data.frame(name=pftnames, numspecies=numspecies)




















sclerophyll_resprout<-rbind(tall_sclerophyll_resprout,med_sclerophyll_resprout,short_sclerophyll_resprout)

ggplot(subset(sclerophyll_resprout,stem_diameter_cm<200) ,aes(x=stem_diameter_cm, y=height_m, col=pft)) +
  geom_point()

ggplot(sclerophyll_resprout$stem_diameter_cm ,aes(x=stem_diameter_cm, y=height_m, col=pft)) +
  geom_point()

ggplot(subset(tall_sclerophyll_resprout,stem_diameter_cm<200) ,aes(x=stem_diameter_cm, y=height_m, col=species)) +
  geom_point()


a<-subset(tall_sclerophyll_seed, species=="Eucalyptus regnans")
par(mfrow=c(1,2)) # similar to subplot in MATLAB
plot(a$stem_diameter_cm,a$height_m, pch=0, cex=0.8)
X<-0:500
Y<-14.9*X^0.31
lines(X,Y,col='red')



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

