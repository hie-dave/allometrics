# explore and plot dbh-height relationships from the shortlisted Aus-Tallo database

rm(list = ls())

library(dplyr)
library(base)
library(stats)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)

setwd("/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom")

##################################### DEFINE FUNCTIONS

# finds the coefficient of a fitted power function and puts them (along with other variables)
# in a new dataframe
df_power_coefficients <- function(input_df) {
  # names of pfts in input_df
  pftnames<-unique(input_df$pft)
  # number of pfts in input_df
  numpfts<-length(unique(input_df$pft))
  
  cutoff<- .9
  
  # define the new dataframe
  output_df <- data.frame(matrix(NA, nrow = numpfts, ncol = 7))
  colnames(output_df)<-c("pft", "numspecies","maxdbh","maxheight","base","power","r_sqr")

  max_stem_diameter_cm<-0
  max_height_m<-0
  numspeciesVec <- 0
  
  for (ii in 1:numpfts){
    # save the maxmun DBH and height
    max_stem_diameter_cm[ii]<-max(subset(input_df,pft==unique(input_df$pft)[ii])$stem_diameter_cm)
    max_height_m[ii]<-max(subset(input_df,pft==unique(input_df$pft)[ii])$height_m)
    
    # X vector (DBH,cm)
    X<-subset(input_df,pft==unique(input_df$pft)[ii])$stem_diameter_cm
    # Y Vvector (hright,m)
    Y<-subset(input_df,pft==unique(input_df$pft)[ii])$height_m
    powerpft = lm(log(Y) ~ log(X))
    # print(paste0('maxdbh= ',max(subset(input_df,pft==unique(input_df$pft)[ii])$stem_diameter_cm))) # debug
    # print(paste0('maxheight= ',max(subset(input_df,pft==unique(input_df$pft)[ii])$height_m))) # debug
    output_df$pft[ii]<-pftnames[ii]
    output_df$numspecies[ii]<-length(unique(subset(input_df,pft==pftnames[ii])$species))
    output_df$base[ii]<-exp(powerpft$coefficients[1]) # exp(Intercept)
    output_df$power[ii]<-powerpft$coefficients[2]
    output_df$r_sqr[ii]<-summary(powerpft)$r.squared
    # print(output_df[1,])
  }
  # save maximum DBH and height vectors into appropriate place in the output_df dataframe
  output_df$maxdbh<-max_stem_diameter_cm
  output_df$maxheight<-max_height_m
  
  return(output_df)
}

# create DBH-height vectors from the coefficients (power functio) for plotting
create_DBH_height_vectors_from_coefficients_for_plotting <-function(input_df, max_DBH){
  
  # names of pfts in input_df
  pftnames<-unique(input_df$pft)
  # number of pfts in input_df
  numpfts<-length(pftnames)
  # define the X axis (DBH range)
  DBH<-seq(0, max_DBH)
  # define the doutput_df
  output_df <- data.frame(matrix(NA, nrow = numpfts*length(DBH), ncol = 3))
  colnames(output_df)<-c("pred_stem_diamter_cm", "pred_height_m","pft")
  
  # pointer to the dataframe location
  pointer<-1
  
  for (iii in 1:numpfts){
    # # define temporary dataframe 
    # temp_height_pred <- data.frame(matrix(NA, nrow = length(DBH), ncol = 3))
    # colnames(temp_height_pred)<-c("pred_stem_diamter_cm", "pred_height_m","pft")
    # start calculating the predicted tree height
    Base<-input_df$base[iii] # from the input dataframe for  pft (iii)
    Power<-input_df$power[iii] # from the input dataframe for  pft (iii)
    for (i in 1:length(DBH)){
      output_df$pft[pointer+i-1]<-pftnames[iii]
      output_df$pred_stem_diamter_cm[pointer+i-1]<-DBH[i]
      output_df$pred_height_m[pointer+i-1]<-Base*(DBH[i]^Power)
    }
    pointer<-pointer+length(DBH)
    rm(Base)
    rm(Power)
  }
  
  return(output_df)
}

# Prepare a dataframe based on Austallo with family data above a
  # user-defined certain minimum threshold of datapoints
prepare_family_df <- function(input_df, obsthreshold) {
  # find how many data points per family
  familynames<-unique(input_df$family)
  Family_info <- data.frame(name=familynames, datarange=NA)
  obsthreshold<-50 # observation threshold (user decides)
  # create empty database
  subset_input_df_family<--data.frame(matrix(ncol=dim(input_df)[2], nrow=0))
  colnames(subset_input_df_family)<-colnames(input_df)
  for (i in 1:length(familynames)) {
    Family_info$datarange[i]<-length(unique(input_df$stem_diameter_cm[input_df$family==familynames[i]]))
    if (Family_info$datarange[i]>obsthreshold){
      #save data for family[i]
      # subset_input_df_family<-rbind(subset_input_df_family,subset(input_df,family==familynames[i]))
      subset_input_df_family<-rbind(subset_input_df_family,input_df[which(input_df$family == familynames[i]),names(input_df)])
    }
  }
  return(subset_input_df_family)
}

############################ END DEFINE FUNCTIONS

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
pft_coeff_df <- df_power_coefficients(Austallo)

###################################### Plot power function models

height_predicted<-create_DBH_height_vectors_from_coefficients_for_plotting (pft_coeff_df, 200)

cool_rainforest_tree<-subset(Austallo,(pft==pftnames[1]))
floodplain_sclerophyll<-subset(Austallo,(pft==pftnames[2]))
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

h_cool_rainforest_tree<-subset(height_predicted,(pft==pftnames[1]))
h_floodplain_sclerophyll<-subset(height_predicted,(pft==pftnames[2]))
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
sclerophyll_coeff<-rbind(h_tall_sclerophyll_resprout,h_med_sclerophyll_resprout,h_short_sclerophyll_resprout,
                   h_tall_sclerophyll_seed,h_med_sclerophyll_seed,h_short_sclerophyll_seed)

sclerophyll_Austallo<-rbind(tall_sclerophyll_resprout, med_sclerophyll_resprout, short_sclerophyll_resprout,
                   tall_sclerophyll_seed, med_sclerophyll_seed, short_sclerophyll_seed)

# create a df that splits all sclerophyll pfts to "short", "med"  and "tall"
sclerophyll_data<-rbind(tall_sclerophyll_resprout,med_sclerophyll_resprout,short_sclerophyll_resprout,
                   tall_sclerophyll_seed,med_sclerophyll_seed,short_sclerophyll_seed)
sclerophyll_data[sclerophyll_data == "tall_sclerophyll_resprout"] <- "tall_sclerophyll"
sclerophyll_data[sclerophyll_data == "tall_sclerophyll_seed"] <- "tall_sclerophyll"
sclerophyll_data[sclerophyll_data == "med_sclerophyll_resprout"] <- "med_sclerophyll"
sclerophyll_data[sclerophyll_data == "med_sclerophyll_seed"] <- "med_sclerophyll"
sclerophyll_data[sclerophyll_data == "short_sclerophyll_resprout"] <- "short_sclerophyll"
sclerophyll_data[sclerophyll_data == "short_sclerophyll_seed"] <- "short_sclerophyll"

sclerophyll_coeff_df <- df_power_coefficients(sclerophyll_data)
sclerophyll_plot_df <- create_DBH_height_vectors_from_coefficients_for_plotting(sclerophyll_coeff_df,200)

#### FAMILY ANALYSIS
min_obs_threshold <- 50
family_Austallo <- prepare_family_df(Austallo, min_obs_threshold)

### PLOT
# pots all the data faceted by PFT
ggplot(Austallo ,aes(x=log(stem_diameter_cm), y=log(height_m))) +
  geom_point()+
  facet_wrap(pft ~ .)+
  stat_smooth(method = "lm", formula = y ~ x , size = 1)

ggplot(Austallo ,aes(x=(stem_diameter_cm), y=(height_m))) +
  geom_point()+
  facet_wrap(pft ~ ., scales="free")

# plot according to family
ggplot(family_Austallo ,aes(x=(stem_diameter_cm), y=(height_m))) +
  geom_point()+
  facet_wrap(family ~ ., scales="free")+
  ylim(0, 100)
ggplot(family_Austallo ,aes(x=(stem_diameter_cm), y=(height_m))) +
  geom_point()+
  facet_wrap(family ~ .)+
  ylim(0, 100)+
  xlim(0, 200)

ggplot(height_predicted ,aes(x=pred_stem_diamter_cm, y=pred_height_m, col=pft)) +
  geom_line()+
  facet_wrap(pft ~ .)+
  ylim(0, 100)

ggplot(sclerophyll_coeff,aes(x=pred_stem_diamter_cm, y=pred_height_m, col=pft)) +
  geom_line()+
  ylim(0, 60)+
  xlim(0,100)

ggplot(sclerophyll_plot_df,aes(x=pred_stem_diamter_cm, y=pred_height_m, col=pft)) +
  geom_line()+
  ylim(0,100)+
  xlim(0,100)

ggplot(sclerophyll_data,aes(x=(stem_diameter_cm), y=(height_m), col=pft)) +
  geom_point()


ggplot(sclerophyll_Austallo,aes(x=stem_diameter_cm, y=height_m, col=pft)) +
  geom_point()+
  facet_wrap(pft ~ ., scales="free")

# Short (data)- color = species
ggplot(rbind(short_sclerophyll_seed, short_sclerophyll_resprout),aes(x=stem_diameter_cm, y=height_m, col=species)) +
  geom_point()+
  facet_wrap(pft ~ .)

# Medium (data) - color = species
ggplot(rbind(med_sclerophyll_seed, med_sclerophyll_resprout),aes(x=stem_diameter_cm, y=height_m, col=species)) +
  geom_point()+
  facet_wrap(pft ~ .)

# Tall (data) - color = species
ggplot(rbind(tall_sclerophyll_seed, tall_sclerophyll_resprout),aes(x=stem_diameter_cm, y=height_m, col=species)) +
  geom_point()+
  facet_wrap(pft ~ .)

# Short (function)
ggplot(rbind(h_short_sclerophyll_seed, h_short_sclerophyll_resprout),aes(x=pred_stem_diamter_cm, y=pred_height_m, col=pft)) +
  geom_line()

# Medium (function)
ggplot(rbind(h_med_sclerophyll_seed, h_med_sclerophyll_resprout),aes(x=pred_stem_diamter_cm, y=pred_height_m, col=pft)) +
  geom_line()

# Tall (function)
ggplot(rbind(h_tall_sclerophyll_seed, h_tall_sclerophyll_resprout),aes(x=pred_stem_diamter_cm, y=pred_height_m, col=pft)) +
  geom_line()


ggplot(sclerophyll_data,aes(x=log(stem_diameter_cm), y=log(height_m), col=pft)) +
  geom_point()+
  stat_smooth(method = "lm", formula = y ~ x , size = 1)


ggplot(sclerophyll_coeff,aes(x=pred_stem_diamter_cm, y=pred_height_m, col=pft)) +
  geom_line()+
  ylim(0, 60)+
  xlim(0,100)

ggplot(sclerophyll_coeff,aes(x=pred_stem_diamter_cm, y=pred_height_m, col=pft)) +
  geom_line()+
  ylim(0, 60)+
  xlim(0,100)

ggplot(data = species_coeff_df)+ 
  geom_point(aes(x = base, y = power, color = pft, size = numberobserv))

ggplot(data = rbind(subset(species_coeff_df,(pft=="tall_sclerophyll_seed")),
                    subset(species_coeff_df,(pft=="tall_sclerophyll_resprout")),
                    subset(species_coeff_df,(pft=="med_sclerophyll_seed")),
                    subset(species_coeff_df,(pft=="med_sclerophyll_resprout")),
                    subset(species_coeff_df,(pft=="short_sclerophyll_seed")),
                    subset(species_coeff_df,(pft=="short_sclerophyll_resprout"))),
       aes(x = base, y = power, color = pft, size = numberobserv))+ 
  geom_point()+
  stat_ellipse(geom = "polygon",
               aes(fill = pft), 
               alpha = 0.1)

ggplot(data = subset(species_coeff_df,(pft=="tall_sclerophyll_resprout"))) +
  geom_point(aes(x = base, y = power, color = species, size = numberobserv))+
  ylim(0, 1)+
  xlim(0,10)

ggplot(data = subset(species_coeff_df,(pft=="tall_sclerophyll_seed"))) +
  geom_point(aes(x = base, y = power, color = species, size = numberobserv))+
  ylim(0, 1)+
  xlim(0,25)

ggplot(data = rbind(subset(species_coeff_df,(pft=="tall_sclerophyll_seed")),
                    subset(species_coeff_df,(pft=="tall_sclerophyll_resprout")),
                    subset(species_coeff_df,(pft=="med_sclerophyll_seed")),
                    subset(species_coeff_df,(pft=="med_sclerophyll_resprout")),
                    subset(species_coeff_df,(pft=="short_sclerophyll_seed")),
                    subset(species_coeff_df,(pft=="short_sclerophyll_resprout")))) +
  geom_point(aes(x = base, y = power, color = pft, size = numberobserv))+
  ylim(0, 1)+
  xlim(0,18)


# sclerophyll_data - plot a boxplot with 5 cm width (DBH) to compare how the
# three sclerophyl PFTs are comparing
sclerophyll_data_binned <- subset(sclerophyll_data,stem_diameter_cm<60) %>% 
  mutate( bin=cut_width(stem_diameter_cm, width=5, boundary=0) ) %>% 
  ggplot( aes(x=bin, y=height_m) ) +
  geom_boxplot(aes(fill = pft)) +
  xlab("DBH")
sclerophyll_data_binned
# sclerophyll_data - plot a boxplot with 5 cm width (DBH) to compare how the
# three sclerophyl PFTs are comparing
sclerophyll_Austallo_binned <- subset(sclerophyll_Austallo,stem_diameter_cm<60) %>% 
  mutate( bin=cut_width(stem_diameter_cm, width=5, boundary=0) ) %>% 
  ggplot( aes(x=bin, y=height_m) ) +
  geom_boxplot(aes(fill = pft)) +
  xlab("DBH")
sclerophyll_Austallo_binned


### plot exponent and bsae for each PFT with X and Y 
plot_scatter_with_error_bars <- function (datain_df){
  nameofpfts <- unique(datain_df$pft)
  numberofpfts <- length(nameofpfts)
  Base <-0
  Power <-0
  zmean<-0
  xsd <-0
  ysd <-0
  for (i in 1:numberofpfts){
    Base[i]<- mean(datain_df$base[datain_df$pft==nameofpfts[i]], na.rm = TRUE)
    Power[i] <- mean(datain_df$power[datain_df$pft==nameofpfts[i]], na.rm = TRUE)  
    zmean[i] <- mean(datain_df$r_sqr[datain_df$pft==nameofpfts[i]], na.rm = TRUE)  
    xsd [i]<- sd(datain_df$base[datain_df$pft==nameofpfts[i]], na.rm = TRUE)
    ysd[i] <- sd(datain_df$power[datain_df$pft==nameofpfts[i]], na.rm = TRUE)
  }
  print("here")
  MaxX <- max(Base+xsd, na.rm = TRUE)
  MinX <- min(Base-xsd, na.rm = TRUE)
  MaxY <- max(Power+ysd, na.rm = TRUE)
  MinY <- min(Power-ysd, na.rm = TRUE)
  # scatter plot
  # plot(datain_df,Power, pch=16, cex=log(zmean), col=rainbow(length(zmean)), 
  #      ylim=c(0,1),xlim=c(0,10))
  plot(Base, Power, pch=16, col=rainbow(length(zmean)), 
       xlim=c(MinX, MaxX), ylim=c(MinY, MaxY), xlab=("Base"), ylab=("Power"))
  # plot(datain_df, Power, pch=16, cex=(zmean)*2, col=rainbow(length(zmean)), xlim=c(-5,20), ylim=c(.2,.8))
  # x error bars
  segments(Base+xsd,Power,Base-xsd,Power, lwd=2)
  # y error bars
  segments(Base,Power+ysd,Base,Power-ysd, lwd=2)  
  
  # legend("right", legend=nameofpfts, inset = c(0.4, 1.2),
  #        col=rainbow(length(zmean)), pt.cex=1.5, pch=16)
  
}

plot_scatter_with_error_bars(species_coeff_df)


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



# 
# # create pft df with details
# PFTs <- data.frame(name=pftnames, numspecies=numspecies)
# 
















# 
# 
# 
# sclerophyll_resprout<-rbind(tall_sclerophyll_resprout,med_sclerophyll_resprout,short_sclerophyll_resprout)
# 
# ggplot(subset(sclerophyll_resprout,stem_diameter_cm<200) ,aes(x=stem_diameter_cm, y=height_m, col=pft)) +
#   geom_point()
# 
# ggplot(sclerophyll_resprout$stem_diameter_cm ,aes(x=stem_diameter_cm, y=height_m, col=pft)) +
#   geom_point()
# 
# ggplot(subset(tall_sclerophyll_resprout,stem_diameter_cm<200) ,aes(x=stem_diameter_cm, y=height_m, col=species)) +
#   geom_point()
# 
# 
# a<-subset(tall_sclerophyll_seed, species=="Eucalyptus regnans")
# par(mfrow=c(1,2)) # similar to subplot in MATLAB
# plot(a$stem_diameter_cm,a$height_m, pch=0, cex=0.8)
# X<-0:500
# Y<-14.9*X^0.31
# lines(X,Y,col='red')
# 
# 
# 
# # plot for all DBH<200cm
# pft_all_col_200 <- ggplot(subset(Austallo,stem_diameter_cm<200), aes(x=stem_diameter_cm, y=height_m, col=pft)) +
#   geom_point() +
#   facet_wrap(pft ~ .) 
# pft_all_nocol_200 <- ggplot(subset(Austallo,stem_diameter_cm<200), aes(x=stem_diameter_cm, y=height_m)) +
#   geom_point() +
#   facet_wrap(pft ~ .) 
# # plot for all data
# pft_all_col <- ggplot(Austallo, aes(x=stem_diameter_cm, y=height_m, col=pft)) +
#   geom_point() +
#   facet_wrap(pft ~ .) 
# pft_all_nocol <- ggplot(Austallo, aes(x=stem_diameter_cm, y=height_m)) +
#   geom_point() +
#   facet_wrap(pft ~ .) 
# 
# 
# ####
# 
# # plot all species for tall-sclerophyl-resprout in one grqph
# tsr_all_col <- ggplot(subset(Austallo,pft=="tall_sclerophyll_resprout"),
#                       aes(x=stem_diameter_cm, y=height_m, col=species)) +
#   geom_point() 

