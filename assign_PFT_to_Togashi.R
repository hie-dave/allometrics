# Togashi et al. Moisture availability constraints on LA:SA dataset
# Units => LA:SA (log10 transformed); Height (m); Alpha (annual unitless); MI (annual unitless)
# Zanne et al. dataset is not included. The dataset is currenlty submitted (Zanne et al.) and it can be requested to the main author (Henrique Furstenau Togashi).
# 
# As a condition for use of the dataset, we request that users agree 
# (1) To notify the main author (Henrique Furstenau Togashi) if the dataset is to be used in any publication;
# (2) To recognise that the researchers who gathered these data with formal recognition that may include co-authorship or acknowledgements on publications

# aim: creates a df that assigns the SA:LA database (Togashi)
# and PFTs to species that exist in Clare's/Laura's lists.

# This code: 1) creates a new df with species and pfts from Laura's csv file; and 
# 2) creates a new df with allometric data for the species found in the df created in (1)
# 3) creates a vector with species names with unmatched allometric data.

# note: the code also changes species name to be made of up to 2 words (to avoid
# species data mssing because species name is not exactly equal)

# outouts saved as files:
# 1) csv file with fata for speceis from Togashi's (latosa) database from species that exist in 
# both Togashi's (latosa) and Clare's/Laura's list
# 2) csv file with species that could not be found in the Togashi's (latosa) database (but exists
# in the other one)

# Last Update: 08/12/2022, Assaf

rm(list = ls())
library(dplyr)
library(base)
library(stats)
library(stringr) 
library(tidyverse) # new

setwd("/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom")

# new - efine a function that converts species names that are larged than 2 words to 
# species names with 2 words.
extract_n_words_from_string <- function(species_name_list){
  library(stringr)
  # species_name<-c('assaf ibbar the 2nd')
  twoWordsSpecies_name_list<-word(species_name_list, 1,2, sep=" ")
  return(twoWordsSpecies_name_list)
}

taxonfilename <- "/Users/30060406/Dropbox/R_CODE/HIE/CFA/species_to_pfts/PFT_classes_taxon_v3_20221114.csv"

# Load Clar's file as df
CSsp2pftdf <- read.csv(file=taxonfilename)
# import Tallo file as dataframe
Togashi <- read.csv(file='Togashi.csv')

## CREATE LIST OF SPECIES ASSOCIATED WITH PFTs FROM Laura'sS CSV
columnnames <- names(CSsp2pftdf)
selected_columns <- seq(1, ncol(CSsp2pftdf), by=3) 
PFTs <- names(CSsp2pftdf[selected_columns]) # vector with the name of the PFTs

# Create an empty data frame with 3 columns
sppftdf <- data.frame(NA, NA, NA)
colnames(sppftdf) <- c('species','count','pft')
# creates a df from the taxon file (Clare/Laura) that aggregates all the 
# species+count+pfts in rows...
n <- 1:length(PFTs)
for (i in n) {
  temppft <- data.frame(CSsp2pftdf[,c(selected_columns[i],(selected_columns[i]+2))])
  temppft$pft <- PFTs[i]
  colnames(temppft) <- c('species','count','pft')
  sppftdf <- rbind(sppftdf, temppft)
}
# new - shorten species name to n=2 words only
sppftdf$species<-extract_n_words_from_string(sppftdf$species)
# new - remove duplicates according to species
sppftdf<-sppftdf[!duplicated(sppftdf$species), ]
# new - remove nans according to first column (species)
sppftdf <- sppftdf[complete.cases(sppftdf[,1]), ] # Omit NAs by columns

# new - shorten species name to n=2 words only (Togashi's (latosa) database)
# ps - no need to remove replicates from this df because each line is data...
Togashi$Specie<-extract_n_words_from_string(Togashi$Specie)
# for Togashi, remove all the lines that were measured from branches
Togashi<-subset(Togashi, Sample=="trunk")
# Create new dataframe for the Australian species
namesnewdf <- names(Togashi)
# newTogashidf <- data.frame(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
newTogashidf <- data.frame(matrix(nrow = 0, ncol = length(names(Togashi))+2))
colnames(newTogashidf) <- c(namesnewdf, "count", "pft") 
# Create new vector for the species that has no data in the Tallo dataset
nonfoundvec <- c()

# go through the pft df and create a new Tallo df with australian species
for (j in 1:nrow(sppftdf)) {
  
  # get one species (j) from the PFt-species list
  CSspecies <- sppftdf$species[j]
  # make vector of all the returns of this species in the Togashi's (latosa) database
  indices <- which(Togashi$Specie == CSspecies)
  
  cond <- (is.integer(indices) && length(indices) == 0) # if the result is integer(0), then cond = TRUE
  if (cond == FALSE) {
    # in case there are some results found in Tallo for this species
    # print("found") # only for testing
    # Add the correct line (with the index found in indices)
    tempTogashi<-Togashi[indices,]
    tempTogashi$count<-sppftdf$count[j]
    tempTogashi$pft<- sppftdf$pft[j]
    newTogashidf <- rbind(newTogashidf, tempTogashi )
  } else {
    # in case there are no results found in Tallo for this species
    # print("not found") # only for testing
    # Make a list of unfound species (from pftdf) to check manually
    nonfoundvec <- append(nonfoundvec, CSspecies)
    
  }
  
}
newTogashidf <- newTogashidf[complete.cases(newTogashidf[,1]), ] # Omit NAs by columns
# order alphabetically acc
newTogashidf <-newTogashidf[order(newTogashidf$pft, newTogashidf$Specie),]
nonfoundvecTogashi <- sort(nonfoundvec)
# save csv
write.csv(newTogashidf,"/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom/Aus-Togashi.csv")
write.csv(nonfoundvecTogashi,"/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom/Aus-Togashinonfoundvec.csv")




