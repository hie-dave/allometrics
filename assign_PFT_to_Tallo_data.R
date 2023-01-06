# This code: 1) creates a new df with species and pfts from Clare's csv file; and 
# 2) creates a new df with allometric data for the species found in the df created in (1)
# 3) creates a vector with species names with unmatched allometric data.

# NEXT STEP
# remove outlayers and create allometric relationships per PFT and then extract parameters for LPJ-GUESS

rm(list = ls())
library(dplyr)
library(base)
library(stats)

setwd("/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom")

# Load Clar's file as df
CSsp2pftdf <- read.csv(file='PFT_classes_taxon_v1.csv')
# import Tallo file as dataframe
Tallodf <- read.csv(file='Tallo.csv')

## CREATE LIST OF SPECIES ASSOCIATED WITH PFTs FROM CLARE'S CSV
col_index <- 1:ncol(CSsp2pftdf)
selected_columns <- col_index[col_index %% 2 != 0] # this is modulus: col_index %% 2
PFTs <- names(CSsp2pftdf[selected_columns]) # vector with the name of the PFTs
# Create an empty data frame with 3 columns
sppftdf <- data.frame(NA, NA, NA)
colnames(sppftdf) <- c('species','count','pft')
# get only the 2 first column
n <- 1:length(PFTs)
for (i in n) {
  temppft <- data.frame(CSsp2pftdf[,selected_columns[i]:(selected_columns[i]+1)])
  temppft$pft <- PFTs[i]
  colnames(temppft) <- c('species','count','pft')
  sppftdf <- rbind(sppftdf, temppft)
}
# remove nans according to first column (species)
# this df is al list of 
sppftdf <- sppftdf[complete.cases(sppftdf[,1]), ] # Omit NAs by columns

# Create new dataframe for the australian species
namesnewdf <- names(Tallodf)
newTallodf <- data.frame(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
colnames(newTallodf) <- c(namesnewdf, "count", "pft") 
# Create new vector for the species that has no dat ain the Tallo dataset
nonfoundvec <- c()

# go through the pft df and create a new Tallo df with australian species
for (j in 1:nrow(sppftdf)) {
  
  CSspecies <- sppftdf$species[j]
  indices <- which(Tallodf$species == CSspecies)
  
  cond <- (is.integer(indices) && length(indices) == 0) # if the result is integer(0), then cond = TRUE
  if (cond == FALSE) {
    # in case there are some results found in Tallo for this species
    print("found") # only for testing
    # Add the correct line (with the index found in indices)
    tempTallo<-Tallodf[indices,]
    tempTallo$count<-sppftdf$count[j]
    tempTallo$pft<- sppftdf$pft[j]
    newTallodf <- rbind(newTallodf, tempTallo )
  } else {
    # in case there are no results found in Tallo for this species
    print("not found") # only for testing
    # Make a list of unfound species (from pftdf) to check manually
    nonfoundvec <- append(nonfoundvec, CSspecies)
    
  }
  
}
newTallodf <- newTallodf[complete.cases(newTallodf[,1]), ] # Omit NAs by columns
# order alphabetically acc
newTallodf <-newTallodf[order(newTallodf$pft, newTallodf$species),]
nonfoundvec <- sort(nonfoundvec)
# save csv
write.csv(newTallodf,"/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom/newTallodf.csv")
write.csv(nonfoundvec,"/Users/30060406/Dropbox/R_CODE/HIE/PFT/allom/nonfoundvec.csv")
