library(tidyverse)
library(ggplot2)
library(dplyr)


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




