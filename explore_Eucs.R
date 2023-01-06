# create a subset of Austallo, but for Eucs
Eucs <- subset(Austallo,genus=="Eucalyptus")

speciesnames<-unique(Eucs$species)
numofspecieseucs <- length(unique(Eucs$species))
# find how many data points per family

# createdata frame that contains allometric data in the species level with data
 # above a certain number of minimum observatiokns
Euc_species_info <- data.frame(name=speciesnames, datarange=NA)
obsthreshold<-20 # observation threshold (user decides)
subset_input_df_euc<-data.frame(matrix(ncol=dim(Eucs)[2], nrow=0))
colnames(subset_input_df_euc)<-colnames(Eucs)
for (i in 1:numofspecieseucs) {
  Euc_species_info$datarange[i]<-length(unique(Eucs$stem_diameter_cm[Eucs$species==speciesnames[i]]))
  if (Euc_species_info$datarange[i]>obsthreshold){
    #save data for family[i]
    # subset_input_df_family<-rbind(subset_input_df_family,subset(input_df,family==speciesnames[i]))
    subset_input_df_euc<-rbind(subset_input_df_euc,Eucs[which(Eucs$species == speciesnames[i]),names(Eucs)])
  }
}

# change the order of the dataframe so it will plot according to pfts (same color)
eucpftnames<-unique(subset_input_df_euc$pft)
subset_input_df_euc$pft_order <- factor(subset_input_df_euc$pft, levels=c("short_sclerophyll_resprout","short_sclerophyll_seed",
                                                                    "med_sclerophyll_resprout","med_sclerophyll_seed",
                                                                    "tall_sclerophyll_resprout","tall_sclerophyll_seed",
                                                                    "floodplain_sclerophyll","savanna_tree","subalpine_tree"))

### PLOT

ggplot(subset_input_df_euc,aes(x=stem_diameter_cm, y=height_m, col=pft_order)) +
  geom_point()+
  facet_wrap(species ~ .)+
  xlim(0, 200)+
  ylim(0, 100)

ggplot(subset(subset_input_df_euc,pft=="tall_sclerophyll_resprout"),aes(x=stem_diameter_cm, y=height_m, col=species)) +
  geom_point()+
  facet_wrap(species ~ .)+
  xlim(0, 200)+
  ylim(0, 100)



