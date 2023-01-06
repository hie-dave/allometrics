# #Make up data
# age<-runif(10, 30,70)
# agesd<-runif(10, 0.1,5)
# bpf<-runif(10, 0,1)
# bpfsd<-runif(10, 0.01,.2)    
# pop.size<-runif(10,5,50)
# 
# #The plot
# plot(age,bpf, pch=16, cex=log(pop.size), col=rainbow(length(pop.size)), 
#      ylim=c(0,1),xlim=c(20,90))
# segments(age+agesd,bpf,age-agesd,bpf, lwd=2)
# segments(age,bpf+bpfsd,age,bpf-bpfsd, lwd=2)    
# legend("topright", legend=paste("Study",1:10), 
#        col=rainbow(length(pop.size)), pt.cex=1.5, pch=16)


sclerophyll_power_coeff = rbind(subset(species_coeff_df,(pft=="tall_sclerophyll_seed")),
             subset(species_coeff_df,(pft=="tall_sclerophyll_resprout")),
             subset(species_coeff_df,(pft=="med_sclerophyll_seed")),
             subset(species_coeff_df,(pft=="med_sclerophyll_resprout")),
             subset(species_coeff_df,(pft=="short_sclerophyll_seed")),
             subset(species_coeff_df,(pft=="short_sclerophyll_resprout")))

nameofpfts <- unique(sclerophyll_power_coeff$pft)
numberofpfts <- length(nameofpfts)
Base <-0
Power <-0
zmean<-0
xsd <-0
ysd <-0
for (i in 1:numberofpfts){
  Base[i]<- mean(sclerophyll_power_coeff$base[sclerophyll_power_coeff$pft==nameofpfts[i]], na.rm = TRUE)
  Power[i] <- mean(sclerophyll_power_coeff$power[sclerophyll_power_coeff$pft==nameofpfts[i]], na.rm = TRUE)  
  zmean[i] <- mean(sclerophyll_power_coeff$r_sqr[sclerophyll_power_coeff$pft==nameofpfts[i]], na.rm = TRUE)  
  xsd [i]<- sd(sclerophyll_power_coeff$base[sclerophyll_power_coeff$pft==nameofpfts[i]], na.rm = TRUE)
  ysd[i] <- sd(sclerophyll_power_coeff$power[sclerophyll_power_coeff$pft==nameofpfts[i]], na.rm = TRUE)
}


plot_scatter_with_error_bars <- function (xmean,ymean,zmean,xsd,ysd){
  MaxX <- max(xmean+xsd, na.rm = TRUE)
  MinX <- min(xmean-xsd, na.rm = TRUE)
  MaxY <- max(ymean+ysd, na.rm = TRUE)
  MinY <- min(ymean-ysd, na.rm = TRUE)
  # scatter plot
  # plot(xmean,ymean, pch=16, cex=log(zmean), col=rainbow(length(zmean)), 
  #      ylim=c(0,1),xlim=c(0,10))
  plot(xmean, ymean, pch=16, col=rainbow(length(zmean)), 
       xlim=c(MinX, MaxX), ylim=c(MinY, MaxY), xlab=("Base"), ylab=("Power"))
  # plot(xmean, ymean, pch=16, cex=(zmean)*2, col=rainbow(length(zmean)), xlim=c(-5,20), ylim=c(.2,.8))
  # x error bars
  segments(xmean+xsd,ymean,xmean-xsd,ymean, lwd=2)
  # y error bars
  segments(xmean,ymean+ysd,xmean,ymean-ysd, lwd=2)  
  
  legend("topright", legend=nameofpfts, inset=c(-.5, 0),
         col=rainbow(length(zmean)), pt.cex=1.5, pch=16)
  
}

plot_scatter_with_error_bars(Base,Power,zmean,xsd,ysd)

