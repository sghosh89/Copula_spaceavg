rm(list=ls())
graphics.off()

#**Data

#Importing the dataset matrix of years by species
d<-readRDS("../Results/hays_results/skewness_results/ts_mat_all_sp_146_hays.RDS")# for all species except "BARE","ROCK","LASP"

#Importing data with "C", "I", "R" values for common, intermediate, and rare species
category_hays<-readRDS("../Results/hays_results/skewness_results/all_sp_146_hays_category.RDS")
all(colnames(d)==category_hays$sp)==T
src<- category_hays$category #src stands for species rarity category

#Keep only the C species, combine the others into a pseudo-species
d<-cbind(d[,src=="C"],pseudo=apply(FUN=sum,X=d[,src!="C"],MARGIN=1))

#**Now get the maps

source("getmap.R")

mapparms<-array(NA,c(dim(d)[2],dim(d)[2],3))
for (iind in 1:(dim(d)[2]-1))
{
  for (jind in (iind+1):(dim(d)[2]))
  {
    print(paste0("iind: ",iind,"; jind: ",jind))
    plotnm<-paste0("Map_",iind,"_",jind,".pdf")
    thisres<-getmap(d[,c(iind,jind)],numpts=25,numsims=100,plotnm=plotnm)
    mapparms[iind,jind,1]<-thisres$fit_parameters["expon"]
    mapparms[iind,jind,2]<-thisres$fit_parameters["prefactor"]
    mapparms[iind,jind,3]<-thisres$fit_parameters["intercept"]
  }
}