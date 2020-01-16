rm(list=ls())
source("PPSurrogObjFun.R")
source("pwlin.R")

#**get some data to work with

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

set.seed(101)
numpts<-18
mapy<-array(NA,c(dim(d)[2],dim(d)[2],numpts))
for (iind in 1:(dim(d)[2]-1))
{
  for (jind in (iind+1):(dim(d)[2]))
  {
    print(paste0("iind: ",iind,"; jind: ",jind))
    plotnm<-paste0("Map_",iind,"_",jind)
    thisres<-getmap(d[,c(iind,jind)],numpts=numpts,numsims=500,plotnm=plotnm)
    mapx<-thisres$fit_parameters$x
    y<-thisres$fit_parameters$y
    if (any(diff(y)<0)){stop("Error: non-monotonic piecewise linear interpolation")}
    mapy[iind,jind,]<-thisres$fit_parameters$y
  }
}

#**now try it out

vijmat<-cor(d)
vij<-vijmat[upper.tri(vijmat)]
pij<-vij[1:(length(vij)-1)]
sumvij<-sum(vij)
PPSurrogObj(pij=pij,sumvij=sumvij,mapx=mapx,mapy=mapy)
PPSurrogObj(c(pij[1]+.001,pij[2:length(pij)]),sumvij,mapx,mapy)
PPSurrogObj(c(pij[1]-.001,pij[2:length(pij)]),sumvij,mapx,mapy) #so it is possible to make improvements

#Note this is not an exacting test, it only shows the function runs and produces a numeric
#value! Maybe add a better test later. 

