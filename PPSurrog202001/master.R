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

if (file.exists("MapsRes.RData"))
{
  load(file="MapsRes.RData")
}else
{
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
  save(mapx,mapy,file="MapsRes.RData")
}

#**now do a constrained optimization

source("PPSurrogObjFun.R")
source("pwlin.R")

#set up the constraints
#getting ready
vijmat<-cor(d)
vij<-vijmat[upper.tri(vijmat)]
pijstart<-vij[1:(length(vij)-1)]
sumvij<-sum(vij)
#first batch
ui1<-diag(length(pijstart))
ci1<-mapy[,,1]
ci1<-ci1[upper.tri(ci1)]
ci1<-ci1[-length(ci1)]
#second batch
ui2<--1*diag(length(pijstart))
ci2<-(-mapy[,,dim(mapy)[3]])
ci2<-ci2[upper.tri(ci2)]
ci2<-ci2[-length(ci2)]
#third and fourth batch
ui34<-matrix(rep(c(-1,1),length(pijstart)),2,length(pijstart))
ci34<-c(mapy[dim(mapy)[1]-1,dim(mapy)[2],1]-sumvij,
        sumvij-mapy[dim(mapy)[1]-1,dim(mapy)[2],dim(mapy)[3]])
#combine them all
ui<-rbind(ui1,ui2,ui34)
ci<-c(ci1,ci2,ci34)

#now do the optimization
res<-constrOptim(theta=pijstart,f=PPSurrogObj,grad=NULL,ui=ui,ci=ci,
                 sumvij=sumvij,mapx=mapx,mapy=mapy,
                 control=list(fnscale=-1,trace=1000,maxit=10000))

#now if it succeeded, then it saved a result, and can pull the result back 
#in from the disk and examine them
pijmat<-readRDS("FirstSuccess_pijmat.RDS")
nijmat<-readRDS("FirstSuccess_nijmat.RDS")
library(matrixcalc)
is.positive.definite(nijmat)
pij<-pijmat[upper.tri(pijmat)]
plot(vij,pij,type='p')
cor(vij,pij)

#next step: make surrogates using nijmat, and compute their CVcom^2 and 
#compare to the actual CVcom2 of the real data 