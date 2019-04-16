#Explorations trying to get toward CV_com-preserving normal surrogates

rm(list=ls())

#**Data

#Importing the dataset matrix of years by species
d<-readRDS("./Results/hays_results/skewness_results/ts_mat_all_sp_146_hays.RDS")# for all species except in "REMOVE" category

#Importing data with "C", "I", "R" values for common, intermediate, and rare species
category_hays<-readRDS("./Results/hays_results/skewness_results/all_sp_146_hays_category.RDS")
all(colnames(d)==category_hays$sp)==T
src<- category_hays$category #a vector with one entry for each species, corresponding to the columns of dHAY
#src stands for species rarity category

#Keep only the C species, combine the others into a pseudo-species
d<-cbind(d[,src=="C"],apply(FUN=sum,X=d[,src!="C"],MARGIN=1))

#**other setup

#***Code imports, packages, and other setup

library(matrixcalc)
library(mvtnorm)

source("copmap.R")
source("getinv.R")
source("getinvrg.R")
source("alignranks.R")

#**study the copmap result, see if I can fit something to it

x<-d[,1]
y<-d[,2]

p<-seq(from=-1,to=1,by=0.02)
numreps<-250
thisres<-copmap(x,y,p,numreps)

mn<-apply(FUN=mean,MARGIN=2,X=thisres)
quants<-apply(FUN=quantile,MARGIN=2,X=thisres,probs=c(.025,.25,.5,.75,.975))
plot(p,mn,ylim=c(-1,1),type='l')
lines(p,quants[3,],type="l",lty="dashed")
lines(p,quants[1,],type="l",lty="dotted")
lines(p,quants[2,],type="l",lty="dotted")
lines(p,quants[4,],type="l",lty="dotted")
lines(p,quants[5,],type="l",lty="dotted")
