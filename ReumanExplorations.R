# Exploration of Pearson-preserving instead of Spearman- or Kendall-preserving surrogates

rm(list=ls())

library(matrixcalc)

### Bring in the data - Shyamolina filled this in

#Importing the two dataset matrices of years by species
dHAY<- readRDS("./Results/hays_results/skewness_results/ts_mat_all_sp_146_hays.RDS")# for all species except in "REMOVE" category
dJRG<- readRDS("./Results/jrg_results/skewness_results/ts_mat_all_sp_39_jrg.RDS")# for all species except "BARE","ROCK","LASP"
  
#Importing data with "C", "I", "R" values for common, intermediate, and rare species
category_hays<-readRDS("./Results/hays_results/skewness_results/all_sp_146_hays_category.RDS")  
category_jrg<-readRDS("./Results/jrg_results/skewness_results/all_sp_39_jrg_category.RDS")

# check :
all(colnames(dHAY)==category_hays$sp)==T
all(colnames(dJRG)==category_jrg$sp)==T

#src stands for species rarity category
srcHAY<- category_hays$category #a vector with one entry for each species, corresponding to the columns of dHAY
srcJRG<- category_jrg$category

### Functions needed
source("copmap.R")
source("getinv.R")

### Now work with it to explore the map, starting with a couple common species
x<-dHAY[,4]
y<-dHAY[,7]

p<-seq(from=-1,to=1,by=0.05)
numreps<-500
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dashed")
lines(p,quants[3,],type='l',lty="dashed")
lines(range(p),rep(cor(x,y),2))
all.equal(res[,1],rep(res[1,1],numreps))
all.equal(res[,length(p)],rep(res[1,length(p)],numreps))
getinv(p,res,cor(x,y))

p<-seq(from=-1,to=1,by=0.2)
numreps<-50
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dashed")
lines(p,quants[3,],type='l',lty="dashed")
lines(range(p),rep(cor(x,y),2))
getinv(p,res,cor(x,y))
#so based on these two results I hypothesize it will be good enough
#to do 50 sims per value of p and to do values of p every .2 or so

### Now use a couple intermediate species
x<-dHAY[,1]
y<-dHAY[,2]

p<-seq(from=-1,to=1,by=0.05)
numreps<-500
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dotted")
lines(p,quants[3,],type='l',lty="dotted")
lines(range(p),rep(cor(x,y),2))
getinv(p,res,cor(x,y))
all.equal(res[,1],rep(res[1,1],numreps))
all.equal(res[,length(p)],rep(res[1,length(p)],numreps))

p<-seq(from=-1,to=1,by=0.2)
numreps<-50
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dotted")
lines(p,quants[3,],type='l',lty="dotted")
lines(range(p),rep(cor(x,y),2))
getinv(p,res,cor(x,y))
#this looks pretty awful

which(x!=0)
which(y!=0)
#So these two species never co-occur! Hard to expect success in what I am doing
#given that! Probably rare species will have to be removed.

#Seems to me we have problems with the skewness analysis for anything other than
#common species, just like we have problems with the Lt-Ut analysis. If there
#are not enough data to do copula fitting, then why would there be enough
#data to do resampling-based approaches to make it like the data came from a normal
#copula when we already know a copula as we consider them (i.e., for continuous
#distributions) cannot produce so many zeros.  

### Try another couple intermediate species, just to get and idea how general this
#problem is

x<-dHAY[,6]
y<-dHAY[,9]

p<-seq(from=-1,to=1,by=0.05)
numreps<-500
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dotted")
lines(p,quants[3,],type='l',lty="dotted")
lines(range(p),rep(cor(x,y),2))
getinv(p,res,cor(x,y))
all.equal(res[,1],rep(res[1,1],numreps))
all.equal(res[,length(p)],rep(res[1,length(p)],numreps))

p<-seq(from=-1,to=1,by=0.2)
numreps<-50
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dotted")
lines(p,quants[3,],type='l',lty="dotted")
lines(range(p),rep(cor(x,y),2))
getinv(p,res,cor(x,y))
#this looks less awful

which(x!=0)
which(y!=0)
#much more co-occurence, which probably explains the reduced awfulness

### Now go through many pairs of species and look for the inverse 
#point in each case

p<-seq(from=-1,to=1,by=0.2)
numreps<-50

numsp<-dim(dHAY)[2]
numsp<-20 #just do the first 20 species, for time reasons - we are 
#just trying to get a sense if this would work, and with which species
allparms<-matrix(NA,numsp,numsp)
for (c1 in 1:(numsp-1))
{
  print(paste0("c1: ",c1))
  for (c2 in (c1+1):numsp)
  {
    x<-dHAY[,c1]
    y<-dHAY[,c2]
    thisres<-copmap(x,y,p,numreps)
    allparms[c1,c2]<-getinv(p,thisres,cor(x,y))
  }
}
allparms[1:numsp,1:numsp]
#problems tend to be with the rare species

### Now go through pairs using only C species

dHAYC<-dHAY[,srcHAY=="C"]
p<-seq(from=-1,to=1,by=0.2)
numreps<-50

numsp<-dim(dHAYC)[2]
allparms<-matrix(NA,numsp,numsp)
for (c1 in 1:(numsp-1))
{
  print(paste0("c1: ",c1))
  for (c2 in (c1+1):numsp)
  {
    x<-dHAYC[,c1]
    y<-dHAYC[,c2]
    thisres<-copmap(x,y,p,numreps)
    allparms[c1,c2]<-getinv(p,thisres,cor(x,y))
  }
}
allparms[1:numsp,1:numsp]
sum(is.infinite(allparms))
#OK, so in this case the function seems always to be ivertible

#do you get a pos semi-def matrix?
allparms[is.na(allparms)]<-0
allparms<-allparms+t(allparms)
diag(allparms)<-1
is.symmetric.matrix(allparms)
is.positive.semi.definite(allparms)
eigen(allparms)$values
#nope, and some evals are meaningfully negative - this is a real problem

#so take the nearest PD matrix and see if you get an acceptable result
allparmsPD<-Matrix::nearPD(allparms,corr=TRUE)
allparmsPD<-matrix(as.numeric(allparmsPD$mat),nrow(allparmsPD$mat),ncol(allparmsPD$mat))
ncop<-normalCopula(param=P2p(allparmsPD),dim=dim(allparmsPD)[1],dispstr="un")
nums<-1000
sims<-rCopula((dim(dHAYC)[1])*nums,ncop)
sims<-aperm(array(sims,c(dim(dHAYC)[1],nums,dim(dHAYC)[2])),c(1,3,2))
remapd<-alignranks(dHAYC,sims)
allcors<-array(apply(FUN=cor,MARGIN=3,X=remapd),c(20,20,1000))
gtres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors>array(cor(dHAYC),c(20,20,1000))))
gtres
ltres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors<array(cor(dHAYC),c(20,20,1000))))
ltres
eqres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors==array(cor(dHAYC),c(20,20,1000))))
eqres
diag(gtres)<-NA
hist(gtres) #this shows the surrogates this way are terrible

#try it another way to guard against mistakes
library(mvtnorm)
nums<-1000
sims<-rmvnorm((dim(dHAYC)[1])*nums,sigma=allparmsPD)
sims<-aperm(array(sims,c(dim(dHAYC)[1],nums,dim(dHAYC)[2])),c(1,3,2))
remapd<-alignranks(dHAYC,sims)
allcors<-array(apply(FUN=cor,MARGIN=3,X=remapd),c(20,20,1000))
gtres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors>array(cor(dHAYC),c(20,20,1000))))
ltres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors<array(cor(dHAYC),c(20,20,1000))))
eqres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors==array(cor(dHAYC),c(20,20,1000))))
diag(gtres)<-NA
hist(gtres) #same result. So we cannot do this.

### Now repeat all the above with JRG, perhaps we will get lucky and 
#same problems will not present

### explore the map, starting with a couple common species
x<-dJRG[,1]
y<-dJRG[,3]

p<-seq(from=-1,to=1,by=0.05)
numreps<-500
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dashed")
lines(p,quants[3,],type='l',lty="dashed")
lines(range(p),rep(cor(x,y),2))
all.equal(res[,1],rep(res[1,1],numreps))
all.equal(res[,length(p)],rep(res[1,length(p)],numreps))
getinv(p,res,cor(x,y))

p<-seq(from=-1,to=1,by=0.2)
numreps<-50
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dashed")
lines(p,quants[3,],type='l',lty="dashed")
lines(range(p),rep(cor(x,y),2))
getinv(p,res,cor(x,y))
#so based on these two results I again hypothesize it will be good enough
#to do 50 sims per value of p and to do values of p every .2 or so

### Now use a couple intermediate species
x<-dJRG[,2]
y<-dJRG[,6]

p<-seq(from=-1,to=1,by=0.05)
numreps<-500
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dotted")
lines(p,quants[3,],type='l',lty="dotted")
lines(range(p),rep(cor(x,y),2))
getinv(p,res,cor(x,y))
all.equal(res[,1],rep(res[1,1],numreps))
all.equal(res[,length(p)],rep(res[1,length(p)],numreps))

p<-seq(from=-1,to=1,by=0.2)
numreps<-50
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dotted")
lines(p,quants[3,],type='l',lty="dotted")
lines(range(p),rep(cor(x,y),2))
getinv(p,res,cor(x,y))
#it basically worked in this case

which(x!=0)
which(y!=0)
#It is odd that it seemed to work because these two species only overlap for one year!

### Try another couple intermediate species
x<-dJRG[,7]
y<-dJRG[,8]

p<-seq(from=-1,to=1,by=0.05)
numreps<-500
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dotted")
lines(p,quants[3,],type='l',lty="dotted")
lines(range(p),rep(cor(x,y),2))
getinv(p,res,cor(x,y))
all.equal(res[,1],rep(res[1,1],numreps))
all.equal(res[,length(p)],rep(res[1,length(p)],numreps))

p<-seq(from=-1,to=1,by=0.2)
numreps<-50
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.025,.5,.975))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dotted")
lines(p,quants[3,],type='l',lty="dotted")
lines(range(p),rep(cor(x,y),2))
getinv(p,res,cor(x,y))

which(x!=0)
which(y!=0)
#much more co-occurence

### Now go through many pairs of species and look for the inverse 
#point in each case

p<-seq(from=-1,to=1,by=0.2)
numreps<-50

numsp<-dim(dJRG)[2]
numsp<-20 #just do the first 20 species, for time reasons - we are 
#just trying to get a sense if this would work, and with which species
allparms<-matrix(NA,numsp,numsp)
for (c1 in 1:(numsp-1))
{
  print(paste0("c1: ",c1))
  for (c2 in (c1+1):numsp)
  {
    x<-dJRG[,c1]
    y<-dJRG[,c2]
    thisres<-copmap(x,y,p,numreps)
    allparms[c1,c2]<-getinv(p,thisres,cor(x,y))
  }
}
allparms[1:numsp,1:numsp]
#problems tend to be with the rare and intermediate species

### Now go through pairs using only C species

dJRGC<-dJRG[,srcJRG=="C"]
p<-seq(from=-1,to=1,by=0.2)
numreps<-50

numsp<-dim(dJRGC)[2]
allparms<-matrix(NA,numsp,numsp)
for (c1 in 1:(numsp-1))
{
  print(paste0("c1: ",c1))
  for (c2 in (c1+1):numsp)
  {
    x<-dJRGC[,c1]
    y<-dJRGC[,c2]
    thisres<-copmap(x,y,p,numreps)
    allparms[c1,c2]<-getinv(p,thisres,cor(x,y))
  }
}
allparms[1:numsp,1:numsp]
sum(is.infinite(allparms))
#OK, so in this case the function seems always to be ivertible

#do you get a pos semi-def matrix?
allparms[is.na(allparms)]<-0
allparms<-allparms+t(allparms)
diag(allparms)<-1
isSymmetric(allparms)
is.positive.semi.definite(allparms)
eigen(allparms)$values
#nope, and some evals are meaningfully negative - this is a real problem, again!

#so take the nearest PD matrix and see if you get an acceptable result
allparmsPD<-Matrix::nearPD(allparms,corr=TRUE)
allparmsPD<-matrix(as.numeric(allparmsPD$mat),nrow(allparmsPD$mat),ncol(allparmsPD$mat))
ncop<-normalCopula(param=P2p(allparmsPD),dim=dim(allparmsPD)[1],dispstr="un")
nums<-1000
sims<-rCopula((dim(dHAYC)[1])*nums,ncop)
sims<-aperm(array(sims,c(dim(dHAYC)[1],nums,dim(dHAYC)[2])),c(1,3,2))
remapd<-alignranks(dHAYC,sims)
allcors<-array(apply(FUN=cor,MARGIN=3,X=remapd),c(20,20,1000))
gtres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors>array(cor(dHAYC),c(20,20,1000))))
gtres
ltres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors<array(cor(dHAYC),c(20,20,1000))))
ltres
eqres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors==array(cor(dHAYC),c(20,20,1000))))
eqres
diag(gtres)<-NA
hist(gtres) #this shows the surrogates this way are again terrible

#try it another way to guard against mistakes
library(mvtnorm)
nums<-1000
sims<-rmvnorm((dim(dHAYC)[1])*nums,sigma=allparmsPD)
sims<-aperm(array(sims,c(dim(dHAYC)[1],nums,dim(dHAYC)[2])),c(1,3,2))
remapd<-alignranks(dHAYC,sims)
allcors<-array(apply(FUN=cor,MARGIN=3,X=remapd),c(20,20,1000))
gtres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors>array(cor(dHAYC),c(20,20,1000))))
ltres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors<array(cor(dHAYC),c(20,20,1000))))
eqres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors==array(cor(dHAYC),c(20,20,1000))))
diag(gtres)<-NA
hist(gtres) #same result. So we cannot do this.




