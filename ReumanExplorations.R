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

# Get the inverse range and see if you can find a PD matrix in there
source("getinvrg.R")

# explore the map for common species first, to make sure getinvrg is giving
#sensible results
x<-dHAY[,4]
y<-dHAY[,7]

p<-seq(from=-1,to=1,by=0.05)
numreps<-500
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.25,.5,.75))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dashed")
lines(p,quants[3,],type='l',lty="dashed")
lines(range(p),rep(cor(x,y),2))
all.equal(res[,1],rep(res[1,1],numreps))
all.equal(res[,length(p)],rep(res[1,length(p)],numreps))
invres<-getinvrg(p,res,cor(x,y))
lines(rep(invres[1],2),c(-1,1))
lines(rep(invres[2],2),c(-1,1))
invres

p<-seq(from=-1,to=1,by=0.1)
numreps<-50
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.25,.5,.75))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dashed")
lines(p,quants[3,],type='l',lty="dashed")
lines(range(p),rep(cor(x,y),2))
all.equal(res[,1],rep(res[1,1],numreps))
all.equal(res[,length(p)],rep(res[1,length(p)],numreps))
invres<-getinvrg(p,res,cor(x,y))
lines(rep(invres[1],2),c(-1,1))
lines(rep(invres[2],2),c(-1,1))
invres
#so based on these two results I hypothesize it will be good enough
#to do 50 sims per value of p and to do values of p every .1 or so

#now do all the common species, and for each pair find the preimage range
dHAYC<-dHAY[,srcHAY=="C"]
p<-seq(from=-1,to=1,by=0.1)
numreps<-50

numsp<-dim(dHAYC)[2]
allparms<-array(NA,c(numsp,numsp,2))
for (c1 in 1:(numsp-1))
{
  print(paste0("c1: ",c1))
  for (c2 in (c1+1):numsp)
  {
    x<-dHAYC[,c1]
    y<-dHAYC[,c2]
    thisres<-copmap(x,y,p,numreps)
    allparms[c1,c2,]<-getinvrg(p,thisres,cor(x,y))
  }
}
allparms[1:numsp,1:numsp,1]
allparms[1:numsp,1:numsp,2]

# now look for a pd matrix with entries subject to these bounds
numtries<-100000
pdres<-NA*numeric(numtries)
h1<-P2p(t(allparms[,,1]))
h2<-P2p(t(allparms[,,2]))
for (counter in 1:numtries)
{
  m<-p2P((h2-h1)*runif(length(h1))+h1)
  pdres[counter]<-is.positive.semi.definite(m)  
}
sum(pdres)
#did not find any pd mats! so we are probably screwed

#but lets do a better search than random guesses
fn<-function(h)
{
  res<-min(eigen(p2P(h))$values)
  counter<-1
  if (res>0)
  {
    save(res,file="pdmat.RData")
    print(paste0("Found one! ",counter))
    counter<-counter+1
  }
  return(res)
}
optim((h2+h1)/2,fn,method="L-BFGS-B",lower=h1,upper=h2,
      control=list(fnscale=-1)) #got -0.008126182
optim(h1+(h2-h1)/3,fn,method="L-BFGS-B",lower=h1,upper=h2,
      control=list(fnscale=-1)) #got -0.005551332
optim(h1+2*(h2-h1)/3,fn,method="L-BFGS-B",lower=h1,upper=h2,
      control=list(fnscale=-1)) #got -0.02355906
set.seed(101)
startvec<-(h2-h1)*runif(length(h1))+h1
optim(startvec,fn,method="L-BFGS-B",lower=h1,upper=h2,
      control=list(fnscale=-1)) #got -0.05817708
startvec<-(h2-h1)*runif(length(h1))+h1
optim(startvec,fn,method="L-BFGS-B",lower=h1,upper=h2,
      control=list(fnscale=-1)) #got -0.06747409
startvec<-(h2-h1)*runif(length(h1))+h1
optim(startvec,fn,method="L-BFGS-B",lower=h1,upper=h2,
      control=list(fnscale=-1)) #got -0.03248601
optim(h1+(h2-h1)/4,fn,method="L-BFGS-B",lower=h1,upper=h2,
      control=list(fnscale=-1)) #got -0.005046021
optim(h1+(h2-h1)/5,fn,method="L-BFGS-B",lower=h1,upper=h2,
      control=list(fnscale=-1)) #got -0.0001480897
optim(h1+(h2-h1)/6,fn,method="L-BFGS-B",lower=h1,upper=h2,
      control=list(fnscale=-1)) #got -0.003144349

#OK, need to come back and figure out a way to search overnight or soemthing.

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

# Get the inverse range and see if you can find a PD matrix in there
source("getinvrg.R")

# explore the map for common species first, to make sure getinvrg is giving
#sensible results
x<-dJRG[,1]
y<-dJRG[,3]

p<-seq(from=-1,to=1,by=0.05)
numreps<-500
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.25,.5,.75))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dashed")
lines(p,quants[3,],type='l',lty="dashed")
lines(range(p),rep(cor(x,y),2))
all.equal(res[,1],rep(res[1,1],numreps))
all.equal(res[,length(p)],rep(res[1,length(p)],numreps))
invres<-getinvrg(p,res,cor(x,y))
lines(rep(invres[1],2),c(-1,1))
lines(rep(invres[2],2),c(-1,1))
invres

p<-seq(from=-1,to=1,by=0.1)
numreps<-50
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.25,.5,.75))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dashed")
lines(p,quants[3,],type='l',lty="dashed")
lines(range(p),rep(cor(x,y),2))
all.equal(res[,1],rep(res[1,1],numreps))
all.equal(res[,length(p)],rep(res[1,length(p)],numreps))
invres<-getinvrg(p,res,cor(x,y))
lines(rep(invres[1],2),c(-1,1))
lines(rep(invres[2],2),c(-1,1))
invres
#so based on these two results I hypothesize it will be good enough
#to do 50 sims per value of p and to do values of p every .1 or so

#now do all the common species, and for each pair find the preimage range
dJRGC<-dJRG[,srcJRG=="C"]
p<-seq(from=-1,to=1,by=0.1)
numreps<-50

numsp<-dim(dJRGC)[2]
allparms<-array(NA,c(numsp,numsp,2))
allparmscent<-matrix(NA,numsp,numsp)
for (c1 in 1:(numsp-1))
{
  print(paste0("c1: ",c1))
  for (c2 in (c1+1):numsp)
  {
    x<-dJRGC[,c1]
    y<-dJRGC[,c2]
    thisres<-copmap(x,y,p,numreps)
    allparms[c1,c2,]<-getinvrg(p,thisres,cor(x,y))
    allparmscent[c1,c2]<-getinv(p,thisres,cor(x,y))
  }
}
#allparms[1:numsp,1:numsp,1]
#allparms[1:numsp,1:numsp,2]
sum(allparms[,,1]<allparmscent,na.rm=TRUE)
sum(!is.na(allparmscent))
sum(allparms[,,2]>allparmscent,na.rm=TRUE)

which(allparms[,,2]<allparmscent,arr.ind=TRUE)
#so check this one out
x<-dJRGC[,7]
y<-dJRGC[,13]
res<-copmap(x,y,p,numreps)
mn<-apply(FUN=mean,MARGIN=2,X=res)
quants<-apply(FUN=quantile,MARGIN=2,X=res,probs=c(.25,.5,.75))
plot(p,mn,type='l',ylim=range(mn,quants))
lines(p,quants[2,],type='l',lty="dashed")
lines(p,quants[1,],type='l',lty="dashed")
lines(p,quants[3,],type='l',lty="dashed")
lines(range(p),rep(cor(x,y),2))
getinv(p,res,cor(x,y))
allparms[7,13,]
allparmscent[9,12]<-



# now look for a pd matrix with entries subject to these bounds
numtries<-100000
pdres<-NA*numeric(numtries)
h1<-P2p(t(allparms[,,1]))
h2<-P2p(t(allparms[,,2]))
for (counter in 1:numtries)
{
  m<-p2P((h2-h1)*runif(length(h1))+h1)
  pdres[counter]<-is.positive.semi.definite(m)  
}
sum(pdres)
#did not find any pd mats! so we are screwed

#but lets do a better search than random guesses
fn<-function(h)
{
  res<-min(eigen(p2P(h))$values)
  counter<-1
  if (res>0)
  {
    save(res,file="pdmat.RData")
    print(paste0("Found one! ",counter))
    counter<-counter+1
  }
  return(res)
}
optim((h2+h1)/2,fn,method="L-BFGS-B",lower=h1,upper=h2,
      control=list(fnscale=-1)) #got 0.07191965
#We seem to have found one!
#We actually want one that is as close as possible to the middle of the ranges 
#we determined. How do we find that?
#Also, the sims for HAY suggest the likelihood surface is rough, so might want to
#use simulated annealing.

#How about we do a bunch of optims from various locations all over the bounding
#box, and we set it up to save all evaluations that have positive eigenvalues.
#Then we just find the closest ones.
fn<-function(h)
{
  res<-min(eigen(p2P(h))$values)
  if (res>0)
  {
    thisline<-matrix(c(h,res),1,1+length(h))
    write.table(thisline,file="AllPosEvals.csv",row.names=FALSE,col.names=FALSE,append=TRUE,sep=",")
    print("Found one!")
  }
  return(res)
}
optim((h2+h1)/2,fn,method="L-BFGS-B",lower=h1,upper=h2,
      control=list(fnscale=-1)) 
set.seed(101)
for (counter in 1:1000)
{
  startvec<-(h2-h1)*runif(length(h1))+h1
  optim(startvec,fn,method="L-BFGS-B",lower=h1,upper=h2,
        control=list(fnscale=-1)) 
}#ran this for 20 mins or so then cut it - the file got big fast!

#now search the results to find the one closest to the middle in all dimensions
pdmats<-read.csv(file="AllPosEvals.csv",header=FALSE)
pdmats<-pdmats[,1:(dim(pdmats)[2]-1)]
pdmats<-as.matrix(pdmats)
cent<-matrix((h1+h2)/2,dim(pdmats)[1],length(h1),byrow = TRUE)
size<-matrix((h2-h1)/2,dim(pdmats)[1],length(h1),byrow = TRUE)
dist<-apply(X=((pdmats-cent)/size)^2,FUN=sum,MARGIN=1)
besth<-unname(pdmats[which(dist==min(dist)),])
bestmat<-p2P(besth)
dim(bestmat)
head(bestmat) #this is supposed to be my parameter matrix for a normal copula!

#see if you get acceptable surrogates from this
ncop<-normalCopula(param=besth,dim=dim(bestmat)[1],dispstr="un")
nums<-1000
sims<-rCopula((dim(dJRGC)[1])*nums,ncop)
sims<-aperm(array(sims,c(dim(dJRGC)[1],nums,dim(dJRGC)[2])),c(1,3,2))
remapd<-alignranks(dJRGC,sims)
allcors<-array(apply(FUN=cor,MARGIN=3,X=remapd),c(23,23,1000))
gtres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors>array(cor(dJRGC),c(23,23,1000))))
gtres
ltres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors<array(cor(dJRGC),c(23,23,1000))))
ltres
eqres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors==array(cor(dJRGC),c(23,23,1000))))
eqres
diag(gtres)<-NA
hist(gtres) #this shows the surrogates this way are again terrible

#try it another way to guard against mistakes
library(mvtnorm)
nums<-1000
sims<-rmvnorm((dim(dJRGC)[1])*nums,sigma=bestmat)
sims<-aperm(array(sims,c(dim(dJRGC)[1],nums,dim(dJRGC)[2])),c(1,3,2))
dim(sims)
remapd<-alignranks(dJRGC,sims)
allcors<-array(apply(FUN=cor,MARGIN=3,X=remapd),c(23,23,1000))
gtres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors>array(cor(dJRGC),c(23,23,1000))))
gtres
ltres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors<array(cor(dJRGC),c(23,23,1000))))
ltres
eqres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors==array(cor(dJRGC),c(23,23,1000))))
diag(gtres)<-NA
hist(gtres) #same result. So we cannot do this.

#why is it so bad?