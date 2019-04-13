#This script looks for Pearson-preserving surrogates for a dataset
#
#The script is currently written to work on the JRG dataset, with
#non-common species lumped into one pseudo-species. To change the 
#dataset on which the script operates, change the first section,
#then the rest should run, perhaps with a bit of hand-holding.
#
#Dan Reuman
#2019 04 12

rm(list=ls())

#***Change this section for different data and other settings

#**Data

#Importing the dataset matrix of years by species
d<-readRDS("./Results/jrg_results/skewness_results/ts_mat_all_sp_39_jrg.RDS")# for all species except "BARE","ROCK","LASP"

#Importing data with "C", "I", "R" values for common, intermediate, and rare species
category_jrg<-readRDS("./Results/jrg_results/skewness_results/all_sp_39_jrg_category.RDS")
all(colnames(d)==category_jrg$sp)==T
src<- category_jrg$category #src stands for species rarity category

#Keep only the C species, combine the others into a pseudo-species
d<-cbind(d[,src=="C"],apply(FUN=sum,X=d[,src!="C"],MARGIN=1))

#**Settings

p<-seq(from=-1,to=1,by=0.1)
numreps<-250
resloc<-"./Results/Results_pps_exist_JRG_CommonAndPseudo/"
theseed1<-103
theseed2<-1001
numpd<-100000
numsimsforcheck<-1000

#***Code imports, packages, and other setup

library(matrixcalc)
library(mvtnorm)

source("copmap.R")
source("getinv.R")
source("getinvrg.R")
source("alignranks.R")

if (!dir.exists(resloc))
{
  dir.create(resloc)
}

#***For each pair of species, find the preimage of the corelation 
#under the map defined by copmap.R

set.seed(theseed1)

numsp<-dim(d)[2]

#receptacle for results, the preimage in allparms[,,2], the range
#specified by the other values
allparms<-array(NA,c(numsp,numsp,3))

#loop through each pair of species
for (c1 in 1:(numsp-1))
{
  print(paste0("c1: ",c1))
  for (c2 in (c1+1):numsp)
  {
    #get the time series for these two species
    x<-d[,c1]
    y<-d[,c2]
    
    #get the (stochastic) map
    thisres<-copmap(x,y,p,numreps)
    
    #find the iverse range
    allparms[c1,c2,c(1,3)]<-getinvrg(p,thisres,cor(x,y))
    
    #find the actual preimage
    allparms[c1,c2,2]<-getinv(p,thisres,cor(x,y))
    
    #make and save a plot
    mn<-apply(FUN=mean,MARGIN=2,X=thisres)
    quants<-apply(FUN=quantile,MARGIN=2,X=thisres,probs=c(.25,.5,.75))
    pdf(paste0(resloc,"mapplot_",c1,"_",c2,".pdf"))
    plot(p,mn,type='l',ylim=range(mn,quants))
    lines(p,quants[2,],type='l',lty="dashed")
    lines(p,quants[1,],type='l',lty="dotted")
    lines(p,quants[3,],type='l',lty="dotted")
    lines(range(p),rep(cor(x,y),2))
    lines(rep(allparms[c1,c2,1],2),c(-1,1))
    lines(rep(allparms[c1,c2,2],2),c(-1,1))
    lines(rep(allparms[c1,c2,3],2),c(-1,1))
    dev.off()
  }
}

#make sure the preimage is within the preimage range in all cases,
#and you got a sensible range and a single-point preimage in all cases
sum(allparms[,,1]<allparms[,,2],na.rm=TRUE)
sum(allparms[,,2]<allparms[,,3],na.rm=TRUE)
(numsp^2-numsp)/2

#***Now look for positive definite

#First check to see if we got lucky and the preimage itself is pos def
m1<-allparms[,,2]
m1[is.na(m1)]<-0
m1<-m1+t(m1)
diag(m1)<-1
m2<-p2P(P2p(m1))
sum(m1==m2)
prod(dim(m1))
is.positive.semi.definite(m1)
eigen(m1)$values
#If m1 is pos semi def, we can skip the next two steps

#now use nearPD and see if you get a result between hlo and hhi
hlo<-P2p(t(allparms[,,1]))
hpre<-P2p(t(allparms[,,2]))
hhi<-P2p(t(allparms[,,3]))
hmid<-(hlo+hhi)/2
bestmatNPD<-matrix(as.numeric(Matrix::nearPD(m2,corr=TRUE)$mat),numsp,numsp)
bestmatNPDvect<-P2p(bestmatNPD)
sum(bestmatNPDvect<=hhi & bestmatNPDvect>=hlo)
length(bestmatNPDvect)
inds<-which(bestmatNPDvect>hhi | bestmatNPDvect<hlo)
cbind(hlo[inds],bestmatNPDvect[inds],hhi[inds])
#so the nearPD matrix is well outside the bounds given by hlo and hhi for
#some of the coords, so I won't bother with it

#prepare to search the space of parameters desribed by allparms for 
#pos def
set.seed(theseed2)
allposdef<-matrix(NA,numpd,length(hlo)+1)
rowcount<-1

fn<-function(h)
{
  res<-min(eigen(p2P(h))$values)
  if (res>0)
  {
    allposdef[rowcount,]<<-c(h,res)
    rowcount<<-rowcount+1 #note the use of global variables, watch out!
    print("Found one!")
  }
  return(res)
}

#The below optimizations will crash out once they have found
#numpd positive definite results. This is intentional.

#do a search starting with the preimage itself
optim(hpre,fn,method="L-BFGS-B",lower=hlo,upper=hhi,
      control=list(fnscale=-1)) 

#do a search starting from the midpoint of the preimage range
optim((hlo+hhi)/2,fn,method="L-BFGS-B",lower=hlo,upper=hhi,
      control=list(fnscale=-1)) 

#now do a bunch of searches from random start points
while (all(is.na(allposdef[numpd,])))
{
  startvec<-(hhi-hlo)*runif(length(hlo))+hlo
  optim(startvec,fn,method="L-BFGS-B",lower=hlo,upper=hhi,
        control=list(fnscale=-1)) 
}

save(allposdef,file=paste0(resloc,"PosDefMats.RData"))

#Now find the one that is closest to hpre
allposdef<-allposdef[,1:(dim(allposdef)[2]-1)]
hpremat<-matrix(hpre,dim(allposdef)[1],length(hpre),byrow = TRUE)
sizemat<-matrix((hhi-hlo)/2,dim(allposdef)[1],length(hlo),byrow = TRUE)
dist<-apply(X=((allposdef-hpremat)/sizemat)^2,FUN=sum,MARGIN=1)
besth<-allposdef[which(dist==min(dist)),]
sum(besth>=hlo)
sum(besth<=hhi)
length(hhi)
bestmat<-p2P(besth) #this is supposed to be my parameter matrix for a normal copula!
isSymmetric(bestmat)
is.positive.semi.definite(bestmat)

#***Now make a plot showing where in the ranges given by hlo to hhi 
#the selected values sit

inds<-order(besth)
plot(1:length(besth),besth[inds],type='l',ylim=range(besth,hlo,hhi))
lines(1:length(besth),hlo[inds],type='l',lty="dashed")
lines(1:length(besth),hhi[inds],type='l',lty="dashed")

#***Now assess whether you get reasonable pearson-preserving surrogates
#using bestmat

sims<-rmvnorm((dim(d)[1])*numsimsforcheck,sigma=bestmat)
sims<-aperm(array(sims,c(dim(d)[1],numsimsforcheck,dim(d)[2])),c(1,3,2))
dim(sims)
sd<-d
for (counter in 1:numsp)
{
  sd[,counter]<-sort(d[,counter])
}#need to sort the columns of d to use aligncall
remapd<-alignranks(sd,sims)
holder<-apply(FUN=cor,MARGIN=3,X=remapd)
allcors<-array(holder,c(numsp,numsp,numsimsforcheck))
gtres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors>array(cor(d),c(numsp,numsp,numsimsforcheck))))
ltres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors<array(cor(d),c(numsp,numsp,numsimsforcheck))))
eqres<-apply(FUN=sum,MARGIN=c(1,2),X=(allcors==array(cor(d),c(numsp,numsp,numsimsforcheck))))

sum(diag(eqres)==numsimsforcheck)
length(diag(eqres))
diag(eqres)<-NA
sum(eqres==0,na.rm=TRUE)
(numsp^2-numsp)
inds<-which(eqres!=0,arr.ind = TRUE)
#so eq res is all 0s except the diagonal which equals numsimforcheck

sum(diag(gtres)==0)
sum(diag(ltres)==0)
diag(gtres)<-NA
diag(ltres)<-NA
hist(gtres)#very few of these should be outside the range 250-750, 
           #and if they are not outside that range you have pretty 
           #good surrogates
hist(ltres)#very few of these should be outside the range 250-750, 
           #and if they are not outside that range you have pretty 
           #good surrogates

save(bestmat,file=paste0(resloc,"FinalParameterMatrix.RData"))
#This is the final result. If you generate data from a normal copula
#with this parameter matrix and then use alignranks, that will give 
#pretty good Person-preserving surrogates.

#***Now produce 100000 surrogates and save for later use

sims<-rmvnorm((dim(d)[1])*100000,sigma=bestmat)
sims<-aperm(array(sims,c(dim(d)[1],100000,dim(d)[2])),c(1,3,2))
dim(sims)
sd<-d
for (counter in 1:numsp)
{
  sd[,counter]<-sort(d[,counter])
}#need to sort the columns of d to use aligncall
surrogs<-alignranks(sd,sims)

save(surrogs,file=paste0(resloc,"SomeSurrogates.RData"))

save.image(file=paste0(resloc,"WholeWorkspace.RData"))
