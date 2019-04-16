#This script looks for CV_com^2-preserving surrogates for a dataset
#
#The script is currently written to work on the JRG dataset, with
#non-common species lumped into one pseudo-species. To change the 
#dataset on which the script operates, change the first section,
#then the rest should more or less run, perhaps with a bit of hand-
#holding.
#
#Dan Reuman
#2019 04 15

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

pfine<-seq(from=-1,to=1,by=0.01) #p stands for a bivariate normal copula parameter
numevalsfine<-1000 #number of times you will evaluate the stochastic map
prough<-seq(from=-1,to=1,by=0.1)
numevalsrough<-250

resloc<-"./Results/Results_pps2_exist_JRG_CommonAndPseudo/"
theseed1<-103
theseed2<-1001

numpd<-10000 #find at least this many positive definite matrices of parameters before stopping

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
if (!dir.exists(paste0(resloc,"fullmap")))
{
  dir.create(paste0(resloc,"fullmap"))
}
if (!dir.exists(paste0(resloc,"partmap")))
{
  dir.create(paste0(resloc,"partmap"))
}

#***For each pair of species, find the preimage of their corelation 
#under the map defined by copmap.R

set.seed(theseed1)

numsp<-dim(d)[2]

#receptacle for results, the preimage in allprg[,,2], the range
#specified by the other values
allprg<-array(NA,c(numsp,numsp,3))

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
    thisres<-copmap(x,y,prough,numevalsrough)
    
    #find the iverse range
    allprg[c1,c2,c(1,3)]<-getinvrg(prough,thisres,cor(x,y))
    
    #find the actual preimage
    allprg[c1,c2,2]<-getinv(prough,thisres,cor(x,y))
    
    #make and save a plot
    mn<-apply(FUN=mean,MARGIN=2,X=thisres)
    quants<-apply(FUN=quantile,MARGIN=2,X=thisres,probs=c(.25,.5,.75))
    pdf(paste0(resloc,"fullmap/mapplot_",c1,"_",c2,".pdf"))
    plot(prough,mn,type='l',ylim=range(mn,quants))
    lines(prough,quants[2,],type='l',lty="dashed")
    lines(prough,quants[1,],type='l',lty="dotted")
    lines(prough,quants[3,],type='l',lty="dotted")
    lines(range(prough),rep(cor(x,y),2))
    lines(rep(allprg[c1,c2,1],2),c(-1,1))
    lines(rep(allprg[c1,c2,2],2),c(-1,1))
    lines(rep(allprg[c1,c2,3],2),c(-1,1))
    dev.off()
  }
}

#make sure the preimage is within the preimage range in all cases,
#and you got a sensible range and a single-point preimage in all cases
sum(allprg[,,1]<allprg[,,2],na.rm=TRUE)
sum(allprg[,,2]<allprg[,,3],na.rm=TRUE)
(numsp^2-numsp)/2

#***Now look for positive definite

#First check to see if we got lucky and the preimage itself is pos def
m1<-allprg[,,2]
m1[is.na(m1)]<-0
m1<-m1+t(m1)
diag(m1)<-1
m2<-p2P(P2p(m1))
sum(m1==m2)
prod(dim(m1))
is.positive.semi.definite(m1)
eigen(m1)$values
#If m1 is pos semi def, we can skip some of the next steps

#now use nearPD and see if you get a result between hlo and hhi
pvlo<-P2p(t(allprg[,,1])) #pv stands for p in vector format
pvpre<-P2p(t(allprg[,,2]))
pvhi<-P2p(t(allprg[,,3]))
pvmid<-(pvlo+pvhi)/2
pmatNPD<-matrix(as.numeric(Matrix::nearPD(m2,corr=TRUE)$mat),numsp,numsp)
pvNPD<-P2p(pmatNPD)
sum(pvNPD<=pvhi & pvNPD>=pvlo)
length(pvNPD)
inds<-which(pvNPD>pvhi | pvNPD<pvlo)
cbind(pvlo[inds],pvNPD[inds],pvhi[inds])
#so the nearPD matrix is well outside the bounds given by hlo and hhi for
#some of the coords, so I won't bother with it

#***Evaluate the stochastic map thoroughly between the bounds obtained 
#above, and perform fitting to nail down the expectation of the map as a
#(deterministic) function, and that function evaluated on allprg[,,1] 
#and allprg[,,3] (which will be used later for bounds for an optimization).

#loop through each pair of species
allcoefs<-array(NA,c(numsp,numsp,3))
vbds<-array(NA,c(numsp,numsp,2))
for (c1 in 1:(numsp-1))
{
  print(paste0("c1: ",c1))
  for (c2 in (c1+1):numsp)
  {
    #get the time series for these two species
    x<-d[,c1]
    y<-d[,c2]
    
    #get the (stochastic) map
    pthis<-seq(from=allprg[c1,c2,1],to=allprg[c1,c2,3],by=0.01)
    thisres<-copmap(x,y,pthis,numevalsfine)
    
    #do a regression
    pthisreg<-as.vector(matrix(pthis,dim(thisres)[1],length(pthis),byrow = TRUE))
    thisresreg<-as.vector(thisres)
    psqthisreg<-pthisreg^2
    mres<-lm(thisresreg~pthisreg+psqthisreg)
    
    #make and save a plot
    mn<-apply(FUN=mean,MARGIN=2,X=thisres)
    quants<-apply(FUN=quantile,MARGIN=2,X=thisres,probs=c(.25,.5,.75))
    pdf(paste0(resloc,"partmap/mapplot_",c1,"_",c2,".pdf"))
    plot(pthis,mn,type='l',ylim=range(mn,quants))
    lines(pthis,quants[2,],type='l',lty="dashed")
    lines(pthis,quants[1,],type='l',lty="dotted")
    lines(pthis,quants[3,],type='l',lty="dotted")
    psqthis<-pthis^2
    y<-coef(mres)[1]+coef(mres)[2]*pthis+coef(mres)[3]*psqthis
    lines(pthis,y,type="l")
    dev.off()
    
    #save the coefficients and v bounds
    allcoefs[c1,c2,]<-unname(coef(mres))
    vbds[c1,c2,]<-c(y[1],y[length(y)])
  }
}

#***Now do the optimization

set.seed(theseed2)

allposdef<-matrix(NA,numpd,length(pvlo)+1)
rowcount<-1


#This function returns a minimum eigenvalue and needs to be maximized
fn<-function(vpart)
{
  #get the remaining coordinate of v - the vs are correlations
    #convert to covariances
    #
  
  #if that coordinate is out of the acceptable range, return -Inf
  
  
  #calculate the inverses to get p
  
  
  #calculate the minimum eigenvalue
  res<-min(eigen(p2P(p))$values)
  if (res>0)
  {
    allposdef[rowcount,]<<-c(p,res)
    rowcount<<-rowcount+1 #note the use of global variables, watch out!
    print("Found one!")
  }
  
  return(res)
}







