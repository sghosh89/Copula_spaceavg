#This script looks for CV_com^2-preserving surrogates for a dataset
#
#The script is currently written to work on the JRG dataset, with
#non-common species lumped into one pseudo-species. To change the 
#dataset on which the script operates, change the first section,
#then the rest should more or less run, perhaps with a bit of hand-
#holding. This is draft code, it is expected that algorithm and 
#workflows will be optemized subsequently.
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
theseed3<-203

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

#***For each pair of species, find the preimage of their covariance 
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
    
    #find the inverse range
    allprg[c1,c2,c(1,3)]<-getinvrg(prough,thisres,cov(x,y))
    
    #find the actual preimage
    allprg[c1,c2,2]<-getinv(prough,thisres,cov(x,y))
    
    #make and save a plot
    mn<-apply(FUN=mean,MARGIN=2,X=thisres)
    quants<-apply(FUN=quantile,MARGIN=2,X=thisres,probs=c(.25,.5,.75))
    pdf(paste0(resloc,"fullmap/mapplot_",c1,"_",c2,".pdf"))
    plot(prough,mn,type='l',ylim=range(mn,quants))
    lines(prough,quants[2,],type='l',lty="dashed")
    lines(prough,quants[1,],type='l',lty="dotted")
    lines(prough,quants[3,],type='l',lty="dotted")
    lines(range(prough),rep(cov(x,y),2))
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
#If m1 is pos semi def, we can skip some or all of the next steps

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

set.seed(theseed2)

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
    pthis<-seq(from=allprg[c1,c2,1],to=allprg[c1,c2,3],length.out=30)
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

set.seed(theseed3)

allresults<-matrix(NA,numpd,2*length(pvlo)+1)
rowcount<-1

covd<-cov(d)
vtotx<-sum(covd)
sumvarsx<-sum(diag(covd))

#For computing the inverse of the fitted expectation of copmap
#
#Args
#threecoefs       The coefficients of the map
#y                Take the inverse of this point under the map
#ylimits, xlimits   Should be the most extreme values
#
#Output: a single number
invval<-function(threecoefs,y,ylimits,xlimits)
{
  c<-threecoefs[1]
  b<-threecoefs[2]
  a<-threecoefs[3]

  #some checks
  h<-a*xlimits^2+b*xlimits+c
  if (!isTRUE(all.equal(h,ylimits)))
  {
    stop("Error in invval: inconsistent arguments")
  }
  #maybe later add a check that the deriv of the function is positive everywhere inside xlimits?
  #if (y<ylimits[1] || y>ylimits[2])
  #{
  #  print(paste0("y: ",y,"; ylimits[1]: ",ylimits[1],"; ylimits[2]: ",ylimits[2]))
  #  stop("Error in invval: out of range y argument")
  #}
  
  c<-(c-y)
  if (isTRUE(all.equal(a,0)))
  {
    res<-(-c/b)
  }else
  {
    res<-(-b+sqrt(b^2-4*a*c))/(2*a)    
  }
  
  #more checks
  #if (res<xlimits[1] || res>xlimits[2])
  #{
  #  stop("Error in invval: res not in range")
  #}
  
  return(res)
}

#This function returns a minimum eigenvalue and needs to be maximized
fn<-function(vpart)
{
  #get the remaining coordinate of v - the vs are covariances
  v12<-(vtotx-2*sum(vpart)-sumvarsx)/2
  
  #if that coordinate is out of the acceptable range, return a low value - can't return
  #-Inf because the optimizer I used needs finite values
  if (v12<vbds[1,2,1] || v12>vbds[1,2,2])
  {
    return(-9)
  }
  v<-p2P(c(v12,vpart))
  
  #calculate the inverses to get p
  p<-matrix(NA,numsp,numsp)
  for (c1 in 1:(numsp-1))
  {
    for (c2 in (c1+1):numsp)
    {
      #print(paste0("c1: ",c1,"; c2: ",c2))
      p[c1,c2]<-invval(threecoefs=allcoefs[c1,c2,],y=v[c1,c2],
                       ylimits=vbds[c1,c2,],xlimits=allprg[c1,c2,c(1,3)])
    }
  }
  p[is.na(p)]<-0
  p<-p+t(p)
  diag(p)<-1
  
  #calculate the minimum eigenvalue
  res<-min(eigen(p)$values)
  if (res>0)
  {
    allresults[rowcount,]<<-c(P2p(p),P2p(v),res)
    rowcount<<-rowcount+1 #note the use of global variables, watch out!
    print("Found one!")
  }
  
  #print(res)
  return(res)
}

vparttest<-P2p(covd)
vparttest<-vparttest[2:length(vparttest)]
fn(vparttest)
fn(vparttest+c(-0.01,rep(0,length(vparttest)-1))) #just to make sure it changes
fn(vparttest+c(0.01,rep(0,length(vparttest)-1))) 

#do a search starting with covariances of the actual data
vpart<-P2p(covd)
vpart<-vpart[2:length(vpart)]
lower<-P2p(t(vbds[,,1]))
lower<-lower[2:length(lower)]
upper<-P2p(t(vbds[,,2]))
upper<-upper[2:length(upper)]
optres<-optim(vpart,fn,method="L-BFGS-B",lower=lower,upper=upper,
      control=list(fnscale=-1)) 

#just eval fn at a bunch of points
lower10<-(99/100)*vpart+(1/100)*lower
upper10<-(99/100)*vpart+(1/100)*upper
for (counter in 1:100)
{
  startv<-(upper10-lower10)*runif(length(lower10))+lower10
  print(fn(startv))
}
#we are getting almost all values that don't work!

hist(log10(upper-lower)) #this is terrible, no wonder the optimizer does not work.
#Makes me think I should go back to cor.

#do some searches from nearby starts
lower10<-(9/10)*vpart+(1/10)*lower
upper10<-(9/10)*vpart+(1/10)*upper
start1<-(upper10-lower10)*runif(length(lower10))+lower10
fn(start1)
optres<-optim(start1,fn,method="L-BFGS-B",lower=lower,upper=upper,
              control=list(fnscale=-1)) 
