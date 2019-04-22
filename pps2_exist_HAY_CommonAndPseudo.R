#This script looks for CV_com^2-approximately-preserving surrogates 
#for a dataset
#
#The script is currently written to work on the HAY dataset, with
#non-common species lumped into one pseudo-species. To change the 
#dataset on which the script operates, change the first section,
#then the rest should more or less run, probably with some hand-
#holding, though. This is draft code, it is expected that algorithm 
#and workflows will be optimized subsequently.
#
#Notation now based on my notes on pp. 36ff.
#
#Dan Reuman

rm(list=ls())

#***Change this section for different data and other settings

#**Data

#Importing the dataset matrix of years by species
d<-readRDS("./Results/hays_results/skewness_results/ts_mat_all_sp_146_hays.RDS")# for all species except "BARE","ROCK","LASP"

#Importing data with "C", "I", "R" values for common, intermediate, and rare species
category_hays<-readRDS("./Results/hays_results/skewness_results/all_sp_146_hays_category.RDS")
all(colnames(d)==category_hays$sp)==T
src<- category_hays$category #src stands for species rarity category

#Keep only the C species, combine the others into a pseudo-species
d<-cbind(d[,src=="C"],apply(FUN=sum,X=d[,src!="C"],MARGIN=1))

#**Settings

pfine<-seq(from=-1,to=1,by=0.01) #p stands for a bivariate normal copula parameter
numevalsfine<-500 #number of times you will evaluate the stochastic map
prough<-seq(from=-1,to=1,by=0.1)
numevalsrough<-250

resloc<-"./Results/Results_pps2_exist_HAY_CommonAndPseudo/"
theseed1<-103
theseed2<-1001
theseed3<-203

numpd<-100

#***Code imports, packages, and other setup

library(matrixcalc)
library(mvtnorm)
library(quantreg)

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
if (!dir.exists(paste0(resloc,"slices")))
{
  dir.create(paste0(resloc,"slices"))
}

#***For each pair of species, find the preimage of their correlation 
#under the map defined by copmap.R with the cflag="cor" option

set.seed(theseed1)

numsp<-dim(d)[2]

#receptacle for results
pijlo<-matrix(NA,numsp,numsp)
pijhi<-matrix(NA,numsp,numsp)
pijpre<-matrix(NA,numsp,numsp)
  
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
    thisres<-copmap(x,y,prough,numevalsrough,"cor")
    
    #find the inverse range
    h<-getinvrg(prough,thisres,cor(x,y))
    pijlo[c1,c2]<-h[1]
    pijhi[c1,c2]<-h[2]
    
    #find the actual preimage
    pijpre[c1,c2]<-getinv(prough,thisres,cor(x,y),"median")
    
    #make and save a plot
    mn<-apply(FUN=mean,MARGIN=2,X=thisres)
    quants<-apply(FUN=quantile,MARGIN=2,X=thisres,probs=c(.25,.5,.75))
    pdf(paste0(resloc,"fullmap/mapplot_",c1,"_",c2,".pdf"))
    plot(prough,mn,type='l',ylim=range(mn,quants),ylab="cor",xlab="norm cop param")
    lines(prough,quants[2,],type='l',lty="dashed")
    lines(prough,quants[1,],type='l',lty="dotted")
    lines(prough,quants[3,],type='l',lty="dotted")
    lines(range(prough),rep(cor(x,y),2))
    lines(rep(pijlo[c1,c2],2),c(-1,1))
    lines(rep(pijpre[c1,c2],2),c(-1,1))
    lines(rep(pijhi[c1,c2],2),c(-1,1))
    dev.off()
  }
}

#make sure the preimage is within the preimage range in all cases,
#and you got a sensible range and a single-point preimage in all cases
sum(pijlo<pijpre,na.rm=TRUE)
sum(pijpre<pijhi,na.rm=TRUE)
(numsp^2-numsp)/2
sum(is.finite(pijlo))
sum(is.finite(pijhi))
sum(is.finite(pijpre))
badindsARR<-rbind(which(pijlo>=pijpre,arr.ind = TRUE),
               which(pijhi<=pijpre,arr.ind = TRUE))
pvlo<-P2p(t(pijlo))
pvhi<-P2p(t(pijhi))
pvpre<-P2p(t(pijpre))
badinds<-which(pvlo>=pvpre | pvhi<=pvpre)
bds<-cbind(pvlo[badinds],pvpre[badinds],pvhi[badinds])
bds

#***Now look for positive definite

#First check to see if the preimage itself is pos def
m1<-pijpre
m1[is.na(m1)]<-0
m1<-m1+t(m1)
diag(m1)<-1
m2<-p2P(P2p(t(pijpre)))
sum(m1==m2)
prod(dim(m1))
is.positive.semi.definite(m1)
eigen(m1)$values
#it is not pd anyway

#now see if the midpoints of the pijlo-pijhi ranges give a pd mat
m1<-(pijlo+pijhi)/2
m1[is.na(m1)]<-0
m1<-m1+t(m1)
diag(m1)<-1
is.positive.semi.definite(m1)
eigen(m1)$values
#still not pd

#now use nearPD and see if you get a result between the bounds given by pijlo and pijhi
pijNPD<-matrix(as.numeric(Matrix::nearPD(m1,corr=TRUE)$mat),numsp,numsp)
pvNPD<-P2p(pijNPD)
inds<-which(pvNPD>pvhi | pvNPD<pvlo)
length(inds)
cbind(pvlo[inds],pvNPD[inds],pvhi[inds])
#so the nearPD matrix is well outside the bounds given by pijlo and pijhi for
#some of the coords, so I won't bother with it

#***Evaluate the stochastic map thoroughly between the bounds obtained 
#above, and perform fitting to nail down the expectation of the map as a
#(deterministic) function, and that function evaluated on pijlo 
#and pijhi (which will be used later for bounds for an optimization).

set.seed(theseed2)

#loop through each pair of species
coefsij<-array(NA,c(numsp,numsp,3))
cijlo<-matrix(NA,numsp,numsp)
cijhi<-matrix(NA,numsp,numsp)
for (c1 in 1:(numsp-1))
{
  print(paste0("c1: ",c1))
  for (c2 in (c1+1):numsp)
  {
    #get the time series for these two species
    x<-d[,c1]
    y<-d[,c2]
    
    #get the (stochastic) map
    pthis<-seq(from=pijlo[c1,c2],to=pijhi[c1,c2],length.out=30)
    thisres<-copmap(x,y,pthis,numevalsfine,"cor")
    
    #do a regression
    pthisreg<-as.vector(matrix(pthis,dim(thisres)[1],length(pthis),byrow = TRUE))
    thisresreg<-as.vector(thisres)
    psqthisreg<-pthisreg^2
    mreslm<-lm(thisresreg~pthisreg+psqthisreg)
    mresrq<-rq(thisresreg~pthisreg+psqthisreg)
      
    #make and save a plot
    mn<-apply(FUN=mean,MARGIN=2,X=thisres)
    quants<-apply(FUN=quantile,MARGIN=2,X=thisres,probs=c(.25,.5,.75))
    pdf(paste0(resloc,"partmap/mapplot_",c1,"_",c2,".pdf"))
    plot(pthis,mn,type='l',ylim=range(mn,quants),ylab="cor",xlab="norm cop param")
    lines(pthis,quants[2,],type='l',lty="dashed")
    lines(pthis,quants[1,],type='l',lty="dotted")
    lines(pthis,quants[3,],type='l',lty="dotted")
    psqthis<-pthis^2
    ylm<-coef(mreslm)[1]+coef(mreslm)[2]*pthis+coef(mreslm)[3]*psqthis
    yrq<-coef(mresrq)[1]+coef(mresrq)[2]*pthis+coef(mresrq)[3]*psqthis
    lines(pthis,ylm,type="l")
    lines(pthis,yrq)
    dev.off()
    
    #save the coefficients and v bounds
    coefsij[c1,c2,]<-unname(coef(mresrq))
    cijlo[c1,c2]<-yrq[1]
    cijhi[c1,c2]<-yrq[length(yrq)]
  }
}

#***Now do the optimization

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
  if (y<ylimits[1] || y>ylimits[2])
  {
    #print(paste0("y: ",y,"; ylimits[1]: ",ylimits[1],"; ylimits[2]: ",ylimits[2]))
    stop("Error in invval: out of range y argument")
  }
  
  c<-(c-y)
  if (isTRUE(all.equal(a,0)))
  {
    res<-(-c/b)
  }else
  {
    res<-(-b+sqrt(b^2-4*a*c))/(2*a)    
  }
  
  #more checks
  if (res<xlimits[1] || res>xlimits[2])
  {
    stop("Error in invval: res not in range")
  }
  
  return(res)
}

vijx<-unname(cov(d))
vtotx<-unname(sum(vijx))
viix<-unname(diag(vijx))
cijx<-unname(cor(d))

#This function returns a minimum eigenvalue and needs to be maximized
#
#Args
#cv      All the correlations, assembled into a vector using the approach of P2p
#
#Output: The minimal eigenvalue
fn<-function(cv)
{
  cij<-p2P(cv)

  #calculate the inverses to get pij
  pij<-matrix(NA,numsp,numsp)
  for (c1 in 1:(numsp-1))
  {
    for (c2 in (c1+1):numsp)
    {
      #print(paste0("c1: ",c1,"; c2: ",c2))
      pij[c1,c2]<-invval(threecoefs=coefsij[c1,c2,],y=cij[c1,c2],
                       ylimits=c(cijlo[c1,c2],cijhi[c1,c2]),xlimits=c(pijlo[c1,c2],pijhi[c1,c2]))
    }
  }
  pij[is.na(pij)]<-0
  pij<-pij+t(pij)
  diag(pij)<-1
  
  #calculate the minimum eigenvalue
  res<-min(eigen(pij)$values)
  if (res>0)
  {
    allresults[rowcount,]<<-c(P2p(p),P2p(v),res)
    rowcount<<-rowcount+1 #note the use of global variables, watch out!
    print("Found one!")
  }
  
  #print(res)
  #return(eigen(pij)$values)
  return(res)
}

#test the function on the values that correspond to the real data
cvx<-P2p(cijx)
fn(cv=cvx)
fn(cvx+c(-0.01,rep(0,length(cvx)-1))) #just to make sure it changes
fn(cvx+c(0.01,rep(0,length(cvx)-1))) 

#set up the constraints for optimizing
h<-2*sqrt(matrix(viix,numsp,numsp)*matrix(viix,numsp,numsp,byrow=TRUE))
h<-P2p(h)
ui<-rbind(diag(length(cvtest)), #for the constraints cij>=cijlo
          -diag(length(cvtest)), #for the constraints cij<=cijhi
          h, #for the complex constraint (see p. 41 of notes) 
          -h #for the complex constraint (see p. 41 of notes) 
)
ci<-c(cvlo, #for the constraints cij>=cijlo
      -cvhi, #for the constraints cij<=cijhi
      .25*vtotx-sum(viix), #for the complex constraint (see p. 41 of notes) 
      sum(viix)-1.75*vtotx #for the complex constraint (see p. 41 of notes) 
) 

#Try doing slices through the original initial condition to get a sense what this 
#objective function looks like. Or maybe do them through the best result I have so 
#far. 
lenplots<-100
y<-NA*numeric(lenplots+1)
maxy<-(-Inf)
miny<-Inf
for (ind in 1:length(cvx))
{
  x<-seq(from=cvlo[ind],to=cvhi[ind],length.out=lenplots)
  x<-c(x,cvx[ind])
  x<-sort(x)
  for (counter in 1:length(x))
  {
    cv<-cvx
    cv[ind]<-x[counter]
    if (all(ui%*%cv-ci>=1e-15))
    {
      y[counter]<-fn(cv)
    }else
    {
      y[counter]<-NA
    }
  }
  miny<-min(y,miny,na.rm=TRUE)
  maxy<-max(y,maxy,na.rm=TRUE)
  pdf(file=paste0(resloc,"slices/SliceAlongVar_",ind,".pdf"))
  myxlims<-range(x)
  myylims<-range(c(y,0),na.rm=TRUE)
  plot(x,y,type='l',xlim=myxlims,ylim=myylims)
  lines(rep(cvx[ind],2),range(c(y,0),na.rm=TRUE),type='l',lty="dashed")
  text(myxlims[2],myylims[2],sum(is.finite(y)),adj=c(1,1),cex=2)
  dev.off()
}

#do a search starting with correlations of the actual data
sum(ui%*%cvx-ci>=0)
length(ci) #check the constraint is satisfied by the initial condition

allresults<-matrix(NA,numpd,2*length(pvlo)+1)
rowcount<-1 #captures all evaluations that result in pd matrix

set.seed(theseed3)
optresNM<-constrOptim(cvx,fn,method="Nelder-Mead",ui=ui,ci=ci,
      control=list(fnscale=-1,maxit=500,trace=1)) 
optresNM$value

optresSANN<-constrOptim(cvx,fn,method="SANN",ui=ui,ci=ci,
      control=list(fnscale=-1,maxit=10000,trace=1,temp=100,tmax=10)) 
optresSANN$value

#now do some searches starting from other places
alloptres<-list()
alloptreslen<-1
for (counter in 1:1000000000)
{
  cvx_alt<-(cvhi-cvlo)*runif(length(cvlo))+cvlo
  if (all(ui%*%cvx_alt-ci>=0) && (fn(cvx_alt)>=fn(cvx)))
  {
    alloptres[[alloptreslen]]<-constrOptim(cvx_alt,fn,method="Nelder-Mead",ui=ui,ci=ci,
                           control=list(fnscale=-1,maxit=500,trace=1)) 
    alloptreslen<-alloptreslen+1
  }
}

for (counter in 1:(alloptreslen-1))
{
  print(alloptres[[counter]]$value)
}

