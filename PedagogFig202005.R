rm(list=ls())

#***Functions that will be used below to generate data for a possibly improved pedagogical figure

#A function for generating data for a pedagogical figure, copula scale
#
#Args
#n        The number of random vectors you want
#sig      A covariance matrix. See the code and examples for how this is used
#
#Output - A matrix of dimensions n by 2*(dim(sig)[1]-1), draws from the random 
#variable corresponding to a 2*(dim(sig)[1]-1)-dimensional copula
#
getdat<-function(n,sig)
{
  coinflip<-runif(n)
  cop<-copula::normalCopula(param=P2p(sig),dim=dim(sig)[1],dispstr="un")
  intres<-rCopula(n,cop)
  d<-dim(sig)[1]-1

  #this will be the result to fill in
  res<-matrix(NA,n,2*d)
  
  #for about half the cases, the first d variables are in the interval (1/2,1)
  #and are the same, and the rest of the variables are in the interval (0,1/2) 
  res[coinflip<1/2,1:d]<-rep(intres[coinflip<1/2,1]/2+1/2,times=d)
  res[coinflip<1/2,(d+1):(2*d)]<-intres[coinflip<1/2,2:(d+1)]/2
  
  #for the other half of the cases, the last d variables are in the interval (1/2,1)
  #and are the same, and the first d variables are in the interval (0,1/2)
  res[coinflip>=1/2,1:d]<-intres[coinflip>=1/2,2:(d+1)]/2
  res[coinflip>=1/2,(d+1):(2*d)]<-rep(intres[coinflip>=1/2,1]/2+1/2,times=d)
  
  return(res)
}

#Skewness ratio. Assumes no NAs in the input matrix. Perhaps Shyamolina has a function
#like this somewhere that should be used instead of this one and this one should be 
#deleted to prevent redundancy.
#
#Args
#x      A matrix, time steps by species
#
source("SkewnessAnd3CentMom.R")
sr<-function(x)
{
  xtot<-apply(FUN=sum,X=x,MARGIN=1)
  
  scom<-myskns(xtot)
  sind<-(sum(apply(FUN=my3cm,MARGIN=2,X=x)))/(sum(apply(FUN=var,MARGIN=2,X=x)))^(3/2)
  
  return(c(scom=scom,sind=sind,sr=scom/sind))
}

#***Now generate the data, copula scale, for the first half of the figure

set.seed(104)

n<-10000
d<-5
#try 1
#syncconst<-0.5 
#compconst<-0
#try 2
#syncconst<-0.5 
#compconst<-0.2
#try 3
syncconst<-0.1 
compconst<-0.2
sig_p<-matrix(syncconst,d,d)
diag(sig_p)<-1
sig<-cbind(rep(-compconst,d),sig_p)
sig<-rbind(c(1,rep(-compconst,d)),sig)
sig
eigen(sig,symmetric=TRUE,only.values=TRUE)

res_p1<-getdat(n=10000,sig=sig)

hist(res_p1[,1],30)
hist(res_p1[,2],30)
hist(res_p1[,3],30)
hist(res_p1[,4],30) 
hist(res_p1[,5],30) 
hist(res_p1[,6],30) 
hist(res_p1[,7],30) 
hist(res_p1[,8],30) 
hist(res_p1[,9],30) 
hist(res_p1[,10],30) #histograms seem uniform, as they should be

plot(res_p1[1:500,1],res_p1[1:500,2],type="p",pch=20)
plot(res_p1[1:500,6],res_p1[1:500,7],type="p",pch=20)

plot(res_p1[1:500,1],res_p1[1:500,6],type="p",pch=20)
plot(res_p1[1:500,1],res_p1[1:500,7],type="p",pch=20)
plot(res_p1[1:500,2],res_p1[1:500,6],type="p",pch=20)
plot(res_p1[1:500,2],res_p1[1:500,7],type="p",pch=20)

#***Now generate the data, copula scale, for the second half of the figure
res_p2<-getdat(n=10000,sig=sig)
res_p2<-(-res_p2+1)

hist(res_p2[,1],30)
hist(res_p2[,2],30)
hist(res_p2[,3],30)
hist(res_p2[,4],30) 
hist(res_p2[,5],30) 
hist(res_p2[,6],30) 
hist(res_p2[,7],30) 
hist(res_p2[,8],30) 
hist(res_p2[,9],30) 
hist(res_p2[,10],30) #histograms seem uniform, as they should be

plot(res_p2[1:500,1],res_p2[1:500,2],type="p",pch=20)
plot(res_p2[1:500,6],res_p2[1:500,7],type="p",pch=20)

plot(res_p2[1:500,1],res_p2[1:500,6],type="p",pch=20)
plot(res_p2[1:500,1],res_p2[1:500,7],type="p",pch=20)
plot(res_p2[1:500,2],res_p2[1:500,6],type="p",pch=20)
plot(res_p2[1:500,2],res_p2[1:500,7],type="p",pch=20)

#***now combine both copulas with standard normal marginals

normres_p1<-qnorm(res_p1)
normres_p2<-qnorm(res_p2)

normres_p1<-normres_p1-min(normres_p1,normres_p2)+1
normres_p2<-normres_p2-min(normres_p2,normres_p2)+1

#***now compare the two resulting datasets in various respects

#very similar species marginal distributions in the two cases
empcdf<-data.frame(x=sort(normres_p1[,1]),y=(1:n)/n)
plot(empcdf[,"x"],empcdf[,"y"],type="l",xlim=c(1,10),ylim=c(0,1))
for (counter in 2:d)
{
  empcdf<-data.frame(x=sort(normres_p1[,counter]),y=(1:n)/n)
  lines(empcdf[,"x"],empcdf[,"y"],type="l")
}
for (counter in (d+1):(2*d))
{
  empcdf<-data.frame(x=sort(normres_p1[,counter]),y=(1:n)/n)
  lines(empcdf[,"x"],empcdf[,"y"],type="l",col="red")
}

#very similar mean of the total
tot_p1<-apply(FUN=sum,MARGIN=1,X=normres_p1)
tot_p2<-apply(FUN=sum,MARGIN=1,X=normres_p2)
mean(tot_p1)
mean(tot_p2)

#very similar variance of the total
var(tot_p1)
var(tot_p2)

#very similar variance ratio
tsvr::vr(t(normres_p1))
tsvr::vr(t(normres_p2))

#very different skewnesses of the total
hist(tot_p1,50)
myskns(tot_p1)
hist(tot_p2,50)
myskns(tot_p2)

#extremely different skewness ratios
sr(normres_p1)
sr(normres_p2)

#***Now make some plots of the time series. The real figure also show some similar plots 
#but will be much nicer than this.

#first set of time series
tim<-1100:1150
plot(tim,normres_p1[tim,1],type="l",col="red",ylim=range(normres_p1[tim,],normres_p2[tim,]))
for (counter in 2:d)
{
  lines(tim,normres_p1[tim,counter],col="red")
}
for (counter in (d+1):(2*d))
{
  lines(tim,normres_p1[tim,counter],col="black")
}

#second set of time series
plot(tim,normres_p2[tim,1],type="l",col="red",ylim=range(normres_p1[tim,],normres_p2[tim,]))
for (counter in 2:d)
{
  lines(tim,normres_p2[tim,counter],col="red")
}
for (counter in (d+1):(2*d))
{
  lines(tim,normres_p2[tim,counter],col="black")
}

#both totals
plot(tim,tot_p1[tim],type="l",ylim=range(tot_p1[tim],tot_p2[tim]))
lines(tim,tot_p2[tim],type="l",col="red")

#The advantage of this approach over the fig currently used in the manuscript is the vr is 
#realistic in this one