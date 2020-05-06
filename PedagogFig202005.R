rm(list=ls())

#Generates data from a multivariate copula I constructed - see math notes dated 2020 05 04,
#entitled "Creating a copula with certain properties", though there are some notational
#differences between there and here, and the written-up math is a more specific case
#(which corresponds to the test below).
#
#Args
#n              Number of random vectors to generate
#sig1           A covariance matrix of dimensions equal to the desired dimension for the copula
#sig2           Another covariance matrix
#vartovar       A vector of length dim(sig1)[1] consisting solely of the integers between 1 and
#                 dim(sig2), inclusive. 
#
#Output - an n by dim(sig1)[1] matrix of copula values
#
rnewcop<-function(n,sig1,sig2,vartovar,betas)
{
  ncop1<-copula::normalCopula(param=P2p(sig1),dim=dim(sig1)[1],dispstr="un")
  ncop2<-copula::normalCopula(param=P2p(sig2),dim=dim(sig2)[1],dispstr="un")
  
  dcop1<-copula::rCopula(n,ncop1)
  dcop2<-copula::rCopula(n,ncop2)
  
  for (counter in 1:dim(sig2)[1])
  {
    inds_var<-which(vartovar==counter)
    indicator1<-0*numeric(n)
    indicator1[dcop2[,counter]>betas[counter]]<-1
    indicator2<-as.numeric(!indicator1)

    dcop1[,inds_var]<-indicator1*dcop2[,counter]+indicator2*betas[counter]*dcop1[,inds_var]
  }
  
  return(dcop1)
}

#Now test it
rho<-.8
n_math<-12 #this is the n from the math notes
eta<-(n_math/2-1)*rho/(n_math/2)
sig11<-matrix(rho,n_math/2,n_math/2)
diag(sig11)<-1
sig22<-sig11
sig12<-matrix(-eta,n_math/2,n_math/2)
sig21<-sig12
sig1<-rbind(cbind(sig11,sig12),cbind(sig21,sig22))
sig1
sum(sig1) #should equal n, making a vr of 1
eigen(sig1,symmetric=TRUE,only.values=TRUE)
sig2<-matrix(c(1,-.8,-.8,1),2,2)
vartovar<-c(rep(1,n_math/2),rep(2,n_math/2))
n_func<-10000 #this is the n argument for the function above
betas<-c(.8,.8)
res<-rnewcop(n_func,sig1,sig2,vartovar,betas)

plot(res[1:500,1],res[1:500,2],type="p",pch=20)
plot(res[1:500,1],res[1:500,3],type="p",pch=20)
plot(res[1:500,2],res[1:500,3],type="p",pch=20)

plot(res[1:500,n_math/2+1],res[1:500,n_math/2+2],type="p",pch=20)
plot(res[1:500,n_math/2+1],res[1:500,n_math/2+3],type="p",pch=20)
plot(res[1:500,n_math/2+2],res[1:500,n_math/2+3],type="p",pch=20)

plot(res[1:500,1],res[1:500,n_math/2+1],type="p",pch=20)
plot(res[1:500,1],res[1:500,n_math/2+2],type="p",pch=20)
plot(res[1:500,1],res[1:500,n_math/2+3],type="p",pch=20)

hist(res[,1],30)
hist(res[,2],30)
hist(res[,3],30)
hist(res[,n_math/2+1],30)
hist(res[,n_math/2+2],30)
hist(res[,n_math/2+3],30) #illustrates the fact that the marginals are uniform


cor(res)
res_norm<-qnorm(res)
cm<-cov(res_norm)
tsvr::vr(t(res_norm)-min(res_norm)+1)
tot<-apply(FUN=sum,MARGIN=1,X=res_norm)
hist(tot,50)
source("SkewnessAnd3CentMom.R")
myskns(tot)

#The above does not seem to be giving me what I was hoping for. Let's try a bone-headed way of doing it.

getdat1<-function(n)
{
  res<-matrix(NA,n,4)
  
  thirds<-runif(n)
  indstop<-which(thirds<1/3)
  indsmid<-which(thirds>=1/3 & thirds<=2/3)
  indsbot<-which(thirds>2/3)
  
  #A third of the time, variables 1 and 2 are the same and uniformly distributed between 2/3 and 1,
  #and variables 3 and 4 are independent and uniformly distributed beyween 0 and 1/3
  res[indstop,1:2]<-rep(runif(length(indstop),min=2/3,max=1),times=2)
  res[indstop,3:4]<-runif(length(indstop)*2,0,1/3)
  
  #A third of the time, variables 3 and 4 are the same and uniformly distributed between 2/3 and 1,
  #and variables 1 and 2 are independent and uniformly distributed beyween 0 and 1/3
  res[indsbot,3:4]<-rep(runif(length(indsbot),min=2/3,max=1),times=2)
  res[indsbot,1:2]<-runif(length(indsbot)*2,0,1/3)
  
  #A third of the time, all four variables are independent and uniformly distributed between 1/3
  #and 2/3
  res[indsmid,]<-runif(length(indsmid)*4,1/3,2/3)
  
  return(res)
}

res<-getdat(10000)
hist(res[,1],30)
hist(res[,2],30)
hist(res[,3],30)
hist(res[,4],30) #histograms seem uniform, as they should be
plot(res[1:500,1],res[1:500,2],type="p",pch=20)
plot(res[1:500,1],res[1:500,3],type="p",pch=20)
#This is probably working, but why use three parts?

#A simpler alternative using 2 parts
getdat2<-function(n)
{
  res<-matrix(NA,n,4)
  
  halves<-runif(n)
  indstop<-which(halves<1/2)
  indsbot<-which(halves>=1/2)
  
  #Half of the time, variables 1 and 2 are the same and uniformly distributed between 1/2 and 1,
  #and variables 3 and 4 are independent and uniformly distributed beyween 0 and 1/2
  res[indstop,1:2]<-rep(runif(length(indstop),min=1/2,max=1),times=2)
  res[indstop,3:4]<-runif(length(indstop)*2,0,1/2)
  
  #The other half of the time, variables 3 and 4 are the same and uniformly distributed between 1/2 and 1,
  #and variables 1 and 2 are independent and uniformly distributed beyween 0 and 1/2
  res[indsbot,3:4]<-rep(runif(length(indsbot),min=1/2,max=1),times=2)
  res[indsbot,1:2]<-runif(length(indsbot)*2,0,1/2)

  return(res)
}

res<-getdat2(10000)
hist(res[,1],30)
hist(res[,2],30)
hist(res[,3],30)
hist(res[,4],30) #histograms seem uniform, as they should be
plot(res[1:500,1],res[1:500,2],type="p",pch=20)
plot(res[1:500,3],res[1:500,4],type="p",pch=20)

plot(res[1:500,1],res[1:500,3],type="p",pch=20)
plot(res[1:500,1],res[1:500,4],type="p",pch=20)
plot(res[1:500,2],res[1:500,3],type="p",pch=20)
plot(res[1:500,2],res[1:500,4],type="p",pch=20)

cov(res)
cor(res)

normres<-qnorm(res)
normres<-normres-min(normres)+1
tsvr::vr(t(normres))



#A maybe more flexbile and controlled alternative, but perhaps not general yet
getdat3<-function(n,sig)
{
  coinflip<-runif(n)
  cop<-copula::normalCopula(param=P2p(sig),dim=dim(sig)[1],dispstr="un")
  intres<-rCopula(n,cop)
  
  res<-matrix(NA,n,4)
  #for about half the cases 
  res[coinflip<1/2,1:2]<-rep(intres[coinflip<1/2,1]/2+1/2,times=2)
  res[coinflip<1/2,3:4]<-intres[coinflip<1/2,2:3]/2
    
  #for the other half
  res[coinflip>=1/2,1:2]<-intres[coinflip>=1/2,2:3]/2
  res[coinflip>=1/2,3:4]<-rep(intres[coinflip>=1/2,1]/2+1/2,times=2)
  
  return(res)
}

n<-10000
#sig<-matrix(c(1,-.8,-.8,-.8,1,.35,-.8,.35,1),3,3)
#sig<-matrix(c(1,-.7,-.7,-.7,1,.25,-.7,.25,1),3,3)
#sig<-matrix(c(1,-0,-0,-0,1,0,-0,0,1),3,3)
sig<-matrix(c(1,-0,-0,-0,1,.5,-0,.5,1),3,3)
sig
eigen(sig,symmetric=TRUE,only.values=TRUE)
res<-getdat3(n=10000,sig=sig)

hist(res[,1],30)
hist(res[,2],30)
hist(res[,3],30)
hist(res[,4],30) #histograms seem uniform, as they should be

plot(res[1:500,1],res[1:500,2],type="p",pch=20)
plot(res[1:500,3],res[1:500,4],type="p",pch=20)

plot(res[1:500,1],res[1:500,3],type="p",pch=20)
plot(res[1:500,1],res[1:500,4],type="p",pch=20)
plot(res[1:500,2],res[1:500,3],type="p",pch=20)
plot(res[1:500,2],res[1:500,4],type="p",pch=20)

normres<-qnorm(res)
cov(normres)
normres<-normres-min(normres)+1
tsvr::vr(t(normres))

tot<-apply(FUN=sum,MARGIN=1,X=normres)
hist(tot,50)
#source("SkewnessAnd3CentMom.R")
myskns(tot)

tim<-1:50
plot(tim,res[tim,1],type="l",col="red")
lines(tim,res[tim,2],type="l",col="pink")
lines(tim,res[tim,3],type="l",col="blue")
lines(tim,res[tim,4],type="l",col="black")

#try to generalize it some

getdat4<-function(n,sig)
{
  coinflip<-runif(n)
  cop<-copula::normalCopula(param=P2p(sig),dim=dim(sig)[1],dispstr="un")
  intres<-rCopula(n,cop)
  d<-dim(sig)[1]-1

  res<-matrix(NA,n,2*d)
  #for about half the cases 
  res[coinflip<1/2,1:d]<-rep(intres[coinflip<1/2,1]/2+1/2,times=d)
  res[coinflip<1/2,(d+1):(2*d)]<-intres[coinflip<1/2,2:(d+1)]/2
  
  #for the other half
  res[coinflip>=1/2,1:d]<-intres[coinflip>=1/2,2:(d+1)]/2
  res[coinflip>=1/2,(d+1):(2*d)]<-rep(intres[coinflip>=1/2,1]/2+1/2,times=d)
  
  return(res)
}

#Skewness ratio. Assumes no NAs.
#
#Args
#x      A matrix, years by species
#
sr<-function(x)
{
  xtot<-apply(FUN=sum,X=x,MARGIN=1)
  
  scom<-myskns(xtot)
  sind<-(sum(apply(FUN=my3cm,MARGIN=2,X=x)))/(sum(apply(FUN=var,MARGIN=2,X=x)))^(3/2)
  
  return(scom/sind)
}

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
res<-getdat4(n=10000,sig=sig)

hist(res[,1],30)
hist(res[,2],30)
hist(res[,3],30)
hist(res[,4],30) 
hist(res[,5],30) 
hist(res[,6],30) 
hist(res[,7],30) 
hist(res[,8],30) 
hist(res[,9],30) 
hist(res[,10],30) #histograms seem uniform, as they should be

plot(res[1:500,1],res[1:500,2],type="p",pch=20)
plot(res[1:500,6],res[1:500,7],type="p",pch=20)

plot(res[1:500,1],res[1:500,6],type="p",pch=20)
plot(res[1:500,1],res[1:500,7],type="p",pch=20)
plot(res[1:500,2],res[1:500,6],type="p",pch=20)
plot(res[1:500,2],res[1:500,7],type="p",pch=20)

normres<-qnorm(res)
normres<-normres-min(normres)+1
tsvr::vr(t(normres))
tot<-apply(FUN=sum,MARGIN=1,X=normres)
hist(tot,50)
myskns(tot)
sr(normres)

lognormres<-exp(qnorm(res))
tsvr::vr(t(lognormres))
tot<-apply(FUN=sum,MARGIN=1,X=lognormres)
hist(tot,50)
myskns(tot)

gamres<-qgamma(res,shape=2,scale=2)
tsvr::vr(t(gamres))
tot<-apply(FUN=sum,MARGIN=1,X=gamres)
hist(tot,50)
myskns(tot)




tim<-1:50
plot(tim,res[tim,1],type="l",col="red")
for (counter in 2:d)
{
  lines(tim,res[tim,counter],col="red")
}
for (counter in (d+1):(2*d))
{
  lines(tim,res[tim,counter],col="blue")
}
