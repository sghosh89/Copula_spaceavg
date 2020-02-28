#***the data to work on

d<-readRDS(file="./Results/hays_results/skewness_results/ts_mat_CP_hays.RDS")

#***setup for the chunk

library(matrixcalc)
library(mvtnorm)

# functions needed
source("PPSurrogObjFun.R")
source("pwlin.R")
source("getmap.R")
source("alignranks.R")

# result folder to save surrogates and function plots
resloc_surrog_hays<-"./Results/hays_results/skewness_results/pp_surrogs_hays_CP/"
if (!dir.exists(resloc_surrog_hays)){
  dir.create(resloc_surrog_hays)
}

#***for each pair of species, create the function from N-copula space
#to covariance space

set.seed(102)
numpts<-18
mapy<-array(NA,c(dim(d)[2],dim(d)[2],numpts))
for (iind in 1:(dim(d)[2]-1))
{
  for (jind in (iind+1):(dim(d)[2]))
  {
    print(paste0("iind: ",iind,"; jind: ",jind))
    plotnm<-paste0(resloc_surrog_hays,"Map_",iind,"_",jind)
    thisres<-getmap(d[,c(iind,jind)],numpts=numpts,numsims=500,plotnm=plotnm)
    mapx<-thisres$fit_parameters$x
    y<-thisres$fit_parameters$y
    if (any(diff(y)<0)){stop("Error in make_surrogs_CP_hays: non-monotonic piecewise linear interpolation")}
    mapy[iind,jind,]<-y
  }
}
save(mapx,mapy,file=paste0(resloc_surrog_hays,"MapsRes.RData"))

#***now do a constrained optimization

#first just test things by calling the objective function on the start vector
pijstart<-cor(d)
pijstart<-pijstart[upper.tri(pijstart)]
PPSurrogObj(pij=pijstart,cijmatd=cijmatd,mapx=mapx,mapy=mapy,saveloc=resloc_surrog_hays)

#set up the constraints - see the math saved in PPSurrog202001
#getting ready
cijmatd<-cov(d)
n<-dim(mapy)[1] #number of species
epsil<-0.1
yijbot<-mapy[,,1]
yijbot<-yijbot[upper.tri(yijbot)]
yijtop<-mapy[,,dim(mapy)[3]]
yijtop<-yijtop[upper.tri(yijtop)]
#first set of constraints
ui1<-diag(length(pijstart))
ci1<-yijbot
#second set of constraints
ui2<-(-1*diag(length(pijstart)))
ci2<-(-yijtop)
#third set of constraints (actually just one constraint)
h<-outer(diag(cijmatd),diag(cijmatd))
h<-sqrt(h[upper.tri(h)])
ui3<-matrix(h,1,length(h))
synccomp<-sum(cijmatd[upper.tri(cijmatd)])
if (synccomp>=0)
{
  ci3<-(1-epsil)*synccomp
} else
{
  ci3<-(1+epsil)*synccomp
}
#fourth set of constraints
ui4<-(-ui3)
if (synccomp>=0)
{
  ci4<-(-(1+epsil))*synccomp
} else
{
  ci4<-(-(1-epsil))*synccomp
}
#combine them all
ui<-rbind(ui1,ui2,ui3,ui4)
ci<-c(ci1,ci2,ci3,ci4)

#make sure the start vector satisfies the constrain
sum(ui%*%pijstart-ci>=0)
length(ci)

#Now just explore the objective function a bit, by looking at slices
for (indcounter in 1:length(pijstart))
{
  print(paste0("indcounter ",indcounter," of ",length(pijstart)))
  bdbot1<-yijbot[indcounter]
  h<-pijstart*ui3
  if (synccomp>=0)
  {
    bdbot2<-((1-epsil)*sum(cijmatd[upper.tri(cijmatd)])-sum(h[-indcounter]))/(ui3[1,indcounter])
  } else
  {
    bdbot2<-((1+epsil)*sum(cijmatd[upper.tri(cijmatd)])-sum(h[-indcounter]))/(ui3[1,indcounter])
  }
  bdbot<-max(bdbot1,bdbot2)
  bdtop1<-yijtop[indcounter]
  if (synccomp>=0)
  {
    bdtop2<-((1+epsil)*sum(cijmatd[upper.tri(cijmatd)])-sum(h[-indcounter]))/(ui3[1,indcounter])
  } else
  {
    bdtop2<-((1-epsil)*sum(cijmatd[upper.tri(cijmatd)])-sum(h[-indcounter]))/(ui3[1,indcounter])
  }
  bdtop<-min(bdtop1,bdtop2)
  x<-seq(from=bdbot,to=bdtop,length.out=100)
  x<-sort(cbind(x,pijstart[indcounter]))
  thisres<-NA*numeric(length(x))
  pij<-pijstart
  for (xcounter in 1:length(x))
  {
    pij[indcounter]<-x[xcounter]
    thisres[xcounter]<-PPSurrogObj(pij=pij,cijmatd=cijmatd,mapx=mapx,mapy=mapy,saveloc=resloc_surrog_hays)
  }
  pdf(paste0(resloc_surrog_hays,"Slice_",indcounter,".pdf"))
  plot(x,thisres,type="b",pch=20,cex=.5)  
  points(x[x==pijstart[indcounter]],thisres[x==pijstart[indcounter]],col="red")
  lines(rep(pijstart[indcounter],2),c(-1000,1000),type="l",col="red")
  dev.off()
}

#Now do the optimization. The goal is not actually to get all the way to
#the optimum, rather it is to get to the first pos def mat you come to.
#The try function is used for that purpose.

#now run the optimization
res<-NA
try(res<-constrOptim(theta=pijstart,f=PPSurrogObj,grad=NULL,ui=ui,ci=ci,
                     cijmatd=cijmatd,mapx=mapx,mapy=mapy,saveloc=resloc_surrog_hays,
                     control=list(fnscale=-1,trace=0,maxit=10000)),silent=TRUE)

#***now if it succeeded in finding a pos def mat, then it saved a result, and 
#we can pull the result back in from the disk and examine it

if (!is.na(res))
{
  stop("Error in make_surrogs_CP_hays: optimization to find good Pearson preserving surrogates failed")
}

pijmat<-readRDS(paste0(resloc_surrog_hays,"FirstSuccess_pijmat.RDS"))
nijmat<-readRDS(paste0(resloc_surrog_hays,"FirstSuccess_nijmat.RDS"))
pij<-pijmat[upper.tri(pijmat)]
pdf(file=paste0(resloc_surrog_hays,"ActualVsOptimizedSurrogatePearsons.pdf"))
plot(pijstart,pij,type='p',xlab="vij",ylab="pij",xlim=c(-1,1),ylim=c(-1,1))
lines(c(-1,1),c(-1,1))
text(1,1,cor(pijstart,pij),adj=c(1,1))
dev.off()

#make surrogates using nijmat, and compute their CV_com^2 and 
#compare to the actual CV_com^2 of the real data 

numsurrog<-100000
sims<-rmvnorm(numsurrog*(dim(d)[1]),sigma=nijmat)
dim(sims)<-c(dim(d)[1],numsurrog,dim(d)[2])
sims<-aperm(sims,c(1,3,2))
dsort<-apply(FUN=sort,X=d,MARGIN=2)
surrogs<-alignranks(dsort,sims)
saveRDS(surrogs,file=paste0(resloc_surrog_hays,"HaysSurrogates.RDS"))

totpops<-apply(FUN=sum,MARGIN=c(1,3),X=surrogs)
survars<-apply(FUN=var,X=totpops,MARGIN=2)
surmns<-apply(FUN=mean,X=totpops,MARGIN=2)
CVcom2_sur<-survars/(surmns^2)

totd<-apply(FUN=sum,X=d,MARGIN=1)
dvar<-var(totd)
dmn<-mean(totd)
CVcom2_d<-dvar/(dmn^2)

pdf(file=paste0(resloc_surrog_hays,"CVcom2_For_Surr_And_Dat.pdf"))
hist(CVcom2_sur,main=sum(CVcom2_sur<CVcom2_d)/numsurrog)
points(CVcom2_d,0,col="red")
dev.off()

surrogs_CP_hays<-surrogs  # this is the pearson preserving surrogs for 
#use in subsequent chunks