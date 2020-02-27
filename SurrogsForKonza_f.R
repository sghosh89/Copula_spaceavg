###***Note from Shya***
#Two different soil types: f, t.
#Each common species time series along each column of the matrix except for the 
#last column which indicates the intermediate and rare ones merged into a pseudo sp.
#
#Matrix data saved as .RDS from this file (Line 257)  with a given soil type (at Line 164)   
#https://github.com/sghosh89/Copula_spaceavg/blob/master/data_cleaning_and_npa_for_KNZ.R
#
#Thanks,
#Shyamolina.
#***End note from Shya***


#***the data to work on

datloc_knz<-"./PPSurrog202001/"
d<-readRDS(file=paste0(datloc_knz,"ts_CP_knz_soiltype_f.RDS"))

#***setup for the chunk

library(matrixcalc)
library(mvtnorm)

# functions needed
source("PPSurrogObjFun.R")
source("pwlin.R")
source("getmap.R")
source("alignranks.R")

# result folder to save surrogates and function plots
resloc_surrog_knz<-"./Results/knz_results/skewness_results/pp_surrogs_knz_f_CP/"
if (!dir.exists(resloc_surrog_knz)){
  dir.create(resloc_surrog_knz)
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
    plotnm<-paste0(resloc_surrog_knz,"Map_",iind,"_",jind)
    thisres<-getmap(d[,c(iind,jind)],numpts=numpts,numsims=500,plotnm=plotnm)
    mapx<-thisres$fit_parameters$x
    y<-thisres$fit_parameters$y
    if (any(diff(y)<0)){stop("Error in make_surrogs_CP_knz: non-monotonic piecewise linear interpolation")}
    mapy[iind,jind,]<-y
  }
}
save(mapx,mapy,file=paste0(resloc_surrog_knz,"MapsRes.RData"))

#***now do a constrained optimization

#first just test things by calling the objective function on the start vector
pijstart<-cor(d)
pijstart<-pijstart[upper.tri(pijstart)]
PPSurrogObj(pij=pijstart,cijmatd=cijmatd,mapx=mapx,mapy=mapy,saveloc=resloc_surrog_knz)

#set up the constraints - see the math saved in PPSurrog202001
#getting ready
cijmatd<-cov(d)
n<-dim(mapy)[1] #number of species
epsil<-0.15
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
    thisres[xcounter]<-PPSurrogObj(pij=pij,cijmatd=cijmatd,mapx=mapx,mapy=mapy,saveloc=resloc_surrog_knz)
  }
  pdf(paste0(resloc_surrog_knz,"Slice_",indcounter,".pdf"))
  plot(x,thisres,type="b",pch=20,cex=.5)  
  points(x[x==pijstart[indcounter]],thisres[x==pijstart[indcounter]],col="red")
  lines(rep(pijstart[indcounter],2),c(-1000,1000),type="l",col="red")
  dev.off()
}

save.image(paste0(resloc_surrog_knz,"SaveBeforeOptimize.RData"))
#resloc_surrog_knz<-"./Results/knz_results/skewness_results/pp_surrogs_knz_f_CP/"
#load(paste0(resloc_surrog_knz,"SaveBeforeOptimize.RData"))

#Now do the optimization. The goal is not actually to get all the way to
#the optimum, rather it is to get to the first pos def mat you come to.
#The try function is used for that purpose.

#now run the optimization
#res<-NA
#try(
res1<-constrOptim(theta=pijstart,f=PPSurrogObj,grad=NULL,ui=ui,ci=ci,
                     cijmatd=cijmatd,mapx=mapx,mapy=mapy,saveloc=resloc_surrog_knz,
                     control=list(fnscale=-1,trace=0,maxit=2500))
#  ,silent=TRUE)

PPSurrogObj(pij=res1$par,cijmatd=cijmatd,mapx=mapx,mapy=mapy,saveloc=resloc_surrog_knz)
res1$value
sum(ui%*%res1$par-ci>=0)
length(ci)
sum(ui%*%res1$par-ci>0)

res2<-constrOptim(theta=res1$par,f=PPSurrogObj,grad=NULL,ui=ui,ci=ci,
                 cijmatd=cijmatd,mapx=mapx,mapy=mapy,saveloc=resloc_surrog_knz,
                 control=list(fnscale=-1,trace=0,maxit=2500))

PPSurrogObj(pij=res2$par,cijmatd=cijmatd,mapx=mapx,mapy=mapy,saveloc=resloc_surrog_knz)
res2$value
sum(ui%*%res2$par-ci>=0)
length(ci)
sum(ui%*%res2$par-ci>0)

res3<-constrOptim(theta=res2$par,f=PPSurrogObj,grad=NULL,ui=ui,ci=ci,
                  cijmatd=cijmatd,mapx=mapx,mapy=mapy,saveloc=resloc_surrog_knz,
                  control=list(fnscale=-1,trace=0,maxit=2500))

PPSurrogObj(pij=res3$par,cijmatd=cijmatd,mapx=mapx,mapy=mapy,saveloc=resloc_surrog_knz)
res3$value
sum(ui%*%res3$par-ci>=0)
length(ci)
sum(ui%*%res3$par-ci>0)

#It did not improve. Not sure what else to do so giving up on the knz_f data for now
#and going to the knx_t data

