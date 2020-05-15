#This is for making a figure to help introducte the theoretical ideas of influence
#of tail associations on variance of a sum and the limits of that.

rm(list=ls())
graphics.off()
set.seed(101)
library(copula)
source("SkewnessAnd3CentMom.R")

#***functions

#Generates points from certain multivariate copulas with tail dependence
#
#Args
#corpar           For the part of the distribution where the variable are
#                   not perfectly related they are generated with a normal
#                   copula with this parameter in all off-diagonal locations
#beta1, beta2     These specify where the tail dependence starts in each tail.
#                   Numbers between 0 and 1 with beta1<beta2.
#numpts           Number of draws from the copula desired
#n                Dimension
#
#Output - a matrix of dimensions numpts by 2
#
taildepcopdat<-function(corpar,beta1,beta2,numpts,n=2)
{
  corpar<-matrix(corpar,n,n)
  diag(corpar)<-1
  corpar<-P2p(corpar)
  ncop<-normalCopula(corpar,n,"un")
  res<-rCopula(numpts,ncop)
  tdornot<-runif(numpts)
  inds<-(tdornot<beta1 | tdornot>beta2)
  res[inds,]<-tdornot[inds]
  res[!inds,]<-(beta2-beta1)*res[!inds,]+beta1
  return(res)
}

#***tests

res<-taildepcopdat(.5,.2,.8,1000)
plot(res[,1],res[,2],type="p",pch=20,cex=.5)
res<-taildepcopdat(.9,.3,.9,1000)
plot(res[,1],res[,2],type="p",pch=20,cex=.5)
res<-taildepcopdat(.8,.1,.9,1000,3)
plot(res[,1],res[,2],type="p",pch=20,cex=.5)
plot(res[,1],res[,3],type="p",pch=20,cex=.5)
plot(res[,2],res[,3],type="p",pch=20,cex=.5)

#***Get the data you need for plotting

numpts<-1000000
numscatpts<-500
dtouse<-20

#Functions which control the marginals
#margfunc<-function(c){qgamma(c,shape=2,scale=2)}
#margfunc<-function(c){qgamma(c,shape=7.5,scale=1)}
#margfunc<-function(c){qgamma(c,shape=9,scale=.5)}
margfunc<-function(c){qnorm(c)}

#get data for example 1 - normal copula case
dcop1<-taildepcopdat(.6,0,1,numpts,dtouse)
d1<-margfunc(dcop1)
cov1<-cov(d1)
diag(cov1)<-NA
cov1<-mean(cov1,na.rm=TRUE)
totX1<-apply(FUN=sum,X=d1,MARGIN=1)
vartotX1<-var(totX1)
sktotX1<-myskns(totX1)
empcdf1<-data.frame(x=sort(totX1),y=(1:numpts)/numpts)

#get data for example 2 - comonotonic
dcop2<-matrix(rep(runif(numpts),times=dtouse),numpts,dtouse)
d2<-margfunc(dcop2)
cov2<-cov(d2)
diag(cov2)<-NA
cov2<-mean(cov2,na.rm=TRUE)
totX2<-apply(FUN=sum,X=d2,MARGIN=1)
vartotX2<-var(totX2)
sktotX2<-myskns(totX2)
empcdf2<-data.frame(x=sort(totX2),y=(1:numpts)/numpts)

#get data for example 3 - upper-tail comonotonic
beta1<-0
beta2<-0.95
objfun<-function(rho)
{
  dat<-margfunc(taildepcopdat(rho,beta1,beta2,numpts))
  h<-cov(dat)[1,2]
  return(c(h,(h-cov1)^2))
}
x<-seq(from=.01,to=0.99,by=0.01)
ycov<-NA*numeric(length(x))
yobj<-NA*numeric(length(x))
for (counter in 1:length(x))
{
  h<-objfun(x[counter])
  ycov[counter]<-h[1]
  yobj[counter]<-h[2]
}
plot(x,yobj,type="l")
lines(range(x),rep(0,2),type="l")
plot(x,ycov,type="l")
lines(range(x),rep(cov1,2),type="l")
mod<-lm(ycov~x)
lines(range(x),coef(mod)[1]+coef(mod)[2]*range(x),type="l")
parmtouse<-(cov1-coef(mod)[1])/coef(mod)[2]
lines(rep(parmtouse,2),range(ycov),type="l")
dcop3<-taildepcopdat(parmtouse,beta1,beta2,numpts,dtouse)
d3<-margfunc(dcop3)
cov3<-cov(d3)
diag(cov3)<-NA
cov3<-mean(cov3,na.rm=TRUE)
totX3<-apply(FUN=sum,X=d3,MARGIN=1)
vartotX3<-var(totX3)
sktotX3<-myskns(totX3)
empcdf3<-data.frame(x=sort(totX3),y=(1:numpts)/numpts)

#get data for example 4 - lower-tail comonotonic
beta1<-0.05
beta2<-1
x<-seq(from=.01,to=0.99,by=0.01)
ycov<-NA*numeric(length(x))
yobj<-NA*numeric(length(x))
for (counter in 1:length(x))
{
  h<-objfun(x[counter])
  ycov[counter]<-h[1]
  yobj[counter]<-h[2]
}
plot(x,yobj,type="l")
lines(range(x),rep(0,2),type="l")
plot(x,ycov,type="l")
lines(range(x),rep(cov1,2),type="l")
mod<-lm(ycov~x)
lines(range(x),coef(mod)[1]+coef(mod)[2]*range(x),type="l")
parmtouse<-(cov1-coef(mod)[1])/coef(mod)[2]
lines(rep(parmtouse,2),range(ycov),type="l")
dcop4<-taildepcopdat(parmtouse,beta1,beta2,numpts,dtouse)
d4<-margfunc(dcop4)
cov4<-cov(d4)
diag(cov4)<-NA
cov4<-mean(cov4,na.rm=TRUE)
totX4<-apply(FUN=sum,X=d4,MARGIN=1)
vartotX4<-var(totX4)
sktotX4<-myskns(totX4)
empcdf4<-data.frame(x=sort(totX4),y=(1:numpts)/numpts)

#get data for example 5
beta1<-0.05
beta2<-0.95
x<-seq(from=.01,to=0.99,by=0.01)
ycov<-NA*numeric(length(x))
yobj<-NA*numeric(length(x))
for (counter in 1:length(x))
{
  h<-objfun(x[counter])
  ycov[counter]<-h[1]
  yobj[counter]<-h[2]
}
plot(x,yobj,type="l")
lines(range(x),rep(0,2),type="l")
plot(x,ycov,type="l")
lines(range(x),rep(cov1,2),type="l")
mod<-lm(ycov~x)
lines(range(x),coef(mod)[1]+coef(mod)[2]*range(x),type="l")
parmtouse<-(cov1-coef(mod)[1])/coef(mod)[2]
lines(rep(parmtouse,2),range(ycov),type="l")
dcop5<-taildepcopdat(parmtouse,beta1,beta2,numpts,dtouse)
d5<-margfunc(dcop5)
cov5<-cov(d5)
diag(cov5)<-NA
cov5<-mean(cov5,na.rm=TRUE)
totX5<-apply(FUN=sum,X=d5,MARGIN=1)
vartotX5<-var(totX5)
sktotX5<-myskns(totX5)
empcdf5<-data.frame(x=sort(totX5),y=(1:numpts)/numpts)

#other prep for consistency among plots
scatylim<-range(d1[,2],d2[,2],d3[,2],d4[,2],d5[,2],na.rm=TRUE)
scatxlim<-range(d1[,1],d2[,1],d3[,1],d4[,1],d5[,1],na.rm=TRUE)
histXbreaks<-seq(from=scatxlim[1],to=scatxlim[2],length.out=50)
histYbreaks<-seq(from=scatylim[1],to=scatylim[2],length.out=50)
histtotXxlim<-range(totX1,totX2,totX3,totX4,totX5,na.rm=TRUE)
histtotXbreaks<-seq(from=histtotXxlim[1],to=histtotXxlim[2],length.out=50)
topthreshrg<-c(20,60)
botthreshrg<-c(-60,-20)

#***plot dimensions, units inches
xmarght<-.5
ymargwd<-.6
totwd<-6
gap<-0.25
spf<-0.45
panwd<-(totwd-2*gap-3*ymargwd)/(3+spf)
smallpan<-panwd*spf
panht<-panwd
cdfpanht<-(panht+gap+smallpan-gap)/2
totht<-xmarght+5*(panht+2*gap+smallpan)
textsz<-.9
pdf(file="TheoryFig.pdf",width=totwd,height=totht)

#***Example 1, top row of panels, normal copula

#scatterplot
panrownum<-5
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
plot(d1[1:numscatpts,1],d1[1:numscatpts,2],xaxt='n',pch=20,cex=.5,col=rgb(red=190/256,green=190/256,blue=190/256,alpha=.3),
     xlim=scatxlim,ylim=scatylim)
axis(side=1,labels=FALSE)
mtext(expression(X[2]),side=2,line=1.2)

#X1 marginal of scatterplot
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap)/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap+smallpan)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d1[,1],breaks=histXbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',yaxt='n',type="l")
axis(side=1,labels=FALSE)
axis(side=2,labels=TRUE)
mtext("Ct.",side=2,line=1.2)

#X2 marginal of scatterplot
par(fig=c((ymargwd+panwd+gap)/totwd,
          (ymargwd+panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d1[,2],breaks=histYbreaks,plot=FALSE)
y<-h$breaks
y<-y[2:length(y)]-diff(y)[1]/2
plot(h$counts,y,yaxt='n',xaxt='n',type="l")
#axis(side=1,at=c(0,400),labels=c("0","400"))
axis(side=1,labels=TRUE)
axis(side=2,labels=FALSE)

#Histogram of totX
par(fig=c((2*ymargwd+panwd+gap+smallpan)/totwd,
          (2*ymargwd+2*panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h1<-hist(totX1,breaks=histtotXbreaks,plot=FALSE)
x1<-h1$breaks
x1<-x1[2:length(x1)]-diff(x1)[1]/2
plot(x1,h1$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Count",side=2,line=1.2)
mtext(bquote(cov(X[i],X[j]) %~~% .(round(cov1,3))),side=3,line=1.6,cex=textsz)
mtext(bquote(sd(Sigma[i] * X[i]) == .(round(sqrt(vartotX1),3))),side=3,line=.8,cex=textsz)
mtext(bquote(sk(Sigma[i] * X[i]) == .(round(sktotX1,3))),side=3,line=0,cex=textsz)

#cdf panel for exceeding large thresholds (top panel of the cdf panels)
par(fig=c((2*ymargwd+2*panwd+gap+smallpan+ymargwd)/totwd,
          (2*ymargwd+3*panwd+gap+smallpan+ymargwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+cdfpanht+gap)/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+2*cdfpanht+gap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
inds<-which(empcdf1$x>topthreshrg[1] & empcdf1$x<topthreshrg[2])
x<-empcdf1$x[inds]
y<-1-empcdf1$y[inds]
plot(x,y,type="l")
mtext("P>th.",side=2,line=1.2)

#cdf panel for falling under small thresholds (bottom panel of the cdf panels)
par(fig=c((2*ymargwd+2*panwd+gap+smallpan+ymargwd)/totwd,
          (2*ymargwd+3*panwd+gap+smallpan+ymargwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+cdfpanht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
inds<-which(empcdf1$x>botthreshrg[1] & empcdf1$x<botthreshrg[2])
x<-empcdf1$x[inds]
y<-empcdf1$y[inds]
plot(x,y,type="l")
mtext("P<th.",side=2,line=1.2)

#***Example 2, second row of panels, comonotonic case

#scatterplot
panrownum<-4
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d2[1:numscatpts,1],d2[1:numscatpts,2],xaxt='n',pch=20,cex=.5,col=rgb(red=255/256,green=0/256,blue=0/256,alpha=.3),
     xlim=scatxlim,ylim=scatylim)
axis(side=1,labels=FALSE)
mtext(expression(X[2]),side=2,line=1.2)

#X1 marginal of scatterplot
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap)/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap+smallpan)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d2[,1],breaks=histXbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Ct.",side=2,line=1.2)

#X2 marginal of scatterplot
par(fig=c((ymargwd+panwd+gap)/totwd,
          (ymargwd+panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d2[,2],breaks=histYbreaks,plot=FALSE)
y<-h$breaks
y<-y[2:length(y)]-diff(y)[1]/2
plot(h$counts,y,yaxt='n',xaxt='n',type="l")
axis(side=1,labels=TRUE)
axis(side=2,labels=FALSE)

#Histogram of totX
par(fig=c((2*ymargwd+panwd+gap+smallpan)/totwd,
          (2*ymargwd+2*panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(totX2,breaks=histtotXbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l",col="red",ylim=range(h$counts,h1$counts))
lines(x1,h1$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Count",side=2,line=1.2)
mtext(bquote(cov(X[i],X[j]) %~~% .(round(cov2,3))),side=3,line=1.6,cex=textsz)
mtext(bquote(sd(Sigma[i] * X[i]) == .(round(sqrt(vartotX2),3))),side=3,line=.8,cex=textsz)
mtext(bquote(sk(Sigma[i] * X[i]) == .(round(sktotX2,3))),side=3,line=0,cex=textsz)

#cdf panel for exceeding large thresholds (top panel of the cdf panels)
par(fig=c((2*ymargwd+2*panwd+gap+smallpan+ymargwd)/totwd,
          (2*ymargwd+3*panwd+gap+smallpan+ymargwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+cdfpanht+gap)/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+2*cdfpanht+gap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
inds<-which(empcdf1$x>topthreshrg[1] & empcdf1$x<topthreshrg[2])
x<-empcdf1$x[inds]
y<-1-empcdf1$y[inds]
inds<-which(empcdf2$x>topthreshrg[1] & empcdf2$x<topthreshrg[2])
xn<-empcdf2$x[inds]
yn<-1-empcdf2$y[inds]
plot(x,y,type="l",lty="dashed",ylim=c(0,max(y,yn)))
lines(xn,yn,type='l',col='red')
mtext("P>th.",side=2,line=1.2)

#cdf panel for falling under small thresholds (bottom panel of the cdf panels)
par(fig=c((2*ymargwd+2*panwd+gap+smallpan+ymargwd)/totwd,
          (2*ymargwd+3*panwd+gap+smallpan+ymargwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+cdfpanht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
inds<-which(empcdf1$x>botthreshrg[1] & empcdf1$x<botthreshrg[2])
x<-empcdf1$x[inds]
y<-empcdf1$y[inds]
inds<-which(empcdf2$x>botthreshrg[1] & empcdf2$x<botthreshrg[2])
xn<-empcdf2$x[inds]
yn<-empcdf2$y[inds]
plot(x,y,type="l",lty="dashed",ylim=c(0,max(y,yn)))
lines(xn,yn,type='l',col='red')
mtext("P<th.",side=2,line=1.2)

#***Example 3, third row of panels, comonotonicity in the right tails only

#scatterplot
panrownum<-3
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d3[1:numscatpts,1],d3[1:numscatpts,2],xaxt='n',pch=20,cex=.5,col=rgb(red=0/256,green=0/256,blue=255/256,alpha=.3),
     xlim=scatxlim,ylim=scatylim)
axis(side=1,labels=FALSE)
mtext(expression(X[2]),side=2,line=1.2)

#X1 marginal of scatterplot
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap)/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap+smallpan)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d3[,1],breaks=histXbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Ct.",side=2,line=1.2)

#X2 marginal of scatterplot
par(fig=c((ymargwd+panwd+gap)/totwd,
          (ymargwd+panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d3[,2],breaks=histYbreaks,plot=FALSE)
y<-h$breaks
y<-y[2:length(y)]-diff(y)[1]/2
plot(h$counts,y,yaxt='n',xaxt='n',type="l")
axis(side=1,labels=TRUE)
axis(side=2,labels=FALSE)

#Histogram of totX
par(fig=c((2*ymargwd+panwd+gap+smallpan)/totwd,
          (2*ymargwd+2*panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(totX3,breaks=histtotXbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l",col="blue",ylim=range(h$counts,h1$counts))
lines(x1,h1$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Count",side=2,line=1.2)
mtext(bquote(cov(X[i],X[j]) %~~% .(round(cov3,3))),side=3,line=1.6,cex=textsz)
mtext(bquote(sd(Sigma[i] * X[i]) == .(round(sqrt(vartotX3),3))),side=3,line=.8,cex=textsz)
mtext(bquote(sk(Sigma[i] * X[i]) == .(round(sktotX3,3))),side=3,line=0,cex=textsz)

#cdf panel for exceeding large thresholds (top panel of the cdf panels)
par(fig=c((2*ymargwd+2*panwd+gap+smallpan+ymargwd)/totwd,
          (2*ymargwd+3*panwd+gap+smallpan+ymargwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+cdfpanht+gap)/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+2*cdfpanht+gap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
inds<-which(empcdf1$x>topthreshrg[1] & empcdf1$x<topthreshrg[2])
x<-empcdf1$x[inds]
y<-1-empcdf1$y[inds]
inds<-which(empcdf2$x>topthreshrg[1] & empcdf2$x<topthreshrg[2])
xr_l<-empcdf2$x[inds]
yr_l<-1-empcdf2$y[inds]
inds<-which(empcdf3$x>topthreshrg[1] & empcdf3$x<topthreshrg[2])
xn<-empcdf3$x[inds]
yn<-1-empcdf3$y[inds]
plot(x,y,type="l",lty="dashed",ylim=c(0,max(y,yn,yr_l)))
lines(xr_l,yr_l,type="l",lty="dashed",col="red")
lines(xn,yn,type='l',col='blue')
mtext("P>th.",side=2,line=1.2)

#cdf panel for falling under small thresholds (bottom panel of the cdf panels)
par(fig=c((2*ymargwd+2*panwd+gap+smallpan+ymargwd)/totwd,
          (2*ymargwd+3*panwd+gap+smallpan+ymargwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+cdfpanht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
inds<-which(empcdf1$x>botthreshrg[1] & empcdf1$x<botthreshrg[2])
x<-empcdf1$x[inds]
y<-empcdf1$y[inds]
inds<-which(empcdf2$x>botthreshrg[1] & empcdf2$x<botthreshrg[2])
xr_s<-empcdf2$x[inds]
yr_s<-empcdf2$y[inds]
inds<-which(empcdf3$x>botthreshrg[1] & empcdf3$x<botthreshrg[2])
xn<-empcdf3$x[inds]
yn<-empcdf3$y[inds]
plot(x,y,type="l",lty="dashed",ylim=c(0,max(y,yn,yr_s)))
lines(xr_s,yr_s,type="l",lty="dashed",col="red")
lines(xn,yn,type='l',col='blue')
mtext("P<th.",side=2,line=1.2)

#***Example 4, fourth row of panels, comonotonicity in the left tails only

#scatterplot
panrownum<-2
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d4[1:numscatpts,1],d4[1:numscatpts,2],xaxt='n',pch=20,cex=.5,col=rgb(red=0/256,green=255/256,blue=0/256,alpha=.3),
     xlim=scatxlim,ylim=scatylim)
axis(side=1,labels=FALSE)
mtext(expression(X[2]),side=2,line=1.2)

#X1 marginal of scatterplot
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap)/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap+smallpan)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d4[,1],breaks=histXbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Ct.",side=2,line=1.2)

#X2 marginal of scatterplot
par(fig=c((ymargwd+panwd+gap)/totwd,
          (ymargwd+panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d4[,2],breaks=histYbreaks,plot=FALSE)
y<-h$breaks
y<-y[2:length(y)]-diff(y)[1]/2
plot(h$counts,y,yaxt='n',xaxt='n',type="l")
axis(side=1,labels=TRUE)
axis(side=2,labels=FALSE)

#Histogram of X+Y
par(fig=c((2*ymargwd+panwd+gap+smallpan)/totwd,
          (2*ymargwd+2*panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(totX4,breaks=histtotXbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l",col="green",ylim=range(h$counts,h1$counts))
lines(x1,h1$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Count",side=2,line=1.2)
mtext(bquote(cov(X[i],X[j]) %~~% .(round(cov4,3))),side=3,line=1.6,cex=textsz)
mtext(bquote(sd(Sigma[i] * X[i]) == .(round(sqrt(vartotX4),3))),side=3,line=.8,cex=textsz)
mtext(bquote(sk(Sigma[i] * X[i]) == .(round(sktotX4,3))),side=3,line=0,cex=textsz)

#cdf panel for exceeding large thresholds (top panel of the cdf panels)
par(fig=c((2*ymargwd+2*panwd+gap+smallpan+ymargwd)/totwd,
          (2*ymargwd+3*panwd+gap+smallpan+ymargwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+cdfpanht+gap)/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+2*cdfpanht+gap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
inds<-which(empcdf1$x>topthreshrg[1] & empcdf1$x<topthreshrg[2])
x<-empcdf1$x[inds]
y<-1-empcdf1$y[inds]
inds<-which(empcdf4$x>topthreshrg[1] & empcdf4$x<topthreshrg[2])
xn<-empcdf4$x[inds]
yn<-1-empcdf4$y[inds]
plot(x,y,type="l",lty="dashed",ylim=c(0,max(y,yn)))
lines(xr_l,yr_l,type="l",lty="dashed",col="red")
lines(xn,yn,type='l',col='green')
mtext("P>th.",side=2,line=1.2)

#cdf panel for falling under small thresholds (bottom panel of the cdf panels)
par(fig=c((2*ymargwd+2*panwd+gap+smallpan+ymargwd)/totwd,
          (2*ymargwd+3*panwd+gap+smallpan+ymargwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+cdfpanht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
inds<-which(empcdf1$x>botthreshrg[1] & empcdf1$x<botthreshrg[2])
x<-empcdf1$x[inds]
y<-empcdf1$y[inds]
inds<-which(empcdf4$x>botthreshrg[1] & empcdf4$x<botthreshrg[2])
xn<-empcdf4$x[inds]
yn<-empcdf4$y[inds]
plot(x,y,type="l",lty="dashed",ylim=c(0,max(y,yn)))
lines(xr_s,yr_s,type="l",lty="dashed",col="red")
lines(xn,yn,type='l',col='green')
mtext("P<th.",side=2,line=1.2)

#***Example 5, fifth row of panels, comonotonicity in the right and left tails only

#scatterplot
panrownum<-1
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d5[1:numscatpts,1],d5[1:numscatpts,2],xaxt='n',pch=20,cex=.5,col=rgb(red=160/256,green=32/256,blue=240/256,alpha=.3),
     xlim=scatxlim,ylim=scatylim)
axis(side=1,labels=TRUE)
mtext(expression(X[2]),side=2,line=1.2)
mtext(expression(X[1]),side=1,line=1.2)

#X1 marginal of scatterplot
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap)/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap+smallpan)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d5[,1],breaks=histXbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Ct.",side=2,line=1.2)

#X2 marginal of scatterplot
par(fig=c((ymargwd+panwd+gap)/totwd,
          (ymargwd+panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d5[,2],breaks=histYbreaks,plot=FALSE)
y<-h$breaks
y<-y[2:length(y)]-diff(y)[1]/2
plot(h$counts,y,yaxt='n',xaxt='n',type="l")
axis(side=1,labels=TRUE)
axis(side=2,labels=FALSE)
mtext("Ct.",side=1,line=1.2)

#Histogram of X+Y
par(fig=c((2*ymargwd+panwd+gap+smallpan)/totwd,
          (2*ymargwd+2*panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(totX5,breaks=histtotXbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l",col="purple",ylim=range(h$counts,h1$counts))
lines(x1,h1$counts,xaxt='n',type="l")
axis(side=1,labels=TRUE)
mtext("Count",side=2,line=1.2)
mtext(bquote(cov(X[i],X[j]) %~~% .(round(cov5,3))),side=3,line=1.6,cex=textsz)
mtext(bquote(sd(Sigma[i] * X[i]) == .(round(sqrt(vartotX5),3))),side=3,line=.8,cex=textsz)
mtext(bquote(sk(Sigma[i] * X[i]) == .(round(sktotX5,3))),side=3,line=0,cex=textsz)
mtext(expression(Sigma[i] * X[i]),side=1,line=1.2)

#cdf panel for exceeding large thresholds (top panel of the cdf panels)
par(fig=c((2*ymargwd+2*panwd+gap+smallpan+ymargwd)/totwd,
          (2*ymargwd+3*panwd+gap+smallpan+ymargwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+cdfpanht+gap)/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+2*cdfpanht+gap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
inds<-which(empcdf1$x>topthreshrg[1] & empcdf1$x<topthreshrg[2])
x<-empcdf1$x[inds]
y<-1-empcdf1$y[inds]
inds<-which(empcdf5$x>topthreshrg[1] & empcdf5$x<topthreshrg[2])
xn<-empcdf5$x[inds]
yn<-1-empcdf5$y[inds]
plot(x,y,type="l",lty="dashed",ylim=c(0,max(y,yn)))
lines(xr_l,yr_l,type="l",lty="dashed",col="red")
lines(xn,yn,type='l',col='purple')
mtext("P>th.",side=2,line=1.2)

#cdf panel for falling under small thresholds (bottom panel of the cdf panels)
par(fig=c((2*ymargwd+2*panwd+gap+smallpan+ymargwd)/totwd,
          (2*ymargwd+3*panwd+gap+smallpan+ymargwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+cdfpanht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
inds<-which(empcdf1$x>botthreshrg[1] & empcdf1$x<botthreshrg[2])
x<-empcdf1$x[inds]
y<-empcdf1$y[inds]
inds<-which(empcdf5$x>botthreshrg[1] & empcdf5$x<botthreshrg[2])
xn<-empcdf5$x[inds]
yn<-empcdf5$y[inds]
plot(x,y,type="l",lty="dashed",ylim=c(0,max(y,yn)))
lines(xr_s,yr_s,type="l",lty="dashed",col="red")
lines(xn,yn,type='l',col='purple')
mtext("P<th.",side=2,line=1.2)
mtext("Threshold (th.)",side=1,line=1.2)

dev.off()

