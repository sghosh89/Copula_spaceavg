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
#                   copula with this parameter
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

numpts<-10000
numscatpts<-500

#Functions which when applied to uniformly distributed data on (0,1) give X and Y.
#To be used to establish the marginals in all 5 examples
#margXfunc<-function(c){qgamma(c,shape=2,scale=2)}
#margYfunc<-function(c){qgamma(c,shape=7.5,scale=1)}
margXfunc<-function(c){qgamma(c,shape=9,scale=.5)}
margYfunc<-function(c){qgamma(c,shape=7.5,scale=1)}
#margXfunc<-function(c){qnorm(c)}
#margYfunc<-function(c){qnorm(c)}

#get data for example 1 - normal copula case
sig<-matrix(c(1,.6,.6,1),2,2)
cop1<-normalCopula(param=P2p(sig),dim=2,dispstr="un")
dcop1<-rCopula(numpts,cop1)
d1<-dcop1
d1[,1]<-margXfunc(d1[,1])
d1[,2]<-margYfunc(d1[,2])
cov1<-cov(d1)[1,2]
XpY1<-d1[,1]+d1[,2]
varXpY1<-var(XpY1)
sk1<-myskns(XpY1)

#get data for example 2 - comonotonic
dcop2<-matrix(rep(runif(numpts),times=2),numpts,2)
d2<-dcop2
d2[,1]<-margXfunc(d2[,1])
d2[,2]<-margYfunc(d2[,2])
cov2<-cov(d2)[1,2]
XpY2<-d2[,1]+d2[,2]
varXpY2<-var(XpY2)
sk2<-myskns(XpY2)

#get data for example 3 - upper-tail comonotonic
d3<-matrix(NA,numpts,2)
beta1<-0
beta2<-0.95
objfun<-function(rho)
{
  cdat<-taildepcopdat(rho,beta1,beta2,numpts)
  dat<-cdat
  dat[,1]<-margXfunc(dat[,1])
  dat[,2]<-margYfunc(dat[,2])
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
dcop3<-taildepcopdat(parmtouse,beta1,beta2,numpts)
d3<-dcop3
d3[,1]<-margXfunc(d3[,1])
d3[,2]<-margYfunc(d3[,2])
cov3<-cov(d3)[1,2]
XpY3<-d3[,1]+d3[,2]
varXpY3<-var(XpY3)
sk3<-myskns(XpY3)

#get data for example 4 - lower-tail comonotonic
d4<-matrix(NA,numpts,2)
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
dcop4<-taildepcopdat(parmtouse,beta1,beta2,numpts)
d4<-dcop4
d4[,1]<-margXfunc(d4[,1])
d4[,2]<-margYfunc(d4[,2])
cov4<-cov(d4)[1,2]
XpY4<-d4[,1]+d4[,2]
varXpY4<-var(XpY4)
sk4<-myskns(XpY4)

#get data for example 5
d5<-matrix(NA,numpts,2)
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
dcop5<-taildepcopdat(parmtouse,beta1,beta2,numpts)
d5<-dcop5
d5[,1]<-margXfunc(d5[,1])
d5[,2]<-margYfunc(d5[,2])
cov5<-cov(d5)[1,2]
XpY5<-d5[,1]+d5[,2]
varXpY5<-var(XpY5)
sk5<-myskns(XpY5)

#other prep for consistency among plots
scatylim<-range(d1[,2],d2[,2],d3[,2],d4[,2],d5[,2],na.rm=TRUE)
scatxlim<-range(d1[,1],d2[,1],d3[,1],d4[,1],d5[,1],na.rm=TRUE)
histXbreaks<-seq(from=floor(scatxlim[1]),to=ceiling(scatxlim[2]),by=0.5)
histYbreaks<-seq(from=floor(scatylim[1]),to=ceiling(scatylim[2]),by=0.5)
histXpYxlim<-range(XpY1,XpY2,XpY3,XpY4,XpY5,na.rm=TRUE)
histXpYbreaks<-seq(from=floor(histXpYxlim[1]),to=ceiling(histXpYxlim[2]),by=1)

#***plot dimensions, units inches
xmarght<-.5
ymargwd<-.5
totwd<-3.5
gap<-0.15
spf<-0.4
panwd<-(totwd-2*gap-2*ymargwd)/(2+spf)
smallpan<-panwd*spf
panht<-panwd
totht<-xmarght+5*(panht+2*gap+smallpan)
pdf(file="TheoryFig.pdf",width=totwd,height=totht)

#***Example 1, top row of panels, normal copula

#scatterplot
panrownum<-5
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
plot(d1[,1],d1[,2],xaxt='n',pch=20,cex=.5,col=rgb(red=190/256,green=190/256,blue=190/256,alpha=.3),
     xlim=scatxlim,ylim=scatylim)
axis(side=1,labels=FALSE)
mtext("Y",side=2,line=1.2)

#X marginal of scatterplot
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap)/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht+gap+smallpan)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d1[,1],breaks=histXbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Ct.",side=2,line=1.2)

#Y marginal of scatterplot
par(fig=c((ymargwd+panwd+gap)/totwd,
          (ymargwd+panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d1[,2],breaks=histYbreaks,plot=FALSE)
y<-h$breaks
y<-y[2:length(y)]-diff(y)[1]/2
plot(h$counts,y,yaxt='n',type="l")
axis(side=2,labels=FALSE)
#mtext("Ct.",side=1,line=1.2)

#Histogram of X+Y
par(fig=c((2*ymargwd+panwd+gap+smallpan)/totwd,
          (2*ymargwd+2*panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h1<-hist(XpY1,breaks=histXpYbreaks,plot=FALSE)
x1<-h1$breaks
x1<-x1[2:length(x1)]-diff(x1)[1]/2
plot(x1,h1$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Count",side=2,line=1.2)
mtext(paste0("cov(X,Y)=",round(cov1,3)),side=3,line=1.6,cex=.75) #adjust the horizontal position,
mtext(paste0("var(X+Y)=",round(varXpY1,3)),side=3,line=.8,cex=.75) #adjust the horizontal position
mtext(paste0("sk(X+Y)=",round(sk1,3)),side=3,line=0,cex=.75) #adjust the horizontal position

#***Example 2, second row of panels, comonotonic case

#scatterplot
panrownum<-4
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d2[,1],d2[,2],xaxt='n',pch=20,cex=.5,col=rgb(red=255/256,green=0/256,blue=0/256,alpha=.3),
     xlim=scatxlim,ylim=scatylim)
axis(side=1,labels=FALSE)
mtext("Y",side=2,line=1.2)

#X marginal of scatterplot
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

#Y marginal of scatterplot
par(fig=c((ymargwd+panwd+gap)/totwd,
          (ymargwd+panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d2[,2],breaks=histYbreaks,plot=FALSE)
y<-h$breaks
y<-y[2:length(y)]-diff(y)[1]/2
plot(h$counts,y,yaxt='n',type="l")
axis(side=2,labels=FALSE)
#mtext("Ct.",side=1,line=1.2)

#Histogram of X+Y
par(fig=c((2*ymargwd+panwd+gap+smallpan)/totwd,
          (2*ymargwd+2*panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(XpY2,breaks=histXpYbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l",col="red",ylim=range(h$counts,h1$counts))
lines(x1,h1$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Count",side=2,line=1.2)
mtext(paste0("cov(X,Y)=",round(cov2,3)),side=3,line=1.6,cex=.75) #adjust the horizontal position,
mtext(paste0("var(X+Y)=",round(varXpY2,3)),side=3,line=.8,cex=.75) #adjust the horizontal position
mtext(paste0("sk(X+Y)=",round(sk2,3)),side=3,line=0,cex=.75) #adjust the horizontal position

#***Example 3, third row of panels, comonotonicity in the right tails only

#scatterplot
panrownum<-3
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d3[,1],d3[,2],xaxt='n',pch=20,cex=.5,col=rgb(red=0/256,green=0/256,blue=255/256,alpha=.3),
     xlim=scatxlim,ylim=scatylim)
axis(side=1,labels=FALSE)
mtext("Y",side=2,line=1.2)

#X marginal of scatterplot
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

#Y marginal of scatterplot
par(fig=c((ymargwd+panwd+gap)/totwd,
          (ymargwd+panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d3[,2],breaks=histYbreaks,plot=FALSE)
y<-h$breaks
y<-y[2:length(y)]-diff(y)[1]/2
plot(h$counts,y,yaxt='n',type="l")
axis(side=2,labels=FALSE)
#mtext("Ct.",side=1,line=1.2)

#Histogram of X+Y
par(fig=c((2*ymargwd+panwd+gap+smallpan)/totwd,
          (2*ymargwd+2*panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(XpY3,breaks=histXpYbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l",col="blue",ylim=range(h$counts,h1$counts))
lines(x1,h1$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Count",side=2,line=1.2)
mtext(paste0("cov(X,Y)=",round(cov3,3)),side=3,line=1.6,cex=.75) #adjust the horizontal position,
mtext(paste0("var(X+Y)=",round(varXpY3,3)),side=3,line=.8,cex=.75) #adjust the horizontal position
mtext(paste0("sk(X+Y)=",round(sk3,3)),side=3,line=0,cex=.75) #adjust the horizontal position

#***Example 4, fourth row of panels, comonotonicity in the left tails only

#scatterplot
panrownum<-2
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d4[,1],d4[,2],xaxt='n',pch=20,cex=.5,col=rgb(red=0/256,green=255/256,blue=0/256,alpha=.3),
     xlim=scatxlim,ylim=scatylim)
axis(side=1,labels=FALSE)
mtext("Y",side=2,line=1.2)

#X marginal of scatterplot
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

#Y marginal of scatterplot
par(fig=c((ymargwd+panwd+gap)/totwd,
          (ymargwd+panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d4[,2],breaks=histYbreaks,plot=FALSE)
y<-h$breaks
y<-y[2:length(y)]-diff(y)[1]/2
plot(h$counts,y,yaxt='n',type="l")
axis(side=2,labels=FALSE)
#mtext("Ct.",side=1,line=1.2)

#Histogram of X+Y
par(fig=c((2*ymargwd+panwd+gap+smallpan)/totwd,
          (2*ymargwd+2*panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(XpY4,breaks=histXpYbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l",col="green",ylim=range(h$counts,h1$counts))
lines(x1,h1$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Count",side=2,line=1.2)
mtext(paste0("cov(X,Y)=",round(cov4,3)),side=3,line=1.6,cex=.75) #adjust the horizontal position, fill in number
mtext(paste0("var(X+Y)=",round(varXpY4,3)),side=3,line=.8,cex=.75) #adjust the horizontal position, fill in number
mtext(paste0("sk(X+Y)=",round(sk4,3)),side=3,line=0,cex=.75) #adjust the horizontal position, fill in number

#***Example 5, fifth row of panels, comonotonicity in the right and left tails only

#scatterplot
panrownum<-1
par(fig=c((ymargwd)/totwd,
          (ymargwd+panwd)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d5[,1],d5[,2],xaxt='n',pch=20,cex=.5,col=rgb(red=160/256,green=32/256,blue=240/256,alpha=.3),
     xlim=scatxlim,ylim=scatylim)
axis(side=1,labels=FALSE)
mtext("Y",side=2,line=1.2)
mtext("X",side=1,line=1.2)

#X marginal of scatterplot
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

#Y marginal of scatterplot
par(fig=c((ymargwd+panwd+gap)/totwd,
          (ymargwd+panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(d5[,2],breaks=histYbreaks,plot=FALSE)
y<-h$breaks
y<-y[2:length(y)]-diff(y)[1]/2
plot(h$counts,y,yaxt='n',type="l")
axis(side=2,labels=FALSE)
mtext("Ct.",side=1,line=1.2)

#Histogram of X+Y
par(fig=c((2*ymargwd+panwd+gap+smallpan)/totwd,
          (2*ymargwd+2*panwd+gap+smallpan)/totwd,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan))/totht,
          (xmarght+(panrownum-1)*(panht+2*gap+smallpan)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
h<-hist(XpY5,breaks=histXpYbreaks,plot=FALSE)
x<-h$breaks
x<-x[2:length(x)]-diff(x)[1]/2
plot(x,h$counts,xaxt='n',type="l",col="purple",ylim=range(h$counts,h1$counts))
lines(x1,h1$counts,xaxt='n',type="l")
axis(side=1,labels=FALSE)
mtext("Count",side=2,line=1.2)
mtext(paste0("cov(X,Y)=",round(cov5,3)),side=3,line=1.6,cex=.75) #adjust the horizontal position, fill in number
mtext(paste0("var(X+Y)=",round(varXpY5,3)),side=3,line=.8,cex=.75) #adjust the horizontal position, fill in number
mtext(paste0("sk(X+Y)=",round(sk5,3)),side=3,line=0,cex=.75) #adjust the horizontal position, fill in number
mtext("X+Y",side=1,line=1.2)

dev.off()