library(mvtnorm)
source("alignranks.R")

#Computes a map from the space of 2D normal-copula parameters
#to the space of Pearson correlation values, based on d
#
#Args
#d          A L by 2 matrix of real numbers
#numpts     The normal-copula parameters 
#             seq(from=-1,to=1,length.out=numpts) are used
#numsims    This many evaluations of the map per value of numpts
#plotnm     NA (default) for no plot, otherwise the file name (without 
#             extension) for a pdf plot to be saved in the working
#             directory
#
#Output is a list consisting of these elements:
#numeric_result     A numsims by numpts matrix, the results (numsims
#                     of them for each domain value) of applying the map
#                     to the points seq(from=-1,to=1,length.out=numpts)
#fit_parameters     A data frame with first column the values 
#                     seq(from=-1,to=1,length.out=numpts) and second column
#                     the mean of the function values over that input. These
#                     specify a piecewise-linear function.
#
#Details
#No error checking of inputs at present
#
getmap<-function(d,numpts,numsims,plotnm=NA)
{
  #***prep
  L<-dim(d)[1]
  d<-apply(FUN=sort,X=d,MARGIN=2)
  
  #***numerically evaluate the map
  
  #make a receptacle for the numeric evaluations
  numcres<-matrix(NA,numsims,numpts)
  
  #now do the evaluations
  mycorfun<-function(x)
  {
    res<-cor(x)
    return(res[1,2])
  }
  parmvals<-seq(from=-1,to=1,length.out=numpts)
  for (counter in 1:length(parmvals))
  {
    sigma<-matrix(c(1,parmvals[counter],parmvals[counter],1),2,2)
    sims<-rmvnorm(numsims*L,sigma=sigma)
    sims<-array(sims,c(L,numsims,2))
    sims<-aperm(sims,c(1,3,2))
    mapres<-alignranks(d,sims)
    numcres[,counter]<-apply(FUN=mycorfun,X=mapres,MARGIN=3)
  }
  
  #***carry out fitting to establish the mathematical version of the map,
  #and its inverse
  mnnumcres<-apply(FUN=mean,X=numcres,MARGIN=2)
  quantres<-apply(FUN=quantile,X=numcres,MARGIN=2,probs=c(.025,.25,.5,.75,.975))
  
  x<-parmvals
  x<-(x+1)/2 #transform so x goes from 0 to 1
  y<-numcres
  ytop<-mnnumcres[length(mnnumcres)]
  ybot<-mnnumcres[1]
  y<-(y-ybot)/(ytop-ybot) #transform so the mean values of y go from 0 to 1
  x<-x[2:(length(x)-1)] #throw out the endpoint because you only want to fit on what's left
  y<-y[,2:(dim(y)[2]-1)]
  x<-matrix(rep(x,each=numsims),numsims,numpts-2)
  myfitfun<-function(expon,x,y)
  {
    return(sum((y-x^expon)^2))
  }
  optres<-optimize(myfitfun,c(.1,10),x=x,y=y)
  expon<-optres$minimum
  prefactor<-0.5^expon*(ytop-ybot)
  intercept<-ybot
  fitparms<-c(expon=expon,prefactor=prefactor,intercept=intercept)

  #***make a plot, if required
  if (!is.na(plotnm))
  {
    pdf(file=paste0(plotnm,".pdf"))
    plot(parmvals,mnnumcres,type='l',ylim=range(mnnumcres,quantres),
         xlab="Normal-copula parameter",ylab="Output Pearson")
    lines(parmvals,quantres[1,],type="l",lty="dotted")
    lines(parmvals,quantres[2,],type="l",lty="dashed")
    lines(parmvals,quantres[3,],type="l",lty="dashed")
    lines(parmvals,quantres[4,],type="l",lty="dashed")
    lines(parmvals,quantres[5,],type="l",lty="dotted")
    newparmvals<-seq(from=parmvals[1],to=parmvals[length(parmvals)],length.out=1000)
    lines(newparmvals,prefactor*(newparmvals+1)^expon+intercept,type='l',col="red")
    dev.off()
  }
  
  #***construct the list that is to be returned and return it
  return(list(numeric_result=numcres,fit_parameters=fitparms))
}