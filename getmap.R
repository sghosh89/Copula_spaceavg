library(mvtnorm)
source("alignranks.R")

#Computes a map from the space of 2D normal-copula parameters
#to the space of Pearson correlation values, based on d
#
#Args
#d          A L by 2 matrix of real numbers
#numpts     The normal-copula parameters 
#             seq(from=-1,to=1,length.out=numpts) are used
#numsims    This many evaluations of the map per value of normal cop params
#plotnm     NA (default) for no plot, otherwise the file name (without 
#             extension) for a pdf plot to be saved 
#
#Output is a list consisting of these elements:
#numeric_result     A numsims by numpts matrix, the results (numsims
#                     of them for each domain value) of applying the map
#                     to the points seq(from=-1,to=1,length.out=numpts)
#fit_parameters     A data frame with first column the values 
#                     seq(from=-1,to=1,length.out=numpts) and second column
#                     the mean of the function values over that input. These
#                     specify a piecewise-linear function which approximates
#                     the map.
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
  
  #***now get the piecewise linear map
  mnnumcres<-apply(FUN=mean,X=numcres,MARGIN=2)
  quantres<-apply(FUN=quantile,X=numcres,MARGIN=2,probs=c(.025,.25,.5,.75,.975))
  fitparms<-data.frame(x=parmvals,y=mnnumcres)  
  
  #***make a plot, if required
  if (!is.na(plotnm))
  {
    pdf(file=paste0(plotnm,".pdf"))
    plot(parmvals,mnnumcres,type='b',ylim=range(mnnumcres,quantres),
         xlab="Normal-copula parameter",ylab="Output Pearson",col="red")
    lines(parmvals,quantres[1,],type="l",lty="dotted")
    lines(parmvals,quantres[2,],type="l",lty="dashed")
    lines(parmvals,quantres[3,],type="l",lty="dashed")
    lines(parmvals,quantres[4,],type="l",lty="dashed")
    lines(parmvals,quantres[5,],type="l",lty="dotted")
    dev.off()
  }
  
  #***construct the list that is to be returned and return it
  return(list(numeric_result=numcres,fit_parameters=fitparms))
}