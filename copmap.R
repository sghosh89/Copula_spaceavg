library(copula)
source("alignranks.R")

#Creates the basis of the map that needs to be inverted
#
#Given two vectors of the same length (x and y), and a vector of values of
#the parameter of the bivariate normal copula family (p), for each entry p[i]
#of p, generates data from the bivariate normal copula with that parameter, 
#maps the ranks of x and y onto the data to acheive permutations of x and y 
#that have the normal copula structure with the parameter p[i], and then 
#computes the covariance. Does so numreps time for each value of p. 
#Returns all the resulting values.
#
#Args
#x, y     two time series of the same length
#p        a vector of values in [-1,1], parameters of a bivariate normal copula
#numreps  number of replicates for each value of p
#
#Output
#A matrix of dimensions numreps by length(p) containing the covariance results
copmap<-function(x,y,p,numreps)
{
  res<-matrix(NA,numreps,length(p))
  
  lenx<-length(x)
  sx<-sort(x)
  sy<-sort(y)
  
  for (pcount in 1:length(p))
  {
    ncop<-normalCopula(p[pcount],2)
    
    #generate a bunch of data from the normal copula, in the end it will be 
    #lenx by 2 by numreps
    ncopdat<-rCopula(lenx*numreps,ncop)
    ncopdat<-aperm(array(ncopdat,c(lenx,numreps,2)),c(1,3,2))
    
    #map ranks
    surrogs<-alignranks(cbind(sx,sy),ncopdat)
    
    #compute the correlations
    res[,pcount]<-apply(FUN=function(x){res<-cov(x);return(res[1,2])},
                        MARGIN=3,X=surrogs)
  }
  
  return(res)
}
