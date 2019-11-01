library(copula)
library(mvtnorm)
source("CopulaFunctions_flexible.R")

#' A function for looking at statistical significance of a measure of 
#' asymmetry of tail dependence, against a normal-copula null
#' hypothesis. The measure is corl-coru for positively associated 
#' variables and the absolute value of this quantity, computed after
#' flipping one of the variables, for negatively associated variables.
#' 
#' @param ub The measure of asymmetry of tail dependence is cor_{0,ub}-
#' cor_{1-ub,1} for positively associated variables, and the absolute 
#' value of this quantity, computed after flipping one of the variables, 
#' for negatively associated variables
#' @param numpts This many points from the null generated for each run
#' @param spcors The distribution of the statistic for these values of
#' the Spearman correlation is computed
#' @param numsims For each value of spcors, this many runs done
#' @param CI  a vector of lower and upper bounds for confidence interval used in significance testing, 95%CI as default
#' @param sigplton logical to have optional plot
#' @param resloc location folder path to save the plot
#' 

# return a matrix with two rows for lower and upper quantiles

tailsignif<-function(ub,numpts,spcors,numsims,CI=c(0.025,0.975),resloc){
  #a small amount of error checking
  if (!is.numeric(ub) || length(ub)!=1 || ub>.5 || ub<=0){
    stop("Error in tailsignif: bad value for ub")  
  }
  if (any(spcors<=-1) || any(spcors>=1)){
    stop("Error in tailsignif: bad value for spcors")
  }
  
  #get the numeric results
  allres<-matrix(NA,numsims,length(spcors))
  dumnorm<-normalCopula(.1,2)
  for (ct_spcors in 1:length(spcors)){
    #find the parameter for which the normal copula has this Spearman
    thisparm<-iRho(dumnorm,spcors[ct_spcors])
    thissig<-matrix(c(1,thisparm,thisparm,1),2,2)
    
    for (ct_sims in 1:numsims){
      #generate numpts data from the normal copula
      d<-rmvnorm(numpts,mean=c(0,0),sig=thissig)
      d<-pobs(d) #note are doing it this way instead of using normalCopula and rCopula
      #because real data will also have had pobs applied to it
      
      #compute the statistic and store
      if (spcors[ct_spcors]>=0){
        allres[ct_sims,ct_spcors]<-Corbds(d[,1],d[,2],0,ub)-Corbds(d[,1],d[,2],1-ub,1)
      }
      
      if (spcors[ct_spcors]<0){ # for negatively correlated cells only absolute value of tail diff. matters
        allres[ct_sims,ct_spcors]<-abs(Corbds(d[,1],-d[,2]+1,0,ub)-Corbds(d[,1],-d[,2]+1,1-ub,1)) 
      }
      
    }
  }
  
  qtl<-apply(FUN=quantile,X=allres,MARGIN=2,prob=CI)
  
  return(qtl)
}
