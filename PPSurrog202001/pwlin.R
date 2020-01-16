#Implements a piecewise linear function
#
#Args
#x          A vector of values to evaluate the function at
#xbp        A vector of strictly increasing values, the x-axis breakpoint values of the function
#ybp        A vector of the same length, the values of the function at the breakpoints
#
#Output
#A vector of values, NA for indices for which x was outside the range min(xbp) to 
#max(xbp)

pwlin<-function(x,xbp,ybp)
{
  #some small amount of error checking
  if (any(diff(xbp)<=0))
  {
    stop("Error in pwlin: xbp values must be strictly increasing")
  }
  
  #receptacle for results
  res<-NA*numeric(length(x))
  
  inds<-which(x>=xbp[1] & x<=xbp[length(xbp)])
  for (counter in 1:length(inds))
  {
    thisx<-x[inds[counter]]
    ilo<-max(which(xbp<=thisx))
    ihi<-min(which(xbp>=thisx))
    if (ilo==ihi)
    {
      res[inds[counter]]<-ybp[ilo]
    }else
    {
      res[inds[counter]]<-((thisx-xbp[ilo])/(xbp[ihi]-xbp[ilo]))*(ybp[ihi]-ybp[ilo])+ybp[ilo]
    }
  }
  
  return(res)
}

