#This is the objective function optimized to search for Pearson preserving
#surrogates for a given dataset. Written to be passed to constrOptim.
#
#Args
#pij          The parameters to optimize
#sumvij       The sum of the (strict) upper triangle of the correlation matrix 
#               of the data, d 
#mapx, mapy   The parameters of the piecewise linear maps obtained by getmap
#saveloc      The location where the first pos def mat is to be saved
#
#Output
#A single number which is a minimum eigenvalue of a symmetric matrix. This should
#be maximized to look for a positive definite matrix.

PPSurrogObj<-function(pij,sumvij,mapx,mapy,saveloc)
{
  #extend pij and put it into the upper traingle
  pij<-c(pij,sumvij-sum(pij)) #The constraints passed to constOptim will be chosen so
  #that all entries of pij, including this newly added one, are not only between -1 
  #and 1, but also are within the range that can be mapped back with the inverse of 
  #the relevant map described in mapx and mapy. 
  pijmat<-matrix(NA,dim(mapy)[1],dim(mapy)[2])
  pijmat[upper.tri(pijmat)]<-pij
  
  #map it back from Pearson space to N-copula parameter space
  nijmat<-matrix(NA,dim(mapy)[1],dim(mapy)[2])
  for (iind in 1:(dim(mapy)[1]-1))
  {
    for (jind in (iind+1):(dim(mapy)[2]))
    {
      nijmat[iind,jind]<-pwlin(pijmat[iind,jind],mapy[iind,jind,],mapx)
    }
  }
  
  #form a symmetric matrix with 1s on the diagonal
  nijmat[is.na(nijmat)]<-0
  nijmat<-nijmat+t(nijmat)
  diag(nijmat)<-1
  
  #calculate the minimum eigenvalue and return, but save and stop if you find
  #a positive value because actually you don't want to maximize you just want the
  #first positive value you can find
  res<-eigen(nijmat,symmetric=TRUE,only.values=TRUE)$values
  if (min(res)>0)
  {
    saveRDS(pijmat,file=paste0(saveloc,"FirstSuccess_pijmat.RDS"))
    saveRDS(nijmat,file=paste0(saveloc,"FirstSuccess_nijmat.RDS"))
    stop("Success!")
  }
  return(min(res))
}