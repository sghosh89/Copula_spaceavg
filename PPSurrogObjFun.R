#This is the objective function optimized to search for Pearson preserving
#surrogates for a given dataset. Written to be passed to constrOptim.
#
#Args
#pij          The parameters to optimize, in "correlation space"
#cijmatd      Covariance matrix of the data
#mapx, mapy   The parameters of the piecewise linear maps obtained by getmap
#saveloc      The location where the first pos def mat is to be saved
#
#Output
#A single number which is a minimum eigenvalue of a symmetric matrix. This should
#be maximized to look for a positive definite matrix. Note the function throws
#an error whenever a positive value is found, since we don't actually want to
#maximize, we just want to find the first pos def matrix we come across.

PPSurrogObj<-function(pij,cijmatd,mapx,mapy,saveloc)
{
  #make pij into an appropriate matrix
  n<-dim(mapy)[1] #number of species
  pijmat<-matrix(NA,n,n)
  pijmat[upper.tri(pijmat)]<-pij
  #cijmat<-pijmat*sqrt(outer(diag(cijmatd),diag(cijmatd)))

  #The constraints passed to constOptim will be chosen so
  #that all entries of pij are within the range that can be 
  #mapped back with the inverse of the relevant map described 
  #in mapx and mapy. 

  #map it back from Pearson space to N-copula parameter space
  nijmat<-matrix(NA,n,n)
  for (iind in 1:(n-1))
  {
    for (jind in (iind+1):n)
    {
      #print(paste0("iind: ",iind,"; jind: ",jind))
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
  print(min(res))
  return(min(res))
}