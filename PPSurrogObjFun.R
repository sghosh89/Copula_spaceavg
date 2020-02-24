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
  #get plast
  n<-dim(mapy)[1] #number of species
  pijmat<-matrix(NA,n,n)
  pijmat[upper.tri(pijmat)]<-c(pij,NA)
  cijmat<-pijmat*sqrt(outer(diag(cijmatd),diag(cijmatd)))
  cijmatd_ut<-cijmatd
  cijmatd_ut[!upper.tri(cijmatd_ut)]<-NA
  plast<-(sum(cijmatd_ut,na.rm=TRUE)-sum(cijmat,na.rm=TRUE))/(sqrt(cijmatd[n-1,n-1]*cijmatd[n,n]))
  
  #insert plast into its place in pijmat
  pijmat[n-1,n]<-plast
  #The constraints passed to constOptim will be chosen so
  #that all entries of pij, including this newly added one, are within the range that 
  #can be mapped back with the inverse of the relevant map described in mapx and mapy. 

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
  return(min(res))
}