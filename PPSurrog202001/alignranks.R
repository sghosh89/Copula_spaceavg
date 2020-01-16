#Rearranges data so the ranks match with simulations.
#
#This is part of various copula surrogate algorithms.
#For each column j of each simulation sims[,,k], the 
#elements of dat[,j] are permuted so the ranks of the
#result are aligned with those of sims[,j,k].
#
#Args
#dat      A T by d matrix, each column sorted
#sims     A T by d by numsims array
#
#Output
#A T by d by numsims array
#
#Notes
#The algorithm assumes the column sims[,j,k] has no 
#ties, for any j, k
#
alignranks<-function(dat,sims)
{
  simsrk<-apply(FUN=rank,MARGIN=c(2,3),X=sims)
  res<-array(NA,dim(sims))
  for (counter in 1:dim(dat)[2])
  {
    res[,counter,]<-dat[simsrk[,counter,],counter]
  }
  
  return(res)
}