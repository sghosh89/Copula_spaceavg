source("ncsurrog.R")
library(e1071)
#--------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Input : ts_matrix : a multivariate dataset (a matrix with dimension = timeseries by species)
#         splist : a vector of selective locations 
#         numsurrog : number of replica you want to create normal cop data
#         ploton : logical
#         corpres : either "spearman" or "kendall"
#         resloc : argument for ncsurrog function call (see ncsurrog.R documentation)

# Output : a list of two : 1. a vector which measures skewness of each surrogate normal cop
#                          2. real skewness of given data
#                         and an optional histogram plot

skewness_testing<-function(ts_matrix,splist,numsurrog,ploton,plotcheckon,corpres,resloc){  
  
  ts_tot<-NA*numeric(length=dim(ts_matrix)[1])
  
  for(counter in c(1:dim(ts_matrix)[1])){
    temp<-sum(ts_matrix[counter,splist],na.rm=T)
    ts_tot[counter]<-temp
  }
  
  ts_matrix_w_tot<-cbind(ts_matrix,ts_tot)
  
  skw<-apply(FUN=skewness,X=ts_matrix_w_tot,MARGIN=2,type=2)
  realskw<-tail(skw,1)
  
  # start surrogating
  #numsurrog<-5
  ts_surrog_array<-ncsurrog(m=ts_matrix[,splist],corpres=corpres,numsurrog,plotcheckon,resloc=resloc)
  
  arr_dim1<-dim(ts_surrog_array)[1]
  ts_surrog_tot<-matrix(NA,nrow=arr_dim1,ncol=numsurrog)
  
  for(carray in c(1:numsurrog)){
    for(jj in c(1:arr_dim1)){
      temp<-sum(ts_surrog_array[,,carray][jj,],na.rm=T)
      ts_surrog_tot[jj,carray]<-temp
    }
  }
  
  surrogskw<-apply(FUN=skewness,X=ts_surrog_tot,MARGIN=2,type=2)
  p_left<-sum(surrogskw < realskw)/numsurrog
  p_right<-sum(surrogskw > realskw)/numsurrog
  
  if(ploton==T){
    #cat("--------hi------","\n")
    hist(surrogskw,breaks=100,main="",xlab="surrogate_skewness",ylab="frequency",cex.lab=2,cex.axis=2,col="grey",border=F)
    abline(v=realskw,col="black")
    mtext(paste0("p_left= ",round(p_left,3),sep=""),
          side = 3, line=-2, adj=0, col="blue")
    #points(x=realskw,y=0,col="red",pch=15,cex=2)
  }
  
  return(list(surrogskw=surrogskw,
              realskw=unname(realskw),
              p_left=p_left,
              p_right=p_right))
  
}  
