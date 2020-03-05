# calculate variance ratio for a matrix where each column represents a time-series for a species

source("mycvsq.R")

get_var_ratio<-function(m){
  
  tot_quantity<-apply(m,MARGIN = 1,FUN = sum)
  var_each_sp<-apply(m,MARGIN = 2,FUN = var)
  #mean_each_sp<-apply(m,MARGIN = 2,FUN = mean)
  
  #cvsq_real<-mycvsq(tot_quantity)
  #cvsq_indep<-sum(var_each_sp)/((sum(mean_each_sp))^2)
  #phi_cvsq<-cvsq_real/cvsq_indep
  
  var_ratio<-var(tot_quantity)/sum(var_each_sp) # same as phi_cvsq
  
  return(var_ratio)
}

#get_var_ratio(m)

