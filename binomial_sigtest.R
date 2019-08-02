# function to get a summary for binomial significance test results 
# Arg:
#     ylist : a list which is the output from NonParamStat_matrixplot function
#     binom_sig : a character; either "LT" or "UT" to be used in binomial testing
# Output:
#     A list of two numbers : one is the overall significance of tail asymmetry
#                             other is the pvalue from binomial test results

binomial_sigtest<-function(ylist,binom_sig){
  
  s<-ylist$res_sig$is_sigmat
  csig<-sum(s==1,na.rm=T) #count on sig. cells (both +ve and -ve correlated)
  cinsig<-sum(s==0,na.rm=T) #count on insig. cells (both +ve and -ve correlated)
  totalcells<-csig+cinsig
  overall_sig<-(csig/totalcells) # overall significance for asymmetry in tail-dep.
  
  y<-ylist$res_sig$mat_tab
  y1<-subset(y,realstat>0 & sprvals>0) # for +ve correlated lower tail dep. cells
  
  sL<-s[cbind(y1$row,y1$col)]
  nL<-sum(is.finite(sL)) # number of positively correlated lower tail dep. cells
  
  y2<-subset(y,realstat<0 & sprvals>0) # for +ve correlated upper tail dep. cells
  sU<-s[cbind(y2$row,y2$col)]
  nU<-sum(is.finite(sU)) # number of positively correlated upper tail dep. cells
  nLU<-nL+nU
  x<-rbinom(10000,nLU,0.5)
  
  if(binom_sig=="LT"){
    p<-sum(x>=nL)/10000
  }else if(binom_sig=="UT"){
    p<-sum(x>=nU)/10000
  }else{
    stop("Error in summary_taildep_sig: binom_sig must be 'LT' or 'UT'")  
  }
  
  res<-list(overall_sig=overall_sig,pval_binom=p)
  
  return(res)
}

#================================================================

# test code for hays
set.seed(seed=101)
ylist<-readRDS("./Results/hays_results/corstat_hays_spaceavg_results/CorlmCoru_hays_spaceavg_nbin_2.RDS")
res<-binomial_sigtest(ylist=ylist,binom_sig="LT") 
res

# test code for knz
set.seed(seed=101)
ylist<-readRDS("./Results/knz_results/corstat_knz_spaceavg_results/CorlmCoru_knz_spaceavg_nbin_2.RDS")
res<-binomial_sigtest(ylist=ylist,binom_sig="UT") 
res



