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
 
  #x<-rbinom(10000,nLU,0.5)
  #if(binom_sig=="LT"){
  #  p<-sum(x>nL)/10000 # note the use of ">" not ">="
  #}else if(binom_sig=="UT"){
  #  p<-sum(x>nU)/10000
  #}else{
  #  stop("Error in summary_taildep_sig: binom_sig must be 'LT' or 'UT'")  
  #}
  
  #if(binom_sig=="LT"){
  #  p<-pbinom(nL,size=nLU,prob=0.5,lower.tail=F) #look up P(X >= nL) when X has the Bin(nLU, 0.5) distribution. 
  #}else if(binom_sig=="UT"){
  #  p<-pbinom(nU,size=nLU,prob=0.5,lower.tail=F) #look up P(X >= nU) when X has the Bin(nLU, 0.5) distribution.
  #}else{
  #  stop("Error in summary_taildep_sig: binom_sig must be 'LT' or 'UT'")  
  #}
  
  bt<-binom.test(x=nL,n=nLU,p=0.5,alternative=c("t"), conf.level=0.95) 
  p<-bt$p.value
  
  res<-list(overall_sig=overall_sig,pval_binom=p)
  
  return(res)
}

#================================================================

#-------------- using rbinom function ----------------

# test code for hays
#set.seed(seed=101)
#ylist<-readRDS("./Results/hays_results/corstat_hays_spaceavg_results/CorlmCoru_hays_spaceavg_nbin_2.RDS")
#res<-binomial_sigtest(ylist=ylist,binom_sig="LT") 
#res$pval_binom #0.0963

# test code for knz
#set.seed(seed=101)
#ylist<-readRDS("./Results/knz_results/corstat_knz_spaceavg_results/CorlmCoru_knz_spaceavg_nbin_2.RDS")
#res<-binomial_sigtest(ylist=ylist,binom_sig="UT") 
#res$pval_binom #0.008

#------------- using pbinom function ----------------

#ylist<-readRDS("./Results/hays_results/corstat_hays_spaceavg_results/CorlmCoru_hays_spaceavg_nbin_2.RDS")
#res<-binomial_sigtest(ylist=ylist,binom_sig="LT")
#res$pval_binom #0.0962

#ylist<-readRDS("./Results/knz_results/corstat_knz_spaceavg_results/CorlmCoru_knz_spaceavg_nbin_2.RDS")
#res<-binomial_sigtest(ylist=ylist,binom_sig="UT")
#res$pval_binom #0.0073

#--------------- can also check similar conclusion with this binom.test() -----------------------------------------
#binom.test(x=34,n=59,p=0.5,alternative=c("g"), conf.level=0.95) #nL=34, nLU=59 for hays data
                                                # gives p-value = 0.1488

#binom.test(x=65,n=106,p=0.5,alternative=c("g"), conf.level=0.95) #nU=65, nLU=106 for knz data
                                                # gives p-value = 0.01251

# if we use ">=" sign in rbinom call, it would give the same results as of binom.test() p-value

#------ so, what to use? ----------------------
# It is reasonable to use a two.tailed binom.test because 
# we should not use a priori reason to get stronger LT or stronger UT
# dep. cells in majority.

#--------------------- summary -----------------------------
# hays data shows dominance of LT dep. cells in the community
# could be by chance (as p>0.05), whereas for knz data UT dep. 
# cells in the community is significantly dominant (as p<0.05).

#------------- web-help-----------------------------------
# https://www.stat.umn.edu/geyer/old/5101/rlook.html
# http://www.endmemo.com/program/R/binomial.php




