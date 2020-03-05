# This code generate some plots just to confiem Pearson preserving surrogates are "good" or not

# This function gives you a dataframe 
# ARGs : 
#           m = a matrix with sp. time series along each column
#           surrogs : (pearson preserving probably) surrogates

source("./get_var_ratio.R")

PPsurrogs_tests<-function(m,surrogs){
  
  #get correlation matrix for real data
  cor_real<-cor(m,method="pearson")
  
  # get variance for aggregated sp. timeseries from real data
  m_aggr_sp<-apply(m,MARGIN = 1,FUN = sum)
  var_aggr_sp_real<-var(m_aggr_sp)
  
  #get covariance matrix for real data
  m_cov<-cov(m,method="pearson")
  
  #get the sum over the covariance matrix for upper triangular part: i>j
  tot_cov_real<-sum(m_cov[lower.tri(m_cov)])
  
  # get variance ratio for real data
  vr_real<-get_var_ratio(m)
  
  #initialize
  cor_surrogs<-array(numeric(),c(nrow(cor_real),ncol(cor_real),0)) # to store pairwise correlation for surrogates
  var_aggr_sp_surrogs<-c() # to store variance of aggregated species timeseries from each surrogate
  tot_cov_surrogs<-c() # to store sum of covariances from each surrogate
  vr_surrogs<-c() # to store variance ratio from each surrogate
  
  for(k in c(1:dim(surrogs)[3])){
    
    ms<-surrogs[,,k] # surrogate matrix
    
    temp<-cor(ms)
    cor_surrogs<-abind::abind(cor_surrogs,temp) # array with correlation for surrogate values
    
    # get variance for aggregated sp. timeseries from each surrogates
    ms_aggr_sp<-apply(ms,MARGIN = 1,FUN = sum)
    get_var_aggr_sp_surrogs<-var(ms_aggr_sp)
    var_aggr_sp_surrogs<-c(var_aggr_sp_surrogs,get_var_aggr_sp_surrogs)
    
    #get covariance matrix for each surrogate
    ms_cov<-cov(ms,method="pearson")
    
    # get the sum over the covariance matrix for upper triangular part: i>j for each surrogates
    get_ms_cov<-sum(ms_cov[lower.tri(ms_cov)])
    tot_cov_surrogs<-c(tot_cov_surrogs,get_ms_cov)
    
    # get variance ratio for surrogate
    gvs<-get_var_ratio(ms)
    vr_surrogs<-c(vr_surrogs,gvs)
  }
  
  return(list(var_aggr_sp_real=var_aggr_sp_real,
              tot_cov_real=tot_cov_real,
              vr_real=vr_real,
              var_aggr_sp_surrogs=var_aggr_sp_surrogs,
              tot_cov_surrogs=tot_cov_surrogs,
              vr_surrogs=vr_surrogs,
              cor_real=cor_real,
              cor_surrogs=cor_surrogs
  ))
  
}

