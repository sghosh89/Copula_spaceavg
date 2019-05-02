source("mycvsq.R")
source("SkewnessAnd3CentMom.R")

# This function gives you a dataframe 
# ARGs : 
#           m = a matrix with sp. time series along each column
#           surrogs : (pearson preserving probably) surrogates
#           surrogs_given : logical

make_tab_stability<-function(m,surrogs,surrogs_given){
  
  #m<-m[,id_sp]
  tot_quantity<-apply(m,MARGIN = 1,FUN = sum)
  var_each_sp<-apply(m,MARGIN = 2,FUN = var)
  mean_each_sp<-apply(m,MARGIN = 2,FUN = mean)
  cm3_each_sp<-apply(m,MARGIN = 2,FUN = my3cm)
  
  df_stability<-as.data.frame(matrix(NA,nrow=1,ncol=6))
 
  colnames(df_stability)<-c("cvsq_real", # square of CV for the real community data
                            "cvsq_indep", # square of CV of the community if all the sp. behave independently
                            "phi_cvsq",  # this is the ratio of cvsq_real/cvsq_indep and compared to 1
                            "skw_real",  # skewness for the real community data
                            "skw_indep", # skewness of the community if all the sp. behave independently
                            "phi_skw" # this is the ratio of skw_real/skw_indep and compared to ??
  )
  
  df_stability$cvsq_real<-mycvsq(tot_quantity)
  df_stability$cvsq_indep<-sum(var_each_sp)/((sum(mean_each_sp))^2)
  df_stability$phi_cvsq<-df_stability$cvsq_real/df_stability$cvsq_indep
  
  df_stability$skw_real<-myskns(tot_quantity)
  df_stability$skw_indep<-sum(cm3_each_sp)/((sum(var_each_sp))^(3/2))
  df_stability$phi_skw<-df_stability$skw_real/df_stability$skw_indep
  
  if(surrogs_given==T){
    # Now start computation for surrogates
    df_stability_surrogs<-as.data.frame(matrix(NA,nrow=1,ncol=6))
    colnames(df_stability_surrogs)<-c("cvsq_ntd_0.025CI",  # 2.5%quantile of square of CVs from the surrogates having no tail dep. 
                                                              #     (possibly with Pearson preserving)
                                      "cvsq_ntd_median",  # median of square of CVs from the surrogates having no tail dep. 
                                      #     (possibly with Pearson preserving)
                                      "cvsq_ntd_0.975CI",  # 97.5% CI of CVs from the surrogates having no tail dep. 
                                      #     (possibly with Pearson preserving)
                                      "skw_ntd_0.025CI",   # 2.5%quantile of skewness from the surrogates having no tail dep. 
                                      #             (possibly with Pearson preserving)
                                      "skw_ntd_median",   # median of skewness from the surrogates having no tail dep. 
                                      #             (possibly with Pearson preserving)
                                      "skw_ntd_0.975CI"   # 97.5%quantile of skewness from the surrogates having no tail dep. 
                                      #             (possibly with Pearson preserving)
                              
    )
    
    #initialize
    cvsq_surrogs<-c()
    skw_surrogs<-c()
    
    for(i in c(1:dim(surrogs)[3])){
      ms<-surrogs[,,i]
      tot_quantity_s<-apply(ms,MARGIN = 1,FUN = sum)
      var_each_sp_s<-apply(ms,MARGIN = 2,FUN = var)
      mean_each_sp_s<-apply(ms,MARGIN = 2,FUN = mean)
      cm3_each_sp_s<-apply(ms,MARGIN = 2,FUN = my3cm)
      cvsq_surrogs<-c(cvsq_surrogs,mycvsq(tot_quantity_s))
      skw_surrogs<-c(skw_surrogs,myskns(tot_quantity_s))
    }
    
    CI_cvsq<-quantile(cvsq_surrogs,probs=c(.025,.975))
    df_stability_surrogs$cvsq_ntd_0.025CI<-CI_cvsq[1]
    df_stability_surrogs$cvsq_ntd_0.975CI<-CI_cvsq[2]
    df_stability_surrogs$cvsq_ntd_median<-median(cvsq_surrogs)
    
    CI_skw<-quantile(skw_surrogs,probs=c(.025,.975))
    df_stability_surrogs$skw_ntd_0.025CI<-CI_skw[1]
    df_stability_surrogs$skw_ntd_0.975CI<-CI_skw[2]
    df_stability_surrogs$skw_ntd_median<-median(skw_surrogs)
    
    df_stability<-cbind(df_stability,df_stability_surrogs)
    return(list(df_stability=df_stability,
                cvsq_surrogs=cvsq_surrogs,
                skw_surrogs=skw_surrogs))
  }else{
    return(df_stability)
  }
}

#----------------------------------------------------------------------------------

















#----------------------------------------this is just to check--------------------------------------
#x<-m[,c(1:2)]
#tot_quantity<-apply(x,MARGIN = 1,FUN = sum)
#v_tot<-var(tot_quantity)

#sum_vij<-cov(x=x[,1],y=x[,1])+cov(x=x[,1],y=x[,2])+cov(x=x[,2],y=x[,1])+cov(x=x[,2],y=x[,2])

#(v_tot-sum_vij)<1e-12 #they are same

#sum_vii<-cov(x=x[,1],y=x[,1])+cov(x=x[,2],y=x[,2])
#vii<-var(x[,1])+var(x[,2])

#sum_vii==vii # they are same











