#rm(list=ls())
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
                                      "cvsq_ntd_mean",  # mean of square of CVs from the surrogates having no tail dep. 
                                      #     (possibly with Pearson preserving)
                                      "cvsq_ntd_0.975CI",  # mean of square of CVs from the surrogates having no tail dep. 
                                      #     (possibly with Pearson preserving)
                                      "skw_ntd_0.025CI",   # 2.5%quantile of skewness from the surrogates having no tail dep. 
                                      #             (possibly with Pearson preserving)
                                      "skw_ntd_mean",   # mean of skewness from the surrogates having no tail dep. 
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
    df_stability_surrogs$cvsq_ntd_mean<-mean(cvsq_surrogs)
    
    CI_skw<-quantile(skw_surrogs,probs=c(.025,.975))
    df_stability_surrogs$skw_ntd_0.025CI<-CI_skw[1]
    df_stability_surrogs$skw_ntd_0.975CI<-CI_skw[2]
    df_stability_surrogs$skw_ntd_mean<-mean(skw_surrogs)
    
    df_stability<-cbind(df_stability,df_stability_surrogs)
    return(list(df_stability=df_stability,
                cvsq_surrogs=cvsq_surrogs,
                skw_surrogs=skw_surrogs))
  }else{
    return(df_stability)
  }
}

#---------------------------------------------------------
### Bring in the data

#Importing the two dataset matrices of years by species
d_jrg<- readRDS("./Results/jrg_results/skewness_results/ts_mat_all_sp_39_jrg.RDS")# for all species except "BARE","ROCK","LASP"
d_hays<- readRDS("./Results/hays_results/skewness_results/ts_mat_all_sp_146_hays.RDS")# for all species except in "REMOVE" category

#Importing data with "C", "I", "R" values for common, intermediate, and rare species
category_jrg<-readRDS("./Results/jrg_results/skewness_results/all_sp_39_jrg_category.RDS")
category_hays<-readRDS("./Results/hays_results/skewness_results/all_sp_146_hays_category.RDS")  

# Now call for jrg data
m<-d_jrg
id_sp<-which(category_jrg$category %in% "C")
m<-m[,id_sp]
jrg_C<-make_tab_stability(m=m,surrogs = NA,surrogs_given = F)

m<-d_jrg
id_sp<-which(category_jrg$category %in% c("C","I"))
m<-m[,id_sp]
jrg_CI<-make_tab_stability(m=m,surrogs = NA,surrogs_given = F)

m<-d_jrg
id_sp<-which(category_jrg$category %in% c("C","I","R"))
m<-m[,id_sp]
jrg_CIR<-make_tab_stability(m=m,surrogs = NA,surrogs_given = F)

jrg_stability<-rbind(jrg_C,jrg_CI,jrg_CIR)
rownames(jrg_stability)<-c("C","C+I","all")

write.csv(jrg_stability,"./Results/jrg_results/skewness_results/jrg_stability_C_CI_CIR.csv")

#--------------------------------------------------------
# Now for pearson preserving surrogates with C+pseudo sp.
m<-d_jrg
id_sp<-which(category_jrg$category %in% "C")
m<-m[,id_sp] # for common sp.
id_sp<-which(!category_jrg$category %in% "C")
m_pseudo<-d_jrg[,id_sp]
pseudo_sp<-apply(m_pseudo,MARGIN = 1,FUN = sum)
m<-cbind(m,pseudo_sp)

jrg_surrogs_pp_CP<-readRDS("./Results/Results_pps_exist_JRG_CommonAndPseudo/SomeSurrogates.RDS")
dim(jrg_surrogs_pp_CP)

numsurrog<-10000
#_________________________________________________________________
set.seed(seed=101) # set seed to sample any surrogates randomly
#_________________________________________________________________
id_surrogs<-sample(c(1:dim(jrg_surrogs_pp_CP)[3]),numsurrog,replace=F)

jrg_surrogs_pp_CP<-jrg_surrogs_pp_CP[,,id_surrogs]

jrg_CP<-make_tab_stability(m=m,surrogs = jrg_surrogs_pp_CP,surrogs_given = T)
ans<-(jrg_CP$df_stability)
rownames(ans)<-"C+P"
class(ans)

write.csv(ans,"./Results/jrg_results/skewness_results/jrg_stability_CP.csv")

#----------------------------------------------------
pdf("./Results/jrg_results/skewness_results/jrg_pearson_preserving_results_cvsq_skw_plots.pdf")
op<-par(mfrow=c(2,1))

#--------------CVsq histogrm-------------------------------------
xlm<-range(ans$cvsq_real,ans$cvsq_indep,jrg_CP$cvsq_surrogs)

hist(jrg_CP$cvsq_surrogs,col="grey",border=F,breaks=100,xaxt="n",xlim=xlm,xlab="10000 PP Surrogs CVsq",main="")
axis(side=1, at=c(xlm[1],0,xlm[2]))
abline(v=ans$cvsq_real,col="red") # actual CVsq from real data

#95% quantiles
abline(v=ans$cvsq_ntd_0.025CI,col="blue") 
abline(v=ans$cvsq_ntd_0.975CI,col="blue")

#50% quantiles
CI_cvsq_50<-quantile(jrg_CP$cvsq_surrogs,probs=c(.25,.75))
abline(v=CI_cvsq_50[1],col="green")
abline(v=CI_cvsq_50[2],col="green")

# Cvsq with no tail dep.
abline(v=ans$cvsq_ntd_mean,col="magenta")

# Cvsq if indep.
abline(v=ans$cvsq_indep,col="black")


#--------------skw histogrm-------------------------------------
xlm<-range(ans$skw_real,ans$skw_indep,jrg_CP$skw_surrogs)

hist(jrg_CP$skw_surrogs,col="grey",border=F,breaks=100,xaxt="n",xlim=xlm,xlab="10000 PP Surrogs Skewness",main="")
axis(side=1, at=c(xlm[1],0,xlm[2]))
abline(v=ans$skw_real,col="red") # actual CVsq from real data

#95% quantiles
abline(v=ans$skw_ntd_0.025CI,col="blue") 
abline(v=ans$skw_ntd_0.975CI,col="blue")

#50% quantiles
CI_skw_50<-quantile(jrg_CP$skw_surrogs,probs=c(.25,.75))
abline(v=CI_skw_50[1],col="green")
abline(v=CI_skw_50[2],col="green")

# Skewness with no tail dep.
abline(v=ans$skw_ntd_mean,col="magenta")

# Skewness if indep.
abline(v=ans$skw_indep,col="black")

# add legend
legend("topright",col=c("red","magenta","blue","green","black"),lty=c(1,1,1,1),
       horiz = F, bty="n",
       legend=c("real value","no Tail-dep. (mean)","95%CI","50%CI","Indep."))


par(op)
dev.off()
#======================================================================================================

# =============================================Now call for hays data=============================
m<-d_hays
id_sp<-which(category_hays$category %in% "C")
hays_C<-make_tab_stability(m=m,surrogs = NA,surrogs_given = F)

m<-d_hays
id_sp<-which(category_hays$category %in% c("C","I"))
hays_CI<-make_tab_stability(m=m,surrogs = NA,surrogs_given = F)

m<-d_hays
id_sp<-which(category_hays$category %in% c("C","I","R"))
hays_CIR<-make_tab_stability(m=m,surrogs = NA,surrogs_given = F)


hays_stability<-rbind(hays_C,hays_CI,hays_CIR)
rownames(hays_stability)<-c("C","C+I","all")

write.csv(hays_stability,"./Results/hays_results/skewness_results/hays_stability_C_CI_CIR.csv")

#--------------------------------------------------------
# Now for pearson preserving surrogates with C+pseudo sp.
m<-d_hays
id_sp<-which(category_hays$category %in% "C")
m<-m[,id_sp] # for common sp.
id_sp<-which(!category_hays$category %in% "C")
m_pseudo<-d_hays[,id_sp]
pseudo_sp<-apply(m_pseudo,MARGIN = 1,FUN = sum)
m<-cbind(m,pseudo_sp)

hays_surrogs_pp_CP<-readRDS("./Results/Results_pps_exist_HAY_CommonAndPseudo/SomeSurrogates.RDS")
dim(hays_surrogs_pp_CP)

numsurrog<-10000
#_________________________________________________________________
set.seed(seed=101) # set seed to sample any surrogates randomly
#_________________________________________________________________
id_surrogs<-sample(c(1:dim(hays_surrogs_pp_CP)[3]),numsurrog,replace=F)

hays_surrogs_pp_CP<-hays_surrogs_pp_CP[,,id_surrogs]

hays_CP<-make_tab_stability(m=m,surrogs = hays_surrogs_pp_CP,surrogs_given = T)
ans<-(hays_CP$df_stability)
rownames(ans)<-"C+P"
class(ans)

write.csv(ans,"./Results/hays_results/skewness_results/hays_stability_CP.csv")

#----------------------------------------------------
pdf("./Results/hays_results/skewness_results/hays_pearson_preserving_results_cvsq_skw_plots.pdf")
op<-par(mfrow=c(2,1))

#--------------CVsq histogrm-------------------------------------
xlm<-range(ans$cvsq_real,ans$cvsq_indep,hays_CP$cvsq_surrogs)

hist(hays_CP$cvsq_surrogs,col="grey",border=F,breaks=100,xaxt="n",xlim=xlm,xlab="10000 PP Surrogs CVsq",main="")
axis(side=1, at=c(xlm[1],0,xlm[2]))
abline(v=ans$cvsq_real,col="red") # actual CVsq from real data

#95% quantiles
abline(v=ans$cvsq_ntd_0.025CI,col="blue") 
abline(v=ans$cvsq_ntd_0.975CI,col="blue")

#50% quantiles
CI_cvsq_50<-quantile(hays_CP$cvsq_surrogs,probs=c(.25,.75))
abline(v=CI_cvsq_50[1],col="green")
abline(v=CI_cvsq_50[2],col="green")

# Cvsq with no tail dep.
abline(v=ans$cvsq_ntd_mean,col="magenta")

# Cvsq if indep.
abline(v=ans$cvsq_indep,col="black")


#--------------skw histogrm-------------------------------------
xlm<-range(ans$skw_real,ans$skw_indep,jrg_CP$skw_surrogs)

hist(hays_CP$skw_surrogs,col="grey",border=F,breaks=100,xaxt="n",xlim=xlm,xlab="10000 PP Surrogs Skewness",main="")
axis(side=1, at=c(xlm[1],0,xlm[2]))
abline(v=ans$skw_real,col="red") # actual CVsq from real data

#95% quantiles
abline(v=ans$skw_ntd_0.025CI,col="blue") 
abline(v=ans$skw_ntd_0.975CI,col="blue")

#50% quantiles
CI_skw_50<-quantile(hays_CP$skw_surrogs,probs=c(.25,.75))
abline(v=CI_skw_50[1],col="green")
abline(v=CI_skw_50[2],col="green")

# Skewness with no tail dep.
abline(v=ans$skw_ntd_mean,col="magenta")

# Skewness if indep.
abline(v=ans$skw_indep,col="black")

# add legend
legend("topright",col=c("red","magenta","blue","green","black"),lty=c(1,1,1,1),
       horiz = F, bty="n",
       legend=c("real value","no Tail-dep. (mean)","95%CI","50%CI","Indep."))

par(op)
dev.off()

























#----------------------------------------this is just to check--------------------------------------
#x<-m[,c(1:2)]
#tot_quantity<-apply(x,MARGIN = 1,FUN = sum)
#v_tot<-var(tot_quantity)

#sum_vij<-cov(x=x[,1],y=x[,1])+cov(x=x[,1],y=x[,2])+cov(x=x[,2],y=x[,1])+cov(x=x[,2],y=x[,2])

#(v_tot-sum_vij)<1e-12 #they are same

#sum_vii<-cov(x=x[,1],y=x[,1])+cov(x=x[,2],y=x[,2])
#vii<-var(x[,1])+var(x[,2])

#sum_vii==vii # they are same











