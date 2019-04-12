#rm(list=ls())
source("mycvsq.R")
source("SkewnessAnd3CentMom.R")

# This function gives you a dataframe 
# ARGs : 
#           m = a matrix with sp. time series along each column
#           id_sp : a vector 

make_tab_stability<-function(m,id_sp){
  
  m<-m[,id_sp]
  tot_quantity<-apply(m,MARGIN = 1,FUN = sum)
  var_each_sp<-apply(m,MARGIN = 2,FUN = var)
  mean_each_sp<-apply(m,MARGIN = 2,FUN = mean)
  cm3_each_sp<-apply(m,MARGIN = 2,FUN = my3cm)
  
  df_stability<-as.data.frame(matrix(NA,nrow=1,ncol=8))
 
  colnames(df_stability)<-c("cvsq_real", # square of CV for the real community data
                            "cvsq_ntd",  # mean of square of CVs from the surrogates having no tail dep. (possibly with Pearson preserving)
                            "cvsq_indep", # square of CV of the community if all the sp. behave independently
                            "phi_cvsq",  # this is the ratio of cvsq_real/cvsq_indep and compared to 1
                            "skw_real",  # skewness for the real community data
                            "skw_ntd",   # mean of skewness from the surrogates having no tail dep. (possibly with Pearson preserving)
                            "skw_indep", # skewness of the community if all the sp. behave independently
                            "phi_skw" # this is the ratio of skw_real/skw_indep and compared to ??
  )
  
  df_stability$cvsq_real<-mycvsq(tot_quantity)
  df_stability$cvsq_ntd<-"???"
  df_stability$cvsq_indep<-sum(var_each_sp)/((sum(mean_each_sp))^2)
  df_stability$phi_cvsq<-df_stability$cvsq_real/df_stability$cvsq_indep
  
  df_stability$skw_real<-myskns(tot_quantity)
  df_stability$skw_ntd<-"???"
  df_stability$skw_indep<-sum(cm3_each_sp)/((sum(var_each_sp))^(3/2))
  df_stability$phi_skw<-df_stability$skw_real/df_stability$skw_indep
  
  return(df_stability)
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
jrg_C<-make_tab_stability(m=m,id_sp=id_sp)

m<-d_jrg
id_sp<-which(category_jrg$category %in% c("C","I"))
jrg_CI<-make_tab_stability(m=m,id_sp=id_sp)

m<-d_jrg
id_sp<-which(category_jrg$category %in% c("C","I","R"))
jrg_CIR<-make_tab_stability(m=m,id_sp=id_sp)


jrg_stability<-rbind(jrg_C,jrg_CI,jrg_CIR)
rownames(jrg_stability)<-c("C","C+I","all")

#=======================================

# Now call for hays data
m<-d_hays
id_sp<-which(category_hays$category %in% "C")
hays_C<-make_tab_stability(m=m,id_sp=id_sp)

m<-d_hays
id_sp<-which(category_hays$category %in% c("C","I"))
hays_CI<-make_tab_stability(m=m,id_sp=id_sp)

m<-d_hays
id_sp<-which(category_hays$category %in% c("C","I","R"))
hays_CIR<-make_tab_stability(m=m,id_sp=id_sp)


hays_stability<-rbind(hays_C,hays_CI,hays_CIR)
rownames(hays_stability)<-c("C","C+I","all")

#----------------------------------------this is just to check--------------------------------------
#x<-m[,c(1:2)]
#tot_quantity<-apply(x,MARGIN = 1,FUN = sum)
#v_tot<-var(tot_quantity)

#sum_vij<-cov(x=x[,1],y=x[,1])+cov(x=x[,1],y=x[,2])+cov(x=x[,2],y=x[,1])+cov(x=x[,2],y=x[,2])

#(v_tot-sum_vij)<1e-12 #they are same

#sum_vii<-cov(x=x[,1],y=x[,1])+cov(x=x[,2],y=x[,2])
#vii<-var(x[,1])+var(x[,2])

#sum_vii==vii # they are same











