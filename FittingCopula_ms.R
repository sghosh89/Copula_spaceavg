# THIS CODE FITS MULTIVARIATE COPULAS (PAIRWISE) BASED ON MODEL SELECTION APPROACH
# This code is written to test synchrony between different species pair at a specific location
#--------------------------------------------------------------------------------------------------------
library(copula)
library(VineCopula)
#-------------------------------------
source("vivj_matrix.R")
source(file="OurBiCopSelect.R")
#---------------------------------------------------------------------------------------
# This function generates result for single location
# Input :
#        loc : location index
#        d_allsp : given data in data[[loc]][[sp]] format to find cross species synchrony
#        good_sp : output from good_splist function or you can give a specified vector
#        families : vector of families of chosen set of copulas
#        level : Indep-test p-value limit
#        timeavg : a logical tag whether you want to use timeavg data (if TRUE) or spatialavg data (if FALSE)

RES_single_loc<-function(loc,d_allsp,good_sp,families,level,timeavg){
  
  len<-length(families)
  
  lengoodsp<-length(good_sp)

  #---------------------------------------------
  # Initialization to store output data
  gfc_p_CvM<-matrix(NA,nrow=lengoodsp,ncol=lengoodsp)
  colnames(gfc_p_CvM) <- paste("sp",good_sp, sep="")
  rownames(gfc_p_CvM) <- paste("sp",good_sp, sep="")
  
  gfc_p_KS<-gfc_p_CvM
  gfc_p_CvM_stat<-gfc_p_CvM
  gfc_p_KS_stat<-gfc_p_CvM
  gfc_normal_p_CvM<-gfc_p_CvM
  gfc_normal_p_KS<-gfc_p_CvM
  gfc_normal_p_CvM_stat<-gfc_p_CvM
  gfc_normal_p_KS_stat<-gfc_p_CvM
  copdata_cor_Kend<-gfc_p_CvM   #to store tauval of Kendall's correlation (either +ve or -ve)
  copdata_cor_Kend_pval<-gfc_p_CvM #to store p-value of Kendall's correlation (sig. cor. means <0.05)
  copdata_indeptest_pval<-gfc_p_CvM #to store p-value of BiCopIndepTest (sig. cor. means <0.05)
  gfc_numBS<-gfc_p_CvM
  gfc_numBS_success<-gfc_p_CvM
  gfc_normal_numBS<-gfc_p_CvM
  gfc_normal_numBS_success<-gfc_p_CvM
  LTdep_AICw<-gfc_p_CvM
  UTdep_AICw<-gfc_p_CvM
  LTmUTdep<-gfc_p_CvM #to store LTdep_AICw - UTdep_AICw for pos corr. and their absolute val for -ve corr.
  neg_sp_mat<-gfc_p_CvM # to store which sp taken as negative if neg cor found btw sp-pair
  
  info_ord_AIC<-array(NA,dim =c(lengoodsp,lengoodsp,len))  # To store all goodfit possibilities
  colnames(info_ord_AIC)<-colnames(gfc_p_CvM)
  rownames(info_ord_AIC)<-rownames(gfc_p_CvM)
  
  info_ord_copcode<-info_ord_AIC
  info_ord_copname<-info_ord_AIC
  info_ord_LTdep<-info_ord_AIC
  info_ord_UTdep<-info_ord_AIC
  info_ord_AICw<-info_ord_AIC
#--------------------------------------- 
  num_indep<-0

  for(i in c(1:lengoodsp)){
    for(j in c(1:lengoodsp)){
      if(i!=j){
        ms<-vivj_matrix(d_allsp=d_allsp,loc=loc,i=good_sp[i],j=good_sp[j],level=level,
                        ploton=F,timeavg=timeavg,tagon=F)
        m<-ms$mat
        u1<-m[,1]
        u2<-m[,2]
        
        copdata_cor_Kend[i,j]<-ms$tauval
        copdata_cor_Kend_pval[i,j]<-ms$pval
        
        #plot(u1,u2)
       #cat("[i, j]=[",i," ,",j,"]","(good_sp[i],good_sp[j])=(",good_sp[i],",",good_sp[j],")","\n")
        
#        print("About to call OurBiCopGofTest")
        ans<-OurBiCopSelect(u1=u1,u2=u2,families,level=0.05,AICBIC="AIC",
                            numBSsmall=100,pthresh=0.2,numBSlarge=1000,
                            gofnormal=FALSE,status=FALSE)
#        print("Finished OurBiCopGofTest")
        copdata_indeptest_pval[i,j]<-ans$IndepTestRes
        #------------------------------------------------------
        if(ans$IndepTestRes<level && ms$tauval<0){
          neg_sp_mat[i,j]<- good_sp[j] # this sp is taken reverse and plot @yaxis in copula
        }
        #------------------------------------------------------
        if(ans$IndepTestRes<level){
          
          gfc_numBS[i,j]<-ans$Numboot
          gfc_numBS_success[i,j]<-ans$Numboot_success
          gfc_p_CvM[i,j]<-ans$GofRes_CvM
          gfc_p_KS[i,j]<-ans$GofRes_KS
          gfc_p_CvM_stat[i,j]<-ans$GofRes_CvM_stat
          gfc_p_KS_stat[i,j]<-ans$GofRes_KS_stat
          
          gfc_normal_numBS[i,j]<-ans$Numboot_Normal
          gfc_normal_numBS_success[i,j]<-ans$Numboot_success_Normal
          gfc_normal_p_CvM[i,j]<-ans$GofRes_Normal_CvM
          gfc_normal_p_KS[i,j]<-ans$GofRes_Normal_KS
          gfc_normal_p_CvM_stat[i,j]<-ans$GofRes_Normal_CvM_stat
          gfc_normal_p_KS_stat[i,j]<-ans$GofRes_Normal_KS_stat
          
          LTdep_AICw[i,j]<-ans$relLTdep_AICw
          UTdep_AICw[i,j]<-ans$relUTdep_AICw
          
          info<-ans$InfCritRes
          info_ord<-info[order(info[,6],decreasing = F),] # ordered according to min to max AIC
          info_ord_AIC[i,j,]<-info_ord$AIC
          info_ord_copcode[i,j,]<-info_ord$copcode
          info_ord_copname[i,j,]<-as.character(info_ord$copname)
          info_ord_LTdep[i,j,]<-info_ord$LTdep
          info_ord_UTdep[i,j,]<-info_ord$UTdep
          info_ord_AICw[i,j,]<-info_ord$AICw
          
          #cat("IndepTestRes=",ans$IndepTestRes,"\n")
          #cat("tauval=",ms$tauval,"\n")
          
          #if(ms$tauval >0){
            LTmUTdep[i,j]<-(LTdep_AICw[i,j] - UTdep_AICw[i,j])
          #}
          
          #if(ms$tauval <0){
          #  LTmUTdep[i,j]<-abs(LTdep_AICw[i,j] - UTdep_AICw[i,j])
          #}
 
        
        }else{
          num_indep<-num_indep+1
          gfc_numBS[i,j]<-Inf
          gfc_numBS_success[i,j]<-Inf
          gfc_p_CvM[i,j]<-Inf         # I just put Inf to see when it's indep?
          gfc_p_KS[i,j]<-Inf
          gfc_p_CvM_stat[i,j]<-Inf         # I just put Inf to see when it's indep?
          gfc_p_KS_stat[i,j]<-Inf
          gfc_normal_numBS[i,j]<-Inf
          gfc_normal_numBS_success[i,j]<-Inf
          gfc_normal_p_CvM[i,j]<-Inf  # If gofnormal==F, then gfc_rmal_p matrices also contains Inf at all off-diagonal entries
          gfc_normal_p_KS[i,j]<-Inf
          gfc_normal_p_CvM_stat[i,j]<-Inf
          gfc_normal_p_KS_stat[i,j]<-Inf
          info_ord_AIC[i,j,]<-Inf
          info_ord_copcode[i,j,]<-Inf
          info_ord_copname[i,j,]<-"Indep"
          info_ord_LTdep[i,j,]<-Inf
          info_ord_UTdep[i,j,]<-Inf
          info_ord_AICw[i,j,]<-Inf
          LTdep_AICw[i,j]<-Inf
          UTdep_AICw[i,j]<-Inf
          neg_sp_mat[i,j]<-Inf
        }
        
      }
    }
    
  }
  #-----------------------------------------------------------------------------------------------
  posnI_ind<-which(copdata_indeptest_pval>=level, arr.ind = T) #indices of indep. pair
  posnP_ind<-which(copdata_indeptest_pval<level & copdata_cor_Kend >0, arr.ind = T) #indices of significantly pos correlated pair
  posnN_ind<-which(copdata_indeptest_pval<level & copdata_cor_Kend <0, arr.ind = T) #indices of significantly neg correlated pair
  
  # Save the results
  RES_loc<-list(num_indep_sp_pair=num_indep/2,
               num_neg_cor_sp_pair=nrow(posnN_ind)/2,
               posnI_ind=posnI_ind,
               posnN_ind=posnN_ind,
               posnP_ind=posnP_ind,
               gfc_numBS=gfc_numBS,
               gfc_numBS_success=gfc_numBS_success,
               gfc_p_CvM=gfc_p_CvM,
               gfc_p_KS=gfc_p_KS,
               gfc_p_CvM_stat=gfc_p_CvM_stat,
               gfc_p_KS_stat=gfc_p_KS_stat,
               gfc_normal_numBS=gfc_normal_numBS,
               gfc_normal_numBS_success=gfc_normal_numBS_success,
               gfc_normal_p_CvM=gfc_normal_p_CvM,
               gfc_normal_p_KS=gfc_normal_p_KS,
               gfc_normal_p_CvM_stat=gfc_normal_p_CvM_stat,
               gfc_normal_p_KS_stat=gfc_normal_p_KS_stat,
               info_ord_copcode=info_ord_copcode,
               info_ord_copname=info_ord_copname,
               info_ord_AIC=info_ord_AIC,
               info_ord_LTdep=info_ord_LTdep,
               info_ord_UTdep=info_ord_UTdep,
               info_ord_AICw=info_ord_AICw,
               LTdep_AICw=LTdep_AICw,
               UTdep_AICw=UTdep_AICw,
               LTmUTdep=LTmUTdep,
               neg_sp_mat=neg_sp_mat,
               copdata_cor_Kend=copdata_cor_Kend,
               copdata_cor_Kend_pval=copdata_cor_Kend_pval,
               copdata_indeptest_pval=copdata_indeptest_pval)
  
  return(RES_loc)
  
} 

#--------------------------------------------------------------------------











