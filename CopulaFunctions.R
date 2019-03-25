#-----------------------------------------------------------------------------------------------------------------------
# THIS CODE CONTAINS STATISTICS FOR LOOKING AT TAIL DEPENDENCE
#-----------------------------------------------------------------------------------------------------------------------

source("CopulaFunctions_flexible.R")

#--------------------------- correlation based Stat ---------------------------------------------------------

CorlCoru<-function(vi,vj,nbin){
  
  lb<-0
  ub<-1/nbin
  
  Corl<-Corbds(vi,vj,lb=lb,ub=ub)
  Coru<-Corbds(vi,vj,lb=1-ub,ub=1-lb)

  return(c(Corl,Coru))
}

#--------------------------- Stat based on counting points a given distance from the diagonal in lower and upper triangle ---------------------------------------------------------

PlPu<-function(vi,vj,nbin){
  
  lb<-0
  ub<-1/nbin
  
  ind_ll<-which(vi+vj<1)
  ind_ur<-which(vi+vj>1)
  
  if(length(ind_ll)!=0 & length(ind_ur)!=0){
    
    res_l<-Pbds(vi,vj,lb=lb,ub=ub)
    Pl<-res_l$abs_res
    dist_Sl_P<-res_l$dist_S
    Sl_P<-res_l$S
    dist_Si_P<-res_l$dist_Si
    Si_P<-res_l$Si
    
    res_u<-Pbds(vi,vj,lb=1-ub,ub=1-lb)
    Pu<-res_u$abs_res
    dist_Su_P<-res_u$dist_S
    Su_P<-res_u$S
    
    Sl_Su_Si_P<-list(dist_Sl_P=dist_Sl_P,Sl_P=Sl_P,dist_Su_P=dist_Su_P,Su_P=Su_P,dist_Si_P=dist_Si_P,Si_P=Si_P)
    
  }else{
    
    Sl_Su_Si_P<-list(dist_Sl_P=NA,Sl_P=NA,dist_Su_P=NA,Su_P=NA,dist_Si_P=NA,Si_P=NA)
    Pl<-NA
    Pu<-NA
    
  }
 
  return(list(Sl_Su_Si_P,Pl,Pu)) 
}

#------------------------------------------------ stat based on distance from the diagonal -----------------------------------------------------------------------

D2lD2u<-function(vi,vj,nbin){
  
  lb<-0
  ub<-1/nbin
  
  D2l<-D2bds(vi,vj,lb=lb,ub=ub)
  D2u<-D2bds(vi,vj,lb=1-ub,ub=1-lb)
  
  return(c(D2l,D2u))
}

#---------------------------------------------------------------------------------------------------------------------
#                                                     CODE ENDS HERE
#---------------------------------------------------------------------------------------------------------------------





















#and the others