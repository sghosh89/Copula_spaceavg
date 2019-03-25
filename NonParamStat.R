# ------------------------------------------------------------------------------------------------
# THIS CODE IS WRITTEN TO TEST NON-PARAMETRIC STATS FOR ANY TWO POSITIVELY CORRELATED SPECIES PAIR
#    ..............to maintain consistency with model selection approach.
#--------------------------------------------------------------------------------------------------
source("CopulaFunctions.R")
source("vivj_matrix.R")
source("good_splist.R")

#---------------------------------------------------------------
#Processing function for a copula approach to synchrony
#
#---Input---------
# m : output from vivj_matrix.R
#
#------Output - A list with these elements-----------
# ranks :       A dataframe with 2 columns, one for 
#               ds1 and one for ds2, corresponding to samples from the copula.
# spear   :     Spearman correlation (single number)
# kend      :   Kendall correlation (single number)
# Corl,Coru  :      covariance based statistics
# Pl,Pu    :  statistics show how the points within a copula scattered from right diagonal of the box within a pentagon
# Tl,Tu    :  statistics show how the points within a copula scattered along right diagonal of the box within a triangle
# D2l,D2u :   measures the squared distance between the points from the right diagonal within a copula 

copsync<-function(m,nbin){

  vi<-m[,1]
  vj<-m[,2]
  
  #get mean and variance
  vi_mean<-mean(vi)
  vj_mean<-mean(vj)
  var_vi<-var(vi)
  var_vj<-var(vj)
  
  if (length(vi)>0){
    #get spear
    spear<-cor(vi,vj, method ="pearson") 
    
    #get kend
    kend<-cor(vi, vj, method ="kendall") 
    
    #----------------------------------------------------STATISTICS :2 ---------------------------------------------------------
    #get Corl, Coru   (covariance based new stat)
    stat2<-CorlCoru(vi,vj,nbin=nbin)
    Corl<-stat2[1]
    Coru<-stat2[2]
    
    #---------------------------------------------------------- STATISTICS : 4 -----------------------------------------------------------------
    #   get New statistics : Pl and Pu   # distance from right diagonal in lower and upper triangle based stat
    stat4<-PlPu(vi,vj,nbin=nbin)
    Sl_Su_Si_P<-stat4[[1]]
    Pl<-stat4[[2]]
    Pu<-stat4[[3]]
    
    #--------------------------------------------------------- STATISTICS : 6 -----------------------------------------------------------------------
    # get D2l : average of squared distance of points from the right diagonal of the box for lower triangle
    # get D2u : average of squared distance of points from the right diagonal of the box for upper triangle
    stat6<-D2lD2u(vi,vj,nbin=nbin)
    D2l<-stat6[1]
    D2u<-stat6[2]
    
  }else{
    
    spear<-NA
    kend<-NA
    Corl<-NA
    Coru<-NA
    Sl_Su_Si_P<-NA
    Pl<-NA
    Pu<-NA
    D2l<-NA
    D2u<-NA
    
  }
  
  
  return(list(ranks=data.frame(Rki=vi,Rkj=vj),
              spear=spear,kend=kend,
              Corl=Corl,Coru=Coru,
              Sl_Su_Si_P=Sl_Su_Si_P,Pl=Pl,Pu=Pu,
              D2l=D2l,D2u=D2u))
}

#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
#Calling the above on all pairs of several time series and plotting
#results and returning all stats. All this function does is call the 
#above function on all pairs of time series, organize the results, 
#and make plots.
#
#-------------Input---------------------------------
# d_allsp  :    A list of data frames, each with columns Year and Dat
#               The years are assumed to be sequential and all included,
#               though there may be NAs in Dat and the years may not
#               be all the same for ds1 and ds2.
# loc  :        location number
# pfname   :    Filename (without extension) prepended to plot files saved.
# good_sp  :    a vector of chosen species

# nvar    :     an integer : the number of drivers for css # omitted at present
# timeavg : a logical tag whether you want to use timeavg data (if TRUE) or spatialavg data (if FALSE)

#----------Output - A list with these elements-----------------
# D                A matrix of geographic distances between sampling locations
# spear            A matrix of spearman results, length(d) by length(d)
# kend             A matrix of kendall results, length(d) by length(d)
# Corl             A matrix of Cl results, length(d) by length(d)
# Coru             A matrix of Cu results, length(d) by length(d)
# Pl               A matrix of Shy_lt results, length(d) by length(d)
# Pu               A matrix of Shy_ut results, length(d) by length(d)
# D2l              A matrix of R_l results, length(d) by length(d)
# D2u              A matrix of R_u results, length(d) by length(d)
# numericdf        A dataframe containing all the statistical data
# summary_nLU_npa  A dataframe counting the numbers of Corl-Coru, Pl-Pu, D2u-D2l positive/negative
#                   only for sp-sp interaction

multcall<-function(d_allsp,loc,pfname,good_sp,nbin,timeavg){
  
  lensp<-length(good_sp)
  
  #first initialize result receptacles for the output
  spear<-matrix(NA,lensp,lensp)
  colnames(spear) <- paste("sp",good_sp, sep="")
  rownames(spear) <-paste("sp",good_sp, sep="")
  
  kend<-matrix(NA,lensp,lensp)
  colnames(kend) <- colnames(spear)
  rownames(kend) <-rownames(spear)
  
  Corl<-matrix(NA,lensp,lensp)
  colnames(Corl) <- colnames(spear)
  rownames(Corl) <-rownames(spear)
  
  Coru<-matrix(NA,lensp,lensp)
  colnames(Coru) <- colnames(spear)
  rownames(Coru) <-rownames(spear)
  
  Pl<-matrix(NA,lensp,lensp)
  colnames(Pl) <- colnames(spear)
  rownames(Pl) <-rownames(spear)
  
  Pu<-matrix(NA,lensp,lensp)
  colnames(Pu) <- colnames(spear)
  rownames(Pu) <-rownames(spear)
  
  D2l<-matrix(NA,lensp,lensp)
  colnames(D2l) <- colnames(spear)
  rownames(D2l) <-rownames(spear)
  
  D2u<-matrix(NA,lensp,lensp)
  colnames(D2u) <- colnames(spear)
  rownames(D2u) <-rownames(spear)
  
  
  #------------------- PLOT :  copula_for_all_sp pair ----------------
  pdf(paste(pfname,"_AllCops.pdf",sep=""),width=6*lensp, height=6*lensp)
  op<-par(mfrow=c(lensp,lensp),mar=c(3,3,3,3), mgp=c(1.5,0.5,0))
  
  for (ii in c(1:lensp)){
    for (jj in c(1:lensp)){
      #compute results
      i<-good_sp[ii]
      j<-good_sp[jj]
      #cat("i,j",i,j,"\n")
      
      #if(i!=j){
      ms<-vivj_matrix(d_allsp=d_allsp,loc=loc,i=i,j=j,level=0.05,ploton=F,timeavg=timeavg,tagon=F)
      m<-ms$mat
      thisres<-copsync(m,nbin=nbin)
      
      spear[ii,jj]<-thisres$spear
      kend[ii,jj]<-thisres$kend
      
      Corl[ii,jj]<-thisres$Corl
      Coru[ii,jj]<-thisres$Coru
      
      Pl[ii,jj]<-thisres$Pl
      Pu[ii,jj]<-thisres$Pu
      
      D2l[ii,jj]<-thisres$D2l
      D2u[ii,jj]<-thisres$D2u
      
      
      plot(thisres$ranks$Rki,thisres$ranks$Rkj,type='p',col=rgb(0,0,0,.2),pch=19,xlim=c(0,1),ylim=c(0,1),xlab=expression(u[i]),ylab=expression(v[j]),cex.lab=2)
      mtext(paste0("[ i, j ] ="," [",i,",",j,"] ", ","," n=",dim(thisres$ranks)[1]),side = 3, line=0.15, adj=0.5, col="red")
    #}
  }
 }
  par(op)
  dev.off()
  
  
  #------------------- PLOT :  (Sl, Su and Si)_P for_all_species pair ----------------
  pdf(paste(pfname,"_Sl_Su_Si_P.pdf",sep=""),width=6*lensp, height=6*lensp)
  op<-par(mfrow=c(lensp,lensp),mar=c(3,3,3,3), mgp=c(1.5,0.5,0))
  
  for (ii in 1:lensp){
    for (jj in 1:lensp){
      
      i<-good_sp[ii]
      j<-good_sp[jj]
      
      #if(i!=j){
        
        ms<-vivj_matrix(d_allsp=d_allsp,loc=loc,i=i,j=j,level=0.05,ploton=F,timeavg=timeavg,tagon=F)
        m<-ms$mat
        thisres<-copsync(m,nbin=nbin)
        
        if(is.na(thisres$Sl_Su_Si_P$Sl_P[1]) == F){    
          plot(thisres$Sl_Su_Si_P$dist_Sl_P,thisres$Sl_Su_Si_P$Sl_P,type='l',col="red",ylim=c(0,1),xlab=" ",ylab=" ")
          lines(thisres$Sl_Su_Si_P$dist_Su_P,thisres$Sl_Su_Si_P$Su_P,type='l',col="blue",ylim=c(0,1))
          lines(thisres$Sl_Su_Si_P$dist_Si_P,thisres$Sl_Su_Si_P$Si_P,type='l',lty="dashed",col="green4",ylim=c(0,1))
          mtext(paste0("[ i, j ] ="," [",i,",",j,"], "," red : Sl_P, "," blue : Su_P, ", " green : Si_P"),side = 3, line=0.15, adj=0.5, col="black")
        }else{
          plot(0,0,type='p',col="white",cex=0,xlim=c(0,sqrt(2)/2),ylim=c(0,1),xlab=" ",ylab=" ")
          mtext(paste0("[ i, j ] ="," [",i,",",j,"]," ),side = 3, line=0.15, adj=0.5, col="black")
          mtext(paste0("NA" ),side = 3, line=-22.5, adj=0.5, col="black", cex=4)
        }
      #}
      
    }
  }
  
  par(op)
  dev.off()
  #-------------------------------------------------------------------------
  #-------------------!!!! ATTENTION !!!!---------------------
  #summary_nLU_npa<-matrix(NA,2,3)
  #rownames(summary_nLU_npa)<-c("num_LT","num_UT")
  #colnames(summary_nLU_npa)<-c("CorlmCoru","PlmPu","D2umD2l")
  #nI<-dim(posnI_list_ln_all[[loc]])[1] #number of independence 
  #nN<-dim(posnN_list_ln_all[[loc]])[1] #number of negcor
  #summary_nILU_npa[1,]<-nI
  #x<-Corl-Coru
  #y<-x[1:(nrow(x)-nvar),1:(ncol(x)-nvar)]
  #summary_nLU_npa[1,1]<-sum(y>0,na.rm=T)
  #summary_nLU_npa[2,1]<-sum(y<0,na.rm=T)
  
  #x<-Pl-Pu
  #y<-x[1:(nrow(x)-nvar),1:(ncol(x)-nvar)]
  #summary_nLU_npa[1,2]<-sum(y>0,na.rm=T)
  #summary_nLU_npa[2,2]<-sum(y<0,na.rm=T)
  
  #x<-D2u-D2l
  #y<-x[1:(nrow(x)-nvar),1:(ncol(x)-nvar)]
  #summary_nLU_npa[1,3]<-sum(y>0,na.rm=T)
  #summary_nLU_npa[2,3]<-sum(y<0,na.rm=T)
  
  #summary_nILU_npa[4,]<-nN
  #----------------------------------------------------------------------
  #-------------------------------------------------------------------------
  
  return(list(spear=spear,kend=kend,
              Corl=Corl,Coru=Coru,
              Pl=Pl,Pu=Pu,
              D2l=D2l,D2u=D2u))
}

#----------------------------------------------------------------------------------------------------------------------------
#                                                       CODE ENDS HERE
#-----------------------------------------------------------------------------------------------------------------------------









