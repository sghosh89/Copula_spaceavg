# This function is written to generate matrix plot for non-parametric stat (Cor stat only) results for each locations
# Input : 
#   data :  output from NonParamStat.R 
#   resloc : location to save the results
#   tagon : logical (argument for vivj_matrix fn call)
#   type : type argument of function mycorrplot : "lower" or "upper" or "full"
#   wd,ht : width and height of genarated plot

#Output
# A (Corl - Coru) matrix
#---------------------------
source("mycorrplot.R")
#---------------------------

NonParamStat_matrixplot<-function(data,resloc,tagon,type,wd,ht){
  
  
  #--------------------------Spearman plot---------------------------
  tempo<-data$spear
  indI<-data$posnI
  tempo[indI]<-NA
  diag(tempo)<-NA
  
  minval<-min(tempo,na.rm=T)
  maxval<-max(tempo,na.rm=T)
  cr<-max(abs(minval),abs(maxval))
  
  pdf(paste(resloc,file="Spearman.pdf",sep=''),width=wd, height=ht)
  z<-tempo
    
  colnames(z)<-rownames(z)
  mycorrplot(z=z,
             posnI_ind=data$posnI,
             posnN_ind=data$posnN,
             colrange=c(0,cr),type=type)
   
  dev.off()
  
  #--------------------------Kendall plot---------------------------
  
  tempo<-data$kend
  indI<-data$posnI
  tempo[indI]<-NA
  diag(tempo)<-NA
  
  minval<-min(tempo,na.rm=T)
  maxval<-max(tempo,na.rm=T)
  cr<-max(abs(minval),abs(maxval))
  
  pdf(paste(resloc,file="Kendall.pdf",sep=''),width=wd, height=ht)
  z<-tempo
  
  colnames(z)<-rownames(z)
  mycorrplot(z=z,
             posnI_ind=data$posnI,
             posnN_ind=data$posnN,
             colrange=c(0,cr))
  
  dev.off()
  
  #========================================= For cor npa stats ===============================================
  
  #if(npa_stats=="cor"){
    #--------------------------Corl plot---------------------------
    tempo<-data$Corl
    indI<-data$posnI
    tempo[indI]<-NA
    diag(tempo)<-NA
    
    minval<-min(tempo,na.rm=T)
    maxval<-max(tempo,na.rm=T)
    cr<-max(abs(minval),abs(maxval))
    
    pdf(paste(resloc,file="Corl.pdf",sep=''),width=wd, height=ht)
    z<-tempo
    
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=data$posnI,
               posnN_ind=data$posnN,
               colrange=c(0,cr))
    
    dev.off()
    
    #--------------------------Coru plot---------------------------
    tempo<-data$Coru
    indI<-data$posnI
    tempo[indI]<-NA
    diag(tempo)<-NA
    
    minval<-min(tempo,na.rm=T)
    maxval<-max(tempo,na.rm=T)
    cr<-max(abs(minval),abs(maxval))
    
    pdf(paste(resloc,file="Coru.pdf",sep=''),width=wd, height=ht)
    z<-tempo
    
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=data$posnI,
               posnN_ind=data$posnN,
               colrange=c(0,cr))
    
    dev.off()
    
    #--------------------------Corl-Coru plot---------------------------
    tempo<-data$Corl-data$Coru
    indI<-data$posnI
    indN<-data$posnN
    tempo[indI]<-NA
    diag(tempo)<-NA
    
    CorlmCoru<-tempo
    
    minval<-min(tempo,na.rm=T)
    maxval<-max(tempo,na.rm=T)
    cr<-max(abs(minval),abs(maxval))
    
    pdf(paste(resloc,file="Corl-Coru.pdf",sep=''),width=wd, height=ht)
    z<-tempo
    
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=data$posnI,
               posnN_ind=data$posnN,
               colrange=c(-cr,cr),type=type)
    
    z[data$posnN]<-NA # this line was added to exclude -vely correlated species pair from nL,nU
                                    # calculation, but it does not matter as for -vely correlated cells [sp_i,sp_j] and 
                                                  # [sp_j,sp_i] nL,nU both will increase by same number
    nL<-sum(z>0,na.rm = T)
    nU<-sum(z<0,na.rm = T)
    
    if(type=="lower" | type=="upper"){
      nL<-nL/2 # as z is a symmetric matrix
      nU<-nU/2
    }
    
    if(tagon == T){
      mtext(paste0("nL =",nL,", nU =",nU),cex=3,side=1,adj=0.6,line=2)
    }
    
    dev.off()
 # }
  
  return(CorlmCoru)
}
