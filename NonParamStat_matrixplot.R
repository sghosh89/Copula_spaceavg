# This function is written to generate matrix plot for non-parametric stat (Cor stat only) results for each locations
# Input : 
#   data :  output from NonParamStat.R 
#   resloc : location to save the results
#   tagon : logical (argument for vivj_matrix fn call)
#   type : type argument of function mycorrplot_with_sig : "lower" or "upper" or "full"
#   wd,ht : width and height of genarated plot
#   for these following arguments: see  arguments from mycorrplot_with_sig.R 
#   sigtest,ub,numpts,numsims=10000,CI=c(0.025,0.975),include_indep

#Output : a list of two:
#   1. A (Corl - Coru) matrix
#   2. A list of two : A data frame having sig. test info and a matrix with significant marking
#---------------------------
source("mycorrplot_with_sig.R")
#---------------------------

NonParamStat_matrixplot<-function(data,resloc,tagon,type,wd,ht,sigtest,ub,numpts,numsims=10000,CI=c(0.025,0.975),include_indep){
  
  #--------------------------Spearman plot---------------------------
  tempo<-data$spear
  indI<-data$posnI
  tempo[indI]<-NA
  diag(tempo)<-NA
  
  minval<-min(tempo,na.rm=T)
  maxval<-max(tempo,na.rm=T)
  cr<-max(abs(minval),abs(maxval))
  
  pdf(paste(resloc,file="Spearman_ub_",ub,".pdf",sep=''),width=wd, height=ht)
  z<-tempo
    
  colnames(z)<-rownames(z)
  mycorrplot_with_sig(z=z,
             posnI_ind=data$posnI,
             posnN_ind=data$posnN,
             colrange=c(0,cr),type=type,
             sigtest=F,spr=NA,realstat=NA,
             ub=NA,numpts=NA,numsims=NA,CI=NA,include_indep=include_indep)
   
  dev.off()
  
  #--------------------------Kendall plot---------------------------
  
  # NOTE: for very week -ve spearman correlated cell can be +vely kendall correlated, however, 
  # on significance test they will show independence relationship.
  
  #tempo<-data$kend
  #indI<-data$posnI
  #tempo[indI]<-NA
  #diag(tempo)<-NA
  
  #minval<-min(tempo,na.rm=T)
  #maxval<-max(tempo,na.rm=T)
  #cr<-max(abs(minval),abs(maxval))
  
  #pdf(paste(resloc,file="Kendall_ub_",ub,".pdf",sep=''),width=wd, height=ht)
  #z<-tempo
  
  #colnames(z)<-rownames(z)
  #mycorrplot_with_sig(z=z,
  #           posnI_ind=data$posnI,
  #           posnN_ind=data$posnN,
  #           colrange=c(0,cr),type=type,
  #          sigtest=F,spr=NA,realstat=NA,
  #          ub=NA,numpts=NA,numsims=NA,CI=NA,include_indep=include_indep)
  
  #dev.off()
  
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
    
    pdf(paste(resloc,file="Corl_ub_",ub,".pdf",sep=''),width=wd, height=ht)
    z<-tempo
    
    colnames(z)<-rownames(z)
    
    if(include_indep==T){
      colrange<-c(0,cr)
    }else{
      colrange<-c(-cr,cr)
    }
    
    mycorrplot_with_sig(z=z,
               posnI_ind=data$posnI,
               posnN_ind=data$posnN,
               colrange=colrange,type=type,
               sigtest=F,spr=NA,realstat=NA,
               ub=NA,numpts=NA,numsims=NA,CI=NA,include_indep=include_indep)
    
    dev.off()
    
    #--------------------------Coru plot---------------------------
    tempo<-data$Coru
    indI<-data$posnI
    tempo[indI]<-NA
    diag(tempo)<-NA
    
    minval<-min(tempo,na.rm=T)
    maxval<-max(tempo,na.rm=T)
    cr<-max(abs(minval),abs(maxval))
    
    pdf(paste(resloc,file="Coru_ub_",ub,".pdf",sep=''),width=wd, height=ht)
    z<-tempo
    
    colnames(z)<-rownames(z)
    
    if(include_indep==T){
      colrange<-c(0,cr)
    }else{
      colrange<-c(-cr,cr)
    }
    
    mycorrplot_with_sig(z=z,
               posnI_ind=data$posnI,
               posnN_ind=data$posnN,
               colrange=colrange,type=type,
               sigtest=F,spr=NA,realstat=NA,
               ub=NA,numpts=NA,numsims=NA,CI=NA,include_indep=include_indep)
    
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
    
    pdf(paste(resloc,file="Corl-Coru_ub_",ub,".pdf",sep=''),width=wd, height=ht)
    z<-tempo
    
    colnames(z)<-rownames(z)
    
    if(sigtest==T){
      res_sig<-mycorrplot_with_sig(z=z,
                          posnI_ind=data$posnI,
                          posnN_ind=data$posnN,
                          colrange=c(-cr,cr),type=type,
                          sigtest=T,spr=data$spear,realstat=CorlmCoru,
                          ub=ub,numpts=numpts,numsims=numsims,CI=CI,include_indep=include_indep)
    }else{
      mycorrplot_with_sig(z=z,
                          posnI_ind=data$posnI,
                          posnN_ind=data$posnN,
                          colrange=c(-cr,cr),type=type,
                          sigtest=F,spr=NA,realstat=NA,
                          ub=NA,numpts=NA,numsims=NA,CI=NA,include_indep=include_indep)
      res_sig<-NA
    }
   
    z[data$posnN]<-NA # this line was added to exclude -vely correlated species pair from nL,nU
                                    # calculation, but it does not matter as for -vely correlated cells [sp_i,sp_j] and 
                                                  # [sp_j,sp_i] nL,nU both will increase by same number
    nL<-sum(z>0,na.rm = T)
    nU<-sum(z<0,na.rm = T)
    total_CorlmCoru<-sum(z,na.rm=T) # for positively correlated cells only
    
    if(isSymmetric(z)==T){ #it's a check
      if(type=="lower" | type=="upper"){
        nL<-nL/2 
        nU<-nU/2
        total_CorlmCoru<-total_CorlmCoru/2
      }
    }
    
    if(tagon == T){
      mtext(paste0("nL = ",nL,", nU = ",nU, ", Total asym. = ",round(total_CorlmCoru,4)),cex=3,side=1,adj=0.6,line=2)
    }
    
    dev.off()
    
    if(sigtest==T){
      
     mat_tab<-res_sig$mat_tab
     
     # generate additional plot
     pdf(paste(resloc,"statistic_vs_spearman_ub_",ub,"_CI_",CI[1],"_",CI[2],".pdf",sep=""),width=8,height=8)
     plot(c(-1,1),c(0,0),ylim=c(-1,1),xlab="Spearman",ylab="Statistic",type='l',col='red')
     lines(c(0,0),c(-1,1),type="l",col="red")
     points(mat_tab$sprvals,mat_tab$realstat,pch=16,col=rgb(1,0,0,0.2))
     lines(mat_tab$sprvals[mat_tab$sprvals>0],mat_tab$lowCI[mat_tab$sprvals>0],type='p',pch=16,col=rgb(0,0,0,0.2))
     lines(mat_tab$sprvals,mat_tab$upCI,type='p',pch=16,col=rgb(0,0,0,0.2))
     dev.off()
  }
  
  return(list(CorlmCoru=CorlmCoru,
              res_sig=res_sig))
}
