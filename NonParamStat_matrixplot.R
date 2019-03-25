# This function is written to generate matrix plot for non-parametric stat results for all locations
# Input : 
#   data_ln_all :        a list of length of all locations (output from NonParamStat.R for all locations)
#   posnI_list_ln_all :  A list of matrices for all locations
#                        each containing the position indices of independent
#   posnN_list_ln_all  : A list of matrices for all locations
#                        each containing the position indices of negatively correlated
#                        vi, vj based on FittingCopula_ms.R code (model selection) results
#   resloc : location to save the results
#   nvar : an integer : number of factors which influence css between considered species
#   nvar_names : a vector of characters containing names of n-variables
#   r : output from FittingCopula_ms.R which will give "nsm" to be used as input in mycorrplot fn call
#   tagon : logical 
#---------------------------
source("mycorrplot.R")
#---------------------------

NonParamStat_matrixplot<-function(data_ln_all,posnI_list_ln_all,posnN_list_ln_all,resloc,nvar,nvar_names,r,tagon){
  
  numloc<-length(data_ln_all)
  selected_loc<-names(data_ln_all)
  
  
  #--------------------------Spearman plot---------------------------
  
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$spear
    indI<-posnI_list_ln_all[[i]]
    tempo[[i]][indI]<-NA
  }
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_Spearman.pdf",sep=''),width=20, height=20)
    z<-tempo[[loc]]
    
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=posnI_list_ln_all[[loc]],
               posnN_ind=posnN_list_ln_all[[loc]],
               colrange=c(-cr,cr),nsm=r[[loc]]$neg_sp_mat)
    #nI<-dim(posnI_list_ln_all[[loc]])[1]
    #mtext(paste0(selected_loc[loc],", #nI =",nI),cex=2,side=1)
    dev.off()
  }
  
  #--------------------------Kendall plot---------------------------
  
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$kend
    indI<-posnI_list_ln_all[[i]]
    tempo[[i]][indI]<-NA
  }
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_Kendall.pdf",sep=''),width=20, height=20)
    z<-tempo[[loc]]
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=posnI_list_ln_all[[loc]],
               posnN_ind=posnN_list_ln_all[[loc]],
               colrange=c(-cr,cr),nsm=r[[loc]]$neg_sp_mat)
    #nI<-dim(posnI_list_ln_all[[loc]])[1]
    #mtext(paste0(selected_loc[loc],", #nI =",nI),cex=2,side=1)
    dev.off()
  }
  
  #--------------------------Corl plot---------------------------
  
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$Corl
    indI<-posnI_list_ln_all[[i]]
    tempo[[i]][indI]<-NA
  }
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_Corl.pdf",sep=''),width=20, height=20)
    z<-tempo[[loc]]
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=posnI_list_ln_all[[loc]],
               posnN_ind=posnN_list_ln_all[[loc]],
               colrange=c(-cr,cr),nsm=r[[loc]]$neg_sp_mat)
    #nI<-dim(posnI_list_ln_all[[loc]])[1]
    #mtext(paste0(selected_loc[loc],", #nI =",nI),cex=2,side=1)
    dev.off()
  }
  
  #--------------------------Coru plot---------------------------
  
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$Coru
    indI<-posnI_list_ln_all[[i]]
    tempo[[i]][indI]<-NA
  }
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_Coru.pdf",sep=''),width=20, height=20)
    z<-tempo[[loc]]
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=posnI_list_ln_all[[loc]],
               posnN_ind=posnN_list_ln_all[[loc]],
               colrange=c(-cr,cr),nsm=r[[loc]]$neg_sp_mat)
    #nI<-dim(posnI_list_ln_all[[loc]])[1]
    #mtext(paste0(selected_loc[loc],", #nI =",nI),cex=2,side=1)
    dev.off()
  }
  
  #--------------------------Corl-Coru plot---------------------------
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$Corl-data_ln_all[[i]]$Coru
    indI<-posnI_list_ln_all[[i]]
    indN<-posnN_list_ln_all[[i]]
    tempo[[i]][indI]<-NA
    #tempo[[i]][indN]<-abs(tempo[[i]][indN])
    diag(tempo[[i]])<-NA
  }
  
  CorlmCoru<-tempo
  names(CorlmCoru)<-selected_loc
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  summary_nLU_CorlmCoru<-matrix(NA,2,numloc) # to keep total count on CorlmCoru for +ve or -ve numbers
  colnames(summary_nLU_CorlmCoru)<-selected_loc
  rownames(summary_nLU_CorlmCoru)<-c("nL","nU")
  
  summary_LU_CorlmCoru<- summary_nLU_CorlmCoru  # to keep sum on CorlmCoru values only for +ve or -ve numbers
  rownames(summary_LU_CorlmCoru)<-c("L","U")
  
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_Corl-Coru.pdf",sep=''),width=20, height=20)
    z<-tempo[[loc]]
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=posnI_list_ln_all[[loc]],
               posnN_ind=posnN_list_ln_all[[loc]],
               colrange=c(-cr,cr),nsm=r[[loc]]$neg_sp_mat)
    #nI<-dim(posnI_list_ln_all[[loc]])[1]
    z[posnN_list_ln_all[[loc]]]<-NA
    dl2<-nrow(z)-nvar
    z1<-z[1:dl2,1:dl2]
    nL<-sum(z1>0,na.rm = T)
    nU<-sum(z1<0,na.rm = T)
    L<-sum(z1[which(z1>0,arr.ind=T)])
    U<-sum(z1[which(z1<0,arr.ind=T)])
    summary_nLU_CorlmCoru[1,loc]<-nL
    summary_nLU_CorlmCoru[2,loc]<-nU
    summary_LU_CorlmCoru[1,loc]<-L
    summary_LU_CorlmCoru[2,loc]<-U
    if(tagon == T){
      mtext(paste0(selected_loc[loc],"  "),cex=5,side=1,col="red",adj=0.3)
    }
    mtext(paste0("nL =",nL,", nU =",nU),cex=5,side=1,adj=0.7)
    dev.off()
  }
  
  #--------------------------Pl plot---------------------------
  
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$Pl
    indI<-posnI_list_ln_all[[i]]
    tempo[[i]][indI]<-NA
  }
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_Pl.pdf",sep=''),width=20, height=20)
    z<-tempo[[loc]]
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=posnI_list_ln_all[[loc]],
               posnN_ind=posnN_list_ln_all[[loc]],
               colrange=c(-cr,cr),nsm=r[[loc]]$neg_sp_mat)
    #nI<-dim(posnI_list_ln_all[[loc]])[1]
    #mtext(paste0(selected_loc[loc],", #nI =",nI),cex=2,side=1)
    dev.off()
  }
  
  
  #--------------------------Pu plot---------------------------
  
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$Pu
    indI<-posnI_list_ln_all[[i]]
    tempo[[i]][indI]<-NA
  }
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_Pu.pdf",sep=''),width=20, height=20)
    z<-tempo[[loc]]
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=posnI_list_ln_all[[loc]],
               posnN_ind=posnN_list_ln_all[[loc]],
               colrange=c(-cr,cr),nsm=r[[loc]]$neg_sp_mat)
    #nI<-dim(posnI_list_ln_all[[loc]])[1]
    #mtext(paste0(selected_loc[loc],", #nI =",nI),cex=2,side=1)
    dev.off()
  }
  
  
  #--------------------------Pl-Pu plot---------------------------
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$Pl-data_ln_all[[i]]$Pu
    indI<-posnI_list_ln_all[[i]]
    indN<-posnN_list_ln_all[[i]]
    tempo[[i]][indI]<-NA
    #tempo[[i]][indN]<-abs(tempo[[i]][indN])
    diag(tempo[[i]])<-NA
  }
  
  PlmPu<-tempo
  names(PlmPu)<-selected_loc
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  summary_nLU_PlmPu<-matrix(NA,2,numloc)
  colnames(summary_nLU_PlmPu)<-selected_loc
  rownames(summary_nLU_PlmPu)<-c("nL","nU")
  
  summary_LU_PlmPu<- summary_nLU_PlmPu  # to keep sum on PlmPu values only for +ve or -ve numbers
  rownames(summary_LU_PlmPu)<-c("L","U")
  
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_Pl-Pu.pdf",sep=''),width=20, height=20)
    z<-tempo[[loc]]
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=posnI_list_ln_all[[loc]],
               posnN_ind=posnN_list_ln_all[[loc]],
               colrange=c(-cr,cr),nsm=r[[loc]]$neg_sp_mat)
    #nI<-dim(posnI_list_ln_all[[loc]])[1]
    z[posnN_list_ln_all[[loc]]]<-NA
    dl2<-nrow(z)-nvar
    z1<-z[1:dl2,1:dl2]
    nL<-sum(z1>0,na.rm = T)
    nU<-sum(z1<0,na.rm = T)
    L<-sum(z1[which(z1>0,arr.ind=T)])
    U<-sum(z1[which(z1<0,arr.ind=T)])
    summary_nLU_PlmPu[1,loc]<-nL
    summary_nLU_PlmPu[2,loc]<-nU
    summary_LU_PlmPu[1,loc]<-L
    summary_LU_PlmPu[2,loc]<-U
    if(tagon == T){
     mtext(paste0(selected_loc[loc],"  "),cex=5,side=1,col="red",adj=0.3)
    }
    mtext(paste0("nL =",nL,", nU =",nU),cex=5,side=1,adj=0.7)
    dev.off()
  }
  
  
  #--------------------------D2l plot---------------------------
  
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$D2l
    indI<-posnI_list_ln_all[[i]]
    tempo[[i]][indI]<-NA
  }
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_D2l.pdf",sep=''),width=20, height=20)
    z<-tempo[[loc]]
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=posnI_list_ln_all[[loc]],
               posnN_ind=posnN_list_ln_all[[loc]],
               colrange=c(-cr,cr),nsm=r[[loc]]$neg_sp_mat)
    #nI<-dim(posnI_list_ln_all[[loc]])[1]
    #mtext(paste0(selected_loc[loc],", #nI =",nI),cex=2,side=1)
    dev.off()
  }
  
  #--------------------------D2u plot---------------------------
  
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$D2u
    indI<-posnI_list_ln_all[[i]]
    tempo[[i]][indI]<-NA
  }
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_D2u.pdf",sep=''),width=20, height=20)
    z<-tempo[[loc]]
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=posnI_list_ln_all[[loc]],
               posnN_ind=posnN_list_ln_all[[loc]],
               colrange=c(-cr,cr),nsm=r[[loc]]$neg_sp_mat)
    #nI<-dim(posnI_list_ln_all[[loc]])[1]
    #mtext(paste0(selected_loc[loc],", #nI =",nI),cex=2,side=1)
    dev.off()
  }
  
  #--------------------------D2u-D2l plot---------------------------
  tempo<-vector("list",numloc)
  for(i in 1:numloc){
    tempo[[i]]<-data_ln_all[[i]]$D2u-data_ln_all[[i]]$D2l
    indI<-posnI_list_ln_all[[i]]
    indN<-posnN_list_ln_all[[i]]
    tempo[[i]][indI]<-NA
    #tempo[[i]][indN]<-abs(tempo[[i]][indN])
    diag(tempo[[i]])<-NA
  }
  
  D2umD2l<-tempo
  names(D2umD2l)<-selected_loc
  
  minval<-min(unlist(lapply(tempo,min,na.rm=T)))
  maxval<-max(unlist(lapply(tempo,max,na.rm=T)))
  
  cr<-max(abs(minval),abs(maxval))
  
  summary_nLU_D2umD2l<-matrix(NA,2,numloc)
  colnames(summary_nLU_D2umD2l)<-selected_loc
  rownames(summary_nLU_D2umD2l)<-c("nL","nU")
  
  summary_LU_D2umD2l<- summary_nLU_D2umD2l  # to keep sum on D2umD2l values only for +ve or -ve numbers
  rownames(summary_LU_D2umD2l)<-c("L","U")
  
  for(loc in 1:numloc){
    resloc2<-paste(resloc,selected_loc[loc],"/",sep="")
    pdf(paste(resloc2,selected_loc[loc],file="_D2u-D2l.pdf",sep=''),width=20, height=20)
    z<-tempo[[loc]]
    if(nvar!=0){
      dl<-nrow(z)-nvar+1
      rownames(z)[c(dl:nrow(z))]<-nvar_names
    }
    colnames(z)<-rownames(z)
    mycorrplot(z=z,
               posnI_ind=posnI_list_ln_all[[loc]],
               posnN_ind=posnN_list_ln_all[[loc]],
               colrange=c(-cr,cr),nsm=r[[loc]]$neg_sp_mat)
    #nI<-dim(posnI_list_ln_all[[loc]])[1]
    z[posnN_list_ln_all[[loc]]]<-NA
    dl2<-nrow(z)-nvar
    z1<-z[1:dl2,1:dl2]
    nL<-sum(z1>0,na.rm = T)
    nU<-sum(z1<0,na.rm = T)
    L<-sum(z1[which(z1>0,arr.ind=T)])
    U<-sum(z1[which(z1<0,arr.ind=T)])
    summary_nLU_D2umD2l[1,loc]<-nL
    summary_nLU_D2umD2l[2,loc]<-nU
    summary_LU_D2umD2l[1,loc]<-L
    summary_LU_D2umD2l[2,loc]<-U
    if(tagon == T){
     mtext(paste0(selected_loc[loc],"  "),cex=5,side=1,col="red",adj=0.3)
    }
     mtext(paste0("nL =",nL,", nU =",nU),cex=5,side=1,adj=0.7)
    dev.off()
  }
  
  return(list(CorlmCoru_all_ln_list=CorlmCoru,
              PlmPu_all_ln_list=PlmPu,
              D2umD2l_all_ln_list=D2umD2l,
              summary_nLU_CorlmCoru=summary_nLU_CorlmCoru,
              summary_nLU_PlmPu=summary_nLU_PlmPu,
              summary_nLU_D2umD2l=summary_nLU_D2umD2l,
              summary_LU_CorlmCoru=summary_LU_CorlmCoru,
              summary_LU_PlmPu=summary_LU_PlmPu,
              summary_LU_D2umD2l=summary_LU_D2umD2l))
}

