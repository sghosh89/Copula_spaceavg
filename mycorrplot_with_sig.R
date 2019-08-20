# This plotter function is to visualize a matrix
# Input:
#     z : a matrix
#     posnI_ind : a matrix containing the row and col indices of z for which z indicates indep values
#     posnN_ind : a matrix containing the row and col indices of z for which z indicates sig. neg cor values
#     colrange : a vector containing min and max value of the color range
#     type : a character (default value "lower" for lower triangular matrix), other options are "upper" and "full" 

#     sigtest : logical (default value F), if TRUE then generated plot with marker on insignificant cell values based on
#                  tail significance test of Normal null hypothesis, with following additional arguments
#     spr: A matrix of spearman correlation; output from NonParamStat function
#     realstat : tail asymmetry difference matrix, e.g., Corl - Coru 
#     ub : 0.5 as default, see tailsignif function argument
#     numpts : see tailsignif function argument (length of timeseries without NA)
#     numsims : see tailsignif function argument

source("./tailsignif.R")
library(corrplot)

mycorrplot_with_sig<-function(z,posnI_ind,posnN_ind,colrange,type="lower",sigtest=F,spr,realstat,
                     ub=0.5,numpts,numsims,CI=c(0.025,0.975)){
  
  col1 <- colorRampPalette(c("blue","white","red")) 
  
  z[is.na(z)]<-mean(colrange)
  diag(z)[1]<-colrange[1] # just to ensure that plot always have specific colorbar range even 
  diag(z)[2]<-colrange[2]     # though all entries are either +ve or -ve
  
  corrplot(z,is.corr = F,col=col1(100),method="color",addgrid.col = "black",type=type,
           diag=F,bg = "white",tl.cex=2,tl.col = "black",tl.offset = 2.4,tl.pos="ld",
           cl.cex = 2,cl.lim = colrange,mar=c(0,0,0,0),cl.length=7,tl.srt=30,
           cl.align.text = "l",cl.ratio = 0.1)
  
  # colorize as black for diagonal indices
  Dg <- matrix(NA,nrow(z),ncol(z))
  diag(Dg)<- 1 
  
  corrplot(Dg, cl.pos = "n", na.label = " ", add = T,addgrid.col = "black",type=type,
           bg = "transparent", tl.col = "transparent",col="black",method="color")
  
  #colorize as yellow for indep posn indices
  if(dim(posnI_ind)[1]!=0){
    
    I <- matrix(NA,nrow(z),ncol(z))
    I[posnI_ind]<- 1 
    
    corrplot(I, cl.pos = "n", na.label = " ", add = T,addgrid.col = "black",type=type,
             bg = "transparent", tl.col = "transparent",col="yellow",method="color")
  } 
  
  #colorize as green for -ve correlated (siginificantly) posn indices
  if(dim(posnN_ind)[1]!=0){
    
    N <- matrix(NA,nrow(z),ncol(z))
    N[posnN_ind]<- -1 
    
    corrplot(N, cl.pos = "n", na.label = " ", add = T,addgrid.col = "transparent",type=type,
             bg = "transparent", tl.col = "transparent",p.mat = N,sig.level = -2,col="transparent",
             pch=20,pch.col="green",pch.cex = 5,number.cex = 2)
    
  }
  
  if(sigtest==T){
    
    diag(spr)<-NA # omit diagonals
    spr[posnI_ind]<-NA  #omit independence
    spr[posnN_ind]<- -spr[posnN_ind] # making absolute vals to -ve values of spearman correlation as it was for real variables
    id_spr<-which(is.finite(spr),arr.ind=T)
    sprvals<-spr[id_spr]
    realstat[posnN_ind]<-abs(realstat[posnN_ind]) # for negatively correlated cells absolute value of tail diff. matters only
    realstat<-realstat[id_spr]
    
    mat_tab<-cbind(id_spr,sprvals,realstat)
    mat_tab<-as.data.frame(mat_tab)
    
    spcors<-mat_tab$sprvals
    qtl<-tailsignif(ub=ub,numpts=numpts,spcors=spcors,numsims=numsims,CI=CI,resloc=NA)
    mat_tab$lowCI<-qtl[1,]
    mat_tab$upCI<-qtl[2,]
    
    # generate additional plot
    #pdf(paste(sigres,"statistic_vs_spearman_ub_",ub,"_CI_",CI[1],"_",CI[2],".pdf",sep=""),width=8,height=8)
    #plot(c(-1,1),c(0,0),ylim=c(-1,1),xlab="Spearman",ylab="Statistic",type='l',col='red')
    #lines(c(0,0),c(-1,1),type="l",col="red")
    #points(mat_tab$sprvals,mat_tab$realstat,pch=16,col=rgb(1,0,0,0.2))
    #lines(mat_tab$sprvals,mat_tab$lowCI,type='p',pch=16,col=rgb(0,0,0,0.2))
    #lines(mat_tab$sprvals,mat_tab$upCI,type='p',pch=16,col=rgb(0,0,0,0.2))
    #dev.off()
    
    # crossmark the insignificant cells
    insig_idposcor<-which(mat_tab$sprvals>0 & mat_tab$realstat>mat_tab$lowCI & mat_tab$realstat<mat_tab$upCI) # two-tailed : +ve correlation
    insig_idnegcor<-which(mat_tab$sprvals<0 & mat_tab$realstat<mat_tab$upCI) # one-tailed for absolute values : -ve correlation
    
    insig_id<-c(insig_idposcor,insig_idnegcor)
    
    mat_tab$is_sig<-1 
    mat_tab$is_sig[insig_id]<-0 # this cells are not significant
    
    ir<-mat_tab$row[insig_id]
    ic<-mat_tab$col[insig_id]
    
    # creating a significance id matrix filled in with 1(for significance) and 0(for insignificance)
    is_sigmat<-matrix(1,nrow(z),ncol(z))
    diag(is_sigmat)<-NA #omit diagonal
    is_sigmat[posnI_ind]<-NA  #omit independence
    is_sigmat[cbind(ir,ic)]<-0
    
    Isg <- matrix(NA,nrow(z),ncol(z))
    Isg[cbind(ir,ic)]<- -1 
    
    corrplot(Isg, cl.pos = "n", na.label = " ", add = T,addgrid.col = "transparent",type=type,
             bg = "transparent", tl.col = "transparent",p.mat = Isg,sig.level = -2,col="transparent",
             pch=4,pch.col="black",pch.cex = 5,number.cex = 2)
    
    if(type=="lower"){
      
      xl<-lower.tri(is_sigmat)
      id_xl<-which(xl==F,arr.ind=T)
      is_sigmat[id_xl]<-NA
      
    }else if(type=="upper"){
      
      xu<-upper.tri(is_sigmat)
      id_xu<-which(xu==F,arr.ind=T)
      is_sigmat[id_xu]<-NA
      
    }else if(type=="full"){
      
      is_sigmat<-is_sigmat
       
    }else{
      stop("Error in mycorrplot_with_sig.R: arg 'type' must be 'lower' or 'upper' or 'full'")
    }
    
    return(list(mat_tab=mat_tab,
                is_sigmat=is_sigmat))
  }
  
}

#------------------------------------------


