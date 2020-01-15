# Select any two species pair [i,j] for a copula of location loc
# This function gives you a matrix with vi and vj as two columns

# Input :
# d_allsp : dataset in data[[loc]][[sp]] format
# loc : location index
# i,j : sp-pair indices
# level : significance level for BiCopIndepTest p-value
# ploton : (optional) logical, if T gives copula plot without transforming j-th variable to it's -ve value
# onbounds : a logical tag (default=FALSE) to get info about the points exactly lying on bounds, if set to TRUE 
#            then the arguments lb and ub must be numeric
# lb : numeric value between [0,1] for lower bound (default =NA)
# ub : numeric value between [0,1] for upper bound (default =NA), ub should be greater than lb 
# include_indep: logical (whether to keep track for indep. test or not)

# Output :
# A list of 4 elements:
#                      mat : a matrix : copula of (vi,vj) with transforming j-th variable to it's -ve value for -ve corr.
#                      corval : Spearman's correlation
#                      pval   : pvalue of Spearman's cor.test
#                      IndepTestRes : BiCopIndepTest p-value

# and an optional plot of the copula

library(VineCopula)
vivj_matrix<-function(d_allsp,loc,i,j,level=0.05,ploton,onbounds=F,lb=NA,ub=NA,include_indep){
  
  ds1<-d_allsp[[loc]][[i]]
  ds2<-d_allsp[[loc]][[j]]
  #----------------------------
  colnames(ds1)<-c("Year","Dat")  # ensuring column names
  colnames(ds2)<-c("Year","Dat")
  
  a1<-ds1$Year[1]
  a2<-ds2$Year[1]
  a3<-ds1$Year[dim(ds1)[1]]
  a4<-ds2$Year[dim(ds2)[1]]
  year_s<-max(a1,a2)
  year_e<-min(a3,a4)
  ind_s1<-which(ds1$Year==year_s)
  ind_s2<-which(ds2$Year==year_s)
  ind_e1<-which(ds1$Year==year_e)
  ind_e2<-which(ds2$Year==year_e)
  ds1<-ds1[ind_s1:ind_e1,]
  ds2<-ds2[ind_s2:ind_e2,]
  # Omitting the years and data containing NA in either d1 or d2 
  #from both d1 and d2
  if(anyNA(ds1$Dat)==T | anyNA(ds2$Dat)==T){
    ind_na1<-which(is.na(ds1$Dat))
    ind_na2<-c(ind_na1,which(is.na(ds2$Dat)))
    ind_na<-unique(ind_na2)
    
    d1Dat<-ds1$Dat[-ind_na]
    d2Dat<-ds2$Dat[-ind_na]
    Years<-ds1$Year[-ind_na]
    d1<-data.frame(Year=Years,Dat=d1Dat)
    d2<-data.frame(Year=Years,Dat=d2Dat)
  }else{
    d1<-ds1
    d2<-ds2
  }
  
  colnames(d1)[2]<-"Dat"  # ensuring column names
  colnames(d2)[2]<-"Dat"
  
  #get ranks modified now
  vi<-VineCopula::pobs(d1$Dat)
  vj<-VineCopula::pobs(d2$Dat)
  
  IndepTestRes<-VineCopula::BiCopIndTest(vi,vj)$p.value
  ct<-cor.test(vi,vj,alternative = "two.sided",method="spearman",exact=F)
  corval<-unname(ct$estimate)
  pval<-ct$p.value
  
  if(ploton==T){
    if(include_indep==T){
      if(IndepTestRes<level && corval>0){ # for significant positive correlation
        plot(vi,vj,type='p',col=rgb(0,0,0,0.3),pch=19,xlim=c(0,1),ylim=c(0,1),
             xlab=names(d_allsp[[loc]])[i],ylab=names(d_allsp[[loc]])[j],cex.lab=1.5)
        
        
        if(j>i){
          if(onbounds==T & identical(vi,vj)==F){
            ind_lb<-which(vi+vj==(2*lb))
            ind_ub<-which(vi+vj==(2*ub))
            onlb<-length(ind_lb)
            onub<-length(ind_ub)
            
            if(onlb!=0 | onub!=0){
              mtext(paste0("onbs = (",onlb," , ",onub,")"),
                    side = 4, line=0.15, adj=0.5, col="red") 
            }
          }
        }
        
        
      }else if(IndepTestRes<level && corval<0){ # for significant negative correlation
        plot(vi,vj,type='p',col=rgb(0,1,0,0.3),pch=19,xlim=c(0,1),ylim=c(0,1),
             xlab=names(d_allsp[[loc]])[i],ylab=names(d_allsp[[loc]])[j],cex.lab=1.5)
        
        if(j>i){
          if(onbounds==T & identical(vi,vj)==F){
            vneg<-VineCopula::pobs(-(d2$Dat)) # see when we count points on bounds we took reverse of second variable
            ind_lb<-which(vi+vneg==(2*lb))
            ind_ub<-which(vi+vneg==(2*ub))
            
            #vneg<-VineCopula::pobs(-(d1$Dat)) # NOTE : onbs will not be same if we consider first variable to be reversed
            #ind_lb<-which(vj+vneg==(2*lb))
            #ind_ub<-which(vj+vneg==(2*ub))
            
            onlb<-length(ind_lb)
            onub<-length(ind_ub)
            
            if(onlb!=0 | onub!=0){ 
              mtext(paste0("onbs = (",onlb," , ",onub,")"),
                    side = 4, line=0.15, adj=0.5, col="red") 
            }
          }
        }
        
      }else{ # independent case
        plot(-1,0,xlim=c(0,1),ylim=c(0,1),xlab=names(d_allsp[[loc]])[i],ylab=names(d_allsp[[loc]])[j],cex.lab=1.5)
        text(0.5,0.5,"Indep.",adj=c(0.5,.5),cex=2)
      }
      mtext(paste0("(sp_x, sp_y) = (",i," , ",j,")"),
            side = 3, line=0.15, adj=0.5, col="black")
      
      if(IndepTestRes<level && corval<0){
        vj<-VineCopula::pobs(-(d2$Dat))
      }
      
    }else{ # when we consider all cells however week their correlations were.
      if(corval>0){ # for all +ve correlations
        plot(vi,vj,type='p',col=rgb(0,0,0,0.3),pch=19,xlim=c(0,1),ylim=c(0,1),
             xlab=names(d_allsp[[loc]])[i],ylab=names(d_allsp[[loc]])[j],cex.lab=1.5)
        
        
        if(j>i){
          if(onbounds==T & identical(vi,vj)==F){
            ind_lb<-which(vi+vj==(2*lb))
            ind_ub<-which(vi+vj==(2*ub))
            onlb<-length(ind_lb)
            onub<-length(ind_ub)
            
            if(onlb!=0 | onub!=0){
              mtext(paste0("onbs = (",onlb," , ",onub,")"),
                    side = 4, line=0.15, adj=0.5, col="red") 
            }
          }
        }
      }else{ # for all -ve correlations
        plot(vi,vj,type='p',col=rgb(0,1,0,0.3),pch=19,xlim=c(0,1),ylim=c(0,1),
             xlab=names(d_allsp[[loc]])[i],ylab=names(d_allsp[[loc]])[j],cex.lab=1.5)
        
        if(j>i){
          if(onbounds==T & identical(vi,vj)==F){
            vneg<-VineCopula::pobs(-(d2$Dat)) # see when we count points on bounds we took reverse of second variable
            ind_lb<-which(vi+vneg==(2*lb))
            ind_ub<-which(vi+vneg==(2*ub))
            
            #vneg<-VineCopula::pobs(-(d1$Dat)) # NOTE : onbs will not be same if we consider first variable to be reversed
            #ind_lb<-which(vj+vneg==(2*lb))
            #ind_ub<-which(vj+vneg==(2*ub))
            
            onlb<-length(ind_lb)
            onub<-length(ind_ub)
            
            if(onlb!=0 | onub!=0){ 
              mtext(paste0("onbs = (",onlb," , ",onub,")"),
                    side = 4, line=0.15, adj=0.5, col="red") 
            }
          }
        }
      }
      
      if(corval<0){ #reverse the variable
        vj<-VineCopula::pobs(-(d2$Dat))
      }
      
    }
    
  }
  
  Years<-d1$Year
  #-------------------------
  #n_datapt<-length(vi)
  #--------------------
  #plot(vi,vj,type="p")
  #-------------------------
  mat<-as.matrix(cbind(vi,vj))
  return(list(mat=mat,   # return reversed mat so that if you plot this mat you get +ve correlation 
              corval=corval,  # but return the actual -ve corr. value  
              pval=pval,
              IndepTestRes=IndepTestRes))  
}



