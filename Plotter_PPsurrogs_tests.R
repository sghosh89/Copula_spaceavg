# Plot to check the pairwise Pearson correlation from real data and surrogs are close enough or not?

library(corrplot)

getqs<-function(x){
  qs<-quantile(x,c(0.025,0.25,0.5,0.75,0.975))
  ql<-unname(qs[1]) #2.5% quantile
  qm<-unname(qs[3]) #50% quantile
  qh<-unname(qs[5]) #97.5% quantile
  return(c(ql,qm,qh))
}


# Args:
# m: a matrix with common and pseudo sp time series along each column
# ans: a list which is the output from PPsurrogs_tests function
# resloc: folder path to save the plots
# tag_legend: a character string of length 3 for comparison hist plot

Plotter_PPsurrogs_tests<-function(m,ans,resloc,tag_legend){
  
  cm<-cor(m)
  s1<-ans$cor_surrogs
  
  cs<-array(numeric(),c(nrow(cm),ncol(cm),0)) #initialize
  
  for(i in c(1:nrow(cm))){
    temp<-cor(s1[,,i])
    cs<-abind::abind(cs,temp) # array with correlation for surrogate values
  }
  
  #dim(cs)
  
  cm_low<-cm_mid<-cm_up<-cm # initialize
  
  for(i in 1:nrow(cs)){
    for(j in 1:ncol(cs)){
      mm<-cs[i,j,]
      tempo<-getqs(mm)
      cm_low[i,j]<-tempo[1]
      cm_mid[i,j]<-tempo[2]
      cm_up[i,j]<-tempo[3]
    }
  }
 
  nspair<-0.5*nrow(cm)*(nrow(cm)-1) # number of unique species pair
  my_dat<-data.frame(real_cor=NA*numeric(nspair),
                     surrogs_cor_lowCI=NA*numeric(nspair),
                     surrogs_cor_midCI=NA*numeric(nspair),
                     surrogs_cor_upCI=NA*numeric(nspair),
                     color_code=NA*numeric(nspair)) 
  
  my_dat$real_cor<-cm[lower.tri(cm)]
  my_dat$surrogs_cor_lowCI<-cm_low[lower.tri(cm_low)]
  my_dat$surrogs_cor_midCI<-cm_mid[lower.tri(cm_mid)]
  my_dat$surrogs_cor_upCI<-cm_up[lower.tri(cm_up)]
  
  id<-which(my_dat$surrogs_cor_lowCI>my_dat$real_cor | my_dat$surrogs_cor_upCI<my_dat$real_cor)
  
  my_dat$color_code<-rgb(0,0,0,0.3) #black
  my_dat$color_code[id]<-rgb(1,0,0,0.6) #red
 
  pdf(paste(resloc,"./comparison_pairwisecor_PPsurrogs.pdf",sep=""),height=5,width=5)
  op<-par(mar=c(5,5,1,1),mgp=c(3,1,0))
  
  plot(my_dat$real_cor,my_dat$surrogs_cor_midCI,
       xlim=c(-1,1),ylim=c(-1,1),
       xlab="pairwise Pearson correlation from real data",ylab="pairwise Pearson correlation from surrogates",
       cex=0.5,pch=20,col=my_dat$color_code)
  segments(x0=-1,y0=-1,x1=1,y1=1,lty="dashed")
  
  segments(x0=my_dat$real_cor,y0=my_dat$surrogs_cor_lowCI,
           x1=my_dat$real_cor,y1=my_dat$surrogs_cor_upCI,col=my_dat$color_code)
  
  par(op)
  dev.off()
  #-----------------------------------------------------------------------------------------------
  #col1 <- colorRampPalette(c("blue","white","red")) 
  #pdf(paste(resloc,"./compare_pairwisecor_PPsurrogs.pdf",sep=""),width=15,height=15)
  #op<-par(mar=c(0.5,0.5,0.5,0.5))
  #rownames(cm)<-colnames(cm)<-paste("sp",c(1:nrow(cm)))
  #diag(cm)<-diag(cm_low)<-diag(cm_high)<-NA
  #corrplot(cm, low = cm_low, upp = cm_high, type='lower',diag=F, plotC = "rect", col=col1(100), 
  #         tl.cex=2,tl.col = "black",tl.offset = 1,tl.pos="ld",addgrid.col = "grey",
  #         cl.cex = 2,cl.length=7,tl.srt=20,
  #         cl.align.text = "l",cl.ratio = 0.1)
  #par(op)
  #dev.off()
  #-----------------------------------------------------------------------------------------------
  
  
  pdf(paste(resloc,"./comparison_histplot_PPsurrogs.pdf",sep=""),height=7,width=9)
  op<-par(mfrow=c(3,1),mar=c(9,9,3,4),mgp=c(5,1,0))
  
  xlm<-sort(c(range(ans$tot_cov_surrogs),ans$tot_cov_real,0))
  pval<-sum(ans$tot_cov_surrogs<ans$tot_cov_real)/length(ans$tot_cov_surrogs)
  hist(ans$tot_cov_surrogs,col="grey",border=F,breaks=1000,xaxt="n",cex.axis=1.6,cex.main=2,
       xlab="Total pairwise-covariance from surrogates",cex.lab=2.4,main=paste("p = ",pval))
  legend("topleft", tag_legend[1], bty="n", cex=3)
  axis(side=1, at=round(xlm,2),cex.axis=1.6)
  abline(v=ans$tot_cov_real,col="black") # actual skw from real data
  
  
  xlm<-sort(c(range(ans$var_aggr_sp_surrogs),ans$var_aggr_sp_real))
  pval<-sum(ans$var_aggr_sp_surrogs<ans$var_aggr_sp_real)/length(ans$var_aggr_sp_surrogs)
  hist(ans$var_aggr_sp_surrogs,col="grey",border=F,breaks=1000,xaxt="n",cex.axis=1.6,cex.main=2,
       xlab="Variance of aggregated species time-series from surrogates",cex.lab=2.4,main=paste("p = ",pval))
  legend("topleft", tag_legend[2], bty="n", cex=3)
  axis(side=1, at=round(xlm,2),cex.axis=1.6)
  abline(v=ans$var_aggr_sp_real,col="black") # actual skw from real data
  
  xlm<-sort(c(range(ans$vr_surrogs),ans$vr_real,1))
  pval<-sum(ans$vr_surrogs<ans$vr_real)/length(ans$vr_surrogs)
  hist(ans$vr_surrogs,col="grey",border=F,breaks=1000,xaxt="n",cex.axis=1.6,cex.main=2,
       xlab="Variance ratio from surrogates",cex.lab=2.4,main=paste("p = ",pval))
  legend("topleft", tag_legend[3], bty="n", cex=3)
  axis(side=1, at=round(xlm,2),cex.axis=1.6)
  abline(v=ans$vr_real,col="black") # actual skw from real data
  
  par(op)
  dev.off()
  
}


# call the function for Hays

#m<-readRDS("./Results/hays_results/skewness_results/ts_mat_CP_hays.RDS")
#ans<-readRDS("./Results/hays_results/skewness_results/pp_surrogs_hays_CP/PPsurrogs_tests_with_HaysSurrogates.RDS")
#resloc<-"./Results/hays_results/skewness_results/pp_surrogs_hays_CP/"
#Plotter_PPsurrogs_tests(m=m,ans=ans,resloc=resloc,tag_legend = c("(A)","(C)","(E)"))


# call the function for KNZ

#m<-readRDS("./Results/knz_results/skewness_results/ts_CP_knz_soiltype_t.RDS")
#ans<-readRDS("./Results/knz_results/skewness_results/pp_surrogs_knz_t_CP/PPsurrogs_tests_with_KNZtSurrogates.RDS")
#resloc<-"./Results/knz_results/skewness_results/pp_surrogs_knz_t_CP/"
#Plotter_PPsurrogs_tests(m=m,ans=ans,resloc=resloc,tag_legend = c("(B)","(D)","(F)"))













