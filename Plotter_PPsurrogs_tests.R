# Plot to check the pairwise Pearson correlation from real data and surrogs are close enough or not?

library(corrplot)

getqs<-function(x){
  qs<-quantile(x,c(0.025,0.25,0.5,0.75,0.975))
  ql<-unname(qs[2]) #25% quantile
  qh<-unname(qs[4]) #75% quantile
  return(c(ql,qh))
}

Plotter_PPsurrogs_tests<-function(m,ans,resloc,wd=15,ht=15){
  
  cm<-cor(m)
  s1<-ans$cor_surrogs
  
  cs<-array(numeric(),c(nrow(cm),ncol(cm),0)) #initialize
  
  for(i in c(1:nrow(cm))){
    temp<-cor(s1[,,i])
    cs<-abind::abind(cs,temp) # array with correlation for surrogate values
  }
  
  #dim(cs)
  
  cm_low<-cm
  cm_high<-cm
  
  for(i in 1:nrow(cs)){
    for(j in 1:ncol(cs)){
      mm<-cs[i,j,]
      tempo<-getqs(mm)
      cm_low[i,j]<-tempo[1]
      cm_high[i,j]<-tempo[2]
    }
  }
  
  #col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
  #                           "cyan", "#007FFF", "blue", "#00007F"))
  #corrplot(cm, type="upper",diag=F,col=col4(10))
  #corrplot(cm, type="upper",diag=F)
  
  col1 <- colorRampPalette(c("blue","white","red")) 
  
  pdf(paste(resloc,"./compare_pairwisecor_PPsurrogs.pdf",sep=""),width=wd,height=ht)
  op<-par(mar=c(0.5,0.5,0.5,0.5))
  rownames(cm)<-colnames(cm)<-paste("sp",c(1:nrow(cm)))
  diag(cm)<-diag(cm_low)<-diag(cm_high)<-NA
  corrplot(cm, low = cm_low, upp = cm_high, type='lower',diag=F, plotC = "rect", col=col1(100), 
           tl.cex=2,tl.col = "black",tl.offset = 1,tl.pos="ld",addgrid.col = "grey",
           cl.cex = 2,cl.length=7,tl.srt=20,
           cl.align.text = "l",cl.ratio = 0.1)
  par(op)
  dev.off()
  
  
  pdf(paste(resloc,"./comparison_histplot_PPsurrogs.pdf",sep=""),height=2.5,width=18)
  op<-par(mfrow=c(1,3),mar=c(5.5,5.5,0.8,0.5),mgp=c(4,0.5,0))
  
  xlm<-sort(c(range(ans$tot_cov_surrogs),ans$tot_cov_real,0))
  hist(ans$tot_cov_surrogs,col="grey",border=F,breaks=1000,xaxt="n",cex.axis=1.1,
       xlab="Total pairwise-covariance \n from surrogates",main="",cex.lab=1.5)
  axis(side=1, at=round(xlm,2),cex.axis=1.1)
  abline(v=ans$tot_cov_real,col="black") # actual skw from real data
  
  
  xlm<-sort(c(range(ans$var_aggr_sp_surrogs),ans$var_aggr_sp_real))
  hist(ans$var_aggr_sp_surrogs,col="grey",border=F,breaks=1000,cex.axis=1.1,xaxt="n",
       xlab="Variance of aggregated species time-series \n from surrogates",main="",cex.lab=1.5)
  axis(side=1, at=xlm,cex.axis=1.1)
  abline(v=ans$var_aggr_sp_real,col="black") # actual skw from real data
  
  xlm<-sort(c(range(ans$vr_surrogs),ans$vr_real,1))
  hist(ans$vr_surrogs,col="grey",border=F,breaks=1000,cex.axis=1.1,xaxt="n",
       xlab="Variance ratio \n from surrogates",main="",cex.lab=1.5)
  axis(side=1, at=round(xlm,2),cex.axis=1.1)
  abline(v=ans$vr_real,col="black") # actual skw from real data
  
  par(op)
  dev.off()
  
}


# call the function for Hays

#m<-readRDS("./Results/hays_results/skewness_results/ts_mat_CP_hays.RDS")
#ans<-readRDS("./Results/hays_results/skewness_results/pp_surrogs_hays_CP/PPsurrogs_tests_with_HaysSurrogates.RDS")
#resloc<-"./Results/hays_results/skewness_results/pp_surrogs_hays_CP/"

#Plotter_PPsurrogs_tests(m=m,ans=ans,resloc=resloc,ht=15,wd=15)


# call the function for KNZ

#m<-readRDS("./Results/knz_results/skewness_results/ts_CP_knz_soiltype_t.RDS")
#ans<-readRDS("./Results/knz_results/skewness_results/pp_surrogs_knz_t_CP/PPsurrogs_tests_with_KNZtSurrogates.RDS")
#resloc<-"./Results/knz_results/skewness_results/pp_surrogs_knz_t_CP/"

#Plotter_PPsurrogs_tests(m=m,ans=ans,resloc=resloc,ht=18,wd=18)













