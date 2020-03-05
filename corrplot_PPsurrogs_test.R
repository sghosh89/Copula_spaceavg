# test the code with Hays data

m<-readRDS("./Results/hays_results/skewness_results/ts_mat_CP_hays.RDS")
cm<-cor(m)

ans<-readRDS("./Results/hays_results/skewness_results/pp_surrogs_hays_CP/PPsurrogs_tests_with_HaysSurrogates.RDS")
s1<-ans$cor_surrogs

cs<-array(numeric(),c(nrow(cm),ncol(cm),0))

for(i in c(1:nrow(cm))){
  temp<-cor(s1[,,i])
  cs<-abind::abind(cs,temp) # array with correlation for surrogate values
}

dim(cs)

getqs<-function(x){
  qs<-quantile(x,c(0.025,0.25,0.5,0.75,0.975))
  ql<-unname(qs[1]) 
  qh<-unname(qs[5]) 
  return(c(ql,qh))
}

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


library(corrplot)
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
col1 <- colorRampPalette(c("blue","white","red")) 
#corrplot(cm, type="upper",diag=F,col=col4(10))
corrplot(cm, type="upper",diag=F)

pdf("./corrplot.pdf",width=19,height=19)
op<-par(mar=c(0.5,0.5,0.5,0.5))
rownames(cm)<-colnames(cm)<-paste("sp",c(1:19))
diag(cm)<-diag(cm_low)<-diag(cm_high)<-NA
corrplot(cm, low = cm_low, upp = cm_high, type='lower',diag=F, plotC = "rect", col=col1(100), 
         tl.cex=2,tl.col = "black",tl.offset = 2.4,tl.pos="ld",addgrid.col = "black",
         cl.cex = 2,cl.length=7,tl.srt=30,
         cl.align.text = "l",cl.ratio = 0.1)
par(op)
dev.off()





