source("./make_tab_stability_assessment.R")

set.seed(seed=104)
library(copula)

#A function for generating data for a pedagogical figure, copula scale
#
#Args
#n        The number of random vectors you want
#sig      A covariance matrix. See the code and examples for how this is used
#
#Output - A matrix of dimensions n by 2*(dim(sig)[1]-1), draws from the random 
#variable corresponding to a 2*(dim(sig)[1]-1)-dimensional copula
#
getdat<-function(n,sig){
  coinflip<-runif(n)
  cop<-copula::normalCopula(param=P2p(sig),dim=dim(sig)[1],dispstr="un")
  intres<-rCopula(n,cop)
  d<-dim(sig)[1]-1
  
  #this will be the result to fill in
  res<-matrix(NA,n,2*d)
  
  #for about half the cases, the first d variables are in the interval (1/2,1)
  #and are the same, and the rest of the variables are in the interval (0,1/2) 
  res[coinflip<1/2,1:d]<-rep(intres[coinflip<1/2,1]/2+1/2,times=d)
  res[coinflip<1/2,(d+1):(2*d)]<-intres[coinflip<1/2,2:(d+1)]/2
  
  #for the other half of the cases, the last d variables are in the interval (1/2,1)
  #and are the same, and the first d variables are in the interval (0,1/2)
  res[coinflip>=1/2,1:d]<-intres[coinflip>=1/2,2:(d+1)]/2
  res[coinflip>=1/2,(d+1):(2*d)]<-rep(intres[coinflip>=1/2,1]/2+1/2,times=d)
  
  return(res)
}


#-------------- Now generate the data, copula scale, for the first half of the figure

n<-10000
d<-5

syncconst<-0.1 
compconst<-0.2
sig_p<-matrix(syncconst,d,d)
diag(sig_p)<-1
sig<-cbind(rep(-compconst,d),sig_p)
sig<-rbind(c(1,rep(-compconst,d)),sig)
sig
eigen(sig,symmetric=TRUE,only.values=TRUE)

res_p1<-getdat(n=10000,sig=sig)

#----------- Now generate the data, copula scale, for the second half of the figure
res_p2<-getdat(n=10000,sig=sig)
res_p2<-1-res_p2

#------------- Now combine both copulas with standard normal marginals

normres_p1<-qnorm(res_p1) # qnorm makes marginals normal from uniform
normres_p2<-qnorm(res_p2)

normres_p1<-normres_p1-min(normres_p1,normres_p2)+1 # we give this normal marginals a shift towards +ve axis 
                                                    # as they are sp. biomass timeseries
normres_p2<-normres_p2-min(normres_p2,normres_p2)+1

#------------------ Stability metrics ------------------
# extremely different skewness ratios though have nearly same variance ratios
(tab1<-make_tab_stability(m=normres_p1,surrogs = NA,surrogs_given = F))
(tab2<-make_tab_stability(m=normres_p2,surrogs = NA,surrogs_given = F))

#--------------- Now plot -----------------------------------------
tim<-c(1:60)
normres_p1<-normres_p1[tim,]
normres_p2<-normres_p2[tim,]
ylm<-max(normres_p1,normres_p2)+4

pdf("./Results/pedagog_figs/pedagog_cv2_skw.pdf",height=3,width=9)
op<-par(mar=c(0.5, 0.5, 0.5, 4), mfcol=c(2,2), oma = c(3.5,4.5, 0.1, 0.1))

# first column for individual biomass of each sp.
plot.ts(normres_p2,col=c(rep(1,d),rep(2,d)),ylab="",xlab="",xaxt="n",plot.type = "single",ylim=c(1,ylm))

legend("topleft","(A) species are synchronously rare \n than abundant within each group of the community",
                      bty="n",cex=0.8)

#get variance ratio
phi_cv_l<-round(tab2$phi_cvsq,2)
legend("topright",as.expression(bquote(phi[cv]==.(phi_cv_l))),bty="n",cex=1,text.col="blue")

plot.ts(normres_p1,col=c(rep(1,d),rep(2,d)),ylab="",xlab="",plot.type = "single",ylim=c(1,ylm))
legend("topleft","(B) species are synchronously abundant \n than rare within each group of the community",
       bty="n",cex=0.8)

#get variance ratio
phi_cv_r<-round(tab1$phi_cvsq,2)
legend("topright",as.expression(bquote(phi[cv]==.(phi_cv_r))),bty="n",cex=1,text.col="blue")

mtext(paste("Individual sp. (i = ",ncol(normres_p1),") biomass \n in the community",sep=""),side=2,line=2.2, adj = 0,cex=0.8)
mtext("Year",side=1,line=2.2, adj = 0.5,cex=0.8)


# second column for aggregated biomass
tot_ts_l<-apply(FUN=sum,MARGIN=1,X=normres_p2)
tot_ts_r<-apply(FUN=sum,MARGIN=1,X=normres_p1)
ylm<-range(tot_ts_l,tot_ts_r)

plot(tim,tot_ts_l,type="l",ylim=c(ylm[1],ylm[2]+5),col="darkgray",ylab="",xaxt="n",xlab="")
legend("topleft",paste("(C) s = ",round(tab2$skw_real,2)),bty="n",cex=0.8)


plot(tim,tot_ts_r,type="l",ylim=c(ylm[1],ylm[2]+5),col="darkgray",xlab="",ylab="")
legend("topleft",paste("(D) s = ",round(tab1$skw_real,2)),bty="n",cex=0.8)

mtext("Total community biomass",side=2,line=2.2, adj = -1.0,cex=0.8)
mtext("Year",side=1,line=2.2, adj = 0.5,cex=0.8)

par(op)
#box(which="figure")
dev.off()







