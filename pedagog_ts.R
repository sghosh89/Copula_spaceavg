library(copula)
library(e1071)
source("./ncsurrog.R")
source("./ExtremeTailDep.R")

nsp<-25
nyr<-200

#------------------calling extreme copula for LTdep.---------------
xlmat<-retd(n=nyr,d=nsp,rl=-1,mn=0,sdev=1) #LTdep.
range(xlmat)
clmat<-xlmat+5 # this line makes cover data positive
any(clmat<0) # This should be FALSE
#plot(xlmat[,1],xlmat[,2])
#plot(clmat[,1],clmat[,2])
#plot(c(1:nyr),clmat[,1])
tot_ts_l<-apply(clmat,MARGIN = 1, FUN=sum)

surrog_array<-ncsurrog(m=clmat,corpres="spearman",numsurrog=1,plotcheckon=F,resloc="./Results")
surrog_mat_l<-surrog_array[,,1]
#plot(surrog_mat_l[,1],surrog_mat_l[,2])
tot_ts_surrogl<-apply(surrog_mat_l,MARGIN = 1, FUN=sum)

skw_surrogl_each<-apply(surrog_mat_l,MARGIN=2,FUN=skewness,type=2)
skw_real_l_each<-apply(clmat,MARGIN=2,FUN=skewness,type=2)

skw_real_l<-skewness(tot_ts_l,type=2)
skw_surrogl<-skewness(tot_ts_surrogl,type=2)

#------------------calling extreme copula for RTdep.---------------
xrmat<-retd(n=nyr,d=nsp,rl=1,mn=0,sdev=1)  #RTdep.  # you can use -xlmat
range(xrmat)
crmat<-xrmat+5 # this line makes cover data positive
any(crmat<0) # This should be FALSE
#plot(crmat[,1],crmat[,2])
#plot(c(1:nyr),crmat[,1])

tot_ts_r<-apply(crmat,MARGIN = 1, FUN=sum)

surrog_array<-ncsurrog(m=crmat,corpres="spearman",numsurrog=1,plotcheckon=F,resloc="./Results")
surrog_mat_r<-surrog_array[,,1]
tot_ts_surrogr<-apply(surrog_mat_r,MARGIN = 1, FUN=sum)

skw_surrogr_each<-apply(surrog_mat_r,MARGIN=2,FUN=skewness,type=2)
skw_real_r_each<-apply(crmat,MARGIN=2,FUN=skewness,type=2)

skw_real_r<-skewness(tot_ts_r,type=2)
skw_surrogr<-skewness(tot_ts_surrogr,type=2)

#---------------Now plot-----------------------------------------
ylm<-range(c(tot_ts_l,tot_ts_surrogl,tot_ts_r,tot_ts_surrogr))
#ped_ts<-as.data.frame(cbind(tot_ts_surrogl,tot_ts_l,tot_ts_surrogr,tot_ts_r))

pdf("./Results/pedagog_figs/pedagog_ts.pdf",height=4,width=8)
par(mar=c(0.5, 0.5, 0.1, 4), mfcol=c(4,2), oma = c(3.5,3.5, 0.1, 0.1))
#op<-par(mfrow=c(4,1),mar=c(2, 4.1, 0.5, 1.1),mgp=c(1,0.3,0))

# first column for individual biomass of each sp.

s_rg<-round(range(skw_surrogl_each),3)
plot.ts(surrog_mat_l,col=rainbow(nsp,alpha=0.3),ylab="",xlab="",xaxt="n",plot.type = "single",ylim=c(1,10))
legend("topleft",paste("(A) No tail association: range of s[i] = (",s_rg[1],",",s_rg[2],")",sep=""),bty="n",cex=1)

s_rg<-round(range(skw_real_l_each),3)
plot.ts(clmat,col=rainbow(nsp,alpha=0.3),ylab="",xlab="",xaxt="n",plot.type = "single",ylim=c(1,10))
legend("topleft",paste("(B) Left tail association: range of s[i] = (",s_rg[1],",",s_rg[2],")",sep=""),bty="n",cex=1)

s_rg<-round(range(skw_surrogr_each),3)
plot.ts(surrog_mat_r,col=rainbow(nsp,alpha=0.3),ylab="",xlab="",xaxt="n",plot.type = "single",ylim=c(1,10))
legend("topleft",paste("(C) No tail association: range of s[i] = (",s_rg[1],",",s_rg[2],")",sep=""),bty="n",cex=1)

s_rg<-round(range(skw_real_r_each),3)
plot.ts(crmat,col=rainbow(nsp,alpha=0.3),ylab="",xlab="",plot.type = "single",ylim=c(1,10))
legend("topleft",paste("(D) Right tail association: range of s[i] = (",s_rg[1],",",s_rg[2],")",sep=""),bty="n",cex=1)

mtext(paste("Biomass for each of ",nsp," species in the community",sep=""),side=2,line=2.2, adj = 0.05)
mtext("Year",side=1,line=2.2, adj = 0.5)


# second column for aggregated biomass

plot(c(1:nyr),tot_ts_surrogl,type="l",ylim=ylm,col="grey",ylab="",xlab="",xaxt="n")
legend("topleft",paste("(E) No tail association: s = ",round(skw_surrogl,3)),bty="n",cex=1)

plot(c(1:nyr),tot_ts_l,type="l",ylim=ylm,col="black",ylab="",xaxt="n",xlab="")
legend("topleft",paste("(F) Left tail association: s = ",round(skw_real_l,3)),bty="n",cex=1)

plot(c(1:nyr),tot_ts_surrogr,type="l",ylim=ylm,col="grey",ylab="",xaxt="n",xlab="")
legend("topleft",paste("(G) No tail association: s = ",round(skw_surrogr,3)),bty="n",cex=1)

plot(c(1:nyr),tot_ts_r,type="l",ylim=ylm,col="black",xlab="",ylab="")
legend("topleft",paste("(H) Right tail association: s = ",round(skw_real_r,3)),bty="n",cex=1)

mtext("Total community biomass",side=2,line=2.2, adj = -0.5)
mtext("Year",side=1,line=2.2, adj = 0.5)

par(op)
dev.off()