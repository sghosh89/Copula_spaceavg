library(copula)
library(e1071)

source("./get_var_ratio.R")
source("./ExtremeTailDep.R")

set.seed(101)
nsp<-15
nyr<-60

#------------------calling copula for major LTdep.---------------
xlmat<-retd(n=nyr,d=nsp,rl=-1,mn=0,sdev=1) #LTdep.
range(xlmat)
clmat<-xlmat+3 # this line makes cover data positive
any(clmat<0) # This should be FALSE

#z<-cor(cc) # z has both +ve and -ve correlation
tot_ts_l<-apply(clmat,MARGIN = 1, FUN=sum)


#------------------calling copula for major RTdep.---------------
#xrmat<-1-xlmat
xrmat<-retd(n=nyr,d=nsp,rl=1,mn=0,sdev=1)
range(xrmat)
crmat<-xrmat+3
any(crmat<0) # This should be FALSE
#crmat[,c(1:6)]<-(1-crmat[,c(1:6)])
tot_ts_r<-apply(crmat,MARGIN = 1, FUN=sum)


var_real_l<-var(tot_ts_l)
var_real_r<-var(tot_ts_r)
skw_real_l<-skewness(tot_ts_l,type=2)
skw_real_r<-skewness(tot_ts_r,type=2)

#---------------Now plot-----------------------------------------
#ylmat<-ceiling(max(clmat,crmat))
ylmat<-10
#ylm<-range(c(tot_ts_l,tot_ts_r))
ylm<-c(0,100)


pdf("./Results/pedagog_figs/pedagog_cv2_skw.pdf",height=3,width=9)
op<-par(mar=c(0.5, 0.5, 0.5, 4), mfcol=c(2,2), oma = c(3.5,4.5, 0.1, 0.1))

# first column for individual biomass of each sp.
plot.ts(clmat,col=rainbow(nsp,alpha=0.3),ylab="",xlab="",xaxt="n",plot.type = "single",ylim=c(0,ylmat))
legend("topleft",paste("(A) species are synchronously \n rare than abundant within the community",sep=""),bty="n",cex=1)

#get variance ratio
phi_cv_l<-get_var_ratio(m=clmat)
phi_cv_l<-round(phi_cv_l,2)
legend("topright",as.expression(bquote(phi[cv]==.(phi_cv_l))),bty="n",cex=1,text.col="blue")

plot.ts(crmat,col=rainbow(nsp,alpha=0.3),ylab="",xlab="",plot.type = "single",ylim=c(0,ylmat))
legend("topleft",paste("(B) species are synchronously \n abundant than rare within the community",sep=""),bty="n",cex=1)

#get variance ratio
phi_cv_r<-get_var_ratio(m=crmat)
phi_cv_r<-round(phi_cv_r,2)
legend("topright",as.expression(bquote(phi[cv]==.(phi_cv_r))),bty="n",cex=1,text.col="blue")

mtext(paste("Individual sp. (i = ",nsp,") biomass \n in the community",sep=""),side=2,line=2.2, adj = 0,cex=0.8)
mtext("Year",side=1,line=2.2, adj = 0.5,cex=0.8)


# second column for aggregated biomass

plot(c(1:nyr),tot_ts_l,type="l",ylim=ylm,col="darkgray",ylab="",xaxt="n",xlab="")
legend("topleft",paste("(C) s = ",round(skw_real_l,3)),bty="n",cex=1)


plot(c(1:nyr),tot_ts_r,type="l",ylim=ylm,col="darkgray",xlab="",ylab="")
legend("topleft",paste("(D) s = ",round(skw_real_r,3)),bty="n",cex=1)

mtext("Total community biomass",side=2,line=2.2, adj = -1.0,cex=0.8)
mtext("Year",side=1,line=2.2, adj = 0.5,cex=0.8)

par(op)
#box(which="figure")
dev.off()


#____________________________________________________
#source("make_tab_stability_assessment.R")
#make_tab_stability(m=clmat,surrogs_given = F) #LT dep. community
#make_tab_stability(m=crmat,surrogs_given = F)  #UT dep. community







