source("./make_tab_stability_assessment.R")
source("./get_var_ratio.R")
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

n<-100000
d<-10

syncconst<-0.1 
compconst<-0.2
sig_p<-matrix(syncconst,d,d)
diag(sig_p)<-1
sig<-cbind(rep(-compconst,d),sig_p)
sig<-rbind(c(1,rep(-compconst,d)),sig)
sig
eigen(sig,symmetric=TRUE,only.values=TRUE)

res_p1<-getdat(n=n,sig=sig) # synchronously right tail dep. within two functional groups, 
                                # but compensatory between groups

#----------- Now generate the data, copula scale, for the second half of the figure
res_p2<-getdat(n=n,sig=sig)
res_p2<-1-res_p2

#------------- Now combine both copulas with standard normal marginals

normres_p1<-qnorm(res_p1) # qnorm makes marginals normal from uniform
normres_p2<-qnorm(res_p2)
range(normres_p1)
normres_p1<-normres_p1-min(normres_p1,normres_p2) 
normres_p2<-normres_p2-mean(normres_p2)+mean(normres_p1) # so that both have same mean
mean(normres_p1)
mean(normres_p2)
range(normres_p2)
#------------------ Stability metrics ------------------
# extremely different skewness ratios though have nearly same variance ratios
(tab1<-make_tab_stability(m=normres_p1,surrogs = NA,surrogs_given = F))
(tab2<-make_tab_stability(m=normres_p2,surrogs = NA,surrogs_given = F))

#--------------- Now plot -----------------------------------------
time_to_show<-c(1:60)
normres_p1_show<-normres_p1[time_to_show,]
normres_p2_show<-normres_p2[time_to_show,]
#ylm<-max(normres_p1_show,normres_p2_show)
ylm<-8
pdf("./Results/pedagog_figs/pedagog_cv2_skw.pdf",height=3,width=9)
op<-par(mar=c(0.5, 0.5, 0.5, 4), mfcol=c(2,3), oma = c(3.5,4.5, 0.1, 0.1))

# first column for individual biomass of each sp.
plot.ts(normres_p2_show,col=c(rep(1,d),rep(2,d)),ylab="",xlab="",xaxt="n",plot.type = "single",ylim=c(0,ylm))
        
text(0,0,"A",adj=c(0,0),cex=1.4)

#get variance ratio
phi_cv_l<-round(tab2$phi_cvsq,2)
cvsq_com<-round(tab2$cvsq_real,4) 
#cvsq_ind<-round(tab2$cvsq_indep,4)
cvsq_ind<-formatC(tab2$cvsq_ind,4,format="f")
phi_loreau<-round(get_var_ratio(m=normres_p2)$loreau_var_ratio,2)
legend(x=38,y=3,c(as.expression(bquote(phi[CV]==.(phi_cv_l))),
                  as.expression(bquote(phi[LdM]==.(phi_loreau)))),
       y.intersp = 1.2,
       bty="n",cex=1,text.col="black")
legend(x=10,y=3.2,
       c(as.expression(bquote(CV[com]^"2"==.(cvsq_com))),
         as.expression(bquote(CV[ind]^"2"==.(cvsq_ind)))),
       y.intersp = 1.2,
       bty="n",cex=1,text.col="black",horiz=F)

plot.ts(normres_p1_show,col=c(rep(1,d),rep(2,d)),ylab="",xlab="",plot.type = "single",ylim=c(0,ylm))
#legend(x=20,y=8,"(B)", bty="n",cex=1.4)
text(0,0,"B",adj=c(0,0),cex=1.4)

#get variance ratio
phi_cv_r<-round(tab1$phi_cvsq,2)
cvsq_com<-round(tab1$cvsq_real,4)
#cvsq_ind<-round(tab1$cvsq_indep,4)
cvsq_ind<-formatC(tab1$cvsq_ind,4,format="f")
phi_loreau<-round(get_var_ratio(m=normres_p1)$loreau_var_ratio,2)
legend(x=38,y=3,c(as.expression(bquote(phi[CV]==.(phi_cv_r))),
                       as.expression(bquote(phi[LdM]==.(phi_loreau)))),
       y.intersp = 1.2,
       bty="n",cex=1,text.col="black")
legend(x=10,y=3.2,
       c(as.expression(bquote(CV[com]^"2"==.(cvsq_com))),
         as.expression(bquote(CV[ind]^"2"==.(cvsq_ind)))),
       y.intersp = 1.2,
       bty="n",cex=1,text.col="black")

mtext(expression("Individual species biomass, x"[i]*"(t)"),side=2,line=2.2, adj = 0,cex=1)
mtext("Year",side=1,line=2.2, adj = 0.5,cex=1)


# second column for aggregated biomass
tot_ts_l<-apply(FUN=sum,MARGIN=1,X=normres_p2)
tot_ts_r<-apply(FUN=sum,MARGIN=1,X=normres_p1)

th_low<-92
th_high<-113
  
# same mean 
mu_l<-mean(tot_ts_l)
mu_r<-mean(tot_ts_r)
mu_l==mu_r # this should be TRUE

# but different variance
v_l<-var(tot_ts_l)
v_r<-var(tot_ts_r)

tot_ts_l_show<-tot_ts_l[time_to_show]
tot_ts_r_show<-tot_ts_r[time_to_show]
ylm<-c(80,145)
#ylm<-range(tot_ts_l_show,tot_ts_r_show)

plot(time_to_show,tot_ts_l_show,type="l",ylim=ylm,col="darkgrey",ylab="",xaxt="n",xlab="")
abline(h=th_low,lty=2)
abline(h=th_high,lty=3)
#legend(x=15,y=64,"(C)", bty="n",cex=1.4)
text(0,ylm[1],"C",adj=c(0,0),cex=1.4)
legend(x=25,y=ylm[2]+5,
       c(as.expression(bquote(skew(x[tot])==.(formatC(tab2$skw_real,2,format="f")))),
         as.expression(bquote(mean(x[tot])==.(round(mu_l,2)))),
                    as.expression(bquote(var(x[tot])==.(round(v_l,2))))),
       bty="n",cex=1,text.col="black",horiz = F)

plot(time_to_show,tot_ts_r_show,type="l",ylim=ylm,col="darkgrey",xlab="",ylab="")
abline(h=th_low,lty=2)
abline(h=th_high,lty=3)
#legend(x=-2.5,y=75,paste("(D) s = ",round(tab1$skw_real,3)),bty="n",cex=1.5)
#legend(x=15,y=64,"(D)", bty="n",cex=1.4)
text(0,ylm[1],"D",adj=c(0,0),cex=1.4)
legend(x=25,y=ylm[2]+5,
       c(as.expression(bquote(skew(x[tot])==.(round(tab1$skw_real,2)))),
         as.expression(bquote(mean(x[tot])==.(round(mu_r,2)))),
         as.expression(bquote(var(x[tot])==.(round(v_r,2))))),
       bty="n",cex=1,text.col="black",horiz = F)
mtext(expression("Total community biomass, x"[tot]*"(t)"),side=2,line=2.2, adj = 0,cex=1)
mtext("Year",side=1,line=2.2, adj = 0.5,cex=1)


# third column for temporal distribution of aggregated biomass
xlm<-c(75,135)
x <- tot_ts_l
y <- hist(x,breaks = 100,plot=F)
pval_low<-sum(x<th_low)/length(x)
pval_hi<-sum(x>th_high)/length(x)
plot(y$breaks,
     c(y$counts,0),
     type="s",col="darkgrey",xlim=xlm,ylim=c(0,max(y$counts)),xaxt="n",
     ylab="Frequency",xlab="Total community Biomass")
abline(v=th_low,lty=2)
abline(v=th_high,lty=3)
#legend(x=38,y=570,"(E)", bty="n",cex=1.4)
text(xlm[1],max(y$counts),"E",adj=c(0,1),cex=1.4)
text(x=th_low,y=1000,paste0("p=",round(pval_low,4)),adj=c(1,0),cex=1.1)
text(x=th_high,y=1000,paste0("p=",round(pval_hi,4)),adj=c(0,0),cex=1.1)

#hist(tot_ts_l,breaks = 100, col="grey",lty="blank")

x <- tot_ts_r
y <- hist(x,breaks = 100,plot=F)
pval_low<-sum(x<th_low)/length(x)
pval_hi<-sum(x>th_high)/length(x)
plot(y$breaks,
     c(y$counts,0),
     type="s",col="darkgrey",xlim=xlm,ylim=c(0,max(y$counts)),
     ylab="Frequency",xlab="Total community Biomass")
abline(v=th_low,lty=2)
abline(v=th_high,lty=3)
mtext("Count",side=2,line=2.2, adj = 1.5,cex=1)
#mtext("Total community biomass",side=1,line=2.2, adj = 0.5,cex=1)
mtext(expression("Total community biomass, x"[tot]*"(t)"),side=1,line=2.4, adj = 0.5,cex=1)
#legend(x=38,y=570,"(F)", bty="n",cex=1.4)
text(xlm[1],max(y$counts),"F",adj=c(0,1),cex=1.4)
#legend(x=20,y=120,paste("p = ",round(pval_low,2)),bty="n",cex=1.1)
#legend(x=50,y=120,paste("p = ",round(pval_hi,2)),bty="n",cex=1.1)
text(x=th_low,y=1000,paste0("p=",round(pval_low,4)),adj=c(1,0),cex=1.1)
text(x=th_high,y=1000,paste0("p=",round(pval_hi,4)),adj=c(0,0),cex=1.1)

par(op)
#box(which="figure")
dev.off()


#x1<-tot_ts_l
#x2<-tot_ts_r
#plot(ecdf(x1))
#plot(ecdf(x2),col="red",add=T)


