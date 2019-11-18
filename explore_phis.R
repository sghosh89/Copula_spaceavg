source("./make_tab_stability_assessment.R")
library(copula)
library(VineCopula)
set.seed(seed=101)

nsp<-25 # number of species in the community
nyr<-1000 # number of years should be very long to see the effect
srho<-0.875 # spearman's correlation

# Now call multivariate Clayton with rho=srho
copC<-claytonCopula(3)
parC<-iRho(copC,rho=srho)

cc<-claytonCopula(par=parC,dim=nsp)
my_com<-rCopula(nyr,cc)  

plot(my_com[,1],my_com[,2])
hist(my_com[,1]) # should be uniform marginal

# Now give a slight right skewed transformation to each marginals of the left-tail dep. community copula
fn_skw_trans <- function(x,gshape,gscale){qgamma(x,shape=gshape,scale=gscale)}
gshape<-10
gscale<-1
my_com_trans<-apply(my_com, MARGIN=2, FUN=fn_skw_trans, gshape=gshape, gscale=gscale)

# Now check if each marginal is slightly right skewed or not
plot(my_com_trans[,1],my_com_trans[,2])
hist(my_com_trans[,1]) # should be slightly right skewed

# calculate cvsq and skewness
#ans0<-make_tab_stability(m=my_com,surrogs=NA,surrogs_given=F)
ans<-make_tab_stability(m=my_com_trans,surrogs=NA,surrogs_given=F)
# see the phi_skw = -ve, means this effects come from tail-dependence (left for Clayton) only not from the 
# right skewed marginals.
















