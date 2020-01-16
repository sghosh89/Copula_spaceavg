library(testthat)
source("pwlin.R")

#test pwlin
set.seed(101)
xbp<-sort(runif(5))
ybp<-runif(5)
x<-seq(from=0,to=1,length.out=100)
res<-pwlin(x,xbp,ybp)

expect_equal(sum(x<xbp[1])+sum(x>xbp[length(xbp)]),sum(is.na(res)))
plot(xbp,ybp,type="l")
points(x,res,type="p",col="red",pch=20,cex=.5)

#test that you can do inverses
set.seed(101)
xbp<-sort(runif(5))
ybp<-sort(runif(5))
x<-seq(from=xbp[1],to=xbp[length(xbp)],length.out=100)
res<-pwlin(x,xbp,ybp)
ires<-pwlin(res,ybp,xbp)
expect_equal(x,ires)

