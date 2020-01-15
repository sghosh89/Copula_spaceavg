library(testthat)
source("alignranks.R")

dat<-matrix(c(1,2,4,2,6,8),3,2) #note each column here is sorted, as required
sims<-array(c(5,3,8,9,1,4,5,4,9,8,0,1),c(3,2,2))
res<-alignranks(dat,sims)
expect_equal(res,array(c(2,1,4,8,2,6,2,1,4,8,2,6),c(3,2,2)))

dat<-matrix(c(0,3,7,9,3,4,6,9),4,2) #note each column here is sorted, as required
sims<-array(c(1,3,1.5,7,2,7,8,9,2,3,4,9,2.5,1,0,2),c(4,2,2))
res<-alignranks(dat,sims)
expect_equal(res,array(c(0,7,3,9,3,4,6,9,0,3,7,9,9,4,3,6),c(4,2,2)))
