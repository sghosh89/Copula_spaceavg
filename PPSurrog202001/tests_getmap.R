library(testthat)
source("getmap.R")

d<-matrix(rnorm(200),100,2)
res<-getmap(d,numpts=20,numsims=100,plotnm="test")
expect_equal(names(res),c("numeric_result","fit_parameters"))
expect_equal(dim(res$numeric_result),c(100,20))

