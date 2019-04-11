context("getinv")
source("getinv.R")

test_that("test a case with one point in the preimage",{
  p<-c(-1,-.5,0,.5,1)
  cors<-matrix(rep(p,each=10),10,length(p))
  imval<-.75
  res<-getinv(p,cors,imval)
  expect_equal(res,.75)
})

test_that("test a case with no preimage",{
  p<-c(-1,-.5,0,.5)
  cors<-matrix(rep(p,each=10),10,length(p))
  imval<-.75
  res<-getinv(p,cors,imval)
  expect_equal(res,-Inf)
})

test_that("test a case with more than one point in the preimage",{
  p<-c(-1,-.5,0,.5,1)
  cors<-matrix(rep(c(-1,-.5,0,-.5,0),each=10),10,length(p))
  imval<-(-.25)
  res<-getinv(p,cors,imval)
  expect_equal(res,Inf)
})

