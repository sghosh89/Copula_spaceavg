context("copmap")
source("copmap.R")

test_that("tests for consistency with what I know should happen",{
  set.seed(101)
  x<-rnorm(100)
  y<-(1/sqrt(2))*x+(1/sqrt(2))*rnorm(100)
  x<-5*x
  y<-5*y
  p<-c(-1,-.5,0,.5,1)  
  numreps<-5

  #test things like the class and the size of the matrix result
  res<-copmap(x,y,p,numreps,"cov")
  expect_equal(class(res),"matrix")
  expect_equal(dim(res),c(numreps,length(p)))

  #test extreme values, for which I know what the result should be
  h<-cov(sort(x),sort(y))
  expect_true(all(res[,5]==h))
  h<-cov(sort(x),rev(sort(y)))
  expect_true(all(res[,1]==h))
  
  #similar tests for "cor" input
  res<-copmap(x,y,p,numreps,"cor")
  expect_equal(class(res),"matrix")
  expect_equal(dim(res),c(numreps,length(p)))
  expect_true(all(res>=-1))
  expect_true(all(res<=1))
  
  h<-cor(sort(x),sort(y))
  expect_true(all(res[,5]==h))
  h<-cor(sort(x),rev(sort(y)))
  expect_true(all(res[,1]==h))
})

#would be good to have some more exacting tests than this