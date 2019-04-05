

#Insert code to import the two datasets
dHAY<- #want a matrix of years by species (all species except those you excluded at the beginning)
dJRG<- #want another such
  
#Insert code to make a vector with "c", "i", "r" values for common, intermediate, and rare species
#src stands for species rarity category
srcHAY<- #a vector with one entry for each species, corresponding to the columns of dHAY
srcJRG<-
  
  
  
#to be deleted when you are done with it Shya
library(copula)
ccop<-claytonCopula(5,2)
d<-rCopula(10000,ccop)
plot(d[,1],d[,2],type="p")
dq<-qnorm(d)
plot(dq[,1],dq[,2],type="p")
hist(dq[,1])
hist(dq[,2])
