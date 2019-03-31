# This plotter function is to visualize a matrix
# Input:
#     z : a matrix
#     posnI_ind : a matrix containing the row and col indices of z for which z indicates indep values
#     posnN_ind : a matrix containing the row and col indices of z for which z indicates sig. neg cor values
#     colrange : a vector containing min and max value of the color range
#     nsm : a matrix containing the label which variable is flipped in case of neg-correlated pair

library(corrplot)

mycorrplot<-function(z,posnI_ind,posnN_ind,colrange,nsm){
  
  col1 <- colorRampPalette(c("blue","white","red")) 
  
  z[is.na(z)]<-mean(colrange)
  diag(z)[1]<-colrange[1] # just to ensure that plot always have specific colorbar range even 
  diag(z)[2]<-colrange[2]     # though all entries are either +ve or -ve
  
  corrplot(z,is.corr = F,col=col1(100),method="color",addgrid.col = "black",
           diag=F,bg = "white",tl.cex=2.5,tl.col = "black",
           cl.cex = 2.2,cl.lim = colrange,
           cl.align.text = "l",cl.ratio = 0.2,p.mat = nsm,insig = "p-value",sig.level = -1)
  
  # colorize as black for diagonal indices
  Dg <- matrix(NA,nrow(z),ncol(z))
  diag(Dg)<- 1 
  
  corrplot(Dg, cl.pos = "n", na.label = " ", add = T,addgrid.col = "black",
           bg = "transparent", tl.col = "transparent",col="black",method="color")
  
  #colorize as yellow for indep posn indices
  if(dim(posnI_ind)[1]!=0){
    
    I <- matrix(NA,nrow(z),ncol(z))
    I[posnI_ind]<- 1 
    
    corrplot(I, cl.pos = "n", na.label = " ", add = T,addgrid.col = "black",
             bg = "transparent", tl.col = "transparent",col="yellow",method="color")
  } 
  
  #colorize as green for -ve correlated (siginificantly) posn indices
  if(dim(posnN_ind)[1]!=0){
    
    N <- matrix(NA,nrow(z),ncol(z))
    N[posnN_ind]<- -1 
    
    corrplot(N, cl.pos = "n", na.label = " ", add = T,addgrid.col = "transparent",
             bg = "transparent", tl.col = "transparent",p.mat = N,sig.level = -2,col="transparent",
             pch=1,pch.col="green",pch.cex = 4,number.cex = 6)
    
  }

}



