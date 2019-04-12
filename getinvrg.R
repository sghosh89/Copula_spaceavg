#Attempts to find a range of acceptable preimages of a point under the fully  
#stochastic version of the map obtained by copmap
#
#Args
#p        An increasing vector of points in [-1,1] - this would have been
#           the p argument of copmap
#cors     The res output of a call to copmap
#imval    The image value to take the preimage of - would be the Pearson 
#           correlation of the x and y inputs to copmap
#
#Output
#The values cors[,i] are samples from a stochastic map evaluated at p[i].
#This function gives the range of values of p for which the 25% and
#75% quantiles of the distribution of values include imval.

getinvrg<-function(p,cors,imval)
{
  #x and y of the function to be inverted
  x<-p
  ylow<-apply(FUN=quantile,X=cors,MARGIN=2,probs=0.25)
  yhi<-apply(FUN=quantile,X=cors,MARGIN=2,probs=0.75)

  #linearly interpolate
  interplow<-approx(x,ylow,n=1000)
  ylow<-interplow$y
  interphi<-approx(x,yhi,n=1000)
  x<-interphi$x
  yhi<-interphi$y
  
  #get result
  res<-c(x[max(which(yhi<=imval))],x[min(which(ylow>=imval))])
  return(res)
}