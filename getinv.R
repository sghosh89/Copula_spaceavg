#Attempts to find the unique preimage of a point under a map  
#obtained by copmap (in expectation)
#
#Args
#p            An increasing vector of points in [-1,1] - this would have been
#               the p argument of copmap
#copmapout    The output of the call to copmap 
#imval        The image value to take the preimage of - would be the correlation 
#               or covariance of the x and y inputs to copmap, depending on the 
#               value of cflag that was used in the call to copmap
#center       Either "mean" or "median" depending on what you want the preimage 
#               of
#
#Output
#Taking the means of each column of copmapout gives a vector of the same length 
#as p. These together determine a function, via linear interpolation between 
#values. The code computes the pre-image under this map, if it exists, of imval. 
#If there is no point in the pre-image, the function returns -Inf. If more than 
#one point, it returns Inf. If exactly one point, it returns that value.

getinv<-function(p,copmapout,imval,center="mean")
{
  #x and y of the function to be inverted
  x<-p
  if (center=="mean")
  {
    y<-apply(FUN=mean,X=copmapout,MARGIN=2)
  }
  if (center="median")
  {
    y<-apply(FUN=median,X=copmapout,MARGIN=2)
  }
  if (!(center %in% c("mean","median")))
  {
    stop("Error in getinv: bad value for center")
  }
  
  res<-numeric(0)
  for (counter in 1:(length(x)-1))
  {
    if (y[counter]==imval && y[counter+1]==imval)
    {
      return(Inf)
    }
    if ((y[counter]<=imval && y[counter+1]>=imval) || 
        (y[counter]>=imval && y[counter+1]<=imval))
    {
      res<-c(res,x[counter]+
               (x[counter+1]-x[counter])*(imval-y[counter])/(y[counter+1]-y[counter]))
    }
  }
  
  if (length(res)==0)
  {
    return(-Inf)
  }
  if (length(res)>1)
  {
    return(Inf)
  }
  return(res)
}