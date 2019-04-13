#Attempts to find the unique preimage of a point under the map  
#obtained by copmap
#
#Args
#p        An increasing vector of points in [-1,1] - this would have been
#           the p argument of copmap
#cors     The res output of a call to copmap
#imval    The image value to take the preimage of - would be the Pearson 
#           correlation of the x and y inputs to copmap
#
#Output
#The means of each column of cors gives a vector of the same length as p.
#These together determine a function, via linear interpolation between 
#values. The code computes the complete pre-image under this map,
#if it exists, of imval. If there is no point in the pre-image, the function
#returns -Inf. If more than one point, it returns Inf. If exactly one point,
#it return that value.

getinv<-function(p,cors,imval)
{
  #x and y of the function to be inverted
  x<-p
  y<-apply(FUN=mean,X=cors,MARGIN=2)
  
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