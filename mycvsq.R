# Function to calculate coefficient of variation^2
# Arg : x is a vector : here used for total timeseries 
mycvsq<-function(x){
  cvsq<-var(x)/(mean(x))^2
  return(cvsq)
}