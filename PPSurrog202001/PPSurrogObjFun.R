#This is the objective function optimized to search for Pearson preserving
#surrogates for a given dataset
#
#Args
#d          The dataset, a matrix with time series in the columns
#fitparms   The parameters of the maps obtained by getmap

PPSurrogObj<-function(d,fitparms)