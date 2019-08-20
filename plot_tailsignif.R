#' A function to plot confidence interval of the statistic Corl - Coru against Spearman correlation
#' 
#' @param spcors A vector of Spearman correlation values : same as given input in tailsignif function
#' @param qtl output matrix of quantiles from tailsignif function
#' @param resloc location folder path to save the plot

plot_tailsignif<-function(spcors,qtl,resloc){
  
  pdf(file=paste0(resloc,"CI_tailsig.pdf"))
  
  plot(c(-1,1),c(0,0),ylim=c(-1,1),xlab="Spearman",ylab="Statistic",type='l',col='red')
  lines(c(0,0),c(-1,1),type="l",col="red")
  
  lines(spcors[spcors>=0],qtl[1,spcors>=0],type='p',pch=16,col=rgb(0,0,0,0.2))
  lines(spcors[spcors>=0],qtl[2,spcors>=0],type='p',pch=16,col=rgb(0,0,0,0.2))
  
  #lines(spcors[spcors<0],qtl[1,spcors<0],type='p',pch=16,col=rgb(0,0,0,0.2))
  lines(spcors[spcors<0],qtl[2,spcors<0],type='p',pch=16,col=rgb(0,0,0,0.2))
  
  dev.off()
}





