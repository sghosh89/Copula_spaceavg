source("./getmap.R")

d<-readRDS(file="./Results/hays_results/skewness_results/ts_mat_CP_hays.RDS")

# result folder to save surrogates and function plots
resloc<-"./Results/hays_results/skewness_results/pp_surrogs_hays_CP/Reform_maps/"
if (!dir.exists(resloc)){
  dir.create(resloc)
}

# I COPIED THE BELOW CHUNK from SurrogsForHays.R and getmaps.R just to save maps with modified axes labels
set.seed(102)
numpts<-18
parmvals<-seq(from=-1,to=1,length.out=numpts)

for (iind in 1:(dim(d)[2]-1)){
  for (jind in (iind+1):(dim(d)[2])){
  
    print(paste0("iind: ",iind,"; jind: ",jind))
    plotnm<-paste0(resloc,"Map_",iind,"_",jind)
    thisres<-getmap(d[,c(iind,jind)],numpts=numpts,numsims=500,plotnm=NA)
    
    numcres<-thisres$numeric_result
    mnnumcres<-apply(FUN=mean,X=numcres,MARGIN=2)
    quantres<-apply(FUN=quantile,X=numcres,MARGIN=2,probs=c(.025,.25,.5,.75,.975))
    
    pdf(file=paste0(plotnm,".pdf"),width=6,height=6)
    op<-par(mar=c(5,5,1,1))
    plot(parmvals,mnnumcres,type='b',ylim=range(mnnumcres,quantres),
         xlab="Normal-copula parameter, p",ylab=expression(varphi(p)),col="red",
         cex.axis=2,cex.lab=2)
    lines(parmvals,quantres[1,],type="l",lty="dotted")
    lines(parmvals,quantres[2,],type="l",lty="dashed")
    lines(parmvals,quantres[3,],type="l",lty="dashed")
    lines(parmvals,quantres[4,],type="l",lty="dashed")
    lines(parmvals,quantres[5,],type="l",lty="dotted")
    par(op)
    dev.off()
    }
}