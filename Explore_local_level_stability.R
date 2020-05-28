# This script is written to explore the variance ratio and skewness ratio on local level, i.e., for each quadrat
# in a given plot which is regional (or community) level 

rm(list=ls())
seed<-101
# First try with Hays data

set.seed(seed)

#====================================== HAYS DATA ======================================================

# basal cover data
cover<-read.csv("./Data/HaysData/allrecords.csv") # this is basal cover

# quadrat sampling ? Yes or not
quadsamp<-read.csv("./Data/HaysData/quadrat_inventory.csv")

quad_info<-read.csv("./Data/HaysData/quadrat_info.csv")
count_quad_grazed<-table(quad_info$Grazing) #36 No grazing, 15 Yes grazed
quad_grazing<-quad_info$quadrat[which(quad_info$Grazing=="Yes")]
quad_grazing<-as.character(quad_grazing) # we will exclude these 15 quadrats where grazing was done

hays_array<-readRDS("./Results/hays_results/hays_array_41yr_51plot_151sp.RDS")
dim(hays_array) # 41 yrs by 51 quad by 151 sp

id_grz<-which(colnames(hays_array)%in%quad_grazing) 
hays_array_nograzing<-hays_array[,-id_grz,]
dim(hays_array_nograzing) # 41 by 36 by 151 only no-grazed plots are included

# Now, we will rearrange the array so that array has dimension = years by species by quadrats
hays_array_nograzing<-aperm(hays_array_nograzing,perm=c(1,3,2))
dim(hays_array_nograzing) # 41 by 151 by 36

# Note for each quadrats for some years there are NA along rows: that means for such years no species 
# have any info as the plot is not surveyed for such years

# so, we have to first omit those rows of NA and then 
# need to pass that species time series matrix for a specific plot 
# in the make_tab_stability_assessment() function

source("./make_tab_stability_assessment.R")

tbl_allsp<-c()
for (i in 1: dim(hays_array_nograzing)[3]){
  z1<-hays_array_nograzing[,,i]
  z1<-na.omit(z1) # omit no survey year from each quadrat
  temp<-make_tab_stability(m=z1,surrogs=NA,surrogs_given=F)
  tbl_allsp<-rbind(tbl_allsp,temp)
}

rownames(tbl_allsp)<-dimnames(hays_array_nograzing)[[3]] # no-grazed quadrat names

hist(tbl_allsp$phi_cvsq,50)
hist(tbl_allsp$phi_skw,50)

tbl_allsp[which(tbl_allsp$phi_cvsq<1 & tbl_allsp$phi_skw>1),]
# only 4 quadrats have similar conclusion that they are more stable than their indep. communities

#=================== Now try to do the above but for different community type ==================================
quad_no_grz<-quad_info[which(quad_info$Grazing=="No"),]
quad_no_grz[,"quadrat"] == rownames(tbl_allsp) # check: this should be all TRUE

(cg<-split(quad_no_grz,quad_no_grz$community)) #4 category: bb, et, lb, sg

tbl_allsp_community_type<-c()
for(j in 1:length(cg)){
  
  id_cg<-which(dimnames(hays_array_nograzing)[3][[1]] %in% as.character(cg[[j]]$quadrat))
  tempo<-hays_array_nograzing[,,id_cg] # array filtering based on community type
  tempo_avg<-apply(tempo, MARGIN=c(1,2), FUN=mean, na.rm=T) # take average over the same category quads
  tempo_avg<-na.omit(tempo_avg) # to omit NaN when certain communities are not suevwyed for a particular year
  ans<-make_tab_stability(m=tempo_avg,surrogs=NA,surrogs_given=F)
  ans$community_type<-names(cg)[j]
  tbl_allsp_community_type<-rbind(tbl_allsp_community_type,ans)
}

tbl_allsp_community_type # for bb and et community it reversed
# i.e. stability from phi_cv and phi_s are opposite, now for "et" I am not surprised as for ecotone, 
# there is higher chance of species reordering but I don't know why that also happens for "bb"

#========================================= KONZA PRAIRIE ================================================



















