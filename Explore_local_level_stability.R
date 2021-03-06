# This script is written to explore the variance ratio and skewness ratio on local level, i.e., for each quadrat
# in a given plot which is regional (or community) level 

rm(list=ls())
seed<-101

library(tidyverse)
#library(ggplot2)

# First try with Hays data

set.seed(seed)

dir.create("./Results/Explore_local_level_stability/")

#====================================== HAYS DATA ======================================================

# basal cover data

# quadrat sampling ? Yes or not
quadsamp<-read.csv("./Data/HaysData/quadrat_inventory.csv")

quad_info<-read.csv("./Data/HaysData/quadrat_info.csv")
count_quad_grazed<-table(quad_info$Grazing) #36 No grazing, 15 Yes grazed
quad_grazing<-quad_info$quadrat[which(quad_info$Grazing=="Yes")]
quad_grazing<-as.character(quad_grazing) # we will exclude these 15 quadrats where grazing was done

hays_array<-readRDS("./Results/hays_results/hays_array_41yr_51plot_151sp.RDS") # total basal cover for each sp. 
                                                                                   # for each quad for each yr
dim(hays_array) # 41 yrs by 51 quad by 151 sp

id_grz<-which(colnames(hays_array)%in%quad_grazing) 
hays_array_nograzing<-hays_array[,-id_grz,]
dim(hays_array_nograzing) # 41 by 36 by 151 only no-grazed plots are included

# Now, we will rearrange the array so that array has dimension = years by species by quadrats
hays_array_nograzing<-aperm(hays_array_nograzing,perm=c(1,3,2))
dim(hays_array_nograzing) # 41 by 151 by 36

# Note for each non-grazed quadrats for some years there are NA along rows: that means for such years no species 
# have any info as the quadrat is not surveyed for such years
apply(quadsamp[,-id_grz], MARGIN = 2, function(x){sum(is.na(x))})

# Here I consider from year 1943 afterwards as 1932-1942 no "et" community plots are surveyed 
# So, for an equal comparison consider only 1943-1972 periods
quadsamp_43to72_nogrz<-quadsamp[12:41,-id_grz]
apply(quadsamp_43to72_nogrz, MARGIN = 2, function(x){sum(is.na(x))})
# Now, we have last 30 years of data where quadrats were surveyed for most of the years 
# if not, then max for 2 years were not surveyed, 
# so it is reasonable to compare

sampled_yr<-which(rownames(hays_array_nograzing)%in%43:72)
hays_array_nograzing<-hays_array_nograzing[sampled_yr,,]
dim(hays_array_nograzing) #30 yr X 151 sp X 36 quadrats

# Now, delete the non-plant species
spinfo<-read.csv("./Data/HaysData/species_list.csv")
nonplant<-as.character(spinfo$species[which(spinfo$type=="remove")]) 
#These are not plant species as per metadata file (page 8 in pdf) for more info: though I doubt on mixed grass and polygonum spp.]
#"Bare ground"    "Fragment"       "Mixed grass"    "Polygonum spp." "Unknown"
id_nonplant<-which(rownames(hays_array_nograzing[1,,])%in%nonplant)
hays_array_nograzing<-hays_array_nograzing[,-id_nonplant,]
dim(hays_array_nograzing) # 30 years by 146 sp. by 36 plots

# so, we have to first omit those (1 or 2) rows of NA and then 
# need to pass that species time series matrix for a specific plot 
# in the make_tab_stability_assessment() function

source("./make_tab_stability_assessment.R")

tbl_allsp<-c()
for (i in 1: dim(hays_array_nograzing)[3]){
  z1<-hays_array_nograzing[,,i]
  z1<-na.omit(z1) # omit no survey year from each quadrat
  
  # ok, now combine all rare species into a single timeseries:
  # otherwise all 0's along each column give you zero variance and cvsq_real = 0 
  
  # rare species means which occurs only twice in a year
  id_rare<-which(apply(z1, MARGIN=2, function(x){length(which(x==0))})>=(nrow(z1)-2))
  z_common<-z1[,-id_rare]
  z_rare<-z1[,id_rare]
  
  rare_sp<-apply(z_rare, MARGIN=1, sum)
  
  # check on rare_sp
  if(sum(rare_sp)==0){print("rare sp. combined still zero, no variance: watch out!")}
  
  m<-cbind(z_common,rare_sp)
  
  temp<-make_tab_stability(m=m,surrogs=NA,surrogs_given=F)
  
  tbl_allsp<-rbind(tbl_allsp,temp)
}

rownames(tbl_allsp)<-dimnames(hays_array_nograzing)[[3]] # no-grazed quadrat names

hist(tbl_allsp$phi_cvsq,50)
hist(tbl_allsp$phi_skw,50)

quad_info_nogrz<-subset(quad_info,quad_info$Grazing=="No")

if(all(quad_info_nogrz$quadrat==rownames(tbl_allsp))==T){
  tbl_allsp$community<-quad_info_nogrz$community
  tbl_allsp$group<-quad_info_nogrz$group
}


p<-ggplot(tbl_allsp,aes(x=phi_cvsq,y=phi_skw,col=community))+
  geom_point(size=0)+xlim(0,1.2)+ylim(-2.5,1.5)+geom_vline(xintercept=1)+geom_hline(yintercept = c(0,1),lty=c(2,3))+
  theme_bw()
p1<-p+geom_text(aes(label = rownames(tbl_allsp)),
                      size = 3.5)
p1

gd<-tbl_allsp%>%group_by(community)%>%summarize(phi_cvsq=mean(phi_cvsq),
                                                phi_skw=mean(phi_skw))
p1<-p1+geom_point(data=gd,size=7,alpha=0.3)
p1
#phi_cv>1 less stable, right side of vertical line, phi_cv<1 more stable, left side of vertical line
#phi_s when in between [0,1] less stable, i.e. within two horizontal lines
# phi_s>1 : more stable
# ph_s<0 : species reordering???

# most of the cases, points fall in between two horizontal lines 
# where two stability indices differ in their interpretation

#======= what happens if we do on regional(plot) level? =================

# first, we will average cover over all quadrats for each sp. and yr
dim(hays_array_nograzing) # 30yrs by 146sp by 36quadrats
bigmat<-apply(hays_array_nograzing, MARGIN=c(1,2),FUN=mean, na.rm=T) # averaged over 36 non-grazed quadrats
dim(bigmat) # 30yrs by 146sp, should have no NAs

# rare sp. means sp. present max for two years
id_rare<-which(apply(bigmat, MARGIN=2, function(x){length(which(x==0))})>=(nrow(bigmat)-2))
bigmat_rare<-bigmat[,id_rare]

rare_sp<-apply(bigmat_rare, MARGIN=1, sum)

# check on rare_sp
if(sum(rare_sp)==0){print("rare sp. combined still zero, no variance: watch out!")}

bigmat<-cbind(bigmat[,-id_rare],rare_sp)
(temp<-make_tab_stability(m=bigmat,surrogs=NA,surrogs_given=F))
p2<-p1+
  annotate(geom="text", x = temp$phi_cvsq, y = temp$phi_skw, color = "red", label="avg.quads",size=4)
p2

#=================== Now try to do the above but for different community type ==================================

quad_info_nogrz[,"quadrat"] == rownames(tbl_allsp) # check: this should be all TRUE

(cg<-split(quad_info_nogrz,quad_info_nogrz$community)) #4 category: bb, et, lb, sg

tbl_allsp_community_type<-c()
for(j in 1:length(cg)){
  
  id_cg<-which(dimnames(hays_array_nograzing)[3][[1]] %in% as.character(cg[[j]]$quadrat))
  tempo<-hays_array_nograzing[,,id_cg] # array filtering based on community type
  tempo_avg<-apply(tempo, MARGIN=c(1,2), FUN=mean, na.rm=T) # take average over the same category quads
  #dim(tempo_avg)
  #tempo_avg<-na.omit(tempo_avg) # to omit NaN when certain communities are not suevwyed for a particular year 1933
  #dim(tempo_avg)
  ans<-make_tab_stability(m=tempo_avg,surrogs=NA,surrogs_given=F)
  ans$community_type<-names(cg)[j]
  tbl_allsp_community_type<-rbind(tbl_allsp_community_type,ans)
}

tbl_allsp_community_type 
# so, I think, we should  concentrate on some continuous years surveyed for all 4 categories

p3<-p2+
  annotate(geom="text", x = tbl_allsp_community_type$phi_cvsq, y = tbl_allsp_community_type$phi_skw, 
           color ="black", 
           label=paste(tbl_allsp_community_type$community_type,"_avg"),size=4)

p3

ggsave(filename="./Results/Explore_local_level_stability/summary_hays.pdf",p3,height=8,width=12,device="pdf")

# what could you conclude from this plot?
# selective harvesting maybe, we can point out which community 
# type makes the whole community less productive, could they be linked with 
# grouped_species-specific response to abiotic factor? I can check this if I get 
# soil depth data from short grass dominated quadrats, maybe Lauren's Jasper Ridge data set could help?

# I am not sure, but my intuition it could be related to functional diversity

# Here, I can see short grass community differs most in phi_cvsq and phi_s interpretation
# No quadrats from blue stem dominated quadrat show species reordering, i.e., phi_s not negative,
# but it could be because of comperatively lower sample size than other 3 communities

# community type wise I can check on JRG, PARKGRASS, 

#================================================================================================================



















