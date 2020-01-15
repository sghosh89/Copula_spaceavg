#============== accessed data on 21-Nov-2019 ===========================
# data source for plant species composition
# http://lter.konza.ksu.edu/content/pvc02-plant-species-composition-selected-watersheds-konza-prairie
# data source for fire history for each watershed
# http://lter.konza.ksu.edu/content/kfh011
rm(list=ls())
seed<-101
library(dplyr)
library(tidyr)
#=============================================================================================
#--------creating results folder for knz data--------------------
if(!dir.exists("./Results/my_knz_results")){
  dir.create("./Results/my_knz_results/")
}
resloc_knz<-"./Results/my_knz_results/"

if(!dir.exists("./Results/my_knz_results/skewness_results")){
  dir.create("./Results/my_knz_results/skewness_results/")
}
resloc_knz_skw<-"./Results/my_knz_results/skewness_results/"
#============================================================================================

# reading raw data 
knz<-read.csv("./Data/KnzData/KNZ_Data_downloaded/knb-lter-knz.69.15/PVC021.csv") 
#unique(knz$Plot)
# [1] 1 2 3 4 5

#unique(knz$Transect)
#[1] A B C D E
#Levels: A B C D E

#unique(knz$SoilType) 
#[1] f t s
#Levels: f s t

#levels(knz$WaterShed)
#[1] "001a" "001c" "001d" "002c" "002d" "004a" "004b" "004f" "00fa" "00fb" "00wa" "00wb" "020a" "020b" "020d" "0spa" "0spb" "0sua"
#[19] "0sub" "n01a" "n01b" "n04a" "n04d" "n20a" "n20b" "r01a" "r01b" "r20a" "r20b"

firehistory<-read.csv("./Data/KnzData/KNZ_Data_downloaded/KFH011.csv")

#subset knz only for watershade "001d" as it is the only one longest available (1983-2018) ungrazed-annually burned watershed
knz_001d<-knz %>% filter(WaterShed=="001d") 

unique(knz_001d$Plot) # 5 plots: 1,2,3,4,5
unique(knz_001d$Transect) # 4 transects: A,B,C,D
unique(knz_001d$SoilType) # 3 soiltypes: f,t,s

my_knz <- knz_001d %>% 
                  # make to small letter in Transect, abbreviated genus and species combined column
                      mutate(Transect = tolower(Transect), 
                             species=paste(AB_genus,AB_species,sep="_"))%>%
                  # make a uniqueID for each observation along with date column
                      mutate(uniqueID=paste(Transect,Plot,SoilType,sep="_"),
                             date=paste(RecYear,RecMonth,RecDay,sep="_"))%>%
                  # choose some required variables
                      select(date,
                             year=RecYear,
                             species_code=SpeCode,
                             species,uniqueID,
                             transect=Transect,
                             plot=Plot,
                             soiltype=SoilType,
                             cover_class=Cover)%>%
                  # add a new column for percent cover value for each class category (see metadata)
                      mutate(midpoint_percent_cover = 
                           recode(cover_class, `1` = 0.5, `2` = 3,`3` = 15,`4` = 37.5,`5` = 62.5,`6` = 85,`7` = 97.5))

#========================================================================================================
# explanation for the next cleaning procedure
# consider a species sampled from a_5_f plot for 1983-2018 
sp_c<-"ambros_psilo"
my_x_sp_a5f<-my_knz%>%filter(species==sp_c,uniqueID=="a_5_f")
dim(my_x_sp_a5f) # it has 39 observations - that means multiple census dates are possible for a year
                     # but we need only one observation per year


     my_knz <- my_knz%>%
                # For species that are sampled on multiple census dates in a year, 
                 # the highest cover class of each plot is used 
                  mutate(date=as.Date(date, "%Y_%m_%d"))%>%
                  tbl_df()%>%
                  group_by(species,year,uniqueID) %>%
                  slice(which.max(cover_class))%>%ungroup()%>%
                # now, drop the date column
                  select(-date)%>%
                # arrange in order
                  arrange(species,year,uniqueID)%>%as.data.frame()


#===================================== difference with Lauren's data ================================================================

# spl1<- unique species list from knz_001d data, spl2<- unique species list from Lauren's data (1983-2015)
# spl1[which((spl1%in%spl2)==F)] # "erecht_hiera" "panicu_capil" two new species observation after 2015, 
                                           # which were not in Lauren's data
#read Lauren's data
#dat<-read.csv("./Data/KnzData/Grassland_group1_master.csv")
#dat_sitewise<-split(dat,dat$site)
#dat_knz<-dat_sitewise$knz

#sort(unique(dat_knz$abundance)) # I have no idea how Lauren got this abundance calculation?
#[1]  0.00  0.10  1.80  3.50  7.55  9.25 15.00 20.50 26.25 33.00 37.50 38.75 50.00 61.25 62.50 67.50 73.75 80.00 85.00 91.25 97.50

#=====================================================================================================================
# check: missing observations

year_c<-1983 #common year
sp_c<-"ambros_psilo"

x1<-my_knz%>%filter(year==year_c)
num_sampled_plot<-length(unique(x1$uniqueID))

x2<-my_knz%>%filter(year==year_c,species==sp_c)
num_reporting<-length(unique(x2$uniqueID))

num_sampled_plot==num_reporting # this is FALSE

# that means though num_sampled_plot are surveyed, sp_c species was not found in all those surveyed plots and 
# reported for plots only with observations

# so, we need to fill in the missing values for sp_c with 0 values for those surveyed plots

my_knz <- my_knz%>%
          # fill in missing observations
             complete(species,nesting(year,uniqueID), fill = list(midpoint_percent_cover = 0))%>%
          # specifying transect,plot,soiltype separately
             mutate(dummyuniqueID=uniqueID)%>%
             separate(col=dummyuniqueID,into=c("transect","plot","soiltype"),sep="_")%>%
          # just rearranging the columns
             select(species_code,
                 species,
                 year,
                 uniqueID,
                 transect,
                 plot,
                 soiltype,
                 cover_class,
                 midpoint_percent_cover)

# check complete: see now they are same 

year_c<-1983 #common year
sp_c<-"ambros_psilo"

x1<-my_knz%>%filter(year==year_c)
num_sampled_plot<-length(unique(x1$uniqueID))

x2<-my_knz%>%filter(year==year_c,species==sp_c)
num_reporting<-length(unique(x2$uniqueID))

num_sampled_plot==num_reporting

#========================================================================================================
# Now calculate the avg. cover for each species and for each year
#       for each 20 transects+plots: (a,b,c,d)+(1:5)combo

as.data.frame(table(my_knz$year,by=my_knz$soiltype))

# But watch out! the above dataframe shows 
# year = 1991 1992 1997 1998 1999 2000 2001 2006 with all 3 soil types: f,s,t
# other years with only 2 soil types: f,t sampled

# so we need to average the 20 plots either with soiltype = "f" or with soiltype = "t"
knz_soiltype<-c("f")

# make a blank dataframe (avg. cover from each species timeseries along each column) to populate
year_list<-sort(unique(my_knz$year)) # 36 years
sp_list<-sort(unique(my_knz$species)) # 125 sp.

ts_all_sp_knz <- matrix(NA,nrow=length(year_list), ncol=length(sp_list))
rownames(ts_all_sp_knz)<- paste("year_",year_list,sep="")
colnames(ts_all_sp_knz)<-sp_list
ts_all_sp_knz<-as.data.frame(ts_all_sp_knz)

for(i_yr in c(1:length(year_list))){
  
  for(i_sp in c(1:length(sp_list))){
    
    # subset the data for a particular year and for a particular species
    x<-my_knz%>%
      filter(year==year_list[i_yr],species==sp_list[i_sp],soiltype%in%knz_soiltype)
    
    # get the average value from all 20 transects+plots combo irrespective of their soiltype
    y<-x%>%summarize(avg_cover=mean(midpoint_percent_cover))%>%as.data.frame()
    
    ts_all_sp_knz[i_yr,i_sp]<-y[1,1]
  }
}


range(ts_all_sp_knz)

#-------saving avg. cover for 36yrs by 125 sp. which has each sp timeseries along each column : knz data-------

saveRDS(ts_all_sp_knz,paste(resloc_knz,"ts_all_sp_knz_soiltype_",knz_soiltype,".RDS",sep=""))


#================ screening for common-intermediate-rare species category ================

ts_all_sp_knz<-readRDS(paste(resloc_knz,"ts_all_sp_knz_soiltype_",knz_soiltype,".RDS",sep=""))

# count on zero values for cover for each species
nyr_0_eachsp<-apply(ts_all_sp_knz,MARGIN = 2,FUN = function(x){sum(x==0)})
nyr_0_eachsp<-as.data.frame(nyr_0_eachsp)

# common sp. for KNZ : these sp present for atleast 35 years, i.e., absent for max 1 years
id_common_sp_knz<-which(nyr_0_eachsp$nyr_0_eachsp<=1) 
ts_common_knz<-ts_all_sp_knz[,id_common_sp_knz]

# rare sp. for KNZ : absent for atleast 34 years, i.e. present max only for 2 years
id_rare_sp_knz<-which(nyr_0_eachsp$nyr_0_eachsp>=34)
ts_rare_knz<-ts_all_sp_knz[,id_rare_sp_knz]

# normal or intermediate sp. for KNZ
id_interm_sp_knz<-which(nyr_0_eachsp$nyr_0_eachsp>1 & nyr_0_eachsp$nyr_0_eachsp<34)
ts_normal_knz<-ts_all_sp_knz[,id_interm_sp_knz]

# sp. category for KNZ
sp_category_knz<-data.frame(sp=rownames(nyr_0_eachsp),category=NA)
sp_category_knz$category[id_common_sp_knz]<-"C"
sp_category_knz$category[id_interm_sp_knz]<-"I"
sp_category_knz$category[id_rare_sp_knz]<-"R"

#---------saving a 125sp by 2 matrix indicating C/I/R category for each knz-sp---------
saveRDS(sp_category_knz,paste(resloc_knz,"all_sp_category_spatialavg_knz_soiltype_",knz_soiltype,".RDS",sep=""))

#===========================Formatting data for tail-asymmetry analysis========================

# Now make the usual format for KNZ data (with common sp and all other merged into a pseudo species) to 
# be used in tail-asymmetry analysis later

knz_spaceavg<-vector("list",1)
names(knz_spaceavg)<-"avg.percent.cover"

sp.screened.data<-vector("list",length(id_common_sp_knz))
names(sp.screened.data)<-rownames(nyr_0_eachsp)[id_common_sp_knz]

for(isp in 1:length(id_common_sp_knz)){
  sp.screened.data[[isp]]<-data.frame(Year=c(1983:2018),Dat=ts_common_knz[,isp])
}

knz_spaceavg$avg.percent.cover<-sp.screened.data

# Append the pseudo species = merged sp. of I & R category
pseudo_knz_IR<-apply(X=ts_all_sp_knz[,which(sp_category_knz$category%in%c("I","R"))],MARGIN = 1,FUN = sum)
pseudo_knz<-data.frame(Year=c(1983:2018),Dat=pseudo_knz_IR)
pseudo_knz<-list(pseudo_knz)

knz_spaceavg$avg.percent.cover<-append(knz_spaceavg$avg.percent.cover,pseudo_knz)
names(knz_spaceavg$avg.percent.cover)[[length(id_common_sp_knz)+1]]<-"pseudo_knz"

#---------saving the spatial avg. data for knz with whigh we will do taildep. analysis later----------------
saveRDS(knz_spaceavg,paste(resloc_knz,"knz_spaceavg_data_CP_soiltype_",knz_soiltype,".RDS",sep=""))

#-----------saving a dataframe with timeseries of common + 1 pseudo (all other merged into 1) sp. for knz------
ts_CP_knz<-cbind(ts_common_knz,pseudo_knz_IR)
saveRDS(ts_CP_knz,paste(resloc_knz_skw,"ts_CP_knz_soiltype_",knz_soiltype,".RDS",sep="")) 

#--------------- time series plot for knz Common sp, Common + Normal(or Intermediate) sp., Common+Normal+Rare sp.-------------------

pdf(paste(resloc_knz_skw,"total_timeseries_soiltype_",knz_soiltype,".pdf",sep=""),height=6,width=6)
op<-par(mar=c(5.1, 5.1, 1.1, 2.1))

total_ts_C<-apply(X=ts_all_sp_knz[,id_common_sp_knz],MARGIN = 1,FUN = sum)
total_ts_CI<-apply(X=ts_all_sp_knz[,sort(c(id_common_sp_knz,id_interm_sp_knz))],MARGIN = 1,FUN = sum)
total_ts_CIR<-apply(X=ts_all_sp_knz[,sort(c(id_common_sp_knz,id_interm_sp_knz,id_rare_sp_knz))],MARGIN = 1,FUN = sum)
total_ts<-cbind(total_ts_C,total_ts_CI,total_ts_CIR)
plot(c(1983:2018),total_ts[,1],ylim=range(total_ts),col=rgb(1,0,0,0.5),type="b",pch=16,xlab="Years",ylab="Total percent cover",xlim=c(1983,2015),cex.lab=2,cex.axis=2)
lines(c(1983:2018),total_ts[,2],type="b",pch=16,col=rgb(0,1,0,0.5))
lines(c(1983:2018),total_ts[,3],type="b",col="black",pch=1)
legend("topright",c("common sp.","common + intermediate sp.","all sp. including rare"),lty=c(1,1,1),
       col=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5),"black"),pch=c(16,16,1),bty="n",cex=1.2)
par(op)
dev.off()

#------plot knz_spaceavg data for all common sp for all yearspan-------------------------------
good_sp<-c(1:length(knz_spaceavg[[1]]))
lensp<-length(good_sp)

summary_knz_commonsp<-data.frame(sp=names(knz_spaceavg$avg.percent.cover),n0=NA,nTies=NA)
pdf(paste(resloc_knz,"rawplot_knz_spaceavg_commonsp_with_pseudosp_avgcover_soiltype_",knz_soiltype,".pdf",sep=""),height = 0.5*lensp,width=0.5*lensp)
op<-par(mfrow=c(6,5),mar=c(5,5,3,3))
for (i in good_sp){
  n0<-sum(knz_spaceavg$avg.percent.cover[[i]]$Dat==0)
  nTies<-sum(duplicated(knz_spaceavg$avg.percent.cover[[i]]$Dat)==T)
  summary_knz_commonsp$n0[i]<-n0
  summary_knz_commonsp$nTies[i]<-nTies
  if(n0==0){
    col1<-rgb(1,0,0,0.3) # these are the sp. present for all years 
  }else{
    col1<-rgb(0,0,1,0.3)
  }
  plot(knz_spaceavg$avg.percent.cover[[i]]$Year,knz_spaceavg$avg.percent.cover[[i]]$Dat,col=col1,pch=19,
       ylim=c(0,max(knz_spaceavg$avg.percent.cover[[i]]$Dat)),xlab="Year (1983-2018)",type="b",ylab="avg. % cover")
  abline(h=0)
  
  mtext(paste0("sp = ",i," : ",names(knz_spaceavg$avg.percent.cover)[i]," ,nT=",nTies,sep=""))
}
par(op)
dev.off()

# ---------------------generate copula plots for all common sp for knz_spaceavg--------------------
source("./vivj_matrix.R")

include_indep<-FALSE

pdf(paste(resloc_knz,"copulaplot_knz_spaceavg_commonsp_with_pseudosp_avgcover_soiltype_",knz_soiltype,".pdf",sep=""),height=2*lensp,width = 2*lensp)
op<-par(mfrow=c(lensp,lensp),mar=c(3,3,3,3), mgp=c(1.5,0.5,0))
for(i in c(1:lensp)){
  for(j in c(1:lensp)){
    vivj_matrix(d_allsp=knz_spaceavg,loc=1,
                i=good_sp[i],j=good_sp[j],level=0.05,
                ploton=T,onbounds=T,lb=0,ub=0.5,include_indep=include_indep)
  }
}
par(op)
dev.off()

#============================== Tail-asymmetry analysis: NPA stats==================================

set.seed(seed) 
source("./NonParamStat.R")

resloc_knz_npa<-paste(resloc_knz,"corstat_knz_spaceavg_results/",sep="")

if(!dir.exists(resloc_knz_npa)){
  dir.create(resloc_knz_npa)
}

resloc2<-paste(resloc_knz_npa,"soiltype_",knz_soiltype,"/",sep="")
if(!dir.exists(resloc2)){
  dir.create(resloc2)
}

resloc<-resloc2
nbin_knz<-2
include_indep<-FALSE 

corstat_knz_spaceavg<-multcall(d_allsp=knz_spaceavg,
                               loc=1,
                               resloc=resloc,
                               good_sp=c(1:length(knz_spaceavg[[1]])),
                               nbin=nbin_knz,include_indep=include_indep)

saveRDS(corstat_knz_spaceavg,paste(resloc,file="corstat_knz_spaceavg_nbin_",nbin_knz,".RDS",sep=''))

#======================= plot NPA stat results for knz data =====================

set.seed(seed)
source("./NonParamStat_matrixplot.R")

resloc<-resloc2

ub<-1/nbin_knz
numpts<-length(knz_spaceavg$avg.percent.cover[[1]]$Year) #for knz data 33years
numsims<-10000
CI<-c(0.025,0.975)
sigtest<-FALSE
include_indep<-FALSE

CorlmCoru_knz_spaceavg<-NonParamStat_matrixplot(data=corstat_knz_spaceavg,
                                                resloc=resloc,tagon=T,
                                                type="lower",wd=15,ht=15,
                                                sigtest=sigtest,ub=ub,numpts=numpts,numsims=numsims,CI=CI,
                                                include_indep=include_indep)
saveRDS(CorlmCoru_knz_spaceavg,paste(resloc,file="CorlmCoru_knz_spaceavg_nbin_",nbin_knz,".RDS",sep=''))
























