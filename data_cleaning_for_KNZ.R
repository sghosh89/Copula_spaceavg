#============== accessed data on 21-Nov-2019 ===========================
# data source for plant species composition
# http://lter.konza.ksu.edu/content/pvc02-plant-species-composition-selected-watersheds-konza-prairie
# data source for fire history for each watershed
# http://lter.konza.ksu.edu/content/kfh011

library(dplyr)
library(tidyr)
#=============================================================================================
#--------creating results folder for knz data--------------------
if(!dir.exists("./Results/knz_results")){
  dir.create("./Results/knz_results/")
}
resloc_knz<-"./Results/knz_results/"

if(!dir.exists("./Results/knz_results/skewness_results")){
  dir.create("./Results/knz_results/skewness_results/")
}
resloc_knz_skw<-"./Results/knz_results/skewness_results/"
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

# we used on purpose the "t" type as DAN did not find "good" surrogates for "f" soil type later

knz_soiltype<-c("t")

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
