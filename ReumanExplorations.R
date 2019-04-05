

#Importing the two dataset matrices of years by species
dHAY<- readRDS("./Results/hays_results/skewness_results/ts_mat_all_sp_146_hays.RDS")# for all species except in "REMOVE" category
dJRG<- readRDS("./Results/jrg_results/skewness_results/ts_mat_all_sp_39_jrg.RDS")# for all species except "BARE","ROCK","LASP"
  
#Importing data with "C", "I", "R" values for common, intermediate, and rare species
category_hays<-readRDS("./Results/hays_results/skewness_results/all_sp_146_hays_category.RDS")  
category_jrg<-readRDS("./Results/jrg_results/skewness_results/all_sp_39_jrg_category.RDS")

# check :
all(colnames(dHAY)==category_hays$sp)==T
all(colnames(dJRG)==category_jrg$sp)==T

#src stands for species rarity category
srcHAY<- category_hays$category #a vector with one entry for each species, corresponding to the columns of dHAY
srcJRG<- category_jrg$category

