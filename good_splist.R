good_splist<-function(d_allsp,loc,data_pt_thrs){
  
  lensp<-length(d_allsp[[1]])
  
  finitedat_ind<-vector(mode = "list", length = lensp)
  num_datapt<-c()
  for(sp in c(1:lensp)){
    finitedat_ind[[sp]]<-which(!is.na(d_allsp[[loc]][[sp]]$Dat))
    num_datapt<-c(num_datapt,length(finitedat_ind[[sp]]))
  }
  
  good_sp<-which(num_datapt>=data_pt_thrs)
  return(good_sp)
  
}