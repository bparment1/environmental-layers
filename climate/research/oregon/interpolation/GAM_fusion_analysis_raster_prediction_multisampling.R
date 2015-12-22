#########################    Raster prediction    ####################################
############################ Interpolation of temperature for given processing region ##########################################
#This script interpolates temperature values using MODIS LST, covariates and GHCND station data.                      
#It requires the text file of stations and a shape file of the study area.           
#Note that the projection for both GHCND and study area is lonlat WGS84.       
#Options to run this program are:
#1) Multisampling: vary the porportions of hold out and use random samples for each run
#2)Constant sampling: use the same sample over the runs
#3)over dates: run over for example 365 dates without mulitsampling
#4)use seed number: use seed if random samples must be repeatable
#5)possibilty of running single and multiple time scale methods:
   # gam_daily, kriging_daily,gwr_daily,gam_CAI,gam_fusion,kriging_fusion,gwr_fusion and other options added.
#For multiple time scale methods, the interpolation is done first at the monthly time scale then delta surfaces are added.
#AUTHOR: Benoit Parmentier                                                                        
#CREATED ON: 04/01/2013  
#MODIFIED ON: 12/29/2015  
#PROJECT: NCEAS INPLANT: Environment and Organisms --TASK#568--     
#
# TO DO:

###################################################################################################

raster_prediction_fun <-function(list_param_raster_prediction){

  ##Function to predict temperature interpolation with 21 input parameters
  #9 parameters used in the data preparation stage and input in the current script
  #1)list_param_data_prep: used in earlier code for the query from the database and extraction for raster brick
  #2)infile_monthly: monthly averages with covariates for GHCND stations obtained after query
  #3)infile_daily: daily GHCND stations with covariates, obtained after query
  #4)infile_locs: vector file with station locations for the processing/study area (ESRI shapefile)
  #5)infile_covariates: raster covariate brick, tif file
  #6)covar_names: covar_names #remove at a later stage...
  #7)var: variable being interpolated-TMIN or TMAX
  #8)out_prefix: output suffix added to files
  #9)CRS_locs_WGS84
  #10)screen_data_training
  #
  #6 parameters for sampling function
  #10)seed_number
  #11)nb_sample
  #12)step
  #13)constant
  #14)prop_minmax
  #15)seed_number_month
  #16)nb_sample_month
  #17)step_month
  #18)constant_month
  #19)prop_minmax_month
  #20)dates_selected
  #
  #13 additional parameters for monthly climatology and more
  #21)list_models: model formulas in character vector
  #22)lst_avg: LST climatology name in the brick of covariate--change later
  #23)in_path
  #24)out_path: path to directory containing daily data
  #25)out_path_clim: path to the directory containing climatology data
  #25)script_path: path to script
  #26)interpolation_method: c("gam_fusion","gam_CAI") #other otpions to be added later
  #27) use_clim_image
  #28) join_daily
  #29)list_models2: models' formulas as string vector for daily devation
  #30)interp_method2: intepolation method for daily devation step
  #31)num_cores: How many cores to use
  #32)max_mem: Max memory to use for predict step
  #33)reg_outline: shapefile with region outline used to create nc output file
  
  ###Loading R library and packages     
  
  library(gtools)                                         # loading some useful tools 
  library(mgcv)                                           # GAM package by Simon Wood
  library(sp)                                             # Spatial pacakge with class definition by Bivand et al.
  library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
  library(rgdal)                               # GDAL wrapper for R, spatial utilities
  library(gstat)                               # Kriging and co-kriging by Pebesma et al.
  library(fields)                             # NCAR Spatial Interpolation methods such as kriging, splines
  library(raster)                              # Hijmans et al. package for raster processing
  library(rasterVis)
  library(parallel)                            # Urbanek S. and Ripley B., package for multi cores & parralel processing
  library(reshape)
  library(plotrix)
  library(maptools)
  library(gdata) #Nesssary to use cbindX
  library(automap)  #autokrige function
  library(spgwr)   #GWR method
  
  ### Parameters and arguments
  #PARSING INPUTS/ARGUMENTS
#   
#   names(list_param_raster_prediction)<-c("list_param_data_prep",
#                                          "seed_number","nb_sample","step","constant","prop_minmax","dates_selected",
#                                          "list_models","lst_avg","in_path","out_path","script_path",
#                                          "interpolation_method")
  
  #9 parameters used in the data preparation stage and input in the current script
  list_param_data_prep<-list_param_raster_prediction$list_param_data_prep
  infile_monthly<-list_param_data_prep$infile_monthly
  infile_daily<-list_param_data_prep$infile_daily
  infile_locs<-list_param_data_prep$infile_locs
  infile_covariates<-list_param_data_prep$infile_covariates #raster covariate brick, tif file
  covar_names<- list_param_data_prep$covar_names #remove at a later stage...
  var<-list_param_data_prep$var
  out_prefix<-list_param_data_prep$out_prefix
  CRS_locs_WGS84<-list_param_data_prep$CRS_locs_WGS84

  
  #6 parameters for sampling function
  seed_number<-list_param_raster_prediction$seed_number
  nb_sample<-list_param_raster_prediction$nb_sample
  step<-list_param_raster_prediction$step
  constant<-list_param_raster_prediction$constant
  prop_minmax<-list_param_raster_prediction$prop_minmax
  dates_selected<-list_param_raster_prediction$dates_selected
  
  seed_number_month <-list_param_raster_prediction$seed_number_month
  nb_sample_month <-list_param_raster_prediction$nb_sample_month
  step_month <-list_param_raster_prediction$step_month
  constant_month <-list_param_raster_prediction$constant_month
  prop_minmax_month <-list_param_raster_prediction$prop_minmax_month
  
  #6 additional parameters for monthly climatology and more
  list_models<-list_param_raster_prediction$list_models
  list_models2<-list_param_raster_prediction$list_models2
  interp_method2 <- list_param_raster_prediction$interp_method2
  
  lst_avg<-list_param_raster_prediction$lst_avg
  out_path<-list_param_raster_prediction$out_path #daily prediction path
  out_path_clim <- list_param_raster_prediction$out_path_clim #clim prediction path
  script_path<-list_param_raster_prediction$script_path
  interpolation_method<-list_param_raster_prediction$interpolation_method
  screen_data_training <-list_param_raster_prediction$screen_data_training
  
  use_clim_image <- list_param_raster_prediction$use_clim_image # use predicted image as a base...rather than average Tmin at the station for delta
  join_daily <- list_param_raster_prediction$join_daily # join monthly and daily station before calucating delta

  #cores and memory usage options
  num_cores <- list_param_raster_prediction$num_cores
  max_mem<- as.numeric(list_param_raster_prediction$max_mem)
 
  #Get the region outline
  reg_outline<-list_param_raster_prediction$reg_outline

  #rasterOptions(maxmemory=max_mem,timer=TRUE,chunksize=1e+08)
  rasterOptions(timer=TRUE,chunksize=5e+05)  
  #rasterOptions(timer=TRUE)

  setwd(out_path) #note that this is now path to daily dir (with the name of the year...)
  
  ###################### START OF THE SCRIPT ########################
   
  #This should not be set here...? master script, modify for precip
  if (var=="TMAX"){
    y_var_name<-"dailyTmax"
    y_var_month<-"TMax"
  }
  if (var=="TMIN"){
    y_var_name<-"dailyTmin"
    y_var_month <-"TMin"
  }
  
  ################# CREATE LOG FILE #####################
  
  #create log file to keep track of details such as processing times and parameters.
  
  #log_fname<-paste("R_log_raster_prediction",out_prefix, ".log",sep="")
  log_fname<-paste("R_log_raster_prediction",out_prefix, ".log",sep="")
  #sink(log_fname) #create new log file
  file.create(file.path(out_path,log_fname)) #create new log file
  
  time1<-proc.time()    #Start stop watch
  
  cat(paste("Starting script at this local Date and Time: ",as.character(Sys.time()),sep=""),
             file=log_fname,sep="\n")
  cat("Starting script process time:",file=log_fname,sep="\n",append=TRUE)
  cat(as.character(time1),file=log_fname,sep="\n",append=TRUE)    
  
  ############### Make nc file from outline ############
  

  ############### READING INPUTS: DAILY STATION DATA AND OTHER DATASETS  #################
  #Takes too long to read shapefiles. Let's try using saveRDS(object,*.rds) and readRDS()
  infile_daily_rds <-sub(".shp",".rds",infile_daily)

  if (file.exists(infile_daily_rds) == TRUE){
     ghcn<-readRDS(infile_daily_rds)
     CRS_interp<-proj4string(ghcn)
  }else{
    ghcn<-readOGR(dsn=dirname(infile_daily),layer=sub(".shp","",basename(infile_daily)))
    CRS_interp<-proj4string(ghcn)                       #Storing projection information (ellipsoid, datum,etc.)
    saveRDS(ghcn,infile_daily_rds)
  }

  infile_locs_rds<-sub(".shp",".rds",infile_locs)
  if (file.exists(infile_locs_rds) == TRUE){
    stat_loc<-readRDS(infile_locs_rds) 
  }else{
    stat_loc<-readOGR(dsn=dirname(infile_locs),layer=sub(".shp","",basename(infile_locs)))
    saveRDS(stat_loc,infile_locs_rds)
  }
  
  #dates2 <-readLines(file.path(in_path,dates_selected)) #dates to be predicted, now read directly from the file  
  
  #Should clean this up, reduce the number of if
  if (dates_selected==""){
    dates<-as.character(sort(unique(ghcn$date))) #dates to be predicted 
  }
  if (dates_selected!=""){
    dates<-dates_selected #dates to be predicted 
  }
  if(class(dates_selected)=="numeric"){ #select n every  observation, may change this later.
    dates<-as.character(sort(unique(ghcn$date))) #dates to be predicted 
    dates <- dates[seq(1, length(dates), dates_selected)]
  }
  #Reading in covariate brickcan be changed...
  
  s_raster<-brick(infile_covariates)                   #read in the data brck
  names(s_raster)<-covar_names               #Assigning names to the raster layers: making sure it is included in the extraction
    
  #Reading monthly data
  infile_monthly_rds<-sub(".shp",".rds",infile_monthly)
  if (file.exists(infile_monthly_rds) == TRUE) {
     dst<-readRDS(infile_monthly_rds)
  }else{
   dst<-readOGR(dsn=dirname(infile_monthly),layer=sub(".shp","",basename(infile_monthly)))
   saveRDS(dst,infile_monthly_rds) 
  }

  #construct date based on input end_year !!!
  day_tmp <- rep("15",length=nrow(dst))
  year_tmp <- rep(as.character(end_year),length=nrow(dst))
  #dates_month <-do.call(paste,c(list(day_tmp,sprintf( "%02d", dst$month ),year_tmp),sep="")) #reformat integer using formatC or sprintf
  dates_month <-do.call(paste,c(list(year_tmp,sprintf( "%02d", dst$month ),day_tmp),sep="")) #reformat integer using formatC or sprintf
  
  dst$date <- dates_month
  
  ########### CREATE SAMPLING -TRAINING AND TESTING STATIONS ###########

  #dates #list of dates for prediction
  #ghcn_name<-"ghcn" #infile daily data 
  
  list_param_sampling<-list(seed_number,nb_sample,step,constant,prop_minmax,dates,ghcn)
  #list_param_sampling<-list(seed_number,nb_sample,step,constant,prop_minmax,dates,ghcn_name)
  names(list_param_sampling)<-c("seed_number","nb_sample","step","constant","prop_minmax","dates","ghcn")
  
  #run function, note that dates must be a character vector!! Daily sampling
  sampling_obj<-sampling_training_testing(list_param_sampling)

  #Now run monthly sampling
  dates_month<-as.character(sort(unique(dst$date)))
  list_param_sampling<-list(seed_number_month,nb_sample_month,step_month,constant_month,prop_minmax_month,dates_month,dst)
  #list_param_sampling<-list(seed_number,nb_sample,step,constant,prop_minmax,dates,ghcn_name)
  names(list_param_sampling)<-c("seed_number","nb_sample","step","constant","prop_minmax","dates","ghcn")
  
  sampling_month_obj <- sampling_training_testing(list_param_sampling)
  
  #save(sampling_month_obj,file="/nobackupp4/aguzman4/climateLayers/output10Deg/reg1/35.0_-115.0/test.RData")

  ########### PREDICT FOR MONTHLY SCALE  ##################
  #First predict at the monthly time scale: climatology
  cat("Predictions at monthly scale:",file=log_fname,sep="\n", append=TRUE)
  cat(paste("Local Date and Time: ",as.character(Sys.time()),sep=""),
      file=log_fname,sep="\n")
  t1<-proc.time()
  j=12
  
  ###Changes 12/21/2015
  #browser() #Missing out_path for gam_fusion list param!!!
  #if (interpolation_method=="gam_fusion"){
  if (interpolation_method %in% c("gam_fusion","kriging_fusion","gwr_fusion")){
    clim_method_mod_obj_file <- file.path(out_path_clim,paste(interpolation_method,"_mod_",y_var_name,out_prefix,".RData",sep=""))
    if(!file.exists(clim_method_mod_obj_file)){
      list_param_runClim_KGFusion<-list(j,s_raster,covar_names,lst_avg,list_models,dst,sampling_month_obj,var,y_var_name, out_prefix,out_path_clim)
      names(list_param_runClim_KGFusion)<-c("list_index","covar_rast","covar_names","lst_avg","list_models","dst","sampling_month_obj","var","y_var_name","out_prefix","out_path")
      #debug(runClim_KGFusion)
      #test<-runClim_KGFusion(1,list_param=list_param_runClim_KGFusion)
      clim_method_mod_obj<-mclapply(1:length(sampling_month_obj$ghcn_data), list_param=list_param_runClim_KGFusion, runClim_KGFusion,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
    
      save(clim_method_mod_obj,file= file.path(out_path_clim,paste(interpolation_method,"_mod_",y_var_name,out_prefix,".RData",sep="")))
      #Use function to extract list
    }else{
      clim_method_mod_obj <- load_obj(clim_method_mod_obj_file) #load the existing file
    }
    #Get relevant data
    list_tmp<-vector("list",length(clim_method_mod_obj))

    for (i in 1:length(clim_method_mod_obj)){
      tmp<-clim_method_mod_obj[[i]]$clim
      list_tmp[[i]]<-tmp
    }
    
    clim_yearlist<-list_tmp
  }
  
  if (interpolation_method %in% c("gam_CAI","kriging_CAI", "gwr_CAI")){
    clim_method_mod_obj_file <- file.path(out_path_clim,paste(interpolation_method,"_mod_",y_var_name,out_prefix,".RData",sep=""))
    if(!file.exists(clim_method_mod_obj_file)){
      num_cores2 = as.integer(num_cores) + 2
      list_param_runClim_KGCAI<-list(j,s_raster,covar_names,lst_avg,list_models,dst,sampling_month_obj,var,y_var_name, out_prefix,out_path_clim)
      names(list_param_runClim_KGCAI)<-c("list_index","covar_rast","covar_names","lst_avg","list_models","dst","sampling_month_obj","var","y_var_name","out_prefix","out_path")
      clim_method_mod_obj<-mclapply(1:length(sampling_month_obj$ghcn_data), list_param=list_param_runClim_KGCAI, runClim_KGCAI,mc.preschedule=FALSE,mc.cores = num_cores2) #This is the end bracket from mclapply(...) statement
      #test<-runClim_KGCAI(1,list_param=list_param_runClim_KGCAI)

      save(clim_method_mod_obj,file= file.path(out_path_clim,paste(interpolation_method,"_mod_",y_var_name,out_prefix,".RData",sep="")))
    }else{
      clim_method_mod_obj <- load_obj(clim_method_mod_obj_file) #load the existing file
    }
    #Now get relevant data
    list_tmp<-vector("list",length(clim_method_mod_obj))
    for (i in 1:length(clim_method_mod_obj)){
      tmp<-clim_method_mod_obj[[i]]$clim
      list_tmp[[i]]<-tmp
    }
    clim_yearlist<-list_tmp
  }
  t2<-proc.time()-t1
  cat(as.character(t2),file=log_fname,sep="\n", append=TRUE)
  
  #Getting rid of raster temp files
  removeTmpFiles(h=0)
  

  ################## PREDICT AT DAILY TIME SCALE #################
  #Predict at daily time scale from single time scale or multiple time scale methods: 2 methods availabe now
  #put together list of clim models per month...
  #rast_clim_yearlist<-list_tmp
  
  #Second predict at the daily time scale: delta
  
  #method_mod_obj<-mclapply(1:1, runGAMFusion,mc.preschedule=FALSE,mc.cores = 1) #This is the end bracket from mclapply(...) statement

  #Set raster options for daily steps
  #rasterOptions(timer=TRUE,chunksize=1e+05)
  #This is for high station count areas
  rasterOptions(timer=TRUE,chunksize=1e+04)
  #rasterOptions(timer=TRUE,chunksize=1e+03)

  cat("Predictions at the daily scale:",file=log_fname,sep="\n", append=TRUE)
  t1<-proc.time()
  cat(paste("Local Date and Time: ",as.character(Sys.time()),sep=""),
      file=log_fname,sep="\n")
  
  #TODO : Same call for all functions!!! Replace by one "if" for all multi time scale methods...
  #The methods could be defined earlier as constant??
  #Create data.frame combining sampling at daily and monthly time scales:
  
  daily_dev_sampling_dat <- combine_sampling_daily_monthly_for_daily_deviation_pred(sampling_obj,sampling_month_obj)
    
  #use_clim_image<- TRUE
  #use_clim_image<-FALSE
  #join_daily <- FALSE
  #join_daily <- TRUE
  
  if (interpolation_method %in% c("gam_CAI","kriging_CAI","gwr_CAI","gam_fusion","kriging_fusion","gwr_fusion")){
    #input a list:note that ghcn.subsets is not sampling_obj$data_day_ghcn
    i<-1
    list_param_run_prediction_daily_deviation <-list(i,clim_yearlist,daily_dev_sampling_dat,sampling_month_obj,sampling_obj,dst,list_models2,interp_method2,
                                                     s_raster,use_clim_image,join_daily,var,y_var_name, interpolation_method,out_prefix,out_path)
    names(list_param_run_prediction_daily_deviation)<-c("list_index","clim_yearlist","daily_dev_sampling_dat","sampling_month_obj","sampling_obj","dst","list_models2","interp_method2",
                                                        "s_raster","use_clim_image","join_daily","var","y_var_name","interpolation_method","out_prefix","out_path")
    #method_mod_obj<-mclapply(1:length(sampling_obj$ghcn_data),list_param=list_param_run_prediction_daily_deviation,run_prediction_daily_deviation,mc.preschedule=FALSE,mc.cores = 9) #This is the end bracket from mclapply(...) statement
    #debug(run_prediction_daily_deviation)
    #test <- run_prediction_daily_deviation(1,list_param=list_param_run_prediction_daily_deviation) #This is the end bracket from mclapply(...) statement
    #test <- mclapply(1:9,list_param=list_param_run_prediction_daily_deviation,run_prediction_daily_deviation,mc.preschedule=FALSE,mc.cores = 9) #This is the end bracket from mclapply(...) statement
    
    method_mod_obj<-mclapply(1:nrow(daily_dev_sampling_dat),list_param=list_param_run_prediction_daily_deviation,run_prediction_daily_deviation,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
    save(method_mod_obj,file= file.path(out_path,paste("method_mod_obj_",interpolation_method,"_",y_var_name,out_prefix,".RData",sep="")))
    
  }
  
  removeTmpFiles(h=0)

  #TODO : Same call for all functions!!! Replace by one "if" for all daily single time scale methods...
  if (interpolation_method=="gam_daily"){
    #input a list:note that ghcn.subsets is not sampling_obj$data_day_ghcn
    i<-1
    list_param_run_prediction_gam_daily <-list(i,s_raster,covar_names,lst_avg,list_models,dst,screen_data_training,var,y_var_name, sampling_obj,interpolation_method,out_prefix,out_path)
    names(list_param_run_prediction_gam_daily)<-c("list_index","covar_rast","covar_names","lst_avg","list_models","dst","screen_data_training","var","y_var_name","sampling_obj","interpolation_method","out_prefix","out_path")
    #test <- runGAM_day_fun(1,list_param_run_prediction_gam_daily)
    
    method_mod_obj<-mclapply(1:length(sampling_obj$ghcn_data),list_param=list_param_run_prediction_gam_daily,runGAM_day_fun,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
    #method_mod_obj<-mclapply(1:22,list_param=list_param_run_prediction_gam_daily,runGAM_day_fun,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
    
    save(method_mod_obj,file= file.path(out_path,paste("method_mod_obj_",interpolation_method,"_",y_var_name,out_prefix,".RData",sep="")))
    
  }
  
  if (interpolation_method=="kriging_daily"){
    #input a list:note that ghcn.subsets is not sampling_obj$data_day_ghcn
    i<-1
    list_param_run_prediction_kriging_daily <-list(i,s_raster,covar_names,lst_avg,list_models,dst,var,y_var_name, sampling_obj,interpolation_method,out_prefix,out_path)
    names(list_param_run_prediction_kriging_daily)<-c("list_index","covar_rast","covar_names","lst_avg","list_models","dst","var","y_var_name","sampling_obj","interpolation_method","out_prefix","out_path")
    #test <- runKriging_day_fun(1,list_param_run_prediction_kriging_daily)
    method_mod_obj<-mclapply(1:length(sampling_obj$ghcn_data),list_param=list_param_run_prediction_kriging_daily,runKriging_day_fun,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
    #method_mod_obj<-mclapply(1:18,list_param=list_param_run_prediction_kriging_daily,runKriging_day_fun,mc.preschedule=FALSE,mc.cores = 9) #This is the end bracket from mclapply(...) statement
    
    save(method_mod_obj,file= file.path(out_path,paste("method_mod_obj_",interpolation_method,"_",y_var_name,out_prefix,".RData",sep="")))
    
  }
  
  if (interpolation_method=="gwr_daily"){
    #input a list:note that ghcn.subsets is not sampling_obj$data_day_ghcn
    i<-1
    list_param_run_prediction_gwr_daily <-list(i,s_raster,covar_names,lst_avg,list_models,dst,var,y_var_name, sampling_obj,interpolation_method,out_prefix,out_path)
    names(list_param_run_prediction_gwr_daily)<-c("list_index","covar_rast","covar_names","lst_avg","list_models","dst","var","y_var_name","sampling_obj","interpolation_method","out_prefix","out_path")
    #test <- run_interp_day_fun(1,list_param_run_prediction_gwr_daily)
    method_mod_obj<-mclapply(1:length(sampling_obj$ghcn_data),list_param=list_param_run_prediction_gwr_daily,run_interp_day_fun,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
    #method_mod_obj<-mclapply(1:22,list_param=list_param_run_prediction_gwr_daily,run_interp_day_fun,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
    #method_mod_obj<-mclapply(1:18,list_param=list_param_run_prediction_kriging_daily,runKriging_day_fun,mc.preschedule=FALSE,mc.cores = 9) #This is the end bracket from mclapply(...) statement
    
    save(method_mod_obj,file= file.path(out_path,paste("method_mod_obj_",interpolation_method,"_",y_var_name,out_prefix,".RData",sep="")))
    
  }
  t2<-proc.time()-t1
  cat(as.character(t2),file=log_fname,sep="\n", append=TRUE)
  #browser()
  
  ############### NOW RUN VALIDATION #########################
  #SIMPLIFY THIS PART: one call
  cat("Validation step:",file=log_fname,sep="\n", append=TRUE)
  t1<-proc.time()
  cat(paste("Local Date and Time: ",as.character(Sys.time()),sep=""),
      file=log_fname,sep="\n")
  
  if (interpolation_method=="gam_daily" | interpolation_method=="kriging_daily" | interpolation_method=="gwr_daily"){
    multi_time_scale <- FALSE
    
    list_data_v <- extract_list_from_list_obj(method_mod_obj,"data_v")
    list_data_s <- extract_list_from_list_obj(method_mod_obj,"data_s")
    rast_day_yearlist <- extract_list_from_list_obj(method_mod_obj,y_var_name) #list_tmp #list of predicted images over full year...
    list_sampling_dat <- extract_list_from_list_obj(method_mod_obj,"sampling_dat")
    
    list_param_validation<-list(i,rast_day_yearlist,list_data_v,list_data_s,list_sampling_dat,y_var_name, multi_time_scale,out_prefix, out_path)
    names(list_param_validation)<-c("list_index","rast_day_year_list",
                                  "list_data_v","list_data_s","list_sampling_dat","y_ref","multi_time_scale","out_prefix", "out_path") #same names for any method
    #debug(calculate_accuracy_metrics)
    #test_val2 <-calculate_accuracy_metrics(1,list_param_validation)
  
    validation_mod_obj <-mclapply(1:length(method_mod_obj), list_param=list_param_validation, calculate_accuracy_metrics,mc.preschedule=FALSE,mc.cores = num_cores) 
    save(validation_mod_obj,file= file.path(out_path,paste(interpolation_method,"_validation_mod_obj_",y_var_name,out_prefix,".RData",sep="")))
    t2<-proc.time()-t1
    cat(as.character(t2),file=log_fname,sep="\n", append=TRUE)
  }
    
  ### Run monthly validation if multi-time scale methods and add information to daily...
  if (interpolation_method %in% c("gam_CAI","kriging_CAI","gwr_CAI","gam_fusion","kriging_fusion","gwr_fusion")){
    multi_time_scale <- TRUE
    i<-1
    
    ## daily time scale
    list_data_v <- extract_list_from_list_obj(method_mod_obj,"data_v")
    list_data_s <- extract_list_from_list_obj(method_mod_obj,"data_s")
    rast_day_yearlist <- extract_list_from_list_obj(method_mod_obj,y_var_name) #list_tmp #list of predicted images over full year...
    list_sampling_dat <- extract_list_from_list_obj(method_mod_obj,"daily_dev_sampling_dat")
    
    list_param_validation<-list(i,rast_day_yearlist,list_data_v,list_data_s,list_sampling_dat,y_var_name, multi_time_scale,out_prefix, out_path)
    names(list_param_validation)<-c("list_index","rast_day_year_list",
                                    "list_data_v","list_data_s","list_sampling_dat","y_ref","multi_time_scale","out_prefix", "out_path") #same names for any method
    #debug(calculate_accuracy_metrics)
    #test_val2 <-calculate_accuracy_metrics(1,list_param_validation)
    
    validation_mod_obj <-mclapply(1:length(method_mod_obj), list_param=list_param_validation, calculate_accuracy_metrics,mc.preschedule=FALSE,mc.cores = num_cores) 
    save(validation_mod_obj,file= file.path(out_path,paste(interpolation_method,"_validation_mod_obj_",y_var_name,out_prefix,".RData",sep="")))
    
    ### monthly time scale
    list_data_v <- extract_list_from_list_obj(clim_method_mod_obj,"data_month_v") #extract monthly testing/validation dataset
    list_data_s <- extract_list_from_list_obj(clim_method_mod_obj,"data_month") #extract monthly training/fitting dataset
    rast_day_yearlist <- extract_list_from_list_obj(clim_method_mod_obj,"clim") #list_tmp #list of predicted images over full year at monthly time scale
    list_sampling_dat <- extract_list_from_list_obj(clim_method_mod_obj,"sampling_month_dat")
    
    #list_param_validation_month <-list(i,clim_yearlist,clim_method_mod_obj,y_var_name, multi_time_scale ,out_prefix, out_path)
    #names(list_param_validation_month)<-c("list_index","rast_day_year_list","method_mod_obj","y_var_name","multi_time_scale","out_prefix", "out_path") #same names for any method
    
    list_param_validation_month <-list(i,rast_day_yearlist,list_data_v,list_data_s,list_sampling_dat,y_var_month, multi_time_scale,out_prefix, out_path)
    names(list_param_validation_month)<-c("list_index","rast_day_year_list",
                                    "list_data_v","list_data_s","list_sampling_dat","y_ref","multi_time_scale","out_prefix", "out_path") #same names for any method
    #debug(calculate_accuracy_metrics)    
    #test_val2 <-calculate_accuracy_metrics(1,list_param_validation_month)
    
    validation_mod_month_obj <- mclapply(1:length(clim_method_mod_obj), list_param=list_param_validation_month, calculate_accuracy_metrics,mc.preschedule=FALSE,mc.cores = num_cores) 
    #test_val<-calculate_accuracy_metrics(1,list_param_validation)
    save(validation_mod_month_obj,file= file.path(out_path,paste(interpolation_method,"_validation_mod_month_obj_",y_var_name,out_prefix,".RData",sep="")))
  
    ##Create data.frame with validation and fit metrics for a full year/full numbe of runs
    tb_month_diagnostic_v <- extract_from_list_obj(validation_mod_month_obj,"metrics_v") 
    #tb_diagnostic_v contains accuracy metrics for models sample and proportion for every run...if full year then 365 rows maximum
    rownames(tb_month_diagnostic_v)<-NULL #remove row names
    tb_month_diagnostic_v$method_interp <- interpolation_method
    tb_month_diagnostic_s<-extract_from_list_obj(validation_mod_month_obj,"metrics_s")
    rownames(tb_month_diagnostic_s)<-NULL #remove row names
    tb_month_diagnostic_s$method_interp <- interpolation_method #add type of interpolation...out_prefix too??
    
  }

  #Cleaning raster temp files
  removeTmpFiles(h=0)

  #################### ASSESSMENT OF PREDICTIONS: PLOTS OF ACCURACY METRICS ###########
  ##Create data.frame with validation and fit metrics for a full year/full numbe of runs
  tb_diagnostic_v<-extract_from_list_obj(validation_mod_obj,"metrics_v") 
  #tb_diagnostic_v contains accuracy metrics for models sample and proportion for every run...if full year then 365 rows maximum
  rownames(tb_diagnostic_v)<-NULL #remove row names
  tb_diagnostic_v$method_interp <- interpolation_method
  tb_diagnostic_s<-extract_from_list_obj(validation_mod_obj,"metrics_s")
  rownames(tb_diagnostic_s)<-NULL #remove row names
  tb_diagnostic_s$method_interp <- interpolation_method #add type of interpolation...out_prefix too??
  
  #Call functions to create plots of metrics for validation dataset
  metric_names<-c("rmse","mae","me","r","m50")
  summary_metrics_v<- boxplot_from_tb(tb_diagnostic_v,metric_names,out_prefix,out_path) #if adding for fit need to change outprefix
  names(summary_metrics_v)<-c("avg","median")
  summary_month_metrics_v<- boxplot_month_from_tb(tb_diagnostic_v,metric_names,out_prefix,out_path)
  
  #################### CLOSE LOG FILE  ####################
  
  #close log_file connection and add meta data
  cat("Finished script process time:",file=log_fname,sep="\n", append=TRUE)
  time2<-proc.time()-time1
  cat(as.character(time2),file=log_fname,sep="\n", append=TRUE)
  #later on add all the parameters used in the script...
  cat(paste("Finished script at this local Date and Time: ",as.character(Sys.time()),sep=""),
             file=log_fname,sep="\n", append=TRUE)
  cat("End of script",file=log_fname,sep="\n", append=TRUE)
  #close(log_fname)
  
  ################### PREPARE RETURN OBJECT ###############
  #Will add more information to be returned
  if (interpolation_method %in% c("gam_CAI","kriging_CAI","gwr_CAI","gam_fusion","kriging_fusion","gwr_fusion")){
    raster_prediction_obj<-list(clim_method_mod_obj,method_mod_obj,validation_mod_obj,validation_mod_month_obj, tb_diagnostic_v,
                                tb_diagnostic_s,tb_month_diagnostic_v,tb_month_diagnostic_s,summary_metrics_v,summary_month_metrics_v)
    names(raster_prediction_obj)<-c("clim_method_mod_obj","method_mod_obj","validation_mod_obj","validation_mod_month_obj","tb_diagnostic_v",
                                    "tb_diagnostic_s","tb_month_diagnostic_v","tb_month_diagnostic_s","summary_metrics_v","summary_month_metrics_v")  
    save(raster_prediction_obj,file= file.path(out_path,paste("raster_prediction_obj_",interpolation_method,"_", y_var_name,out_prefix,".RData",sep="")))
    
  }
  
  #use %in% instead of "|" operator
  if (interpolation_method=="gam_daily" | interpolation_method=="kriging_daily" | interpolation_method=="gwr_daily"){
    raster_prediction_obj<-list(method_mod_obj,validation_mod_obj,tb_diagnostic_v,
                                tb_diagnostic_s,summary_metrics_v,summary_month_metrics_v)
    names(raster_prediction_obj)<-c("method_mod_obj","validation_mod_obj","tb_diagnostic_v",
                                    "tb_diagnostic_s","summary_metrics_v","summary_month_metrics_v")  
    save(raster_prediction_obj,file= file.path(out_path,paste("raster_prediction_obj_",interpolation_method,"_", y_var_name,out_prefix,".RData",sep="")))
    
  }

  return(raster_prediction_obj)
}

####################################################################
######################## END OF SCRIPT/FUNCTION #####################
