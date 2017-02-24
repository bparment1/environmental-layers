####################################  Cimate Layers Predctions  #######################################
############################  Script for assessment of scaling up on NEX: part 1a ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#The purpose is to create as set of functions to diagnose and assess quickly a set of predictd tiles.
#Part 1 create summary tables and inputs files for figure in part 2 and part 3.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 03/23/2014  
#MODIFIED ON: 01/04/2016            
#Version: 5
#PROJECT: Environmental Layers project  
#TO DO:
# - 
# - Separate call in a master script for assessment
# - add second stage in the master script for assessment
# - add mosaicing in the master script for assessment
# - clean up the code by making two function to clarify the code and remove repetition 

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh  
#MODULEPATH=$MODULEPATH:/nex/modules/files
#module load pythonkits/gdal_1.10.0_python_2.7.3_nex

# These are the names and number for the current subset regions used for global runs:
#reg1 - North America (NAM)
#reg2 - Europe (WE)
#reg3 - Asia 
#reg4 - South America (SAM)
#reg5 - Africa (AF)
#reg6 - East Asia and Australia 

###################################################################################################


run_assessment_prediction_fun <-function(i,list_param_run_assessment_prediction){
  #This function assesses results from prediction of climate variables. 
  #Predictions are run on the NEX-NASA computer by tiles for a given year.
  #The function collects information by region and tiles and generate tables of accuracy by tiles/regions.
  #There are currently five processing:
  #reg1 (North Am), reg2(Europe),reg3(Asia), reg4 (South Am), reg5 (Africa), reg6 (Australia-Asia)
  ##There are 20 inputs to this function.
  #
  ##INPUTS:
  #1) in_dir1 : location of root directory containging tiles
  #2) region_names : region_names #e.g. c("reg23","reg4") 
  #3) y_var_name : list_param_run_assessment_prediction$y_var_name # e.g. dailyTmax" #PARAM3
  #4) interpolation_method :  e.g. #c("gam_CAI") #PARAM4
  #5) out_prefix : e.g. "run_global_analyses_pred_12282015" #PARAM5
  #6) out_dir <- list_param_run_assessment_prediction$out_dir #<- "/nobackupp8/bparmen1/" #PARAM6
  #7) create_out_dir_param: if true a new dirrectory  is craeted for the outputs 
  #8) proj_str: coordinates  CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84, #PARAM8
  #9) list_year_predicted : list of years predicted e.g. 1984:2004
  #10) file_format : format for mosaiced files #PARAM10
  #11) NA_flag_val : No data value, #PARAM11
  #12) num_cores : number of cores used #PARAM13
  #13) plotting_figures: if TRUE, figures are generated from tables using assessment part2 script 
  ##Parameters related to assessment part 2: plot functions
  #14) mosaic_plot <- FALSE #PARAM6
  #15) day_to_mosaic <- c("19920101","19920102","19920103") #PARAM7,#if daily mosaics NULL then mosaic all days of the year
  #16) multiple_region <- TRUE #PARAM 12
  #17) countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #PARAM 13, copy this on NEX too
  #18) plot_region <- TRUE
  #20) threshold_missing_day <- c(367,365,300,200) #PARM18
  ##OUTPUTS
  #1)
  #2)
  #3)
  #
  #
  
  ###Loading R library and packages   
  
  ### Loading R library and packages        
  #library used in the workflow production:
  library(gtools)                              # loading some useful tools 
  library(mgcv)                                # GAM package by Simon Wood
  library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
  library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
  library(rgdal)                               # GDAL wrapper for R, spatial utilities
  library(gstat)                               # Kriging and co-kriging by Pebesma et al.
  library(fields)                              # NCAR Spatial Interpolation methods such as kriging, splines
  library(raster)                              # Hijmans et al. package for raster processing
  library(gdata)                               # various tools with xls reading, cbindX
  library(rasterVis)                           # Raster plotting functions
  library(parallel)                            # Parallelization of processes with multiple cores
  library(maptools)                            # Tools and functions for sp and other spatial objects e.g. spCbind
  library(maps)                                # Tools and data for spatial/geographic objects
  library(reshape)                             # Change shape of object, summarize results 
  library(plotrix)                             # Additional plotting functions
  library(plyr)                                # Various tools including rbind.fill
  library(spgwr)                               # GWR method
  library(automap)                             # Kriging automatic fitting of variogram using gstat
  library(rgeos)                               # Geometric, topologic library of functions
  #RPostgreSQL                                 # Interface R and Postgres, not used in this script
  library(gridExtra)
  #Additional libraries not used in workflow
  library(pgirmess)                            # Krusall Wallis test with mulitple options, Kruskalmc {pgirmess}  
  library(colorRamps)
  
  #### FUNCTION USED IN SCRIPT
  
  #This is read in in the master script
  #function_analyses_paper1 <- "global_run_scalingup_assessment_part1_functions_02112015.R" #PARAM12
  #script_path <- "/nobackupp8/bparmen1/env_layers_scripts" #path to script
  #source(file.path(script_path,function_analyses_paper1)) #source all functions used in this script 
  
  ##############################
  #### Parameters and constants  
  

  in_dir1 <- list_param_run_assessment_prediction$in_dir1 
  region_name <- list_param_run_assessment_prediction$region_name #e.g. c("reg23","reg4") #run only for one region
  y_var_name <- list_param_run_assessment_prediction$y_var_name # e.g. dailyTmax" #PARAM3
  interpolation_method <- list_param_run_assessment_prediction$interpolation_method #c("gam_CAI") #PARAM4
  out_prefix <- list_param_run_assessment_prediction$out_prefix #output suffix e.g."run_global_analyses_pred_12282015" #PARAM5
  out_dir <- list_param_run_assessment_prediction$out_dir #<- "/nobackupp8/bparmen1/" #PARAM6
  create_out_dir_param <-list_param_run_assessment_prediction$create_out_dir_param #if TRUE output dir created #PARAM7
  proj_str <- list_param_run_assessment_prediction$proj_str # CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84, #PARAM8
  list_year_predicted <- list_param_run_assessment_prediction$list_year_predicted # 1984:2004
  file_format <- list_param_run_assessment_prediction$file_format #<- ".tif" #format for mosaiced files #PARAM10
  NA_flag_val <- list_param_run_assessment_prediction$NA_flag_val #<- -9999  #No data value, #PARAM11
  num_cores <- list_param_run_assessment_prediction$num_cores #<- 6 #number of cores used #PARAM13
  plotting_figures <- list_param_run_assessment_prediction$plotting_figures #if true run part2 of assessment
  
  ##for plotting assessment function
  
  mosaic_plot <- list_param_run_assessment_prediction$mosaic_plot  #PARAM14
  day_to_mosaic <- list_param_run_assessment_prediction$day_to_mosaic #PARAM15
  multiple_region <- list_param_run_assessment_prediction$multiple_region #PARAM16
  countries_shp <- list_param_run_assessment_prediction$countries_shp #PARAM17
  plot_region <- list_param_run_assessment_prediction$plot_region #PARAM18
  threshold_missing_day <- list_param_run_assessment_prediction$threshold_missing_day #PARM20

  ########################## START SCRIPT #########################################
  
  #Need to make this a function to run as a job...
  
  ######################## PART0: Read content of predictions first.... #####
  #function looped over i, correspoding to year predicted
  
  year_predicted <- list_param_run_assessment_prediction$list_year_predicted[i] 
  #region_name is not null then restrict the assessment to a specific region
  #if(!is.null(region_name)){
  #  in_dir1 <- file.path(in_dir1,region_name)
  #}
  in_dir1 <- file.path(in_dir1,region_name)
  
  list_outfiles <- vector("list", length=14) #collect names of output files
  
  in_dir_list <- list.dirs(path=in_dir1,recursive=FALSE) #get the list regions processed for this run
  #basename(in_dir_list)
  #                       y=in_dir_list) 
  
  #in_dir_list_all  <- unlist(lapply(in_dir_list,function(x){list.dirs(path=x,recursive=F)}))
  in_dir_list_all <- in_dir_list
  #in_dir_list <- in_dir_list_all
  #in_dir_list <- in_dir_list[grep("bak",basename(basename(in_dir_list)),invert=TRUE)] #the first one is the in_dir1
  
  #this was changed on 10052015 because the shapefiles were not matching!!!
  #in_dir_list_tmp <- list.dirs(path=in_dir1,recursive=FALSE) #get the list regions processed for this run
  #in_dir_subset <- in_dir_list_tmp[grep("subset",basename(in_dir_list_tmp),invert=FALSE)] #select directory with shapefiles...
  #in_dir_shp <- file.path(in_dir_subset,"shapefiles")
  #in_dir_list_tmp <- list.dirs(path=in_dir1,recursive=FALSE) #get the list regions processed for this run
  in_dir_subset <- in_dir_list_all[grep("subset",basename(in_dir_list_all),invert=FALSE)] #select directory with shapefiles...
  in_dir_shp <- file.path(in_dir_subset,"shapefiles")
  
  
  #select only directories used for predictions
  #nested structure, we need to go to higher level to obtain the tiles...
  in_dir_reg <- in_dir_list[grep(".*._.*.",basename(in_dir_list),invert=FALSE)] #select directory with shapefiles...
  #in_dir_reg <- in_dir_list[grep("july_tiffs",basename(in_dir_reg),invert=TRUE)] #select directory with shapefiles...
  in_dir_list <- in_dir_reg
  
  
  in_dir_list <- in_dir_list[grep("bak",basename(basename(in_dir_list)),invert=TRUE)] #the first one is the in_dir1
  #list of shapefiles used to define tiles
  in_dir_shp_list <- list.files(in_dir_shp,".shp",full.names=T)
  
  ## load problematic tiles or additional runs
  #modify later...
  
  #system("ls /nobackup/bparmen1")
  
  if(create_out_dir_param==TRUE){
    out_dir <- create_dir_fun(out_dir,out_prefix)
    setwd(out_dir)
  }else{
    setwd(out_dir) #use previoulsy defined directory
  }
  
  setwd(out_dir)
  
  ##raster_prediction object : contains testing and training stations with RMSE and model object
  in_dir_list_tmp <- file.path(in_dir_list,year_predicted)
  list_raster_obj_files <- try(lapply(in_dir_list_tmp,FUN=function(x){list.files(path=x,pattern="^raster_prediction_obj.*.RData",full.names=T)}))
  #Add stop message here...if no raster object in any tiles then break from the function
  
  list_names_tile_coord <- lapply(list_raster_obj_files,FUN=function(x){basename(dirname(x))})
  list_names_tile_id <- paste("tile",1:length(list_raster_obj_files),sep="_")
  names(list_raster_obj_files)<- list_names_tile_id
  
  #one level up
  lf_covar_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar_obj.*.RData",full.names=T)})
  lf_covar_tif <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar.*.tif",full.names=T)})
  
  lf_sub_sampling_obj_files <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern=paste("^sub_sampling_obj_",interpolation_method,".*.RData",sep=""),full.names=T)})
  lf_sub_sampling_obj_files <- lapply(in_dir_list_tmp,FUN=function(x){list.files(path=x,pattern=paste("^sub_sampling_obj_",interpolation_method,".*.RData",sep=""),full.names=T)})
  
  
  lf_sub_sampling_obj_daily_files <- lapply(in_dir_list_tmp,FUN=function(x){list.files(path=x,pattern="^sub_sampling_obj_daily.*.RData",full.names=T)})
  
  ################################################################
  ######## PART 1: Generate tables to collect information:
  ######## over all tiles in North America 
  
  ##Function to collect all the tables from tiles into a table
  ###Table 1: Average accuracy metrics
  ###Table 2: daily accuracy metrics for all tiles
  
  #First create table of tiles under analysis and their coord
  df_tile_processed <- data.frame(tile_coord=basename(in_dir_list))
  df_tile_processed$tile_id <- unlist(list_names_tile_id) #Arbitrary tiling number!!
  df_tile_processed$path_NEX <- in_dir_list
  df_tile_processed$year_predicted <- year_predicted
  #Deal with the abscence of subsampling object for specific tiles
  lf_sub_sampling_obj_files_tmp <- lapply(1:length(lf_sub_sampling_obj_files),FUN=function(i,x){val <- x[[i]];if(length(val)==0){val<-0};val},x=lf_sub_sampling_obj_files)
  lf_sub_sampling_obj_daily_files_tmp <- lapply(1:length(lf_sub_sampling_obj_daily_files),FUN=function(i,x){val <- x[[i]];if(length(val)==0){val<-0};val},x=lf_sub_sampling_obj_daily_files)
  df_tile_processed$sub_sampling_clim  <- unlist(lf_sub_sampling_obj_files_tmp)
  df_tile_processed$sub_sampling_daily  <- unlist(lf_sub_sampling_obj_daily_files_tmp)
  ##review this part!!
  df_tile_processed$reg  <- region_name #eg c("reg4) #should only be one char string  
  df_tile_processed$year_predicted  <- year_predicted #eg c("reg4) #should only be one char string    
  #lf_sub_sampling_obj_files
 
  ################
  #### Table 1: Average accuracy metrics per tile and predictions
  
  #can use a maximum of 6 cores on the NEX Bridge
  #For 28 tiles and  28 RData boject it takes 15-16 min
  #summary_metrics_v_list <- mclapply(list_raster_obj_files[5:6],FUN=function(x){try( x<- load_obj(x)); try(x[["summary_metrics_v"]]$avg)},mc.preschedule=FALSE,mc.cores = 2)                           
  
  summary_metrics_v_list <- mclapply(list_raster_obj_files,FUN=function(x){try( x<- load_obj(x)); try(x[["summary_metrics_v"]]$avg)},mc.preschedule=FALSE,mc.cores = num_cores)                         
  #summary_metrics_v_list <- lapply(summary_metrics_v_list,FUN=function(x){try(x$avg)})
  names(summary_metrics_v_list) <- list_names_tile_id
  
  summary_metrics_v_tmp <- remove_from_list_fun(summary_metrics_v_list)$list
  df_tile_processed$metrics_v <- as.integer(remove_from_list_fun(summary_metrics_v_list)$valid)
  #Now remove "try-error" from list of accuracy)
  
  summary_metrics_v_NA <- do.call(rbind.fill,summary_metrics_v_tmp) #create a df for NA tiles with all accuracy metrics
  #tile_coord <- lapply(1:length(summary_metrics_v_list),FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=summary_metrics_v_list)
  #add the tile id identifier
  tile_id_tmp <- lapply(1:length(summary_metrics_v_tmp),
                        FUN=function(i,x,y){rep(y[i],nrow(x[[i]]))},x=summary_metrics_v_tmp,y=names(summary_metrics_v_tmp))
  #adding tile id summary data.frame
  summary_metrics_v_NA$tile_id <-unlist(tile_id_tmp)
  summary_metrics_v_NA$n <- as.integer(summary_metrics_v_NA$n)
  
  summary_metrics_v_NA <- merge(summary_metrics_v_NA,df_tile_processed[,1:2],by="tile_id")
  
  tx<-strsplit(as.character(summary_metrics_v_NA$tile_coord),"_")
  lat<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][1]},x=tx))
  long<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][2]},x=tx))
  summary_metrics_v_NA$lat <- lat
  summary_metrics_v_NA$lon <- long
  summary_metrics_v_NA$reg <- region_name #add region name
  summary_metrics_v_NA$year_predicted <- year_predicted #add year predicted
  
  write.table(as.data.frame(summary_metrics_v_NA),
              file=file.path(out_dir,paste("summary_metrics_v2_NA_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
  
  list_outfiles[[1]] <- file.path(out_dir,paste("summary_metrics_v2_NA_",year_predicted,"_",out_prefix,".txt",sep=""))
    
  #################
  ###Table 2: daily validation/testing accuracy metrics for all tiles
  #this takes about 15min for 28 tiles (reg4)
  #tb_diagnostic_v_list <- lapply(list_raster_obj_files,FUN=function(x){x<-load_obj(x);x[["tb_diagnostic_v"]]})                           
  tb_diagnostic_v_list <- mclapply(list_raster_obj_files,FUN=function(x){try(x<-load_obj(x));try(x[["tb_diagnostic_v"]])},mc.preschedule=FALSE,mc.cores = num_cores)                           
  
  names(tb_diagnostic_v_list) <- list_names_tile_id
  tb_diagnostic_v_tmp <- remove_from_list_fun(tb_diagnostic_v_list)$list
  #df_tile_processed$tb_diag <- remove_from_list_fun(tb_diagnostic_v_list)$valid
  
  tb_diagnostic_v_NA <- do.call(rbind.fill,tb_diagnostic_v_tmp) #create a df for NA tiles with all accuracy metrics
  tile_id_tmp <- lapply(1:length(tb_diagnostic_v_tmp),
                        FUN=function(i,x,y){rep(y[i],nrow(x[[i]]))},x=tb_diagnostic_v_tmp,y=names(tb_diagnostic_v_tmp))
  
  tb_diagnostic_v_NA$tile_id <- unlist(tile_id_tmp) #adding identifier for tile
  tb_diagnostic_v_NA <- merge(tb_diagnostic_v_NA,df_tile_processed[,1:2],by="tile_id")
  tb_diagnostic_v_NA$reg <- region_name #add region name
  tb_diagnostic_v_NA$year_predicted <- year_predicted #add year
  
  write.table((tb_diagnostic_v_NA),
              file=file.path(out_dir,paste("tb_diagnostic_v_NA_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
 
  list_outfiles[[2]] <- file.path(out_dir,paste("tb_diagnostic_v_NA_",year_predicted,"_",out_prefix,".txt",sep=""))
 
  #################
  ###Table 3: monthly fit/training accuracy information for all tiles
  
  ## Monthly fitting information
  tb_month_diagnostic_s_list <- mclapply(list_raster_obj_files,FUN=function(x){try(x<-load_obj(x));try(x[["tb_month_diagnostic_s"]])},mc.preschedule=FALSE,mc.cores = num_cores)                           
  
  names(tb_month_diagnostic_s_list) <- list_names_tile_id
  tb_month_diagnostic_s_tmp <- remove_from_list_fun(tb_month_diagnostic_s_list)$list
  #df_tile_processed$tb_diag <- remove_from_list_fun(tb_diagnostic_v_list)$valid
  
  tb_month_diagnostic_s_NA <- do.call(rbind.fill,tb_month_diagnostic_s_tmp) #create a df for NA tiles with all accuracy metrics
  tile_id_tmp <- lapply(1:length(tb_month_diagnostic_s_tmp),
                        FUN=function(i,x,y){rep(y[i],nrow(x[[i]]))},x=tb_month_diagnostic_s_tmp,y=names(tb_month_diagnostic_s_tmp))
  
  tb_month_diagnostic_s_NA$tile_id <- unlist(tile_id_tmp) #adding identifier for tile
  
  tb_month_diagnostic_s_NA <- merge(tb_month_diagnostic_s_NA,df_tile_processed[,1:2],by="tile_id")
  
  date_f<-strptime(tb_month_diagnostic_s_NA$date, "%Y%m%d")   # interpolation date being processed
  tb_month_diagnostic_s_NA$month <-strftime(date_f, "%m")          # current month of the date being processed
  tb_month_diagnostic_s_NA$reg <- region_name #add region name
  tb_month_diagnostic_s_NA$year_predicted <- year_predicted #add year

  write.table((tb_month_diagnostic_s_NA),
              file=file.path(out_dir,paste("tb_month_diagnostic_s_NA_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")

  list_outfiles[[3]] <- file.path(out_dir,paste("tb_month_diagnostic_s_NA_",year_predicted,"_",out_prefix,".txt",sep=""))
  
  #################
  ###Table 4: daily fit/training accuracy information with predictions for all tiles
  
  ## daily fit info:
  
  tb_diagnostic_s_list <- mclapply(list_raster_obj_files,FUN=function(x){try(x<-load_obj(x));try(x[["tb_diagnostic_s"]])},mc.preschedule=FALSE,mc.cores = num_cores)                           
  
  names(tb_diagnostic_s_list) <- list_names_tile_id
  tb_diagnostic_s_tmp <- remove_from_list_fun(tb_diagnostic_s_list)$list
  #df_tile_processed$tb_diag <- remove_from_list_fun(tb_diagnostic_v_list)$valid
  
  tb_diagnostic_s_NA <- do.call(rbind.fill,tb_diagnostic_s_tmp) #create a df for NA tiles with all accuracy metrics
  tile_id_tmp <- lapply(1:length(tb_diagnostic_s_tmp),
                        FUN=function(i,x,y){rep(y[i],nrow(x[[i]]))},x=tb_diagnostic_s_tmp,y=names(tb_diagnostic_s_tmp))
  
  tb_diagnostic_s_NA$tile_id <- unlist(tile_id_tmp) #adding identifier for tile
  
  tb_diagnostic_s_NA <- merge(tb_diagnostic_s_NA,df_tile_processed[,1:2],by="tile_id")
  tb_diagnostic_s_NA$reg <- region_name #add region name
  tb_diagnostic_s_NA$year_predicted <- year_predicted #add year

  write.table((tb_diagnostic_s_NA),
              file=file.path(out_dir,paste("tb_diagnostic_s_NA_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
  list_outfiles[[4]] <- file.path(out_dir,paste("tb_diagnostic_s_NA_",year_predicted,"_",out_prefix,".txt",sep=""))
  
  ##### Table 5: monthly station information used in training

  ### Make this part a function...this is repetitive
  #load data_month for specific tiles
  #10.45pm
  #data_month <- extract_from_list_obj(robj1$clim_method_mod_obj,"data_month")
  #names(data_month) #this contains LST means (mm_1, mm_2 etc.) as well as TMax and other info
  
  #data_month_s_list <- mclapply(list_raster_obj_files,FUN=function(x){try(x<-load_obj(x));try(x$validation_mod_month_obj[["data_s"]])},mc.preschedule=FALSE,mc.cores = 6)                           
  data_month_s_list <- mclapply(list_raster_obj_files,FUN=function(x){try(x<-load_obj(x));try(extract_from_list_obj(x$validation_mod_month_obj,"data_s"))},mc.preschedule=FALSE,mc.cores = 6)                           
  #test <- mclapply(list_raster_obj_files[1:6],FUN=function(x){try(x<-load_obj(x));try(extract_from_list_obj(x$validation_mod_month_obj,"data_s"))},mc.preschedule=FALSE,mc.cores = 6)                           
  
  names(data_month_s_list) <- list_names_tile_id
  
  data_month_tmp <- remove_from_list_fun(data_month_s_list)$list
  #df_tile_processed$metrics_v <- remove_from_list_fun(data_month_s_list)$valid
  
  tile_id <- lapply(1:length(data_month_tmp),
                    FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=data_month_tmp)
  data_month_NAM <- do.call(rbind.fill,data_month_tmp) #combined data_month for "NAM" North America
  data_month_NAM$tile_id <- unlist(tile_id)
  data_month_NAM$reg <- region_name #add region name
  data_month_NAM$year_predicted <- year_predicted #add year

  write.table((data_month_NAM),
              file=file.path(out_dir,paste("data_month_s_NAM_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
  list_outfiles[[5]] <- file.path(out_dir,paste("data_month_s_NAM_",year_predicted,"_",out_prefix,".txt",sep=""))
  browser()
  
  ##### Table 6 and table 7: stations for daily predictions
  
  data_day_s_list <- mclapply(list_raster_obj_files,FUN=function(x){try(x<-load_obj(x));try(extract_from_list_obj(x$validation_mod_obj,"data_s"))},mc.preschedule=FALSE,mc.cores = num_cores)    
  data_day_v_list <- mclapply(list_raster_obj_files,FUN=function(x){try(x<-load_obj(x));try(extract_from_list_obj(x$validation_mod_obj,"data_v"))},mc.preschedule=FALSE,mc.cores = num_cores)    
  
  names(data_day_s_list) <- list_names_tile_id
  names(data_day_v_list) <- list_names_tile_id
  
  data_day_s_tmp <- remove_from_list_fun(data_day_s_list)$list
  data_day_v_tmp <- remove_from_list_fun(data_day_v_list)$list
  
  #df_tile_processed$metrics_v <- remove_from_list_fun(data_month_s_list)$valid
  
  tile_id <- lapply(1:length(data_day_s_tmp),
                    FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=data_day_s_tmp)
  data_day_s_NAM <- do.call(rbind.fill,data_day_s_tmp) #combined data_month for "NAM" North America
  data_day_s_NAM$tile_id <- unlist(tile_id)
  data_day_s_NAM$reg <- region_name #add region name
  data_day_s_NAM$year_predicted <- year_predicted #add year predicted
  
  tile_id <- lapply(1:length(data_day_v_tmp),
                    FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=data_day_v_tmp)
  data_day_v_NAM <- do.call(rbind.fill,data_day_v_tmp) #combined data_month for "NAM" North America
  data_day_v_NAM$tile_id <- unlist(tile_id)
  data_day_v_NAM$reg <- region_name #add region name
  data_day_v_NAM$year_predicted <- year_predicted #add year predicted
  
  write.table((data_day_s_NAM),
              file=file.path(out_dir,paste("data_day_s_NAM_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
  write.table((data_day_v_NAM),
              file=file.path(out_dir,paste("data_day_v_NAM_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
  list_outfiles[[6]] <- file.path(out_dir,paste("data_day_s_NAM_",year_predicted,"_",out_prefix,".txt",sep=""))
  list_outfiles[[7]] <- file.path(out_dir,paste("data_day_v_NAM_",year_predicted,"_",out_prefix,".txt",sep=""))
  
  ##### Table 8: validation stations for monthly predictions
  
  #### Recover subsampling data
  #For tiles with many stations, there is a subsampling done in terms of distance (spatial pruning) and 
  #in terms of station numbers if there are still too many stations to consider. This is done at the 
  #daily and monthly stages.
  
  #lf_sub_sampling_obj_files <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern=paste("^sub_sampling_obj_",interpolation_method,".*.RData",sep=""),full.names=T)})
  #lf_sub_sampling_obj_daily_files <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="^sub_sampling_obj_daily.*.RData",full.names=T)})
  #sub_sampling_obj <- try(load_obj(lf_sub_sampling_obj_files[[3]])) #This is an example tile
  #data_removed contains the validation data...
  #this data can be used for validation of the product. Note that it may be missing for some tiles
  #as no stations are removed if there are too enough stations in the tile
  #this will need to be checked later on...
  
  ### This is not working on 02/24/2017, check if data input is present
  data_month_v_subsampling_list <- mclapply(lf_sub_sampling_obj_files,FUN=function(x){try(x<-load_obj(x));try(extract_from_list_obj(x$validation_mod_month_obj,"data_removed"))},mc.preschedule=FALSE,mc.cores = 6)                           
  #test <- mclapply(list_raster_obj_files[1:6],FUN=function(x){try(x<-load_obj(x));try(extract_from_list_obj(x$validation_mod_month_obj,"data_s"))},mc.preschedule=FALSE,mc.cores = 6)                           
  names(data_month_v_subsampling_list) <- list_names_tile_id
  data_month_v_subsampling_tmp <- remove_from_list_fun(data_month_v_subsampling_list)$list
  #df_tile_processed$metrics_v <- remove_from_list_fun(data_month_s_list)$valid
  #if no stations have been removed then there are no validation stations !!!
  if(length(data_month_v_subsampling_tmp)!=0){
    tile_id <- lapply(1:length(data_month_v_subsampling_tmp),
                      FUN=function(i,x){try(rep(names(x)[i],nrow(x[[i]])))},x=data_month_v_subsampling_tmp)
    data_month_v_subsmapling_NAM <- do.call(rbind.fill,ddata_month_v_subsampling_tmp) #combined data_month for "NAM" North America
    data_month_v_subsampling_NAM$tile_id <- unlist(tile_id)
    data_month_v_subsampling_NAM$reg <- region_name
    data_month_v_subsampling_NAM$year_predicted <- year_predicted
    
    write.table((data_month_v_subsampling_NAM),
                file=file.path(out_dir,paste("data_month_v_subsampling_NAM_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
    list_outfiles[[8]] <- file.path(out_dir,paste("data_month_v_subsampling_NAM_",year_predicted,"_",out_prefix,".txt",sep=""))
  }else{
    list_outfiles[[8]] <- NA
  }
  
  ####
  data_day_v_subsampling_list <- mclapply(lf_sub_sampling_obj_daily_files,FUN=function(x){try(x<-load_obj(x));try(extract_from_list_obj(x$validation_mod_month_obj,"data_removed"))},mc.preschedule=FALSE,mc.cores = 6)                           
  data_day_v_subsampling_list <- mclapply(lf_sub_sampling_obj_daily_files[1:6],FUN=function(x){try(x<-load_obj(x));try(extract_from_list_obj(x$validation_mod_month_obj,"data_removed"))},mc.preschedule=FALSE,mc.cores = 6)                           

  lf_sub_sampling_obj_daily_files
  
  ##### Table 9: validation accuracy metrics for monthly predictions
  
  #Get validation data?? Find other object from within the dir 
  #Som region don't have validation data at monthly time scale.

  #### To be changed later...there is no validation data at this stage
  ## Monthly fitting information
  #tb_month_diagnostic_v_list <- mclapply(list_raster_obj_files,FUN=function(x){try(x<-load_obj(x));try(x[["tb_month_diagnostic_v"]])},mc.preschedule=FALSE,mc.cores = num_cores)                           
  #names(tb_month_diagnostic_v_list) <- list_names_tile_id
  #tb_month_diagnostic_v_tmp <- remove_from_list_fun(tb_month_diagnostic_v_list)$list
  #tb_month_diagnostic_v_NA <- do.call(rbind.fill,tb_month_diagnostic_v_tmp) #create a df for NA tiles with all accuracy metrics
  #tile_id_tmp <- lapply(1:length(tb_month_diagnostic_v_tmp),
  #                      FUN=function(i,x,y){rep(y[i],nrow(x[[i]]))},x=tb_month_diagnostic_v_tmp,y=names(tb_month_diagnostic_v_tmp))
  #tb_month_diagnostic_v_NA$tile_id <- unlist(tile_id_tmp) #adding identifier for tile
  #tb_month_diagnostic_v_NA <- merge(tb_month_diagnostic_v_NA,df_tile_processed[,1:2],by="tile_id")
  #date_f<-strptime(tb_month_diagnostic_v_NA$date, "%Y%m%d")   # interpolation date being processed
  #tb_month_diagnostic_v_NA$month<-strftime(date_f, "%m")          # current month of the date being processed
  #write.table((tb_month_diagnostic_v_NA),
  #            file=file.path(out_dir,paste("tb_month_diagnostic_v_NA_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")

  #list_outfiles[[9]] <- file.path(out_dir,paste("tb_month_diagnostic_v_NA_",year_predicted,"_",out_prefix,".txt",sep=""))
  list_outfiles[[9]] <- NA
  
  ######################################################
  ####### PART 3: EXAMINE STATIONS AND MODEL FITTING ###
  
  ##### Table 10 and Table 11: extracting accuracy information from daily and monthly predictions
  
  ### Stations and model fitting ###
  #summarize location and number of training and testing used by tiles
  #names(data_month) #this contains LST means (mm_1, mm_2 etc.) as well as TMax and other info
  
  use_day=TRUE
  use_month=TRUE
  
  list_names_tile_id <- df_tile_processed$tile_id
  list_raster_obj_files[list_names_tile_id]
  #list_names_tile_id <- c("tile_1","tile_2")
  list_param_training_testing_info <- list(list_raster_obj_files[list_names_tile_id],use_month,use_day,list_names_tile_id)
  names(list_param_training_testing_info) <- c("list_raster_obj_files","use_month","use_day","list_names_tile_id")
  
  list_param <- list_param_training_testing_info
  #debug(extract_daily_training_testing_info)
  #pred_data_info <- extract_daily_training_testing_info(1,list_param=list_param_training_testing_info)
  pred_data_info <- mclapply(1:length(list_raster_obj_files[list_names_tile_id]),FUN=extract_daily_training_testing_info,list_param=list_param_training_testing_info,mc.preschedule=FALSE,mc.cores = num_cores)

  pred_data_info_tmp <- remove_from_list_fun(pred_data_info)$list #remove data not predicted
  ##Add tile nanmes?? it is alreaready there
  #names(pred_data_info)<-list_names_tile_id
  pred_data_month_info <- do.call(rbind,lapply(pred_data_info_tmp,function(x){x$pred_data_month_info}))
  pred_data_day_info <- do.call(rbind,lapply(pred_data_info_tmp,function(x){x$pred_data_day_info}))
  
  pred_data_month_info$reg <- region_name
  pred_data_day_info$reg <- region_name
  pred_data_month_info$year_predicted <- year_predicted
  pred_data_day_info$year_predicted <- year_predicted
    
  #putput inforamtion in csv !!
  write.table(pred_data_month_info,
              file=file.path(out_dir,paste("pred_data_month_info_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
  write.table(pred_data_day_info,
              file=file.path(out_dir,paste("pred_data_day_info_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
  list_outfiles[[10]] <- file.path(out_dir,paste("pred_data_month_info_",year_predicted,"_",out_prefix,".txt",sep=""))
  list_outfiles[[11]] <- file.path(out_dir,paste("pred_data_day_info_",year_predicted,"_",out_prefix,".txt",sep=""))
  
  #browser() #debugging on 01/042016
  
  ######################################################
  ####### PART 4: Get shapefiles defining region tiling with centroids ###
  
  ##### Table 12, Table 13, Table 14: collect location of predictions from shapefiles

  #get shape files for the region being assessed:
  list_shp_world <- list.files(path=in_dir_shp,pattern=".*.shp",full.names=T)
  l_shp <- gsub(".shp","",basename(list_shp_world))
  l_shp <- sub("shp_","",l_shp)
  l_shp <- unlist(lapply(1:length(l_shp),
                         FUN=function(i){paste(strsplit(l_shp[i],"_")[[1]][1:2],collapse="_")}))
  
  df_tiles_all <- as.data.frame(as.character(unlist(list_shp_world)))
  df_tiles_all$tile_coord <- l_shp
  #names(df_tiles_all) <- "list_shp_world"
  names(df_tiles_all) <- c("shp_files","tile_coord")
  #add tiles id
  df_tiles_all$reg <- region_name
  df_tiles_all$year_predicted <- year_predicted
  
  matching_index <- match(basename(in_dir_list),l_shp)
  list_shp_reg_files <- list_shp_world[matching_index]
  df_tile_processed$shp_files <-list_shp_reg_files

  tx<-strsplit(as.character(df_tile_processed$tile_coord),"_")
  lat<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][1]},x=tx))
  long<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][2]},x=tx))
  df_tile_processed$lat <- lat
  df_tile_processed$lon <- long
  
  #put that list in the df_processed and also the centroids!!
  write.table(df_tile_processed,
              file=file.path(out_dir,paste("df_tile_processed_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
  
  write.table(df_tiles_all,
              file=file.path(out_dir,paste("df_tiles_all_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
  list_outfiles[[12]] <- file.path(out_dir,paste("df_tile_processed_",year_predicted,"_",out_prefix,".txt",sep=""))
  list_outfiles[[13]] <- file.path(out_dir,paste("df_tiles_all_",year_predicted,"_",out_prefix,".txt",sep=""))
  
  #Copy to local home directory on NAS-NEX
  dir.create(file.path(out_dir,"shapefiles"))
  #list_shp_world_dbf <- gsub(".shp",".*",basename(list_shp_world)) #remove shp extension, did not owork with file.copy
  #list_shp_world_dbf <- gsub(".shp",".dbf",basename(list_shp_world)) #remove shp extension
  #list_shp_world_prj <- gsub(".shp",".dbf",basename(list_shp_world)) #remove shp extension
  #list_shp_world_tmp <- file.path(dirname(list_shp_world), list_shp_world_tmp)
  list_shp_all <- list.files(path=in_dir_shp,full.names=T)#need all the files!!! including .dbf etc.
  file.copy(list_shp_all,file.path(out_dir,"shapefiles"))
  
  #save a list of all files...
  write.table(df_tiles_all,
              file=file.path(out_dir,"shapefiles",paste("df_tiles_all_",year_predicted,"_",out_prefix,".txt",sep="")),sep=",")
  list_outfiles[[14]] <- file.path(out_dir,"shapefiles",paste("df_tiles_all_",year_predicted,"_",out_prefix,".txt",sep=""))

  ##############################################################
  ############## Collect information from assessment ##########
  
  outfiles_names <- c("summary_metrics_v_names","tb_v_accuracy_name","tb_month_s_name","tb_s_accuracy_name", 
  "data_month_s_name","data_day_v_name","data_day_s_name","data_month_v_name", "tb_month_v_name",
  "pred_data_month_info_name","pred_data_day_info_name","df_tile_processed_name","df_tiles_all_name", 
  "df_tiles_all_name") 
  names(list_outfiles) <- outfiles_names
  
  #This data.frame contains all the files from the assessment
  df_assessment_files <- data.frame(filename=outfiles_names,files=unlist(list_outfiles),
                                    reg=region_name,year=year_predicted)
  ###Prepare files for copying back?
  df_assessment_files_name <- file.path(out_dir,paste("df_assessment_files_",region_name,"_",year_predicted,"_",out_prefix,".txt",sep=""))
  write.table(df_assessment_files,
              file=df_assessment_files_name,sep=",")

  #with reg4 prediction 2014, it took 1h36minutes to reach this point in the code.
  #this was processed using the bridge1 with 6 cores...
  
  ######################################################
  ####### PART 5: run plotting functions to produce figures

  #browser() #debugging on 01/042016
  #out_dir <- "/nobackupp8/bparmen1/output_run_global_analyses_pred_12282015"
  in_dir <- out_dir #PARAM 0
  #y_var_name <- "dailyTmax" #PARAM1 , already set
  #interpolation_method <- c("gam_CAI") #PARAM2, already set
  out_suffix <- out_prefix #PARAM3
  #out_dir <-  #PARAM4, already set
  create_out_dir_param <- FALSE #PARAM 5, already created and set
  #mosaic_plot <- FALSE #PARAM6
  #if daily mosaics NULL then mosaicas all days of the year
  #day_to_mosaic <- c("19920101","19920102","19920103") #PARAM7
  #CRS_locs_WGS84 already set
  proj_str<- CRS_locs_WGS84 #PARAM 8 #check this parameter
  #file_format <- ".rst" #PARAM 9, already set
  #NA_flag_val <- -9999 #PARAM 11, already set
  #multiple_region <- TRUE #PARAM 12
  #countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #PARAM 13, copy this on NEX too
  #plot_region <- TRUE
  #num_cores <- 6 #PARAM 14, already set
  #region_name <- c("reg4") #reference region to merge if necessary, if world all the regions are together #PARAM 16
  #use previous files produced in step 1a and stored in a data.frame
  #df_assessment_files_name <- "df_assessment_files_reg4_2014_run_global_analyses_pred_12282015.txt"# #PARAM 17, set in the script
  #df_assessment_files <- read.table(df_assessment_files_name,stringsAsFactors=F,sep=",")
  #threshold_missing_day <- c(367,365,300,200) #PARM18

  list_param_run_assessment_plotting <- list(in_dir,y_var_name, interpolation_method, out_suffix, 
                      out_dir, create_out_dir_param, mosaic_plot, proj_str, file_format, NA_flag_val,
                      multiple_region, countries_shp, plot_region, num_cores, 
                      region_name, df_assessment_files_name, threshold_missing_day,year_predicted) 

  names(list_param_run_assessment_plotting) <- c("in_dir","y_var_name","interpolation_method","out_suffix", 
                      "out_dir","create_out_dir_param","mosaic_plot","proj_str","file_format","NA_flag_val",
                      "multiple_region","countries_shp","plot_region","num_cores", 
                      "region_name","df_assessment_files_name","threshold_missing_day","year_predicted") 
  
  #function_assessment_part2 <- "global_run_scalingup_assessment_part2_01032016.R"
  #source(file.path(script_path,function_assessment_part2)) #source all functions used in this script 

  #debug(run_assessment_plotting_prediction_fun)
  df_assessment_figures_files <- run_assessment_plotting_prediction_fun(list_param_run_assessment_plotting) 
  
  ######################################################
  ##### Prepare objet to return ####

  assessment_obj <- list(df_assessment_files, df_assessment_figures_files)
  names(assessment_obj) <- c("df_assessment_files", "df_assessment_figures_files")
  ## Prepare list of files to return...
  return(assessment_obj)
}

##################### END OF SCRIPT ######################


