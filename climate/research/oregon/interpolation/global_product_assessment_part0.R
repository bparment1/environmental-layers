##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part 0 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This script checks the number of predictions by tiles and years.
#with the goal of predicting potential gaps or missing predictions in fugure mosaics by region.
#The general logic is to check the number of overlap by shapefile polyon tiles
#along with the predicitons for every day of the year (*.tif)
#Summary tables and data are also produced in the script.
#
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/27/2016  
#MODIFIED ON: 10/30/2016            
#Version: 1
#PROJECT: Environmental Layers project     
#COMMENTS: Major update of code with changes to the listing of tiles files 
#TODO:
#1) 
#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rm u:aguzman4:rwx /nobackupp6/aguzman4/climateLayers/LST_tempSpline/
#COMMIT: generating animation for region 4 for multiple years sequences

#################################################################################################


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
library(zoo)
library(xts)
library(lubridate)
library(mosaic)

###### Function used in the script #######
  
script_path <- "/nobackupp8/bparmen1/env_layers_scripts" #path to script, NEX
#script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts" #path to script

## NASA poster and paper related
#source(file.path(script_path,"NASA2016_conference_temperature_predictions_function_05032016b.R"))

#Mosaic related on NEX
#script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts"
function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_09282016.R" #Functions used to mosaic predicted tiles
function_mosaicing <-"global_run_scalingup_mosaicing_09282016.R" #main scripts for mosaicing predicted tiles

source(file.path(script_path,function_mosaicing)) #source all functions used in this script 
source(file.path(script_path,function_mosaicing_functions)) #source all functions used in this script 

#Assessment on NEX
function_assessment_part1_functions <- "global_run_scalingup_assessment_part1_functions_12282015.R" #PARAM12
function_assessment_part1a <-"global_run_scalingup_assessment_part1a_01042016.R"
function_assessment_part2 <- "global_run_scalingup_assessment_part2_02092016.R"
function_assessment_part2_functions <- "global_run_scalingup_assessment_part2_functions_01032016.R"
function_assessment_part3 <- "global_run_scalingup_assessment_part3_07292016.R"

source(file.path(script_path,function_assessment_part1_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part1a)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part3)) #source all functions used in this script 

#Product assessment
function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_09192016b.R"
source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 
function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_10222016.R"
source(file.path(script_path,function_product_assessment_part2_functions)) #source all functions used in this script 

###############################
####### Parameters, constants and arguments ###

CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #constant 1

var<-"TMAX" # variable being interpolated #param 1, arg 1


#interpolation_method<-c("gam_fusion") #other otpions to be added later
interpolation_method<-c("gam_CAI") #param 2
CRS_interp <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #param 3
#CRS_interp <-"+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs";

#out_region_name<-""
#list_models<-c("y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)") #param 4
metric_name <- "var_pred" #use RMSE if accuracy

item_no <- 13
day_start <- "2000101" #PARAM 12 arg 12
day_end <- "20001231" #PARAM 13 arg 13
#date_start <- day_start
#date_end <- day_end
date_start <- day_start
date_end <- day_end
#day_start <- "1984101" #PARAM 12 arg 12
#day_end <- "20141231" #PARAM 13 arg 13
day_to_mosaic_range <- NULL

in_dir <- "/nobackupp6/aguzman4/climateLayers/out/reg6/assessment"
#in_dir <- "/nobackupp8/bparmen1/climateLayers/out/reg6/assessment"
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/mosaics/mosaic" #predicted mosaic

region_name <- c("reg6") #param 6, arg 3
out_suffix <- "global_assessment_reg6_10232016"

create_out_dir_param <- TRUE #param 9, arg 6

out_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/assessment"

#run_figure_by_year <- TRUE # param 10, arg 7

file_format <- ".tif" #format for mosaiced files # param 11
NA_flag_val <- -32768  #No data value, # param 12

#num_cores <- 6 #number of cores used # param 13, arg 8
plotting_figures <- TRUE #running part2 of assessment to generate figures... # param 14
num_cores <- 11 #number of cores used # param 13, arg 8
#python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin" #PARAM 30
python_bin <- "/usr/bin" #PARAM 30

#NA_flag_val_mosaic <- -3399999901438340239948148078125514752.000
NA_flag_val_mosaic <- -32768
in_dir_list_filename <- NULL #if NULL, use the in_dir directory to search for info
countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #Atlas
lf_raster <- NULL #list of raster to consider
item_no <- 13

#On NEX
#contains all data from the run by Alberto
in_dir1 <- "/nobackupp6/aguzman4/climateLayers/out" #On NEX
#parent output dir for the current script analyes

y_var_name <- "dailyTmax" #PARAM1
interpolation_method <- c("gam_CAI") #PARAM2
out_suffix <- "predictions_assessment_reg6_10302016"
#out_suffix <- "output_run10_1000x3000_global_analyses_02102015"
#out_suffix <- "run10_1500x4500_global_analyses_pred_1992_10052015" #PARAM3
#out_dir <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_10052015" #PARAM4
#create_out_dir_param <- TRUE #PARAM 5
#mosaic_plot <- FALSE #PARAM6
#if daily mosaics NULL then mosaicas all days of the year
#day_to_mosaic <- c("19920101","19920102","19920103") #PARAM7
#CRS_WGS84 <-    CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84 #CONSTANT1
#CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
#proj_str<- CRS_WGS84 #PARAM 8 #check this parameter
#file_format <- ".rst" #PARAM 9
NA_value <- -9999 #PARAM10
NA_flag_val <- NA_value
#multiple_region <- TRUE #PARAM 12
region_name <- "reg6" #PARAM 13
countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #PARAM 13, copy this on NEX too
#plot_region <- TRUE
num_cores <- 6 #PARAM 14
#/nobackupp6/aguzman4/climateLayers/out/reg6/subset/shapefiles
list_year_predicted <- c(2000,2012,2013) #year still on disk for reg6
  

############################################
#### Parameters and constants  

##Add for precip later...
if (var == "TMAX") {
  y_var_name <- "dailyTmax"
  y_var_month <- "TMax"
}
if (var == "TMIN") {
  y_var_name <- "dailyTmin"
  y_var_month <- "TMin"
}

##Add for precip later...
if (var == "TMAX") {
  variable_name <- "maximum temperature"
}
if (var == "TMIN") {
  variable_name <- "minimum temperature"
}

i<-1

predictions_tiles_missing_fun <- function(list_param,i){
  
  #from script:
  #interpolation/global_run_scalingup_assessment_part1a.R

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
  
  #list_outfiles <- vector("list", length=35) #collect names of output files, this should be dynamic?
  #list_outfiles_names <- vector("list", length=35) #collect names of output files

  year_predicted <- list_param_run_assessment_prediction$list_year_predicted[i] 

  in_dir1_reg <- file.path(in_dir1,region_name)
  
  list_outfiles <- vector("list", length=14) #collect names of output files
  
  in_dir_list <- list.dirs(path=in_dir1_reg,recursive=FALSE) #get the list regions processed for this run
  #basename(in_dir_list)
  #                       y=in_dir_list) 
  
  #in_dir_list_all  <- unlist(lapply(in_dir_list,function(x){list.dirs(path=x,recursive=F)}))
  in_dir_list_all <- in_dir_list
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
  out_dir <- in_dir
  if(create_out_dir_param==TRUE){
    out_dir <- create_dir_fun(out_dir,out_suffix)
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
  lf_sub_sampling_obj_daily_files <- lapply(in_dir_list_tmp,FUN=function(x){list.files(path=x,pattern="^sub_sampling_obj_daily.*.RData",full.names=T)})

  if(is.null(day_to_mosaic_range)){
  #  start_date <- #first date
     start_date <- paste0(year_processed,"0101") #change this later!!
     end_date <-   paste0(year_processed,"1231") #change this later!!
     day_to_mosaic <- seq(as.Date(strptime(start_date,"%Y%m%d")), as.Date(strptime(end_date,"%Y%m%d")), 'day')
     day_to_mosaic <- format(day_to_mosaic,"%Y%m%d") #format back to the relevant date format for files
  }else{
    start_date <- day_to_mosaic_range[1]
    end_date <- day_to_mosaic_range[2]
    day_to_mosaic <- seq(as.Date(strptime(start_date,"%Y%m%d")), as.Date(strptime(end_date,"%Y%m%d")), 'day')
    day_to_mosaic <- format(day_to_mosaic,"%Y%m%d") #format back to the relevant date format for files
  }
  
  in_dir_tiles_tmp <- in_dir1 #
  #in_dir_tiles_tmp <- in_dir_reg
  
  ### Do this by tile!!!
  
  # Making listing of files faster with multicores use
  lf_mosaic <- mclapply(1:length(day_to_mosaic),FUN=function(i){
    searchStr = paste(in_dir_tiles_tmp,"/*/",year_processed,"/gam_CAI_dailyTmax_predicted_",pred_mod_name,"*",day_to_mosaic[i],"*.tif",sep="")
    Sys.glob(searchStr)},mc.preschedule=FALSE,mc.cores = num_cores)
  
  lf_mosaic <- mclapply(1:length(day_to_mosaic),FUN=function(i){
    searchStr = paste(in_dir_tiles_tmp,"/*/",year_processed,"/gam_CAI_dailyTmax_predicted_",pred_mod_name,"*",day_to_mosaic[i],"*.tif",sep="")
    Sys.glob(searchStr)},mc.preschedule=FALSE,mc.cores = num_cores)

  ##Run this on reg4 and reg5 after
  #Add report by year in text file?
  #Using specified values for parameters
  test_missing <- check_missing(lf=lf_mosaic[[1]], 
                              pattern_str=NULL,
                              in_dir=in_dir_mosaic,
                              date_start="1984101",
                              date_end="20141231",
                              item_no=13,
                              out_suffix=out_suffix,
                              num_cores=num_cores,
                              out_dir=".")

  df_time_series <- test_missing$df_time_series
  head(df_time_series)

  table(df_time_series$missing)
  table(df_time_series$year)
  
  ########################
  #### Step 2: examine overlap
  
  #combine polygon
  #http://gis.stackexchange.com/questions/155328/merging-multiple-spatialpolygondataframes-into-1-spdf-in-r

  #http://gis.stackexchange.com/questions/116388/count-overlapping-polygons-in-single-shape-file-with-r

  #### Use the predictions directory
  #By region
  #For each polygon/tile find polygon overlapping with count and ID (like list w)
  #for each polygon/tile and date find if there is a prediction using the tif (multiply number)
  #for each date of year report data in table.

  #go through table and hsow if there are missing data (no prediction) or report min predictions for tile set?
    
  #for each polygon find you overlap!!
  #plot number of overlap
  #for specific each find prediction...
  
  ########################
  #### Step 3: combine overlap information and number of predictions by day
  
  
  

  return()
}

############################ END OF SCRIPT ##################################
