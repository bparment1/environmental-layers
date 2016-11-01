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
#MODIFIED ON: 11/01/2016            
#Version: 1
#PROJECT: Environmental Layers project     
#COMMENTS: testing of files by tiles and combining listing 
#TODO:
#1) 
#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rm u:aguzman4:rwx /nobackupp6/aguzman4/climateLayers/LST_tempSpline/
#COMMIT: combining tif tiles and shapefiles to examine potential gaps

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
function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_10312016.R"
source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 
##Don't load part 1 and part2, mosaic package does not work on NEX
#function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_09192016b.R"
#source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 
#function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_10222016.R"
#source(file.path(script_path,function_product_assessment_part2_functions)) #source all functions used in this script 

###############################
####### Parameters, constants and arguments ###

CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #constant 1
var<-"TMAX" # variable being interpolated #param 1, arg 1
interpolation_method<-c("gam_CAI") #param 2
CRS_interp <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #param 3
#list_models<-c("y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)") #param 4
metric_name <- "var_pred" #use RMSE if accuracy
pred_mod_name <- "mod1"
item_no <- 9
day_start <- "2000101" #PARAM 12 arg 12
day_end <- "20001231" #PARAM 13 arg 13
#date_start <- day_start
#date_end <- day_end
date_start <- day_start
date_end <- day_end
#day_start <- "1984101" #PARAM 12 arg 12
#day_end <- "20141231" #PARAM 13 arg 13
day_to_mosaic_range <- NULL
#infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg6.tif"
infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg6.tif"

in_dir <- "/nobackupp6/aguzman4/climateLayers/out/reg6/assessment"
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/mosaics/mosaic" #predicted mosaic
region_name <- c("reg6") #param 6, arg 3
out_suffix <- "global_assessment_reg6_10232016"
create_out_dir_param <- TRUE #param 9, arg 6
out_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/assessment"
file_format <- ".tif" #format for mosaiced files # param 11
NA_flag_val <- -32768  #No data value, # param 12
plotting_figures <- TRUE #running part2 of assessment to generate figures... # param 14
num_cores <- 6 #number of cores used # param 13, arg 8
#python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin" #PARAM 30
python_bin <- "/usr/bin" #PARAM 30
NA_flag_val_mosaic <- -32768
in_dir_list_filename <- NULL #if NULL, use the in_dir directory to search for info
#countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #Atlas
countries_hsp <- "/nobackupp8/bparmen1/NEX_data/countries.shp"
lf_raster <- NULL #list of raster to consider
#On NEX
#contains all data from the run by Alberto
in_dir1 <- "/nobackupp6/aguzman4/climateLayers/out" #On NEX
#parent output dir for the current script analyes
y_var_name <- "dailyTmax" #PARAM1
out_suffix <- "predictions_assessment_reg6_10302016"
file_format <- ".rst" #PARAM 9
NA_value <- -9999 #PARAM10
NA_flag_val <- NA_value
#multiple_region <- TRUE #PARAM 12
region_name <- "reg6" #PARAM 13
countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #PARAM 13, copy this on NEX too
#/nobackupp6/aguzman4/climateLayers/out/reg6/subset/shapefiles
list_year_predicted <- c(2000,2012,2013) #year still on disk for reg6
  
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

#### Function to check missing tiles and estimate potential gaps
predictions_tiles_missing_fun <- function(list_param,i){
  

  ##############################
  #### Parameters and constants  
  

  in_dir1 <- list_param$in_dir1 
  region_name <- list_param$region_name #e.g. c("reg23","reg4") #run only for one region
  y_var_name <- list_param$y_var_name # e.g. dailyTmax" #PARAM3
  interpolation_method <- list_param_run_assessment_prediction$interpolation_method #c("gam_CAI") #PARAM4
  out_suffix <- list_param_run$out_suffix #output suffix e.g."run_global_analyses_pred_12282015" #PARAM5
  out_dir <- list_param$out_dir #<- "/nobackupp8/bparmen1/" #PARAM6
  create_out_dir_param <-list_param$create_out_dir_param #if TRUE output dir created #PARAM7
  proj_str <- list_param$proj_str # CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84, #PARAM8
  list_year_predicted <- list_param$list_year_predicted # 1984:2004
  file_format <- list_param$file_format #<- ".tif" #format for mosaiced files #PARAM10
  NA_flag_val <- list_param$NA_flag_val #<- -9999  #No data value, #PARAM11
  num_cores <- list_param$num_cores #<- 6 #number of cores used #PARAM13
  plotting_figures <- list_param$plotting_figures #if true run part2 of assessment
  
  ##for plotting assessment function
  
  item_no <- list_param_run_assessment_prediction$mosaic_plot  #PARAM14
  day_to_mosaic <- list_param$day_to_mosaic #PARAM15
  countries_shp <- list_param$countries_shp #PARAM17
  plot_region <- list_param$plot_region #PARAM18
  threshold_missing_day <- list_param$threshold_missing_day #PARM20
  pred_mod_name <- list_param$pred_mod_name
  metric_name <- list_param$metric_name
  
  ########################## START SCRIPT #########################################
  
  #system("ls /nobackup/bparmen1")
  out_dir <- in_dir
  if(create_out_dir_param==TRUE){
    out_dir <- create_dir_fun(out_dir,out_suffix)
    setwd(out_dir)
  }else{
    setwd(out_dir) #use previoulsy defined directory
  }
  
  setwd(out_dir)
  
  #list_outfiles <- vector("list", length=35) #collect names of output files, this should be dynamic?
  #list_outfiles_names <- vector("list", length=35) #collect names of output files

  year_predicted <- list_param_run_assessment_prediction$list_year_predicted[i] 

  in_dir1_reg <- file.path(in_dir1,region_name)
  
  list_outfiles <- vector("list", length=14) #collect names of output files
  
  in_dir_list <- list.dirs(path=in_dir1_reg,recursive=FALSE) #get the list regions processed for this run

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

  ##raster_prediction object : contains testing and training stations with RMSE and model object
  in_dir_list_tmp <- file.path(in_dir_list,year_predicted)
  list_raster_obj_files <- mclapply(in_dir_list_tmp,
                                    FUN=function(x){list.files(path=x,pattern="^raster_prediction_obj.*.RData",full.names=T)},
                                    mc.preschedule=FALSE,mc.cores = num_cores)
  
  #list_raster_obj_files <- try(lapply(in_dir_list_tmp,FUN=function(x){list.files(path=x,pattern="^raster_prediction_obj.*.RData",full.names=T)}))
  #Add stop message here...if no raster object in any tiles then break from the function
  
  list_names_tile_coord <- lapply(list_raster_obj_files,FUN=function(x){basename(dirname(x))})
  list_names_tile_id <- paste("tile",1:length(list_raster_obj_files),sep="_")
  names(list_raster_obj_files)<- list_names_tile_id
  
  #pred_mod_name <- "mod1"
  list_lf_raster_tif_tiles <- mclapply(in_dir_list_tmp,
                                    FUN=function(x){list.files(path=x,pattern=paste0("gam_CAI_dailyTmax_predicted_",pred_mod_name,".*.tif"),full.names=T)},
                                    mc.preschedule=FALSE,mc.cores = num_cores)
  list_names_tile_coord <- lapply(list_lf_raster_tif_tiles,FUN=function(x){basename(dirname(dirname(x)))})
  list_names_tile_id <- paste("tile",1:length(list_lf_raster_tif_tiles),sep="_")
  names(list_lf_raster_tif_tiles)<- list_names_tile_id
  
  #one level up
  #lf_covar_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar_obj.*.RData",full.names=T)})
  #lf_covar_tif <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar.*.tif",full.names=T)})
  
  #lf_sub_sampling_obj_files <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern=paste("^sub_sampling_obj_",interpolation_method,".*.RData",sep=""),full.names=T)})
  #lf_sub_sampling_obj_daily_files <- lapply(in_dir_list_tmp,FUN=function(x){list.files(path=x,pattern="^sub_sampling_obj_daily.*.RData",full.names=T)})
  year_processed <- year_predicted
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
  
  #gam_CAI_dailyTmax_predicted_mod1_0_1_20001231_30_1-39.7_165.1.tif
  
  #undebug(check_missing)
  test_missing <- try(lapply(1:length(list_lf_raster_tif_tiles),function(i){check_missing(lf=list_lf_raster_tif_tiles[[i]], 
                                                                      pattern_str=NULL,
                                                                      in_dir=out_dir,
                                                                      date_start=start_date,
                                                                      date_end=end_date,
                                                                      item_no=item_no, #9 for predicted tiles
                                                                      out_suffix=out_suffix,
                                                                      num_cores=num_cores,
                                                                      out_dir=".")}))

 
  df_time_series <- test_missing[[1]]$df_time_series
  head(df_time_series)

  table(df_time_series$missing)
  table(df_time_series$year)
  
  ###Now combined in one table?
  
  list_missing <- lapply(1:length(test_missing),FUN=function(i){df_time_series <- test_missing[[i]]$df_time_series$missing})
  
  df_missing <- as.data.frame(do.call(cbind,list_missing))
  names(df_missing) <- unlist(basename(in_dir_reg))
  df_missing$tot_missing <- rowSums (df_missing, na.rm = FALSE, dims = 1)
  df_missing$reg <- region_name
  df_missing$date <- day_to_mosaic

  filename_df_missing <- paste0("df_missing_by_tiles_predicted_tif_",region_name,"_",pred_mod_name,"_",out_suffix)
  write.table(df_missing,file=filename_df_missing)
  
  ########################
  #### Step 2: examine overlap
  
  path_to_shp <- dirname(countries_shp)
  layer_name <- sub(".shp","",basename(countries_shp))
  reg_layer <- readOGR(path_to_shp, layer_name)
  
  #collect info: read in all shapefiles
  
  centroids_shp_fun <- function(i,list_shp_reg_files){
    #
    shp_filename <- list_shp_reg_files[[i]]
    layer_name <- sub(".shp","",basename(shp_filename))
    path_to_shp <- dirname(shp_filename)
    shp1 <- try(readOGR(path_to_shp, layer_name)) #use try to resolve error below
    #shp_61.0_-160.0
    #Geographical CRS given to non-conformant data: -186.331747678
    
    #shp1<-readOGR(dirname(list_shp_reg_files[[i]]),sub(".shp","",basename(list_shp_reg_files[[i]])))
    if (!inherits(shp1,"try-error")) {
      pt <- gCentroid(shp1)
      #centroids_pts[[i]] <- pt
    }else{
      pt <- shp1
      #centroids_pts[[i]] <- pt
    }
    #shps_tiles[[i]] <- shp1
    #centroids_pts[[i]] <- centroids
    
    shp_obj <- list(shp1,pt)
    names(shp_obj) <- c("spdf","centroid")
    return(shp_obj)
  }

  obj_centroids_shp <- centroids_shp_fun(1,list_shp_reg_files=in_dir_shp_list)
                                         

  obj_centroids_shp <- mclapply(1:length(in_dir_shp_list),
                                FUN=centroids_shp_fun,
                                list_shp_reg_files=in_dir_shp_list,
                                mc.preschedule=FALSE,
                                mc.cores = num_cores)

  centroids_pts <- lapply(obj_centroids_shp, FUN=function(x){x$centroid})
  shps_tiles <-   lapply(obj_centroids_shp, FUN=function(x){x$spdf})

  #remove try-error polygons...we loose three tiles because they extend beyond -180 deg
  tmp <- shps_tiles
  shps_tiles <- remove_errors_list(shps_tiles) #[[!inherits(shps_tiles,"try-error")]]
  #shps_tiles <- tmp
  length(tmp)-length(shps_tiles) #number of tiles with error message
  
  tmp_pts <- centroids_pts 
  centroids_pts <- remove_errors_list(centroids_pts) #[[!inherits(shps_tiles,"try-error")]]
  #centroids_pts <- tmp_pts 
  
  r <- raster(infile_mask)
  plot(r)
  plot(shp1,add=T,border="blue",usePolypath = FALSE) #added usePolypath following error on brige and NEX

  ## find overlap
  #http://gis.stackexchange.com/questions/156660/loop-to-check-multiple-polygons-overlap-in-r
  
  matrix_overlap <- matrix(data=NA,nrow=length(shps_tiles),ncol=length(shps_tiles))
  for(i in 1:length(shps_tiles)){
     for(j in 2:length(shps_tiles)){
      overlap_val <- as.numeric(over(shps_tiles[[i]],shps_tiles[[j]]))
      matrix_overlap[i,j]<- overlap_val
    }
    #
  }
  
  matrix_overlap%*%df_missing[1,1:26]
  
  ## For each day can do overalp matrix* prediction
  ## if prediction and overlap then 1 else 0, if no-overlap is then NA
  ## then for each tile compute the number of excepted predictions taken into account in a tile
  
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
