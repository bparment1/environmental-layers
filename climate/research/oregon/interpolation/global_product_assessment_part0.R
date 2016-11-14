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
#MODIFIED ON: 11/14/2016            
#Version: 1
#PROJECT: Environmental Layers project     
#COMMENTS: 
#TODO:
#1) 
#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rm u:aguzman4:rwx /nobackupp6/aguzman4/climateLayers/LST_tempSpline/
#COMMIT: adding and debugging function to generate raster of number of predictions for day with missing tiles   

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
#library(mosaic) #not installed on NEX

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
function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_11052016.R"
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
#python_bin <- "/usr/bin" #PARAM 30, NCEAS
python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin" #PARAM 30"
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
#file_format <- ".rst" #PARAM 9
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
  
  #####
  #Now combine df_time_series in one table
  
  dim(test_missing[[1]]$df_time_series)
  list_lf <- lapply(1:length(test_missing),FUN=function(i){df_time_series <- as.character(test_missing[[i]]$df_time_series$lf)})
  df_lf_tiles_time_series <- as.data.frame(do.call(cbind,list_lf))
  #http://stackoverflow.com/questions/26220913/replace-na-with-na
  #Use dfr[dfr=="<NA>"]=NA where dfr is your dataframe.
  names(df_lf_tiles_time_series) <- unlist(basename(in_dir_reg))
  filename_df_lf_tiles <- paste0("df_files_by_tiles_predicted_tif_",region_name,"_",pred_mod_name,"_",out_suffix)
  write.table(df_lf_tiles_time_series,file=filename_df_lf_tiles)

  ###Now combined missing in one table?
  
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
  
  r_mask <- raster(infile_mask)
  plot(r)
  plot(shps_tiles[[1]],add=T,border="blue",usePolypath = FALSE) #added usePolypath following error on brige and NEX

  ## find overlap
  #http://gis.stackexchange.com/questions/156660/loop-to-check-multiple-polygons-overlap-in-r
  
  matrix_overlap <- matrix(data=NA,nrow=length(shps_tiles),ncol=length(shps_tiles))
  for(i in 1:length(shps_tiles)){
     for(j in 1:length(shps_tiles)){
      overlap_val <- as.numeric(over(shps_tiles[[i]],shps_tiles[[j]]))
      matrix_overlap[i,j]<- overlap_val
    }
    #
  }
  
  names(shps_tiles) <- basename(unlist(in_dir_reg))
  r_ref <- raster(list_lf_raster_tif_tiles[[1]][1])
  list_r_ref <- lapply(1:length(in_dir_reg), function(i){raster(list_lf_raster_tif_tiles[[i]][1])})
  tile_spdf <- shps_tiles[[1]]
  tile_coord <- basename(in_dir_reg[1])
  date_val <- df_missing$date[1]
  
  ### use rasterize
  spdf_tiles <- do.call(bind, shps_tiles) #bind all tiles together in one shapefile

  #function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_11052016.R"
  #source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 

  #undebug(rasterize_tile_day)
  list_predicted <- rasterize_tile_day(1,
           list_spdf=shps_tiles,
           df_missing=df_missing,
           list_r_ref=list_r_ref,
           col_name="overlap",
           date_val=df_missing$date[1])
  #list_predicted <- mclapply(1:6,
  #         FUN=rasterize_tile_day,
  #         list_spdf=shps_tiles,
  #         df_missing=df_missing,
  #         list_r_ref=list_r_ref,
  #         col_name = "overlap",
  #         date_val=df_missing$date[1],
  #          mc.preschedule=FALSE,
  #         mc.cores = num_cores)
  
  list_predicted <- mclapply(1:length(shps_tiles),
           FUN=rasterize_tile_day,
           list_spdf=shps_tiles,
           df_missing=df_missing,
           list_r_ref=list_r_ref,
           col_name = "overlap",
           date_val=df_missing$date[1],
            mc.preschedule=FALSE,
           mc.cores = num_cores)

  ##check that everything is correct:
  plot(r_mask)
  plot(raster(list_predicted[[1]]),add=T)
  plot(spdf_tiles_test,add=T,border="green",usePolypath = FALSE) #added usePolypath following error on brige and NEX

  ### Make a list of file
  out_suffix_str_tmp <- paste0(region_name,"_",out_suffix)
  out_dir_str <- out_dir
  filename_list_predicted <- file.path(out_dir_str,paste("list_to_mosaics_",out_suffix_str_tmp,".txt",sep=""))

  #writeLines(unlist(list_weights_m),con=filename_list_mosaics_weights_m) #weights files to mosaic 
  #writeLines(unlist(list_weights_prod_m),con=filename_list_mosaics_prod_weights_m) #prod weights files to mosaic
      
  writeLines(unlist(list_predicted),con=filename_list_predicted) #weights files to mosaic 

  #out_mosaic_name_weights_m <- r_weights_sum_raster_name <- file.path(out_dir,paste("r_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix,".tif",sep=""))
  #out_mosaic_name_prod_weights_m <- r_weights_sum_raster_name <- file.path(out_dir,paste("r_prod_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix,".tif",sep=""))
  out_mosaic_name_predicted_m  <- file.path(out_dir_str,paste("r_overlap_sum_m_",out_suffix_str_tmp,".tif",sep=""))

  rast_ref_name <- infile_mask
  mosaic_python <- "/nobackupp6/aguzman4/climateLayers/sharedCode/"
  rast_ref_name <- infile_mask
  #python /nobackupp6/aguzman4/climateLayers/sharedCode//gdal_merge_sum.py --config GDAL_CACHEMAX=1500 --overwrite=TRUE -o /nobackupp8/bparmen1/climateLayers/out
  mosaic_overlap_tiles_obj <- mosaic_python_merge(NA_flag_val=NA_flag_val,
                                                module_path=mosaic_python,
                                                module_name="gdal_merge_sum.py",
                                                input_file=filename_list_predicted,
                                                out_mosaic_name=out_mosaic_name_predicted_m,
                                                raster_ref_name = rast_ref_name) ##if NA, not take into account
  r_overlap_raster_name <- mosaic_overlap_tiles_obj$out_mosaic_name
  cmd_str1 <-   mosaic_overlap_tiles_obj$cmd_str

  r_overlap <- raster(r_overlap_raster_name)
  r_mask <- raster(infile_mask)
  
  out_mosaic_name_overlap_masked  <- file.path(out_dir_str,paste("r_overlap_sum_masked_",out_suffix_str_tmp,".tif",sep=""))

  r_overlap_m <- mask(r_overlap,r_mask,filename=out_mosaic_name_overlap_masked)
  #plot(r_overlap_m)
  #plot(spdf_tiles_test,add=T,border="green",usePolypath = FALSE) #added usePolypath following error on brige and NEX
  
  r_table <- ratify(r_overlap_m) # build the Raster Attibute table
  rat <- levels(r_table)[[1]]#get the values of the unique cell frot the attribute table
  #rat$legend <- paste0("tile_",1:26)
  tb_freq <- as.data.frame(freq(r_table))
  rat$legend <- tb_freq$value
  levels(r_table) <- rat
  
  
  res_pix <- 800
  #res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1
  
  png_filename <-  file.path(out_dir,paste("Figure_maximum_overlap_",region_name,"_",out_suffix,".png",sep =""))
    
  title_str <-  paste("Maximum overlap: Number of predicted pixels for ",variable_name, sep = "")
  
  png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
    #my_col=c('blue','red','green')
  my_col <- rainbow(length(tb_freq$value))

  plot(r_table,col=my_col,legend=F,box=F,axes=F,main=title_str)
  legend(x='topright', legend =rat$legend,fill = my_col,cex=0.8)
  
  dev.off()

  #r <- raster(matrix(runif(20),5,4))
  #r[r>.5] <- NA
  #g <- as(r, 'SpatialGridDataFrame')
  #p <- as(r, 'SpatialPixels')
  #p_spdf <- as(r_overlap_m,"SpatialPointsDataFrame")
  #plot(r)
  #points(p)
  
  ### now assign id and match extent for tiles
  
  lf_files <- unlist(list_predicted)
  rast_ref_name <- infile_mask
  rast_ref <- rast_ref_name
  
  ##Maching resolution is probably only necessary for the r mosaic function
  #Modify later to take into account option R or python...
  list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix_str_tmp,out_dir_str)
  names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")

  #undebug(raster_match)
  #r_test <- raster_match(1,list_param_raster_match)
  #r_test <- raster(raster_match(1,list_param_raster_match))

  list_tiles_predicted_m <- unlist(mclapply(1:length(lf_files),
                                            FUN=raster_match,list_param=list_param_raster_match,
                                            mc.preschedule=FALSE,mc.cores = num_cores))                           

  extension_str <- extension(lf_files)
  raster_name_tmp <- gsub(extension_str,"",basename(lf_files))
  out_suffix_str <- paste0(region_name,"_",out_suffix)
  raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_","masked_",out_suffix_str,file_format,sep=""))
  
  #writeRaster(r, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  

  #r_stack <- stack(list_tiles_predicted_m)
  list_mask_out_file_name <- raster_name
  list_tiles_predicted_masked <- unlist(mclapply(1:length(list_tiles_predicted_m),
                                                 FUN=function(i){mask(raster(list_tiles_predicted_m[i]),r_mask,filename=list_mask_out_file_name[i])},
                                                       mc.preschedule=FALSE,mc.cores = num_cores))                         
  #r_stack_masked <- mask(r, m2) #, maskvalue=TRUE)
  
  ########################
  #### Step 3: combine overlap information and number of predictions by day
  ##Now loop through every day if missing then generate are raster showing map of number of prediction
  
  #r_tiles_stack <- stack(list_tiles_predicted_masked)
  #names(r_tiles_stack) <- basename(in_dir_reg) #this does not work, X. is added to the name, use list instead
  
  names(list_tiles_predicted_masked) <- basename(in_dir_reg)
  df_missing_tiles_day <- subset(df_missing,tot_missing > 0)
  #r_tiles_s <- r_tiles_stack
  names_tiles <- basename(in_dir_reg)
  
  generate_raster_number_of_prediction_by_day <- function(i,list_param){
    
    list_names_tile_coord
    df_time_series
    missing_tiles <- df_missing_tiles_day[i,]
    date_str <- missing_tiles$date
    
    #r_tiles_s <- list_param$r_tiles_s
    list_param$list_tiles_predicted_masked
    
    selected_col <- names(list_tiles_predicted_masked)
    missing_tiles_subset <- subset(missing_tiles,select=selected_col)
    selected_missing <- missing_tiles_subset==1
    #names(df_missing_tiles_day)[selected_missing]
    
    #names(list_tiles_predicted_masked)[selected_missing]
    list_missing_tiles_raster <- list_tiles_predicted_masked[selected_missing]
    r_tiles_s <- stack(list_missing_tiles_raster)
    
    ### first sum missing
    datasum <- stackApply(r_tiles_s, 1:nlayers(r_tiles_s), fun = sum)
     
    ### then substract missing tiles...
    r_day_predicted <- r_overlap_m - datasum
    
    r_table <- ratify(r_day_predicted) # build the Raster Attibute table
    rat <- levels(r_table)[[1]]#get the values of the unique cell frot the attribute table
    #rat$legend <- paste0("tile_",1:26)
    tb_freq <- as.data.frame(freq(r_table))
    rat$legend <- tb_freq$value
    levels(r_table) <- rat

    res_pix <- 800
    #res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
  
    png_filename <-  file.path(out_dir,paste("Figure_number_of_predictionds_by_pixel_",date_str,"_",region_name,"_",out_suffix,".png",sep =""))
    
    title_str <-  paste("Number of predicted pixels for ",variable_name," on ",date_str, sep = "")
  
    png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
    #my_col=c('blue','red','green')
    my_col <- rainbow(length(tb_freq$value))
    plot(r_table,col=my_col,legend=F,box=F,axes=F,main=title_str)
    legend(x='topright', legend =rat$legend,fill = my_col,cex=0.8)
  
    dev.off()

    ### Day missing reclass above
    
    ## do this in gdalcalc?
    r_missing_day <- r_day_predicted == 0
    
    res_pix <- 800
    #res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
  
    png_filename <-  file.path(out_dir,paste("Figure_missing_predictionds_by_pixel_",date_str,"_",region_name,"_",out_suffix,".png",sep =""))
    
    title_str <-  paste("Number of predicted pixels for ",variable_name," on ",date_str, sep = "")
  
    png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
    #my_col=c('blue','red','green')
    my_col <- c("black","red")
    plot(r_missing_day,col=my_col,legend=F,box=F,axes=F,main=title_str)
    legend(x='topright', legend =c("prediced","missing"),fill = my_col,cex=0.8)
  
    dev.off()

    ### generate retunr object
    returm(r_day_predicted)
  }
  
  ## Make a function,
  #for specifi i (day) select tile with missing info, load it and subsetract to overlap raster, save it.
  
  #http://stackoverflow.com/questions/19586945/how-to-legend-a-raster-using-directly-the-raster-attribute-table-and-displaying
  #
  #http://gis.stackexchange.com/questions/148398/how-does-spatial-polygon-over-polygon-work-to-when-aggregating-values-in-r
  #ok other way of doing this:
  #1. find overlap assuming all predictions!
  #2. Us raster image with number of overlaps in the mosaic tile
  #3. for every pixel generate and ID (tile ID) as integer, there should  be 26 layers at the mosaic extent
  #4. generate a table? for each pixel it can say if is part of a specific tile
  #5. workout a formula to generate the number of predictions for each pixel based on tile predicted for each date!!

  

  
  
  

  return()
}

############################ END OF SCRIPT ##################################

  # spdf_tiles <- do.call(intersect, shps_tiles)
  # 
  # ### Now use intersect to retain actual overlap
  # 
  # for(i in 1:length(shps_tiles)){
  #   overlap_intersect <- intersect(shps_tiles[[1]],shps_tiles[[i]])
  # }
  # 
  # overlap_intersect <- lapply(1:length(shps_tiles),FUN=function(i){intersect(shps_tiles[[1]],shps_tiles[[i]])})
  # #test <- overlap_intersect <- intersect(shps_tiles[[1]],shps_tiles[[2]]))
  # 
  # names(overlap_intersect) <- basename(in_dir_reg)
  # shp_selected <- unlist(lapply(1:length(overlap_intersect),function(i){!is.null(overlap_intersect[[i]])}))
  # test_list <- overlap_intersect[shp_selected]
  # spdf_tiles_test <- do.call(bind, test_list) #combines all intersect!!
  # #ll <- ll[ ! sapply(ll, is.null) ]
  # test <- overlap_intersect[!lapply(overlap_intersect,is.null)]
  # spdf_tiles_test <- do.call(bind, test) #combines all intersect!!
  # #ll <- ll[ ! sapply(ll, is.null) ]
  # spdf_tiles <- do.call(bind, overlap_intersect[1:4]) #combines all intersect!!
  # spdf_tiles_test <- do.call(bind, test) #combines all intersect!!
  # 
  # plot(spdf_tiles_test,add=T,border="green",usePolypath = FALSE) #added usePolypath following error on brige and NEX
  # 
  # matrix_overlap%*%df_missing[1,1:26]
  # 
  # 
  # ## For each day can do overalp matrix* prediction
  # ## if prediction and overlap then 1 else 0, if no-overlap is then NA
  # ## then for each tile compute the number of excepted predictions taken into account in a tile
  # 
  # #combine polygon
  # #http://gis.stackexchange.com/questions/155328/merging-multiple-spatialpolygondataframes-into-1-spdf-in-r
  # 
  # #http://gis.stackexchange.com/questions/116388/count-overlapping-polygons-in-single-shape-file-with-r
  # 
  # #### Use the predictions directory
  # #By region
  # #For each polygon/tile find polygon overlapping with count and ID (like list w)
  # #for each polygon/tile and date find if there is a prediction using the tif (multiply number)
  # #for each date of year report data in table.
  # 
  # #go through table and hsow if there are missing data (no prediction) or report min predictions for tile set?
  #   
  # #for each polygon find you overlap!!
  # #plot number of overlap
  # #for specific each find prediction...
