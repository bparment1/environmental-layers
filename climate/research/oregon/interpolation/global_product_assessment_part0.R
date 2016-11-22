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
#MODIFIED ON: 11/22/2016            
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
#COMMIT: moving function of number of predictions for day with missing tiles to functon script

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
#library(mosaic) #not installed on NEX, drop for the workflow

###### Function used in the script #######
  
#Adding command line arguments to use mpiexec
args<-commandArgs(TRUE)

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
function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_11222016b.R"
source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 
##Don't load part 1 and part2, mosaic package does not work on NEX
#function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_09192016b.R"
#source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 
#function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_10222016.R"
#source(file.path(script_path,function_product_assessment_part2_functions)) #source all functions used in this script 

###############################
####### Parameters, constants and arguments ###

#Find number of param
#var <- args[1] # variable being interpolated #param 1, arg 1
#var<-"TMAX" # variable being interpolated #param 1, arg 1

#CRS_locs_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #constant 1
proj_str <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"

var<-"TMAX" # variable being interpolated #param 1, arg 1
interpolation_method<-c("gam_CAI") #param 2
#CRS_interp <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #param 3
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
out_dir <- "/nobackupp6/aguzman4/climateLayers/out/reg6/assessment"
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
countries_shp <- "/nobackupp8/bparmen1/NEX_data/countries.shp"
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


##### prepare list of parameters for call of function
function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_11222016b.R"
source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 

list_param_predictions_tiles_missing <- list(in_dir1,region_name,y_var_name,interpolation_method,out_suffix,out_dir,
                                             create_out_dir_param,proj_str,list_year_predicted,file_format,NA_flag_val,
                                             num_cores,plotting_figures,item_no,day_to_mosaic_range,countries_shp,plotting_figures,
                                             #threshold_missing_day,
                                             pred_mod_name,metric_name)

names(list_param_predictions_tiles_missing) <- c("in_dir1","region_name","y_var_name","interpolation_method","out_suffix","out_dir",
                                             "create_out_dir_param","proj_str","list_year_predicted","file_format","NA_flag_val",
                                             "num_cores","plotting_figures","item_no","day_to_mosaic_range","countries_shp","plotting_figures",
                                             #"threshold_missing_day",
                                             "pred_mod_name","metric_name")

#debug(predictions_tiles_missing_fun)
#Started at 11.06pm
predictions_tiles_missing_fun(1,list_param=list_param_predictions_tiles_missing)

#### Function to check missing tiles and estimate potential gaps

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
