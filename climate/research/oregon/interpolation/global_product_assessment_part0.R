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
#MODIFIED ON: 11/28/2016            
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
#COMMIT: making callable from shel for function of number of predictions for day with missing tiles 

### Testing several years on the bridge before running jobs on nodes with qsub
#Use the following command to run as script via the shell on the bridge 
#run mosaic of predictions

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
function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_11292016b.R"
source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 
##Don't load part 1 and part2, mosaic package does not work on NEX
#function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_09192016b.R"
#source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 
#function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_10222016.R"
#source(file.path(script_path,function_product_assessment_part2_functions)) #source all functions used in this script 

###############################
####### Parameters, constants and arguments ###

#Rscript /nobackupp8/bparmen1/env_layers_scripts/global_product_assessment_part0_11272016.R TMAX /nobackupp6/aguzman4/climateLayers/out/reg6/assessment reg6 predictions_assessment_reg6_10302016 /nobackupp8/bparmen1/climateLayers/out/reg6/assessment TRUE 2000 6 1e+07 9 rmse 20000101 20001231 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg6.tif /nobackupp6/aguzman4/climateLayers/out var_pred

### ARGUMENTS: inputs parameters set from the command line

var <- args[1] # variable being interpolated #param 1, arg 1
in_dir <- args[2] #PARAM2,#region_name <- "reg4" #PARAM 3 #reg4 South America, Africa reg5,Europe reg2, North America reg1, Asia reg3
region_name <- args[3] #PARAM3
out_suffix <- args[4] #PARAM 4
out_suffix_str <- region_name #PARAM 4, CONST 3
out_dir <- args[5] #PARAM 5
create_out_dir_param <- args[6] #PARAM 6
year_predicted <- args[7] #PARAM 7
num_cores <- args[8] #PARAM 8
max_mem<-args[9] #PARAM 9
##mosaicing_method <- c("unweighted","use_edge_weights") #PARAM10
item_no <- args[10] #PARAM10
metric_name <- args[11]
day_start <- args[12] #PARAM 12
day_end <- args[13] #PARAM 13
infile_mask <- args[14]
in_dir1 <- args[15] #PARAM 15, files containing assessment information
layers_option <- args[16] # PARAM 17 options are:
tmp_files <- args[17]

#### values used for testing
var <- "TMAX" # variable being interpolated #PARAM 1, arg 1
in_dir <- "/nobackupp6/aguzman4/climateLayers/out/reg6/assessment" #PARAM2
region_name <- c("reg6") #PARAM 3, arg 3
out_suffix <- "predictions_assessment_reg6_10302016" #PARAM 4
#out_suffix_str <- region_name #PARAM 4, CONST 3
out_dir <- "/nobackupp6/aguzman4/climateLayers/out/reg6/assessment" #PARAM 5
create_out_dir_param <- TRUE #PARAM 12, arg 6
year_predicted <- c(2000) #PARAM 7, arg7
num_cores <- 6 #number of cores used # PARAM 8, arg 8
max_mem <- 1e+07 #PARAM 9
#mosaicing_method <- args[10] #PARAM10
item_no <- 9 #PARAM 10, arg 10
metric_name <- "rmse" # "mae", "r" for MAE, R etc.; can also be ns or nv? #PARAM 11, arg 11
day_start <- "20000101" #PARAM 12, arg 12
day_end <- "20001231" #PARAM 13, arg 13
infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg6.tif" #PARAM 14, arg 14
in_dir1 <- "/nobackupp6/aguzman4/climateLayers/out" # PARAM 15 On NEX
layers_option <- c("var_pred") #PARAM 16, arg 16
tmp_files <- FALSE #PARAM 17, arg 17

###################
### CONSTANT: not set from command line

pred_mod_name <- "mod1"
date_start <- day_start
date_end <- day_end
NA_value <- -32768 #PARAM 26
NA_flag_val <- NA_value #PARAM 26
proj_str <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"
interpolation_method<-c("gam_CAI") #PARAM 2
day_to_mosaic_range <- NULL #PARAM 7
file_format <- ".tif" #format for mosaiced files # PARAM 14
plotting_figures <- TRUE #running part2 of assessment to generate figures... # PARAM 13
#python_bin <- "/usr/bin" #PARAM 15, NCEAS
python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin" #PARAM 15"
in_dir_list_filename <- NULL # PARAM 16, if NULL, use the in_dir directory to search for info
#countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #Atlas
countries_shp <- "/nobackupp8/bparmen1/NEX_data/countries.shp" #PARAM 17
lf_raster <- NULL #list of raster to consider #PARAM 18
scaling <- 1
data_type <- "Int16"
#CRS_interp <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #param 3
#list_models<-c("y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)") #param 4

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

if(!(is.null(day_start)) & !(is.null(day_end))){
  day_to_mosaic_range <- c(day_start,day_end) #if null run all year #PARAM 12
}else{
  day_to_mosaic_range <- NULL
}

#parse input value range
#values_range <- as.numeric(unlist(strsplit(values_range,",")))
#scaling <- 1
scaling <- as.numeric(scaling)
#data_type <- "Int16"
#i<-1

##### prepare list of parameters for call of function

list_param_predictions_tiles_missing <- list(in_dir1,region_name,y_var_name,interpolation_method,out_suffix,out_dir,
                                             create_out_dir_param,proj_str,year_predicted,file_format,NA_flag_val,
                                             num_cores,plotting_figures,item_no,day_to_mosaic_range,countries_shp,plotting_figures,
                                             scaling, data_type, python_bin,tmp_files,
                                             pred_mod_name,metric_name)

names(list_param_predictions_tiles_missing) <- c("in_dir1","region_name","y_var_name","interpolation_method","out_suffix","out_dir",
                                             "create_out_dir_param","proj_str","year_predicted","file_format","NA_flag_val",
                                             "num_cores","plotting_figures","item_no","day_to_mosaic_range","countries_shp","plotting_figures",
                                             "scaling", "data_type", "python_bin","tmp_files",
                                             "pred_mod_name","metric_name")

#Product assessment
function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_11292016b.R"
source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 

#debug(predictions_tiles_missing_fun)
#Started at 11.06pm
obj_predictions_tiles_missing_fun <- predictions_tiles_missing_fun(list_param=list_param_predictions_tiles_missing)


############################ END OF SCRIPT ##################################


