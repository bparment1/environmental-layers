##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part 0 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This script checks the number of predictions by tiles and years.
#with the goal of predicting potential gaps or missing predictions in fugure mosaics by region.
#The general logic is to check the number of overlap by shapefile polyon tiles
#along with the predictions for every day of the year (*.tif)
#Summary tables and data are also produced in the script.
#
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/27/2016  
#MODIFIED ON: 05/28/2017            
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
#COMMIT: testing command line for reg 6 after debugging

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
function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_04142017.R" #Functions used to mosaic predicted tiles
function_mosaicing <-"global_run_scalingup_mosaicing_04172017.R" #main scripts for mosaicing predicted tiles

source(file.path(script_path,function_mosaicing)) #source all functions used in this script 
source(file.path(script_path,function_mosaicing_functions)) #source all functions used in this script 

#Assessment on NEX
function_assessment_part1_functions <- "global_run_scalingup_assessment_part1_functions_12282015.R" #PARAM12
function_assessment_part1a <-"global_run_scalingup_assessment_part1a_02242017.R"
function_assessment_part2 <- "global_run_scalingup_assessment_part2_02092016.R"
function_assessment_part2_functions <- "global_run_scalingup_assessment_part2_functions_01032016.R"
function_assessment_part3 <- "global_run_scalingup_assessment_part3_07292016.R"

source(file.path(script_path,function_assessment_part1_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part1a)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part3)) #source all functions used in this script 

#Product assessment
function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_05282017.R"
source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 
##Don't load part 1 and part2, mosaic package does not work on NEX
#function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_09192016b.R"
#source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 
#function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_10222016.R"
#source(file.path(script_path,function_product_assessment_part2_functions)) #source all functions used in this script 

###############################
####### Parameters, constants and arguments ###

#Rscript /nobackupp8/bparmen1/env_layers_scripts/global_product_assessment_part0_12182016b.R TMAX /nobackupp6/aguzman4/climateLayers/out/reg6/assessment reg6 predictions_assessment_reg6_2000_test2 /nobackupp8/bparmen1/climateLayers/out/reg6/assessment TRUE 2000 6 1e+07 9 rmse 20000101 20001231 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg6.tif /nobackupp6/aguzman4/climateLayers/out var_pred FALSE FALSE
#Rscript /nobackupp8/bparmen1/env_layers_scripts/global_product_assessment_part0_12182016b.R TMAX /nobackupp6/aguzman4/climateLayers/out/reg6/assessment reg6 predictions_tiles_assessment_reg6_2000_test3 /nobackupp8/bparmen1/climateLayers/out/reg6/assessment TRUE 2000 6 1e+07 9 rmse 20000101 20001231 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg6.tif /nobackupp6/aguzman4/climateLayers/out var_pred FALSE FALSE
#Rscript /nobackupp8/bparmen1/env_layers_scripts/global_product_assessment_part0_04292017.R TMIN /nobackupp6/aguzman4/climateLayers/tMinOut/reg1/assessment reg1 predictions_gaps_tiles_assessment_reg1_1985 /nobackupp8/bparmen1/climateLayers/tMinOut/reg1/assessment TRUE 1985 6 1e+07 9 rmse 19850101 19851231 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1.tif /nobackupp6/aguzman4/climateLayers/tMinOut var_pred FALSE FALSE FALSE FALSE
 
#Rscript /nobackupp8/bparmen1/env_layers_scripts/global_product_assessment_part0_05282017.R TMIN /nobackupp6/aguzman4/climateLayers/tMinOut/reg6/assessment reg6 predictions_gaps_tiles_assessment_reg6_1985 /nobackupp8/bparmen1/climateLayers/tMinOut/reg6/assessment TRUE 1985 6 1e+07 9 rmse 19850101 19851231 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg6.tif /nobackupp6/aguzman4/climateLayers/tMinOut var_pred FALSE FALSE FALSE FALSE                                                                                                
#Rscript /nobackupp8/bparmen1/env_layers_scripts/global_product_assessment_part0_05282017.R TMIN /nobackupp6/aguzman4/climateLayers/tMinOut/reg5/assessment reg5 predictions_gaps_tiles_assessment_reg5_1985 /nobackupp8/bparmen1/climateLayers/tMinOut/reg5/assessment TRUE 1985 6 1e+07 9 rmse 19850101 19851231 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg5.tif /nobackupp6/aguzman4/climateLayers/tMinOut var_pred FALSE FALSE FALSE FALSE                                                                                                

### ARGUMENTS: inputs parameters set from the command line

var <- args[1] # variable being interpolated #param 1, arg 1
in_dir <- args[2] #input dir containing tiles predictions from stage 4 workflow
region_name <- args[3] #PARAM3  #reg4 South America, Africa reg5,Europe reg2, North America reg1, Asia reg3
out_suffix <- args[4] #PARAM 4 # output suffix, add region and year of assessment
out_dir <- args[5] #PARAM 5, parent output dir, a new dir is generated using the "output_"+out_suffix
create_out_dir_param <- args[6] #PARAM 6, if true create out_dir otherwise use given out_dir
year_predicted <- args[7] #PARAM 7, year being assessed
num_cores <- args[8] #PARAM 8, number of cores used in the parraleliation
max_mem<-args[9] #PARAM 9, maximum memory used in raster package
item_no <- args[10] #PARAM10, string position of date in tile tif prediction, use 9 as default
metric_name <- args[11] #PARAM 11, prediction or accuracy: rmse, mae
day_start <- args[12] #PARAM 12, start of day to process
day_end <- args[13] #PARAM 13, end of day to process
infile_mask <- args[14]#PARAM 14, input mask file for the region
in_dir1 <- args[15] #PARAM 15, files containing assessment information
layers_option <- args[16] # PARAM 16 options are: prediction or accuracy
tmp_files <- args[17] # PARAM 17, if FALSE, temporary files are removed
plotting_figures <- args[18]# PARAM 18, if TRUE, png files are produced for missing tiles and day predicted
raster_overlap <- args[19] # PARAM 19, if TRUE, raster overlap is generated
raster_pred <- args[20] # PARAM 20, if TRUE, raster prediction is generated

#### values used for testing
# var <- "TMIN" # variable being interpolated #PARAM 1, arg 1
# in_dir <- "/nobackupp6/aguzman4/climateLayers/tMinOut/reg6/assessment" #PARAM2
# region_name <- c("reg6") #PARAM 3, arg 3
# out_suffix <- "predictions_gaps_tiles_assessment_reg6_1985" #PARAM 4
# #out_suffix_str <- region_name #PARAM 4, CONST 3
# out_dir <- "/nobackupp8/bparmen1/climateLayers/tMinOut/reg6/assessment"
# create_out_dir_param <- TRUE #PARAM 12, arg 6
# year_predicted <- c(1985) #PARAM 7, arg7
# num_cores <- 6 #number of cores used # PARAM 8, arg 8
# max_mem <- 1e+07 #PARAM 9
# item_no <- 9 #PARAM 10, arg 10
# metric_name <- "rmse" # "mae", "r" for MAE, R etc.; can also be ns or nv? #PARAM 11, arg 11
# day_start <- "19850101" #PARAM 12, arg 12
# day_end <- "19851231" #PARAM 13, arg 13
# infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg6.tif" #PARAM 14, arg 14
# in_dir1 <- "/nobackupp6/aguzman4/climateLayers/tMinOut" # PARAM 15 On NEX
# layers_option <- c("var_pred") #PARAM 16, arg 16
# tmp_files <- FALSE #PARAM 17, arg 17
# plotting_figures <- TRUE #PARAm 18, arg 18
# raster_overlap <- FALSE # PARAM 19, if TRUE, raster overlap is generated
# raster_pred <- FALSE # PARAM 20, if TRUE, raster prediction is generated

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
#python_bin <- "/usr/bin" #PARAM 15, NCEAS
python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin" #PARAM 15"
in_dir_list_filename <- NULL # PARAM 16, if NULL, use the in_dir directory to search for info
#countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #Atlas
countries_shp <- "/nobackupp8/bparmen1/NEX_data/countries.shp" #PARAM 17
lf_raster <- NULL #list of raster to consider #PARAM 18
scaling <- 1
data_type <- "Int16"
out_suffix_str <- region_name #PARAM 4, CONST 3

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
item_no <- as.numeric(item_no) #PARAM 10, arg 10, args does not take numeric inputs!!

#data_type <- "Int16"

##### prepare list of parameters for call of function

list_param_predictions_tiles_missing <- list(in_dir1,region_name,y_var_name,interpolation_method,out_suffix,out_dir,
                                             create_out_dir_param,proj_str,year_predicted,file_format,NA_flag_val,
                                             num_cores,plotting_figures,item_no,day_to_mosaic_range,countries_shp,plotting_figures,
                                             scaling, data_type, python_bin,tmp_files,
                                             pred_mod_name,metric_name,raster_overlap,raster_pred)

names(list_param_predictions_tiles_missing) <- c("in_dir1","region_name","y_var_name","interpolation_method","out_suffix","out_dir",
                                             "create_out_dir_param","proj_str","year_predicted","file_format","NA_flag_val",
                                             "num_cores","plotting_figures","item_no","day_to_mosaic_range","countries_shp","plotting_figures",
                                             "scaling", "data_type", "python_bin","tmp_files",
                                             "pred_mod_name","metric_name","raster_overlap","raster_pred")

#Product assessment
#function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_05282017.R"
#source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 

#debug(predictions_tiles_missing_fun)
#Started at 9.35am
obj_predictions_tiles_missing_fun <- predictions_tiles_missing_fun(list_param=list_param_predictions_tiles_missing)
 
###Generate summary from object here to simplify output?

############################ END OF SCRIPT ##################################


