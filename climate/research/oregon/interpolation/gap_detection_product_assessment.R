####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part 0 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This script checks the number of predictions by tiles and years.
#with the goal of predicting potential gaps or missing predictions in fugure mosaics by region.
#The general logic is to check the number of overlap by shapefile polyon tiles
#along with the predictions for every day of the year (*.tif)
#Summary tables and data are also produced in the script.
#
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/31/2016  
#MODIFIED ON: 07/04/2017            
#Version: 1
#PROJECT: Environmental Layers project     
#COMMENTS: removing unused functions and clean up for part0 global product assessment part0 
#TODO:#PROJECT: Environmental Layers project     
#COMMENTS:
#TODO:
#1) Add plot broken down by year and region 
#2) Modify code for overall assessment accross all regions and year
#3) Clean up

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rmd user:aguzman4:rwx /nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015

##COMMIT: modified gap detection function for overlap plot and return object

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
#library(mosaic)

###### Function used in the script #######
#Product assessment
function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_06072017.R"
function_product_assessment_gap_detection_functions <- "gap_detection_product_assessment_functions_07042017.R"

script_path <- "/nobackupp8/bparmen1/env_layers_scripts"
source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 
source(file.path(script_path,function_product_assessment_gap_detection_functions)) #source all functions used in this script 

###############################
####### Parameters, constants and arguments ###


#################### Begin Script ##################################

var <- "TMIN" # variable being interpolated #PARAM 1, arg 1
in_dir <- "/nobackupp6/aguzman4/climateLayers/tMinOut/testGaps" #PARAM2
region_name <- c("reg1") #PARAM 3, arg 3
out_suffix <- "mosaic_gaps_tiles_assessment_reg1_combined" #PARAM 4
#out_suffix_str <- region_name #PARAM 4, CONST 3
out_dir <- "/nobackupp8/bparmen1/climateLayers/tMinOut/testGaps" #PARAM 5
create_out_dir_param <- TRUE #PARAM 12, arg 6
num_cores <- 6 #number of cores used # PARAM 8, arg 8
max_mem <- 1e+07 #PARAM 9
metric_name <- "rmse" # "mae", "r" for MAE, R etc.; can also be ns or nv? #PARAM 11, arg 11
#day_start <- "19930101" #PARAM 12, arg 12
#day_end <- "19931231" #PARAM 13, arg 13
infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1.tif" #PARAM 14, arg 14
in_dir1 <- "/nobackupp6/aguzman4/climateLayers/tMinOut" # PARAM 15 On NEX
layers_option <- c("var_pred") #PARAM 16, arg 16
tmp_files <- FALSE #PARAM 17, arg 17
plotting_figures <- TRUE #PARAm 18, arg 18
raster_overlap <- FALSE # PARAM 19, if TRUE, raster overlap is generated
#raster_pred <- FALSE # PARAM 20, if TRUE, raster prediction is generated
create_out_dir_param =T
#in_dir_shp <- file.path(in_dir_subset,"shapefiles")
in_dir_shp <- NULL # if NULL look in a predetermined place (see below)
#in_dir_shp <- /nobackupp6/aguzman4/climateLayers/tMinOut/reg*/subset/shapefiles/
  
### constant

pred_mod_name <- "mod1"
#pred_mod_name <- "mod1"
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
#mosaic_python <- "/nobackupp6/aguzman4/climateLayers/sharedCode/"
mosaic_python_script <- "/nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum.py"

#Product assessment
function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_06072017.R"
function_product_assessment_gap_detection_functions <- "gap_detection_product_assessment_functions_07042017d.R"

script_path <- "/nobackupp8/bparmen1/env_layers_scripts"
source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 
source(file.path(script_path,function_product_assessment_gap_detection_functions)) #source all functions used in this script 

debug(gap_tiles_assessment_fun)

#gap_tiles_assessment_fun(in_dir,y_var_name,region_name,num_cores,NA_flag_val,data_type_str,
#                                     shps_tiles,list_lf_raster_tif_tiles,infile_mask,countries_shp,
#                                     moscaic_python_script,out_dir,out_suffix)
  
gap_tiles_obj <- gap_tiles_assessment_fun(in_dir,y_var_name,region_name,num_cores,NA_flag_val,
                                     data_type_str,list_lf_raster_tif_tiles,
                                     infile_mask,countries_shp,moscaic_python_script,
                                     in_dir_shp,out_dir,out_suffix)
  
############################# END OF SCRIPT ###################################