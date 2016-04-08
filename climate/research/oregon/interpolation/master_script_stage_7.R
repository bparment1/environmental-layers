##################    Master script for climate predictions  #######################################
############################ TMIN AND TMAX predictions ##########################################
#                           
##This script produces intperpolated surface of TMIN and TMAX for specified processing region(s) given sets 
#of inputs and parameters.
#STAGE 1: LST climatology downloading and/or calculation
#STAGE 2: Covariates preparation for study/processing area: calculation of covariates (spect,land cover,etc.) and reprojection
#STAGE 3: Data preparation: meteorological station database query and extraction of covariates values from raster brick
#STAGE 4: Raster prediction: run interpolation method (-- gam fusion, gam CAI, ...) and perform validation 
#STAGE 5: Output analyses: assessment of results for specific dates and tiles
#STAGE 6: Assessement of predictions by tiles and regions 
#STAGE 7: Mosaicing of predicted surfaces and accuracy metrics (RMSE,MAE) by regions 
#STAGE 8: Comparison of predictions across regions and years with figures generation.

#AUTHOR: Benoit Parmentier                                                                        
#CREATED ON: 01/01/2016  
#MODIFIED ON: 04/08/2016  
#PROJECT: NCEAS INPLANT: Environment and Organisms                                                                           

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh  
#MODULEPATH=$MODULEPATH:/nex/modules/files
#module load pythonkits/gdal_1.10.0_python_2.7.3_nex

## TODO:
# 
# Clean up temporary files
# 
##################################################################################################

###Loading R library and packages  ou 
library(RPostgreSQL)
library(maps)
library(maptools)
library(parallel)
library(gtools)                              # loading some useful tools 
library(mgcv)                                # GAM package by Simon Wood
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gstat)                               # Kriging and co-kriging by Pebesma et al.
library(fields)                              # NCAR Spatial Interpolation methods such as kriging, splines
library(raster)                              # Hijmans et al. package for raster processing
library(rasterVis)
library(spgwr)
library(reshape)
library(plotrix)

######## PARAMETERS FOR WORK FLOW #########################
### Need to add documentation ###

#Adding command line arguments to use mpiexec
args<-commandArgs(TRUE)
#script_path<-"/nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/oregon/interpolation"
#dataHome<-"/nobackupp6/aguzman4/climateLayers/interp/testdata/"
#script_path2<-"/nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation"

#CALLED FROM MASTER SCRIPT:

#script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts"
script_path <- "/nobackupp8/bparmen1/env_layers_scripts" #path to script
function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_04102016.R" #PARAM12
function_mosaicing <-"global_run_scalingup_mosaicing_04082016.R"
source(file.path(script_path,function_mosaicing)) #source all functions used in this script 
source(file.path(script_path,function_mosaicing_functions)) #source all functions used in this script 

function_assessment_part1_functions <- "global_run_scalingup_assessment_part1_functions_02112015.R" #PARAM12
function_assessment_part1a <-"global_run_scalingup_assessment_part1a_01042016.R"
function_assessment_part2 <- "global_run_scalingup_assessment_part2_02092016.R"
function_assessment_part2_functions <- "global_run_scalingup_assessment_part2_functions_01032016.R"
function_assessment_part3 <- "global_run_scalingup_assessment_part3_02102016.R"
source(file.path(script_path,function_assessment_part1_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part1a)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part3)) #source all functions used in this script 

### Parameters, constants and arguments ###
stages_to_run<-c(0,0,0,0,0,0,7) #this is stage, other stages are stored in files.
  
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

#####mosaicing parameters

#Data is on ATLAS or NASA NEX

var <- args[1] # variable being interpolated #param 1, arg 1
var<-"TMAX" # variable being interpolated #param 1, arg 1

var<-"TMAX" # variable being interpolated
if (var == "TMAX") {
  y_var_name <- "dailyTmax"
  y_var_month <- "TMax"
}
if (var == "TMIN") {
  y_var_name <- "dailyTmin"
  y_var_month <- "TMin"
}

#PARAM 2
#in_dir <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_12072015" #NCEAS
#in_dir <- "/nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_12072015" #NEX
in_dir <- "/nobackupp6/aguzman4/climateLayers/out/"

interpolation_method <- c("gam_CAI") #PARAM3

region_name <- "reg4" #PARAM 4 #reg4 South America, Africa reg5,Europe reg2, North America reg1, Asia reg3
mosaicing_method <- c("unweighted","use_edge_weights") #PARAM5
#out_suffix <- paste(region_name,"_","run10_1500x4500_global_analyses_pred_1991_04052016",sep="") #PARAM 6
#out_suffix_str <- "run10_1500x4500_global_analyses_pred_1991_04052016" #PARAM 7

out_suffix <- region_name #PARAM 6
out_suffix_str <- region_name #PARAM 7

metric_name <- "rmse" #RMSE, MAE etc. #PARAM 8
pred_mod_name <- "mod1" #PARAM 9
var_pred <- "res_mod1" #used in residuals mapping #PARAM 10

#out_dir <- in_dir #PARAM 11
out_dir <- "/nobackupp8/bparmen1/climateLayers/out/reg4" #PARAM 11, use this location for now
create_out_dir_param <- TRUE #PARAM 12

#if daily mosaics NULL then mosaicas all days of the year #PARAM 13
#day_to_mosaic <- c("19910101","19910102","19910103") #,"19920104","19920105") #PARAM9, two dates note in /tiles for now on NEX
day_to_mosaic_range <- c("19910101","19910103") #if null run all year
#CRS_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84 #CONSTANT1
#CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
#proj_str<- CRS_WGS84 #PARAM 8 #check this parameter
 
file_format <- ".tif" #PARAM 14
NA_value <- -9999 #PARAM 15
NA_flag_val <- NA_value #PARAM 16
     
num_cores <- 6 #PARAM 17                 
#region_names <- c("reg23","reg4") #selected region names, ##PARAM 18 
use_autokrige <- F #PARAM 19
proj_str <- CRS_locs_WGS84

###Separate folder for masks by regions, should be listed as just the dir!!... #PARAM 20
infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_reg4.tif"
#infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_reg4.tif"
## All of this is interesting so use df_assessment!!

year_processed <- 1991 #PARAM 31
#path_assessment <- "/nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1991"
#path_assessment <- file.path(in_dir,region_name,"assessment",paste("output_",region_name,year_processed,sep=""))
df_assessment_files_name <- "/nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1991/df_assessment_files_reg4_1991_reg4_1991.txt" # data.frame with all files used in assessmnet, PARAM 21

#in_dir can be on NEX or Atlas

#python script and gdal on NEX NASA:
mosaic_python <- "/nobackupp6/aguzman4/climateLayers/sharedCode/"
python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin"
#python script and gdal on Atlas NCEAS
#mosaic_python <- "/data/project/layers/commons/NEX_data/sharedCode" #PARAM 26
#python_bin <- "/usr/bin" #PARAM 27

algorithm <- "python" #PARAM 28 #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
#algorithm <- "R" #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
match_extent <- "FALSE" #PARAM 29 #try without matching!!!

#for residuals...
list_models <- NULL #PARAM 30
#list_models <- paste(var_pred,"~","1",sep=" ") #if null then this is the default...
  
#max number of cells to read in memory
max_mem<-args[11]
#in_dir_tiles <- file.path(in_dir,"tiles") #this is valid both for Atlas and NEX
layers_option <- c("var_pred") #options are:
#res_training, res_testing,ac_training, ac_testing, var_pred
tmp_files <- FALSE

#rasterOptions(maxmemory=1e+07,timer=TRUE)
list_param_run_mosaicing_prediction <- list(in_dir,y_var_name,interpolation_method,region_name,
                 mosaicing_method,out_suffix,out_suffix_str,metric_name,pred_mod_name,var_pred,
                 create_out_dir_param,day_to_mosaic_range,proj_str,file_format,NA_value,num_cores,
                 use_autokrige,infile_mask,df_assessment_files_name,mosaic_python,
                 python_bin,algorithm,match_extent,list_models,layers_option,tmp_files)
param_names <- c("in_dir","y_var_name","interpolation_method","region_name",
                 "mosaicing_method","out_suffix","out_suffix_str","metric_name","pred_mod_name","var_pred",
                 "create_out_dir_param","day_to_mosaic_range","proj_str","file_format","NA_value","num_cores",
                 "use_autokrige","infile_mask","df_assessment_files_name","mosaic_python",
                 "python_bin","algorithm","match_extent","list_models","layers_option","tmp_files")
names(list_param_run_mosaicing_prediction) <- param_names
#list_param_run_mosaicing_prediction
#debug(run_mosaicing_prediction_fun)
#debug(debug_fun_test)
#debug_fun_test(list_param_raster_prediction)
i <- 1 #this select the first year of list_year_predicted
#function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_04072016.R" #PARAM12
#function_mosaicing <-"global_run_scalingup_mosaicing_04072016.R"
#source(file.path(script_path,function_mosaicing)) #source all functions used in this script 
#source(file.path(script_path,function_mosaicing_functions)) #source all functions used in this script 

if (stages_to_run[7]==7){
  assessment_prediction_obj <- run_mosaicing_prediction_fun(i,list_param_run_mosaicing_prediction)
  #list_param_run_mosaicing_prediction
}

###############   END OF SCRIPT   ###################
#####################################################

# #LAND COVER INFORMATION
# LC1: Evergreen/deciduous needleleaf trees
# LC2: Evergreen broadleaf trees
# LC3: Deciduous broadleaf trees
# LC4: Mixed/other trees
# LC5: Shrubs
# LC6: Herbaceous vegetation
# LC7: Cultivated and managed vegetation
# LC8: Regularly flooded shrub/herbaceous vegetation
# LC9: Urban/built-up
# LC10: Snow/ice
# LC11: Barren lands/sparse vegetation
# LC12: Open water

