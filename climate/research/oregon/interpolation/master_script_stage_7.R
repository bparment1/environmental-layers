##################    Master script for climate predictions  #######################################
############################ TMIN AND TMAX predictions ##########################################
#                           
##This script produces intperpolated surface of TMIN and TMAX for specified processing region(s) given sets 
#of inputs and parameters.
#STAGE 1: LST climatology downloading and/or calculation
#STAGE 2: Covariates preparation for study/processing area: calculation of covariates (spect,land cover,etc.) and reprojection
#STAGE 3: Data preparation: meteorological station database query and extraction of covariates values from raster brick
#STAGE 4: Raster prediction: run interpolation method (-- gam fusion, gam CAI, ...) and perform validation 
#STAGE 5: Output analyses: assessment of results for specific dates for individual tiles
#STAGE 6: Assessement of predictions by tiles, analyses and figures generation, summarized by regions 
#STAGE 7: Mosaicing of predicted surfaces and accuracy metrics (RMSE,MAE) by regions 
#STAGE 8: Comparison of predictions across regions and years with additional figures generation.

#AUTHOR: Benoit Parmentier                                                                        
#CREATED ON: 01/01/2016  
#MODIFIED ON: 02/09/2017
#PROJECT: NCEAS INPLANT: Environment and Organisms                                                                           

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh  
#MODULEPATH=$MODULEPATH:/nex/modules/files
#module load pythonkits/gdal_1.10.0_python_2.7.3_nex
#Set permissions:
#setfacl -Rm u:aguzman4:rwx /nobackupp6/aguzman4/climateLayers/LST_tempSpline/
#
## TODO:
# 
## Comments: Fixed accuracy bugs and tested command line script for jobs
#
## Commit: adding variable to deal with world mosaic

### Testing several years on the bridge before running jobs on nodes with qsub
#Use the following command to run as script via the shell on the bridge 
#run mosaic of predictions
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_06192016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 reg4_1999 /nobackupp8/bparmen1/climateLayers/out/reg4 TRUE 1999 6 1e+07 use_edge_weights rmse 19990107 19990108 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg4.tif /nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1999/df_assessment_files_reg4_1999_reg4_1999.txt python var_pred FALSE Int16 100 -100,100
#run mosaic for kriged residuals testing
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_06192016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 reg4_1999 /nobackupp8/bparmen1/climateLayers/out/reg4 TRUE 1999 6 1e+07 use_edge_weights rmse 19990107 19990108 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg4.tif /nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1999/df_assessment_files_reg4_1999_reg4_1999.txt python res_testing FALSE Int16 100 -100,100
#run mosaic for ac testing, rmse
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_06192016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 reg4_1999 /nobackupp8/bparmen1/climateLayers/out/reg4 TRUE 1999 6 1e+07 use_edge_weights rmse 19990107 19990108 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg4.tif /nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1999/df_assessment_files_reg4_1999_reg4_1999.txt python ac_testing FALSE Int16 100 -100,100
#run mosaic for ac testing, n
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_06192016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 reg4_1999 /nobackupp8/bparmen1/climateLayers/out/reg4 TRUE 1999 6 1e+07 use_edge_weights n 19990107 19990108 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg4.tif /nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1999/df_assessment_files_reg4_1999_reg4_1999.txt python ac_testing FALSE Int16 100 -100,100
#Region 5: test
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_07052016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg5 reg5_1991 /nobackupp6/aguzman4/climateLayers/out/reg5/mosaicsAc/ TRUE 1991 6 1e+07 use_edge_weights rmse 19910101 19910103 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg5.tif /nobackupp6/aguzman4/climateLayers/out/reg5/assessment/output_reg5_1991/df_assessment_files_reg5_1991_reg5_1991.txt python ac_testing FALSE Int16 100 -100,100

### Region 6
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_06192016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg6 reg6_2000 /nobackupp8/bparmen1/climateLayers/out/reg6/mosaics/ TRUE 2000 6 1e+07 use_edge_weights rmse 20000115 20000117 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg4.tif /nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1999/df_assessment_files_reg6_2000_reg6_2000.txt python var_pred FALSE Int16 100 -100,100

#Region 1: test

#Alberto run
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_04232016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg1 reg1_1984 /nobackupp6/aguzman4/climateLayers/out/reg1/mosaics/ TRUE 1984 40 1e+07 use_edge_weights rmse 19840101 19841231 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1.tif /nobackupp6/aguzman4/climateLayers/out/reg1/assessment//output_reg1_1984/df_assessment_files_reg1_1984_reg1_1984.txt python var_pred FALSE Int16 100 -100,100

#### Note it is mosaic not "mosaics"
#Run reg1 for reg1 and number of station n:
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_08232016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg1 reg1_1984 /nobackupp8/bparmen1/climateLayers/out/reg1/mosaic/ TRUE 1984 6 1e+07 use_edge_weights n 19840101 19840101 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1.tif /nobackupp6/aguzman4/climateLayers/out/reg1/assessment//output_reg1_1984/df_assessment_files_reg1_1984_reg1_1984.txt python ac_testing FALSE Int16 1 0,32767

#Run reg1 for reg1 and residuals testing:
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_08232016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg1 reg1_1984 /nobackupp8/bparmen1/climateLayers/out/reg1/mosaic/ TRUE 1984 6 1e+07 use_edge_weights rmse 19840101 19840101 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1.tif /nobackupp6/aguzman4/climateLayers/out/reg1/assessment//output_reg1_1984/df_assessment_files_reg1_1984_reg1_1984.txt python res_testing FALSE Int16 100 -100,100

#Run reg1 for reg1 and rmse testing:
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_08232016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg1 reg1_1984 /nobackupp8/bparmen1/climateLayers/out/reg1/mosaic/ TRUE 1984 6 1e+07 use_edge_weights rmse 19840101 19840101 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1.tif /nobackupp6/aguzman4/climateLayers/out/reg1/assessment//output_reg1_1984/df_assessment_files_reg1_1984_reg1_1984.txt python ac_testing FALSE Int16 100 -100,100

#Run reg1 for reg1 and  var_pred (tmax):
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_08232016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg1 reg1_1984 /nobackupp8/bparmen1/climateLayers/out/reg1/mosaic/ TRUE 1984 6 1e+07 use_edge_weights rmse 19840101 19840101 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1.tif /nobackupp6/aguzman4/climateLayers/out/reg1/assessment//output_reg1_1984/df_assessment_files_reg1_1984_reg1_1984.txt python var_pred FALSE Int16 100 -100,100

### For year 1992

#Run reg1 for reg1 and ac testing rmse:
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_09282016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg1 reg1_1992 /nobackupp8/bparmen1/climateLayers/out/reg1/mosaic/ TRUE 1992 6 1e+07 use_edge_weights rmse 19920101 19920101 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1.tif /nobackupp6/aguzman4/climateLayers/out/reg1/assessment//output_reg1_1992/df_assessment_files_reg1_1992_reg1_1992.txt python ac_testing FALSE Int16 100 -100,100

#Run reg1 for reg1 and var_pred (tmax):
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_09282016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg1 reg1_1992 /nobackupp8/bparmen1/climateLayers/out/reg1/mosaic/ TRUE 1992 6 1e+07 use_edge_weights rmse 19920101 19920101 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1.tif /nobackupp6/aguzman4/climateLayers/out/reg1/assessment//output_reg1_1992/df_assessment_files_reg1_1992_reg1_1992.txt python var_pred FALSE Int16 100 -100,100

##################################################################################################

### PARAMETERS DEFINED IN THE SCRIPT


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
function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_02012017.R"
#function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_07052016.R" #PARAM12
function_mosaicing <-"global_run_scalingup_mosaicing_02092017.R"
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
stages_to_run<-c(0,0,0,0,0,0,7) #this is stage, other stages are stored in files.#CONSTANT1
  
#CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84, CONST2
CRS_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84, CONST2

#####mosaicing parameters

#Data is on ATLAS or NASA NEX

### PARAMETERS DEFINED IN THE SCRIPT
#There are 33 parameters, 3 constants and 20 arguments (drawn from the parameters) for the Rscript shell call.
#The 20 arguments are passed directly from Rscript:

#/nobackupp6/aguzman4/climateLayers/out/mosaics/reg4/mosaic/output_reg4_1984

## Use the following values to run code from the shell:
var <- "TMAX" # variable being interpolated #param 1, arg 1
in_dir <- "/nobackupp6/aguzman4/climateLayers/out/" #PARAM2,arg 2
region_name <- "world"
#region_name <- "reg1" #PARAM 3, arg 3 #reg4 South America, Africa reg5,Europe reg2, North America reg1, Asia reg3
#if region name is world then run full mosaic with data available for the corresponding year
#out_suffix <- "reg1_1992" #PARAM 4, arg 4
out_suffix <- "world_1984"
out_suffix_str <- region_name #PARAM 4, CONST 3
out_dir <- "/nobackupp8/bparmen1/climateLayers/out/world" #PARAM 5,arg 5 use this location for now
#out_dir <- "/nobackupp8/bparmen1/climateLayers/out/reg5/mosaicsAc" #PARAM 5,arg 5 use this location for now
# out_dir <- "/nobackupp8/bparmen1/climateLayers/out/reg1/mosaic"
create_out_dir_param <- TRUE #PARAM 6, arg 6
year_predicted <- 1984 #PARAM 7, arg 7
num_cores <- 6 #PARAM 8, arg 8
max_mem = 1e+07 #param 9, arg 9
mosaicing_method <- "use_edge_weights" #PARAM10, arg 10
metric_name <- "rmse" # "mae", "r" for MAE, R etc.; can also be ns or nv? #PARAM 11, arg 11
# #metric_name <- "n"
# #metric_name <- "mae"
# 
day_start <- "19840701" #PARAM 12 arg 12
day_end <- "19840705" #PARAM 13 arg 13
# #day_start <- "19920102" #PARAM 12 arg 12
# #day_end <- "19920104" #PARAM 13 arg 13
# 
# #infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg5.tif" #PARAM 14, arg 14
#infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1.tif" #PARAM 14, arg 14
infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1_reg4.tif" #PARAM 14, arg 14
#infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_world.tif" #PARAM 14, arg 14

df_assessment_files_name <- NULL
# #df_assessment_files_name <- "/nobackupp6/aguzman4/climateLayers/out/reg5/assessment/output_reg5_1991/df_assessment_files_reg5_1991_reg5_1991.txt"  # data.frame with all files used in assessmnet, PARAM 15
# #df_assessment_files_name <- "/nobackupp6/aguzman4/climateLayers/out/reg5/assessment/output_reg5_1985/df_assessment_files_reg5_1985_reg5_1985.txt"
algorithm <- "python" #PARAM 16 #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
layers_option <- c("var_pred") #arg 17 ,param 17, options are:#res_training, res_testing,ac_training, ac_testing, var_pred
# #layers_option <- c("ac_training") #arg 17 ,param 17, options are:#res_training, res_testing,ac_training, ac_testing, var_pred
# #layers_option <- c("res_training") # #arg 17 ,param 17, options are:#res_training, res_testing,ac_training, ac_testing, var_pred
# #layers_option <- c("res_testing") #arg 17 ,param 17, options are:#res_training, res_testing,ac_training, ac_testing, var_pred
# #layers_option <- c("ac_testing") #arg 17 ,param 17, options are:#res_training, res_testing,ac_training, ac_testing, var_pred
# 
tmp_files <- FALSE #arg 18, param 18, keep temp files if TRUE
data_type <- "Int16" #, param 19, use int32 for output layers mosaiced
scaling <- 100 #, param 20, if null use 1
# #scaling <- 1 #use this if predicting n rather than other variables
values_range <- "-100,100"
#this could be a list of folder or file with location of region mosaics to list...
#list_reg <- "reg1,reg4" # if NULL then use other information, use this if using world mosaicing #param 22
#infile_reg_mosaics <- "/data/project/layers/commons/NEX_data/regions_input_files/world_input_mosaics_tmax_var_02102017.csv"
infile_reg_mosaics <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/world_input_mosaics_tmax_var_02102017.csv"

NA_flag_val <- -32768 #should be here

# values_range <- "0,32767" #this is for n variable

#path_assessment <- NOT USED "/nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1991" #PARAM 14a, arg 14


############################

var <- args[1] # variable being interpolated #param 1, arg 1
#var<-"TMAX" # variable being interpolated #param 1, arg 1

#PARAM 2
#in_dir <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_12072015" #NCEAS
#in_dir <- "/nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_12072015" #NEX
#in_dir <- "/nobackupp6/aguzman4/climateLayers/out/" #PARAM2
in_dir <- args[2] #PARAM2

region_name <- args[3] #PARAM3
#region_name <- "reg4" #PARAM 3 #reg4 South America, Africa reg5,Europe reg2, North America reg1, Asia reg3

#out_suffix <- paste(region_name,"_","run10_1500x4500_global_analyses_pred_1991_04052016",sep="") #PARAM 6
#out_suffix_str <- "run10_1500x4500_global_analyses_pred_1991_04052016" #PARAM 7

out_suffix <- args[4] #PARAM 4
#out_suffix <- region_name #PARAM 4
out_suffix_str <- region_name #PARAM 4, CONST 3
out_dir <- args[5] #PARAM 5
#out_dir <- "/nobackupp8/bparmen1/climateLayers/out/reg4" #PARAM 5, use this location for now
create_out_dir_param <- args[6] #PARAM 6
#create_out_dir_param <- TRUE #PARAM 6

year_predicted <- args[7] #PARAM 7
#year_predicted <- 1991 #PARAM 7

num_cores <- args[8] #PARAM 8
#num_cores <- 6 #PARAM 8         

#max number of cells to read in memory
max_mem<-args[9] #PARAM 9

#mosaicing_method <- c("unweighted","use_edge_weights") #PARAM10
mosaicing_method <- args[10] #PARAM10

metric_name <- args[11]
#metric_name <- "rmse" #RMSE, MAE etc. #PARAM 11
#if daily mosaics NULL then mosaicas all days of the year #PARAM 12
#day_to_mosaic_range <- c("19910101","19910103") #if null run all year #PARAM 12
#day_to_mosaic_range <- c("19910101","19910101") #if null run all year #PARAM 12
#day_to_mosaic_range <- args[12] #PARAM 12
day_start <- args[12] #PARAM 12
day_end <- args[13] #PARAM 13

###Separate folder for masks by regions, should be listed as just the dir!!... #PARAM 14
#infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_reg4.tif" #PARAM 14
infile_mask <- args[14]
#infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_reg4.tif"
## All of this is interesting so use df_assessment!!

#path_assessment <- "/nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1991"
#path_assessment <- file.path(in_dir,region_name,"assessment",paste("output_",region_name,year_processed,sep=""))

df_assessment_files_name <- args[15] #PARAM 15, files containing assessment information
#df_assessment_files_name <- "/nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1991/df_assessment_files_reg4_1991_reg4_1991.txt" # data.frame with all files used in assessmnet, PARAM 14
algorithm <- args[16] #PARAM 16
#algorithm <- "python" #PARAM 15 #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
#algorithm <- "R" #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann

layers_option <- args[17] # PARAM 17 options are:
#layers_option <- c("var_pred") #options are:
#res_training, res_testing,ac_training, ac_testing, var_pred
tmp_files <- args[18] #PARAM 18
data_type <- args[19] #PARAM 19 #use integer 32 for layers outputs
scaling <- args[20] #PARAM 19 #use integer 32 for layers outputs
values_range <-  args[21] #Param 21, args 21
infile_reg_mosaics  <- args[22]


#tmp_files <- FALSE
interpolation_method <- c("gam_CAI") #PARAM21
pred_mod_name <- "mod1" #PARAM 22
var_pred <- "res_mod1" #used in residuals mapping #PARAM 23
proj_str<- CRS_WGS84 #PARAM 24 #check this parameter
file_format <- ".tif" #PARAM 25
#NA_value <- -9999 #PARAM 26
NA_value <- -32768 #PARAM 26
NA_flag_val <- NA_value #PARAM 26
use_autokrige <- F #PARAM 28
#proj_str <- CRS_locs_WGS84 #PARAM 29
#python script and gdal on NEX NASA:
mosaic_python <- "/nobackupp6/aguzman4/climateLayers/sharedCode/" #PARAM 29
python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin" #PARAM 30
#python script and gdal on Atlas NCEAS
#mosaic_python <- "/data/project/layers/commons/NEX_data/sharedCode" #PARAM 30
#python_bin <- "/usr/bin" #PARAM 30 #Atlas
match_extent <- "FALSE" #PARAM 31 #try without matching!!!
#for residuals...
list_models <- NULL #PARAM 32
#list_models <- paste(var_pred,"~","1",sep=" ") #if null then this is the default...
  
if (var == "TMAX") {
  y_var_name <- "dailyTmax"
  y_var_month <- "TMax"
}
if (var == "TMIN") {
  y_var_name <- "dailyTmin"
  y_var_month <- "TMin"
}

if(!(is.null(day_start)) & !(is.null(day_end))){
  day_to_mosaic_range <- c(day_start,day_end) #if null run all year #PARAM 12
}else{
  day_to_mosaic_range <- NULL
}

#parse input value range
values_range <- as.numeric(unlist(strsplit(values_range,",")))
scaling <- as.numeric(scaling)

#browser()

#rasterOptions(maxmemory=1e+07,timer=TRUE)
##30 parameters passed
list_param_run_mosaicing_prediction <- list(in_dir,y_var_name,interpolation_method,region_name,
                 mosaicing_method,out_suffix,out_suffix_str,metric_name,pred_mod_name,var_pred, out_dir,
                 create_out_dir_param,day_to_mosaic_range,year_predicted,proj_str,file_format,NA_value,
                 num_cores,use_autokrige,infile_mask,df_assessment_files_name,mosaic_python,
                 python_bin,algorithm,match_extent,list_models,layers_option,tmp_files,data_type,scaling,
                 values_range,infile_reg_mosaics)
param_names <- c("in_dir","y_var_name","interpolation_method","region_name",
                 "mosaicing_method","out_suffix","out_suffix_str","metric_name","pred_mod_name","var_pred","out_dir",
                 "create_out_dir_param","day_to_mosaic_range","year_predicted","proj_str","file_format","NA_value",
                 "num_cores","use_autokrige","infile_mask","df_assessment_files_name","mosaic_python",
                 "python_bin","algorithm","match_extent","list_models","layers_option","tmp_files","data_type","scaling",
                 "values_range","infile_reg_mosaics")
names(list_param_run_mosaicing_prediction) <- param_names
#list_param_run_mosaicing_prediction
#debug(run_mosaicing_prediction_fun)
#debug(debug_fun_test)
#debug_fun_test(list_param_raster_prediction)
i <- 1 #this select the first year of list_year_predicted
#function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_07012016.R" #PARAM12
#function_mosaicing <-"global_run_scalingup_mosaicing_02012017.R"
#source(file.path(script_path,function_mosaicing)) #source all functions used in this script 
#source(file.path(script_path,function_mosaicing_functions)) #source all functions used in this script 

if (stages_to_run[7]==7){
  assessment_prediction_obj <- run_mosaicing_prediction_fun(i,list_param_run_mosaicing_prediction)
  #list_param_run_mosaicing_prediction
}

#runs in 42 minutes for 3 dates but note that beyond date 1, the process is about 11 minutes or so.
#19h09
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

