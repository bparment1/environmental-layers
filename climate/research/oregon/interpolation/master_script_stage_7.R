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
#MODIFIED ON: 04/11/2016  
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
function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_04102016.R" #PARAM12
function_mosaicing <-"global_run_scalingup_mosaicing_04112016.R"
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

### PARAMETERS DEFINED IN THE SCRIPT
#There are 31 parameters, 1 constant and 17 arguments (drawn from the parameters) for the Rscript call.
#The arguments are passed directly from Rscript:
#var <- args[1] # variable being interpolated #param 1, arg 1
#in_dir <- args[2] # This is the output directory containing global prediction e.g./nobackupp6/aguzman4/climateLayers/out/ param 5, arg 2
#region_name <- args[3] # region e.g. "reg4" param 3, arg 3
#out_suffix <- args[4] # formely out_prefix, this is used in creating an output directory, it is suggested to use "reg4" or same as region_name
#out_dir <- args[5] # output parent dir, can be home dir or other, param 5, arg 5
#create_out_dir_param <- args[6] # if TRUE create a output from "output"+out_prefix param 6, arg 6
#year_predicted <- args[7] # enter as list but currently runs on the first element of the list, param 7, arg 7
#num_cores <- args[8] #number of cores used # param 8, arg 8
#max_mem <- args[9] # maximum memory, param 9
#mosaicing_method <- arg[10] #PARAM10
#metric_name <- arg[11] #"rmse" #RMSE, MAE etc. #PARAM 11
#day_to_mosaic_range <- arg[12] #c("19910101","19910103") #if null run all year, param 12
#infile_mask <- arg[13] # "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_reg4.tif"
#df_assessment_files_name <- arg[14] #"/nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1991/df_assessment_files_reg4_1991_reg4_1991.txt" # data.frame with all files used in assessmnet, PARAM 21
#algorithm <- arg[15] #"python" #PARAM 28 #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
#layers_option <- arg[16] #c("var_pred") #options are:res_training, res_testing,ac_training, ac_testing, var_pred
#tmp_files <- arg[17] #FALSE

### Use the following values to run code from the shell:
#var<-"TMAX" # variable being interpolated #param 1, arg 1
#in_dir <- "/nobackupp6/aguzman4/climateLayers/out/" #PARAM2,arg 2
#region_name <- "reg4" #PARAM 3, arg 3 #reg4 South America, Africa reg5,Europe reg2, North America reg1, Asia reg3
#out_suffix <- "reg4" #PARAM 4, arg 4
#out_dir <- "/nobackupp8/bparmen1/climateLayers/out/reg4" #PARAM 5,arg 5 use this location for now
#create_out_dir_param <- TRUE #PARAM 6, arg 6
#year_predicted <- 1991 #PARAM 7, arg 7
#num_cores <- 6 #PARAM 8, arg 8         
#max_mem = 1e+07 #param 9, arg 9
#mosaicing_method <- use_edge_weights" #PARAM10, arg 10
#metric_name <- "rmse" #RMSE, MAE etc. #PARAM 11, arg 11
#day_start <- "19910101" #PARAM 12
#day_end <- "19910101" #PARAM 13
#infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_reg4.tif" #PARAM 14, arg 14
#df_assessment_files_name <- "/nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1991/df_assessment_files_reg4_1991_reg4_1991.txt"  # data.frame with all files used in assessmnet, PARAM 15
#algorithm <- "python" #PARAM 16 #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
#layers_option <- c("var_pred") #arg 17 ,param 17, options are:#res_training, res_testing,ac_training, ac_testing, var_pred
#tmp_files <- FALSE #arg 18, param 18

#path_assessment <- NOT USED "/nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1991" #PARAM 14a, arg 14

### Testing several years on the bridge before running jobs on nodes with qsub
#Use the following command to run as script via the shell on the bridge 
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_7_04112016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 reg4 /nobackupp8/bparmen1/climateLayers/out/reg4 TRUE 1991 6 1e+07 use_edge_weights rmse 19910101 19910103 /nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_reg4.tif /nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1991/df_assessment_files_reg4_1991_reg4_1991.txt python var_pred FALSE 

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
out_suffix_str <- region_name #PARAM 4
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

df_assessment_files_name <- args[15] #PARAM 15
#df_assessment_files_name <- "/nobackupp6/aguzman4/climateLayers/out/reg4/assessment/output_reg4_1991/df_assessment_files_reg4_1991_reg4_1991.txt" # data.frame with all files used in assessmnet, PARAM 14
algorithm <- args[16] #PARAM 16
#algorithm <- "python" #PARAM 15 #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
#algorithm <- "R" #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann

layers_option <- args[17] # PARAM 17 options are:
#layers_option <- c("var_pred") #options are:
#res_training, res_testing,ac_training, ac_testing, var_pred
tmp_files <- args[18] #PARAM 18
#tmp_files <- FALSE
interpolation_method <- c("gam_CAI") #PARAM19
pred_mod_name <- "mod1" #PARAM 20
var_pred <- "res_mod1" #used in residuals mapping #PARAM 21
proj_str<- CRS_WGS84 #PARAM 22 #check this parameter
file_format <- ".tif" #PARAM 23
NA_value <- -9999 #PARAM 24
NA_flag_val <- NA_value #PARAM 24
use_autokrige <- F #PARAM 25
proj_str <- CRS_locs_WGS84 #PARAM 26
#python script and gdal on NEX NASA:
mosaic_python <- "/nobackupp6/aguzman4/climateLayers/sharedCode/" #PARAM 27
python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin" #PARAM 28
#python script and gdal on Atlas NCEAS
#mosaic_python <- "/data/project/layers/commons/NEX_data/sharedCode" #PARAM 29
#python_bin <- "/usr/bin" #PARAM 30
match_extent <- "FALSE" #PARAM 29 #try without matching!!!
#for residuals...
list_models <- NULL #PARAM 30
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

#browser()

#rasterOptions(maxmemory=1e+07,timer=TRUE)
list_param_run_mosaicing_prediction <- list(in_dir,y_var_name,interpolation_method,region_name,
                 mosaicing_method,out_suffix,out_suffix_str,metric_name,pred_mod_name,var_pred, out_dir,
                 create_out_dir_param,day_to_mosaic_range,year_predicted,proj_str,file_format,NA_value,num_cores,
                 use_autokrige,infile_mask,df_assessment_files_name,mosaic_python,
                 python_bin,algorithm,match_extent,list_models,layers_option,tmp_files)
param_names <- c("in_dir","y_var_name","interpolation_method","region_name",
                 "mosaicing_method","out_suffix","out_suffix_str","metric_name","pred_mod_name","var_pred","out_dir",
                 "create_out_dir_param","day_to_mosaic_range","year_predicted","proj_str","file_format","NA_value","num_cores",
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

#runs in 42 minutes for 3 dates but note that beyond date 1, the process is about 11 minutes or so.

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

