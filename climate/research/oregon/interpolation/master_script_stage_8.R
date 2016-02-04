##################    Master script for climate predictions  #######################################
############################ TMIN AND TMAX predictions ##########################################
#                           
##This script produces intperpolated surface of TMIN and TMAX for specified processing region(s) given sets 
#of inputs and parameters.
#STAGE 1: LST climatology downloading and/or calculation
#STAGE 2: Covariates preparation for study/processing area: calculation of covariates (spect,land cover,etc.) and reprojection
#STAGE 3: Data preparation: meteorological station database query and extraction of covariates values from raster brick
#STAGE 4: Raster prediction: run interpolation method (-- gam fusion, gam CAI, ...) and perform validation 
#STAGE 5: Output analyses: assessment of results for specific dates... (tile based)
#STAGE 6: Assessement of predictions by tiles and region: summary tables and figures 
#STAGE 7: Mosaicing of predictions and accuracy layer productions
#STAGE 8: Comparison of predictions across regions and years with figure generation.
#AUTHOR: Benoit Parmentier                                                                        
#CREATED ON: 12/29/2015  
#MODIFIED ON: 02/04/2016  
#PROJECT: NCEAS-IPLANT-NASA: Environment Layers                                                                           

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh  
#MODULEPATH=$MODULEPATH:/nex/modules/files
#module load pythonkits/gdal_1.10.0_python_2.7.3_nex

## TODO:
# Add  assessment part 2 (figures): still testing
# Make this callable from the shell
# Adapt for python script
# Fix figure for part2
# Call also use library(optparse)

##################################################################################################
#
### PARAMETERS DEFINED IN THE SCRIPT
#There are 21 parameters, 1 constant and 8 arguments (drawn from the parameters) for the Rscript call.
#The arguments are passed directly from Rscript:
#var <- args[1] # variable being interpolated #param 1, arg 1
#in_dir1 <- args[2] # This is the output directory containing global prediction e.g./nobackupp6/aguzman4/climateLayers/out/ param 5, arg 2
#region_name <- args[3] # region e.g. "reg4" param 6, arg 3
#out_prefix <- args[4] # this is used in creating an output directory,include region name? param 7, arg 4
#out_dir <- args[5] # output parent dir, can be home dir or other, param 8, arg 5)
#create_out_dir_param <- args[6] # if TRUE create a output from "output"+out_prefix param 9, arg 6
#list_year_predicted <- args[7] # enter as list but currently runs on the first element of the list, param 10, arg 7
#num_cores <- args[8] #number of cores used # param 13, arg 8
#max_mem <- args[9] # maximum memory, param 21

#var = "TMAX" # variable being interpolated #param 1, arg 1
#in_dir1 = "/nobackupp6/aguzman4/climateLayers/out/" #param 5, arg 2
#region_name = "reg4" #param 6, arg 3
#out_prefix = "run_global_analyses_pred_12282015" #param 7, arg 4
#out_dir = "/nobackupp8/bparmen1/" #param 8, arg 5
#create_out_dir_param = "TRUE" #param 9, arg 6
#list_year_predicted = c(2010) # param 10, arg 7
#num_cores = 6 #number of cores used # param 13, arg 8
#max_mem = 1e+07 #param 21, arg 9

### Testing several years on the bridge before running jobs on nodes with qsub
#This can be made in a data.frame to run through...  
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_6_01222016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 run_global_analyses_pred_12282015 /nobackupp8/bparmen1/ TRUE 2010 6 1e+07
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_6_01222016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 run_global_analyses_pred_2011_reg4 /nobackupp8/bparmen1/ TRUE 2011 6 1e+07
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_6_01222016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 run_global_analyses_pred_2012_reg4 /nobackupp8/bparmen1/ TRUE 2012 6 1e+07
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_6_01222016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 run_global_analyses_pred_2013_reg4 /nobackupp8/bparmen1/ TRUE 2013 6 1e+07
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_6_01222016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 run_global_analyses_pred_2014_reg4 /nobackupp8/bparmen1/ TRUE 2014 6 1e+07
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_6_01222016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 run_global_analyses_pred_2009_reg4 /nobackupp8/bparmen1/ TRUE 2009 6 1e+07
#Rscript /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_6_01222016.R TMAX /nobackupp6/aguzman4/climateLayers/out/ reg4 run_global_analyses_pred_2010_reg4 /nobackupp8/bparmen1/ TRUE 2010 6 1e+07

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
args <- commandArgs(TRUE)
#script_path<-"/nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/oregon/interpolation"
#dataHome<-"/nobackupp6/aguzman4/climateLayers/interp/testdata/"
#script_path2<-"/nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation"
#var <- args[1] # variable being interpolated #param 1, arg 1
#in_dir1 <- args[2] #param 5, arg 2
#region_name <- args[3] #param 6, arg 3
#out_prefix <- args[4] #param 7, arg 4
#out_dir <- args[5] #param 8, arg 5
#out_dir <-paste(out_dir,"_",out_prefix,sep="")
#create_out_dir_param <- args[6] #param 9, arg 6
#list_year_predicted <- args[7] # param 10, arg 7
#num_cores <- args[8] #number of cores used # param 13, arg 8
#max_mem <- args[9] #param 21

#CALLED FROM MASTER SCRIPT:

script_path <- "/nobackupp8/bparmen1/env_layers_scripts" #path to script
function_assessment_part1_functions <- "global_run_scalingup_assessment_part1_functions_02112015.R" #PARAM12
function_assessment_part1a <-"global_run_scalingup_assessment_part1a_01042016.R"
function_assessment_part2 <- "global_run_scalingup_assessment_part2_02032016.R"
function_assessment_part2_functions <- "global_run_scalingup_assessment_part2_functions_01032016.R"
source(file.path(script_path,function_assessment_part1_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part1a)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2_functions)) #source all functions used in this script 

### Parameters, constants and arguments ###
stages_to_run<-c(0,0,0,0,0,6) #this stage 2 to 5 currently stored in another file.

CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #constant 1

#var<-"TMAX" # variable being interpolated #param 1, arg 1
var <- args[1] # variable being interpolated #param 1, arg 1

##Add for precip later...
if (var == "TMAX") {
  y_var_name <- "dailyTmax"
  y_var_month <- "TMax"
}
if (var == "TMIN") {
  y_var_name <- "dailyTmin"
  y_var_month <- "TMin"
}

#interpolation_method<-c("gam_fusion") #other otpions to be added later
interpolation_method<-c("gam_CAI") #param 2
CRS_interp <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #param 3
#CRS_interp <-"+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs";

out_region_name<-""
list_models<-c("y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)") #param 4

#reg1 (North Am), reg2(Europe),reg3(Asia), reg4 (South Am), reg5 (Africa), reg6 (Australia-Asia)
#master directory containing the definition of tile size and tiles predicted
#in_dir1 <- "/nobackupp6/aguzman4/climateLayers/out/" #param 5, arg 2
in_dir1 <- args[2] #param 5, arg 2

#region_names <- c("reg1","reg2","reg3","reg4","reg5","reg6") #selected region names, #PARAM2
#region_names <- c("reg23","reg4") #selected region names,
#run assessment by region, this is a unique region only 
#region_name <- c("reg4") #param 6, arg 3
region_name <- args[3] #param 6, arg 3

#out_prefix <- "run_global_analyses_pred_12282015" #param 7, arg 4
#out_dir <- "/nobackupp8/bparmen1/" #param 8, arg 5
#out_dir <-paste(out_dir,"_",out_prefix,sep="")
create_out_dir_param <- TRUE #param 9, arg 6
out_prefix <- args[4] #param 7, arg 4
out_dir <- args[5] #param 8, arg 5
#out_dir <-paste(out_dir,"_",out_prefix,sep="")
create_out_dir_param <- args[6] #param 9, arg 6

#list_year_predicted <- 1984:2004
#list_year_predicted <- c("2014") # param 10, arg 7
#year_predicted <- list_year_predicted[1]
list_year_predicted <- args[7] # param 10, arg 7

file_format <- ".tif" #format for mosaiced files # param 11
NA_flag_val <- -9999  #No data value, # param 12
#num_cores <- 6 #number of cores used # param 13, arg 8
plotting_figures <- TRUE #running part2 of assessment to generate figures... # param 14
num_cores <- args[8] #number of cores used # param 13, arg 8
  
##Additional parameters used in part 2, some these may be removed as code is simplified
mosaic_plot <- FALSE #param 15
day_to_mosaic <- c("19920102","19920103","19920103") #param 16, not in use...
multiple_region <- TRUE #param 17
countries_shp <- "/nobackupp8/bparmen1/NEX_data/countries.shp" #param 18
#countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #Atlas
plot_region <- TRUE  # param 19
threshold_missing_day <- c(367,365,300,200) # param 20

#max number of cells to read in memory
max_mem <- args[9] #param 21
in_dir_list_filename <- args[10] #param 22
#rasterOptions(maxmemory=1e+07,timer=TRUE)


list_param_run_assessment_part2_plotting <-list(  
    in_dir,y_var_name, interpolation_method, out_suffix,
    out_dir, create_out_dir_param, mosaic_plot, proj_str, file_format, NA_flag_val,
    multiple_region, countries_shp, plot_region, num_cores,
    region_name, df_assessment_files_name, threshold_missing_day,year_predicted
  )

names(list_param_run_assessment_part2_plotting) <- c(
  "in_dir","y_var_name","interpolation_method","out_suffix",
    "out_dir","create_out_dir_param","mosaic_plot","proj_str","file_format","NA_flag_val",
    "multiple_region","countries_shp","plot_region","num_cores",
    "region_name","df_assessment_files_name","threshold_missing_day","year_predicted"
  )

list_param_run_assessment_combined_region_plotting_prediction <-list(  
    in_dir_list_filename,
    in_dir,y_var_name, interpolation_method, out_suffix,
    out_dir, create_out_dir_param, mosaic_plot, proj_str, file_format, NA_flag_val,
    multiple_region, countries_shp, plot_region, num_cores,
    region_name, df_assessment_files_name, threshold_missing_day,year_predicted)

names(list_param_run_assessment_combined_region_plotting_prediction) <- c(
    "in_dir_list_filename",
    "in_dir","y_var_name","interpolation_method","out_suffix",
    "out_dir","create_out_dir_param","mosaic_plot","proj_str","file_format","NA_flag_val",
    "multiple_region","countries_shp","plot_region","num_cores",
    "region_name","df_assessment_files_name","threshold_missing_day","year_predicted")

i <- 1 #this select the first year of list_year_predicted
#Step 1: run figures production by region using table (part2 assessment script)
#Step 2: run figures and tables generation across region and years
#Step 3: latex/slidify presentation?
#Step 4: latex/slidify presentation?

if (stages_to_run[8]==8){
  
  #Step 1: run individual figure production if needed:
  if(run_figure_by_year==TRUE){
    #debug(run_assessment_plotting_prediction_fun)
    df_assessment_figures_files <-
    run_assessment_plotting_prediction_fun(list_param_run_assessment_plotting)
  }

  #Step 2: run combination of all files...
  
  #function_assessment_part2 <- "global_run_scalingup_assessment_part2_01032016.R"
  #source(file.path(script_path,function_assessment_part2)) #source all functions used in this script

  #debug(run_assessment_combined_region_plotting_prediction_fun)
  df_assessment_combined_figures_files <-
  run_assessment_combined_region_plotting_prediction_fun(list_param_run_assessment_combined_region_plotting_prediction)
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

#"Rscript %s %s wgs84Grid %s %s %s %s/subset/mean_LST_%s_jan_6_wgs84.tif FALSE %s/%s/covar_obj_%s.RData %s/%s/%s/met_stations_outfiles_obj_gam_CAI_%s.RData 10 4800  %s %s > %s/outLogs/%s_stage4_%s.log 2>  %s/outLogs/%s_stage4_err_%s.log" % (scriptFn,ll,ll,outputFolder,b[0],outputFolder,ll,outputFolder,ll,ll,outputFolder,ll,year,ll,year,yearInt,outputFolder,ll,year,outputFolder,ll,year)
      #print outSt

