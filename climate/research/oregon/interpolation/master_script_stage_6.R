##################    Master script for climate predictions  #######################################
############################ TMIN AND TMAX predictions ##########################################
#                           
##This script produces intperpolated surface of TMIN and TMAX for specified processing region(s) given sets 
#of inputs and parameters.
#STAGE 1: LST climatology downloading and/or calculation
#STAGE 2: Covariates preparation for study/processing area: calculation of covariates (spect,land cover,etc.) and reprojection
#STAGE 3: Data preparation: meteorological station database query and extraction of covariates values from raster brick
#STAGE 4: Raster prediction: run interpolation method (-- gam fusion, gam CAI, ...) and perform validation 
#STAGE 5: Output analyses: assessment of results for specific dates...
#STAGE 6: Assessement of predictions by tiles and regions with mosaicing of predictions and accuracy
#AUTHOR: Benoit Parmentier                                                                        
#CREATED ON: 12/29/2015  
#MODIFIED ON: 01/03/2015  
#PROJECT: NCEAS INPLANT: Environment and Organisms                                                                           

## TODO:
# Add  assessment part 2 (figures)
# Add mosaicing part 2
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

script_path <- "/nobackupp8/bparmen1/env_layers_scripts" #path to script
function_assessment_part1_functions <- "global_run_scalingup_assessment_part1_functions_02112015.R" #PARAM12
function_assessment_part1a <-"global_run_scalingup_assessment_part1a_01042016.R"
function_assessment_part2 <- "global_run_scalingup_assessment_part2_01042016.R"
function_assessment_part2_functions <- "global_run_scalingup_assessment_part2_functions_01032016.R"
source(file.path(script_path,function_assessment_part1_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part1a)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2_functions)) #source all functions used in this script 

### Parameters and arguments ###
  
var<-"TMAX" # variable being interpolated
if (var == "TMAX") {
  y_var_name <- "dailyTmax"
  y_var_month <- "TMax"
}
if (var == "TMIN") {
  y_var_name <- "dailyTmin"
  y_var_month <- "TMin"
}

#interpolation_method<-c("gam_fusion") #other otpions to be added later
interpolation_method<-c("gam_CAI")
CRS_interp<-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs";
#CRS_interp <-"+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs";
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

out_region_name<-""
list_models<-c("y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)")

#reg1 (North Am), reg2(Europe),reg3(Asia), reg4 (South Am), reg5 (Africa), reg6 (Australia-Asia)
#master directory containing the definition of tile size and tiles predicted
in_dir1 <- "/nobackupp6/aguzman4/climateLayers/out/"
#/nobackupp6/aguzman4/climateLayers/out_15x45/1982

#region_names <- c("reg23","reg4") #selected region names, #PARAM2
region_name <- c("reg4") #run assessment by region, this is a unique region only
#region_names <- c("reg1","reg2","reg3","reg4","reg5","reg6") #selected region names, #PARAM2
interpolation_method <- c("gam_CAI") #PARAM4
out_prefix <- "run_global_analyses_pred_12282015" #PARAM5
out_dir <- "/nobackupp8/bparmen1/" #PARAM6
#out_dir <-paste(out_dir,"_",out_prefix,sep="")
create_out_dir_param <- TRUE #PARAM7

#CRS_interp<-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs";
#CRS_interp <-"+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs";
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

#list_year_predicted <- 1984:2004
list_year_predicted <- c("2014")
#year_predicted <- list_year_predicted[1]

file_format <- ".tif" #format for mosaiced files #PARAM10
NA_flag_val <- -9999  #No data value, #PARAM11
num_cores <- 6 #number of cores used #PARAM13
plotting_figures <- TRUE #running part2 of assessment to generate figures...
  
##Additional parameters used in part 2, some these may be removed as code is simplified
mosaic_plot <- FALSE #PARAM14
day_to_mosaic <- c("19920102","19920103","19920103") #PARAM15
multiple_region <- TRUE #PARAM16
countries_shp <- "/nobackupp8/bparmen1/NEX_data/countries.shp" #PARAM17
#countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #Atlas
plot_region <- TRUE  #PARAM18
threshold_missing_day <- c(367,365,300,200)#PARAM19

list_param_run_assessment_prediction <- list(in_dir1,region_name,y_var_name,interpolation_method,out_prefix,
                                  out_dir,create_out_dir_param,CRS_locs_WGS84,list_year_predicted,
                                  file_format,NA_flag_val,num_cores,plotting_figures,
                                  mosaic_plot,day_to_mosaic,multiple_region,countries_shp,plot_region)
list_names <- c("in_dir1","region_name","y_var_name","interpolation_method","out_prefix",
                                  "out_dir","create_out_dir_param","CRS_locs_WGS84","list_year_predicted",
                                  "file_format","NA_flag_val","num_cores","plotting_figures",
                                  "mosaic_plot","day_to_mosaic","multiple_region","countries_shp","plot_region")
                                      

names(list_param_run_assessment_prediction)<-list_names
  
#max number of cells to read in memory
max_mem<-args[11]
#rasterOptions(maxmemory=1e+07,timer=TRUE)

#debug(run_assessment_prediction_fun)
#debug(debug_fun_test)
#debug_fun_test(list_param_raster_prediction)
i <- 1 #this select the first year of list_year_predicted
if (stages_to_run[6]==6){
  assessment_prediction_obj <- run_assessment_prediction_fun(i,list_param_run_assessment_prediction)
}

## Add stage 7 (mosaicing) here??
#i <- 1 #this select the first year of list_year_predicted
#if (stages_to_run[7]==7){
#  assessment_prediction_obj <- run_assessment_prediction_fun(i,list_param_run_assessment_prediction)
#}

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

