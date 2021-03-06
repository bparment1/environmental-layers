##################    Master script for temperature predictions  #######################################
############################ TMIN AND TMAX predictions ##########################################
#                           
##This script produces intperpolated surface of TMIN and TMAX for specified processing region(s) given sets 
#of inputs and parameters.
#STAGE 1: LST climatology downloading and/or calculation
#STAGE 2: Covariates preparation for study/processing area: calculation of covariates (spect,land cover,etc.) and reprojection
#STAGE 3: Data preparation: meteorological station database query and extraction of covariates values from raster brick
#STAGE 4: Raster prediction: run interpolation method (-- gam fusion, gam CAI, ...) and perform validation 
#STAGE 5: Output analyses: assessment of results for specific dates...
#
#AUTHOR: Benoit Parmentier                                                                       
#DATE: 11/29/2013                                                                                 

#PROJECT: NCEAS INPLANT: Environment and Organisms --TASK#363, TASK$568--   

## TODO:
# Modify code for stage 1 and call python script from R in parallel
# Add options to run only specific stage + additional out_suffix?
# Make master script a function?
# Add log file for master script,add function to collect inputs and outputs
##################################################################################################

###Loading R library and packages   
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

#/data/project/layers/commons/data_workflow : this directory contains the input data and scripts

script_path<-"/data/project/layers/commons/data_workflow/env_layers_scripts/"

##SCRIPT USED FOR THE PREDICTIONS: Source or list all scripts here to avoid confusion on versions being run!!!!

#source(file.path(script_path,"master_script_temp_08122013.R")) #Master script can be run directly...

#CALLED FROM MASTER SCRIPT:

modis_download_script <- file.path(script_path,"modis_download_05142013.py") # LST modis download python script
clim_script <- file.path(script_path,"climatology_05312013.py") # LST climatology python script
grass_setting_script <- file.path(script_path,"grass-setup.R") #Set up system shell environment for python+GRASS
#source(file.path(script_path,"download_and_produce_MODIS_LST_climatology_06112013.R"))
source(file.path(script_path,"covariates_production_temperatures_08052013.R"))
source(file.path(script_path,"Database_stations_covariates_processing_function_06112013.R"))
source(file.path(script_path,"GAM_fusion_analysis_raster_prediction_multisampling_11032013.R"))
source(file.path(script_path,"results_interpolation_date_output_analyses_08052013.R"))
#source(file.path(script_path,"results_covariates_database_stations_output_analyses_04012013.R")) #to be completed

#FUNCTIONS CALLED FROM GAM ANALYSIS RASTER PREDICTION ARE FOUND IN...

source(file.path(script_path,"sampling_script_functions_08252013.R"))
source(file.path(script_path,"GAM_fusion_function_multisampling_11032013.R")) #Includes Fusion and CAI methods
source(file.path(script_path,"interpolation_method_day_function_multisampling_10112013.R")) #Include GAM_day
source(file.path(script_path,"GAM_fusion_function_multisampling_validation_metrics_10102013.R"))

#stages_to_run<-c(1,2,3,4,5) #May decide on antoher strategy later on...
#stages_to_run<-c(0,2,3,4,5) #May decide on antoher strategy later on...
stages_to_run<-c(0,2,3,4,5) #MRun only raster fitting, prediction and assessemnt (providing lst averages, covar brick and met stations)
#If stage 2 is skipped then use previous covar object
covar_obj_file<-"/data/project/layers/commons/data_workflow/output_data_365d_gam_fus_lst_test_run_07172013/covar_obj__365d_gam_fus_lst_test_run_07172013.RData"
#covar_obj_file<-"covar_obj__365d_gam_CAI_lst_comb3_08252013.RData"
#If stage 3 is skipped then use previous met_stations object
met_stations_outfiles_obj_file<-"/data/project/layers/commons/data_workflow/output_data_365d_gam_fus_lst_test_run_07172013/met_stations_outfiles_obj_gam_fusion__365d_gam_fus_lst_test_run_07172013.RData"
#met_stations_outfiles_obj_file<-"met_stations_outfiles_obj_gam_CAI__365d_gam_CAI_lst_comb3_08252013.RData"

var<-"TMAX" # variable being interpolated
#out_prefix<-"_365d_gam_cai_lst_comb3_10102013"                #User defined output prefix
out_prefix<-"_365d_gwr_daily_lst_comb5p4_7_11292013"                #User defined output prefix

out_suffix<-"_OR_11292013"                                       #Regional suffix
out_suffix_modis <-"_05302013"                       #pattern to find tiles produced previously     

#interpolation_method<-c("gam_fusion","gam_CAI","gam_daily") #other otpions to be added later
#interpolation_method<-c("gam_CAI") #other otpions to be added later
#interpolation_method<-c("gam_fusion") #other otpions to be added later
#interpolation_method<-c("kriging_fusion") #other otpions to be added later
#interpolation_method<-c("gwr_fusion") #other otpions to be added later
#interpolation_method<-c("gwr_CAI") #other otpions to be added later
#interpolation_method<-c("kriging_CAI") 

#interpolation_method<-c("gam_daily") #other otpions to be added later
#interpolation_method<-c("kriging_daily") #other otpions to be added later
interpolation_method<-c("gwr_daily") #other otpions to be added later

#out_path<-"/home/parmentier/Data/IPLANT_project/Oregon_interpolation/Oregon_03142013/output_data"
out_path <- "/data/project/layers/commons/Oregon_interpolation/output_data"

#out_path<-"/data/project/layers/commons/data_workflow/output_data"
out_path <-paste(out_path,out_prefix,sep="")

if (!file.exists(out_path)){
  dir.create(out_path)
  #} else{
  #  out_path <-paste(out_path..)
}
  
lc_path<-"/data/project/layers/commons/data_workflow/inputs/lc-consensus-global"
infile_modis_grid<-"/data/project/layers/commons/data_workflow/inputs/modis_grid/modis_sinusoidal_grid_world.shp" #modis grid tiling system, global
infile_elev<-"/data/project/layers/commons/data_workflow/inputs/dem-cgiar-srtm-1km-tif/srtm_1km.tif"  #elevation at 1km, global extent to be replaced by the new fused product 
#infile_canheight<-"/data/project/layers/commons/data_workflow/inputs/treeheight-simard2011/Simard_Pinto_3DGlobalVeg_JGR.tif"         #Canopy height, global extent
infile_canheight<-"/data/project/layers/commons/data_workflow/inputs/treeheight-simard2011/canheight_int.tif" #changed to INT4S
infile_distoc <- "/data/project/layers/commons/data_workflow/inputs/distance_to_coast/GMT_intermediate_coast_distance_01d_rev.tif" #distance to coast, global extent at 0.01 deg
#infile_covariates<- "/home/parmentier/Data/IPLANT_project/Venezuela_interpolation/Venezuela_01142013/covariates_Oregon_region_TMAX__OR_04052013.tif" #Oregon covar TMAX from earlier codes...for continuity
#infile_reg_outline=""  #input region outline defined by polygon: none for Venezuela
#This is the shape file of outline of the study area                                                      #It is an input/output of the covariate script
infile_reg_outline <- "/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/OR83M_state_outline.shp"  #input region outline defined by polygon: Oregon
#ref_rast_name<-""  #local raster name defining resolution, exent, local projection--. set on the fly?? 
#this may be redundant with infile_reg_outline
ref_rast_name<-"/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/mean_day244_rescaled.rst"  #local raster name defining resolution, exent: oregon
buffer_dist<-0 #not in use yet, must change climatology step to make sure additional tiles are downloaded and LST averages
               #must also be calculated for neighbouring tiles.

#list_tiles_modis <- c("h11v08,h11v07,h12v07,h12v08,h10v07,h10v08") #tile for Venezuela and surrounding area
list_tiles_modis <- c("h08v04,h09v04") #tiles for Oregon
  
#CRS_interp<-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs";
CRS_interp <-"+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs";

#"+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80"
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
#out_region_name<-"_venezuela_region" #generated on the fly
out_region_name<-"_oregon_region" #generated on the fly
  
#The names of covariates can be changed...these names should be output/input from covar script!!!
rnames<-c("x","y","lon","lat","N","E","N_w","E_w","elev_s","slope","aspect","CANHGHT","DISTOC")
lc_names<-c("LC1","LC2","LC3","LC4","LC5","LC6","LC7","LC8","LC9","LC10","LC11","LC12")
lst_names<-c("mm_01","mm_02","mm_03","mm_04","mm_05","mm_06","mm_07","mm_08","mm_09","mm_10","mm_11","mm_12",
               "nobs_01","nobs_02","nobs_03","nobs_04","nobs_05","nobs_06","nobs_07","nobs_08",
               "nobs_09","nobs_10","nobs_11","nobs_12")
covar_names<-c(rnames,lc_names,lst_names)
  
list_val_range <-c("lon,-180,180","lat,-90,90","N,-1,1","E,-1,1","N_w,-1,1","E_w,-1,1","elev_s,0,6000","slope,0,90",
                   "aspect,0,360","DISTOC,-0,10000000","CANHGHT,0,255","LC1,0,100","LC5,0,100","mm_01,-15,50",
                   "mm_02,-15,50","mm_03,-15,50","mm_04,-15,50","mm_05,-15,50","mm_06,-15,50","mm_07,-15,50",
                   "mm_08,-15,50","mm_09,-15,50","mm_10,-15,50","mm_11,-15,50","mm_12,-15,50")

############ STAGE 1: LST Climatology ###############

#Parameters,Inputs from R to Python??
start_year = "2001"
end_year = "2010"
hdfdir <- "/data/project/layers/commons/data_workflow/inputs/MOD11A1_tiles" #path directory to MODIS data
download=0 #download MODIS product if 1
clim_calc=1 #calculate lst averages/climatology if 1

list_param_download_clim_LST_script <- list(list_tiles_modis,start_year,end_year,hdfdir,
                                            var,grass_setting_script,modis_download_script, clim_script,
                                            download,clim_calc,out_suffix_modis)
names(list_param_download_clim_LST_script)<-c("list_tiles_modis","start_year","end_year","hdfdir",
                                              "var","grass_setting_script","modis_download_script","clim_script",
                                              "download","clim_calc","out_suffix_modis")
no_tiles <- length(unlist(strsplit(list_tiles_modis,",")))  # transform string into separate element in char vector

if (stages_to_run[1]==1){
  #clim_production_obj <-mclapply(1:2, list_param=list_param_download_clim_LST_script, download_calculate_MODIS_LST_climatology,mc.preschedule=FALSE,mc.cores = 2) #This is the end bracket from mclapply(...) statement
  clim_production_obj <-lapply(1:no_tiles, list_param=list_param_download_clim_LST_script, download_calculate_MODIS_LST_climatology) #,mc.preschedule=FALSE,mc.cores = 2) #This is the end bracket from mclapply(...) statement
  
}
#Collect LST climatology list as output???

############ STAGE 2: Covariate production ################

#list of 18 parameters
list_param_covar_production<-list(var,out_path,lc_path,infile_modis_grid,infile_elev,infile_canheight,
                                  infile_distoc,list_tiles_modis,infile_reg_outline,CRS_interp,CRS_locs_WGS84,out_region_name,
                                  buffer_dist,list_val_range,out_suffix,out_suffix_modis,ref_rast_name,hdfdir,covar_names) 

names(list_param_covar_production)<-c("var","out_path","lc_path","infile_modis_grid","infile_elev","infile_canheight",
                                      "infile_distoc","list_tiles_modis","infile_reg_outline","CRS_interp","CRS_locs_WGS84","out_region_name",
                                      "buffer_dist","list_val_range","out_suffix","out_suffix_modis","ref_rast_name","hdfdir","covar_names") 

## Modify to store infile_covar_brick in output folder!!!
if (stages_to_run[2]==2){
  covar_obj <- covariates_production_temperature(list_param_covar_production)
  infile_covariates <- covar_obj$infile_covariates
  infile_reg_outline <- covar_obj$infile_reg_outline
  covar_names<- covar_obj$covar_names
}else{
  covar_obj <-load_obj(covar_obj_file)
  infile_covariates <- covar_obj$infile_covariates
  infile_reg_outline <- covar_obj$infile_reg_outline
  covar_names<- covar_obj$covar_names
  
}

#Note that if stages_to_run[2]!=2, then use values defined at the beginning of the script for infile_covariates and infile_reg_outline

############# STAGE 3: Data preparation ###############

#specific to this stage
db.name <- "ghcn"       # name of the Postgres database
range_years<-c("2010","2011") #right bound not included in the range!!
range_years_clim<-c("2000","2011") #right bound not included in the range!!
infile_ghncd_data <-"/data/project/layers/commons/data_workflow/inputs/ghcn/v2.92-upd-2012052822/ghcnd-stations.txt"                              #This is the textfile of station locations from GHCND
#qc_flags_stations<-c("0","S")    #flags allowed for screening after the query from the GHCND??
qc_flags_stations<-c("0")    #flags allowed for screening after the query from the GHCND??

#infile_covariates and infile_reg_outline defined in stage 2 or at the start of script...

#list of 12 parameters for input in the function...

list_param_prep<-list(db.name,var,range_years,range_years_clim,infile_reg_outline,infile_ghncd_data,infile_covariates,CRS_locs_WGS84,out_path,covar_names,qc_flags_stations,out_prefix)
cnames<-c("db.name","var","range_years","range_years_clim","infile_reg_outline","infile_ghncd_data","infile_covariates","CRS_locs_WGS84","out_path","covar_names","qc_flags_stations","out_prefix")
names(list_param_prep)<-cnames

##### RUN SCRIPT TO GET STATION DATA WITH COVARIATES #####

if (stages_to_run[3]==3){
  list_outfiles<-database_covariates_preparation(list_param_prep)
}else{
  list_outfiles <-load_obj(met_stations_outfiles_obj_file) 
}

############### STAGE 4: RASTER PREDICTION #################

#Prepare parameters for for raster prediction... 

#Collect parameters from the previous stage: data preparation stage

#3 parameters from output
infile_monthly<-list_outfiles$monthly_covar_ghcn_data #outile4 from database_covar script
infile_daily<-list_outfiles$daily_covar_ghcn_data  #outfile3 from database_covar script
infile_locs<- list_outfiles$loc_stations_ghcn #outfile2? from database covar script

#names(outfiles_obj)<- c("loc_stations","loc_stations_ghcn","daily_covar_ghcn_data","monthly_covar_ghcn_data")

list_param_data_prep <- list(infile_monthly,infile_daily,infile_locs,infile_covariates,covar_names,var,out_prefix,CRS_locs_WGS84)
names(list_param_data_prep) <- c("infile_monthly","infile_daily","infile_locs","infile_covariates","covar_names","var","out_prefix","CRS_locs_WGS84")

#Set additional parameters
#Input for sampling function...need to reorganize inputs!!!
seed_number<- 100  #if seed zero then no seed?     

nb_sample<-1           #number of time random sampling must be repeated for every hold out proportion
#step<- 0.1         
step<- 0         
constant<-0             #if value 1 then use the same samples as date one for the all set of dates
prop_minmax<-c(0.3,0.3)  #if prop_min=prop_max and step=0 then predictions are done for the number of dates...
#prop_minmax<-c(0.1,0.7)  #if prop_min=prop_max and step=0 then predictions are done for the number of dates...

seed_number_month <- 100
nb_sample_month <-1           #number of time random sampling must be repeated for every hold out proportion
step_month <-0         
#step_month <-0.1
constant_month <- 0             #if value 1 then use the same samples as date one for the all set of dates
prop_minmax_month <-c(0,0)  #if prop_min=prop_max and step=0 then predictions are done for the number of dates...

#dates_selected<-c("20100101","20100102","20100103","20100901") # Note that the dates set must have a specific format: yyymmdd
#dates_selected<-c("20100101","20100102","20100301","20100302","20100501","20100502","20100701","20100702","20100901","20100902","20101101","20101102")
dates_selected<-"" # if empty string then predict for the full year specified earlier
#dates_selected <- 2 # if integer then predict for the evert n dat in the year specified earlier

screen_data_training<- FALSE #screen training data for NA and use same input training for all models fitted
use_clim_image <- TRUE # use predicted image as a base...rather than average Tmin at the station for delta
join_daily <- FALSE # join monthly and daily station before calucating delta

#Models to run...this can be changed for each run
#LC1: Evergreen/deciduous needleleaf trees

#Combination 5: for paper multi-timescale  paper
#list_models<-c("y_var ~ s(lat,lon)",
#               "y_var ~ s(lat,lon) + s(LST)",
#               "y_var ~ s(lat,lon) + s(elev_s)",
#               "y_var ~ s(lat,lon) + s(elev_s) + s(N_w,E_w)",
#               "y_var ~ s(lat,lon) + s(elev_s) + s(DISTOC)",
#               "y_var ~ s(lat,lon) + s(elev_s) + s(LST)",
#               "y_var ~ s(lat,lon) + s(elev_s) + s(LST) + ti(LST,LC1)")

#Combination 5: for paper multi-timescale  paper
#list_models<-c("y_var ~ lat*lon",
#               "y_var ~ lat*lon + LST",
#               "y_var ~ lat*lon + elev_s")
list_models<-c("y_var ~ lat*lon + elev_s + N_w*E_w",
               "y_var ~ lat*lon + elev_s + DISTOC",
               "y_var ~ lat*lon + elev_s + LST",
               "y_var ~ lat*lon + elev_s + LST + I(LST*LC1)")

#Combination 3: for paper baseline=s(lat,lon)+s(elev)
#list_models<-c("y_var ~ s(lat,lon) + s(elev_s)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(N_w)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(E_w)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(LST)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(DISTOC)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(LC1)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(CANHGHT)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(LST) + ti(LST,LC1)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(LST) + ti(LST,CANHGHT)")

#Combination 4: for paper baseline=s(lat,lon)
#list_models<-c("y_var ~ s(lat,lon)",
#              "y_var ~ s(lat,lon) + s(elev_s)",
#                "y_var ~ s(lat,lon) + s(N_w)",
#                "y_var ~ s(lat,lon) + s(E_w)",
#               "y_var ~ s(lat,lon) + s(LST)",
#               "y_var ~ s(lat,lon) + s(DISTOC)",
#              "y_var ~ s(lat,lon) + s(LC1)",
#             "y_var ~ s(lat,lon) + s(CANHGHT)",
#            "y_var ~ s(lat,lon) + s(LST) + ti(LST,LC1)",
#           "y_var ~ s(lat,lon) + s(LST) + ti(LST,CANHGHT)")

#Combination 5: for paper baseline=s(lat,lon)+s(elev)
#list_models<-c("y_var ~ s(lat,lon) + s(elev_s)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(N_w,E_w)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(LST)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(DISTOC)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(LC1)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(CANHGHT)",
#                "y_var ~ s(lat,lon) + s(elev_s) + s(LST) + ti(LST,LC1)")#,
#                #"y_var ~ s(lat,lon) + s(elev_s) + s(LST) + ti(LST,CANHGHT)")

#list_models2<-c("y_var ~ s(lat,lon) + s(DISTOC)")
list_models2 <- NULL
interp_method2 <- NULL #other options are "gwr" and "kriging"

#interp_method2 <- "gam" #other options are "gwr" and "kriging"
#list_models <- c("y_var ~ lat*lon + elev_s")

#list_models<-c("y_var ~ s(lat,lon) + s(elev_s)")

#list_models<-c("y_var ~ lat*lon + elev_s") #,
#               "y_var ~ lat*lon + elev_s + N_w",
#               "y_var ~ lat*lon + elev_s + E_w",
#               "y_var ~ lat*lon + elev_s + LST",
#               "y_var ~ lat*lon + elev_s + DISTOC",
#               "y_var ~ lat*lon + elev_s + LC1",
#               "y_var ~ lat*lon + elev_s + CANHGHT",
#               "y_var ~ lat*lon + elev_s + LST + I(LST*LC1)",
#               "y_var ~ lat*lon + elev_s + LST + I(LST*CANHGHT)")

#Default name of LST avg to be matched               
lst_avg<-c("mm_01","mm_02","mm_03","mm_04","mm_05","mm_06","mm_07","mm_08","mm_09","mm_10","mm_11","mm_12")  

#Collect all parameters in a list
list_param_raster_prediction<-list(list_param_data_prep,screen_data_training,
                                seed_number,nb_sample,step,constant,prop_minmax,dates_selected,
                                seed_number_month,nb_sample_month,step_month,constant_month,prop_minmax_month,
                                list_models,list_models2,interp_method2,lst_avg,out_path,script_path,use_clim_image,join_daily,
                                interpolation_method)
names(list_param_raster_prediction)<-c("list_param_data_prep","screen_data_training",
                                "seed_number","nb_sample","step","constant","prop_minmax","dates_selected",
                                "seed_number_month","nb_sample_month","step_month","constant_month","prop_minmax_month",
                                "list_models","list_models2","interp_method2","lst_avg","out_path","script_path","use_clim_image","join_daily",
                                "interpolation_method")

#debug(raster_prediction_fun)
#debug(debug_fun_test)
#debug_fun_test(list_param_raster_prediction)

raster_prediction_obj <- raster_prediction_fun(list_param_raster_prediction)

############## STAGE 5: OUTPUT ANALYSES ##################

date_selected_results<-c("20100101") 
list_param_results_analyses<-list(out_path,script_path,raster_prediction_obj,interpolation_method,
                                  covar_obj,date_selected_results,var,out_prefix)
names(list_param_results_analyses)<-c("out_path","script_path","raster_prediction_obj","interpolation_method",
                     "covar_obj","date_selected_results","var","out_prefix")
#plots_assessment_by_date<-function(j,list_param){
if (stages_to_run[5]==5){
  #source(file.path(script_path,"results_interpolation_date_output_analyses_08052013.R"))
  #Use lapply or mclapply
  summary_v_day <-plots_assessment_by_date(1,list_param_results_analyses)
  #Call as function...
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

