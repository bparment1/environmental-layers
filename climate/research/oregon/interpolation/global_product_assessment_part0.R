####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Assessment of product part 1: mosaic and accuracy ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This part 2 of the assessment focuses on graphics to explore the spatial patterns of raster times series as figures and movie
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/03/2016  
#MODIFIED ON: 10/24/2016            
#Version: 1
#PROJECT: Environmental Layers project     
#COMMENTS: Initial commit, script based on part NASA biodiversity conferenc 
#TODO:
#1) Add plot broken down by year and region 
#2) Modify code for overall assessment accross all regions and year
#3) Clean up

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rm u:aguzman4:rwx /nobackupp6/aguzman4/climateLayers/LST_tempSpline/
#COMMIT: generating animation for region 4 for multiple years sequences

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
library(mosaic)

###### Function used in the script #######
  
#script_path <- "/nobackupp8/bparmen1/env_layers_scripts" #path to script
script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts" #path to script

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
function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_09192016b.R"
source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 
function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_10222016.R"
source(file.path(script_path,function_product_assessment_part2_functions)) #source all functions used in this script 

###############################
####### Parameters, constants and arguments ###

CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #constant 1

var<-"TMAX" # variable being interpolated #param 1, arg 1

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

#interpolation_method<-c("gam_fusion") #other otpions to be added later
interpolation_method<-c("gam_CAI") #param 2
CRS_interp <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #param 3
#CRS_interp <-"+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs";

#out_region_name<-""
#list_models<-c("y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)") #param 4
metric_name <- "var_pred" #use RMSE if accuracy

#reg1 (North Am), reg2(Europe),reg3(Asia), reg4 (South Am), reg5 (Africa), reg6 (Australia-Asia)
#master directory containing the definition of tile size and tiles predicted
#in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/assessment"
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/mosaic"
in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/assessment"
in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/mosaics/mosaic" #predicted mosaic
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaics/mosaic"
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaics/mosaic"
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/mosaic" #note dropped the s in mosaics

region_name <- c("reg6") #param 6, arg 3
out_suffix <- "global_assessment_reg6_10232016"

create_out_dir_param <- TRUE #param 9, arg 6


out_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/assessment"

#run_figure_by_year <- TRUE # param 10, arg 7

file_format <- ".tif" #format for mosaiced files # param 11
NA_flag_val <- -32768  #No data value, # param 12

#num_cores <- 6 #number of cores used # param 13, arg 8
plotting_figures <- TRUE #running part2 of assessment to generate figures... # param 14
num_cores <- 11 #number of cores used # param 13, arg 8
#python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin" #PARAM 30
python_bin <- "/usr/bin" #PARAM 30

day_start <- "1984101" #PARAM 12 arg 12
day_end <- "20141231" #PARAM 13 arg 13
#date_start <- day_start
#date_end <- day_end

#infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg4.tif"
#infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg5.tif"
#infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg4.tif"
#infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg6.tif"

#run_figure_by_year <- TRUE # param 10, arg 7
#list_year_predicted <- "1984,2014"
scaling <- 0.01 #was scaled on 100 
#if scaling is null then perform no scaling!!

#df_centroids_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/output_reg5_1999/df_centroids_19990701_reg5_1999.txt"
#df_centroids_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/output_reg4_1999/df_centroids_19990701_reg4_1999.txt"
#df_centroids_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/mosaic/output_reg6_1984/df_centroids_19840101_reg6_1984.txt"
#/nobackupp6/aguzman4/climateLayers/out/reg1/assessment//output_reg1_1984/df_assessment_files_reg1_1984_reg1_1984.txt

#dates to plot and analyze

#l_dates <- c("19990101","19990102","19990103","19990701","19990702","19990703")
#l_dates <- c("19990101","19990102","19990103","19990104","19990105") 
#df_points_extracted_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/int_mosaics/data_points_extracted.txt"
#df_points_extracted_fname <- NULL #if null extract on the fly
#r_mosaic_fname <- "r_mosaic.RData"
#r_mosaic_fname <- NULL #if null create a stack from input dir

#NA_flag_val_mosaic <- -3399999901438340239948148078125514752.000
NA_flag_val_mosaic <- -32768
in_dir_list_filename <- NULL #if NULL, use the in_dir directory to search for info
countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #Atlas
lf_raster <- NULL #list of raster to consider
item_no <- 13

##################### START SCRIPT #################

####### PART 1: Read in data ########
out_dir <- in_dir
if (create_out_dir_param == TRUE) {
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

#setwd(out_dir)

###########  ####################

############ Using predicting first ##########

## using predictions
#pattern_str <- ".*.tif"
pattern_str <-"*.tif"
lf_raster <- list.files(path=in_dir_mosaic,pattern=pattern_str,recursive=F,full.names=T)
r_stack <- stack(lf_raster,quick=T) #this is very fast now with the quick option!
#save(r_mosaic,file="r_mosaic.RData")

#### check for missing dates from list of tif

###This should be a function!!

#####  This could be moved in a separate file!!
###############  PART4: Checking for mosaic produced for given region ##############
## From list of mosaic files predicted extract dates
## Check dates predicted against date range for a given date range
## Join file information to centroids of tiles data.frame
#list_dates_produced <- unlist(mclapply(1:length(lf_mosaic_list),FUN=extract_date,x=lf_mosaic_list,item_no=13,mc.preschedule=FALSE,mc.cores = num_cores))                         
#list_dates_produced <-  mclapply(1:2,FUN=extract_date,x=lf_mosaic_list,item_no=13,mc.preschedule=FALSE,mc.cores = 2)                         
item_no <- 13
date_start <- day_start
date_end <- day_end
#day_start <- "1984101" #PARAM 12 arg 12
#day_end <- "20141231" #PARAM 13 arg 13

#Using default values for parameters exectpt for num_cores=11 instead of 1
#debug(check_missing)
#test_missing <- check_missing(lf=lf_raster, 
#                              pattern_str=NULL,
#                              in_dir=".", #this is not used if lf is given
#                              date_start="1984101",
#                              date_end="20141231",
#                              item_no=13,
#                              out_suffix="",
#                              num_cores=num_cores,
#                              out_dir=".")
  
##Run this on reg4 and reg5 after
#Add report by year in text file?
#Using specified values for parameters
test_missing <- check_missing(lf=lf_raster, 
                              pattern_str=NULL,
                              in_dir=in_dir_mosaic,
                              date_start="1984101",
                              date_end="20141231",
                              item_no=13,
                              out_suffix=out_suffix,
                              num_cores=num_cores,
                              out_dir=".")

df_time_series <- test_missing$df_time_series
head(df_time_series)

table(df_time_series$missing)
table(df_time_series$year)

#############################
##### Creating animation based on prediction

#####
NAvalue(r_stack)
plot(r_stack,y=6,zlim=c(-10000,10000)) #this is not rescaled
#plot(r_stack,zlim=c(-50,50),col=matlab.like(255))
var_name <- "dailyTmax"

#debug(plot_and_animate_raster_time_series)

#metric_name <- "var_pred" #use RMSE if accuracy
#df_raster <- read.table("df_raster_global_assessment_reg6_10102016.txt",sep=",",header=T)
#plot_figure <- 
#function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_10222016.R"
#source(file.path(script_path,function_product_assessment_part2_functions)) #source all functions used in this script 

#undebug(plot_and_animate_raster_time_series)
range_year <- c(1984,1985)
subset_df_time_series <- subset(df_time_series,year%in% range_year)
subset_df_time_series <- subset_df_time_series[!is.na(subset_df_time_series$lf),]

lf_subset <- file.path(subset_df_time_series$dir,subset_df_time_series$lf)
range_year_str <- paste(range_year, sep = "_", collapse = "_")

out_suffix_str <- paste(range_year_str,out_suffix,sep="_")

#started on 10/22/2016 at 9.57
animation_obj <- plot_and_animate_raster_time_series(lf_subset, 
                                                     item_no,
                                                     region_name,
                                                     var_name,
                                                     metric_name,
                                                     NA_flag_val,
                                                     filenames_figures=NULL,
                                                     frame_speed=60,
                                                     animation_format=".gif",
                                                     zlim_val=NULL,
                                                     plot_figure=T,
                                                     generate_animation=T,
                                                     num_cores=num_cores,
                                                     out_suffix=out_suffix_str,
                                                     out_dir=out_dir)
  
zlim_val <- c(-2000,5000)
animation_obj <- plot_and_animate_raster_time_series(lf_subset, 
                                                     item_no,
                                                     region_name,
                                                     var_name,
                                                     metric_name,
                                                     NA_flag_val,
                                                     filenames_figures=NULL,
                                                     frame_speed=60,
                                                     animation_format=".gif",
                                                     zlim_val=zlim_val,
                                                     plot_figure=T,
                                                     generate_animation=T,
                                                     num_cores=num_cores,
                                                     out_suffix=out_suffix_str,
                                                     out_dir=out_dir)

#ffmpeg -i yeay.gif outyeay.mp4

#/Applications/ffmpeg -r 25 -i input%3d.png -vcodec libx264 -x264opts keyint=25 -pix_fmt yuv420p -r 25 ../output.mp4

#ffmpeg -f gif -i file.gif -c:v libx264 outfile.mp4

#ffmpeg -i animation_frame_60_-2500_6000_.gif animation_frame_60_-2500_6000_.mp4

#ffmpeg -i animation_frame_60_-2500_6000_.gif animation_frame_60_-2500_6000_.mp4

#ffmpeg -f gif -i animation_frame_60_-2500_6000_.gif -vcodec libx264 -x264opts keyint=25 -pix_fmt yuv420p -r 25 outfile.mp4
#ffmpeg -f gif -i animation_frame_60_-2500_6000_.gif -vcodec libx264 -x264opts keyint=11 -pix_fmt yuv420p -r 11 outfile.mp4

#ffmpeg -r 10 -i animation_frame_60_-2500_6000_.gif animation.avi

#ffmpeg -f gif -i animation_frame_60_-2500_6000_.gif -vcodec libx264 -x264opts -pix_fmt yuv420p outfile.mp4

#combine polygon
#http://gis.stackexchange.com/questions/155328/merging-multiple-spatialpolygondataframes-into-1-spdf-in-r

#http://gis.stackexchange.com/questions/116388/count-overlapping-polygons-in-single-shape-file-with-r

#### Use the predictions directory
#By region
#For each polygon/tile find polygon overlapping with count and ID (like list w)
#for each polygon/tile and date find if there is a prediction using the tif (multiply number)
#for each date of year report data in table.

#go through table and hsow if there are missing data (no prediction) or report min predictions for tile set?

############################################
#### Parameters and constants  

#on ATLAS
#in_dir1 <- "/data/project/layers/commons/NEX_data/test_run1_03232014/output" #On Atlas
#parent output dir : contains subset of the data produced on NEX
#in_dir1 <- "/data/project/layers/commons/NEX_data/output_run6_global_analyses_09162014/output20Deg2"
# parent output dir for the curent script analyes
#out_dir <-"/data/project/layers/commons/NEX_data/output_run3_global_analyses_06192014/" #On NCEAS Atlas
# input dir containing shapefiles defining tiles
#in_dir_shp <- "/data/project/layers/commons/NEX_data/output_run5_global_analyses_08252014/output/subset/shapefiles"

#On NEX
#contains all data from the run by Alberto
#in_dir1 <- " /nobackupp6/aguzman4/climateLayers/out_15x45/" #On NEX
#parent output dir for the current script analyes
#out_dir <- "/nobackup/bparmen1/" #on NEX
#in_dir_shp <- "/nobackupp4/aguzman4/climateLayers/output4/subset/shapefiles/"

#in_dir <- "/data/project/layers/commons/NEX_data/reg4_assessment"
#list_in_dir_run <-
#in_dir_list <-  c("output_run_global_analyses_pred_2009_reg4","output_run_global_analyses_pred_2010_reg4",
#                  "output_run_global_analyses_pred_2011_reg4","output_run_global_analyses_pred_2012_reg4",
#                  "output_run_global_analyses_pred_2013_reg4","output_run_global_analyses_pred_2014_reg4")
#in_dir_list_filename <- "/data/project/layers/commons/NEX_data/reg4_assessment/stage6_reg4_in_dir_list_02072016.txt"
#in_dir <- "" #PARAM 0
#y_var_name <- "dailyTmax" #PARAM1
#interpolation_method <- c("gam_CAI") #PARAM2
#out_suffix <- "global_analyses_overall_assessment_reg4_02072016"
#out_suffix <- "output_run10_1000x3000_global_analyses_02102015"
#out_suffix <- "run10_1500x4500_global_analyses_pred_1992_10052015" #PARAM3
#out_dir <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_10052015" #PARAM4
#create_out_dir_param <- TRUE #PARAM 5
#mosaic_plot <- FALSE #PARAM6
#if daily mosaics NULL then mosaicas all days of the year
#day_to_mosaic <- c("19920101","19920102","19920103") #PARAM7
#CRS_WGS84 <-    CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84 #CONSTANT1
#CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
#proj_str<- CRS_WGS84 #PARAM 8 #check this parameter
#file_format <- ".rst" #PARAM 9
#NA_value <- -9999 #PARAM10
#NA_flag_val <- NA_value
#multiple_region <- TRUE #PARAM 12
#region_name <- "world" #PARAM 13
#countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #PARAM 13, copy this on NEX too
#plot_region <- TRUE
#num_cores <- 6 #PARAM 14
#region_name <- c("reg4") #reference region to merge if necessary, if world all the regions are together #PARAM 16
#use previous files produced in step 1a and stored in a data.frame
#df_assessment_files <- "df_assessment_files_reg4_1984_run_global_analyses_pred_12282015.txt" #PARAM 17
#threshold_missing_day <- c(367,365,300,200) #PARM18

#list_param_run_assessment_plottingin_dir <- list(in_dir,y_var_name, interpolation_method, out_suffix, 
#                      out_dir, create_out_dir_param, mosaic_plot, proj_str, file_format, NA_value,
#                      multiple_region, countries_shp, plot_region, num_cores, 
#                      region_name, df_assessment_files, threshold_missing_day) 

#names(list_param_run_assessment_plottingin_dir) <- c("in_dir","y_var_name","interpolation_method","out_suffix", 
#                      "out_dir","create_out_dir_param","mosaic_plot","proj_str","file_format","NA_value",
#                      "multiple_region","countries_shp","plot_region","num_cores", 
#                      "region_name","df_assessment_files","threshold_missing_day") 

#run_assessment_plotting_prediction_fun(list_param_run_assessment_plottingin_dir) 

  ####### PARSE INPUT ARGUMENTS/PARAMETERS #####
  in_dir_list_filename <- list_param_run_assessment_plotting$in_dir_list_filename #PARAM 0
  in_dir <- list_param_run_assessment_plotting$in_dir #PARAM 1
  y_var_name <- list_param_run_assessment_plotting$y_var_name #PARAM2
  interpolation_method <- list_param_run_assessment_plotting$interpolation_method #c("gam_CAI") #PARAM3
  out_suffix <- list_param_run_assessment_plotting$out_suffix #PARAM4
  out_dir <- list_param_run_assessment_plotting$out_dir # PARAM5
  create_out_dir_param <- list_param_run_assessment_plotting$create_out_dir_param # FALSE #PARAM 6
  mosaic_plot <- list_param_run_assessment_plotting$mosaic_plot #FALSE #PARAM7
  proj_str<- list_param_run_assessment_plotting$proj_str #CRS_WGS84 #PARAM 8 #check this parameter
  file_format <- list_param_run_assessment_plotting$file_format #".rst" #PARAM 9
  NA_flag_val <- list_param_run_assessment_plotting$NA_flag_val #-9999 #PARAM10
  multiple_region <- list_param_run_assessment_plotting$multiple_region # <- TRUE #PARAM 11
  countries_shp <- list_param_run_assessment_plotting$countries_shp #<- "world" #PARAM 12
  plot_region <- list_param_run_assessment_plotting$plot_region # PARAM13 
  num_cores <- list_param_run_assessment_plotting$num_cores # 6 #PARAM 14
  region_name <- list_param_run_assessment_plotting$region_name #<- "world" #PARAM 15
  #df_assessment_files_name <- list_param_run_assessment_plotting$df_assessment_files_name #PARAM 16
  threshold_missing_day <- list_param_run_assessment_plotting$threshold_missing_day #PARM17
  year_predicted <- list_param_run_assessment_plotting$year_predicted
 
  NA_value <- NA_flag_val 
  metric_name <- "rmse" #to be added to the code later...
  
  ##################### START SCRIPT #################
  
  ####### PART 1: Read in data ########
  out_dir <- in_dir
  if(create_out_dir_param==TRUE){
    out_dir <- create_dir_fun(out_dir,out_suffix)
    setwd(out_dir)
  }else{
    setwd(out_dir) #use previoulsy defined directory
  }

  setwd(out_dir)
  
  list_outfiles <- vector("list", length=35) #collect names of output files, this should be dynamic?
  list_outfiles_names <- vector("list", length=35) #collect names of output files
  counter_fig <- 0 #index of figure to collect outputs
  
  #i <- year_predicted
  ###Table 1: Average accuracy metrics
  ###Table 2: daily accuracy metrics for all tiles

  if(!is.null(in_dir_list_filename)){
    in_dir_list <- as.list(read.table(in_dir_list_filename,stringsAsFactors=F)[,1])
  }else{
    pattern_str <- paste0("^output_",region_name,".*.")
    in_dir_list_all <- list.dirs(path=in_dir,recursive = T)
    in_dir_list <- in_dir_list_all[grep(pattern_str,basename(in_dir_list_all),invert=FALSE)] #select directory with shapefiles...
    #in_dir_shp <- file.path(in_dir_list_all,"shapefiles")
  }
  #pattern_str <- file.path(in_dir,paste0("output_",region_name,".*."))
  #test <- Sys.glob(pattern_str,FALSE)
  #  searchStr = paste(in_dir_tiles_tmp,"/*/",year_processed,"/gam_CAI_dailyTmax_predicted_",pred_mod_name,"*",day_to_mosaic[i],"*.tif",sep="")
  #  #print(searchStr)
  #  Sys.glob(searchStr)})

  #lf_mosaic <- lapply(1:length(day_to_mosaic),FUN=function(i){
  #  searchStr = paste(in_dir_tiles_tmp,"/*/",year_processed,"/gam_CAI_dailyTmax_predicted_",pred_mod_name,"*",day_to_mosaic[i],"*.tif",sep="")
  #  #print(searchStr)
  #  Sys.glob(searchStr)})

  ##Read in data list from in_dir_list
  #list_tb_fname <- list.files(path=file.path(in_dir,in_dir_list),"tb_diagnostic_v_NA_.*.txt",full.names=T)
  #list_df_fname <- list.files(path=file.path(in_dir,in_dir_list),"df_tile_processed_.*..txt",full.names=T)
  #list_summary_metrics_v_fname <- list.files(path=file.path(in_dir,in_dir_list),"summary_metrics_v2_NA_.*.txt",full.names=T)
  #list_tb_s_fname <- list.files(path=file.path(in_dir,in_dir_list),"tb_diagnostic_s_NA.*.txt",full.names=T)
  #list_tb_month_s_fname <- list.files(path=file.path(in_dir,in_dir_list),"tb_month_diagnostic_s.*.txt",full.names=T)
  #list_data_month_s_fname <- list.files(path=file.path(in_dir,in_dir_list),"data_month_s.*.txt",full.names=T)
  #list_data_s_fname <- list.files(path=file.path(in_dir,in_dir_list),"data_day_s.*.txt",full.names=T)
  #list_data_v_fname <- list.files(path=file.path(in_dir,in_dir_list),"data_day_v.*.txt",full.names=T)
  #list_pred_data_month_info_fname <- list.files(path=file.path(in_dir,in_dir_list),"pred_data_month_info.*.txt",full.names=T)
  #list_pred_data_day_info_fname <- list.files(path=file.path(in_dir,in_dir_list),"pred_data_day_info.*.txt",full.names=T)
  
  list_tb_fname <- list.files(path=in_dir_list,"tb_diagnostic_v_NA_.*.txt",full.names=T)
  list_df_fname <- list.files(path=in_dir_list,"df_tile_processed_.*..txt",full.names=T)
  list_summary_metrics_v_fname <- list.files(path=in_dir_list,"summary_metrics_v2_NA_.*.txt",full.names=T)
  list_tb_s_fname <- list.files(path=in_dir_list,"tb_diagnostic_s_NA.*.txt",full.names=T)
  list_tb_month_s_fname <- list.files(path=in_dir_list,"tb_month_diagnostic_s.*.txt",full.names=T)
  list_data_month_s_fname <- list.files(path=in_dir_list,"data_month_s.*.txt",full.names=T)
  list_data_s_fname <- list.files(path=in_dir_list,"data_day_s.*.txt",full.names=T)
  list_data_v_fname <- list.files(path=in_dir_list,"data_day_v.*.txt",full.names=T)
  list_pred_data_month_info_fname <- list.files(path=in_dir_list,"pred_data_month_info.*.txt",full.names=T)
  list_pred_data_day_info_fname <- list.files(path=in_dir_list,"pred_data_day_info.*.txt",full.names=T)
  
  #need to fix this !! has all of the files in one list (for a region)
  #list_shp <- list.files(path=file.path(in_dir,file.path(in_dir_list,"shapefiles")),"*.shp",full.names=T)

############################ END OF SCRIPT ##################################
