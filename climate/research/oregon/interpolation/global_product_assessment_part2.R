####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Assessment of product part 1: mosaic and accuracy ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This part 2 of the assessment focuses on graphics to explore the spatial patterns of raster times series as figures and movie
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/03/2016  
#MODIFIED ON: 01/14/2018            
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
#COMMIT: clean up and testing animatio for reg6

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
  
#script_path <- "/nobackupp8/bparmen1/env_layers_scripts" #path to script
script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts" #path to script

## NASA poster and paper related
#source(file.path(script_path,"NASA2016_conference_temperature_predictions_function_05032016b.R"))

#Mosaic related on NEX
script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts"
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
#function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_09192016b.R"
#source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 
function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_01142018.R"
source(file.path(script_path,function_product_assessment_part2_functions)) #source all functions used in this script 

###############################
####### Parameters, constants and arguments ###

CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #constant 1
var<-"TMIN" # variable being interpolated #param 1, arg 1
interpolation_method<-c("gam_CAI") #param 2
CRS_interp <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #param 3
#list_models<-c("y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)") #param 4
metric_name <- "var_pred" #use RMSE if accuracy
in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/tMinOut/reg6/assessment2"
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/tMinOut/reg6/mosaics/mosaic/output_reg6_1984"
in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/tMinOut/reg6/mosaics/mosaic/"

#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaics/mosaic"
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaics/mosaic"
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/mosaic" #note dropped the s in mosaics

#reg1 (North Am), reg2(Europe),reg3(Asia), reg4 (South Am), reg5 (Africa), reg6 (Australia-Asia)
#master directory containing the definition of tile size and tiles predicted
region_name <- c("reg6") #param 6, arg 3
out_suffix <- "global_assessment_mosaic_reg6_01142018"
create_out_dir_param <- TRUE #param 9, arg 6
#out_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/assessment"
out_dir <- "/data/project/layers/commons/NEX_data/climateLayers/tMinOut/reg6/assessment2/"
#run_figure_by_year <- TRUE # param 10, arg 7

file_format <- ".tif" #format for mosaiced files # param 11
NA_flag_val <- -32768  #No data value, # param 12
plotting_figures <- TRUE #running part2 of assessment to generate figures... # param 14
num_cores <- 11 #number of cores used # param 13, arg 8
#python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin" #PARAM 30
python_bin <- "/usr/bin" #PARAM 30
date_start <- "19840101" #PARAM 12 arg 12
date_end <- "20141231" #PARAM 13 arg 13
#date_end <- "19841231" #PARAM 13 arg 13

#infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg5.tif"
#infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg4.tif"
#infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg6.tif"
#run_figure_by_year <- TRUE # param 10, arg 7
#list_year_predicted <- "1984,2014"
scaling <- 0.01 #was scaled on 100 
#if scaling is null then perform no scaling!!
NA_flag_val_mosaic <- -32768
in_dir_list_filename <- NULL #if NULL, use the in_dir directory to search for info
countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #Atlas
lf_raster <- NULL #list of raster to consider
#item_no <- 13 #was this value in the previous version in 2016
item_no <- 12
stat_opt <- TRUE
frame_speed <- 50
animation_format <- ".mp4" #options are mp4,avi,gif at this stage 
  
##################### START SCRIPT #################

## Setting up of parameters and constant

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

out_dir <- in_dir
if (create_out_dir_param == TRUE) {
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

#setwd(out_dir)

####### PART 1: Read in data ########

## List prediction mosaiced by region 
#pattern_str <- ".*.tif"
pattern_str <-"r_m_use_edge_weights_weighted_mean_mask_gam_CAI_dailyTmin_.*.tif"
#> lf_raster[c(1,367,733)]
#[1] "/data/project/layers/commons/NEX_data/climateLayers/tMinOut/reg6/mosaics/mosaic/output_reg6_1984/r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmin_19840101_reg6_1984.tif"       
#[2] "/data/project/layers/commons/NEX_data/climateLayers/tMinOut/reg6/mosaics/mosaic/output_reg6_1984/r_m_use_edge_weights_weighted_mean_masked_gam_CAI_dailyTmin_19840101_reg6_1984.tif"
#[3] "/data/project/layers/commons/NEX_data/climateLayers/tMinOut/reg6/mosaics/mosaic/output_reg6_1984/r_m_use_edge_weights_weighted_mean_mask_gam_CAI_dailyTmin_19840101_reg6_1984.tif"  

lf_raster <- list.files(path=in_dir_mosaic,
                        pattern=pattern_str,
                        recursive=T,
                        full.names=T)
r_stack <- stack(lf_raster,quick=T) #this is very fast now with the quick option!
#save(r_mosaic,file="r_mosaic.RData")

###############  PART2: Checking for mosaic produced for given region ##############
## From list of mosaic files predicted extract dates
## Check dates predicted against date range for a given date range

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
  
#Add report by year in text file?
#Using specified values for parameters
#debug(check_missing)
#This takes about 6 minutes for full time series (11204)
test_missing <- check_missing(lf=lf_raster, 
                              pattern_str=NULL,
                              in_dir=in_dir_mosaic,
                              date_start=date_start,
                              date_end=date_end,
                              item_no=item_no,
                              out_suffix=out_suffix,
                              num_cores=num_cores,
                              out_dir=out_dir)

df_time_series <- test_missing$df_time_series
head(df_time_series)

table(df_time_series$missing) 
count_by_year<- t(table(df_time_series$year)) #missing by year
count_by_year
#names(count_by_year)

#############################
##### PART 3: Creating animation based on prediction

#####
NAvalue(r_stack)
#plot(r_stack,y=6,zlim=c(-10000,10000)) #this is not rescaled
#plot(r_stack,zlim=c(-50,50),col=matlab.like(255))
var_name <- y_var_name

#df_raster <- read.table("df_raster_global_assessment_reg6_10102016.txt",sep=",",header=T)
#range_year <- c(1984,2014)
range_year <- range(as.numeric(as.character(df_time_series$year)),na.rm=T)
#subset_df_time_series <- subset(df_time_series,year%in% range_year)
#subset_df_time_series <- subset_df_time_series[!is.na(subset_df_time_series$lf),]

#lf_subset <- file.path(subset_df_time_series$dir,subset_df_time_series$lf)
range_year_str <- paste(range_year, sep = "_", collapse = "_")
out_suffix_str <- paste(range_year_str,out_suffix,sep="_")

#debug(plot_and_animate_raster_time_series)
#animation_frame_60_min_max_1984_1985_global_assessment_reg6_12112017.gif
#function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_01142018.R"
#source(file.path(script_path,function_product_assessment_part2_functions)) #source all functions used in this script 
## 12/13 at 21:01 to 23.42 for 613
#12/15:9:56
animation_obj <- plot_and_animate_raster_time_series(lf_raster, 
                                                     item_no=item_no,
                                                     region_name,
                                                     var_name,
                                                     metric_name,
                                                     NA_flag_val,
                                                     filenames_figures=NULL,
                                                     frame_speed=frame_speed,
                                                     animation_format=animation_format,
                                                     zlim_val=NULL,
                                                     plot_figure=T,
                                                     generate_animation=T,
                                                     stat_opt=stat_opt,
                                                     num_cores=num_cores,
                                                     out_suffix=out_suffix_str,
                                                     out_dir=out_dir)
  
stat_df_fname <- animation_obj$filenames_figures_mosaic
stat_df_fname <- list.files(pattern="stat_df.*.txt",full.names = F)

#stat_df <- read.table(file.path(out_dir,"stat_df_dailyTmin_var_pred_1984_2014_global_assessment_reg6_12132017.txt"),sep=",",header=T)
stat_df <- read.table(file.path(out_dir,stat_df_fname),sep=",",header=T)


View(stat_df)
plot(stat_df$mean,type="l",main="Mean tmin")
plot(stat_df$min, type="l",col="blue",main="Min tmin")
plot(stat_df$max,type="l",col="red",main="Max tmin")

plot(stat_df$mean,type="l",col="black")
lines(stat_df$min,type="l",col="blue")
lines(stat_df$max,type="l",col="red")

zlim_val <- c(-2000,5000)
animation_obj <- plot_and_animate_raster_time_series(basename(lf_subset), 
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
                                                     stat_opt=stat_opt,
                                                     num_cores=num_cores,
                                                     out_suffix=out_suffix_str,
                                                     out_dir=out_dir)

names(animation_obj)

stat_df_fname <- "stat_df_dailyTmin_var_pred_1984_2014_global_assessment_mosaic_reg6_12152017.txt"
stat_df <- read.table(
  stat_df_fname,
           sep=",",
  header=T)
############################ END OF SCRIPT ##################################

#ffmpeg -i yeay.gif outyeay.mp4

#/Applications/ffmpeg -r 25 -i input%3d.png -vcodec libx264 -x264opts keyint=25 -pix_fmt yuv420p -r 25 ../output.mp4

#ffmpeg -f gif -i file.gif -c:v libx264 outfile.mp4

#ffmpeg -i animation_frame_60_-2500_6000_.gif animation_frame_60_-2500_6000_.mp4
#http://user.astro.columbia.edu/~robyn/ffmpeghowto.html
#ffmpeg -f image2 -r 10 -i ./img%d.gif -b 600k ./out.mp4
#/data/project/layers/commons/NEX_data/climateLayers/tMinOut/reg6/assessment2/output_global_assessment_reg6_12112017/animation_frame_60_min_max_1984_1985_global_assessment_reg6_12112017.gif

#r: how many frames per seconds

#ffmpeg -i animation_frame_60_-2500_6000_.gif animation_frame_60_-2500_6000_.mp4

#ffmpeg -f gif -i animation_frame_60_-2500_6000_.gif -vcodec libx264 -x264opts keyint=25 -pix_fmt yuv420p -r 25 outfile.mp4
#ffmpeg -f gif -i animation_frame_60_-2500_6000_.gif -vcodec libx264 -x264opts keyint=11 -pix_fmt yuv420p -r 11 outfile.mp4

#ffmpeg -r 10 -i animation_frame_60_-2500_6000_.gif animation.avi

#This shrinks by 25 Mb
#ffmpeg -f gif -i animation_frame_60_-2500_6000_.gif -vcodec libx264 -x264opts -pix_fmt yuv420p outfile.mp4

#rate of one per second:
#ffmpeg -f image2 -r 1 -pattern_type glob -i '*.png' out.mp4
#crf: used for compression count rate factor is between 18 to 24, the lowest is the highest quality
#ffmpeg -f image2 -r 1 -vcodec libx264 -crf 24 -pattern_type glob -i '*.png' out.mp4
