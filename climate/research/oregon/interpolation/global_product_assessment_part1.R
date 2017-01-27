#################################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Assessment of product part 1: extraction of values and matching to testing/traing ######################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This part 1 of the porduct assessment focuses on extraction from raster times series of mosaics.
#Observed data contained in testing and training data.frame by tiles are aggregated.
#Aggregated observed data are combiend with extracted predictions.
#There are three options for each step:
#step 1: combining testing and training from observed stations
#step 2: extraction from predicted raster mosaic
#step 3: combining extraction predicted and observed data
#
#AUTHOR: Benoit Parmentier 
#CREATED ON: 05/15/2016  
#MODIFIED ON: 01/27/2017            
#Version: 1
#PROJECT: Environmental Layers project     
#COMMENTS: clean up and moving function to function script
#TODO:
#1) Add plot broken down by year and region 
#2) Modify code for overall assessment accross all regions and year
#3) Clean up

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rmd user:aguzman4:rwx /nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015

#COMMIT: major clean up- visualization of extraction and mosacing from script and splitting code

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
function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_08232016.R" #Functions used to mosaic predicted tiles
function_mosaicing <-"global_run_scalingup_mosaicing_08222016.R" #main scripts for mosaicing predicted tiles

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
#global_product_assessment_part1_functions_01242017.R
function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_12182016b.R"
source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 
function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_01272017b.R"
source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 
##Don't load part 1 and part2, mosaic package does not work on NEX
#function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_09192016b.R"
#source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 
#function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_10222016.R"
#source(file.path(script_path,function_product_assessment_part2_functions)) #source all functions used in this script 
###############################
####### Parameters, constants and arguments ###

##Running steps
run_steps <- c(FALSE,FALSE,TRUE)
#step 1: combining testing and training from observed stations
#step 2: extraction from predicted raster mosaic
#step 3: combining extraction predicted and observed data

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

out_region_name<-""
list_models<-c("y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)") #param 4

#reg1 (North Am), reg2(Europe),reg3(Asia), reg4 (South Am), reg5 (Africa), reg6 (Australia-Asia)
#master directory containing the definition of tile size and tiles predicted
#in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/assessment"
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/mosaic"
in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/assessment"
in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaic/mosaic"

region_name <- c("reg1") #param 6, arg 3
out_suffix <- "_global_assessment_reg1_01112017"

create_out_dir_param <- TRUE #param 9, arg 6

out_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/assessment"

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
date_start <- day_start
date_end <- day_end

#infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg4.tif"
#infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg5.tif"
infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg1.tif"

#run_figure_by_year <- TRUE # param 10, arg 7
list_year_predicted <- "1984,2014"
scaling <- 0.01 #was scaled on 100 
#if scaling is null then perform no scaling!!

#df_centroids_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/output_reg5_1999/df_centroids_19990701_reg5_1999.txt"
df_centroids_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaic/output_reg1_1984/df_centroids_19840101_reg1_1984.txt"
#/nobackupp6/aguzman4/climateLayers/out/reg1/assessment//output_reg1_1984/df_assessment_files_reg1_1984_reg1_1984.txt

#dates to plot and analyze

#l_dates <- c("19990101","19990102","19990103","19990701","19990702","19990703")
l_dates <- c("19990101","19990102","19990103","19990104","19990105") 
#df_points_extracted_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/int_mosaics/data_points_extracted.txt"
df_points_extracted_fname <- NULL #if null extract on the fly
#r_mosaic_fname <- "r_mosaic.RData"
r_mosaic_fname <- NULL #if null create a stack from input dir

#NA_flag_val_mosaic <- -3399999901438340239948148078125514752.000
NA_flag_val_mosaic <- -32768
in_dir_list_filename <- NULL #if NULL, use the in_dir directory to search for info
countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #Atlas


##################### START SCRIPT #################

####### PART 1: Read in data ########
out_dir <- in_dir
if (create_out_dir_param == TRUE) {
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

setwd(out_dir)

###########  ####################

##Get the assessment information for every year for the matching region
if(!is.null(in_dir_list_filename)){
  in_dir_list <- as.list(read.table(in_dir_list_filename,stringsAsFactors=F)[,1])
}else{
  pattern_str <- paste0("^output_",region_name,".*.")
  in_dir_list_all <- list.dirs(path=in_dir,recursive = T)
  in_dir_list <- in_dir_list_all[grep(pattern_str,basename(in_dir_list_all),invert=FALSE)] #select directory with shapefiles...
  #in_dir_shp <- file.path(in_dir_list_all,"shapefiles")
}

## Now get the list of files for assessment of global product
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

### Get station data information: create a data.frame with stations info
data_date <- unlist(lapply(list_data_s_fname, function(x){unlist(strsplit(basename(x),"_"))[5]}))
#df_data_s_fname <- data.frame(data_date,region_name,dataset="data_s",file=list_data_s_fname)
year_str <- as.character(unlist(lapply(list_data_s_fname, function(x){unlist(strsplit(basename(x),"_"))[5]})))
df_data_s_fname <- data.frame(year_str,region_name,dataset="data_s",file=list_data_s_fname)

year_str <- as.character(unlist(lapply(list_data_v_fname, function(x){unlist(strsplit(basename(x),"_"))[5]})))
df_data_v_fname <- data.frame(year_str,region_name,dataset="data_v",file=list_data_v_fname)

df_data_v_fname$year <- df_data_v_fname$year_str
df_data_v_fname$year <- as.character(df_data_v_fname$year_str)
df_data_v_fname$file <- as.character(df_data_v_fname$file)
df_data_s_fname$file <- as.character(df_data_s_fname$file)
df_data_s_fname$year <- as.character(df_data_s_fname$year_str)

list_dates_val <- as.Date(strptime(l_dates,"%Y%m%d"))
l_dates_month_str <- format(list_dates_val, "%b") ## Month, char, abbreviated
l_dates_year_str <- format(list_dates_val, "%Y") ## Year with century
l_dates_day_str <- as.numeric(format(list_dates_val, "%d")) ## numeric month

##start new function
### Needs to make this a function!!!
#Step 1 for a list of dates, load relevant files with year matching,
#Step 2 for giving dates subset the data.frame
#Step 3: find duplicates, create return object and return given station data for the date

#### Before aggregation and extraction by stations combine data_s and data_v?
## retain all fields? this must be tested
list_data_v_fname
list_data_s_fname

###########################################################################
############### PART1: Combine testing and training stations #############
#######

if(run_steps[1]==TRUE){
  #selected_var <- c("elev","LC1","LC2","LC3","LC4","LC5","LC6","LC7",           
  #                    "LC8", "LC9", "LC10","LC11","LC12","nobs_01","nobs_02","nobs_03","nobs_04",     
  #                    "nobs_05","nobs_06","nobs_07","nobs_08","nobs_09","nobs_10","nobs_11",     
  #                    "nobs_12","lon","lat","N","E","N_w","E_w","elev_s","slope","aspect",     
  #                    "DISTOC","CANHGHT","mm_01","mm_02","mm_03","mm_04","mm_05",        
  #                    "mm_06","mm_07","mm_08","mm_09","mm_10","mm_11","mm_12")     #

  #selected_var <- c("mflag","qflag","sflag","elev","LC1","LC2","LC3","LC4","LC5","LC6","LC7",           
  #                  "LC8", "LC9", "LC10","LC11","LC12","nobs_01","nobs_02","nobs_03","nobs_04",     
  #                  "nobs_05","nobs_06","nobs_07","nobs_08","nobs_09","nobs_10","nobs_11",     
  #                  "nobs_12","lon","lat","N","E","N_w","E_w","elev_s","slope","aspect",     
  #                  "DISTOC","CANHGHT","mm_01","mm_02","mm_03","mm_04","mm_05",        
  #                  "mm_06","mm_07","mm_08","mm_09","mm_10","mm_11","mm_12")     

  #read from 
  #list_year_predicted <- "1984,2014"

  l_dates_year <- 1984:2014
  list_year_str <- unique(l_dates_year)
  list_out_suffix <- paste(region_name,list_year_str,sep="_")

  list_data_df_training <- list_data_s_fname
  list_data_df_testing <- list_data_v_fname

  #function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_01242017.R"
  #source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 
  #16:24
  #17:21
  #For 31 years and region 1, it took 2h27 minutes to produce the combine data
  #combine_and_aggregate_df_data_fun <- function(i,list_data_df_training,list_data_df_testing,selected_var=NULL,fun_selected_var="mean",list_out_suffix=NULL,out_dir="."){
  list_combine_agg_fname <- mclapply(1:length(list_year_str),
                                  FUN=combine_and_aggregate_df_data_fun,
                                  list_data_df_training=list_data_s_fname,
                                  list_data_df_testing=list_data_v_fname,
                                  selected_var=NULL,
                                  fun_selected_var="mean",
                                  list_out_suffix=list_out_suffix,
                                  out_dir=out_dir,
                                  mc.preschedule=FALSE,
                                  mc.cores = 11)
  #16:47
  #20:48


  #data_stations_var_pred1 <- list_combine_agg_fanme[[1]]$data_stations_var_pred
  #list_combine_agg_fanme[[1]]$data_stations_combined_v_s

  data_stations_var_pred1 <- read.table(list_combine_agg_fname[[31]]$data_stations_var_pred,header=T,stringsAsFactors = F,sep=",")
  #data_stations <- read.table(list_combine_agg_fname[[31]]$data_stations_combined_v_s,header=T,stringsAsFactors = F,sep=",")

  list_data_stations_var_pred_df_filename <- list.files(path=out_dir,pattern=paste0("data_stations_var_pred_",region_name,"_",".*.","txt"))

  ######## 
  ### Use combine training and testing data!!!

  #debug(aggregate_by_id_and_coord)
  data_name_point <- c("combined_data_v_s") #move this up, this can be any data.frame

  #for(k in 1:length(data_name_point)){
  
  out_suffix_str <- paste0(y_var_name,"_",interpolation_method,"_",region_name,out_suffix)
  list_out_suffix <- paste0(data_name_point[k],"_",1984:2014)
  #test<- mclapply(1:length(list_data_v_fname[1:1]),
  #                FUN=aggregate_by_id_and_coord,
  #         list_df_data = list_data_v_fname,
  #         list_out_suffix = list_out_suffix, 
  #         out_dir=out_dir,
  #         mc.preschedule=FALSE,
  #         mc.cores = 1
  #         )

  #21:34
  aggregated_data_points <- mclapply(1:length(list_data_stations_var_pred_df_filename),
                                     FUN=aggregate_by_id_and_coord,
                                     list_df_data = list_data_stations_var_pred_df_filename,
                                     list_out_suffix = list_out_suffix, 
                                     out_dir=out_dir,
                                     mc.preschedule=FALSE,
                                     mc.cores = 11)
  
  df_tmp <- do.call(rbind,aggregated_data_points)         
  df_data_points <- aggregate(id  ~ x + y,data=df_tmp,FUN=mean)
  filename_df_data_points <- file.path(out_dir,
                                       paste0("df_data_points_id_stations","_",region_name,"_",list_year_str[1],"_",list_year_str[length(list_year_str)],".txt"))
  write.table(df_data_points,filename_df_data_points)
}

##############################################################
######## #  Part 2: perform extraction ###########################################
## Do the extraction with training and test combined station data!!

if(run_steps[2]==TRUE){
#if(is.null(df_points_extracted_fname)){
  
  #df_data_points <- df_data_points 
  coords<- df_data_points[,c("x","y")]
  coordinates(df_data_points)<- coords #maybe change this to data_stations? #these are the lcoations of climate stations

  #now extract from mosaic
  #t<-unique(test$id)
  #md <- melt(pix_ts, id=(c("date")),measure.vars=c(var_pred_tmp, "missing")) #c("x","y","dailyTmax","mod1","res_mod1"))
  #formula_str <- "id + date ~ x + y + dailyTmax + mod1 + res_mod1"
  #pix_ts <- cast(md, date ~ variable, fun.aggregate = mean, 
  #function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_01142017.R"
  #source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script

  in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaics/mosaic"
  #in_dir_mosaic_rmse <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaicsRMSE/mosaic"
  pattern_str <- ".*.tif"
  lf_raster <- list.files(pattern=pattern_str,path = in_dir_mosaic)
  out_suffix_str <- paste0(region_name,"_",list_year_str[1],"_",list_year_str[length(list_year_str)])
  #started at 22:10pm
  #For region 1, with 4,803 stations for 31 years but missing some stations
  #undebug(extract_from_time_series_raster_stack)
  #extract from var pred mosaic, tmax in this case:
  extract_obj_var_pred <- extract_from_time_series_raster_stack(df_points=df_data_points,
                                                                date_start,
                                                                date_end,
                                                                lf_raster=NULL,
                                                                item_no=13,
                                                                num_cores=11,
                                                                pattern_str=NULL,
                                                                in_dir=in_dir_mosaic,
                                                                out_dir=out_dir,
                                                                out_suffix=out_suffix_str)
}else{
  extract_obj_fname <- file.path(out_dir,paste("raster_extract_obj_",out_suffix,".RData",sep=""))
  extract_obj_var_pred <-load_obj(extract_obj_fname)
  #df_points_day_extracted <- read.table(df_points_extracted_fname,sep=",")
  #df_time_series <- read.table( df_time_series,sep=",")
  #df_points_extracted <- read.table(df_points_extracted_fname,sep=",")
  #Ended at 16:50
  #so for 10,289 mosaics and 4,803 stations it took about 18h40 minutes
}

df_points_extracted_fname <-extract_obj_var_pred$df_points_extracted_fname
df_raster_fname <- extract_obj_var_pred$df_raster_fname
df_time_series_fname <- extract_obj_var_pred$df_time_series_fname

df_raster <- read.table(df_raster_fname,sep=",",header=T,stringsAsFactors =F)
df_time_series <- read.table( df_time_series_fname,sep=",",header=T,stringsAsFactors =F)
df_points_extracted <- read.table(df_points_extracted_fname,sep=",",header=T,stringsAsFactors = F)

#dim(df_time_series)

#################################
###### STEP 3: Combine extracted and observation with testing and training information ####
## Need to combine back with original data, this is done station by station... and year?

if(run_steps[3]==TRUE){

  
  k <- 1
  #data_var<- list_df_v_stations[[i]] 
  data_var <- read.table(list_data_stations_var_pred_df_filename[[k]],header=T, stringsAsFactors = F,sep=",")

  ### prepare arguments to combine stations
  list_selected_ID <- unique(df_points_extracted$id) #4800 stations selected
  #test_list<- table(df_points_extracted$id)
  #> test_list[test_list >1]
  #
  #71534099999 71689099999 71706099999 
  #          2           2           2 
          
  #data_var# This contains testing or training data with variable being modeled and covariates
  #list_data_stations_var_pred_df_filename

  in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaics/mosaic"
  #in_dir_mosaic_rmse <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaicsRMSE/mosaic"
  pattern_str <- ".*.tif"
  lf_raster <- list.files(pattern=pattern_str,path = in_dir_mosaic,full.names = T)
  #df_ts_pix <- df_points_extracted #this is extracted points with rows stations, columns are x,y, id and predictions from raster 
  r_ts_name <- sub(extension(lf_raster),"",basename(lf_raster))
  var_name <- "dailyTmax" #observed measurements, y_var_name
  var_pred <- "mod1" #predictions
  #dates_str <-
  #dates_val <-
  #df_raster #contains dates of raster mosaic produced
  df_time_series #contains de date , name, and dir for raster time series from date_star to end including missing dates field 
  plot_fig <- FALSE
  i<-1

  #Product assessment
  function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_01272017.R"
  source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script
  #out_suffix_str <- paste0(region_name,"_",out_suffix)
  #debug(combine_measurements_and_predictions_df)
  #this can be run with mclapply, very fast right now:
  #18:00pm
  station_summary_obj <- combine_measurements_and_predictions_df(
                                       i=2400,
                                       df_raster=df_raster,
                                       df_time_series=df_time_series,
                                       df_points_extracted=df_points_extracted,
                                       data_var=list_data_stations_var_pred_df_filename, #this can be a list
                                       list_selected_ID=list_selected_ID,
                                       r_ts_name=r_ts_name,
                                       var_name=var_name,
                                       var_pred = var_pred,
                                       scaling= scaling,
                                       out_dir =out_dir,
                                       out_suffix=out_suffix_str,
                                       plot_fig=F)

  #18:02pm
  ##Quick look at the data:
  df_pix_ts <- station_summary_obj$df_pix_ts
  plot(df_pix_ts$mod1_mosaic,type="l",col="blue")
  lines(df_pix_ts$dailyTmax,type="l",col="red")
  station_summary_obj$metric_stat_df

  #####
  ##combine information for the 1458 stations
  #started at 17.12pm
  list_station_summary_obj <- mclapply(1:length(list_selected_ID[2400:2410]),
                                        FUN=combine_measurements_and_predictions_df,
                                        df_raster=df_raster,
                                        df_time_series=df_time_series,
                                        df_points_extracted=df_points_extracted,
                                        data_var=list_data_stations_var_pred_df_filename, #this can be a list
                                        list_selected_ID=list_selected_ID[2400:2410],
                                        r_ts_name=r_ts_name,
                                        var_name=var_name,
                                        var_pred = var_pred,
                                        scaling=scaling,
                                        out_dir =out_dir,
                                        out_suffix=out_suffix_str,
                                        plot_fig=F,
                                        mc.preschedule=FALSE,
                                        mc.cores = num_cores)

  #station_summary_obj <- list(nb_zero,nb_NA, df_pix_ts)
  #check ID:70162026508
  test<- unlist(lapply(list_station_summary_obj, function(x) !inherits(x, "try-error")))
  #problem with 1379
  ####
}

##################################  END OF SCRIPT ########################################################
