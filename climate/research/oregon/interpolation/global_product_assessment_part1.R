####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Assessment of product part 1: mosaic and accuracy ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Combining tables and figures for individual runs for years and tiles.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 05/15/2016  
#MODIFIED ON: 09/19/2016            
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
#setfacl -Rmd user:aguzman4:rwx /nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015

#COMMIT: plotting extracted predicted values and measured tmax 

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
function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_09192016b.R"
source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 

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

out_region_name<-""
list_models<-c("y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)") #param 4

#reg1 (North Am), reg2(Europe),reg3(Asia), reg4 (South Am), reg5 (Africa), reg6 (Australia-Asia)
#master directory containing the definition of tile size and tiles predicted
#in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/assessment"
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/mosaic"
in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/assessment"
in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaic/mosaic"

region_name <- c("reg1") #param 6, arg 3
out_suffix <- "_global_assessment_reg1_08282016"

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
day_end <- "19991231" #PARAM 13 arg 13
#date_start <- day_start
#date_end <- day_end

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

############### PART1: Select stations for specific list of dates #############
#### Extract corresponding stationsfor given dates and/or plot stations used
## Use station from specific year and date?

list_dates_val <- as.Date(strptime(l_dates,"%Y%m%d"))
l_dates_month_str <- format(list_dates_val, "%b") ## Month, char, abbreviated
l_dates_year_str <- format(list_dates_val, "%Y") ## Year with century
l_dates_day_str <- as.numeric(format(list_dates_val, "%d")) ## numeric month


##start new function
### Needs to make this a function!!!
#Step 1 for a list of dates, load relevant files with year matching,
#Step 2 for giving dates subset the data.frame
#Step 3: find duplicates, create return object and return given station data for the date

list_year_str <- unique(l_dates_year_str)
list_df_v_stations <- vector("list",length=length(list_year_str))
list_df_s_stations <- vector("list",length=length(list_year_str))

for(i in 1:length(list_year_str)){
  filename_tmp<- df_data_v_fname$file[df_data_v_fname$year==list_year_str[i]]
  list_df_v_stations[[i]] <- read.table(filename_tmp,stringsAsFactors=F,sep=",")
  filename_tmp<- df_data_s_fname$file[df_data_s_fname$year==list_year_str[i]]
  list_df_s_stations[[i]] <- read.table(filename_tmp,stringsAsFactors=F,sep=",")
}

df_combined_data_v <- do.call(rbind,list_df_v_stations) #reading only the years related to the the dates e.g. 1999
df_combined_data_s <- do.call(rbind,list_df_s_stations)

#### Get df points for specific dates
#lapply(1:length(l_dates)list_df_v_stations,function(x){x$date==l_dates})
#dim(list_df_v_stations[[1]])
list_df_points_data_v_dates <- vector("list",length=length(l_dates))
list_df_points_data_s_dates <- vector("list",length=length(l_dates))
for(i in 1:length(l_dates)){
  #year_str <- list_year_str[[i]]
  list_df_points_data_v_dates[[i]] <- subset(df_combined_data_v,df_combined_data_v$date==l_dates[i])
  list_df_points_data_s_dates[[i]] <- subset(df_combined_data_s,df_combined_data_s$date==l_dates[i])

}

df_combined_data_v_dates <- do.call(rbind,list_df_points_data_v_dates)
df_combined_data_s_dates <- do.call(rbind,list_df_points_data_s_dates)

####


#function(x){x$date==}
#df_points <- subset(df_stations,df_stations_tmp$date==l_dates[1])
reg_layer <- readOGR(dsn=dirname(countries_shp),sub(".shp","",basename(countries_shp)))
countries_shp_tmp <- reg_layer

## Now plot by dates:

model_name <- "res_mod1"
var_selected <- "res_mod1"
variable_name
#8 minutes for 18 dates and reg1 ?
station_type_name <- "testing"
proj_str <- CRS_locs_WGS84
list_param_plot_stations_val_by_date <- list(l_dates,df_combined_data_v_dates,countries_shp_tmp,CRS_locs_WGS84,
                                             station_type_name,model_name,
                                             var_selected,y_var_name,out_suffix,out_dir)
names(list_param_plot_stations_val_by_date) <- c("l_dates","df_combined_data_points","countries_shp","proj_str",
                                                 "station_type_name","model_name",
                                                 "var_selected","y_var_name","out_suffix","out_dir")
  

list_plot_obj_data_v <- mclapply(1:length(l_dates),FUN=plot_stations_val_by_date,list_param=list_param_plot_stations_val_by_date,
                           mc.preschedule=FALSE,
                           mc.cores = num_cores)       

#### run for training now
station_type_name <- "training"
list_param_plot_stations_val_by_date <- list(l_dates,df_combined_data_s_dates,countries_shp_tmp,CRS_locs_WGS84,
                                             station_type_name,model_name,
                                             var_selected,y_var_name,out_suffix,out_dir)
names(list_param_plot_stations_val_by_date) <- c("l_dates","df_combined_data_points","countries_shp","proj_str",
                                                 "station_type_name","model_name",
                                                 "var_selected","y_var_name","out_suffix","out_dir")
#debug(plot_stations_val_by_date)
#test <- plot_stations_val_by_date(2,list_param = list_param_plot_stations_val_by_date)

list_plot_obj_data_s <- mclapply(1:length(l_dates),FUN=plot_stations_val_by_date,list_param=list_param_plot_stations_val_by_date,
                           mc.preschedule=FALSE,
                           mc.cores = num_cores)   
## do training

##plot summary of staistics
list_plot_obj_data_v[[1]]$df_basic_stat
names(list_plot_obj_data_v[[1]])
list_data_stations_var_pred_data_v <- lapply(list_plot_obj_data_v,FUN=function(x){x$data_stations_var_pred})
list_data_stations_var_pred_data_s <- lapply(list_plot_obj_data_s,FUN=function(x){x$data_stations_var_pred})

data_stations_var_pred_data_v <- do.call(rbind,list_data_stations_var_pred_data_v)
data_stations_var_pred_data_s <- do.call(rbind,list_data_stations_var_pred_data_s)
dim(data_stations_var_pred_data_v)
dim(data_stations_var_pred_data_s)

write.table(data_stations_var_pred_data_v,"data_stations_var_pred_data_v")
write.table(data_stations_var_pred_data_s,"data_stations_var_pred_data_s")


############### PART2: Select stations by ID to build a time series #############
#### Extract corresponding stationsfor given dates and plot stations used
## Use station from specific year and date?

###################### make this part a function!!!
#select one station based on id or coordinates and find that in the list of data.frame??

#Make a function to find the closest stations to a givine coordinates?
#42.262003, -71.965866 #this is near Spencer MA

#id_selected <- "82111099999"
#dim(df_points)
############### PART5: Extraction of temporal profile #############
#### Extract time series from raster stack
### Add dates to mkae it a date object?
#-72,42.24 #near Spencer MA
#48.281724, -71.401815 near Saguenay Quebec
#39.805293, -89.872626 near Springfield Illinois
#32.009676, -106.990266 near El Paso Us Texas
#39.955529, -105.263851: near Boulder Colorado
#45.384051, -121.859146 : near lost Creek near Mount Hood Oregon
#53.283685, -113.341702: leduc Canada near Edmonton
#39.009052, -77.026875: Silver Spring Sligo Creek, MD
#36.627806, -119.928901 : Raisin City, Fresno CA
#36.677139, -115.330122: Las Vegas
#35.696143, -78.396011: near Raleigh NC

list_lat_long <- list(
c( -72, 42.24),
c( -71.401815, 48.281724),
c( -89.872626, 39.805293),
c( -106.990266, 32.009676),
c( -105.263851, 39.955529), 
c( -121.859146,45.384051),
c( -113.341702,53.283685),
c( -77.026875,39.009052),
c( -119.928901,36.627806),
c(  -115.330122,36.677139),
c( -78.396011,35.696143)) 


test_day_query2 <- lapply(list_lat_long,FUN=query_for_station_lat_long,df_points_spdf=df_point_day,step_x=1,step_y=1)
#test_day_query <-query_for_station_lat_long(c(-72,42),df_points_spdf=df_point_day,step_x=1,step_y=0.25)
df_stations_selected <- do.call(rbind,test_day_query2)
proj4string(df_stations_selected) <- proj_str
#debug(query_for_station_lat_long)

##Next use the day results and select a mean station, quartile and min and max?

#list_id_data_v <- unique(data_stations_var_pred_data_v$id)
#list_id_data_s <- unique(data_stations_var_pred_data_s$id)

#Started at 4pm: on sept 9
list_id_data_v <- df_stations_selected$id
list_id_data_s <- df_stations_selected$id

### loop through all files and get the time series

lf_data_s_subset <- mclapply(list_data_s_fname,
                           FUN=extract_from_df,
                           col_selected="id",
                           val_selected=list_id_data_s,
                           mc.preschedule=FALSE,
                           mc.cores = num_cores)   
#took less than 8 minutes for 1839 stations

lf_data_v_subset <- mclapply(list_data_v_fname,
                           FUN=extract_from_df,
                           col_selected="id",
                           val_selected=list_id_data_v,
                           mc.preschedule=FALSE,
                           mc.cores = num_cores)   

data_v_subset <- do.call(rbind,lf_data_v_subset)
data_s_subset <- do.call(rbind,lf_data_s_subset)

data_s_subset$training <- 1
data_v_subset$training <- 0

## Need a testing variable to  count  later the use of  a station
data_s_subset$testing <- 0
data_v_subset$testing <- 1
# a station can be used multipel times as trainin gor testing within a day because of the overlap of tiles.
write.table(data_v_subset,"data_v_subset_test.txt")
write.table(data_s_subset,"data_s_subset_test.txt")
##finish at 16.10pm on 09/09


#data_stations <- rbind(data_s_subset,data_v_subset)
dim(data_s_subset)
#[1] 21991826       9
dim(data_v_subset)
#[1] 9319967      85

##36 minutes to get here
#rbind.fill(mtcars[c("mpg", "wt")], mtcars[c("wt", "cyl")])
data_stations <- rbind.fill(data_v_subset, data_s_subset)
#dim(data_stations)#one station only but repetition of records because of tiles and dates!!!
#[1] 31311793       90
dim(data_stations)#one station only but repetition of records because of tiles and dates!!!
#[1] 30202891       91

#coordinates(data_stations) <- cbind(data_stations$x,data_stations$y)
#proj4string(data_stations) <- CRS_locs_WGS84

#data_stations_var_pred <- aggregate(id ~ date, data = data_stations, min)
#data_stations_var_pred <- aggregate(id ~ x + y + date + dailyTmax + mod1 + res_mod1 , data = data_stations, min)

##Add tile id here...and check if data stations was training or testing.

#16.30 pm on 09/09
#data_stations_var_pred <- aggregate(id2 + date ~ x + y + dailyTmax + mod1 + res_mod1 ,data = data_stations, FUN=mean ) #+ mod1 + res_mod1 , data = data_stations, min)
dim(data_stations_var_pred)
#md <- melt(mydata, id=(c("id", "time")),)
md <- melt(data_stations, id=(c("id", "date")),measure.vars=c("x","y","dailyTmax","mod1","res_mod1"))
#formula_str <- "id + date ~ x + y + dailyTmax + mod1 + res_mod1"
data_stations_var_pred <- cast(md, id + date ~ variable, fun.aggregate = mean, 
  na.rm = TRUE)

#write.table(data_stations_var_pred,
#            file=file.path(out_dir,paste0("data_stations_var_pred_tmp_",out_suffix,".txt",
#                                                                 sep=",")))
write.table(data_stations_var_pred,
            file=file.path(out_dir,paste0("data_stations_var_pred_tmp_",out_suffix,".txt")),
            sep=",")

#data_stations_var_pred <- read.table(
#            file=file.path(out_dir,paste0("data_stations_var_pred_tmp_",out_suffix,".txt",
#                                                                 sep=",")))

md <- melt(data_stations, id=(c("id", "date")),measure.vars=c("training","testing"))
data_stations_training_testing <- cast(md, id + date ~ variable, fun.aggregate = sum, 
  na.rm = TRUE)

#write.table(data_stations_training_testing,
#            file=file.path(out_dir,paste0("data_stations_training_testing_",out_suffix,".txt",
#                                                                 sep=",")))
write.table(data_stations_training_testing,
            file=file.path(out_dir,paste0("data_stations_training_testing_",out_suffix,".txt")),
            sep=",")
#data_stations_var_pred <- aggregate(id2 ~ x + y + date + dailyTmax + mod1 + res_mod1 ,data = data_stations, mean ) #+ mod1 + res_mod1 , data = data_stations, min)
#data_stations$id2 <- as.numeric(data_stations$id)
#data_stations$date <- as.character(data_stations$date)

dim(data_stations_var_pred)
#> dim(data_stations_var_pred)
#[1] 57154     7
unique(data_stations_var_pred$id)
dim(data_stations_training_testing)
#[1] 57154     4

data_stations_var_pred$date_str <- data_stations_var_pred$date
data_stations_var_pred$date <- as.Date(strptime(data_stations_var_pred$date_str,"%Y%m%d"))

##Find stations used for training and testing
#data_stations_var_pred2 <- aggregate(id ~ training,data = data_stations, sum ) #+ mod1 + res_mod1 , data = data_stations, min)
#data_stations_var_pred2 <- aggregate(date ~ training,data = data_stations, sum ) #+ mod1 + res_mod1 , data = data_stations, min)

#data_stations_var_pred <- merge(data_stations_var_pred, data_stations_training_testing , by="id") #this is slow maybe do cbind?
#data_stations_var_pred_test <- data_stations_var_pred 

data_stations_var_pred <- cbind(data_stations_var_pred,data_stations_training_testing)

write.table(data_stations_var_pred,
            file=file.path(out_dir,paste0("data_stations_var_pred_",out_suffix,".txt")),
            sep=",")

#started at 16.51, 09/07


#####  This could be moved in a separate file!!
###############  PART4: Checking for mosaic produced for given region ##############
## From list of mosaic files predicted extract dates
## Check dates predicted against date range for a given date range
## Join file information to centroids of tiles data.frame

##Now grab year year 1992 or matching year...maybe put this in a data.frame??

#pattern_str <- "r_m_use_edge_weights_weighted_mean_mask_gam_CAI_dailyTmax_.*._reg4_.*.tif"
#searchStr = paste(in_dir_mosaic,"/output_reg4_2014",year_processed,"/gam_CAI_dailyTmax_predicted_",pred_mod_name,"*",day_to_mosaic[i],"*.tif",sep="")

#debug(extract_date)
#test <- extract_date(6431,lf_mosaic_list,12) #extract item number 12 from the name of files to get the data
#list_dates_produced <- unlist(mclapply(1:length(lf_mosaic_list),FUN=extract_date,x=lf_mosaic_list,item_no=13,mc.preschedule=FALSE,mc.cores = num_cores))                         
#list_dates_produced <-  mclapply(1:2,FUN=extract_date,x=lf_mosaic_list,item_no=13,mc.preschedule=FALSE,mc.cores = 2)                         

if(is.null(df_points_extracted_fname)){
  
  in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaics/mosaic"
  #in_dir_mosaic_rmse <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaicsRMSE/mosaic"
  pattern_str <- ".*.tif"
  #df_points <- #this contains the location of points to be used for extraction
  #df_points <- data_stations_var_pred_data_s #collected previously
    
  #extract from var pred mosaic, tmax in this case:
  extract_obj_var_pred <- extract_from_time_series_raster_stack(df_points,date_start,date_end,lf_raster=NULL,item_no=13,
                                      num_cores=11,pattern_str=NULL,in_dir=in_dir_mosaic,out_dir=out_dir,out_suffix=out_suffix)
  df_points_extracted_fname <-extract_obj_var_pred$df_points_extracted_fname
  df_raster_fname <- extract_obj_var_pred$df_raster_fname
  df_time_series_fname <- extract_obj_var_pred$df_time_series_fname

  df_raster <- read.table(df_points_extracted_fname,sep=",")
  df_time_series <- read.table( df_time_series,sep=",")
  df_points_extracted <- read.table(df_points_extracted_fname,sep=",")


}else{
  df_points_day_extracted <- read.table(df_points_extracted_fname,sep=",")
  df_time_series <- read.table( df_time_series,sep=",")
  df_points_extracted <- read.table(df_points_extracted_fname,sep=",")
}

#################################################################################################
###################### PART 5: combine stations data and extracted values ############################
#### Now combined with the station data extracted from the assessment stage

#combine
head(data_stations_var_pred)


#           id     date        x      y dailyTmax     mod1  res_mod1          id     date training testing
#1 71238099999 19930703 -112.867 53.683      20.2 18.19909 -2.000915 71238099999 19930703        4       0
#2 71238099999 19930704 -112.867 53.683      23.4 17.74476 -5.655237 71238099999 19930704        3       1
#3 71238099999 19930705 -112.867 53.683      21.7 19.22313 -2.476870 71238099999 19930705        3       1
#4 71238099999 19930706 -112.867 53.683      22.0 19.53294 -2.467061 71238099999 19930706        3       1
#5 71238099999 19930707 -112.867 53.683      20.9 17.84168 -3.058324 71238099999 19930707        2       2
#6 71238099999 19930708 -112.867 53.683      21.2 16.50887 -4.691128 71238099999 19930708        1       3

#need to combine:
#data_stations_var_pred: constains id, x, y, dailyTmax, mod1, res_mod1, date, training, testing
#df_time_series
#df_points_extracted #contains id,x,y of stations and extracted values from raster stack
dim(df_time_series)
#dim(data_stations_var_pred)
df_points_extracted_tmp <- df_points_extracted
#df_points_extracted <- cbind(df_points,df_points_extracted)

df_points_extracted$id <- df_points$id #this should have been done earlier in the extraction function
df_points_extracted$x <- df_points$x #this should have been done earlier in the extraction function
df_points_extracted$y <- df_points$y #this should have been done earlier in the extraction function

dim(df_time_series)

list_selected_ID <- unique(data_stations_var_pred$id) #11 stations selected
data_var <- data_stations_var_pred
df_ts_pix <- df_points_extracted
r_ts_name <- sub(extension(lf_raster),"",basename(lf_raster))
var_name <- "dailyTmax" #observed measurements
var_pred <- "mod1" #predictions
#dates_str <-
#dates_val <-
df_raster
df_time_series
plot_fig <- FALSE
i<-1

#Product assessment
#function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_09192016c.R"
#source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 
#undebug(combine_measurements_and_predictions_df)
out_suffix_str <- paste0(region_name,"_",out_suffix)
#this can be run with mclapply, very fast right now:
station_summary_obj <- combine_measurements_and_predictions_df(i=i,
                                        df_raster=df_raster,
                                        df_time_series=df_time_series,
                                        df_ts_pix=df_ts_pix,
                                        data_var=data_var,
                                        list_selected_ID=list_selected_ID,
                                        r_ts_name=r_ts_name,
                                        var_name=var_name,
                                        var_pred = var_pred,
                                        out_dir =out_dir,
                                        out_suffix=out_suffix_str,
                                        plot_fig=F)
df_pix_ts <- station_summary_obj$df_pix_ts
#station_summary_obj <- list(nb_zero,nb_NA, df_pix_ts)

#id_name <- list_selected_ID[i]
#df_pix_ts_filename <- file.path(out_dir,paste0("df_pix_ts_",id_name,out_suffix,y_var_name,".txt"))
#write.table(df_pix_ts,df_pix_ts_filename,sep=",")


##################### plot time series: make this a function!!! ####
###start of plotting
### makes this a function:

#pix_ts$date <- as.Date(strptime(pix_ts$date_str,"%Y%m%d"))
#head(pix_ts)
var_pred_mosaic <- paste0(var_pred,"_mosaic")
id_name <- list_selected_ID[i]
station_id <- id_name
#> myDF[ is.na(myDF) ] <- NA 

plot_observation_predictions_time_series <- function(df_pix_time_series,var_name,df_pix_time_series,station_id,scaling=NULL,out_dir=".",out_suffix=""){
  #This function plots time series for observed and predicted points.
  #This assumes that time series with residuals were extracted from raster stacks.
  #
  #Inputs:
  #1) df_pix_time_series
  #2) var_name1: name of variables use
  #3) var_name2:
  #4) var_name3: residuals time series (var_name2 -var_name 1)
  #5) station_id
  #6) scaling: if null no scaling
  #7) out_dir
  #8) out_suffix
  #
  
  ########## Start script #########
  
  ## Screen of NaN and make sure we have NA
  df_pix_ts[is.na(df_pix_ts)] <- NA
  
  #Scale if necessary
  if(!is.null(scaling)){
    
    df_pix_ts[[var_pred_mosaic]] <-  df_pix_ts[[var_pred_mosaic]]*scaling 
  }

  #### var_name 1
  d_z_obs <- zoo(as.numeric(df_pix_ts[[y_var_name]]),as.Date(df_pix_ts$date))
  plot(d_z_obs)  
  
  #### var_name 2
  d_z_var <- zoo(as.numeric(df_pix_ts[[var_pred_mosaic]]),as.Date(df_pix_ts$date)) #make sure date is a date object !!!
  #names(d_z_var) <- var_pred_mosaic
  plot(d_z_var)  

  d_z_diff <- d_z_var - d_z_obs 
  #d_z_res <- d_z_var - d_z_obs 
  
  #d_z_res <- zoo(as.numeric(df_pix_ts[[paste0("res_",var_pred_mosaic)]]),as.Date(df_pix_ts$date))
  plot(d_z_res)  

  d_z_all <- merge(d_z_var,d_z_obs,d_z_res)
  names(d_z_all) <- c("pred","obs","res")
  d_z <-  merge(d_z_var,d_z_obs)
  names(d_z) <- c("pred","obs")
  
  range(d_z_res,na.rm=T)
  mean(d_z_res,na.rm=T)
  quantile(d_z_res,na.rm=T)

  range_dates <- range(as.Date(df_pix_time_series$date))
  #day_start <- "1984-01-01" #PARAM 12 arg 12
  #day_end <- "2014-12-31" #PARAM 13 arg 13
  #start_date <- as.Date(day_start)
  #end_date <- as.Date(day_end)
  start_year <- year(range_dates[1])
  end_year <- year(range_dates[2])
  #range_year <- range(df_pix_ts$year)
  #start_year <- range_year[1]
  #end_year <- range_year[2]
  var="tmax"

  res_pix <- 1000
  #res_pix <- 480
  col_mfrow <- 2
  row_mfrow <- 1
  
  png_filename <-  file.path(out_dir,paste("Figure_5a_time_series_profile_",region_name,"_",id_name,y_var_name,out_suffix,".png",sep =""))
  title_str <- paste("Daily", var_name1,"and", var_name2,"for", station_id," for the ", start_year,"-",end_year," time period",sep=" ")
  
  png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
  #this is the whole time series
  #plot(d_z_var,ylab="tmax in deg C",xlab="Daily time steps",
  #     main=title_str,cex=3,font=2,
  #     cex.main=1.5,cex.lab=1.5,font.lab=2,
  #     lty=3)
  #range(df_pix_ts[[var_pred_mosaic]],na.rm=T)
   
  #plot(d_z,plot.type="single",col=c("blue","red"))
  col_str <- c("blue","red")
  plot(d_z,plot.type="single",col=col_str,
       ylab="tmax in deg C",xlab="Daily time steps",
       main=title_str,cex=1,font=2,type="l",
       cex.main=1.5,cex.lab=1.5,font.lab=2)
  
  #plot(df_pix_ts[[var_pred_mosaic]] ~ as.Date(df_pix_ts$date),ylab="tmax in deg C",xlab="Daily time steps",
  #   main=title_str,cex=1,font=2,type="l",
  #   cex.main=1.5,cex.lab=1.5,font.lab=2,
  #   lty=3)
  #lines(as.numeric(df_pix_ts[[y_var_name]]) ~ as.Date(df_pix_ts$date),ylab="tmax in deg C",xlab="Daily time steps",
  #   main=title_str,cex=1,font=2,col="red",type="l",
  #   cex.main=1.5,cex.lab=1.5,font.lab=2,
  #   lty=3)

  dev.off()

  #for list_windows go in a loop
  if(!is.null(list_windows)){
    
    for(i in 1:length(list_windows)){
      
      
      ### get smaller window

      window_range <- list_windows[i]
      #day_start <- "1994-01-01" #PARAM 12 arg 12
      #day_end <- "1999-12-31" #PARAM 13 arg 13
      day_start <- window_range[1]
      day_end <- windwow_range[2]
      
      start_date <- as.Date(day_start)
      end_date <- as.Date(day_end)
      start_year <- year(start_date)
      end_year <- year(end_date)
      
      ### now select from the original time series
      #
      d_z_w <- window(d_z_tmp,start=start_date,end=end_date)
      names(d_z) <- var_pred_mosaic
      df_var <- as.data.frame(d_z)
 
      d_z_tmp <- d_z_obs
      d_z_tmp <- d_z
      d_z_tmp <- window(d_z_tmp,start=start_date,end=end_date)
      plot(d_z_tmp,col=c("blue","red"),plot.type="single")
      df_z <- as.data.frame(d_z)

      names(d_z) <- y_var_name
      df_obs <- as.data.frame(d_z)
      names(d_z) <- y_var_name

      res_pix <- 1000
      #res_pix <- 480
      col_mfrow <- 2
      row_mfrow <- 1
  
      png_filename <-  file.path(out_dir,paste("Figure5b_time_series_profile_",region_name,"_",out_suffix,".png",sep =""))
      title_str <- paste("Predicted daily ", station_id," ",var," for the ", start_year,"-",end_year," time period",sep="")
  
      png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)

      plot(df_pix_ts[[var_pred_mosaic]] ~ as.Date(df_pix_ts$date),ylab="tmax in deg C",xlab="Daily time steps",
      main=title_str,cex=1,font=2,type="l",
      cex.main=1.5,cex.lab=1.5,font.lab=2,
      lty=3)
      lines(as.numeric(df_pix_ts[[y_var_name]]) ~ as.Date(df_pix_ts$date),ylab="tmax in deg C",xlab="Daily time steps",
      main=title_str,cex=1,font=2,col="red",type="l",
      cex.main=1.5,cex.lab=1.5,font.lab=2,
      lty=3)

      dev.off()

    }
    
    
  }
  
}

############### PART5: Make raster stack and display maps #############
#### Extract corresponding raster for given dates and plot stations used


#start_date <- day_to_mosaic_range[1]
#end_date <- day_to_mosaic_range[2]
#start_date <- day_start #PARAM 12 arg 12
#end_date <- day_end #PARAM 13 arg 13

#date_to_plot <- seq(as.Date(strptime(start_date,"%Y%m%d")), as.Date(strptime(end_date,"%Y%m%d")), 'day')
#l_dates <- format(date_to_plot,"%Y%m%d") #format back to the relevant date format for files
#mask_pred <- FALSE
#matching <- FALSE #to be added after mask_pred option
#list_param_pre_process <- list(raster_name_lf,python_bin,infile_mask,scaling,mask_pred,NA_flag_val,out_suffix,out_dir) 
#names(list_param_pre_process) <- c("lf","python_bin","infile_mask","scaling","mask_pred","NA_flag_val","out_suffix","out_dir") 
  
#debug(pre_process_raster_mosaic_fun)

#lf_mosaic_scaled <- mclapply(1:length(raster_name_lf),FUN=pre_process_raster_mosaic_fun,list_param=list_param_pre_process,mc.preschedule=FALSE,mc.cores = num_cores)                         
#lf_mosaic_scaled <- mclapply(1:length(raster_name_lf),FUN=pre_process_raster_mosaic_fun,list_param=list_param_pre_process,mc.preschedule=FALSE,mc.cores = num_cores)                         

#test <- pre_process_raster_mosaic_fun(2,list_param_pre_process)
#lf_mosaic_scaled <- unlist(lf_mosaic_scaled)

##################################### PART 5  ######
##### Plotting specific days for the mosaics

r_mosaic_scaled <- stack(lf_mosaic_scaled)
NAvalue(r_mosaic_scaled)<- -3399999901438340239948148078125514752.000
plot(r_mosaic_scaled,y=6,zlim=c(-50,50))
plot(r_mosaic_scaled,zlim=c(-50,50),col=matlab.like(255))

#layout_m<-c(1,3) #one row two columns
#levelplot(r_mosaic_scaled,zlim=c(-50,50),col.regions=matlab.like(255))
#levelplot(r_mosaic_scaled,zlim=c(-50,50),col.regions=matlab.like(255))

#png(paste("Figure7a__spatial_pattern_tmax_prediction_levelplot_",date_selected,out_prefix,".png", sep=""),
#    height=480*layout_m[1],width=480*layout_m[2])
#plot(r_pred,col=temp.colors(255),zlim=c(-3500,4500))
#plot(r_pred,col=matlab.like(255),zlim=c(-40,50))
#paste(raster_name[1:7],collapse="_")
#add filename option later

#NA_flag_val_mosaic <- -3399999901438340239948148078125514752.000

list_param_plot_raster_mosaic <- list(l_dates,r_mosaic_scaled,NA_flag_val_mosaic,out_dir,out_suffix,
                                      region_name,variable_name)
names(list_param_plot_raster_mosaic) <- c("l_dates","r_mosaic_scaled","NA_flag_val_mosaic","out_dir","out_suffix",
                                          "region_name","variable_name")

lf_mosaic_plot_fig <- mclapply(1:length(lf_mosaic_scaled),FUN=plot_raster_mosaic,list_param=list_param_plot_raster_mosaic,mc.preschedule=FALSE,mc.cores = num_cores)                         

####################
###### Now add figures with additional met stations?

#df_selected2 <- data_stations_var_pred

#d_z_tmp <- zoo(df_selected[[station_id]],df_selected$date)
#d_z_tmp2 <- zoo(df_selected2$dailyTmax,df_selected2$date)

#d_z_tmp <- zoo(df[[station_id]],df$date)
#names(d_z_tmp)<-station_id
#names(d_z_tmp2)<-station_id

#min(d_z_tmp$ID_8)
#max(d_z_tmp$ID_8)
#day_start <- "1984-01-01" #PARAM 12 arg 12
#day_end <- "2014-12-31" #PARAM 13 arg 13
#day_start <- "1991-01-01" #PARAM 12 arg 12
#day_end <- "1992-12-31" #PARAM 13 arg 13

#start_date <- as.Date(day_start)
#end_date <- as.Date(day_end)
#start_year <- year(start_date)
#end_year <- year(end_date)

#res_pix <- 1000
#res_pix <- 480
#col_mfrow <- 2
#row_mfrow <- 1
  
#png_filename <-  file.path(out_dir,paste("Figure5a_time_series_profile_",region_name,"_",out_suffix,".png",sep =""))
#title_str <- paste("Predicted daily ", station_id," ",var," for the ", start_year,"-",end_year," time period",sep="")
  
#png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
#this is the whole time series

#plot(d_z_tmp,ylab="tmax in deg C",xlab="Daily time steps",
#     main=title_str,cex=3,font=2,
#     cex.main=1.5,cex.lab=1.5,font.lab=2,
#     lty=3)
#lines(d_z_tmp2,ylab="tmax in deg C",xlab="Daily time steps",
#     main=title_str,cex=3,font=2,
#     col="red",
#     cex.main=1.5,cex.lab=1.5,font.lab=2,
#     lty=3)
#
#dev.off()

#This is a lot of replication!!! okay cut it down
data_stations$date <- as.Date(strptime(data_stations_var_pred$date_str,"%Y%m%d"))
data_stations_tmp <- data_stations
data_stations_tmp <- data_stations

data_stations_tmp$date <- as.Date(strptime(data_stations_tmp$date,"%Y%m%d"))
#data_stations_tmp$date <- as.Date(strptime(data_stations_tmp$date_str,"%Y%m%d"))
#d_z4 <- 
d_z_tmp4 <- zoo(data_stations_tmp$dailyTmax,data_stations_tmp$date)
plot(d_z_tmp,cex.main=1.5,cex.lab=1.5,font.lab=2,
     lty=3)
lines(d_z_tmp4,ylab="tmax in deg C",xlab="Daily time steps",
     main=title_str,cex=3,font=2,
     col="red",
     cex.main=1.5,cex.lab=1.5,font.lab=2,
     lty=3)

day_start <- "1991-01-01" #PARAM 12 arg 12
day_end <- "1992-12-31" #PARAM 13 arg 13

start_date <- as.Date(day_start)
end_date <- as.Date(day_end)
start_year <- year(start_date)
end_year <- year(end_date)
d_z4 <- window(d_z_tmp4,start=start_date,end=end_date)

data_stations_tmp$date[7190]

###TO DO:
#1) compute correlation r, RMSE, MAE etc by station and profile
#2) Compute temporal autocorrelation by station and profile
#3) 

#### PLOT ACCURACY METRICS: First test ####
##this will be cleaned up later:

#dir_ac_mosaics <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/output_reg4_1999"
lf_tmp <-list.files(path=dir_ac_mosaics,pattern="r_m_use_edge_weights_weighted_mean_mask_gam_CAI_.*.ac.*._reg4_1999.tif")

#lf_tmp1 <- lf_tmp[21:24]
#list_param_plot_raster_mosaic
lf_tmp <-list.files(path=dir_ac_mosaics,pattern="r_m_use_edge_weights_weighted_mean_mask_gam_CAI_.*.ac.*._reg4_1999.tif",full.names=T)
#Product assessment
#function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_06142016b.R"
#source(file.path(script_path,function_product_assessment_part1_functions)) #source all functions used in this script 

r_mosaiced_ac <- stack(lf_tmp)
l_dates <- unlist(lapply(1:length(lf_tmp),FUN=extract_date,x=basename(lf_tmp),item_no=14))
variable_name
zlim_val <- NULL
list_param_plot_raster_mosaic <- list(l_dates,r_mosaiced_ac,NA_flag_val_mosaic,out_dir,out_suffix,
                                      region_name,variable_name, zlim_val)
names(list_param_plot_raster_mosaic) <- c("l_dates","r_mosaiced_scaled","NA_flag_val_mosaic","out_dir","out_suffix",
                                          "region_name","variable_name","zlim_val")
#debug(plot_raster_mosaic)
plot_raster_mosaic(1,list_param_plot_raster_mosaic)
lf_mosaic_plot_fig <- mclapply(1:length(lf_tmp),
                               FUN=plot_raster_mosaic,
                               list_param=list_param_plot_raster_mosaic,
                               mc.preschedule=FALSE,
                               mc.cores = num_cores)                         

#### Now plot kriged residuals from mosaiced surfaces

lf_tmp_res <-list.files(path=dir_ac_mosaics,pattern="r_m_use_edge_weights_weighted_mean_mask_gam_CAI_.*.residuals.*._reg4_1999.tif",full.names=T)

l_dates <- unlist(lapply(1:length(lf_tmp_res),FUN=extract_date,x=basename(lf_tmp),item_no=14))
variable_name
zlim_val <- NULL
r_mosaiced_res <- stack(lf_tmp_res)
list_param_plot_raster_mosaic <- list(l_dates,r_mosaiced_res,NA_flag_val_mosaic,out_dir,out_suffix,
                                      region_name,variable_name, zlim_val)
names(list_param_plot_raster_mosaic) <- c("l_dates","r_mosaiced_scaled","NA_flag_val_mosaic","out_dir","out_suffix",
                                          "region_name","variable_name","zlim_val")
#debug(plot_raster_mosaic)
plot_raster_mosaic(1,list_param_plot_raster_mosaic)
lf_mosaic_plot_fig_res <- mclapply(1:length(lf_tmp_res),
                               FUN=plot_raster_mosaic,
                               list_param=list_param_plot_raster_mosaic,
                               mc.preschedule=FALSE,
                               mc.cores = num_cores)                         

### New plot of residuals surface with zlim
zlim_val <- c(-60,60)
#r_mosaiced_res <- stack(lf_tmp_res)
list_param_plot_raster_mosaic <- list(l_dates,r_mosaiced_res,NA_flag_val_mosaic,out_dir,out_suffix,
                                      region_name,variable_name, zlim_val)
names(list_param_plot_raster_mosaic) <- c("l_dates","r_mosaiced_scaled","NA_flag_val_mosaic","out_dir","out_suffix",
                                          "region_name","variable_name","zlim_val")
#debug(plot_raster_mosaic)
#plot_raster_mosaic(1,list_param_plot_raster_mosaic)
lf_mosaic_plot_fig_res <- mclapply(1:length(lf_tmp_res),
                               FUN=plot_raster_mosaic,
                               list_param=list_param_plot_raster_mosaic,
                               mc.preschedule=FALSE,
                               mc.cores = num_cores)                         


############################ END OF SCRIPT ##################################
