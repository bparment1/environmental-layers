####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  NASA 2016 Meeting: biodiversity and ecological forecasting ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Combining tables and figures for individual runs for years and tiles.
#Figures and data for the AAG conference are also produced.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 05/01/2016  
#MODIFIED ON: 05/03/2016            
#Version: 1
#PROJECT: Environmental Layers project     
#COMMENTS: Initial commit, script based on part 2 of assessment, will be modified further for overall assessment 
#TODO:
#1) Add plot broken down by year and region 
#2) Modify code for overall assessment accross all regions and year
#3) Clean up

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rmd user:aguzman4:rwx /nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015

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

###### Function used in the script #######
  
#script_path <- "/nobackupp8/bparmen1/env_layers_scripts" #path to script
script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts" #path to script

## NASA poster and paper related
source(file.path(script_path,"NASA2016_conference_temperature_predictions_function_05032016b.R"))

#Mosaic related on NEX
#script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts"
function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_04232016.R" #PARAM12
function_mosaicing <-"global_run_scalingup_mosaicing_05012016.R"
source(file.path(script_path,function_mosaicing)) #source all functions used in this script 
source(file.path(script_path,function_mosaicing_functions)) #source all functions used in this script 

#Assessment on NEX
function_assessment_part1_functions <- "global_run_scalingup_assessment_part1_functions_02112015.R" #PARAM12
function_assessment_part1a <-"global_run_scalingup_assessment_part1a_01042016.R"
function_assessment_part2 <- "global_run_scalingup_assessment_part2_02092016.R"
function_assessment_part2_functions <- "global_run_scalingup_assessment_part2_functions_01032016.R"
function_assessment_part3 <- "global_run_scalingup_assessment_part3_04292016b.R"
source(file.path(script_path,function_assessment_part1_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part1a)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part3)) #source all functions used in this script 

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
in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/assessment"
in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics"

region_name <- c("reg4") #param 6, arg 3

create_out_dir_param <- TRUE #param 9, arg 6
out_suffix <- "_meeting_NASA_reg4_04292016"

out_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/assessment"

#run_figure_by_year <- TRUE # param 10, arg 7
list_year_predicted <- "1984,2014"

file_format <- ".tif" #format for mosaiced files # param 11
NA_flag_val <- -32768  #No data value, # param 12

#num_cores <- 6 #number of cores used # param 13, arg 8
plotting_figures <- TRUE #running part2 of assessment to generate figures... # param 14
num_cores <- 11 #number of cores used # param 13, arg 8
#python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin" #PARAM 30
python_bin <- "/usr/bin" #PARAM 30

day_start <- "1986101" #PARAM 12 arg 12
day_end <- "19981231" #PARAM 13 arg 13

#infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg4.tif"
infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg4.tif"

#run_figure_by_year <- TRUE # param 10, arg 7
list_year_predicted <- "1984,2014"
scaling <- 0.01 #was scaled on 100 

df_centroids_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/output_reg4_1999/df_centroids_19990701_reg4_1999.txt"

raster_name_lf <- c("/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920101_reg4_1992_m_gam_CAI_dailyTmax_19920101_reg4_1992.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920102_reg4_1992_m_gam_CAI_dailyTmax_19920102_reg4_1992.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920103_reg4_1992_m_gam_CAI_dailyTmax_19920103_reg4_1992.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920701_reg4_1992_m_gam_CAI_dailyTmax_19920701_reg4_1992.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920702_reg4_1992_m_gam_CAI_dailyTmax_19920702_reg4_1992.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920703_reg4_1992_m_gam_CAI_dailyTmax_19920703_reg4_1992.tif")

#l_dates <- c("19990101","19990102","19990103","19990701","19990702","19990703")
l_dates <- c("19920101","19920102","19920103","19920701","19920702","19990703")

df_points_extracted_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics/data_points_extracted.txt"
NA_flag_val_mosaic <- -3399999901438340239948148078125514752.000

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

#start_date <- day_to_mosaic_range[1]
#end_date <- day_to_mosaic_range[2]
#start_date <- day_start #PARAM 12 arg 12
#end_date <- day_end #PARAM 13 arg 13

#date_to_plot <- seq(as.Date(strptime(start_date,"%Y%m%d")), as.Date(strptime(end_date,"%Y%m%d")), 'day')
#l_dates <- format(date_to_plot,"%Y%m%d") #format back to the relevant date format for files
mask_pred <- TRUE
list_param_pre_process <- list(raster_name_lf,python_bin,infile_mask,scaling,mask_pred,NA_flag_val,out_suffix,out_dir) 
names(list_param_pre_process) <- c("lf","python_bin","infile_mask","scaling","mask_pred","NA_flag_val","out_suffix","out_dir") 
  
#debug(pre_process_raster_mosaic_fun)

#lf_mosaic_scaled <- mclapply(1:length(raster_name_lf),FUN=pre_process_raster_mosaic_fun,list_param=list_param_pre_process,mc.preschedule=FALSE,mc.cores = num_cores)                         
lf_mosaic_scaled <- mclapply(1:length(raster_name_lf),FUN=pre_process_raster_mosaic_fun,list_param=list_param_pre_process,mc.preschedule=FALSE,mc.cores = num_cores)                         

#test <- pre_process_raster_mosaic_fun(2,list_param_pre_process)
lf_mosaic_scaled <- unlist(lf_mosaic_scaled)

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

############### PART2: temporal profile #############
#### Extract time series
###
#-65,-22

df_points <- read.table(df_points_extracted_fname,sep=",") 
df_points_tmp <- df_points
df_points <- as.data.frame(t(df_points))
names(df_points) <- paste0("ID_",1:ncol(df_points))

#df_centroids <- read.table(df_centroids_fname,sep=",")

coordinates(df_centroids)<- c("long","lat")
proj4string(df_centroids) <- CRS_locs_WGS84

lf_mosaic_list <- list.files(path=in_dir_mosaic,pattern="*.tif")
#r_mosaic_ts <- stack(lf_mosaic_list)
#df_centroids <- extract(df_centroids,r_mosaic_ts)

df_points$files <- lf_mosaic_list

#debug(extract_date)
#test <- extract_date(6431,lf_mosaic_list,12) #extract item number 12 from the name of files to get the data

list_dates_produced <- unlist(mclapply(1:length(lf_mosaic_list),FUN=extract_date,x=lf_mosaic_list,item_no=12,mc.preschedule=FALSE,mc.cores = num_cores))                         
#list_dates_produced <-  mclapply(6400:6431,FUN=extract_date,x=lf_mosaic_list,item_no=12,mc.preschedule=FALSE,mc.cores = num_cores)                         

list_dates_produced_date_val <- as.Date(strptime(list_dates_produced,"%Y%m%d"))
month_str <- format(list_dates_produced_date_val, "%b") ## Month, char, abbreviated
year_str <- format(list_dates_produced_date_val, "%Y") ## Year with century
day_str <- as.numeric(format(list_dates_produced_date_val, "%d")) ## numeric month

df_points$date <- list_dates_produced_date_val
df_points$month <- month_str
df_points$year <- year_str
df_points$day <- day_str

unique_date_tb <-table(df_points$date)
unique_date <- unique(df_points$date)

station_id <- 8
var_name <-paste0("ID_",station_id)

##Screen for unique date values
if(max(unique_date_tb)>1){
#  formula_str <- paste(var_name," ~ ","TRIP_START_DATE_f",sep="")
   var_pix <- aggregate(ID_8 ~ date, data = df_points, mean) #aggregate by date
}

var_pix$ID_8 <- var_pix$ID_8*scaling

d_z_tmp <- zoo(var_pix$ID_8,var_pix$date)
names(d_z_tmp)<-"ID_8"
min(d_z_tmp$ID_8)
max(d_z_tmp$ID_8)

plot(d_z_tmp)

day_start <- "1986-01-01" #PARAM 12 arg 12
day_end <- "1998-12-31" #PARAM 13 arg 13

start_date <- as.Date(day_start)
end_date <- as.Date(day_end)
start_year <- year(start_date)
end_year <- year(end_date)

d_z <- window(d_z_tmp,start=start_date,end=end_date)
#d_z2 <- window(d_z_tmp2,start=start_date,end=end_date)

title_str <- paste("Predicted daily ",variable_name," for the ", start_year,"-",end_year," time period",sep="")
plot(d_z,ylab="tmax in deg C",xlab="daily time steps",
     main=title_str,
     lty=3)


#data_pixel <- data_df[id_selected,]
#data_pixel$rainfall <- as.numeric(data_pixel$rainfall)
#d_z_tmp <-zoo(data_pixel$rainfall,as.Date(data_pixel$date))
#names(d_z_tmp)<- "rainfall"
#data_pixel <- as.data.frame(data_pixel)
#d_z_tmp2 <- zoo(data_pixel[[var_name]],as.Date(data_pixel$date))
    
#df_tmp <- subset(data_var,data_var$ID_stat==id_name)
#if(da)


############################ END OF SCRIPT ##################################