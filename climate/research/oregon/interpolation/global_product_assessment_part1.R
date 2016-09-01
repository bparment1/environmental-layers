####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Assessment of product part 1: mosaic and accuracy ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Combining tables and figures for individual runs for years and tiles.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 05/15/2016  
#MODIFIED ON: 08/31/2016            
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
function_assessment_part3 <- "global_run_scalingup_assessment_part3_05162016.R"
source(file.path(script_path,function_assessment_part1_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part1a)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part3)) #source all functions used in this script 

#Product assessment
function_product_assessment_part1_functions <- "global_product_assessment_part1_functions_08312016.R"
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
day_end <- "19981231" #PARAM 13 arg 13

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

raster_name_lf <- c("/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920101_reg4_1992_m_gam_CAI_dailyTmax_19920101_reg4_1992.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920102_reg4_1992_m_gam_CAI_dailyTmax_19920102_reg4_1992.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920103_reg4_1992_m_gam_CAI_dailyTmax_19920103_reg4_1992.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920701_reg4_1992_m_gam_CAI_dailyTmax_19920701_reg4_1992.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920702_reg4_1992_m_gam_CAI_dailyTmax_19920702_reg4_1992.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920703_reg4_1992_m_gam_CAI_dailyTmax_19920703_reg4_1992.tif")

#l_dates <- c("19990101","19990102","19990103","19990701","19990702","19990703")
l_dates <- c("19840110","19840111","19840112","19840113","19840114") #dates to plot and analze

#df_points_extracted_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/int_mosaics/data_points_extracted.txt"
df_points_extracted_fname <- NULL #if null compute on the fly
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
df_data_v_fname <- data.frame(year_str,region_name,dataset="data_v",file=list_data_s_fname)

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

df_combined <- do.call(rbind,list_df_v_stations)

#### Get df points for specific dates
#lapply(1:length(l_dates)list_df_v_stations,function(x){x$date==l_dates})
#dim(list_df_v_stations[[1]])
list_df_points_dates <- vector("list",length=length(l_dates))
for(i in 1:length(l_dates)){
  #year_str <- list_year_str[[i]]
  list_df_points_dates[[i]] <- subset(df_combined,df_combined$date==l_dates[i])
}

df_combined_dates <- do.call(rbind,list_df_points_dates)

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

for(i in 1:length(l_dates)){
  #d

    
  date_processed <- l_dates[i]
  
  df_points <- list_df_points_dates[[i]]
  #plot(df_points)
  table(df_points$tile_id)
  #plot(df_points,add=T)
  coordinates(df_points) <- cbind(df_points$x,df_points$y)
  proj4string(df_points) <- CRS_locs_WGS84
  # # No spatial duplicates
  df_points_day <- remove.duplicates(df_points) #remove duplicates...
  #plot(df_points_day)
  dim(df_points_day)
  dim(df_points)
  #plot(df_points)
  
  ##layer to be clipped
  if(class(countries_shp)=!"SpatialPolygonsDataFrame"){
    reg_layer <- readOGR(dsn=dirname(countries_shp),sub(".shp","",basename(countries_shp)))
  }

  #extent_df_points_day <- extent(df_points_day)
  
  reg_layer_clipped <- gClip(shp=reg_layer, bb=df_points_day , keep.attribs=TRUE,outDir=NULL,outSuffix=NULL)
  
  #data_stations_temp <- aggregate(id ~ date, data = data_stations, min)
  #data_stations_temp <- aggregate(id ~ x + y + date + dailyTmax + mod1 + res_mod1 , data = data_stations, min)
  data_stations_temp <- aggregate(id ~ x + y + date + dailyTmax + res_mod1,data = df_points, mean ) #+ mod1 + res_mod1 , data = data_stations, min)
  dim(data_stations_temp)
  #> dim(data_stations_temp)
  #[1] 11171     5

  data_stations_temp$date_str <- data_stations_temp$date
  data_stations_temp$date <- as.Date(strptime(data_stations_temp$date_str,"%Y%m%d"))
  coordinates(data_stations_temp) <- cbind(data_stations_temp$x,data_stations_temp$y)
  data_stations_temp$constant <- c(1000,2000)
  #plot(data_stations_temp)
  #plot(reg_layer)
  #res_pix <- 1200
  res_pix <- 400
  col_mfrow <- 1
  row_mfrow <- 1
  
  png_filename <- paste("Figure_ac_metrics_map_stations_l_",station_type_name,"_",model_name,"_",y_var_name,"_",date_processed,
                       "_",out_suffix,".png",sep="")
  png(filename=png_filename,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  #plot(data_stations_temp
  #p_shp <- spplot(reg_layer_clipped,"ISO3" ,col.regions=NA, col="black") #ok problem solved!!
  #title("(a) Mean for 1 January")
  #p <- bubble(data_stations_temp,"constant",main=paste0("Average Residuals by validation stations ",
  #                                                      date_processed,
  #                                                      "for ",var_selected))
  #p <- spplot(data_stations_temp,"constant",col.regions=NA, col="black",
  #            main=paste0("Average Residuals by validation stations ",pch=3,cex=10,
  #            date_processed, "for ",var_selected))

  #p1 <- p+p_shp
  #try(print(p1)) #error raised if number of missing values below a threshold does not exist
  
  title_str <- paste0("Average Residuals by ",station_type_name," stations ",date_processed, " for ",y_var_name," ", var_selected)
  plot(reg_layer_clipped,main=title_str)
  plot(data_stations_temp,add=T,cex=0.5,col="blue")
  legend("topleft",legend=paste("n= ", nrow(data_stations_temp)),bty = "n")
  
  dev.off()
  
  res_pix <- 800
  col_mfrow <- 1
  row_mfrow <- 1
  png_filename <- paste("Figure_ac_metrics_map_stations",station_type_name,"averaged",model_name,y_var_name,date_processed,
                       out_suffix,".png",sep="_")
  png(filename=png_filename,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
  #model_name[j]
  title_str <- paste0("Average Residuals by ",station_type_name," stations ",date_processed, " for ",y_var_name," ", var_selected)
  #p_shp <- layer(sp.polygons(reg_layer, lwd=1, col='black'))
  p_shp <- spplot(reg_layer,"ISO3" ,col.regions=NA, col="black") #ok problem solved!!
  #title("(a) Mean for 1 January")
  p <- bubble(data_stations_temp,"res_mod1",main=title_str)
  p1 <- p+p_shp
  try(print(p1)) #error raised if number of missing values below a threshold does not exist
  
  dev.off()

  #### histogram of values
  res_pix <- 800
  col_mfrow <- 1
  row_mfrow <- 1
  png_filename <- paste("Figure_ac_metrics_histograms_stations",station_type_name,"averaged",model_name,y_var_name,date_processed,
                       out_suffix,".png",sep="_")
  png(filename=png_filename,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  h<- histogram(data_stations_temp$res_mod1,breaks=seq(-50,50,5))
  
  print(h)
  dev.off()
  
  ##Make data.frame with dates for later use!!
  #from libary mosaic
  df_basic_stat <- fav_stats(data_stations_temp$res_mod1)
  #quantile(data_stations_temp$res_mod1,c(1,5,10,90,95,99))
  df_quantile_val <- (quantile(data_stations_temp$res_mod1,c(0.01,0.05,0.10,0.45,0.50,0.55,0.90,0.95,0.99)))
  #quantile(c(1,5,10,90,95,99),data_stations_temp$res_mod1,)
  #rmse(data_stations_temp$res_mod1)
}



############### PART2: Select stations by ID to build a time series #############
#### Extract corresponding stationsfor given dates and plot stations used
## Use station from specific year and date?

#####################
#select one station based on id or coordinates and find that in the list of data.frame??

id_selected <- "82111099999"
dim(df_points)

### loop through all files and get the time series
extract_from_df <- function(x,col_selected,val_selected){
  df_tmp <- read.table(x,stringsAsFactors=F,sep=",")
  #data_subset <- subset(data_stations,col_selected==val_selected)
  data_subset <- subset(df_tmp,df_tmp$id%in%val_selected)
  return(data_subset)
}

lf_data_s_subset <- mclapply(list_data_s_fname,
                           FUN=extract_from_df,
                           col_selected="id",
                           val_selected=id_selected,
                           mc.preschedule=FALSE,
                           mc.cores = num_cores)                         
lf_data_v_subset <- mclapply(list_data_v_fname,
                           FUN=extract_from_df,
                           col_selected="id",
                           val_selected=id_selected,
                           mc.preschedule=FALSE,
                           mc.cores = num_cores)   

data_v_subset <- do.call(rbind,lf_data_v_subset)
data_s_subset <- do.call(rbind,lf_data_s_subset)

data_s_subset$training <- 1
data_v_subset$training <- 0

data_stations <- rbind(data_s_subset,data_v_subset)
dim(data_s_subset)
dim(data_v_subset)

#rbind.fill(mtcars[c("mpg", "wt")], mtcars[c("wt", "cyl")])
data_stations <- rbind.fill(data_v_subset, data_s_subset)

coordinates(data_stations) <- cbind(data_stations$x,data_stations$y)
proj4string(data_stations) <- CRS_locs_WGS84

dim(data_stations) #one station only but repitition of records because of tiles and dates!!!
#> dim(data_stations)
#[1] 100458     90

#data_stations_temp <- aggregate(id ~ date, data = data_stations, min)
#data_stations_temp <- aggregate(id ~ x + y + date + dailyTmax + mod1 + res_mod1 , data = data_stations, min)
data_stations_temp <- aggregate(id ~ x + y + date + dailyTmax,data = data_stations, min ) #+ mod1 + res_mod1 , data = data_stations, min)
dim(data_stations_temp)
#> dim(data_stations_temp)
#[1] 11171     5

data_stations_temp$date_str <- data_stations_temp$date
data_stations_temp$date <- as.Date(strptime(data_stations_temp$date_str,"%Y%m%d"))

##Find stations used for training and testing
#data_stations_temp2 <- aggregate(id ~ training,data = data_stations, sum ) #+ mod1 + res_mod1 , data = data_stations, min)
#data_stations_temp2 <- aggregate(date ~ training,data = data_stations, sum ) #+ mod1 + res_mod1 , data = data_stations, min)

############### PART3: Make raster stack and display maps #############
#### Extract corresponding raster for given dates and plot stations used

##Now grab year year 1992 or matching year...maybe put this in a data.frame??

if(is.null(r_mosaic_fname)){
  pattern_str <-"*.tif"
  lf_mosaic_list <- list.files(path=in_dir_mosaic,pattern=pattern_str,recursive=F,full.names=T)
  r_mosaic <- stack(lf_mosaic_list)
  save(r_mosaic,file="r_mosaic.RData")
}else{
  r_mosaic <- load_obj(r_mosaic_fname) #load raster stack of images
}

#start_date <- day_to_mosaic_range[1]
#end_date <- day_to_mosaic_range[2]
#start_date <- day_start #PARAM 12 arg 12
#end_date <- day_end #PARAM 13 arg 13

#date_to_plot <- seq(as.Date(strptime(start_date,"%Y%m%d")), as.Date(strptime(end_date,"%Y%m%d")), 'day')
#l_dates <- format(date_to_plot,"%Y%m%d") #format back to the relevant date format for files
mask_pred <- FALSE
matching <- FALSE #to be added after mask_pred option
list_param_pre_process <- list(raster_name_lf,python_bin,infile_mask,scaling,mask_pred,NA_flag_val,out_suffix,out_dir) 
names(list_param_pre_process) <- c("lf","python_bin","infile_mask","scaling","mask_pred","NA_flag_val","out_suffix","out_dir") 
  
#debug(pre_process_raster_mosaic_fun)

lf_mosaic_scaled <- mclapply(1:length(raster_name_lf),FUN=pre_process_raster_mosaic_fun,list_param=list_param_pre_process,mc.preschedule=FALSE,mc.cores = num_cores)                         
#lf_mosaic_scaled <- mclapply(1:length(raster_name_lf),FUN=pre_process_raster_mosaic_fun,list_param=list_param_pre_process,mc.preschedule=FALSE,mc.cores = num_cores)                         

#test <- pre_process_raster_mosaic_fun(2,list_param_pre_process)
#lf_mosaic_scaled <- unlist(lf_mosaic_scaled)

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

###############  PART4: Checking for mosaic produced for given region ##############
## From list of mosaic files predicted extract dates
## Check dates predicted against date range for a given date range
## Join file information to centroids of tiles data.frame

## Checking new files:
#in_dir_mosaic <- "/nobackupp6/aguzman4/climateLayers/out/reg4/mosaics2/mosaic"
#/nobackupp6/aguzman4/climateLayers/out/reg4/mosaics2/mosaic/output_reg4_*/r_m_use_edge_weights_weighted_mean_mask_gam_CAI_dailyTmax_*_reg4_*.tif
pattern_str <- "r_m_use_edge_weights_weighted_mean_mask_gam_CAI_dailyTmax_.*._reg4_.*.tif"
searchStr = paste(in_dir_mosaic,"/output_reg4_2014",year_processed,"/gam_CAI_dailyTmax_predicted_",pred_mod_name,"*",day_to_mosaic[i],"*.tif",sep="")

#lf_mosaic_list <- list.files(path=in_dir_mosaic,pattern="*.tif",recursive=T)
lf_mosaic_list <- list.files(path=in_dir_mosaic,pattern=pattern_str,recursive=T)
lf_mosaic_list <- lapply(1:length(day_to_mosaic),
                         FUN=function(i){
                           searchStr = paste(in_dir_tiles_tmp,"/*/",year_processed,"/gam_CAI_dailyTmax_predicted_",pred_mod_name,"*",day_to_mosaic[i],"*.tif",sep="")
                           Sys.glob(searchStr)})

#r_mosaic_ts <- stack(lf_mosaic_list)
#df_centroids <- extract(df_centroids,r_mosaic_ts)

df_points$files <- lf_mosaic_list

#debug(extract_date)
#test <- extract_date(6431,lf_mosaic_list,12) #extract item number 12 from the name of files to get the data
list_dates_produced <- unlist(mclapply(1:length(lf_mosaic_list),FUN=extract_date,x=lf_mosaic_list,item_no=14,mc.preschedule=FALSE,mc.cores = num_cores))                         
#list_dates_produced <-  mclapply(6400:6431,FUN=extract_date,x=lf_mosaic_list,item_no=12,mc.preschedule=FALSE,mc.cores = num_cores)                         

list_dates_produced_date_val <- as.Date(strptime(list_dates_produced,"%Y%m%d"))
month_str <- format(list_dates_produced_date_val, "%b") ## Month, char, abbreviated
year_str <- format(list_dates_produced_date_val, "%Y") ## Year with century
day_str <- as.numeric(format(list_dates_produced_date_val, "%d")) ## numeric month

df_produced <- data.frame(lf_mosaic_list,list_dates_produced_date_val,month_str,year_str,day_str)

date_start <- "19840101"
date_end <- "20141231"
date1 <- as.Date(strptime(date_start,"%Y%m%d"))
date2 <- as.Date(strptime(date_end,"%Y%m%d"))
dates_range <- seq.Date(date1, date2, by="1 day") #sequence of dates

missing_dates <- setdiff(as.character(dates_range),as.character(list_dates_produced_date_val))

month_str <- format(list_dates_produced_date_val, "%b") ## Month, char, abbreviated
year_str <- format(list_dates_produced_date_val, "%Y") ## Year with century
day_str <- as.numeric(format(list_dates_produced_date_val, "%d")) ## numeric month

df_points$date <- list_dates_produced_date_val
df_points$month <- month_str
df_points$year <- year_str
df_points$day <- day_str

############### PART5: Extraction of temporal profile #############
#### Extract time series from raster stack
### Add dates to mkae it a date object?
#-65,-22

#Use the global output??

##23.09 (on 05/22)
#df_points_day_extracted <- extract(r_mosaic,data_stations,df=T)
#df_points_day_extracted_fname <- paste0("df_points_day_extracted_fname",".txt") 
#write.table(df_points_day_extracted,file=df_points_day_extracted_fname) #5.1 Giga
#4.51 (on 05/23)

if(is.null(df_points_extracted_fname)){
  
  ##10.41 (on 05/22)
  df_points_day_extracted <- extract(r_mosaic,df_points_day,df=T)
  #17.19 (on 05/23)
  df_points_day_extracted_fname <- paste0("df_points_day_extracted_fname2",".txt")
  #17.27 (on 05/23)
  write.table(df_points_day_extracted,file=df_points_day_extracted_fname,sep=",",row.names = F) 
  #17.19 (on 05/23)

}else{
  df_points_day_extracted <- read.table(df_points_extracted_fname,sep=",")
}

head(df_points_day) #contains ID, dates and coordinates
df_points_day$id

pix_ts <- as.data.frame(t(df_points_day_extracted))
pix_ts <- pix_ts[-1,]
#var_names <- rownames(df_points_day_extracted) #same as lf_mosaic_list but without (*.tif)
var_names <- rownames(pix_ts) #same as lf_mosaic_list but without (*.tif)
var_id <- df_points_day$id
df_points_day_extracted$id <- var_id
 
#lf_var <- names(r_mosaic)
### Not get the data from the time series
#data_pixel <- df_ts_pix[id_selected,]
#data_pixel <- as.data.frame(data_pixel)
#pix_ts <- t(as.data.frame(subset(data_pixel,select=r_ts_name))) #can subset to range later
#pix_ts <- subset(as.data.frame(pix_ts),select=r_ts_name)

list_dates_produced <- unlist(mclapply(1:length(var_names),
                                       FUN=extract_date,
                                       x=var_names,
                                       item_no=12,
                                       mc.preschedule=FALSE,
                                       mc.cores = num_cores))                         

list_dates_produced_date_val <- as.Date(strptime(list_dates_produced,"%Y%m%d"))
pix_ts$date <- list_dates_produced_date_val
pix_ts[,1]*scaling #scale?

names(pix_ts) <- var_id

pix_ts <- read.table("pix_ts2.txt",sep=",")
names(pix_ts) <- c(var_id,"files","date_str")
pix_ts$date <- as.Date(strptime(pix_ts$date_str,"%Y%m%d"))
  
#head(pix_ts)

###start of plotting
### makes this a function:
id_selected <- "82111099999"
#station_id <- 8
station_id <- id_selected
df <- pix_ts
#scaling
#start_date: subset dates
#end_date: subset dates
df2 <- data_stations_temp
  
  
df_selected <- subset(df,select=station_id)
#df_selected <- subset(pix_ts,select=station_id)
names(df_selected) <- station_id
df_selected$date <- df$date

if(!is.null(scaling)){
  df_selected[[station_id]]<- df_selected[[station_id]]*scaling
}
if(!is.null(df2)){
  df_selected2 <- df2
  rm(df2)
  d_z_tmp2 <- zoo(df_selected2$dailyTmax,df_selected2$date)
  names(d_z_tmp2)<-station_id
}else{
  df_selected2 <- NULL
}

d_z_tmp <- zoo(df_selected[[station_id]],df_selected$date)
#d_z_tmp <- zoo(df[[station_id]],df$date)
names(d_z_tmp)<-station_id
#min(d_z_tmp$ID_8)
#max(d_z_tmp$ID_8)
day_start <- "1984-01-01" #PARAM 12 arg 12
day_end <- "2014-12-31" #PARAM 13 arg 13
start_date <- as.Date(day_start)
end_date <- as.Date(day_end)
start_year <- year(start_date)
end_year <- year(end_date)

res_pix <- 1000
#res_pix <- 480
col_mfrow <- 2
row_mfrow <- 1
  
png_filename <-  file.path(out_dir,paste("Figure5a_time_series_profile_",region_name,"_",out_suffix,".png",sep =""))
title_str <- paste("Predicted daily ", station_id," ",var," for the ", start_year,"-",end_year," time period",sep="")
  
png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
#this is the whole time series
plot(d_z_tmp,ylab="tmax in deg C",xlab="Daily time steps",
     main=title_str,cex=3,font=2,
     cex.main=1.5,cex.lab=1.5,font.lab=2,
     lty=3)
if(!is.null(df_selected2)){
  lines(d_z_tmp2,ylab="tmax in deg C",xlab="Daily time steps",
     main=title_str,cex=3,font=2,
     col="red",
     cex.main=1.5,cex.lab=1.5,font.lab=2,
     lty=3)
}
dev.off()

day_start <- "1984-01-01" #PARAM 12 arg 12
day_end <- "1998-12-31" #PARAM 13 arg 13

start_date <- as.Date(day_start)
end_date <- as.Date(day_end)
start_year <- year(start_date)
end_year <- year(end_date)

d_z <- window(d_z_tmp,start=start_date,end=end_date)
if(!is.null(df_selected2)){
  d_z2 <- window(d_z_tmp2,start=start_date,end=end_date)
}

res_pix <- 1000
#res_pix <- 480
col_mfrow <- 2
row_mfrow <- 1
  
png_filename <-  file.path(out_dir,paste("Figure5b_time_series_profile_",region_name,"_",out_suffix,".png",sep =""))
title_str <- paste("Predicted daily ", station_id," ",var," for the ", start_year,"-",end_year," time period",sep="")
  
png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)

plot(d_z,ylab="tmax in deg C",xlab="Daily time steps",
     main=title_str,cex=3,font=2,
     cex.main=1.5,cex.lab=1.5,font.lab=2,
     lty=3)
if(!is.null(df_selected2)){
  lines(d_z2,ylab="tmax in deg C",xlab="Daily time steps",
        main=title_str,cex=3,font=2,
        col="red",
        cex.main=1.5,cex.lab=1.5,font.lab=2,
        lty=3)
}

dev.off()

#### Subset for 5c: zoom in

day_start <- "1991-01-01" #PARAM 12 arg 12
day_end <- "1992-12-31" #PARAM 13 arg 13

start_date <- as.Date(day_start)
end_date <- as.Date(day_end)
start_year <- year(start_date)
end_year <- year(end_date)
d_z <- window(d_z_tmp,start=start_date,end=end_date)
if(!is.null(df_selected2)){
  d_z2 <- window(d_z_tmp2,start=start_date,end=end_date)
}

res_pix <- 1000
#res_pix <- 480
col_mfrow <- 2
row_mfrow <- 1
  
png_filename <-  file.path(out_dir,paste("Figure5c_subset_time_series_profile_",region_name,"_",out_suffix,".png",sep =""))
title_str <- paste("Predicted daily ", station_id," ",var," for the ", start_year,"-",end_year," time period",sep="")
  
png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)

plot(d_z,ylab="tmax in deg C",xlab="Daily time steps",
     main=title_str,cex=3,font=2,
     cex.main=1.5,cex.lab=1.5,font.lab=2,
     lty=3)
if(!is.null(df_selected2)){
  lines(d_z2,ylab="tmax in deg C",xlab="Daily time steps",
        main=title_str,cex=3,font=2,
        col="red",
        cex.main=1.5,cex.lab=1.5,font.lab=2,
        lty=3)
}

dev.off()

####################
###### Now add figures with additional met stations?

#df_selected2 <- data_stations_temp

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
data_stations$date <- as.Date(strptime(data_stations_temp$date_str,"%Y%m%d"))
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

dir_ac_mosaics <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/output_reg4_1999"
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
