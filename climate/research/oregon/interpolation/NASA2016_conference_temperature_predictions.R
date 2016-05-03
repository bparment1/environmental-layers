####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  NASA 2016 Meeting: biodiversity and ecological forecasting ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Combining tables and figures for individual runs for years and tiles.
#Figures and data for the AAG conference are also produced.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 05/01/2016  
#MODIFIED ON: 05/02/2016            
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

###### Function used in the script #######
  
#script_path <- "/nobackupp8/bparmen1/env_layers_scripts" #path to script
script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts" #path to script

#Mosaic related
#script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts"
function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_04232016.R" #PARAM12
function_mosaicing <-"global_run_scalingup_mosaicing_05012016.R"
source(file.path(script_path,function_mosaicing)) #source all functions used in this script 
source(file.path(script_path,function_mosaicing_functions)) #source all functions used in this script 

#Assessment
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

pre_process_raster_mosaic_fun <- function(i,list_param){
  
  
  ## Extract parameters
  
  lf <- list_param$lf
  python_bin <- list_param$python_bin
  infile_mask <- list_param$infile_mask
  scaling <- list_param$scaling
  mask_pred <- list_param$mask_pred
  NA_flag_val <- list_param$NA_flag_val
  out_suffix <- list_param$out_suffix
  out_dir <- list_param$out_dir
  
  raster_name_in <- lf[i]
  
  #Step 1: match extent and resolution
  
  lf_files <- c(raster_name_in) #match to mask
  rast_ref <- infile_mask
  ##Maching resolution is probably only necessary for the r mosaic function
  #Modify later to take into account option R or python...
  list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix,out_dir)
  names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")
  r_pred_matched <- raster_match(1,list_param_raster_match)
  raster_name_in <- c(r_pred_matched)

  #Step 2: mask
  if(mask_pred==TRUE){
    r_mask <- raster(infile_mask)
    extension_str <- extension(raster_name_in)
    raster_name_tmp <- gsub(extension_str,"",basename(raster_name_in))
    raster_name <- file.path(out_dir,paste(raster_name_tmp,"_masked.tif",sep = ""))
    r_pred <- mask(raster(raster_name_in),r_mask,filename = raster_name,overwrite = TRUE)
  }
  
  NAvalue(r_pred) <- NA_flag_val
  #r_pred <- setMinMax(r_pred)

  #Step 3: remove scaling factor
  raster_name_in <- filename(r_pred)
  extension_str <- extension(raster_name_in)
  raster_name_tmp <- gsub(extension_str,"",basename(filename(r_pred)))
  raster_name_out <- file.path(out_dir,paste(raster_name_tmp,"_rescaled.tif",sep = ""))
  #r_pred <- overlay(r_pred, fun=function(x){x*1/scaling},filename=raster_name,overwrite=T)

  #raster_name_in <- "comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19990102_reg4_1999_m_gam_CAI_dailyTmax_19990102_reg4_1999_m__meeting_NASA_reg4_04292016_masked.tif"
  python_cmd <- file.path(python_bin,"gdal_calc.py")
  cmd_str3 <- paste(python_cmd, 
                     paste("-A ", raster_name_in,sep=""),
                     paste("--outfile=",raster_name_out,sep=""),
                     #paste("--type=","Int32",sep=""),
                     "--co='COMPRESS=LZW'",
                     paste("--NoDataValue=",NA_flag_val,sep=""),
                     paste("--calc='(A)*",scaling,"'",sep=""),
                     "--overwrite",sep=" ") #division by zero is problematic...
  #system(cmd_str3)
  system(cmd_str3)
  #NAvalue(r_pred) <- NA_flag_val
  #r_pred <- 
  r_pred <- setMinMax(raster(raster_name_out))
  
  return(raster_name_out)
}

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

day_start <- "19840101" #PARAM 12 arg 12
day_end <- "20021231" #PARAM 13 arg 13

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

NA_flag_val_mosaic <- -3399999901438340239948148078125514752.000

list_param_plot_raster_mosaic <- list(l_dates,r_mosaic_scaled,NA_flag_val_mosaic,out_dir,out_suffix,
                                      region_name,variable_name)
names(list_param_plot_raster_mosaic) <- c("l_dates","r_mosaic_scaled","NA_flag_val_mosaic","out_dir","out_suffix",
                                          "region_name","variable_name")

lf_mosaic_plot_fig <- mclapply(1:length(lf_mosaic_scaled),FUN=plot_raster_mosaic,list_param=list_param_plot_raster_mosaic,mc.preschedule=FALSE,mc.cores = num_cores)                         

plot_raster_mosaic <- function(i,list_param){
  #Function to plot mosaic for poster
  #
  l_dates <- list_param$l_dates
  r_mosaiced_scaled <- list_param$r_mosaiced_scaled
  NA_flag_val <- list_param$NA_flag_val
  out_dir <- list_param$out_dir
  out_suffix <- list_param$out_suffix
  region_name <- list_param$region_name
  variable_name <- list_param$variable_name

#for (i in 1:length(nlayers(r_mosaic_scaled))){
  
  date_proc <- l_dates[i]
  r_pred <- subset(r_mosaic_scaled,i)
  NAvalue(r_pred)<- NA_flag_val 
 
  date_proc <- l_dates[i]
  date_val <- as.Date(strptime(date_proc,"%Y%m%d"))
  #month_name <- month.name(date_val)
  month_str <- format(date_val, "%b") ## Month, char, abbreviated
  year_str <- format(date_val, "%Y") ## Year with century
  day_str <- as.numeric(format(date_val, "%d")) ## numeric month
  date_str <- paste(month_str," ",day_str,", ",year_str,sep="")
  
  res_pix <- 1200
  #res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1
  
  png_filename <-  file.path(out_dir,paste("Figure4_clim_mosaics_day_","_",date_proc,"_",region_name,"_",out_suffix,".png",sep =""))
  title_str <-  paste("Predicted ",variable_name, " on ",date_str , " ", sep = "")
  
  png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
  plot(r_pred,main =title_str,cex.main =1.5,col=matlab.like(255),zlim=c(-50,50),
       legend.shrink=0.8,legend.width=0.8)
       #axis.args = list(cex.axis = 1.6), #control size of legend z
       #legend.args=list(text='dNBR', side=4, line=2.5, cex=2.2))
       #legend.args=list(text='dNBR', side=4, line=2.49, cex=1.6))
  dev.off()

  return(png_filename)
}

############### PART2: temporal profile #############
#### Extract time series
###
#-65,-22

df_points <- read.table(df_points_extracted_fname,sep=",") 
df_centroids <- read.table(df_centroids_fname,sep=",")

coordinates(df_centroids)<- c("long","lat")
proj4string(df_centroids) <- CRS_locs_WGS84

lf_mosaic_list <- list.files(path=in_dir_mosaic,pattern="*.tif")
#r_mosaic_ts <- stack(lf_mosaic_list)
#df_centroids <- extract(df_centroids,r_mosaic_ts)

df_points$files <- lf_mosaic_list

extract_date <- function(i,x){
  y <- unlist(strsplit(x[[i]],"_"))
  date_str <- y[length(y)-2]
}
#debug(extract_date)
#test <- (extract_date(1,lf_mosaic_list))

list_dates_produced <- unlist(mclapply(1:length(lf_mosaic_list),FUN=extract_date,mc.preschedule=FALSE,mc.cores = num_cores))                         
#list_dates_produced <- mclapply(1:11,FUN=extract_date,mc.preschedule=FALSE,x=lf_mosaic_list,mc.cores = num_cores)                         


list_dates_produced_date_val <- as.Date(strptime(list_dates_produced,"%Y%m%d"))
month_str <- format(list_dates_produced_date_val, "%b") ## Month, char, abbreviated
year_str <- format(list_dates_produced_date_val, "%Y") ## Year with century
day_str <- as.numeric(format(list_dates_produced_date_val, "%d")) ## numeric month

df_points$date <- list_dates_produced_date_val
df_points$month <- month_str
df_points$year <- year_str
df_points$day <- day_str

#data_pixel <- data_df[id_selected,]
#data_pixel$rainfall <- as.numeric(data_pixel$rainfall)
#d_z_tmp <-zoo(data_pixel$rainfall,as.Date(data_pixel$date))
#names(d_z_tmp)<- "rainfall"
#data_pixel <- as.data.frame(data_pixel)
#d_z_tmp2 <- zoo(data_pixel[[var_name]],as.Date(data_pixel$date))
    
#start_date <-input$dates[1]
#end_date <-input$dates[2]
    
#d_z <- window(d_z_tmp,start=start_date,end=end_date)
#d_z2 <- window(d_z_tmp2,start=start_date,end=end_date)
#df2 <- as.data.frame(d_z2)
#names(df2)<- var_name


############################ END OF SCRIPT ##################################