####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Assessment of product part 1: mosaic and accuracy ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This part 2 of the assessment focuses on graphics to explore the spatial patterns of raster times series as figures and movie
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/03/2016  
#MODIFIED ON: 10/03/2016            
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
function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_10032016b.R"
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

out_region_name<-""
list_models<-c("y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)") #param 4

#reg1 (North Am), reg2(Europe),reg3(Asia), reg4 (South Am), reg5 (Africa), reg6 (Australia-Asia)
#master directory containing the definition of tile size and tiles predicted
#in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/assessment"
#in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg5/mosaic/mosaic"
in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/assessment"
in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg1/mosaics/mosaic"

region_name <- c("reg1") #param 6, arg 3
out_suffix <- "_global_assessment_reg1_10032016"

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

#setwd(out_dir)

###########  ####################


#https://www.r-bloggers.com/animated-plots-with-r/


if(is.null(lf_raster)){
  
  #pattern_str <- ".*.tif"
  pattern_str <-"*.tif"
  lf_raster <- list.files(path=in_dir_mosaic,pattern=pattern_str,recursive=F,full.names=T)
  r_stack <- stack(lf_raster,quick=T) #this is very fast now with the quick option!
  #save(r_mosaic,file="r_mosaic.RData")
    
}else{
  r_stack <- stack(lf_raster,quick=T) #this is very fast now with the quick option!
}

NAvalue(r_stack)
plot(r_stack,y=6,zlim=c(-10000,10000)) #this is not rescaled
plot(r_stack,zlim=c(-50,50),col=matlab.like(255))

#plot(r_mosaic_scaled,y=6,zlim=c(-50,50))
#plot(r_mosaic_scaled,zlim=c(-50,50),col=matlab.like(255))

#debug(extract_date)
#test <- extract_date(6431,lf_mosaic_list,12) #extract item number 12 from the name of files to get the data
#list_dates_produced <- unlist(mclapply(1:length(lf_raster),FUN=extract_date,x=lf_raster,item_no=13,mc.preschedule=FALSE,mc.cores = num_cores))                         
lf_mosaic_list <- lf_raster
list_dates_produced <-  mclapply(1:2,
                                 FUN=extract_date,
                                 x=lf_mosaic_list,
                                 item_no=13,
                                 mc.preschedule=FALSE,
                                 mc.cores = 2)  
item_no <-13
list_dates_produced <- unlist(mclapply(1:length(lf_raster),
                                       FUN=extract_date,
                                       x=lf_raster,
                                       item_no=item_no,
                                       mc.preschedule=FALSE,
                                       mc.cores = num_cores))                         

list_dates_produced_date_val <- as.Date(strptime(list_dates_produced,"%Y%m%d"))
month_str <- format(list_dates_produced_date_val, "%b") ## Month, char, abbreviated
year_str <- format(list_dates_produced_date_val, "%Y") ## Year with century
day_str <- as.numeric(format(list_dates_produced_date_val, "%d")) ## numeric month

df_raster <- data.frame(lf=basename(lf_raster),
                          date=list_dates_produced_date_val,
                          month_str=month_str,
                          year=year_str,
                          day=day_str,
                          dir=dirname(lf_mosaic_list))

df_raster_fname <- file.path(out_dir,paste0("df_raster_",out_suffix,".txt"))
write.table(df_raster,file= df_raster_fname,sep=",",row.names = F) 

############### PART5: Make raster stack and display maps #############
#### Extract corresponding raster for given dates and plot stations used

## TODO: make movies from time series in png

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

function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_10032016b.R"
source(file.path(script_path,function_product_assessment_part2_functions)) #source all functions used in this script 

#NA_flag_val_mosaic <- -3399999901438340239948148078125514752.000
r_stack_subset <- subset(r_stack,1:11)
l_dates <- list_dates_produced_date_val[1:11]

#undebug(plot_raster_mosaic)
zlim_val <- NULL
list_param_plot_raster_mosaic <- list(l_dates,r_stack_subset,NA_flag_val,out_dir,out_suffix,
                                      region_name,variable_name, zlim_val)
names(list_param_plot_raster_mosaic) <- c("l_dates","r_mosaiced_scaled","NA_flag_val_mosaic","out_dir","out_suffix",
                                          "region_name","variable_name","zlim_val")

lf_mosaic_plot_fig <- lapply(1:2,
                               FUN=plot_raster_mosaic,
                               list_param=list_param_plot_raster_mosaic)         

lf_mosaic_plot_fig <- mclapply(1:length(l_dates),
                               FUN=plot_raster_mosaic,
                               list_param=list_param_plot_raster_mosaic,
                               mc.preschedule=FALSE,
                               mc.cores = num_cores)    
if(is.null(zlim_val)){
  out_suffix_movie <- paste("min_max_",out_suffix,sep="")
}else{
  zlim_val_str <- paste(zlim_val,sep="_",collapse="_")
  out_suffix_movie <- paste(zlim_val_str,"_",out_suffix,sep="")
}

#filename_figures_mosaic <- file.path(out_dir,"mosaic_plot_fig.txt")
filename_figures_mosaic <- file.path(out_dir,paste("list_mosaic_plot_fig",out_suffix_movie,".txt",sep=""))
out_figure_movies <- file.path(out_dir,paste("mosaic_movie_",out_suffix_movie,".gif",sep=""))

write.table(unlist(lf_mosaic_plot_fig),filename_figures_mosaic,row.names = F,col.names = F,quote = F)
#now generate movie with imageMagick

#-delay 20
delay_option <- 60
cmd_str <- paste("convert",paste("-delay",delay_option),paste0("@",filename_figures_mosaic),out_figure_movies)
#convert @myimages.txt mymovie.gif

system(cmd_str)

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
