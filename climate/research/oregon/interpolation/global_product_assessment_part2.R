####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Assessment of product part 1: mosaic and accuracy ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This part 2 of the assessment focuses on graphics to explore the spatial patterns of raster times series as figures and movie
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/03/2016  
#MODIFIED ON: 10/10/2016            
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

#COMMIT: generating animation for reg6 (Australia and South East Asia)

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
function_product_assessment_part2_functions <- "global_product_assessment_part2_functions_10102016b.R"
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
in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/assessment"
in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/mosaics/mosaic" #predicted mosaic

region_name <- c("reg6") #param 6, arg 3
out_suffix <- "global_assessment_reg6_10102016"

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
day_end <- "19991231" #PARAM 13 arg 13
#date_start <- day_start
#date_end <- day_end

#infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg4.tif"
#infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg5.tif"
infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg6.tif"

#run_figure_by_year <- TRUE # param 10, arg 7
list_year_predicted <- "1984,2014"
scaling <- 0.01 #was scaled on 100 
#if scaling is null then perform no scaling!!

#df_centroids_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/output_reg5_1999/df_centroids_19990701_reg5_1999.txt"
df_centroids_fname <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/mosaic/output_reg6_1984/df_centroids_19840101_reg6_1984.txt"
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
lf_raster <- NULL #list of raster to consider

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

in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/mosaics/mosaic"
## using predictions
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
#plot(r_stack,zlim=c(-50,50),col=matlab.like(255))
var_name <- "dailyTmax"
debug(plot_and_animate_raster_time_series)

metric_name <- "var_pred" #use RMSE if accuracy
#df_raster <- read.table("df_raster_global_assessment_reg6_10102016.txt",sep=",",header=T)
#plot_figure <- 
plot_and_animate_raster_time_series(lf_raster=lf_raster,
                                    NAvalue=NA_flag_val, 
                                    item_no=13,
                                    region_name=region_name,
                                    var_name=var_name,
                                    metric_name=metric_name,
                                    frame_speed=frame_speed,
                                    animation_format=animation_format,
                                    zlim_val=NULL,
                                    plot_figure=F,
                                    generate_animation=T,
                                    num_cores=num_cores,
                                    out_suffix=out_suffix,
                                    out_dir=out_dir)

## Create function here:

plot_and_animate_raster_time_series <- function(lf_raster, item_no,region_name,var_name,metric_name,NA_flag_val,filenames_figures=NULL,frame_speed=60,animation_format=".gif",zlim_val=NULL,plot_figure=T,generate_animation=T,num_cores=2,out_suffix="",out_dir="."){
  #Function to generate figures and animation for a list of raster
  #
  #
  #INPUTS
  #1) lf_raster
  #2) filenames_figures
  #2) NAvalue
  #3) item_no
  #4) region_name,
  #5) var_name
  #6) metric_name
  #7) frame_speed
  #8) animation_format
  #9) zlim_val
  #10) plot_figure
  #11) generate_animation
  #12) num_cores
  #13) out_suffix
  #14) out_dir
  #OUTPUTS
  #
  #
  

  
  #lf_mosaic_list <- lf_raster
  variable_name <- var_name
  
  if(!is.null(plot_figure)){
    #item_no <- 13
    list_dates_produced <- unlist(mclapply(1:length(lf_raster),
                                           FUN = extract_date,
                                           x = lf_raster,
                                           item_no = item_no,
                                           mc.preschedule = FALSE,
                                           mc.cores = num_cores))

    list_dates_produced_date_val <- as.Date(strptime(list_dates_produced, "%Y%m%d"))
    month_str <- format(list_dates_produced_date_val, "%b") ## Month, char, abbreviated
    year_str <- format(list_dates_produced_date_val, "%Y") ## Year with century
    day_str <- as.numeric(format(list_dates_produced_date_val, "%d")) ## numeric month
    df_raster <- data.frame(lf = basename(lf_raster),
                          date = list_dates_produced_date_val,
                          month_str = month_str,
                          year = year_str,
                          day = day_str,
                          dir = dirname(lf_raster))

    df_raster_fname <- file.path(out_dir, paste0("df_raster_", out_suffix, ".txt"))
    write.table(df_raster,file = df_raster_fname,sep = ",",row.names = F)
    
    ############### PART5: Make raster stack and display maps #############
    #### Extract corresponding raster for given dates and plot

    r_stack <- stack(lf_raster,quick=T)
    l_dates <- list_dates_produced_date_val #[1:11]

    #undebug(plot_raster_mosaic)
    zlim_val <- zlim_val

    ### Now run for the full time series
    #13.26 Western time: start
    #l_dates <- list_dates_produced_date_val
    #r_stack_subset <- r_stack
    #zlim_val <- NULL
    out_suffix_str <- paste0(var_name,"_",metric_name,"_",out_suffix)
    list_param_plot_raster_mosaic <- list(l_dates,r_stack,NA_flag_val,out_dir,
                                          out_suffix_str,region_name,variable_name,zlim_val)
    names(list_param_plot_raster_mosaic) <- c("l_dates","r_mosaiced_scaled","NA_flag_val_mosaic","out_dir",
                                              "out_suffix", "region_name","variable_name","zlim_val")
  
    #lf_mosaic_plot_fig <- mclapply(1:length(l_dates[1:11]),
    #                               FUN = plot_raster_mosaic,
    #                               list_param = list_param_plot_raster_mosaic,
    #                               mc.preschedule = FALSE,
    #                               mc.cores = num_cores)
    
    ##start at 12.29
    ##finished at 15.23 (for reg 6 with 2,991 figures)
    lf_mosaic_plot_fig <- mclapply(1:length(l_dates),
                                   FUN = plot_raster_mosaic,
                                   list_param = list_param_plot_raster_mosaic,
                                   mc.preschedule = FALSE,
                                   mc.cores = num_cores)
    
    if (is.null(zlim_val)) {
      out_suffix_movie <- paste("min_max_", out_suffix_str, sep = "")
    } else{
      zlim_val_str <- paste(zlim_val, sep = "_", collapse = "_")
      out_suffix_movie <- paste(zlim_val_str, "_", out_suffix, sep = "")
    }
    filenames_figures_mosaic <- paste0("list_figures_animation_", out_suffix_movie, ".txt")

    write.table(unlist(lf_mosaic_plot_fig),filenames_figures_mosaic,row.names = F,col.names = F,quote = F)
    
  }
  
  ### Part 2 generate movie
  
  if(generate_animation==TRUE){
    
    out_suffix_str <- paste0(var_name,"_",metric_name,"_",out_suffix)

    if (is.null(zlim_val)) {
      out_suffix_movie <- paste("min_max_", out_suffix, sep = "")
    } else{
      zlim_val_str <- paste(zlim_val, sep = "_", collapse = "_")
      out_suffix_movie <- paste(zlim_val_str, "_", out_suffix, sep = "")
    }

    #already provided as a parameter
    #filenames_figures_mosaic <- paste0("list_figures_animation_", out_suffix_movie, ".txt")
    #write.table(unlist(lf_mosaic_plot_fig),filenames_figures_mosaic,row.names = F,col.names = F,quote = F)

    #now generate movie with imageMagick
    #frame_speed <- 60
    #animation_format <- ".gif"
    #out_suffix_str <- out_suffix
    #started
    #debug(generate_animation_from_figures_fun)
    #lf_mosaic_plot_fig <- read.table(filenames_figure,sep=",")
    
    #out_filename_figure_animation <- generate_animation_from_figures_fun(filenames_figures = unlist(lf_mosaic_plot_fig[1:11]),
    #                                    frame_speed = frame_speed,
    #                                    format_file = animation_format,
    #                                    out_suffix = out_suffix_str,
    #                                    out_dir = out_dir,
    #                                    out_filename_figure_animation = "test2_reg6_animation.gif")

    #started 17.36 Western time on Oct 10 and 18.18
    #        15.58                 oct 11 for 
    out_filename_figure_animation <- generate_animation_from_figures_fun(filenames_figures = filenames_figures_mosaic,
                                        frame_speed = frame_speed,
                                        format_file = animation_format,
                                        out_suffix = out_suffix_movie,
                                        out_dir = out_dir,
                                        out_filename_figure_animation = NULL)
  }
  
  ## prepare object to return
  
  figure_animation_obj <- list(filenames_figures_mosaic,out_filename_figure_animation)
  names(figure_animation_obj) <- c("filenames_figures_mosaic","out_filename_figure_animation")
  return(figure_animation_obj)

}

############ Now accuracy
#### PLOT ACCURACY METRICS: First test ####
##this will be cleaned up later:

in_dir_mosaic_RMSE <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg6/mosaicsRMSE/mosaic"
pattern_str <-"*.tif"
in_dir_mosaic <- in_dir_mosaic_RMSE
lf_raster_rmse <- list.files(path=in_dir_mosaic,pattern=pattern_str,recursive=F,full.names=T)
lf_raster <- lf_raster_rmse
r_stack <- stack(lf_raster,quick=T) #this is very fast now with the quick option!
#save(r_mosaic,file="r_mosaic.RData")

NAvalue(r_stack)
plot(r_stack,y=6,zlim=c(0,8000)) #this is not rescaled

lf_mosaic_list <- lf_raster
list_dates_produced_RMSE <-  mclapply(1:2,
                                 FUN=extract_date,
                                 x=lf_mosaic_list,
                                 item_no=15,
                                 mc.preschedule=FALSE,
                                 mc.cores = 2)  
item_no <-15
list_dates_produced_RMSE <- unlist(mclapply(1:length(lf_raster),
                                       FUN=extract_date,
                                       x=lf_raster,
                                       item_no=item_no,
                                       mc.preschedule=FALSE,
                                       mc.cores = num_cores))                         

list_dates_produced_date_val <- as.Date(strptime(list_dates_produced_RMSE,"%Y%m%d"))
month_str <- format(list_dates_produced_date_val, "%b") ## Month, char, abbreviated
year_str <- format(list_dates_produced_date_val, "%Y") ## Year with century
day_str <- as.numeric(format(list_dates_produced_date_val, "%d")) ## numeric month

df_raster_rmse <- data.frame(lf=basename(lf_raster),
                          date=list_dates_produced_date_val,
                          month_str=month_str,
                          year=year_str,
                          day=day_str,
                          dir=dirname(lf_mosaic_list))

df_raster_fname <- file.path(out_dir,paste0("df_raster_rmse",out_suffix,".txt"))
write.table(df_raster,file= df_raster_fname,sep=",",row.names = F) 


r_stack_subset <- subset(r_stack,1:11)
l_dates <- list_dates_produced_date_val[1:11]

#undebug(plot_raster_mosaic)
out_suffix_str <- paste0("rmse_",out_suffix)

zlim_val <- NULL
##Need to add title option!!
list_param_plot_raster_mosaic <- list(l_dates,r_stack_subset,NA_flag_val,out_dir,out_suffix_str,
                                      region_name,variable_name, zlim_val)
names(list_param_plot_raster_mosaic) <- c("l_dates","r_mosaiced_scaled","NA_flag_val_mosaic","out_dir","out_suffix",
                                          "region_name","variable_name","zlim_val")
lf_mosaic_plot_fig <- lapply(1:2,
                               FUN=plot_raster_mosaic,
                               list_param=list_param_plot_raster_mosaic)         

### Now run for the full time series
#13.26 Western time: start
l_dates <- list_dates_produced_date_val
r_stack_subset <- r_stack
zlim_val <- NULL
list_param_plot_raster_mosaic <- list(l_dates,r_stack_subset,NA_flag_val,out_dir,out_suffix,
                                      region_name,variable_name, zlim_val)
names(list_param_plot_raster_mosaic) <- c("l_dates","r_mosaiced_scaled","NA_flag_val_mosaic","out_dir","out_suffix",
                                          "region_name","variable_name","zlim_val")
#started at 20.16 on 10/10/2016, finished 22.39
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
#r_stack_subset <- subset(r_stack,1:11)
#l_dates <- list_dates_produced_date_val[1:11]

filenames_figures_mosaic_test <- "list_figures_animation_test_reg6.txt"

write.table(unlist(lf_mosaic_plot_fig[1:11]),filenames_figures_mosaic_test,row.names = F,col.names = F,quote = F)

filenames_figures_mosaic <- paste0("list_figures_animation_",out_suffix_movie,".txt")

write.table(unlist(lf_mosaic_plot_fig),filenames_figures_mosaic,row.names = F,col.names = F,quote = F)

#now generate movie with imageMagick
frame_speed <- 60
animation_format <- ".gif"
out_suffix_str <- out_suffix
out_suffix_movie <- paste0("rmse_",out_suffix_movie)
#started: 

generate_animation_from_figures_fun(filenames_figures= unlist(lf_mosaic_plot_fig[1:11]),
                                    frame_speed=frame_speed,
                                    format_file=animation_format,
                                    out_suffix=out_suffix_str,
                                    out_dir=out_dir,
                                    out_filename_figure_animation="test2_reg6_animation.gif")
#started 6.43 Western time on Oct 11 and 7.19
generate_animation_from_figures_fun(filenames_figures= filenames_figures_mosaic,
                                    frame_speed=frame_speed,
                                    format_file=animation_format,
                                    out_suffix=out_suffix_movie,
                                    out_dir=out_dir,
                                    out_filename_figure_animation=NULL)

############################ END OF SCRIPT ##################################
