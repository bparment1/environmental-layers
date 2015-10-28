##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Different options to explore mosaicing are tested.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 04/14/2015  
#MODIFIED ON: 10/28/2015            
#Version: 5
#PROJECT: Environmental Layers project     
#COMMENTS: analyses run for reg4 1992 for test of mosaicing using 1500x4500km and other tiles
#TODO:
#1) Make this is a script/function callable from the shell/bash
#2) generalize to run dates and region fast (use python mosaic Alberto code)
#3) clean up temporary files, it builds currently on the disk
#4) fix output folder for some of output files
#

### Before running, the gdal modules and other environment parameters need to be set if on NEX-NASA.
### This can be done by running the following commands:
#
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#MODULEPATH=$MODULEPATH:/nex/modules/files
#module load pythonkits/gdal_1.10.0_python_2.7.3_nex
#
#
#################################################################################################

### Loading R library and packages        
#library used in the workflow production:
library(gtools)                              # loading some useful tools 
library(mgcv)                                # GAM package by Simon Wood
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
#library(gstat)                               # Kriging and co-kriging by Pebesma et al.
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
#library(spgwr)                               # GWR method
#library(automap)                             # Kriging automatic fitting of variogram using gstat
#library(rgeos)                               # Geometric, topologic library of functions
#RPostgreSQL                                 # Interface R and Postgres, not used in this script
library(gridExtra)
#Additional libraries not used in workflow
library(pgirmess)                            # Krusall Wallis test with mulitple options, Kruskalmc {pgirmess}  
library(colorRamps)
library(zoo)
library(xts)

#### FUNCTION USED IN SCRIPT

function_mosaicing <-"global_run_scalingup_mosaicing_function_10282015.R"

#in_dir_script <-"/home/parmentier/Data/IPLANT_project/env_layers_scripts" #NCEAS UCSB
in_dir_script <- "/nobackupp8/bparmen1/env_layers_scripts" #NASA NEX
source(file.path(in_dir_script,function_mosaicing))

############################################
#### Parameters and constants  

#Data is on ATLAS: reg4 (South America)

#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test" #PARAM1
#in_dir <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_10052015" #PARAM4
in_dir <- "/nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015" #NEX

in_dir_tiles <- file.path(in_dir,"tiles")
#in_dir_tiles <- "/nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015/tiles" #North America
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg2" #Europe
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg4" #South America
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg5" #Africa

y_var_name <- "dailyTmax" #PARAM2
interpolation_method <- c("gam_CAI") #PARAM3
region_name <- "reg4" #PARAM 4 #reg4 South America, Africa reg5,Europe reg2, North America reg1, Asia reg3
mosaicing_method <- c("unweighted","use_edge_weights") #PARAM5
out_suffix <- paste(region_name,"_","run10_1500x4500_global_analyses_pred_1992_10052015",sep="") #PARAM 6
out_suffix_str <- "run10_1500x4500_global_analyses_pred_1992_10052015"
metric_name <- "rmse" #RMSE, MAE etc.
pred_mod_name <- "mod1"

#PARAM3
out_dir <- in_dir #PARAM 7
create_out_dir_param <- FALSE #PARAM 8

#if daily mosaics NULL then mosaicas all days of the year
day_to_mosaic <- c("19920101","19920102","19920103") #,"19920104","19920105") #PARAM9, two dates note in /tiles for now on NEX

#CRS_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84 #CONSTANT1
#CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
#proj_str<- CRS_WGS84 #PARAM 8 #check this parameter
 
file_format <- ".tif" #PARAM 10
NA_value <- -9999 #PARAM 11
NA_flag_val <- NA_value
     
num_cores <- 6 #PARAM 12                  

#infile_mask <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_10052015/r_mask_reg4.tif"
infile_mask <- "/nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015/r_mask_reg4.tif"

#tb_accuracy_name <- file.path(in_dir,paste("tb_diagnostic_v_NA","_",out_suffix_str,".txt",sep=""))
#tb_accuracy_name <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_10052015/tb_diagnostic_v_NA_run10_1500x4500_global_analyses_pred_1992_10052015.txt"
tb_accuracy_name <- "/nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015/tb_diagnostic_v_NA_run10_1500x4500_global_analyses_pred_1992_10052015.txt"

#python script and gdal on NEX NASA:
mosaic_python <- "/nobackupp6/aguzman4/climateLayers/sharedCode/"
python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin"
#python script and gdal on Atlas NCEAS
#mosaic_python <- "/data/project/layers/commons/NEX_data/sharedCode"
#python_bin <- "/usr/bin"

algorithm <- "python" #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
#algorithm <- "R" #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
 
########################## START SCRIPT ##############################


####### PART 1: Read in data and process data ########

#in_dir <- file.path(in_dir,region_name)
out_dir <- in_dir
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

setwd(out_dir)


tb <- read.table(file=tb_accuracy_name,sep=",")

#list all files to mosaic for a list of day
lf_mosaic <- lapply(1:length(day_to_mosaic),FUN=function(i){list.files(path=file.path(in_dir_tiles),    
           pattern=paste(".*.",day_to_mosaic[i],".*.tif$",sep=""),full.names=T)})
       
##################### PART 2: produce the mosaic ##################

##########################
#### First generate rmse images for each date and tile for the region

#lf <- lf_mosaic1 #list of files to mosaic for date 1, there a 28 files for reg4, South America
lf <- lf_mosaic #list of files to mosaic for date 1, there a 28 files for reg4, South America

#tb <- list_param$tb #fitting or validation table with all days
#metric_name <- "rmse" #RMSE, MAE etc.
#pred_mod_name <- "mod1"
#y_var_name 
#interpolation_method #c("gam_CAI") #PARAM3
days_to_process <- day_to_mosaic
#NA_flag_val <- list_param$NA_flag_val
#file_format <- list_param$file_format
out_dir_str <- out_dir
out_suffix_str <- out_suffix

#Improve by adding multicores option
list_param_accuracy_metric_raster <- list(lf,tb,metric_name,pred_mod_name,y_var_name,interpolation_method,
                    days_to_process,NA_flag_val,file_format,out_dir_str,out_suffix_str) 
names(list_param_accuracy_metric_raster) <- c("lf","tb","metric_name","pred_mod_name","y_var_name","interpolation_method",
                       "days_to_process","NA_flag_val","file_format","out_dir_str","out_suffix_str") 
list_raster_created_obj <- lapply(1:length(day_to_mosaic),FUN=create_accuracy_metric_raster,
                                  list_param=list_param_accuracy_metric_raster)
#debug(create_accuracy_metric_raster)
#list_raster_created_obj <- lapply(1:1,FUN=create_accuracy_metric_raster,
#                                  list_param=list_param_accuracy_metric_raster)

#raster_created_obj <- create_accuracy_metric_raster(1, list_param_accuracy_metric_raster)

#Extract list of files for rmse and date 1 (19920101), there should be 28 raster images
lf_accuracy_raster <- lapply(1:length(list_raster_created_obj),FUN=function(i){unlist(list_raster_created_obj[[i]]$list_raster_name)}) 

#Plot as quick check
r1 <- raster(lf_mosaic[[1]][1]) 
#r2 <- raster(lf_mosaic2[2]) 

#plot(r1)
#plot(r2)
#r1_ac <- raster(lf_accuracy_raster[[1]][1]) 
#plot(r1_ac)

#################################
#### Second mosaic tiles for the variable and accuracy metrid

#methods availbable:use_sine_weights,use_edge,use_linear_weights
#only use edge method for now
#loop to dates..., make this a function...
list_mosaic_obj <- vector("list",length=length(day_to_mosaic))
for(i in 1:length(day_to_mosaic)){
  
  mosaic_method <- "use_edge_weights" #this is distance from edge
  out_suffix_tmp <- paste(interpolation_method,y_var_name,day_to_mosaic[i],out_suffix,sep="_")
  debug(mosaicFiles)
  #can also loop through methods!!!
  #python_bin <- "/usr/bin/" #python gdal bin, on Atlas NCEAS
  #python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules/bin" #on NEX
  mosaic_edge_obj_prediction <- mosaicFiles(lf_mosaic[[i]],
                                        mosaic_method="use_edge_weights",
                                        num_cores=num_cores,
                                        r_mask_raster_name=infile_mask,
                                        python_bin=python_bin,
                                        mosaic_python=mosaic_python,
                                        algorithm=algorithm,
                                        df_points=NULL,
                                        NA_flag=NA_flag_val,
                                        file_format=file_format,
                                        out_suffix=out_suffix_tmp,
                                        out_dir=out_dir)
  
  mosaic_method <- "use_edge_weights" #this is distance from edge
  out_suffix_tmp <- paste(interpolation_method,metric_name,day_to_mosaic[i],out_suffix,sep="_")

  #debug(mosaicFiles)
  #can also loop through methods!!!
  mosaic_edge_obj_accuracy <- mosaicFiles(lf_accuracy_raster[[i]],
                                        mosaic_method="use_edge_weights",
                                        num_cores=num_cores,
                                        r_mask_raster_name=infile_mask,
                                        python_bin=python_bin,
                                        mosaic_python=mosaic_python,
                                        algorithm=algorithm,
                                        df_points=NULL,
                                        NA_flag=NA_flag_val,
                                        file_format=file_format,
                                        out_suffix=out_suffix_tmp,
                                        out_dir=out_dir)
  
  list_mosaic_obj[[i]] <- list(prediction=mosaic_edge_obj_prediction,accuracy=mosaic_edge_obj_accuracy)
}

#####################
###### PART 2: Analysis and figures for the outputs of mosaic function #####

#### compute and aspect and slope with figures
#list_lf_mosaic_obj <- vector("list",length(day_to_mosaic))
#lf_mean_mosaic <- vector("list",length(mosaicing_method))#2methods only
#l_method_mosaic <- vector("list",length(mosaicing_method))
#list_out_suffix <- vector("list",length(mosaicing_method))

#for(i in 1:length(day_to_mosaic)){
#  list_lf_mosaic_obj[[i]] <- list.files(path=out_dir,pattern=paste("*",day_to_mosaic[i],
#                                                                   "_.*.RData",sep=""))
#  lf_mean_mosaic[[i]] <- unlist(lapply(list_lf_mosaic_obj[[i]],function(x){load_obj(x)[["mean_mosaic"]]}))
#  l_method_mosaic[[i]] <- paste(unlist(lapply(list_lf_mosaic_obj[[i]],function(x){load_obj(x)[["method"]]})),day_to_mosaic[i],sep="_")
#  list_out_suffix[[i]] <- unlist(paste(l_method_mosaic[[i]],day_to_mosaic[[i]],out_suffix,sep="_"))
#}


#list_param_plot_mosaic <- list(lf_mosaic=unlist(lf_mean_mosaic),
#                               method=unlist(l_method_mosaic),
#                               out_suffix=unlist(list_out_suffix))

#plot and produce png movie...
#plot_mosaic(1,list_param=list_param_plot_mosaic)
num_cores <- 1
l_png_files <- mclapply(1:length(unlist(lf_mean_mosaic)),FUN=plot_mosaic,
                        list_param= list_param_plot_mosaic,
                        mc.preschedule=FALSE,mc.cores = num_cores)

lf_plot<- list.files(pattern="r_m_use.*.mask.*.tif$")
lf_mean_mosaic <- lf_plot


list_param_plot_mosaic <- list(lf_raster_fname=unlist(lf_mean_mosaic[1:3]),
                               screenRast=FALSE,
                               l_dates=day_to_mosaic,
                               out_dir_str=out_dir,
                               out_prefix_str <- "dailyTmax_",
                               out_suffix_str=out_suffix)
#plot_screen_raster_val(3,list_param_plot_mosaic)
#debug(plot_screen_raster_val)

num_cores <- 3
l_png_files <- mclapply(1:length(unlist(lf_mean_mosaic)[1:3]),FUN=plot_screen_raster_val,
                        list_param= list_param_plot_mosaic,
                        mc.preschedule=FALSE,mc.cores = num_cores)


list_param_plot_mosaic <- list(lf_raster_fname=unlist(lf_mean_mosaic[4:6]),
                               screenRast=FALSE,
                               l_dates=day_to_mosaic,
                               out_dir_str=out_dir,
                               out_prefix_str <- paste("rmse_",out_suffix,sep=""),
                               out_suffix_str=paste("rmse_",out_suffix,sep=""))
num_cores <- 3
l_png_files_rmse <- mclapply(1:length(unlist(lf_mean_mosaic)[4:6]),FUN=plot_screen_raster_val,
                        list_param= list_param_plot_mosaic,
                        mc.preschedule=FALSE,mc.cores = num_cores)

###############

##################### END OF SCRIPT ######################

