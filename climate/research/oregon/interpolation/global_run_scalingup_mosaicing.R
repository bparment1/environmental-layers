##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Different options to explore mosaicing are tested.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 04/14/2015  
#MODIFIED ON: 12/20/2015            
#Version: 5
#PROJECT: Environmental Layers project     
#COMMENTS: analyses run for reg4 1992 for test of mosaicing using 1500x4500km and other tiles
#TODO:
#1) Make this is a script/function callable from the shell/bas
#2) clean up temporary files, it builds currently on the disk
#3) fix output folder for some of output files: create a mosaic output folder if doesn't exist?
#4) create a helper function for inputs/arguments to automate...?? Could also be in the assessment stage

### Before running, the gdal modules and other environment parameters need to be set if on NEX-NASA.
### This can be done by running the following commands:
#
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#MODULEPATH=$MODULEPATH:/nex/modules/files
#module load pythonkits/gdal_1.10.0_python_2.7.3_nex
#
#setfacl -Rmd user:aguzman4:rwx /nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015
#
#reg1   : North America
#reg23" : Europe + Asia
#reg4"  : South America
#reg5   : Africa
#reg6   : Oceania+ South East Asia
#

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
#library(spgwr)                               # GWR method
library(automap)                             # Kriging automatic fitting of variogram using gstat
library(rgeos)                               # Geometric, topologic library of functions
#RPostgreSQL                                 # Interface R and Postgres, not used in this script
library(gridExtra)
#Additional libraries not used in workflow
library(pgirmess)                            # Krusall Wallis test with mulitple options, Kruskalmc {pgirmess}  
library(colorRamps)
library(zoo)
library(xts)

#### FUNCTION USED IN SCRIPT

function_mosaicing <-"global_run_scalingup_mosaicing_function_12192015.R"

in_dir_script <-"/home/parmentier/Data/IPLANT_project/env_layers_scripts" #NCEAS UCSB
#in_dir_script <- "/nobackupp8/bparmen1/env_layers_scripts" #NASA NEX
source(file.path(in_dir_script,function_mosaicing))

load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

create_dir_fun <- function(out_dir,out_suffix){
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    out_dir <- file.path(out_dir,out_name)
  }
  #create if does not exists
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  return(out_dir)
}

############################################
#### Parameters and constants  

#Data is on ATLAS or NASA NEX
#PARAM 1
in_dir <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_12072015" #NCEAS
#in_dir <- "/nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_12072015" #NEX

in_dir_tiles <- file.path(in_dir,"tiles") #this is valid both for Atlas and NEX

y_var_name <- "dailyTmax" #PARAM2
interpolation_method <- c("gam_CAI") #PARAM3
region_name <- "reg4" #PARAM 4 #reg4 South America, Africa reg5,Europe reg2, North America reg1, Asia reg3
mosaicing_method <- c("unweighted","use_edge_weights") #PARAM5
out_suffix <- paste(region_name,"_","run10_1500x4500_global_analyses_pred_1992_12072015",sep="") #PARAM 6
out_suffix_str <- "run10_1500x4500_global_analyses_pred_1992_12072015" #PARAM 7
metric_name <- "rmse" #RMSE, MAE etc. #PARAM 8
pred_mod_name <- "mod1" #PARAM 9
var_pred <- "res_mod1" #used in residuals mapping #PARAM 10

out_dir <- in_dir #PARAM 11
create_out_dir_param <- FALSE #PARAM 12

#if daily mosaics NULL then mosaicas all days of the year #PARAM 13
day_to_mosaic <- c("19920101","19920102","19920103") #,"19920104","19920105") #PARAM9, two dates note in /tiles for now on NEX

#CRS_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84 #CONSTANT1
#CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
#proj_str<- CRS_WGS84 #PARAM 8 #check this parameter
 
file_format <- ".tif" #PARAM 14
NA_value <- -9999 #PARAM 15
NA_flag_val <- NA_value #PARAM 16
     
num_cores <- 6 #PARAM 17                 
region_names <- c("reg23","reg4") #selected region names, ##PARAM 18 
use_autokrige <- F #PARAM 19

###Separate folder for masks by regions, should be listed as just the dir!!... #PARAM 20
#infile_mask <- "/nobackupp8/bparmen1/regions_input_files/r_mask_reg4.tif"
infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_reg4.tif"

#tb_accuracy_name <- file.path(in_dir,paste("tb_diagnostic_v_NA","_",out_suffix_str,".txt",sep=""))
#tb_accuracy_name <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_12072015/tb_diagnostic_v_NA_run10_1500x4500_global_analyses_pred_1992_12072015.txt" #PARAM 21
#data_month_s_name <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_12072015/data_month_s_NAM_run10_1500x4500_global_analyses_pred_1992_12072015.txt" #PARAM 22
#data_day_v_name <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_12072015/data_day_v_NAM_run10_1500x4500_global_analyses_pred_1992_12072015.txt" #PARAM 23
#data_day_s_name <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_12072015/data_day_s_NAM_run10_1500x4500_global_analyses_pred_1992_12072015.txt" ##PARAM 24
#df_tile_processed_name <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_12072015/df_tile_processed_run10_1500x4500_global_analyses_pred_1992_12072015.txt" ##PARAM 25

#in_dir can be on NEX or Atlas
tb_accuracy_name <- file.path(in_dir,"tb_diagnostic_v_NA_run10_1500x4500_global_analyses_pred_1992_12072015.txt") #PARAM 21
data_month_s_name <- file.path(in_dir,"data_month_s_NAM_run10_1500x4500_global_analyses_pred_1992_12072015.txt") #PARAM 22
data_day_v_name <- file.path(in_dir,"data_day_v_NAM_run10_1500x4500_global_analyses_pred_1992_12072015.txt") #PARAM 23
data_day_s_name <- file.path(in_dir,"data_day_s_NAM_run10_1500x4500_global_analyses_pred_1992_12072015.txt") ##PARAM 24
df_tile_processed_name <- file.path(in_dir,"df_tile_processed_run10_1500x4500_global_analyses_pred_1992_12072015.txt") ##PARAM 25

#python script and gdal on NEX NASA:
#mosaic_python <- "/nobackupp6/aguzman4/climateLayers/sharedCode/"
#python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin"
#python script and gdal on Atlas NCEAS
mosaic_python <- "/data/project/layers/commons/NEX_data/sharedCode" #PARAM 26
python_bin <- "/usr/bin" #PARAM 27

algorithm <- "python" #PARAM 28 #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
#algorithm <- "R" #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
 
match_extent <- "FALSE" #PARAM 29 #try without matching!!!

#for residuals...
list_models <- NULL #PARAM 30
#list_models <- paste(var_pred,"~","1",sep=" ") #if null then this is the default...

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

# accuracy table by tiles
tb <- read.table(tb_accuracy_name,sep=",")
# textfiles of stations by month
data_month_s <- read.table(file.path(data_month_s_name),sep=",")
data_day_s <- read.table(file.path(data_day_s_name),sep=",") #daily testing/validation stations by dates and tiles
data_day_v <- read.table(file.path(data_day_v_name),sep=",") #daily training stations by dates and tiles

df_tile_processed <- read.table( df_tile_processed_name,sep=",")

#list all files to mosaic for a list of day
#Take into account multiple region in some cases!!!  
reg_lf_mosaic <- vector("list",length=length(region_names))
for(k in 1:length(region_names)){
  in_dir_tiles_tmp <- file.path(in_dir_tiles, region_names[k])
  reg_lf_mosaic[[k]] <- lapply(1:length(day_to_mosaic),FUN=function(i){list.files(path=file.path(in_dir_tiles_tmp),    
           pattern=paste(".*.",day_to_mosaic[i],".*.tif$",sep=""),full.names=T,recursive=F)})
}

##################### PART 2: produce the mosaic ##################

#This is is assuming a list of file for a region!! 
#this is where the main function for mosaicing region starts!!
#use reg4 to test the code for now, redo later for any regions!!!
k<-2
for(k in 1:length(region_names)){
  region_selected <- region_names[k]
  ##########################
  #### First generate rmse images for each date and tile for the region


  lf_mosaic <- reg_lf_mosaic[[k]] #list of files to mosaic by regions
                                #There a 28 files for reg4, South America

  #######################################
  ################### PART I: Accuracy layers by tiles #############
  #first generate accuracy layers using tiles definitions and output from the accuracy assessment
  
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
  lf <- lf_mosaic

  #Improved by adding multicores option
  num_cores_tmp <- num_cores
  list_param_accuracy_metric_raster <- list(lf,tb,metric_name,pred_mod_name,y_var_name,interpolation_method,
                    days_to_process,num_cores_tmp,NA_flag_val,file_format,out_dir_str,out_suffix_str) 
  names(list_param_accuracy_metric_raster) <- c("lf","tb","metric_name","pred_mod_name","y_var_name","interpolation_method",
                       "days_to_process","num_cores","NA_flag_val","file_format","out_dir_str","out_suffix_str") 
  list_raster_created_obj <- lapply(1:length(day_to_mosaic),FUN=create_accuracy_metric_raster,
                                  list_param=list_param_accuracy_metric_raster)

  #debug(create_accuracy_metric_raster)
  #list_raster_created_obj <- lapply(1:1,FUN=create_accuracy_metric_raster,
  #                                  list_param=list_param_accuracy_metric_raster)
  #raster_created_obj <- create_accuracy_metric_raster(1, list_param_accuracy_metric_raster)

  #Extract list of files for rmse and date 1 (19920101), there should be 28 raster images
  lf_accuracy_raster <- lapply(1:length(list_raster_created_obj),FUN=function(i){unlist(list_raster_created_obj[[i]]$list_raster_name)}) 

  #Plot as quick check
  #r1 <- raster(lf_mosaic[[1]][1]) 
  #plot(r1)

  ####################################
  ### Now create accuracy surfaces from residuals...

  ## Create accuracy surface by kriging

  num_cores_tmp <-num_cores
  lf_day_tiles  <- lf_mosaic #list of raster files by dates
  data_df <- data_day_v # data.frame table/spdf containing stations with residuals and variable
  #df_tile_processed  #tiles processed during assessment usually by region
  #var_pred  #variable being modeled
  #if not list of models is provided generate one
  if(is.null(list_models)){
    list_models <- paste(var_pred,"~","1",sep=" ")
  }

  #use_autokrige #if TRUE use automap/gstat package
  #y_var_name  #"dailyTmax" #PARAM2
  #interpolation_method #c("gam_CAI") #PARAM3, need to select reg!!
  #date_processed #can be a monthly layer
  #num_cores #number of cores used
  #NA_flag_val 
  #file_format 
  out_dir_str <- out_dir
  out_suffix_str <- out_suffix 
  days_to_process <- day_to_mosaic
  df_tile_processed$path_NEX <- as.character(df_tile_processed$path_NEX) 
  df_tile_processed$reg <- basename(dirname(df_tile_processed$path_NEX))

  ##By regions, selected earlier
  #for(k in 1:length(region_names)){
  df_tile_processed_reg <- subset(df_tile_processed,reg==region_selected)#use reg4
  #i<-1 #loop by days/date to process!!
  #test on the first day 
  list_param_create_accuracy_residuals_raster <- list(lf_day_tiles,data_df,df_tile_processed_reg,
                    var_pred,list_models,use_autokrige,y_var_name,interpolation_method,
                    days_to_process,num_cores_tmp,NA_flag_val,file_format,out_dir_str,
                    out_suffix_str) 
  names(list_param_create_accuracy_residuals_raster) <- c("lf_day_tiles","data_df","df_tile_processed_reg",
                    "var_pred","list_models","use_autokrige","y_var_name","interpolation_method",
                    "days_to_process","num_cores_tmp","NA_flag_val","file_format","out_dir_str",
                    "out_suffix_str") 

  list_create_accuracy_residuals_raster_obj <- lapply(1:length(day_to_mosaic),FUN=create_accuracy_residuals_raster,
                                  list_param=list_param_create_accuracy_residuals_raster)

  #undebug(create_accuracy_residuals_raster)
  #list_create_accuracy_residuals_raster_obj <- lapply(1:1,FUN=create_accuracy_residuals_raster,
  #                                list_param=list_param_create_accuracy_residuals_raster)

  #create_accuracy_residuals_raster_obj <- create_accuracy_metric_raster(1, list_param_create_accuracy_residuals_raster_obj)

  #note that three tiles did not produce a residuals surface!!! find out more later, join the output
  #to df_raste_tile to keep track of which one did not work...
  #lf_accuracy_residuals_raster <- as.character(unlist(lapply(1:length(list_create_accuracy_residuals_raster_obj),FUN=function(i,x){unlist(extract_from_list_obj(x[[i]]$list_pred_res_obj,"raster_name"))},x=list_create_accuracy_residuals_raster_obj))) 
  lf_accuracy_residuals_raster <- lapply(1:length(list_create_accuracy_residuals_raster_obj),FUN=function(i,x){as.character(unlist(extract_from_list_obj(x[[i]]$list_pred_res_obj,"raster_name")))},x=list_create_accuracy_residuals_raster_obj)

  #Plot as quick check
  #r1 <- raster(lf_mosaic[[1]][1]) 
  #list_create_accuracy_residuals_raster_obj
  
  ##Run for data_day_s
  ##
  data_df <- data_day_s # data.frame table/spdf containing stations with residuals and variable

  num_cores_tmp <-num_cores
  lf_day_tiles  <- lf_mosaic #list of raster files by dates
  #data_df <- data_day_v # data.frame table/spdf containing stations with residuals and variable
  #df_tile_processed  #tiles processed during assessment usually by region
  #var_pred  #variable being modeled
  #if not list of models is provided generate one
  if(is.null(list_models)){
    list_models <- paste(var_pred,"~","1",sep=" ")
  }

  #use_autokrige #if TRUE use automap/gstat package
  #y_var_name  #"dailyTmax" #PARAM2
  #interpolation_method #c("gam_CAI") #PARAM3, need to select reg!!
  #date_processed #can be a monthly layer
  #num_cores #number of cores used
  #NA_flag_val 
  #file_format 
  out_dir_str <- out_dir
  out_suffix_str <- paste("data_day_s_",out_suffix,sep="") 
  days_to_process <- day_to_mosaic
  df_tile_processed$path_NEX <- as.character(df_tile_processed$path_NEX) 
  df_tile_processed$reg <- basename(dirname(df_tile_processed$path_NEX))

  ##By regions, selected earlier
  #for(k in 1:length(region_names)){
  df_tile_processed_reg <- subset(df_tile_processed,reg==region_selected)#use reg4
  #i<-1 #loop by days/date to process!!
  #test on the first day 
  list_param_create_accuracy_residuals_raster <- list(lf_day_tiles,data_df,df_tile_processed_reg,
                    var_pred,list_models,use_autokrige,y_var_name,interpolation_method,
                    days_to_process,num_cores_tmp,NA_flag_val,file_format,out_dir_str,
                    out_suffix_str) 
  names(list_param_create_accuracy_residuals_raster) <- c("lf_day_tiles","data_df","df_tile_processed_reg",
                    "var_pred","list_models","use_autokrige","y_var_name","interpolation_method",
                    "days_to_process","num_cores_tmp","NA_flag_val","file_format","out_dir_str",
                    "out_suffix_str") 

  list_create_accuracy_residuals_raster_obj <- lapply(1:length(day_to_mosaic),FUN=create_accuracy_residuals_raster,
                                  list_param=list_param_create_accuracy_residuals_raster)

  #undebug(create_accuracy_residuals_raster)
  #list_create_accuracy_residuals_raster_obj <- lapply(1:1,FUN=create_accuracy_residuals_raster,
  #                                list_param=list_param_create_accuracy_residuals_raster)

  #create_accuracy_residuals_raster_obj <- create_accuracy_metric_raster(1, list_param_create_accuracy_residuals_raster_obj)

  #note that three tiles did not produce a residuals surface!!! find out more later, join the output
  #to df_raste_tile to keep track of which one did not work...
  #lf_accuracy_residuals_raster <- as.character(unlist(lapply(1:length(list_create_accuracy_residuals_raster_obj),FUN=function(i,x){unlist(extract_from_list_obj(x[[i]]$list_pred_res_obj,"raster_name"))},x=list_create_accuracy_residuals_raster_obj))) 
  lf_accuracy_residuals_data_s_raster <- lapply(1:length(list_create_accuracy_residuals_raster_obj),FUN=function(i,x){as.character(unlist(extract_from_list_obj(x[[i]]$list_pred_res_obj,"raster_name")))},x=list_create_accuracy_residuals_raster_obj)

  ##took 31 minutes to generate the residuals maps for each tiles (28) for region 4
  
  ######################################################
  #### PART 2: GENETATE MOSAIC FOR LIST OF FILES #####
  #################################
  #### Mosaic tiles for the variable predicted and accuracy metric

  #methods availbable:use_sine_weights,use_edge,use_linear_weights
  #only use edge method for now
  #loop to dates..., make this a function...
  list_mosaic_obj <- vector("list",length=length(day_to_mosaic))
  for(i in 1:length(day_to_mosaic)){
    #
    mosaic_method <- "use_edge_weights" #this is distance from edge
    out_suffix_tmp <- paste(interpolation_method,y_var_name,day_to_mosaic[i],out_suffix,sep="_")
    #debug(mosaicFiles)
    #can also loop through methods!!!
    #python_bin <- "/usr/bin/" #python gdal bin, on Atlas NCEAS
    #python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules/bin" #on NEX
    #gdal_merge_sum_noDataTest.py

    mosaic_edge_obj_prediction <- mosaicFiles(lf_mosaic[[i]],
                                        mosaic_method="use_edge_weights",
                                        num_cores=num_cores,
                                        r_mask_raster_name=infile_mask,
                                        python_bin=python_bin,
                                        mosaic_python=mosaic_python,
                                        algorithm=algorithm,
                                        match_extent=match_extent,
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

    ### produce residuals mosaics
    #for now add data_day_s in the name!!
    mosaic_method <- "use_edge_weights" #this is distance from edge
    out_suffix_tmp <- paste(interpolation_method,"kriged_residuals","data_day_s",day_to_mosaic[i],out_suffix,sep="_")
    #lf_tmp<-list.files(pattern="*kriged_residuals.*.tif",full.names=T)
    lf_tmp <- lf_accuracy_residuals_raster[[i]]
    #lf_accuracy_residuals_raster[[i]]
    #debug(mosaicFiles)
    mosaic_edge_obj_residuals <- mosaicFiles(lf_tmp,
                                        mosaic_method="use_edge_weights",
                                        num_cores=num_cores,
                                        r_mask_raster_name=infile_mask,
                                        python_bin=python_bin,
                                        mosaic_python=mosaic_python,
                                        algorithm=algorithm,
                                        match_extent=match_extent,
                                        df_points=NULL,
                                        NA_flag=NA_flag_val,
                                        file_format=file_format,
                                        out_suffix=out_suffix_tmp,
                                        out_dir=out_dir)
  
    list_mosaic_obj[[i]] <- list(prediction=mosaic_edge_obj_prediction,accuracy=mosaic_edge_obj_accuracy,mosaic_edge_obj_residuals)
  }

  ##End of mosaicing function for region predictions
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


list_param_plot_mosaic <- list(lf_raster_fname=unlist(lf_mean_mosaic[1:2]),
                               screenRast=TRUE,
                               l_dates=day_to_mosaic,
                               out_dir_str=out_dir,
                               out_prefix_str <- "dailyTmax_",
                               out_suffix_str=out_suffix)
debug(plot_screen_raster_val)
plot_screen_raster_val(1,list_param_plot_mosaic)


num_cores <- 2
l_png_files <- mclapply(1:length(unlist(lf_mean_mosaic)[1:2]),FUN=plot_screen_raster_val,
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

