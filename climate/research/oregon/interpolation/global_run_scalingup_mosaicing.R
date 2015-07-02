##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Different options to explore mosaicing are tested.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 04/14/2015  
#MODIFIED ON: 06/29/2015            
#Version: 5
#PROJECT: Environmental Layers project     
#COMMENTS: analyses run for reg5 for test of mosaicing using 1500x4500km and other tiles
#TODO:
#1) Make this is a script/function callable from the shell/bash
#2) generalize to run dates and region fast

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

#### FUNCTION USED IN SCRIPT

function_mosaicing <-"global_run_scalingup_mosaicing_function_07012015.R"

in_dir_script <-"/home/parmentier/Data/IPLANT_project/env_layers_scripts"
source(file.path(in_dir_script,function_mosaicing))

############################################
#### Parameters and constants  

#Data is on ATLAS: reg4 (South America)

in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test" #PARAM1
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg1" #North America
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg2" #Europe
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg4" #South America
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg5" #Africa

y_var_name <- "dailyTmax" #PARAM2
interpolation_method <- c("gam_CAI") #PARAM3
region_name <- "reg23" #PARAM 4 #reg4 South America, Africa reg5,Europe reg2, North America reg1, Asia reg3
mosaicing_method <- c("unweighted","use_edge_weights") #PARAM5
out_suffix <- paste(region_name,"_","mosaic_run10_1500x4500_global_analyses_06292015",sep="") #PARAM 6
#PARAM3
out_dir <- in_dir #PARAM 7
create_out_dir_param <- TRUE #PARAM 8

#if daily mosaics NULL then mosaicas all days of the year
day_to_mosaic <- c("20100831",
                   "20100901") #PARAM 9
  
file_format <- ".tif" #PARAM 10
NA_value <- -9999 #PARAM 11
NA_flag_val <- NA_value
     
num_cores <- 11 #PARAM 12                  

########################## START SCRIPT ##############################


####### PART 1: Read in data and process data ########

in_dir <- file.path(in_dir,region_name)
out_dir <- in_dir
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

setwd(out_dir)

lf_mosaic1 <-list.files(path=file.path(in_dir),    
           pattern=paste(".*.",day_to_mosaic[1],".*.tif$",sep=""),full.names=T) #choosing date 2...20100901

lf_mosaic2 <-list.files(path=file.path(in_dir),    
           pattern=paste(".*.",day_to_mosaic[2],".*.tif$",sep=""),full.names=T) #choosing date 2...20100901
#lf_mosaic <- lf_mosaic[1:20]
r1 <- raster(lf_mosaic1[1]) 
r2 <- raster(lf_mosaic2[2]) 

plot(r1)
plot(r2)

#methods availbable:use_sine_weights,use_edge,use_linear_weights
#only use edge method for now
#loop to dates...
list_mosaic_obj <- vector("list",length=length(day_to_mosaic))
for(i in 1:length(day_to_mosaic)){
  
  mosaic_method <- "use_edge_weights"
  out_suffix_str <- paste(day_to_mosaic[i],out_suffix,sep="_")
  #undebug(mosaicFiles)
  #can also loop through methods!!!
  mosaic_edge_obj <- mosaicFiles(lf_mosaic1,mosaic_method="use_edge_weights",
                                        num_cores=num_cores,
                                        python_bin=NULL,
                                        df_points=NULL,NA_flag=NA_flag_val,
                                        file_format=file_format,out_suffix=out_suffix_str,
                                        out_dir=out_dir)

  mosaic_unweighted_obj <- mosaicFiles(lf_mosaic1,mosaic_method="unweighted",
                                        num_cores=num_cores,
                                        python_bin=NULL,
                                        df_points=NULL,NA_flag=NA_flag_val,
                                        file_format=file_format,out_suffix=out_suffix_str,
                                        out_dir=out_dir)

  list_mosaic_obj[[i]] <- list(unweighted=mosaic_unweighted_obj,edge=mosaic_edge_obj)
}

#####################
###### PART 2: Analysis and figures for the outputs of mosaic function #####

#### compute and aspect and slope with figures
list_lf_mosaic_obj <- vector("list",length(day_to_mosaic))
lf_mean_mosaic <- vector("list",length(mosaicing_method))#2methods only
l_method_mosaic <- vector("list",length(mosaicing_method))
list_out_suffix <- vector("list",length(mosaicing_method))

for(i in 1:length(day_to_mosaic)){
  list_lf_mosaic_obj[[i]] <- list.files(path=out_dir,pattern=paste("*",day_to_mosaic[i],
                                                                   "_.*.RData",sep=""))
  lf_mean_mosaic[[i]] <- unlist(lapply(list_lf_mosaic_obj[[i]],function(x){load_obj(x)[["mean_mosaic"]]}))
  l_method_mosaic[[i]] <- paste(unlist(lapply(list_lf_mosaic_obj[[i]],function(x){load_obj(x)[["method"]]})),day_to_mosaic[i],sep="_")
  list_out_suffix[[i]] <- unlist(paste(l_method_mosaic[[i]],day_to_mosaic[[i]],out_suffix,sep="_"))
}


list_param_plot_mosaic <- list(lf_mosaic=unlist(lf_mean_mosaic),
                               method=unlist(l_method_mosaic),
                               out_suffix=unlist(list_out_suffix))
#undebug(plot_mosaic)
#plot_mosaic(1,list_param=list_param_plot_mosaic)
num_cores <- 4
l_png_files <- mclapply(1:length(unlist(lf_mean_mosaic)),FUN=plot_mosaic,
                        list_param= list_param_plot_mosaic,
                        mc.preschedule=FALSE,mc.cores = num_cores)

####################
#### Now difference figures...

lf_obj1 <- list.files(path=out_dir,pattern="*unweighted.*.RData")
lf_obj2 <- list.files(path=out_dir,pattern="*edge_.*.RData")

lf1 <- unlist(lapply(lf_obj1,function(x){load_obj(x)[["mean_mosaic"]]}))
lf2 <- unlist(lapply(lf_obj2,function(x){load_obj(x)[["mean_mosaic"]]}))

out_suffix_str <- paste(paste(mosaicing_method,collapse="_"),day_to_mosaic,out_suffix,sep="_")

list_param_plot_diff <- list(lf1=lf1,lf2=lf2,out_suffix=out_suffix_str)

#debug(plot_diff_raster)
#plot_diff_raster(1,list_param=list_param_plot_diff)

num_cores <- 2
l_diff_png_files <- mclapply(1:length(lf1),FUN=plot_diff_raster,list_param= list_param_plot_diff,
                        mc.preschedule=FALSE,mc.cores = num_cores)


###############
##### Get all the tiles togheter merged

#ls -ltr ./reg*/*/*mean*.tif | wc
in_dir1 <- "/data/project/layers/commons/NEX_data/mosaicing_data_test"
lf_unweighted_20100831<- list.files(path=in_dir1,pattern="r_m_mean_20100831.*.tif",recursive=T,full.names=T)
lf_edge_weighted_20100831<- list.files(path=in_dir1,pattern="r_.*.edge.*._mean_20100831.*.tif",recursive=T,full.names=T)
lf_unweighted_20100901<- list.files(path=in_dir1,pattern="r_m_mean_20100901.*.tif",recursive=T,full.names=T)
lf_edge_weighted_20100901<- list.files(path=in_dir1,pattern="r_.*.edge.*.mean_20100901.*.tif",recursive=T,full.names=T)

output_fnames <- c("mean_unweighted_world_20100831_global_analyses_07012015.tif",
                   "mean_eged_weighted_world_20100831_global_analyses_07012015.tif",
                   "mean_unweighted_world_20100901_global_analyses_07012015.tif",
                   "mean_edge_weighted_world_20100901_global_analyses_07012015.tif"
                   )

list_lf <- list(lf_unweighted_20100831,lf_edge_weighted_20100831,lf_unweighted_20100901,lf_edge_weighted_20100901)
out_dir_str <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/mosaic_world_07012015"
for(i in 1:length(output_fnames)){
  rast_ref <- file.path(out_dir_str,output_fnames[i]) #this is a the ref ouput file
  lf_to_mosaic <- list_lf[[i]]
  cmd_str <- paste("python","/usr/bin/gdal_merge.py","-o ",rast_ref,paste(lf_to_mosaic,collapse=" ")) 
  system(cmd_str)
}
 
list_lf_m <- list.files(path=out_dir_str,pattern="mean.*.world.*.global_analyses_07012015.tif",full.names=T)

reg_name <- "world"
l_dates <- c("unweighted_20100831","edge_weighted_20100831","unweighted_20100901","edge_weighted_20100901")
list_param_plot_daily_mosaics <- list(list_lf_m,reg_name,out_dir_str,out_suffix,l_dates)
names(list_param_plot_daily_mosaics) <- c("lf_m","reg_name","out_dir_str","out_suffix","l_dates")

#debug(plot_daily_mosaics)
#test<- plot_daily_mosaics(1,list_param_plot_daily_mosaics)
num_cores <- 4
lf_plot <- mclapply(1:length(l_dates),FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics,
                    mc.preschedule=FALSE,mc.cores = num_cores)

##################### END OF SCRIPT ######################

