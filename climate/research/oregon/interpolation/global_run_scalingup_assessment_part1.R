####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Accuracy methods are added in the the function scripts to evaluate results.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 03/23/2014  
#MODIFIED ON: 03/23/2014            
#Version: 1
#PROJECT: Environmental Layers project                                     
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

#### FUNCTION USED IN SCRIPT

function_analyses_paper1 <-"contribution_of_covariates_paper_interpolation_functions_10222013.R" #first interp paper
function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_03182014.R"

load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

##############################
#### Parameters and constants  

script_path<-"/home/parmentier/Data/IPLANT_project/env_layers_scripts/" #path to script
source(file.path(script_path,function_analyses_paper1)) #source all functions used in this script 1.
source(file.path(script_path,function_analyses_paper2)) #source all functions used in this script 2.


#in_dir1 <- "/data/project/layers/commons/NEX_data/test_run1_03232014/output" #On Atlas
in_dir1 <- "/nobackupp4/aguzman4/climateLayers/output" #On NEX
in_dir_list <- list.dirs(path=in_dir1) #get the list of directories with resutls by 10x10 degree tiles
#in_dir_list <- as.list(in_dir_list[-1])
in_dir_list <- in_dir_list[-1] #the first one is the in_dir1
in_dir_list <- in_dir_list[-25] # the last directory contains shapefiles #]in_dir_list[-1]


##raster_prediction object : contains testing and training stations with RMSE and model object

list_raster_obj_files <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="^raster_prediction_obj.*.RData",full.names=T)})

names(list_raster_obj_files)<- paste("tile",1:length(list_raster_obj_files),sep="_")
                                     
y_var_name <- "dailyTmax"
out_prefix<-"run1_NA_analyses_03232013"
#out_dir<-"/data/project/layers/commons/NEX_data/"
out_dir <- "/nobackup/bparmen1/"
out_dir <-paste(out_dir,"_",out_prefix,sep="")

#system("ls /nobackup/bparmen1")

if (!file.exists(out_dir)){
  dir.create(out_dir)
  #} else{
  #  out_path <-paste(out_path..)
}
setwd(out_dir)
                                   
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
                                     
###################### PART I: Generate tables to collect information over all tiles in North America ##########
#Table 1: Average accuracy metrics
#Table 2: daily accuracy metrics for all tiles

##Quick exploration of raster object
robj1 <- load_obj(list_raster_obj_files[[1]])
names(robj1)
names(robj1$clim_method_mod_obj[[1]]$data_month) #for January
names(robj1$validation_mod_month_obj[[1]]$data_s) #for January with predictions

data_month_list <- lapply(list_raster_obj_files,FUN=function(x){x<-load_obj(x);x[["summary_metrics_v"]]$avg})                           

robj1$tb_diagnostic_v[1:10,] #first 10 rows of accuarcy metrics per day and model (for specific tile)
robj1$summary_metrics_v #first 10 rows of accuarcy metrics per day and model (for specific tile)

#summary_metrics_v_list <- lapply(list_raster_obj_files,FUN=function(x){x<-load_obj(x);x[["summary_metrics_v"]]$avg$rmse})                           
summary_metrics_v_list <- lapply(list_raster_obj_files,FUN=function(x){x<-load_obj(x);x[["summary_metrics_v"]]$avg})                           
#summary_metrics_v_NA <- do.call(rbind,summary_metrics_v_list) #create a df for NA tiles with all accuracy metrics
summary_metrics_v_NA <- do.call(rbind.fill,summary_metrics_v_list) #create a df for NA tiles with all accuracy metrics
tile_id <- lapply(1:length(summary_metrics_v_list),FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=summary_metrics_v_list)
summary_metrics_v_NA$tile_id <-unlist(tile_id)
summary_metrics_v_NA$n <- as.integer(summary_metrics_v_NA$n)
write.table(as.data.frame(summary_metrics_v_NA),file=paste("summary_metrics_v_NA_",out_prefix,".txt",sep=""),sep=",")
#Function to collect all the tables from tiles into a table

tb_diagnostic_v_list <- lapply(list_raster_obj_files,FUN=function(x){x<-load_obj(x);x[["tb_diagnostic_v"]]})                           
names(tb_diagnostic_v_list)

tb_diagnostic_v_NA <- do.call(rbind.fill,tb_diagnostic_v_list) #create a df for NA tiles with all accuracy metrics
tile_id <- lapply(1:length(tb_diagnostic_v_list),FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=tb_diagnostic_v_list)
tb_diagnostic_v_NA$tile_id <- unlist(tile_id) #adding identifier for tile
write.table((tb_diagnostic_v_NA),file=paste("tb_diagnostic_v_NA","_",out_prefix,".txt",sep=""),sep=",")

### Create a combined shape file for all region?

#get centroid
#plot the average RMSE at the centroid??
#quick kriging for every tile?
                    
### Create a combined boxplot for every tile (can also do that in pannel)
### Create a quick mosaic (mean method...)
#mean predicitons
#mean of kriging error?
                    
tb <- read.table("/data/project/layers/commons/NEX_data/test_run1_03232014/tb_diagnostic_v_NA_run1_NA_analyses_03232013.txt",sep=",")

boxplot(rmse~tile_id,data=subset(tb,tb$pred_mod=="mod1"))
bwplot(rmse~as.factor(tile_id), data=subset(tb,tb$pred_mod=="mod1"))
