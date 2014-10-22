####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for manuscript analyses,tables and figures #######################################
#This script uses the worklfow code applied to the Oregon case study. Daily methods (GAM,GWR, Kriging) are tested with
#different covariates using two baselines. Accuracy methods are added in the the function script to evaluate results.
#Figures, tables and data for the contribution of covariate paper are also produced in the script.
#AUTHOR: Benoit Parmentier                                                                      
#MODIFIED ON: 09/11/2014            
#Version: 5
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
library(gridExtra)                           # Combining lattice plots
library(colorRamps)                          # Palette/color ramps for symbology
#Additional libraries not used in workflow
library(pgirmess)                            # Krusall Wallis test with mulitple options, Kruskalmc {pgirmess}  
library(ncf)                                 # No paramtric covariance function

#### FUNCTION USED IN SCRIPT

function_analyses_paper1 <-"contribution_of_covariates_paper_interpolation_functions_07182014.R" #first interp paper
function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_10062014.R"

##############################
#### Parameters and constants  

script_path<-"/home/parmentier/Data/IPLANT_project/env_layers_scripts/" #path to script
source(file.path(script_path,function_analyses_paper1)) #source all functions used in this script 1.
source(file.path(script_path,function_analyses_paper2)) #source all functions used in this script 2.

#Multi time scale - CAI: gam, kriging, gwr
in_dir4 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gam_cai_lst_comb5_11032013"
raster_obj_file_4 <- "raster_prediction_obj_gam_CAI_dailyTmax_365d_gam_cai_lst_comb5_11032013.RData"

raster_prediction_obj_4 <- load_obj(file.path(in_dir4,raster_obj_file_4))

names(raster_prediction_obj_4)
names(raster_prediction_obj_4$method_mod_obj[[1]])

(raster_prediction_obj_4$method_mod_obj[[1]]$data_month_v) #no holdout...
data_month <- raster_prediction_obj_4$method_mod_obj[[1]]$data_month_s #no holdout...
dim(data_month)
test <- zerodist(data_month, zero = 0.0, unique.ID = FALSE)

target_max_nb <- 200
target_min_nb <- 100
max_dist <- 10
min_dist <- 0
step_dist <- 1000
target_range_nb <- c(target_min_nb,target_max_nb)
#debug(sub_sampling_by_dist)
#First increase distance till 5km
#then use random sampling...to get the extact target?
test <- sub_sampling_by_dist(target_range_nb,dist=min_dist,max_dist=max_dist,step=step_dist,data_in=data_month)
test <- sub_sampling_by_dist(target_range_nb,dist=min_dist,max_dist=NULL,step=step_dist,data_in=data_month)

dist_range <- c(0,5000) 

sub_sampling_by_dist_nb_stat(target_range_nb=c(10000,10000),dist_range,data_in,sampling=T){

  
sub_sampling_by_dist <- function(target_range_nb=c(10000,10000),dist=0.0,max_dist=NULL,step,data_in){
  data <- data_in
  target_min_nb <- target_range_nb[1]
  station_nb <- nrow(data_in)
  if(is.null(max_dist)){
    while(station_nb > target_min_nb){
      data <- remove.duplicates(data, zero = dist) #spatially sub sample...
      dist <- dist + step
      station_nb <- nrow(data)
    }
  }
  if(!is.null(max_dist)){
    while(station_nb > target_min_nb || dist < max_dist){ 
      d#} #|| nrow > 0){
      #test <- zerodist(data, zero = 0.0, unique.ID = FALSE)
      #test <- remove.duplicates(data_month, zero = 5000)
      data <- remove.duplicates(data, zero = dist) #spatially sub sample...
      dist <- dist + step
      station_nb <- nrow(data)
    }
  }
  
  obj_sub_sampling <- list(data,dist)
  names(obj_sub_sampling) <- c("data","dist")
  return(obj_sub_sampling)
}

sub_sampling_by_dist_nb_stat <- function(target_range_nb=c(10000,10000),dist_range,data_in,sampling=T){
  
  data <- data_in
  min_dist <- dist_range[1]
  max_dist <- dist_range[2]
  
  if(sampling==T){
    dat <- sub_sampling_by_dist(target_range_nb,dist=min_dist,max_dist=max_dist,step=step_dist,data_in=data_month)
    data_out <-sample(dat$data, target_range_nb[2], replace = FALSE, prob = NULL)
    data_out <- list(data_out,dat$dist,dat$data)
    data_out <- c("data","dist","data_dist")
  }
  if(sampling!=T){
    data_out <- sub_sampling_by_dist(target_range_nb,dist=min_dist,max_dist=NULL,step=step_dist,data_in=data_month)
  }
  return(data_out)
}

#sub_sampling_by_dist <- function(target_range_nb=c(10000,10000),dist=0.0,step,data_in){
#  data <- data_in
#  target_min_nb <- target_range_nb[1]
#  station_nb <- nrow(data_in)
#  while(station_nb > target_min_nb){  #} #|| nrow > 0){
#      #test <- zerodist(data, zero = 0.0, unique.ID = FALSE)
#      #test <- remove.duplicates(data_month, zero = 5000)
#      data <- remove.duplicates(data, zero = dist) #spatially sub sample...
#      dist <- dist + step
#      station_nb <- nrow(data)
#  }
#  obj_sub_sampling <- list(data,dist)
#  names(obj_sub_sampling) <- c("data","dist")
#  return(obj_sub_sampling)
#}

############ END OF SCRIPT #########