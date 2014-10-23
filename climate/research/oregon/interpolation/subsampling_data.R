####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script to prune and sample stations data #######################################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA and NCEAS UCSB.
#The purpose is sub sample stations data in areas where there are a very high number of stations.
#
#AUTHOR: Benoit Parmentier                                                                      
#CREATED ON: 10/16/2014            
#MODIFIED ON: 10/23/2014            
#Version: 1
#
#PROJECT: Environmental Layers project  NCEAS-NASA
#TO DO:
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

function_analyses_paper1 <- "contribution_of_covariates_paper_interpolation_functions_07182014.R" #first interp paper
function_analyses_paper2 <- "multi_timescales_paper_interpolation_functions_10062014.R"

sub_sampling_by_dist <- function(target_range_nb=c(10000,10000),dist=0.0,max_dist=NULL,step,data_in){
  #Function to select stations data that are outside a specific spatial range from each other
  #Parameters:
  #max_dist: maximum spatial distance at which to stop the pruning
  #min_dist: minimum distance to start pruning the data
  #step: spatial distance increment

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
    
    while(station_nb > target_min_nb & dist < max_dist){ 
      data <- remove.duplicates(data, zero = dist) #spatially sub sample...
      id_rm <- zerodist(data, zero = dist, unique.ID = FALSE)
      data_rm <- data[id_rm,]
      dist <- dist + step
      station_nb <- nrow(data)
    }
  }
  
  obj_sub_sampling <- list(data,dist)
  names(obj_sub_sampling) <- c("data","dist")
  return(obj_sub_sampling)
}

sub_sampling_by_dist_nb_stat <- function(target_range_nb,dist_range,data_in,sampling=T,combined=F){
  
  data <- data_in
  min_dist <- dist_range[1]
  max_dist <- dist_range[2]
  
  if(sampling==T){
    dat <- sub_sampling_by_dist(target_range_nb,dist=min_dist,max_dist=max_dist,step=step_dist,data_in=data_month)
    ind_s1  <- sample(nrow(dat$data), size=target_range_nb[1], replace = FALSE, prob = NULL)
    ind_s2 <- setdiff(1:nrow(dat$data), ind_s1)
    data_out <- dat$data[ind_s1,] #selected the randomly sampled stations
    data_removed <- dat[ind_s2,]
    
    #Find the corresponding 
    #data_sampled<-ghcn.subsets[[i]][ind.training,] #selected the randomly sampled stations

    data_out <- list(data_out,dat$dist,data_removed,dat$data)
    data_out <- c("data","dist","data_removed","data_dist")
  }
  if(sampling!=T){
    dat <- sub_sampling_by_dist(target_range_nb,dist=min_dist,max_dist=NULL,step=step_dist,data_in=data_month)
    #
    data_out <- list(dat$data,dat$dist,data_removed)
    data_out <- c("data","dist","data_removed")
    
  }
  return(data_out)
}

debug(sub_sampling_by_dist_nb_stat)
test3 <- sub_sampling_by_dist_nb_stat(target_range_nb=c(100,200),dist_range=c(0,1000),data_in=data_month,sampling=T,combined=F)


#n<-nrow(ghcn.subsets[[i]])
#prop<-(sampling_dat$prop[i])/100
#ns<-n-round(n*prop)   #Create a sample from the data frame with 70% of the rows
#nv<-n-ns              #create a sample for validation with prop of the rows
#ind.training <- sample(nrow(ghcn.subsets[[i]]), size=ns, replace=FALSE) #This selects the index position for 70% of the rows taken randomly
#ind.testing <- setdiff(1:nrow(ghcn.subsets[[i]]), ind.training)
#Find the corresponding 
#data_sampled<-ghcn.subsets[[i]][ind.training,] #selected the randomly sampled stations
#station_id.training<-data_sampled$station     #selected id for the randomly sampled stations (115)
#Save the information
#sampling[[i]]<-ind.training #index of training sample from data.frame
#sampling_station_id[[i]]<- station_id.training #station ID for traning samples

##############################
#### Parameters and constants  

script_path<-"/home/parmentier/Data/IPLANT_project/env_layers_scripts/" #path to script
source(file.path(script_path,function_analyses_paper1)) #source all functions used in this script 1.
source(file.path(script_path,function_analyses_paper2)) #source all functions used in this script 2.

#Multi time scale - CAI: gam, kriging, gwr
in_dir4 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gam_cai_lst_comb5_11032013" #data used in multi-timescale paper
raster_obj_file_4 <- "raster_prediction_obj_gam_CAI_dailyTmax_365d_gam_cai_lst_comb5_11032013.RData"

###################################################################
############### BEGIN SCRIPT ##################################

raster_prediction_obj_4 <- load_obj(file.path(in_dir4,raster_obj_file_4)) #load raster object to get monthly data

names(raster_prediction_obj_4)
names(raster_prediction_obj_4$method_mod_obj[[1]])

raster_prediction_obj_4$method_mod_obj[[1]]$data_month_v #no holdout...
data_month <- raster_prediction_obj_4$method_mod_obj[[1]]$data_month_s #no holdout...
dim(data_month) #193x70, there are 70 stations in this particular case

plot(data_month)

#set up input parameters

target_max_nb <- 200 #this is not actually used yet in the current implementation
target_min_nb <- 100 #this is the target number of stations we would like
max_dist <- 1000 # the maximum distance used for pruning ie removes stations that are closer than 1000m 
min_dist <- 0    #minimum distance to start with
step_dist <- 1000 #iteration step to remove the stations
target_range_nb <- c(target_min_nb,target_max_nb) #target range of number of stations
#First increase distance till 5km
#then use random sampling...to get the extact target?
#First test with maximum distance of 100m
test1 <- sub_sampling_by_dist(target_range_nb,dist=min_dist,max_dist=max_dist,step=step_dist,data_in=data_month)
#Second test with no maximum distance: the process of spatial pruning only stops when the number of stations is below
#the minimum target number
test2 <- sub_sampling_by_dist(target_range_nb,dist=min_dist,max_dist=NULL,step=step_dist,data_in=data_month)
#set max distance to NULL

dim(test1$data) #192 stations selected
test1$dist # for distance of 1000 m (max_dist)
dim(test2$data) #97 stations selected 
test2$dist # for distance of 31,000 m (no max_dist is set)

dist_range <- c(0,5000) 

#Now use the other function to sample the station data points:

sub_sampling_by_dist_nb_stat(target_range_nb=c(10000,10000),dist_range,data_in,sampling=T)


#n<-nrow(ghcn.subsets[[i]])
#prop<-(sampling_dat$prop[i])/100
#ns<-n-round(n*prop)   #Create a sample from the data frame with 70% of the rows
#nv<-n-ns              #create a sample for validation with prop of the rows
#ind.training <- sample(nrow(ghcn.subsets[[i]]), size=ns, replace=FALSE) #This selects the index position for 70% of the rows taken randomly
#ind.testing <- setdiff(1:nrow(ghcn.subsets[[i]]), ind.training)
#Find the corresponding 
#data_sampled<-ghcn.subsets[[i]][ind.training,] #selected the randomly sampled stations
#station_id.training<-data_sampled$station     #selected id for the randomly sampled stations (115)
#Save the information
#sampling[[i]]<-ind.training #index of training sample from data.frame
#sampling_station_id[[i]]<- station_id.training #station ID for traning samples
  
############ END OF SCRIPT #########