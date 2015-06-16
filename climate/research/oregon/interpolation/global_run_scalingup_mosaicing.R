##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Different options to explore mosaicing are tested.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 04/14/2015  
#MODIFIED ON: 06/16/2015            
#Version: 4
#PROJECT: Environmental Layers project     
#COMMENTS: analyses run for reg5 for test of mosaicing using 1500x4500km and other tiles
#TODO:
#1) Split functions and master script
#2) Make this is a script/function callable from the shell/bash
#3) Check image format for tif
#4) generalize to run dates and region fast

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

#function_analyses_paper1 <-"contribution_of_covariates_paper_interpolation_functions_07182014.R" #first interp paper
#function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_08132014.R"

function_mosaicing <-"multi_timescales_paper_interpolation_functions_08132014.R"

in_dir_script <-"/home/parmentier/Data/IPLANT_project/env_layers_scripts"
source(file.path(in_dir_script,function_mosaicing))

############################################
#### Parameters and constants  

#Data is on ATLAS: reg4 (South America)

in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test" #reg1 is North America and reg5 is Africa
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg1"
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg2" #Europe
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg5" #Africa

y_var_name <- "dailyTmax" #PARAM1
interpolation_method <- c("gam_CAI") #PARAM2
region_name <- "reg2" #PARAM 13 #reg4 South America, Africa reg5,Europe reg2, North America reg1, Asia reg3

out_suffix <- paste(region_name,"_","mosaic_run10_1500x4500_global_analyses_06152015",sep="") 
#PARAM3
out_dir <- in_dir #PARAM4
create_out_dir_param <- TRUE #PARAM 5

mosaic_plot <- FALSE #PARAM6

#if daily mosaics NULL then mosaicas all days of the year
day_to_mosaic <- c("20100831",
                   "20100901") #PARAM7
  
#CRS_locs_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
CRS_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84 #CONSTANT1
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

proj_str<- CRS_WGS84 #PARAM 8 #check this parameter
file_format <- ".tif" #PARAM 9
NA_value <- -9999 #PARAM10
NA_flag_val <- NA_value
     
num_cores <- 11                              
tile_size <- "1500x4500" #PARAM 11
mulitple_region <- TRUE #PARAM 12

plot_region <- FALSE

########################## START SCRIPT ##############################


####### PART 1: Read in data and process data ########

#make this a loop?, fistt use sept 1, 2010 data
#out_suffix <- paste(day_to_mosaic[2],out_suffix,sep="_")

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

lf <- sub(".tif","",lf_mosaic2)
tx<-strsplit(as.character(lf),"_")

lat<- as.character(lapply(1:length(tx),function(i,x){x[[i]][13]},x=tx))
long<- as.character(lapply(1:length(tx),function(i,x){x[[i]][14]},x=tx))
lat <- as.character(lapply(1:length(lat),function(i,x){substr(x[[i]],2,nchar(x[i]))},x=lat)) #first number not in the coordinates

#Produce data.frame with centroids of each tiles...

df_centroids <- data.frame(long=as.numeric(long),lat=as.numeric(lat))
df_centroids$ID <- as.numeric(1:nrow(df_centroids))
coordinates(df_centroids) <- cbind(df_centroids$long,df_centroids$lat)
proj4string(df_centroids) <- projection(r1)
df_points <- df_centroids
#methods availbable:use_sine_weights,use_edge,use_linear_weights
out_suffix_str <- paste(day_to_mosaic[2],out_suffix,sep="_")

#debug(mosaicFiles)
mosaic_edge_20100901_obj <- mosaicFiles(lf_mosaic2,mosaic_method="use_edge_weights",
                                        num_cores=num_cores,
                                        python_bin=NULL,
                                        df_points=NULL,NA_flag=NA_flag_val,
                                        file_format=file_format,out_suffix=out_suffix_str,
                                        out_dir=out_dir)
#debug(mosaicFiles)
save(mosaic_unweighted_20100901_obj,file=file.path(out_dir,
                                                   paste(mosaic_method,"_","mosaic_obj_","20100901_",out_suffix,".RData",sep="")))

mosaic_unweighted_20100901_obj <- mosaicFiles(lf_mosaic2,mosaic_method="unweighted",
                                        num_cores=num_cores,
                                        python_bin=NULL,
                                        df_points=NULL,NA_flag=NA_flag_val,
                                        file_format=file_format,out_suffix=out_suffix_str,
                                        out_dir=out_dir)
mosaic_method <- "unweighted"
save(mosaic_edge_20100901_obj,file=file.path(out_dir,paste(mosaic_method,"_","mosaic_obj_","20100901",out_suffix,".RData")))

mosaic_edge_20100831_obj <- mosaicFiles(lf_mosaic1,mosaic_method="use_edge_weights",
                                        num_cores=num_cores,
                                        python_bin=NULL,
                                        df_points=NULL,NA_flag=NA_flag_val,
                                        file_format=file_format,out_suffix=out_suffix_str,
                                        out_dir=out_dir)

mosaic_unweighted_20100831_obj <- mosaicFiles(lf_mosaic1,mosaic_method="unweighted",
                                        num_cores=num_cores,
                                        python_bin=NULL,
                                        df_points=NULL,NA_flag=NA_flag_val,
                                        file_format=file_format,out_suffix=out_suffix_str,
                                        out_dir=out_dir)

r2_unweighted <-raster(mosaic_unweighted_20100901_obj$mean_mosaic)
r2_edge <-raster(mosaic_edge_20100901_obj$mean_mosaic)

r1_unweighted <-raster(mosaic_unweighted_20100831_obj$mean_mosaic)
r1_edge <-raster(mosaic_edge_20100831_obj$mean_mosaic)
plot(r1_edge)

#####################
###### PART 2: Analysis and figures for the outputs of mosaic function #####

#### compute and aspect and slope with figures
list_mosaic_unweighted <- list(mosaic_unweighted_20100831_obj,mosaic_edge_20100831_obj)
list_mosaic_edge <- list(mosaic_unweighted_20100901_obj,mosaic_edge_20100901_obj)

list_mosaiced_files <- c(list_mosaiced_files,r_m_mean_unweighted)
names(list_mosaiced_files2) <- c(names(list_mosaiced_files),"unweighted")

plot_mosaic <- function(f_mosaic,method,out_dir,out_stuffix){

  method_str <- method
  r_mosaic <- raster(f_mosaiced)

  r_mosaic_terrain <- terrain(r_mosaic,opt=c("slope","aspect"),unit="degrees")

  res_pix <- 1200
  col_mfrow <- 1 
  row_mfrow <- 0.8

  png(filename=paste("Figure2_mosaic_mean_",method_str,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_mosaic,main=paste("mosaic mean ",method_str,sep=""))

  dev.off()
  
  #### plot terrain to emphasize possible edges..
  res_pix <- 1200
  col_mfrow <- 1 
  row_mfrow <- 0.8

  png(filename=paste("Figure2_slope_mean_",method_str,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_mosaic_terrain,y=1,main=paste("slope mosaic mean ",method_str,sep=""))

  dev.off()

  png(filename=paste("Figure2_aspect_mean_",method_str,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_mosaic_terrain,y=2,main=paste("aspect mean ",method_str,sep=""))

  dev.off()
}

####################
#### Now difference figures...
r_m_edge_weighted_mean <- raster(list_mosaiced_files2[1])#edge
r_m_linear_weighted_mean <- raster(list_mosaiced_files2[2])#linear
r_m_sine_weighted_mean <- raster(list_mosaiced_files2[3])#sine  
r_m_unweighted_mean <- raster(list_mosaiced_files2[4])#unweighted

r_diff_linear_sine_weighted_mean <- r_m_linear_weighted_mean - r_m_sine_weighted_mean

res_pix <- 1200
col_mfrow <- 1 
row_mfrow <- 0.8

png(filename=paste("Figure2_diff_linear_sine_weigthed_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_diff_linear_sine_weighted_mean)

dev.off()

r_diff_linear_edge_weighted_mean <- r_m_linear_weighted_mean - r_m_edge_weighted_mean

png(filename=paste("Figure2_diff_linear_edge_weigthed_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_diff_linear_edge_weighted_mean)

dev.off()


#r_diff_linear_edge_weighted_mean <- r_m_linear_weighted_mean - r_m_edge_weighted_mean
r_diff_edge_sine_weighted_mean <- r_m_edge_weighted_mean - r_m_sine_weighted_mean

png(filename=paste("Figure2_diff_edge_sine_weigthed_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_diff_edge_sine_weighted_mean)

dev.off()

###############
##### Now compare to unweighted values

r_diff_unweighted_linear_weighted_mean <- r_m_mean - r_m_linear_weighted_mean 
r_diff_unweighted_sine_weighted_mean <- r_m_mean - r_m_sine_weighted_mean 
r_diff_unweighted_edge_weighted_mean <- r_m_mean - r_m_edge_weighted_mean 

png(filename=paste("Figure2_diff_unweighted_edge_weigthed_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_diff_unweighted_edge_weighted_mean)

dev.off()

png(filename=paste("Figure2_diff_unweighted_linear_weighted_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_diff_unweighted_linear_weighted_mean)

dev.off()

png(filename=paste("Figure2_diff_unweighted_sine_weigthed_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_diff_unweighted_sine_weighted_mean)

dev.off()

##################### END OF SCRIPT ######################

