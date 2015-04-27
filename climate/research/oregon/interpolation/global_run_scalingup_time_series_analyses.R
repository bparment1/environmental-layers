##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX ##############################
#This script uses outputs from the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#We examine outputs from time series mosaiced regionally and globablly.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 04/27/2015  
#MODIFIED ON: 04/30/2015            
#Version: 1
#PROJECT: NCEAS/IPLANT/NASA Environmental Layers project     
#COMMENTS: analyses for run 10 global analyses,all regions 1500x4500km and other tiles
#TODO:
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
library(gridExtra)
#Additional libraries not used in workflow
library(pgirmess)                            # Krusall Wallis test with mulitple options, Kruskalmc {pgirmess}  
library(colorRamps)                          # Palettes for raster maps
library(zoo)                                 # time series analysis objects and tools/functionalities
library(xts)                                 # time series analysis extended tools

#### FUNCTION USED IN SCRIPT

function_analyses_paper1 <-"contribution_of_covariates_paper_interpolation_functions_07182014.R" #first interp paper
function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_08132014.R"

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

 #Remove models that were not fitted from the list
#All modesl that are "try-error" are removed
remove_errors_list<-function(list_items){
  
  #This function removes "error" items in a list
  list_tmp<-list_items
    if(is.null(names(list_tmp))){
    names(list_tmp) <- paste("l",1:length(list_tmp),sep="_")
    names(list_items) <- paste("l",1:length(list_tmp),sep="_")
  }

  for(i in 1:length(list_items)){
    if(inherits(list_items[[i]],"try-error")){
      list_tmp[[i]]<-0
    }else{
    list_tmp[[i]]<-1
   }
  }
  cnames<-names(list_tmp[list_tmp>0])
  x <- list_items[match(cnames,names(list_items))]
  return(x)
}


############################################
#### Parameters and constants  

#on ATLAS

in_dir <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_04172015/mosaics"

y_var_name <- "dailyTmax" #PARAM1
interpolation_method <- c("gam_CAI") #PARAM2
out_suffix <- "ts_run10_1500x4500_global_analyses_04172015" #PARAM3
out_dir <- in_dir #PARAM4
create_out_dir_param <- TRUE #PARAM 5
  
#CRS_locs_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
CRS_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84 #CONSTANT1
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

proj_str<- CRS_WGS84 #PARAM 8 #check this parameter
file_format <- ".rst" #PARAM 9
NA_value <- -9999 #PARAM10
NA_flag_val <- NA_value
                                   
tile_size <- "1500x4500" #PARAM 11

region_name <- "North_America" #PARAM 13

########################## START SCRIPT ##############################


####### PART 1: Read in data ########

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

setwd(out_dir)

################## PLOTTING  MOSAICS AND TIME SERIES ################

lf_mosaic_pred_1500x4500 <-list.files(path=file.path(in_dir),    
           pattern=paste("reg1.*.tif$",sep=""),full.names=T) 
r_stack <- stack(lf_mosaic_pred_1500x4500)#

summary_metrics_v <- read.table(file=file.path("/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_04172015",paste("summary_metrics_v2_NA_","run10_1500x4500_global_analyses_04172015",".txt",sep="")),sep=",")
coordinates(summary_metrics_v) <- cbind(summary_metrics_v$lon,summary_metrics_v$lat)

r1 <- raster(lf_mosaic_pred_1500x4500[11]) 

idx <- seq(as.Date('2010-01-01'), as.Date('2010-12-31'), 'day')
date_l <- strptime(idx[1], "%Y%m%d") # interpolation date being processed
dates_l <- format(idx, "%Y%m%d") # interpolation date being processed

## Plot mosaics for North America (region1) for daily predictions in 2010

res_pix <- 480

col_mfrow <- 2
row_mfrow <- 1
  
#  png(filename=paste("Figure10_clim_world_mosaics_day_","_",date_proc,"_",tile_size,"_",out_suffix,".png",sep=""),
#    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
png(filename=paste("Figure1_","time_series_step_in_raster_mosaics",dates_l[11],"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r1)
plot(summary_metrics_v,add=T)
text(summary_metrics_v,summary_metrics_v$tile_id,cex=1.4)

dev.off()

## Get pixel time series at centroids of tiles used in the predictions

df_ts_pixel <- extract(r_stack,summary_metrics_v,df=T,sp=F)

df_ts_pixel <- cbind(summary_metrics_v,df_ts_pixel)

#make a function later on?

#inputs
list_selected_pix <- c("tile_4","tile_6","tile_8","tile_11","tile_14","tile_3","tile_5","tile_7","tile_38","tile_12")
list_pix <- vector("list",length=length(list_selected_pix))
#idx <- seq(as.Date('2010-01-01'), as.Date('2010-12-31'), 'day')
#df_ts_pix

#Select one pix to profile/plot

df_ts_pix <- subset(df_ts_pixel,pred_mod=="mod1")

for(i in 1:length(list_selected_pix)){
  
  selected_pix <- list_selected_pix[i]
  data_pixel <- subset(df_ts_pix,tile_id==selected_pix)
  pix <- t(data_pixel[1,24:388])
  d_z <- zoo(pix,idx) #make a time series ...
  list_pix[[i]] <- pix

  res_pix <- 480

  col_mfrow <- 2
  row_mfrow <- 1
  
  #  png(filename=paste("Figure10_clim_world_mosaics_day_","_",date_proc,"_",tile_size,"_",out_suffix,".png",sep=""),
  #    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  png(filename=paste("Figure2_","pixel_profile_",selected_pix,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)


  plot(d_z,type="b")
  title(paste("Pixel time series",selected_pix,sep=" "))

 dev.off()
}

data_dz <- do.call(cbind,list_pix)
colnames(data_dz) <- list_selected_pix
data_dz <- zoo(data_dz,idx)

list_selected_pix <- c("tile_4","tile_6","tile_8","tile_11","tile_14","tile_3","tile_5","tile_7","tile_38","tile_12")
df_ts_pix2 <- subset(df_ts_pix,tile_id%in% list_selected_pix)

pix_data <- t(df_ts_pix2[,24:388])

d_z2 <- zoo(pix_data,idx)
names(d_z2)<-

##################### END OF SCRIPT ######################
