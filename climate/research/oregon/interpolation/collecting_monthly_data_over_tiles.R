####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#It shows how to load the monthly data at station locations.

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

load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

extract_list_from_list_obj<-function(obj_list,list_name){
  #Create a list of an object from a given list of object using a name prodived as input
  
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp<-obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]]<-tmp
  }
  return(list_tmp) #this is  a data.frame
}

#This extract a data.frame object from raster prediction obj and combine them in one data.frame 
extract_from_list_obj<-function(obj_list,list_name){
  #extract object from list of list. This useful for raster_prediction_obj
  library(plyr)
  
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp<-obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]]<- as.data.frame(tmp) #if spdf
  }
  tb_list_tmp<-do.call(rbind.fill,list_tmp) #long rownames
  #tb_list_tmp<-do.call(rbind,list_tmp) #long rownames
  
  return(tb_list_tmp) #this is  a data.frame
}


##############################
#### Parameters and constants  

#in_dir1 <- "/data/project/layers/commons/NEX_data/test_run1_03232014/output" #On Atlas
in_dir1 <- "/nobackupp4/aguzman4/climateLayers/output" #On NEX
in_dir_list <- list.dirs(path=in_dir1) #get the list of directories with resutls by 10x10 degree tiles
#in_dir_list <- as.list(in_dir_list[-1])
in_dir_list <- in_dir_list[grep("output",basename(in_dir_list),invert=TRUE)] #the first one is the in_dir1
in_dir_list <- in_dir_list[grep("shapefiles",basename(in_dir_list),invert=TRUE)] 
#the first one is the in_dir1
# the last directory contains shapefiles 

##raster_prediction object : contains testing and training stations with RMSE and model object

list_raster_obj_files <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="^raster_prediction_obj.*.RData",full.names=T)})

names(list_raster_obj_files)<- paste("tile",1:length(list_raster_obj_files),sep="_")
                                                                        
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
                                     
###################### PART I: Generate tables to collect information over all tiles in North America ##########

##Quick exploration of raster object

robj1 <- load_obj(list_raster_obj_files[[12]]) #load raster object
names(robj1)  #list the content of raster object, it is a R list
names(robj1$clim_method_mod_obj[[1]]$data_month) # monthly data for January
#names(robj1$validation_mod_month_obj[[1]]$data_s) # monthly for January with predictions
robj1$tb_diagnostic_v[1:10,] #first 10 rows of accuarcy metrics per day and model (for specific tile)
robj1$summary_metrics_v # accuracy averages per model (for specific tile)

#load data_month for specific tiles
data_month <- extract_from_list_obj(robj1$clim_method_mod_obj,"data_month")

names(data_month) #this contains LST means (mm_1, mm_2 etc.) as well as TMax and other info

#problem with tile 12...the raster ojbect has missing sub object
#data_month_list <- lapply(1:length(list_raster_obj_files),x=list_raster_obj_files,
#                          FUN=function(i,x){x<-load_obj(x[[i]]);
#                                            extract_from_list_obj(x$validation_mod_month_obj,"data_s")})                           

data_month_list <- lapply(1:length(list_raster_obj_files),x=list_raster_obj_files,
                          FUN=function(i,x){x<-load_obj(x[[i]]);
                                            extract_from_list_obj(x$clim_method_mod_obj,"data_month")})                           

names(data_month_list) <- paste("tile","_",1:length(data_month_list),sep="")


#names(data_month_list) <- basename(in_dir_list) #use folder id instead

tile_id <- lapply(1:length(data_month_list),
                  FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=data_month_list)
data_month_NAM <- do.call(rbind.fill,data_month_list) #combined data_month for "NAM" North America
data_month_NAM$tile_id <- unlist(tile_id)

#plot(mm_01 ~ elev_s,data=data_month_NAM) #Jan across all tiles
#plot(mm_06 ~ elev_s,data=data_month_NAM) #June across all tiles
#plot(TMax ~ mm_01,data=data_month_NAM) #monthly tmax averages and LST across all tiles

###### Additional information ######

# #LAND COVER INFORMATION (Tunamu et al., Jetz lab)

# LC1: Evergreen/deciduous needleleaf trees
# LC2: Evergreen broadleaf trees
# LC3: Deciduous broadleaf trees
# LC4: Mixed/other trees
# LC5: Shrubs
# LC6: Herbaceous vegetation
# LC7: Cultivated and managed vegetation
# LC8: Regularly flooded shrub/herbaceous vegetation
# LC9: Urban/built-up
# LC10: Snow/ice
# LC11: Barren lands/sparse vegetation
# LC12: Open water

## LST information: mm_01, mm_02 ...to mm_12 are monthly mean LST at station locaitons
## LST information: nobs_01, nobs_02 ... to nobs_12 number of valid obs used in mean LST averages
## TMax : monthly mean tmax at meteorological stations
## nbs_stt: number of stations used in the monthly mean tmax at stations
###

