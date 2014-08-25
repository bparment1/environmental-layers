####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for assessment of scaling up on NEX: part 1 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#The purpose is to create as set of functions to diagnose and assess quickly a set of predictd tiles.
#Part 1 create summary tables and inputs for figure in part 2 and part 3.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 03/23/2014  
#MODIFIED ON: 08/25/2014            
#Version: 3
#PROJECT: Environmental Layers project  
#TO DO:
# - generate delta and clim mosaic
# - generate monthly inputs data_month
# - generate table of number of observations per tile for use in map part 2
# - generate data_s and data_v inputs as giant table
# - generate accuracy for mosaic (part 2 and part3)
# - clean up

#First source file:
#source /nobackupp4/aguzman4/climateLayers/sharedModules/etc/environ.sh
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
  #create if does not exists: create the output dir as defined 
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  return(out_dir)
}

extract_list_from_list_obj<-function(obj_list,list_name){
  library(plyr)
  
  #Create a list of an object from a given list of object using a name prodived as input
  
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp <- obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]] <- as.data.frame(tmp) #deal with spdf...
  }
  return(list_tmp) #this is  a data.frame
}

#This extract a data.frame object from raster prediction obj and combine them in one data.frame 
extract_from_list_obj<-function(obj_list,list_name){
  #extract object from list of list. This useful for raster_prediction_obj
  library(plyr)
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp <- obj_list[[i]]
    if(inherits(tmp,"try-error")){     
      print(paste("no model generated or error in list",sep=" ")) #change message for any model type...
      list_tmp[[i]] <- NULL #double bracket to return data.frame
    }else{
      #tmp<-obj_list[[i]][[list_name]] #double bracket to return data.frame
      list_tmp[[i]] <- as.data.frame(tmp[[list_name]])
    }
    #
    #tmp<-obj_list[[i]][[list_name]] #double bracket to return data.frame
    #list_tmp[[i]]<- as.data.frame(tmp) #if spdf
  }
  #
  #list_tmp <-list_tmp[!is.null(list_tmp)]
  list_tmp <- list_tmp[unlist(lapply(list_tmp,FUN=function(x){!is.null(x)}))]

  tb_list_tmp<-do.call(rbind.fill,list_tmp) #long rownames
  #tb_list_tmp<-do.call(rbind,list_tmp) #long rownames 
  return(tb_list_tmp) #this is  a data.frame
}

## Function to mosaic modis or other raster images

mosaic_m_raster_list<-function(j,list_param){
  #This functions returns a subset of tiles from the modis grid.
  #Arguments: modies grid tile,list of tiles
  #Output: spatial grid data frame of the subset of tiles
  #Note that rasters are assumed to be in the same projection system!!
  
  #rast_list<-vector("list",length(mosaic_list))
  #for (i in 1:length(mosaic_list)){  
  # read the individual rasters into a list of RasterLayer objects
  # this may be changed so that it is not read in the memory!!!
  
  #parse output...
  
  #j<-list_param$j
  mosaic_list<-list_param$mosaic_list
  out_path<-list_param$out_path
  out_names<-list_param$out_rastnames
  file_format <- list_param$file_format
  NA_flag_val <- list_param$NA_flag_val
  out_suffix <- list_param$out_suffix
  ## Start
  
  if(class(mosaic_list[[j]])=="list"){
    m_list <- unlist(mosaic_list[[j]])
  }
  input.rasters <- lapply(m_list, raster)
  mosaiced_rast<-input.rasters[[1]]
  
  for (k in 2:length(input.rasters)){
    mosaiced_rast<-mosaic(mosaiced_rast,input.rasters[[k]], fun=mean)
    #mosaiced_rast<-mosaic(mosaiced_rast,raster(input.rasters[[k]]), fun=mean)
  }
  
  data_name<-paste("mosaiced_",sep="") #can add more later...
  #raster_name<-paste(data_name,out_names[j],".tif", sep="")
  raster_name<-paste(data_name,out_names[j],file_format, sep="")
  
  writeRaster(mosaiced_rast, NAflag=NA_flag_val,filename=file.path(out_path,raster_name),overwrite=TRUE)  
  #Writing the data in a raster file format...  
  rast_list<-file.path(out_path,raster_name)
  
  ## The Raster and rgdal packages write temporary files on the disk when memory is an issue. This can potential build up
  ## in long  loops and can fill up hard drives resulting in errors. The following  sections removes these files 
  ## as they are created in the loop. This code section  can be transformed into a "clean-up function later on
  ## Start remove
  #tempfiles<-list.files(tempdir(),full.names=T) #GDAL transient files are not removed
  #files_to_remove<-grep(out_suffix,tempfiles,value=T) #list files to remove
  #if(length(files_to_remove)>0){
  #  file.remove(files_to_remove)
  #}
  #now remove temp files from raster package located in rasterTmpDir
  removeTmpFiles(h=0) #did not work if h is not set to 0
  ## end of remove section
  
  return(rast_list)
}

### Function:
pred_data_info_fun <- function(k,list_data,pred_mod,sampling_dat_info){
  #Summarizing input info from sampling and df used in training/testing
    
  data <- list_data[[k]]
  sampling_dat <- sampling_dat_info[[k]]
  if(data!="try-error"){
    n <- nrow(data)
    n_mod <- vector("numeric",length(pred_mod))
    for(j in 1:length(pred_mod)){
      n_mod[j] <- sum(!is.na(data[[pred_mod[j]]]))
    }
    n <- rep(n,length(pred_mod))
    sampling_dat <- sampling_dat[rep(seq_len(nrow(sampling_dat)), each=length(pred_mod)),]
    row.names(sampling_dat) <- NULL
    df_n <- data.frame(n,n_mod,pred_mod)
    df_n <- cbind(df_n,sampling_dat)
  }else{        
    n <- rep(NA,length(pred_mod))
    n_mod <- vector("numeric",length(pred_mod))
    n_mod <- rep(NA,length(pred_mod))
    df_n <- data.frame(n,n_mod,pred_mod)
    sampling_dat <- sampling_dat[rep(seq_len(nrow(sampling_dat)), each=length(pred_mod)),]
    row.names(sampling_dat) <- NULL
    df_n <- data.frame(n,n_mod,pred_mod)
    df_n <- cbind(df_n,sampling_dat)

  }
  
  return(df_n)
}

extract_daily_training_testing_info <- function(i,list_param){
  #This function extracts training and testing information from the raster object produced for each tile
  #This is looping through tiles...
  ##Functions used
  #Defined outside this script:
  #pred_data_info_fun
  
  ##### Parse parameters and arguments ####
  
  raster_obj_file <- list_param$list_raster_obj_files[i]
  use_month <- list_param$use_month
  use_day <- list_param$use_day
  tile_id <- list_param$list_names_tile_id[i]
  
  ### Start script ##
  
  raster_obj <- load_obj(unlist(raster_obj_file)) #may not need unlist
  nb_models <- length((raster_obj$clim_method_mod_obj[[1]]$formulas))
  pred_mod <- paste("mod",c(1:nb_models,"_kr"),sep="")
  #we are assuming no monthly hold out...
  #we are assuming only one specific daily prop?
  nb_models <- length(pred_mod)
  #names(raster_obj$method_mod_obj[[1]])
  var_interp <- unique(raster_obj$tb_diagnostic_s$var_interp)
  method_interp <- unique(raster_obj$tb_diagnostic_s$method_interp)
  
  if(use_day==TRUE){
    
    list_data_day_v <- try(extract_list_from_list_obj(raster_obj$validation_mod_obj,"data_v"))
    list_data_day_s <- try(extract_list_from_list_obj(raster_obj$validation_mod_obj,"data_s"))
    sampling_dat_day <- extract_list_from_list_obj(raster_obj$method_mod_obj,"daily_dev_sampling_dat")
    #debug(pred_data_info_fun)
    #list_pred_data_day_s_info <- pred_data_info_fun(1,list_data=list_data_day_s,pred_mod=pred_mod,sampling_dat_info=sampling_dat_day)
    list_pred_data_day_s_info <- lapply(1:length(sampling_dat_day),FUN=pred_data_info_fun,
           list_data=list_data_day_s,pred_mod=pred_mod,sampling_dat_info=sampling_dat_day)
    list_pred_data_day_v_info <- lapply(1:length(sampling_dat_day),FUN=pred_data_info_fun,
           list_data=list_data_day_v,pred_mod=pred_mod,sampling_dat_info=sampling_dat_day)
    pred_data_day_s_info <- do.call(rbind,list_pred_data_day_s_info)
    pred_data_day_v_info <- do.call(rbind,list_pred_data_day_v_info)
    pred_data_day_s_info$training <- rep(1,nrow(pred_data_day_s_info)) 
    pred_data_day_v_info$training <- rep(0,nrow(pred_data_day_v_info)) 
    pred_data_day_info <-rbind(pred_data_day_v_info,pred_data_day_s_info)
    
    pred_data_day_info$method_interp <- rep(method_interp,nrow(pred_data_day_info)) 
    pred_data_day_info$var_interp <- rep(var_interp,nrow(pred_data_day_info)) 
    pred_data_day_info$tile_id <- rep(tile_id,nrow(pred_data_day_info))                                       
  }
  
  if(use_month==TRUE){
    
    list_data_month_s <- try(extract_list_from_list_obj(raster_obj$validation_mod_month_obj,"data_s"))
    list_data_month_v <- try(extract_list_from_list_obj(raster_obj$validation_mod_month_obj,"data_v"))
    sampling_dat_month <- extract_list_from_list_obj(raster_obj$clim_method_mod_obj,"sampling_month_dat")
    list_pred_data_month_s_info <- lapply(1:length(sampling_dat_month),FUN=pred_data_info_fun,
           list_data=list_data_month_s,pred_mod=pred_mod,sampling_dat_info=sampling_dat_month)
    list_pred_data_month_v_info <- mclapply(1:length(sampling_dat_month),FUN=pred_data_info_fun,
           list_data=list_data_month_v,pred_mod=pred_mod,sampling_dat_info=sampling_dat_month,mc.preschedule=FALSE,mc.cores = 1)
    #Note for list_pred_data_month_v_info it will be try-error when there is no holdout.
    #combine training and testing later? also combined with accuracy
    pred_data_month_s_info <- do.call(rbind,list_pred_data_month_s_info)
    #Adding a column for training: if data item used as training then it is 1
    pred_data_month_s_info$training <- try(rep(1,nrow(pred_data_month_s_info)))

    #Remove try-error??
    list_pred_data_month_v_info_tmp <- remove_from_list_fun(list_pred_data_month_v_info)
    list_pred_data_month_v_info <- list_pred_data_month_v_info_tmp$list
    if (length(list_pred_data_month_v_info)>0){
      pred_data_month_v_info <- do.call(rbind,list_pred_data_month_v_info)
      pred_data_month_v_info$training <- try(rep(0,nrow(pred_data_month_v_info)))
      pred_data_month_info <- rbind(pred_data_month_v_info,pred_data_month_s_info)
    }else{
      pred_data_month_info <- pred_data_month_s_info
    }
    
    pred_data_month_info$method_interp <- rep(method_interp,nrow(pred_data_month_info)) 
    pred_data_month_info$var_interp <- rep(var_interp,nrow(pred_data_month_info)) 
    pred_data_month_info$tile_id <- rep(tile_id,nrow(pred_data_month_info)) 
  }    
    
  if(use_month==FALSE){
    pred_data_month_info <- NULL
  }
  if(use_day==FALSE){
    pred_data_day_info <- NULL
  }  
  #prepare object to return
  pred_data_info_obj <- list(pred_data_month_info,pred_data_day_info)
  names(pred_data_info_obj) <- c("pred_data_month_info","pred_data_day_info")
  #could add data.frame data_s and data_v later...
  ###
  return(pred_data_info_obj)
}

remove_from_list_fun <- function(l_x,condition_class ="try-error"){
  index <- vector("list",length(l_x))
  for (i in 1:length(l_x)){
    if (inherits(l_x[[i]],condition_class)){
      index[[i]] <- FALSE #remove from list
    }else{
      index[[i]] <- TRUE
    }
  }
  l_x<-l_x[unlist(index)] #remove from list all elements using subset
  
  obj <- list(l_x,index)
  names(obj) <- c("list","valid")
  return(obj)
}

##Function to list predicted tif
list_tif_fun <- function(i,in_dir_list,pattern_str){
  #list.files(path=x,pattern=".*predicted_mod1_0_1.*20100101.*.tif",full.names=T)})
  pat_str<- pattern_str[i]
  list_tif_files_dates <-lapply(in_dir_list,
         FUN=function(x,pat_str){list.files(path=x,pattern=pat_str,full.names=T)},pat_str=pat_str)
  return(list_tif_files_dates)
} 

##############################
#### Parameters and constants  

#in_dir1 <- "/data/project/layers/commons/NEX_data/test_run1_03232014/output" #On Atlas
#in_dir1 <- "/nobackupp4/aguzman4/climateLayers/output10Deg/reg1/" #On NEX
in_dir1 <- "/nobackupp4/aguzman4/climateLayers/output20Deg/"

#/nobackupp4/aguzman4/climateLayers/output10Deg/reg1/finished.txt
#in_dir_list <- list.dirs(path=in_dir1,recursive=FALSE) #get the list of directories with resutls by 10x10 degree tiles
in_dir_list <- c(
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg2//-10.0_-70.0/",
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg4//40.0_0.0/",
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg4//50.0_0.0/",
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg6//60.0_40.0/",
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg6//30.0_40.0/",
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg8//40.0_130.0/")

#Models used.
#list_models<-c("y_var ~ s(lat,lon,k=4) + s(elev_s,k=3) + s(LST,k=3)",
#               "y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)",
#               "y_var ~ s(lat,lon,k=8) + s(elev_s,k=4) + s(LST,k=4)",
  
#use subset for now:

#in_dir_list <- c(
#"/nobackupp4/aguzman4/climateLayers/output10Deg/reg1/40.0_-120.0/",
#"/nobackupp4/aguzman4/climateLayers/output10Deg/reg1/35.0_-115.0/")
  
#in_dir_list <- file.path(in_dir1,read.table(file.path(in_dir1,"processed.txt"))$V1)
#in_dir_list <- as.list(in_dir_list[-1])
#in_dir_list <- in_dir_list[grep("bak",basename(basename(in_dir_list)),invert=TRUE)] #the first one is the in_dir1
#in_dir_shp <- in_dir_list[grep("shapefiles",basename(in_dir_list),invert=FALSE)] #select directory with shapefiles...
in_dir_shp <- c(
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg2/subset/shapefiles/",
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg4/subset/shapefiles/",
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg4/subset/shapefiles/",
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg6/subset/shapefiles/",
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg6/subset/shapefiles/",
"/nobackupp4/aguzman4/climateLayers/output20Deg/reg8/subset/shapefiles/")

#in_dir_shp <- "/nobackupp4/aguzman4/climateLayers/output10Deg/reg1/subset/shapefiles/"
#in_dir_shp <- "/nobackupp4/aguzman4/climateLayers/output20Deg/reg2/subset/shapefiles"
in_dir_shp_list <- list.files(in_dir_shp,".shp",full.names=T)

#in_dir_list <- in_dir_list[grep("shapefiles",basename(in_dir_list),invert=TRUE)] 
#the first one is the in_dir1
# the last directory contains shapefiles 
y_var_name <- "dailyTmax"
interpolation_method <- c("gam_CAI")
out_prefix<-"run5_global_analyses_08252014"

#out_dir<-"/data/project/layers/commons/NEX_data/" #On NCEAS Atlas
out_dir <- "/nobackup/bparmen1/" #on NEX
#out_dir <-paste(out_dir,"_",out_prefix,sep="")
create_out_dir_param <- TRUE

#system("ls /nobackup/bparmen1")

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_prefix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

setwd(out_dir)
                                   
CRS_locs_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

day_to_mosaic <- c("20100101","20100901")
file_format <- ".tif" #format for mosaiced files
NA_flag_val <- -9999  #No data value

#day_to_mosaic <- NULL #if day to mosaic is null then mosaic all dates

##raster_prediction object : contains testing and training stations with RMSE and model object

#l_shp <- lapply(1:length(in_dir_shp_list),FUN=function(i){paste(strsplit(in_dir_shp_list[i],"_")[[1]][2:3],collapse="_")})
#match(l_shp,in_dir_list)
#in_dir_list[match(in_dir_list,l_shp]

list_raster_obj_files <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="^raster_prediction_obj.*.RData",full.names=T)})
basename(dirname(list_raster_obj_files[[1]]))
list_names_tile_coord <- lapply(list_raster_obj_files,FUN=function(x){basename(dirname(x))})
list_names_tile_id <- paste("tile",1:length(list_raster_obj_files),sep="_")
names(list_raster_obj_files)<- list_names_tile_id

lf_covar_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar_obj.*.RData",full.names=T)})
lf_covar_tif <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar.*.tif",full.names=T)})
#diagnostics_obj_gam_fitting_dailyTmax7__08062014.RData
lf_diagnostic_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="diagnostics_.*.RData",full.names=T)})
lf_diagnostic_obj <- lf_diagnostic_obj[grep("lk_min",lf_diagnostic_obj,invert=T)] #remove object that have lk_min...

lf_validation_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="gam_CAI_validation_mod_obj_dailyTmax.*.RData",full.names=T)})
#validation_mod_obj <-load_obj("/nobackupp4/aguzman4/climateLayers/output20Deg/reg6/60.0_40.0/gam_CAI_validation_mod_obj_dailyTmax60.0_40.0.RData")
debug(extract_from_list_obj)
tb_diagnostic_v<-extract_from_list_obj(validation_mod_obj,"metrics_v") 

########################## START SCRIPT ##############################

################################################################
######## PART 1: Generate tables to collect information 
######## over all tiles in North America 

##Function to collect all the tables from tiles into a table
###Table 1: Average accuracy metrics
###Table 2: daily accuracy metrics for all tiles

#First create table of tiles under analysis and their coord
df_tile_processed <- data.frame(tile_coord=basename(in_dir_list))
df_tile_processed$tile_id <- unlist(list_names_tile_id) #Arbitrary tiling number!!
df_tile_processed$path_NEX <- in_dir_list
  
##Quick exploration of raster object
robj1 <- load_obj(list_raster_obj_files[[2]]) #This is tile corresponding to Oregon

names(robj1)
names(robj1$method_mod_obj[[1]]) #for January 1, 2010
names(robj1$method_mod_obj[[1]]$dailyTmax) #for January

names(robj1$clim_method_mod_obj[[1]]$data_month) #for January
names(robj1$validation_mod_month_obj[[1]]$data_s) #for January with predictions
#Get the number of models predicted
nb_mod <- length(unique(robj1$tb_diagnostic_v$pred_mod))

validation_mod_obj <-load_obj("/nobackupp4/aguzman4/climateLayers/output20Deg/reg6/60.0_40.0/gam_CAI_validation_mod_obj_dailyTmax60.0_40.0.RData")
tb_diagnostic_v<-extract_from_list_obj(validation_mod_obj,"metrics_v") 


#clim_method_mod_obj <- load_obj("/nobackupp4/aguzman4/climateLayers/output20Deg/reg6/60.0_40.0/gam_CAI_mod_dailyTmax60.0_40.0.RData")
#list_data_v <- extract_list_from_list_obj(clim_method_mod_obj,"data_month_v") #extract monthly testing/validation dataset
#list_data_s <- extract_list_from_list_obj(clim_method_mod_obj,"data_month") #extract monthly training/fitting dataset
#rast_day_yearlist <- extract_list_from_list_obj(clim_method_mod_obj,"clim") #list_tmp #list of predicted images over full year at monthly time scale
#list_sampling_dat <- extract_list_from_list_obj(clim_method_mod_obj,"sampling_month_dat")

list_tb_diagnostic_v <- mclapply(lf_validation_obj,FUN=function(x){try( x<- load_obj(x)); try(extract_from_list_obj(x,"metrics_v"))},mc.preschedule=FALSE,mc.cores = 6)                           
names(list_tb_diagnostic_v) <- list_names_tile_id

#undebug(extract_from_list_obj)
#validation_mod_month_obj <-load_obj("/nobackupp4/aguzman4/climateLayers/output20Deg/reg6/60.0_40.0/gam_CAI_validation_mod_month_obj_dailyTmax60.0_40.0.RData")
#tb_diagnostic_v<-extract_from_list_obj(validation_mod_obj,"metrics_v") 


################
#### Table 1: Average accuracy metrics

#can use a maximum of 6 cores on the NEX Bridge
#summary_metrics_v_list <- mclapply(list_raster_obj_files[5:6],FUN=function(x){try( x<- load_obj(x)); try(x[["summary_metrics_v"]]$avg)},mc.preschedule=FALSE,mc.cores = 2)                           

summary_metrics_v_list <- mclapply(list_raster_obj_files,FUN=function(x){try( x<- load_obj(x)); try(x[["summary_metrics_v"]]$avg)},mc.preschedule=FALSE,mc.cores = 6)                           
names(summary_metrics_v_list) <- list_names_tile_id

summary_metrics_v_tmp <- remove_from_list_fun(summary_metrics_v_list)$list
df_tile_processed$metrics_v <- as.integer(remove_from_list_fun(summary_metrics_v_list)$valid)
#Now remove "try-error" from list of accuracy)

summary_metrics_v_NA <- do.call(rbind.fill,summary_metrics_v_tmp) #create a df for NA tiles with all accuracy metrics
#tile_coord <- lapply(1:length(summary_metrics_v_list),FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=summary_metrics_v_list)
#add the tile id identifier
tile_id_tmp <- lapply(1:length(summary_metrics_v_tmp),
                     FUN=function(i,x,y){rep(y[i],nrow(x[[i]]))},x=summary_metrics_v_tmp,y=names(summary_metrics_v_tmp))
#adding tile id summary data.frame
summary_metrics_v_NA$tile_id <-unlist(tile_id_tmp)
summary_metrics_v_NA$n <- as.integer(summary_metrics_v_NA$n)

summary_metrics_v_NA <- merge(summary_metrics_v_NA,df_tile_processed[,1:2],by="tile_id")

tx<-strsplit(as.character(summary_metrics_v_NA$tile_coord),"_")
lat<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][1]},x=tx))
long<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][2]},x=tx))
summary_metrics_v_NA$lat <- lat
summary_metrics_v_NA$lon <- long

write.table(as.data.frame(summary_metrics_v_NA),
            file=file.path(out_dir,paste("summary_metrics_v2_NA_",out_prefix,".txt",sep="")),sep=",")

#################
###Table 2: daily accuracy metrics for all tiles
#this takes about 25min
#tb_diagnostic_v_list <- lapply(list_raster_obj_files,FUN=function(x){x<-load_obj(x);x[["tb_diagnostic_v"]]})                           
tb_diagnostic_v_list <- mclapply(list_raster_obj_files,FUN=function(x){try(x<-load_obj(x));try(x[["tb_diagnostic_v"]])},mc.preschedule=FALSE,mc.cores = 6)                           

names(tb_diagnostic_v_list) <- list_names_tile_id
tb_diagnostic_v_tmp <- remove_from_list_fun(tb_diagnostic_v_list)$list
#df_tile_processed$tb_diag <- remove_from_list_fun(tb_diagnostic_v_list)$valid

tb_diagnostic_v_NA <- do.call(rbind.fill,tb_diagnostic_v_tmp) #create a df for NA tiles with all accuracy metrics
tile_id_tmp <- lapply(1:length(tb_diagnostic_v_tmp),
                     FUN=function(i,x,y){rep(y[i],nrow(x[[i]]))},x=tb_diagnostic_v_tmp,y=names(tb_diagnostic_v_tmp))

tb_diagnostic_v_NA$tile_id <- unlist(tile_id_tmp) #adding identifier for tile

tb_diagnostic_v_NA <- merge(tb_diagnostic_v_NA,df_tile_processed[,1:2],by="tile_id")

write.table((tb_diagnostic_v_NA),
            file=file.path(out_dir,paste("tb_diagnostic_v_NA","_",out_prefix,".txt",sep="")),sep=",")

####### process gam fitting diagnostic info

#/nobackupp4/aguzman4/climateLayers/output20Deg/reg5/20.0_30.0//diagnostics_obj_gam_fitting_TMax_9_mod2_08062014.RData
#lf_diagnostic_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="diagnostics_.*.RData",full.names=T)})
#lf_diagnostic_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="diagnostics_obj_gam_fitting_TMax_*_mod*_08062014.RData",full.names=T)})
lf_diagnostic_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="diagnostics_obj_gam_fitting_TMax_.*._mod.*._08062014.RData",full.names=T)})

#lf_diagnostic_obj <- lf_diagnostic_obj[grep("lk_min",lf_diagnostic_obj,invert=T)] #remove object that have lk_min...

names(lf_diagnostic_obj) <- list_names_tile_id
lf_diagnostic_obj_tmp <- remove_from_list_fun(lf_diagnostic_obj)$list
#df_tile_processed$tb_diag <- remove_from_list_fun(tb_diagnostic_v_list)$valid

gam_diagnostic_tb_list <- vector("list",length=length(lf_diagnostic_obj_tmp))
for(i in 1:length(lf_diagnostic_obj_tmp)){
  l_diagnostic_obj_tmp <- lf_diagnostic_obj_tmp[[i]]
  tile_id_name <-  names(lf_diagnostic_obj_tmp)[i]
  #l_diagnostic_obj_tmp <- l_diagnostic_obj_tmp[grep("lk_min",l_diagnostic_obj_tmp,invert=T)] #remove object that have lk_min...
  l_diagnostic_obj_tmp_list <- lapply(l_diagnostic_obj_tmp,FUN=function(x){try(x<-load_obj(x));try(x[["df_diagnostics"]])})#,mc.preschedule=FALSE,mc.cores = 6)                            
  gam_diagnostic_tb <- do.call(rbind.fill,l_diagnostic_obj_tmp_list)#create a df for NA tiles with all accuracy metrics
  gam_diagnostic_tb$tile_id <- tile_id_name
  gam_diagnostic_tb_list[[i]] <- gam_diagnostic_tb    
}

gam_diagnostic_df <- do.call(rbind.fill,gam_diagnostic_tb_list) #create a df for NA tiles with all accuracy metrics

write.table(gam_diagnostic_df,
            file=file.path(out_dir,paste("gam_diagnostic_df_",out_prefix,".txt",sep="")),sep=",")


#Now look at the 100 tiles of 10x10
#lf_test<-list.files("/nobackupp4/aguzman4/climateLayers/output10Deg/*/*/","diagnostics_obj_gam_fitting*")  
lf_test <-list.files("/nobackupp4/aguzman4/climateLayers/output10Deg/","diagnostics_obj_gam_fitting.*.RData",recursive=T,full.names=T)

gam_diagnostic_10x10tb_list <- vector("list",length=length(lf_test))
lf_diagnostic_obj_tmp <- lf_test  
for(i in 1:length( lf_diagnostic_obj_tmp)){
  l_diagnostic_obj_tmp <-  lf_diagnostic_obj_tmp[[i]]
  tile_coord <-  basename(dirname(lf_diagnostic_obj_tmp[i]))
  #l_diagnostic_obj_tmp <- l_diagnostic_obj_tmp[grep("lk_min",l_diagnostic_obj_tmp,invert=T)] #remove object that have lk_min...
  l_diagnostic_obj_tmp_list <- lapply(l_diagnostic_obj_tmp,FUN=function(x){try(x<-load_obj(x));try(x[["df_diagnostics"]])})#,mc.preschedule=FALSE,mc.cores = 6)                            
  gam_diagnostic_tb <- do.call(rbind.fill,l_diagnostic_obj_tmp_list)#create a df for NA tiles with all accuracy metrics
  gam_diagnostic_tb$tile_coord <- tile_coord
  gam_diagnostic_10x10tb_list[[i]] <- gam_diagnostic_tb    
}

gam_diagnostic_10x10_df <- do.call(rbind.fill,gam_diagnostic_10x10tb_list) #create a df for NA tiles with all accuracy metrics

list_tile_coord <- unique(gam_diagnostic_10x10_df$tile_coord)
list_tile_id <- paste("tile_",1:length(list_tile_coord),sep="")

tile_id_df <- data.frame(tile_coord=list_tile_coord,tile_id=list_tile_id)
gam_diagnostic_10x10_df <- merge(gam_diagnostic_10x10_df,tile_id_df,by="tile_coord")

# write.table(gam_diagnostic_10x10_df,
#             file=file.path(out_dir,paste("gam_diagnostic_10x10_df_",out_prefix,".txt",sep="")),sep=",")

#################
###Table 3: monthly station information with predictions for all tiles

#load data_month for specific tiles
# data_month <- extract_from_list_obj(robj1$clim_method_mod_obj,"data_month")
# names(data_month) #this contains LST means (mm_1, mm_2 etc.) as well as TMax and other info
# 
# data_month_s_list <- mclapply(list_raster_obj_files,FUN=function(x){try(x<-load_obj(x));try(x$validation_mod_month_obj[["data_s"]])},mc.preschedule=FALSE,mc.cores = 6)                           
# 
# names(data_month_s_list) <- list_names_tile_id
# 
# data_month_tmp <- remove_from_list_fun(data_month_s_list)$list
# #df_tile_processed$metrics_v <- remove_from_list_fun(data_month_s_list)$valid
# 
# tile_id <- lapply(1:length(data_month_tmp),
#                   FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=data_month_tmp)
# data_month_NAM <- do.call(rbind.fill,data_month_list) #combined data_month for "NAM" North America
# data_month_NAM$tile_id <- unlist(tile_id)
# 
# write.table((data_month_NAM),
#             file=file.path(out_dir,paste("data_month_s_NAM","_",out_prefix,".txt",sep="")),sep=",")

##### Table 4: Add later on: daily info
### with also data_s and data_v saved!!!

#Insert here...compute input and predicted ranges to spot potential errors?

######################################################
####### PART 4: Get shapefile tiling with centroids ###

#get shape files for the region being assessed:

list_shp_world <- list.files(path=in_dir_shp,pattern=".*.shp",full.names=T)
l_shp <- unlist(lapply(1:length(list_shp_world),FUN=function(i){paste(strsplit(list_shp_world[i],"_")[[1]][2:3],collapse="_")}))
matching_index <- match(basename(in_dir_list),l_shp)
list_shp_reg_files <- list_shp_world[matching_index]
df_tile_processed$shp_files <-list_shp_world[matching_index]

tx<-strsplit(as.character(df_tile_processed$tile_coord),"_")
lat<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][1]},x=tx))
long<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][2]},x=tx))
df_tile_processed$lat <- lat
df_tile_processed$lon <- long

#put that list in the df_processed and also the centroids!!
write.table(df_tile_processed,
            file=file.path(out_dir,paste("df_tile_processed_",out_prefix,".txt",sep="")),sep=",")

######################################################
####### PART 2 CREATE MOSAIC OF PREDICTIONS PER DAY, Delta surfaces and clim ###

dates_l <- unique(robj1$tb_diagnostic_s$date) #list of dates to query tif

## make this a function? report on number of tiles used for mosaic...

#inputs: build a pattern to find files
y_var_name <- "dailyTmax" #set up in parameters of this script
interpolation_method <- c("gam_CAI") #set up in parameters of the script
name_method <- paste(interpolation_method,"_",y_var_name,"_",sep="")
#make it general using nb_mod!!
#could be set up at the begining?

mod_id <- c(1:(nb_mod-1),"_kr")
pred_pattern_str <- paste(".*predicted_mod",mod_id,"_0_1.*",sep="")
#,".*predicted_mod2_0_1.*",".*predicted_mod3_0_1.*",".*predicted_mod_kr_0_1.*")
#l_pattern_models <- lapply(c(".*predicted_mod1_0_1.*",".*predicted_mod2_0_1.*",".*predicted_mod3_0_1.*",".*predicted_mod_kr_0_1.*"),
#                           FUN=function(x){paste(x,dates_l,".*.tif",sep="")})
l_pattern_models <- lapply(pred_pattern_str,
                           FUN=function(x){paste(x,dates_l,".*.tif",sep="")})
#gam_CAI_dailyTmax_predicted_mod_kr_0_1_20101231_30_145.0_-120.0.tif
#gam_CAI_dailyTmax_predicted_mod_kr_0_1_20101231_30_145.0_-120.0.tif                    

##Get list of predicted tif across all tiles, models and dates...
#this takes time, use mclapply!!
lf_pred_tif <- vector("list",length=length(l_pattern_models)) #number of models is 3
for (i in 1:length(l_pattern_models)){
  l_pattern_mod <- l_pattern_models[[i]] #365 dates
  #list_tif_files_dates <-lapply(1:length(l_pattern_mod),FUN=list_tif_fun, 
  #                            in_dir_list=in_dir_list,pattern_str=l_pattern_models[[i]])
  list_tif_files_dates <-mclapply(1:length(l_pattern_mod),FUN=list_tif_fun, 
                              in_dir_list=in_dir_list,pattern_str=l_pattern_models[[i]],mc.preschedule=FALSE,mc.cores = 6)
  
  lf_pred_tif[[i]] <- list_tif_files_dates
}

#Now get the clim surfaces:
month_l <- paste("clim_month_",1:12,sep="")
l_pattern_models <- lapply(c("_mod1_0_1.*","_mod2_0_1.*","_mod3_0_1.*","_mod_kr_0_1.*"),
                           FUN=function(x){paste("*.",month_l,x,".*.tif",sep="")})
#"CAI_TMAX_clim_month_11_mod2_0_145.0_-120.0.tif"
lf_clim_tif <- vector("list",length=nb_mod) #number of models is 3
for (i in 1:length(l_pattern_models)){
  l_pattern_mod <- l_pattern_models[[i]] #12 dates
  list_tif_files_dates <- mclapply(1:length(l_pattern_mod),FUN=list_tif_fun, 
                              in_dir_list=in_dir_list,pattern_str=l_pattern_models[[i]],mc.preschedule=FALSE,mc.cores = 6)
  lf_clim_tif[[i]] <- list_tif_files_dates
}

#Now get delta surfaces:
date_l# <- paste("clim_month_",1:12,sep="")
#l_pattern_models <- lapply(c("_mod1_0_1.*","_mod2_0_1.*","_mod3_0_1.*","_mod_kr_0_1.*"),
#                           FUN=function(x){paste("*.",month_l,x,".*.tif",sep="")})
l_pattern_models <- lapply(c(".*delta_dailyTmax_mod1_del_0_1.*",".*delta_dailyTmax_mod2_del_0_1.*",".*delta_dailyTmax_mod3_del_0_1.*",".*delta_dailyTmax_mod_kr_del_0_1.*"),
                           FUN=function(x){paste(x,dates_l,".*.tif",sep="")})

lf_delta_tif <- vector("list",length=nb_mod) #number of models is 3
for (i in 1:length(l_pattern_models)){
  l_pattern_mod <- l_pattern_models[[i]]
  list_tif_files_dates <- mclapply(1:length(l_pattern_mod),FUN=list_tif_fun, 
                              in_dir_list=in_dir_list,pattern_str=l_pattern_models[[i]],mc.preschedule=FALSE,mc.cores = 6)
  lf_delta_tif[[i]] <- list_tif_files_dates
}


#### NOW create mosaic images for daily prediction

out_prefix_s <- paste(name_method,c("predicted_mod1_0_01","predicted_mod2_0_01","predicted_mod3_0_01","predicted_mod_kr_0_1"),sep="")
dates_l #list of predicted dates
#l_out_rastnames_var <- paste(name_method,"predicted_mod1_0_01_",dates_l,sep="")
l_out_rastnames_var <- lapply(out_prefix_s,
                              FUN=function(x){paste(x,"_",dates_l,sep="")})

#nb_mod <- 4 #this is set up earlier
##Add option to specify wich dates to mosaic??
day_to_mosaic <- c("20100101","20100901")
if (!is.null(day_to_mosaic)){
  list_days <-match(day_to_mosaic,dates_l)
}else{
  list_days <- 1:365 #should check for year in case it has 366, add later!!
}
###Make this a function later??
for (i in 1:nb_mod){
  
  list_tif_files_dates <- lf_pred_tif[[i]] 

  mosaic_list_var <- list_tif_files_dates  
  out_rastnames_var <- l_out_rastnames_var[[i]]

  file_format <- ".tif"
  NA_flag_val <- -9999

  j<-1 #date index for loop
  list_param_mosaic<-list(j,mosaic_list_var,out_rastnames_var,out_dir,file_format,NA_flag_val)
  names(list_param_mosaic)<-c("j","mosaic_list","out_rastnames","out_path","file_format","NA_flag_val")
  #list_var_mosaiced <- mclapply(1:2,FUN=mosaic_m_raster_list,list_param=list_param_mosaic,mc.preschedule=FALSE,mc.cores = 2)
  list_var_mosaiced <- mclapply(list_days,FUN=mosaic_m_raster_list,list_param=list_param_mosaic,mc.preschedule=FALSE,mc.cores = 2)
  #list_var_mosaiced <- mclapply(1:1,FUN=mosaic_m_raster_list,list_param=list_param_mosaic,mc.preschedule=FALSE,mc.cores = 1)
  #list_var_mosaiced <- mclapply(1:365,FUN=mosaic_m_raster_list,list_param=list_param_mosaic,mc.preschedule=FALSE,mc.cores = 2)
  
  #mosaic for delt sufaces?
  #mosoaic for clim months?
  
}

######################
### mosaic clim monthly data...this will be a function later...

#Now get the clim surfaces:
month_l <- paste("clim_month_",1:12,sep="")
#l_pattern_models <- lapply(c("_mod1_0_1","_mod2_0_1","_mod3_0_1","_mod_kr_0_1"),
#                           FUN=function(x){paste(x,"_",month_l,sep="")})

out_prefix_s <- paste(name_method,c("_mod1_0_01","_mod2_0_01","_mod3_0_01","_mod_kr_0_1"),sep="")
dates_l #list of predicted dates
#l_out_rastnames_var <- paste(name_method,"predicted_mod1_0_01_",dates_l,sep="")
l_out_rastnames_var <- lapply(out_prefix_s,
                              FUN=function(x){paste(x,"_",month_l,sep="")})

for (i in 1:nb_mod){
  
  #this should be the input param for the new function generated automatically...
  list_tif_files_dates <- lf_clim_tif[[i]] 
  mosaic_list_var <- list_tif_files_dates  
  out_rastnames_var <- l_out_rastnames_var[[i]]
  #file_format <- ".tif"
  #NA_flag_val <- -9999

  j<-1 #date index for loop
  list_param_mosaic<-list(j,mosaic_list_var,out_rastnames_var,out_dir,file_format,NA_flag_val)
  names(list_param_mosaic)<-c("j","mosaic_list","out_rastnames","out_path","file_format","NA_flag_val")
  #list_var_mosaiced <- mclapply(1:2,FUN=mosaic_m_raster_list,list_param=list_param_mosaic,mc.preschedule=FALSE,mc.cores = 2)
  list_var_mosaiced <- mclapply(1:12,FUN=mosaic_m_raster_list,list_param=list_param_mosaic,mc.preschedule=FALSE,mc.cores = 4)
}

######################
#### NOW create mosaic images for daily delta prediction
#This should be a function!!!
date_l# <- paste("clim_month_",1:12,sep="")
#l_pattern_models <- lapply(c("_mod1_0_1.*","_mod2_0_1.*","_mod3_0_1.*","_mod_kr_0_1.*"),
#                           FUN=function(x){paste("*.",month_l,x,".*.tif",sep="")})
#l_pattern_models <- lapply(c(".*delta_dailyTmax_mod1_del_0_1.*",".*delta_dailyTmax_mod2_del_0_1.*",".*delta_dailyTmax_mod3_del_0_1.*",".*delta_dailyTmax_mod_kr_del_0_1.*"),
#                           FUN=function(x){paste(x,dates_l,".*.tif",sep="")})

out_prefix_s <- paste(name_method,c("delta_mod1_0_01","delta_mod2_0_01","delta_mod3_0_01","delta_mod_kr_0_1"),sep="")
dates_l #list of predicted dates
#l_out_rastnames_var <- paste(name_method,"predicted_mod1_0_01_",dates_l,sep="")
l_out_rastnames_var <- lapply(out_prefix_s,
                              FUN=function(x){paste(x,"_",dates_l,sep="")})

#nb_mod <- 4 #this is set up earlier
##Add option to specify wich dates to mosaic??
day_to_mosaic <- c("20100101","20100901")
if (!is.null(day_to_mosaic)){
  list_days <-match(day_to_mosaic,dates_l)
}else{
  list_days <- 1:365 #should check for year in case it has 366, add later!!
}
###Make this a function later??
for (i in 1:nb_mod){
  
  list_tif_files_dates <- lf_pred_tif[[i]] 
  mosaic_list_var <- list_tif_files_dates  
  out_rastnames_var <- l_out_rastnames_var[[i]]
  #this is be set up earlier...
  #file_format <- ".tif"
  #NA_flag_val <- -9999

  j<-1 #date index for loop
  list_param_mosaic<-list(j,mosaic_list_var,out_rastnames_var,out_dir,file_format,NA_flag_val)
  names(list_param_mosaic)<-c("j","mosaic_list","out_rastnames","out_path","file_format","NA_flag_val")
  #list_var_mosaiced <- mclapply(1:2,FUN=mosaic_m_raster_list,list_param=list_param_mosaic,mc.preschedule=FALSE,mc.cores = 2)
  list_var_mosaiced <- mclapply(list_days,FUN=mosaic_m_raster_list,list_param=list_param_mosaic,mc.preschedule=FALSE,mc.cores = 2)
}

### Now find out how many files were predicted


######################################################
####### PART 3: EXAMINE STATIONS AND MODEL FITTING ###

### Stations and model fitting ###
#summarize location and number of training and testing used by tiles

names(robj1$clim_method_mod_obj[[1]]$data_month) # monthly data for January
#names(robj1$validation_mod_month_obj[[1]]$data_s) # daily for January with predictions
#note that there is no holdout in the current run at the monthly time scale:

robj1$clim_method_mod_obj[[1]]$data_month_v #zero rows for testing stations at monthly timescale
#load data_month for specific tiles
data_month <- extract_from_list_obj(robj1$clim_method_mod_obj,"data_month")

names(data_month) #this contains LST means (mm_1, mm_2 etc.) as well as TMax and other info

#problem with tile 12...the raster ojbect has missing sub object
#data_month_list <- lapply(1:length(list_raster_obj_files),x=list_raster_obj_files,
#                          FUN=function(i,x){x<-load_obj(x[[i]]);
#                                            extract_from_list_obj(x$validation_mod_month_obj,"data_s")})                           

### make this part a function:

#create a table for every month, day and tiles...
# data_month_list <- lapply(1:length(list_raster_obj_files),x=list_raster_obj_files,
#                           FUN=function(i,x){x<-load_obj(x[[i]]);
#                                             extract_from_list_obj(x$clim_method_mod_obj,"data_month")})                           
# 
# names(data_month_list) <- paste("tile","_",1:length(data_month_list),sep="")
# 
# #names(data_month_list) <- basename(in_dir_list) #use folder id instead
# 
# list_names_tile_id <- paste("tile",1:length(list_raster_obj_files),sep="_")
# 
# #tile_id <- lapply(1:length(data_month_list),
# #                  FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=data_month_list)
# 
# data_month_NAM <- do.call(rbind.fill,data_month_list) #combined data_month for "NAM" North America
# data_month_NAM$tile_id <- unlist(tile_id)
# 
# names(robj1$validation_mod_day_obj[[1]]$data_s) # daily for January with predictions
# dim(robj1$validation_mod_month_obj[[1]]$data_s) # daily for January with predictions
# 

use_day=TRUE
use_month=TRUE
 
#list_raster_obj_files <- c("/data/project/layers/commons/NEX_data/output_run3_global_analyses_06192014/output10Deg/reg1//30.0_-100.0/raster_prediction_obj_gam_CAI_dailyTmax30.0_-100.0.RData",
#                    "/data/project/layers/commons/NEX_data/output_run3_global_analyses_06192014/output10Deg/reg1//30.0_-105.0/raster_prediction_obj_gam_CAI_dailyTmax30.0_-105.0.RData")

list_names_tile_id <- df_tile_processed$tile_id
list_raster_obj_files[list_names_tile_id]
#list_names_tile_id <- c("tile_1","tile_2")
list_param_training_testing_info <- list(list_raster_obj_files[list_names_tile_id],use_month,use_day,list_names_tile_id)
names(list_param_training_testing_info) <- c("list_raster_obj_files","use_month","use_day","list_names_tile_id")
 
list_param <- list_param_training_testing_info
#debug(extract_daily_training_testing_info)
#pred_data_info <- extract_daily_training_testing_info(1,list_param=list_param_training_testing_info)
pred_data_info <- mclapply(1:length(list_raster_obj_files[list_names_tile_id]),FUN=extract_daily_training_testing_info,list_param=list_param_training_testing_info,mc.preschedule=FALSE,mc.cores = 6)
#pred_data_info <- mclapply(1:length(list_raster_obj_files[list_names_tile_id][1:6]),FUN=extract_daily_training_testing_info,list_param=list_param_training_testing_info,mc.preschedule=FALSE,mc.cores = 6)
#pred_data_info <- lapply(1:length(list_raster_obj_files),FUN=extract_daily_training_testing_info,list_param=list_param_training_testing_info)
#pred_data_info <- lapply(1:length(list_raster_obj_files[1]),FUN=extract_daily_training_testing_info,list_param=list_param_training_testing_info)

pred_data_info_tmp <- remove_from_list_fun(pred_data_info)$list #remove data not predicted
##Add tile nanmes?? it is alreaready there
#names(pred_data_info)<-list_names_tile_id
pred_data_month_info <- do.call(rbind,lapply(pred_data_info_tmp,function(x){x$pred_data_month_info}))
pred_data_day_info <- do.call(rbind,lapply(pred_data_info_tmp,function(x){x$pred_data_day_info}))

#putput inforamtion in csv !!
write.table(pred_data_month_info,
            file=file.path(out_dir,paste("pred_data_month_info_",out_prefix,".txt",sep="")),sep=",")
write.table(pred_data_day_info,
            file=file.path(out_dir,paste("pred_data_day_info_",out_prefix,".txt",sep="")),sep=",")

########### LAST PART: COPY SOME DATA BACK TO ATLAS #####

### This assumes the tree structure has been replicated on Atlas:
#for i in 1:length(df_tiled_processed$tile_coord)
#output_atlas_dir <- "/data/project/layers/commons/NEX_data/output_run3_global_analyses_06192014/output10Deg/reg1"
output_atlas_dir <- "/data/project/layers/commons/NEX_data/output_run4_global_analyses_08142014/output20Deg"

#Make directories on ATLAS
#for (i in 1:length(df_tiled_processed$tile_coord)){
#  create_dir_fun(file.path(output_atlas_dir,as.character(df_tiled_processed$tile_coord[i])),out_suffix=NULL)
#}  

#Make directories on ATLAS for shapefiles
#for (i in 1:length(df_tiled_processed$tile_coord)){
#  create_dir_fun(file.path(output_atlas_dir,as.character(df_tiled_processed$tile_coord[i]),"/shapefiles"),out_suffix=NULL)
#}  


#Copy summary textfiles and mosaic back to atlas

Atlas_dir <- file.path("/data/project/layers/commons/NEX_data/",basename(out_dir))#,"output/subset/shapefiles")
Atlas_hostname <- "parmentier@atlas.nceas.ucsb.edu"
lf_cp_f <- list.files(out_dir,full.names=T)#copy all files can filter later
filenames_NEX <- paste(lf_cp_f,collapse=" ")  #copy raster prediction object
cmd_str <- paste("scp -p",filenames_NEX,paste(Atlas_hostname,Atlas_dir,sep=":"), sep=" ")
system(cmd_str)

#system("scp -p ./*.txt parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/output_run2_global_analyses_05122014")
#system("scp -p ./*.txt ./*.tif parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/output_run2_global_analyses_05122014")

#### COPY SHAPEFILES, TIF MOSAIC, COMBINED TEXT FILES etc...

#copy shapefiles defining regions
Atlas_dir <- file.path("/data/project/layers/commons/NEX_data/",basename(out_dir),"output/subset/shapefiles")
Atlas_hostname <- "parmentier@atlas.nceas.ucsb.edu"
lf_cp_shp <- df_tile_processed$shp_files #get all the files...

#lf_cp_shp <- list.files(in_dir_shp, ".shp",full.names=T)
list_tile_scp <- 1:8

for (j in 1:length(list_tile_scp)){
  tile_nb <- list_tile_scp[j]
  
  in_dir_tile <-dirname(df_tile_processed$shp_files[tile_nb])
  #/data/project/layers/commons/NEX_data/output_run2_05122014/output
  #output_atlas_dir
  #Atlas_dir <- file.path(file.path("/data/project/layers/commons/NEX_data/",basename(out_dir),"output"),in_dir_tile)
  Atlas_dir <- file.path(output_atlas_dir,as.character(df_tile_processed$tile_coord[j]),"/shapefiles")

  Atlas_hostname <- "parmentier@atlas.nceas.ucsb.edu"

  filenames_NEX <- paste(lf_cp_shp,collapse=" ")  #copy raster prediction object
  cmd_str <- paste("scp -p",filenames_NEX,paste(Atlas_hostname,Atlas_dir,sep=":"), sep=" ")
  system(cmd_str)
}

#### FIRST COPY DATA FOR SPECIFIC TILES #####
#Copy specific tiles info back...This assumes that the tree structre 
#has been created on ATLAS:
#../$out_dir/ouput/tile_coord

#list_tile_scp <- c(1,2)
list_tile_scp <- 1:8

for (j in 1:length(list_tile_scp)){
  tile_nb <- list_tile_scp[j]
  #nb_mod <- 3+1 #set up earlier
  date_selected <- c("20100101","20100901") #should be set up earlier
  date_index <- c(1,244) #list_day??
  #tile_nb <- 1

  in_dir_tile <- basename(df_tile_processed$path_NEX[tile_nb])
  #/data/project/layers/commons/NEX_data/output_run2_05122014/output
  #output_atlas_dir
  #Atlas_dir <- file.path(file.path("/data/project/layers/commons/NEX_data/",basename(out_dir),"output"),in_dir_tile)
  Atlas_dir <- file.path(output_atlas_dir,in_dir_tile)
  Atlas_hostname <- "parmentier@atlas.nceas.ucsb.edu"
  #filenames_NEX <- list_raster_obj_files[tile_nb] #copy raster prediction object
  #cmd_str <- paste("scp -p",filenames_NEX,paste(Atlas_hostname,Atlas_dir,sep=":"), sep=" ")
  #system(cmd_str)

  #Now copy back tif for specific dates and tile (date 1 and date 244)
  #nb_mod <- 3+1
  lf_cp_day <- vector("list",length=length(date_selected))
  #Get relevant daily info
  for(i in 1:length(date_selected)){
    #d
    index <- date_index[i]  
    #get all predicted tmax files for all models and specific date, tile
    lf_cp_pred_tif  <- unlist(lapply(1:nb_mod,FUN=function(x){lf_pred_tif[[x]][[index]][[tile_nb]]}))
    lf_cp_delta_tif <- unlist(lapply(1:nb_mod,FUN=function(x){lf_delta_tif[[x]][[index]][[tile_nb]]}))
    lf_cp_day[[i]] <- unlist(c(unlist(lf_cp_pred_tif),unlist(lf_cp_delta_tif)))
  }
  #get the monthly info...
  month_index <- 1:12 #can subset
  #month_index <- c(1,9) #can subset
  lf_cp_month <- vector("list",length=length(month_index))

  for(i in 1:length(month_index)){
    #d
    index <- month_index[i]  
    #get all predicted tmax files for all models and specific date, tile
    lf_cp_month[[i]]  <- unlist(lapply(1:nb_mod,FUN=function(x){lf_clim_tif[[x]][[index]][[tile_nb]]}))
  }
  ##Add RData object for specified tile...
  lf_cp_RData_tif <- c(lf_covar_obj[tile_nb],lf_covar_tif[tile_nb],list_raster_obj_files[[tile_nb]],lf_diagnostic_obj[[tile_nb]])
  #unlist(lf_cp_RData_tif)
  lf_cp <- unlist(c(lf_cp_day,lf_cp_month,lf_cp_RData_tif))
  #lf_cp <- c(unlist(c(lf_cp_day,lf_cp_month)),list_raster_obj_files[tile_nb])
  filenames_NEX <- paste(lf_cp,collapse=" ")
  #filenames_NEX <- paste(list_tif_files_dates[[1]][[6]],list_tif_files_dates[[244]][[6]],lf_covar_tif[6]) #to get first date and tile 6 prediction mod1
  cmd_str <- paste("scp -p",filenames_NEX,paste(Atlas_hostname,Atlas_dir,sep=":"), sep=" ")
  system(cmd_str)
}


##################### END OF SCRIPT ######################
