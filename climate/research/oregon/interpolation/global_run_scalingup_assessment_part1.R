####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Accuracy methods are added in the the function scripts to evaluate results.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 03/23/2014  
#MODIFIED ON: 05/15/2014            
#Version: 3
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

function_analyses_paper1 <-"contribution_of_covariates_paper_interpolation_functions_05212014.R" #first interp paper
function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_05052014.R"

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
    tmp<-obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]]<- as.data.frame(tmp) #if spdf
  }
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

extract_daily_training_testing_info<- function(i,list_param){
  #This function extracts training and testing information from the raster object produced for each tile
  #This is looping through tiles...
  
  ### Function:
  pred_data_info_fun <- function(k,list_data,pred_mod,sampling_dat_info){
    
    data <- list_data[[k]]
    sampling_dat <- sampling_dat_info[[k]]
    if(data_day!="try-error"){
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

    #pred_data_day_s_info$method_interp <- rep(method_interp,nrow(pred_data_day_s_info)) 
    #pred_data_day_s_info$var_interp <- rep(var_interp,nrow(pred_data_day_s_info)) 
    #pred_data_day_s_info$tile_id <- rep(tile_id,nrow(pred_data_day_s_info)) 
    #pred_data_day_v_info <- do.call(rbind,list_pred_data_day_v_info)
    #pred_data_day_v_info$method_interp <- rep(method_interp,nrow(pred_data_day_v_info)) 
    #pred_data_day_v_info$var_interp <- rep(var_interp,nrow(pred_data_day_v_info)) 
    #pred_data_day_v_info$tile_id <- rep(tile_id,nrow(pred_data_day_v_info)) 
                                      
  }
  if(use_month==TRUE){
    
    list_data_month_s <- try(extract_list_from_list_obj(raster_obj$validation_mod_month_obj,"data_s"))
    list_data_month_v <- try(extract_list_from_list_obj(raster_obj$validation_mod_month_obj,"data_v"))
    sampling_dat_month <- extract_list_from_list_obj(raster_obj$clim_method_mod_obj,"sampling_month_dat")
    list_pred_data_month_s_info <- lapply(1:length(sampling_dat_month),FUN=pred_data_info_fun,
           list_data=list_data_month_s,pred_mod=pred_mod,sampling_dat_info=sampling_dat_month)
    list_pred_data_month_v_info <- lapply(1:length(sampling_dat_month),FUN=pred_data_info_fun,
           list_data=list_data_month_v,pred_mod=pred_mod,sampling_dat_info=sampling_dat_month)
    
    #combine training and testing later? also combined with accuracy
    pred_data_month_s_info <- do.call(rbind,list_pred_data_month_s_info)
    pred_data_month_v_info <- do.call(rbind,list_pred_data_month_v_info)
    
    pred_data_month_v_info$training <- rep(0,nrow(pred_data_month_v_info))
    pred_data_month_s_info$training <- rep(1,nrow(pred_data_month_v_info))
    pred_data_month_info <- rbind(pred_data_month_v_info,pred_data_month_s_info)
    
    pred_data_month_info$method_interp <- rep(method_interp,nrow(pred_data_month_info)) 
    pred_data_month_info$var_interp <- rep(var_interp,nrow(pred_data_month_info)) 
    pred_data_month_info$tile_id <- rep(tile_id,nrow(pred_data_month_info)) 

    #pred_data_month_s_info$method_interp <- rep(method_interp,nrow(pred_data_month_s_info)) 
    #pred_data_month_s_info$var_interp <- rep(var_interp,nrow(pred_data_month_s_info)) 
    #pred_data_month_s_info$tile_id <- rep(tile_id,nrow(pred_data_month_s_info)) 
    #pred_data_month_v_info$method_interp <- rep(method_interp,nrow(pred_data_month_v_info)) 
    #pred_data_month_v_info$var_interp <- rep(var_interp,nrow(pred_data_month_v_info)) 
    #pred_data_month_v_info$tile_id <- rep(tile_id,nrow(pred_data_month_v_info)) 

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
in_dir1 <- "/nobackupp4/aguzman4/climateLayers/output4" #On NEX

#in_dir_list <- list.dirs(path=in_dir1) #get the list of directories with resutls by 10x10 degree tiles
#use subset for now:

in_dir_list <- file.path(in_dir1,read.table(file.path(in_dir1,"processed.txt"))$V1)
#in_dir_list <- as.list(in_dir_list[-1])
#in_dir_list <- in_dir_list[grep("bak",basename(basename(in_dir_list)),invert=TRUE)] #the first one is the in_dir1
#in_dir_shp <- in_dir_list[grep("shapefiles",basename(in_dir_list),invert=FALSE)] #select directory with shapefiles...
in_dir_shp <- "/nobackupp4/aguzman4/climateLayers/output4/subset/shapefiles/"
#in_dir_list <- in_dir_list[grep("shapefiles",basename(in_dir_list),invert=TRUE)] 
#the first one is the in_dir1
# the last directory contains shapefiles 
y_var_name <- "dailyTmax"
interpolation_method <- c("gam_CAI")
out_prefix<-"run2_global_analyses_05122014"

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
                                   
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

##raster_prediction object : contains testing and training stations with RMSE and model object

list_raster_obj_files <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="^raster_prediction_obj.*.RData",full.names=T)})
basename(dirname(list_raster_obj_files[[1]]))
list_names_tile_coord <- lapply(list_raster_obj_files,FUN=function(x){basename(dirname(x))})
list_names_tile_id <- paste("tile",1:length(list_raster_obj_files),sep="_")
names(list_raster_obj_files)<- list_names_tile_id

lf_covar_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar_obj.*.RData",full.names=T)})
lf_covar_tif <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar.*.tif",full.names=T)})

########################## START SCRIPT ##############################

################################################################
######## PART 1: Generate tables to collect information 
######## over all tiles in North America 

##Function to collect all the tables from tiles into a table
###Table 1: Average accuracy metrics
###Table 2: daily accuracy metrics for all tiles

#First create table of tiles under analysis and their coord
df_tile_processed <- data.frame(tile_coord=basename(in_dir_list))
df_tile_processed$tile_id <- unlist(list_names_tile_id)
df_tile_processed$path_NEX <- in_dir_list
  
##Quick exploration of raster object
robj1 <- load_obj(list_raster_obj_files[[1]])
names(robj1)
names(robj1$method_mod_obj[[1]]) #for January 1, 2010
names(robj1$method_mod_obj[[1]]$dailyTmax) #for January

names(robj1$clim_method_mod_obj[[1]]$data_month) #for January
names(robj1$validation_mod_month_obj[[1]]$data_s) #for January with predictions

################
#### Table 1: Average accuracy metrics

#can use a maximum of 6 cores on the NEX Bridge
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


######################################################
####### PART 2 CREATE MOSAIC OF PREDICTIONS PER DAY ###

dates_l <- unique(robj1$tb_diagnostic_s$date) #list of dates to query tif

## make this a function? report on number of tiles used for mosaic...

#inputs: build a pattern to find files
y_var_name <- "dailyTmax"
interpolation_method <- c("gam_CAI")
name_method <- paste(interpolation_method,"_",y_var_name,"_",sep="")
l_pattern_models <- lapply(c(".*predicted_mod1_0_1.*",".*predicted_mod2_0_1.*",".*predicted_mod3_0_1.*",".*predicted_mod_kr_0_1.*"),
                           FUN=function(x){paste(x,dates_l,".*.tif",sep="")})
#gam_CAI_dailyTmax_predicted_mod_kr_0_1_20101231_30_145.0_-120.0.tif
out_prefix_s <- paste(name_method,c("predicted_mod1_0_01","predicted_mod2_0_01","predicted_mod3_0_01","predicted_mod_kr_0_1"),sep="")
dates_l #list of predicted dates
#l_out_rastnames_var <- paste(name_method,"predicted_mod1_0_01_",dates_l,sep="")
l_out_rastnames_var <- lapply(out_prefix_s,
                              FUN=function(x){paste(x,"_",dates_l,sep="")})
#gam_CAI_dailyTmax_predicted_mod_kr_0_1_20101231_30_145.0_-120.0.tif                    

##Get list of predicted tif across all tiles, models and dates...
#this takes time, use mclapply!!
lf_pred_tif <- vector("list",length=length(l_pattern_models)) #number of models is 3
for (i in 1:length(l_pattern_models)){
  l_pattern_mod <- l_pattern_models[[i]]
  list_tif_files_dates <-lapply(1:length(l_pattern_mod),FUN=list_tif_fun, 
                              in_dir_list=in_dir_list,pattern_str=l_pattern_models[[i]])
  lf_pred_tif[[i]] <- list_tif_files_dates
}

#Now get the clim surfaces:
month_l <- paste("clim_month_",1:12,sep="")
l_pattern_models <- lapply(c("_mod1_0_1.*","_mod2_0_1.*","_mod3_0_1.*","_mod_kr_0_1.*"),
                           FUN=function(x){paste("*.",month_l,x,".*.tif",sep="")})
#"CAI_TMAX_clim_month_11_mod2_0_145.0_-120.0.tif"
lf_clim_tif <- vector("list",length=nb_mod) #number of models is 3
for (i in 1:length(l_pattern_models)){
  l_pattern_mod <- l_pattern_models[[i]]
  list_tif_files_dates <- lapply(1:length(l_pattern_mod),FUN=list_tif_fun, 
                              in_dir_list=in_dir_list,pattern_str=l_pattern_models[[i]])
  lf_clim_tif[[i]] <- list_tif_files_dates
}



#### NOW create mosaic images
nb_mod <- 4

for (i in 1:nb_mod){
  
  #l_pattern_mod <- l_pattern_models[[i]]
  #out_prefix_s <-    
    
  #list_tif_files_dates <- list_tif_fun(1,in_dir_list,l_pattern_mod)

  ##List of predicted tif ...
  #list_tif_files_dates <-lapply(1:length(l_pattern_mod),FUN=list_tif_fun, 
  #                            in_dir_list=in_dir_list,pattern_str=l_pattern_mod)
  list_tif_files_dates <- lf_pred_tif[[i]] 
  #save(list_tif_files_dates, file=paste("list_tif_files_dates","_",out_prefix,".RData",sep=""))

  mosaic_list_var <- list_tif_files_dates
#  l_out_rastnames_var <- paste(name_method,"predicted_mod1_0_01_",dates_l,sep="")
#  out_rastnames_var <- l_out_rastnames_var[i]  
  #l_out_rastnames_var <- paste(name_method,"predicted_mod2_0_01_",dates_l,sep="")
  l_out_rastnames_var <- paste(name_method,"predicted_mod3_0_01_",dates_l,sep="")
  out_rastnames_var <- l_out_rastnames_var  

  file_format <- ".tif"
  NA_flag_val <- -9999

  j<-1
  list_param_mosaic<-list(j,mosaic_list_var,out_rastnames_var,out_dir,file_format,NA_flag_val)
  names(list_param_mosaic)<-c("j","mosaic_list","out_rastnames","out_path","file_format","NA_flag_val")
  list_var_mosaiced <- mclapply(1:2,FUN=mosaic_m_raster_list,list_param=list_param_mosaic,mc.preschedule=FALSE,mc.cores = 2)
  #list_var_mosaiced <- mclapply(1:365,FUN=mosaic_m_raster_list,list_param=list_param_mosaic,mc.preschedule=FALSE,mc.cores = 2)

  
}

### Now find out how many files were predicted

l_pattern_mod1<-paste(".*predicted_mod1_0_1.*",dates_l,".*.tif",sep="")

l_f_t12 <- list.files(path=in_dir_list[12],".*predicted_mod1_0_1.*")


l_f_bytiles<-lapply(in_dir_list,function(x){list.files(path=x,pattern=".*predicted_mod1_0_1.*")})
l_f_bytiles<-lapply(in_dir_list,function(x){list.files(path=x,pattern=".*predicted_mod1_0_1.*")})
l_f_bytiles<-lapply(in_dir_list,function(x){list.files(path=x,pattern=".*predicted_mod1_0_1.*")})


#system("scp -p ./*.tif parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/output_run2_global_analyses_05122014")

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
# use_day=TRUE
# use_month=TRUE
# 
# list_param_training_testing_info <- list(list_raster_obj_files,use_month,use_day,list_names_tile_id)
# names(list_param_training_testing_info) <- c("list_raster_obj_files","use_month","use_day","list_names_tile_id")
# 
# list_param <- list_param_training_testing_info
# 
# pred_data_info <- lapply(1:length(list_raster_obj_files),FUN=extract_daily_training_testing_info,list_param=list_param_training_testing_info)
# 
# pred_data_month_info <- do.call(rbind,lapply(pred_data_info,function(x){x$pred_data_month_info}))
# pred_data_day_info <- do.call(rbind,lapply(pred_data_info,function(x){x$pred_data_day_info}))

######################################################
####### PART 4: Get shapefile tiling with centroids ###

system("scp -p /nobackupp4/aguzman4/climateLayers/output4/subset/shapefiles/* parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/shapefiles")
#in_dir_shp <- "/nobackupp4/aguzman4/climateLayers/output4/subset/shapefiles/"

#get shape files for the region being assessed:
list_shp_global_tiles_files <- list.files(path=in_dir_shp,pattern="*.shp")
pattern_str <- basename(in_dir_list)
list_shp_reg_files <- lapply(pattern_str,function(x){list_shp_global_tiles_files[grep(x,invert=FALSE,list_shp_global_tiles_files)]}) #select directory with shapefiles...
df_tile_processed$shp_files <- unlist(list_shp_reg_files)

tx<-strsplit(as.character(df_tile_processed$tile_coord),"_")
lat<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][1]},x=tx))
long<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][2]},x=tx))
df_tile_processed$lat <- lat
df_tile_processed$lon <- long

#put that list in the df_processed and also the centroids!!
write.table(df_tile_processed,
            file=file.path(out_dir,paste("df_tile_processed_",out_prefix,".txt",sep="")),sep=",")

########### LAST PART: COPY SOME DATA BACK TO ATLAS #####

#Copy summary and mosaic back to atlas
system("scp -p ./*.txt parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/output_run2_global_analyses_05122014")
system("scp -p ./*.txt ./*.tif parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/output_run2_global_analyses_05122014")

#Copy specific tiles info back...
#tile_6: Oregon (45.0_-120.0)
system("scp -p ./*.txt ./*.tif parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/output_run2_05122014/output/45.0_-120.0")
#tile_8: Californi-Arizona
system("scp -p ./*.txt ./*.tif parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/output_run2_05122014/output/45.0_-120.0")
df_tile_processed$path_NEX
list_raster_obj_files[6]
list_raster_obj_files[8]
Atlas_dir <- "/data/project/layers/commons/NEX_data/output_run2_05122014/output/45.0_-120.0"
Atlas_hostname <- "parmentier@atlas.nceas.ucsb.edu"

#Oregon tile
filenames_NEX <- list_raster_obj_files[6] #copy raster prediction object
cmd_str <- paste("scp -p",filenames_NEX,paste(Atlas_hostname,Atlas_dir,sep=":"), sep=" ")
system(cmd_str)
#Now copy back tif for specific dates and tile (date 1 and date 244)
date_selected <- c("20100101","20100901")
date_index <- c(1,244)
tile_nb <- 6
nb_mod <- 3+1
lf_cp <- vector("list",length=length(date_selected))
for(i in 1:length(date_selected)){
  index <- date_index[i]  
  #get all predicted tmax files for all models and specific date, tile
  lf_cp[[i]] <- unlist(lapply(1:nb_mod,FUN=function(x){lf_pred_tif[[x]][[index]][[tile_nb]]}))
}

lf_clim_tiff[[]]

filenames_NEX <- paste(list_tif_files_dates[[1]][[6]],list_tif_files_dates[[244]][[6]],lf_covar_tif[6]) #to get first date and tile 6 prediction mod1
cmd_str <- paste("scp -p",filenames_NEX,paste(Atlas_hostname,Atlas_dir,sep=":"), sep=" ")
system(cmd_str)

#California-Arizona tile
Atlas_dir <- "/data/project/layers/commons/NEX_data/output_run2_05122014/output/35.0_-115.0"
filenames_NEX <- paste(list_raster_obj_files[8],lf_covar_obj[8]) #copy raster prediction object
cmd_str <- paste("scp -p",filenames_NEX,paste(Atlas_hostname,Atlas_dir,sep=":"), sep=" ")
system(cmd_str)
#Now copy back tif for specific dates and tile (date 1 and date 244)
filenames_NEX <- paste(list_tif_files_dates[[1]][[8]],list_tif_files_dates[[244]][[8]],lf_covar_tif[8]) #to get first date and tile 6 prediction mod1
cmd_str <- paste("scp -p",filenames_NEX,paste(Atlas_hostname,Atlas_dir,sep=":"), sep=" ")
system(cmd_str)

##################### END OF SCRIPT ######################
