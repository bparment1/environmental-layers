##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Different options to explore mosaicing are tested.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 04/14/2015  
#MODIFIED ON: 05/07/2015            
#Version: 4
#PROJECT: Environmental Layers project     
#COMMENTS: analyses for run 10 global analyses,all regions 1500x4500km and other tiles
#TODO:
#1) Split functions and master script
#2) Make this is a script/function callable from the shell/bash
#3) Check image format for tif

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


plot_daily_mosaics <- function(i,list_param_plot_daily_mosaics){
  #Purpose:
  #This functions mask mosaics files for a default range (-100,100 deg).
  #It produces a masked tif in a given dataType format (FLT4S)
  #It creates a figure of mosaiced region being interpolated.
  #Author: Benoit Parmentier
  #Parameters:
  #lf_m: list of files 
  #reg_name:region name with tile size included
  #To do:
  #Add filenames
  #Add range
  #Add output dir
  #Add dataType_val option
  
  ##### BEGIN ########
  
  #Parse the list of parameters
  lf_m <- list_param_plot_daily_mosaics$lf_m
  reg_name <- list_param_plot_daily_mosaics$reg_name
  out_dir_str <- list_param_plot_daily_mosaics$out_dir_str
  out_suffix <- list_param_plot_daily_mosaics$out_suffix
  l_dates <- list_param_plot_daily_mosaics$l_dates


  #list_param_plot_daily_mosaics <- list(lf_m=lf_m,reg_name=reg_name,out_dir_str=out_dir_str)

  
  #rast_list <- vector("list",length=length(lf_m))
  r_test<- raster(lf_m[i])

  m <- c(-Inf, -100, NA,  
         -100, 100, 1, 
         100, Inf, NA) #can change the thresholds later
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  rc <- reclassify(r_test, rclmat,filename=paste("rc_tmp_",i,".tif",sep=""),dataType="FLT4S",overwrite=T)
  file_name <- unlist(strsplit(basename(lf_m[i]),"_"))
  
  #date_proc <- file_name[7] #specific tot he current naming of files
  date_proc <- l_dates[i]
  #paste(raster_name[1:7],collapse="_")
  #add filename option later
  extension_str <- extension(filename(r_test))
  raster_name_tmp <- gsub(extension_str,"",basename(filename(r_test)))
  raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_masked.tif",sep=""))
  r_pred <- mask(r_test,rc,filename=raster_name,overwrite=TRUE)
  
  res_pix <- 1200
  #res_pix <- 480

  col_mfrow <- 1
  row_mfrow <- 1
  
  png(filename=paste("Figure9_clim_mosaics_day_test","_",date_proc,"_",reg_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_pred,main=paste("Predicted on ",date_proc , " ", reg_name,sep=""),cex.main=1.5)
  dev.off()
  
  return(raster_name)
  
}

plot_screen_raster_val<-function(i,list_param){
  ##USAGE###
  #Screen raster list and produced plot as png.
  fname <-list_param$lf_raster_fname[i]
  screenRast <- list_param$screenRast
  l_dates<- list_param$l_dates
  out_dir_str <- list_param$out_dir_str
  prefix_str <-list_param$prefix_str
  out_suffix_str <- list_param$out_suffix_str
  
  ### START SCRIPT ####
  date_proc <- l_dates[i]
  
  if(screenRast==TRUE){
    r_test <- raster(fname)

    m <- c(-Inf, -100, NA,  
         -100, 100, 1, 
         100, Inf, NA) #can change the thresholds later
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    rc <- reclassify(r_test, rclmat,filename=paste("rc_tmp_",i,".tif",sep=""),dataType="FLT4S",overwrite=T)
    #file_name <- unlist(strsplit(basename(lf_m[i]),"_"))
    extension_str <- extension(filename(r_test))
    raster_name_tmp <- gsub(extension_str,"",basename(filename(r_test)))
    raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_masked.tif",sep=""))
  
    r_pred <- mask(r_test,rc,filename=raster_name,overwrite=TRUE)
  }else{
    r_pred <- raster(fname)
  }
  
  #date_proc <- file_name[7] #specific tot he current naming of files
  #date_proc<- "2010010101"

  #paste(raster_name[1:7],collapse="_")
  #add filename option later

  res_pix <- 960
  #res_pix <- 480

  col_mfrow <- 1
  row_mfrow <- 1
  
#  png(filename=paste("Figure10_clim_world_mosaics_day_","_",date_proc,"_",tile_size,"_",out_suffix,".png",sep=""),
#    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  png(filename=paste(prefix_str,"_",date_proc,"_",tile_size,"_",out_suffix_str,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_pred,main=paste("Predicted on ",date_proc , " ", tile_size,sep=""),cex.main=1.5)
  dev.off()

}

create_weights_fun <- function(i, list_param){
  #This function generates weights from a point location on a raster layer.
  #Note that the weights are normatlized on 0-1 scale using max and min values.
  #Inputs:
  #lf: list of raster files
  #df_points: 
  #Outputs:
  #
  ############
  
  lf <- list_param$lf
  df_centroids <- list_param$df_points
  out_dir_str <- list_param$out_dir_str
    
  ####### START SCRIPT #####
  
  r1 <- raster(lf[i]) #input image

  set1f <- function(x){rep(NA, x)}
  r_init <- init(r1, fun=set1f)

  cell_ID <- cellFromXY(r_init,xy=df_centroids[i,])
  r_init[cell_ID] <- df_centroids$ID[i]

  r_dist <- distance(r_init)
  min_val <- cellStats(r_dist,min) 
  max_val <- cellStats(r_dist,max)
  r_dist_normalized <- abs(r_dist - max_val)/ (max_val - min_val)
  
  extension_str <- extension(lf_mosaic_pred_1500x4500[i])
  raster_name_tmp <- gsub(extension_str,"",basename(lf_mosaic_pred_1500x4500[i]))
  raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_weights.tif",sep=""))
  writeRaster(r_dist_normalized, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  
  
  r_var_prod <- r1*r_dist_normalized
  raster_name_prod <- file.path(out_dir_str, paste(raster_name_tmp,"_prod_weights.tif"))
  writeRaster(r_var_prod, NAflag=NA_flag_val,filename=raster_name_prod,overwrite=TRUE)  
    
  weights_obj <- list(raster_name,raster_name_prod)
  names(weights_obj) <- c("r_weights","r_weights_prod")
  return(weights_obj)
}

############################################
#### Parameters and constants  

#on ATLAS

in_dir <- "~/Data/IPLANT_project/mosaicing_data_test/reg2"

y_var_name <- "dailyTmax" #PARAM1
interpolation_method <- c("gam_CAI") #PARAM2
out_suffix <- "mosaic_run10_1500x4500_global_analyses_03252015" #PARAM3
out_dir <- in_dir #PARAM4
create_out_dir_param <- TRUE #PARAM 5

mosaic_plot <- FALSE #PARAM6

#if daily mosaics NULL then mosaicas all days of the year
day_to_mosaic <- c("20100101","20100102","20100103","20100104","20100105",
                   "20100301","20100302","20100303","20100304","20100305",
                   "20100501","20100502","20100503","20100504","20100505",
                   "20100701","20100702","20100703","20100704","20100705",
                   "20100901","20100902","20100903","20100904","20100905",
                   "20101101","20101102","20101103","20101104","20101105") #PARAM7
  
#CRS_locs_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
CRS_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84 #CONSTANT1
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

proj_str<- CRS_WGS84 #PARAM 8 #check this parameter
file_format <- ".rst" #PARAM 9
NA_value <- -9999 #PARAM10
NA_flag_val <- NA_value
                                   
tile_size <- "1500x4500" #PARAM 11
mulitple_region <- TRUE #PARAM 12

region_name <- "world" #PARAM 13
plot_region <- FALSE

########################## START SCRIPT ##############################


####### PART 1: Read in data ########

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

setwd(out_dir)

###Table 1: Average accuracy metrics
################## PLOTTING WORLD MOSAICS ################

lf_mosaic_pred_1500x4500 <-list.files(path=file.path(in_dir),    
           pattern=paste(".*.tif$",sep=""),full.names=T) 

r1 <- raster(lf_mosaic_pred_1500x4500[1]) 
r2 <- raster(lf_mosaic_pred_1500x4500[2]) 


lf <- sub(".tif","",lf_mosaic_pred_1500x4500)
tx<-strsplit(as.character(lf),"_")

lat<- as.character(lapply(1:length(tx),function(i,x){x[[i]][13]},x=tx))
long<- as.character(lapply(1:length(tx),function(i,x){x[[i]][14]},x=tx))
lat <- as.character(lapply(1:length(lat),function(i,x){substr(x[[i]],2,nchar(x[i]))},x=lat)) #first number not in the coordinates

df_centroids <- data.frame(long=as.numeric(long),lat=as.numeric(lat))
df_centroids$ID <- as.numeric(1:nrow(df_centroids))
#
extract(r1,)
coordinates(df_centroids) <- cbind(df_centroids$long,df_centroids$lat)
proj4string(df_centroids) <- projection(r1)
#centroid distance
#c1 <- gCentroid(x,byid=TRUE) 
#pt <- gCentroid(shp1)

#Make this a function later...
#then distance...
out_dir_str <- out_dir
lf_r_weights <- vector("list",length=length(lf_mosaic_pred_1500x4500))

list_param_create_weights <- list(lf_mosaic_pred_1500x4500,df_centroids,out_dir_str) 
names(list_param_create_weights) <- c("lf","df_points","out_dir_str") 

create_weights_fun
num_cores <- 6

#debug(create_weights_fun)
weights_obj <- create_weights_fun(1,list_param=list_param_create_weights)

weights_obj_list <- mclapply(1:length(lf_mosaic_pred_1500x4500),FUN=create_weights_fun,list_param=list_param_create_weights,mc.preschedule=FALSE,mc.cores = num_cores)                           

list_args_weights_prod <- lapply(1:length(weights_obj_list), FUN=function(x){raster(x[[i]]$r_weights_prod})}
list_args_weights_prod$fun <- "sum"


#"r_weights","r_weights_prod"

list_args_weights <- lapply(1:length(weights_obj_list), FUN=function(i,x){raster(x[[i]]$r_weights)},x=weights_obj_list)
list_args_weights_prod <- lapply(1:length(weights_obj_list), FUN=function(i,x){raster(x[[i]]$r_weights_prod)},x=weights_obj_list)

list_args_weights$fun <- "sum"
#list_args_weights$fun <- "mean"

list_args_weights$na.rm <- TRUE
list_args_weights$tolerance <- 1

list_args_weights_prod$fun <- "sum"
list_args_weights_prod$na.rm <- TRUE
list_args_weights_prod$na.rm <- TRUE

r_weights_sum <- do.call(mosaic,list_args_weights) #sum
r_prod_sum <- do.call(mosaic,list_args_weights_prod) #sum

r_weighted_mean <- r_weights_sum/r_prod_sum


#r_weights_sum <- do.call(overlay,list_args_weights) #sum
#r_weights_sum <- do.call(overlay,list_args_weights) #sum

#r_test_w <-do.call(overlay,list_args_w) #prod


###Mosaic with do.call...
#rasters1.mosaicargs <- rasters1
#rasters1.mosaicargs$fun <- mean
#mos2 <- do.call(mosaic, rasters1.mosaicargs)

#Then scale on 1 to zero? or 0 to 1000
#e.g. for a specific pixel
#weight_sum=0.2 +0.4 +0.4+0.2=1.2
#val_w_sum= (0.2*val1+0.4*val2+0.4*val3+0.2*val4)
#no_valid = 
#m_val= sum_val/weight_sum #mean value
#

#in raster term
#r_weights_sum <- ...
#r_val_w_sum <-  
#

#################################################
#Ok testing on fake data:

##Quick function to generate test dataset
make_raster_from_lf <- function(i,list_lf,r_ref){
  vect_val <- list_lf[[i]]
  r <-  r_ref
  values(r) <-vect_val
  #writeRaster...
  return(r)
}

vect_pred1 <- c(9,4,1,3,5,9,9,9,2)
vect_pred2 <- c(10,3,1,2,4,8,7,8,2)
vect_pred3 <- c(10,NA,NA,3,5,9,8,9,2)
vect_pred4 <- c(9,3,2,NA,5,8,9,9,2)
lf_vect_pred <- list(vect_pred1,vect_pred2,vect_pred3,vect_pred4)

vect_w1 <- c(0.2,0.5,0.1,0.3,0.4,0.5,0.5,0.3,0.2)
vect_w2 <- c(0.3,0.4,0.1,0.3,0.4,0.5,0.7,0.1,0.2)
vect_w3 <- c(0.5,0.3,0.2,0.2,0.3,0.6,0.7,0.3,0.2)
vect_w4 <- c(0.2,0.5,0.3,0.3,0.4,0.5,0.5,0.2,0.2)
lf_vect_w <- list(vect_w1,vect_w2,vect_w3,vect_w4)
df_vect_w <-do.call(cbind,lf_vect_w)
df_vect_pred <-do.call(cbind,lf_vect_pred)

tr_ref <- raster(nrows=3,ncols=3)

r_pred_l <- lapply(1:length(lf_vect_pred),FUN=make_raster_from_lf,list_lf=lf_vect_pred,r_ref=r_ref)
r_w_l <- lapply(1:length(lf_vect_w),FUN=make_raster_from_lf,list_lf=lf_vect_w,r_ref=r_ref)

#r_w1<- make_raster_from_lf(2,list_lf=lf_vect_w,r_ref)

list_args_pred <- r_pred_l
list_args_pred$fun <- "sum"

list_args_w <- r_w_l
list_args_w$fun <- prod

r_test_val <-do.call(overlay,list_args) #sum
r_test_w <-do.call(overlay,list_args_w) #prod

#need to do sumprod
r1<- r_w_l[[1]]*r_pred_l[[1]]
r2<- r_w_l[[2]]*r_pred_l[[2]]
r3<- r_w_l[[3]]*r_pred_l[[3]]
r4<- r_w_l[[4]]*r_pred_l[[4]]

r_pred <- stack(r_pred_l)
r_w <- stack(r_w_l)

list_args_pred <- r_pred_l
list_args_pred$fun <- mean
list_args_pred$na.rm <- TRUE
#r_sum_pred <-do.call(overlay,list_args_pred) #prod

#r_sum_pred <-do.call(mean,list_args_pred) #prod
r_sum_pred <-do.call(mosaic,list_args_pred) #prod

list_args_pred$na.rm <- FALSE
r_sum_pred <-do.call(overlay,list_args_pred) #prod

r_sum_pred <-do.call(overlay,list_args_w) #prod

list_args_w$fun <- sum
r_sum_w <-do.call(overlay,list_args_w) #prod

r_m_w <- ((r1+r2+r3+r4)/(r_sum_w)) #mean weiated
#n33e to check the result!! especially the nubmer of valid pix val

#r_test_val <-do.call(overlay,list_args) #sum

#can do mosaic with sum?? for both weighted sum and val
#
#can use gdal calc...

#r_m <- r1 + r2
#name_method <- paste(interpolation_method,"_",y_var_name,"_",sep="")
##Use python code written by Alberto Guzman

#system("MODULEPATH=$MODULEPATH:/nex/modules/files")
#system("module load /nex/modules/files/pythonkits/gdal_1.10.0_python_2.7.3_nex")
#lf1 <- lf_world_pred_1000x3000
#lf2 <- lf_world_pred_1500x4500

#module_path <- ""
#module_path <- "/nobackupp6/aguzman4/climateLayers/sharedCode/"
#sh /nobackupp6/aguzman4/climateLayers/sharedCode/gdalCalDiff.sh file1.tif file2.tif output.tif
#/nobackupp6/aguzman4/climateLayers/sharedCode/mosaicUsingGdalMerge.py
#l_dates <- paste(day_to_mosaic,collapse=",",sep=" ")
#l_dates <- paste(day_to_mosaic,collapse=",")
## use region 2 first
#lf_out <- paste("diff_world_","1000_3000","by1500_4500_","mod1","_",l_dates,out_suffix,"_",file_format,sep="")


#for (i in 1:length(lf_out)){
#  out_file <- lf_out[i]
#  in_file1 <- lf1[i]
#  in_file2 <- lf2[i]
#    
#  cmd_str <- paste("sh", file.path(module_path,"gdalCalDiff.sh"),
#                 in_file1,
#                 in_file2,
#                 out_file,sep=" ")
#  system(cmd_str)
#
#}

##################### END OF SCRIPT ######################
