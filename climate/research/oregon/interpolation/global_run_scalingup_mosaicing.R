##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Different options to explore mosaicing are tested.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 04/14/2015  
#MODIFIED ON: 06/05/2015            
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

sine_structure_fun <-function(x,T,phase,a,b,use_cos=FALSE){
  
  #Create sine for a one dimensional series
  #Note that sine function uses radian unit.
  #a=amplitude
  #b=mean or amplitude 0 of the series
  #T= stands for period definition
  #phase=phase angle (in radian!!)
  #cos: use cosine instead of sine if TRUE
  
  if(use_cos==FALSE){
    y <- a*sin((x*pi/T)+ phase) + b
  }else{
    y <- a*cos((x*pi/T)+ phase) + b
  }
  return(y)
}

create_weights_fun <- function(i, list_param){
  #This function generates weights from a point location on a raster layer.
  #Note that the weights are normatlized on 0-1 scale using max and min values.
  #Inputs:
  #lf: list of raster files
  #df_points: reference points from which to compute distance
  #r_feature: reference features as raster image from which to compute distance from
  #methods: options available: use_sine_weights,use_edge,use_linear_weights
  #Outputs:
  #raster list of weights and product of wegihts and inuts
  #TODO: 
  # -use gdal proximity for large files and use_edge option
  # - add raster options
  # - improve efficiency
  # - change name options
  #
  ############
  
  lf <- list_param$lf
  df_points <- list_param$df_points
  r_feature <- list_param$r_feature #this should be change to a list
  padding <- TRUE #if padding true then make buffer around edges??
  method <- list_param$method #differnt methods available to create weights
  #NAflag,file_format,out_suffix etc...
  out_dir_str <- list_param$out_dir_str
    
  ####### START SCRIPT #####
  
  r_in <- raster(lf[i]) #input image

  set1f <- function(x){rep(NA, x)}
  r_init <- init(r_in, fun=set1f)

  if(!is.null(r_feature)){
    r_init <- r_feature
  }

  if(!is.null(df_points)){ #reference points as SPDF object
    cell_ID <- cellFromXY(r_init,xy=df_points[i,])
    r_init[cell_ID] <- df_points$ID[i]
  }

  if(method=="use_sine_weights"){
    #Generate spatial pattern 5:     
    n_col <- ncol(r_init)
    n_row <- nrow(r_init)

    #u <- xFromCol(r_init,col=1:n_col)
    #add padding option later...buffer from a specific distance and tailling of at 0.1
    u <- 1:n_col
    a<- 1 #amplitude in this case
    b<- 0
    T<- n_col
    phase <- 0
    use_cos <- FALSE
    ux <- sine_structure_fun(u,T,phase,a,b,use_cos)
    ux_rep <-rep(ux,time=n_row)  
    r1 <-setValues(r_init,ux_rep)  #note efficient in memory might need to revise this
    #plot(r)
  
    v <- 1:n_row
    a<- 1 #amplitude in this case
    b<- 0
    T<- n_row
    phase <- 0
    use_cos <- FALSE
    vx <- sine_structure_fun(v,T,phase,a,b,use_cos)
    vx_rep <- unlist((lapply(1:n_row,FUN=function(i){rep(vx[i],time=n_col)})))  
  
    r2 <-setValues(r_init,vx_rep)  
    #plot(r2)
 
    r <- (r1+r2)*0.5  #combine patterns to make a elliptic surface, -0.5 could be modified
    #plot(r)
  }

  #change here to distance from edges..
  if(method=="use_edge"){ #does not work with large images
      #change...use gdal
      n_col <- ncol(r_init)
      n_row <- nrow(r_init)
    
      #xfrom
      r_init[1,1:n_col] <- 1
      r_init[n_row,1:n_col] <- 1
      r_init[1:n_row,1] <- 1
      r_init[1:n_row,n_col] <- 1
      r_dist <- distance(r_init)
      min_val <- cellStats(r_dist,min) 
      max_val <- cellStats(r_dist,max)
      r <- abs(r_dist - max_val)/ (max_val - min_val)
  } #too slow with R use http://www.gdal.org/gdal_proximity.html
  
  #plot(r_init)
  
  if(method=="use_linear_weights"){
    #
    r_dist <- distance(r_init)
    min_val <- cellStats(r_dist,min) 
    max_val <- cellStats(r_dist,max)
    r <- abs(r_dist - max_val)/ (max_val - min_val)
  }
  
  extension_str <- extension(lf[i])
  raster_name_tmp <- gsub(extension_str,"",basename(lf[i]))
  raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_",method,"_weights.tif",sep=""))
  writeRaster(r, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  
  
  r_var_prod <- r_in*r
  raster_name_prod <- file.path(out_dir_str, paste(raster_name_tmp,"_",method,"_prod_weights.tif",sep=""))
  writeRaster(r_var_prod, NAflag=NA_flag_val,filename=raster_name_prod,overwrite=TRUE)  
    
  weights_obj <- list(raster_name,raster_name_prod)
  names(weights_obj) <- c("r_weights","r_weights_prod")
  return(weights_obj)
}

mosaic_m_raster_list<-function(j,list_param){
  #This functions returns a subset of tiles from the modis grid.
  #Arguments: modies grid tile,list of tiles
  #Output: spatial grid data frame of the subset of tiles
  #Note that rasters are assumed to be in the same projection system!!
  #modified for global mosaic...still not working right now...
  
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
  }else{
    m_list <- mosaic_list[[j]]
  }
  input.rasters <- lapply(m_list, raster) #create raster image for each element of the list
  #inMemory(input.rasters[[1]])
  #note that input.rasters are not stored in memory!!
  mosaiced_rast<-input.rasters[[1]]
  
  for (k in 2:length(input.rasters)){
    mosaiced_rast<-mosaic(mosaiced_rast,input.rasters[[k]], tolerance=1,fun=mean)
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

raster_match <- function(i,list_param){
  ### Read in parameters/arguments
  lf_files <- list_param$lf_files
  rast_ref <- list_param$rast_ref #name of reference file
  file_format <- list_param$file_format #".tif",".rst" or others
  out_suffix <- list_param$out_suffix
  out_dir_str <- list_param$out_dir_str
    
  ### START SCRIPT ##
  
  r_m <- raster(rast_ref) #ref image with resolution and extent to match

  set1f <- function(x){rep(NA, x)}
  
  inFilename <- lf_files[i]
  
  extension_str <- extension(inFilename)
  raster_name_tmp <- gsub(extension_str,"",basename(inFilename))
  #outFilename <- file.path(out_dir,paste(raster_name_tmp,"_","m_",out_suffix,file_format,sep="")) #for use in function later...
  
  raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_","m_",out_suffix,file_format,sep=""))#output file
  r_ref <- init(r_m, fun=set1f, filename=raster_name, overwrite=TRUE)
  #NAvalue(r_ref) <- -9999

  cmd_str <- paste("/usr/bin/gdalwarp",inFilename,raster_name,sep=" ") 
  #gdalwarp -t_srs '+proj=utm +zone=11 +datum=WGS84' raw_spot.tif utm11.tif

  system(cmd_str)
  
  ##return name of file created
  return(raster_name)
}


############################################
#### Parameters and constants  

#Data is on ATLAS

in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test" #reg1 is North America and reg5 is Africa
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg1"
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg2" #Europe
#in_dir <- "/data/project/layers/commons/NEX_data/mosaicing_data_test/reg5" #Africa

y_var_name <- "dailyTmax" #PARAM1
interpolation_method <- c("gam_CAI") #PARAM2
region_name <- "reg2" #PARAM 13 #reg1 is North America, Africa Region 5

out_suffix <- paste(region_name,"_","mosaic_run10_1500x4500_global_analyses_06052015",sep="") 
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
                                   
tile_size <- "1500x4500" #PARAM 11
mulitple_region <- TRUE #PARAM 12

plot_region <- FALSE

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

lf_mosaic <-list.files(path=file.path(in_dir),    
           pattern=paste(".*.",day_to_mosaic[2],".*.tif$",sep=""),full.names=T) #choosing date 2...20100901
#lf_mosaic <- lf_mosaic[1:20]
r1 <- raster(lf_mosaic[1]) 
r2 <- raster(lf_mosaic[2]) 

plot(r1)
plot(r2)

lf <- sub(".tif","",lf_mosaic)
tx<-strsplit(as.character(lf),"_")

lat<- as.character(lapply(1:length(tx),function(i,x){x[[i]][13]},x=tx))
long<- as.character(lapply(1:length(tx),function(i,x){x[[i]][14]},x=tx))
lat <- as.character(lapply(1:length(lat),function(i,x){substr(x[[i]],2,nchar(x[i]))},x=lat)) #first number not in the coordinates

#Produce data.frame with centroids of each tiles...

df_centroids <- data.frame(long=as.numeric(long),lat=as.numeric(lat))
df_centroids$ID <- as.numeric(1:nrow(df_centroids))
coordinates(df_centroids) <- cbind(df_centroids$long,df_centroids$lat)
proj4string(df_centroids) <- projection(r1)

###############
### PART 2: prepare weights using tile rasters ############

out_dir_str <- out_dir
lf_r_weights <- vector("list",length=length(lf_mosaic))

#methods availbable:use_sine_weights,use_edge,use_linear_weights

### First use linear weights
method <- "use_linear_weights"
df_points <- df_centroids
#df_points <- NULL
r_feature <- NULL

list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,out_dir_str) 
names(list_param_create_weights) <- c("lf","df_points","r_feature","method","out_dir_str") 
num_cores <- 11

#debug(create_weights_fun)
weights_obj <- create_weights_fun(1,list_param=list_param_create_weights)

#This is the function creating the weights by tile. Distance from the centroids needs to be change from distance to
#the edges...can use rows and columsn to set edges to 1 and 0 for the others.
linear_weights_obj_list <- mclapply(1:length(lf_mosaic),FUN=create_weights_fun,list_param=list_param_create_weights,mc.preschedule=FALSE,mc.cores = num_cores)                           

list_linear_r_weights <- lapply(1:length(linear_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights},x=linear_weights_obj_list)
list_linear_r_weights_prod <- lapply(1:length(linear_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights_prod},x=linear_weights_obj_list)

### Second use sine weights
method <- "use_sine_weights"
#df_points <- df_centroids
df_points <- NULL
r_feature <- NULL

list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,out_dir_str) 
names(list_param_create_weights) <- c("lf","df_points","r_feature","method","out_dir_str") 
num_cores <- 11

#debug(create_weights_fun)
weights_obj <- create_weights_fun(1,list_param=list_param_create_weights)

#This is the function creating the weights by tile. Distance from the centroids needs to be change from distance to
#the edges...can use rows and columsn to set edges to 1 and 0 for the others.
sine_weights_obj_list <- mclapply(1:length(lf_mosaic),FUN=create_weights_fun,list_param=list_param_create_weights,mc.preschedule=FALSE,mc.cores = num_cores)                           

list_sine_r_weights <- lapply(1:length(sine_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights},x=sine_weights_obj_list)
list_sine_r_weights_prod <- lapply(1:length(sine_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights_prod},x=sine_weights_obj_list)

##

#Then scale on 1 to zero? or 0 to 1000
#e.g. for a specific pixel
#weight_sum=0.2 +0.4 +0.4+0.2=1.2
#val_w_sum= (0.2*val1+0.4*val2+0.4*val3+0.2*val4)
#no_valid = number of valid observation, needs to be added in the function
#m_val= sum_val/weight_sum #mean value
#

###############
### PART 3: prepare weightsfor mosaicing by matching extent ############

## Rasters tiles vary slightly in resolution, they need to be matched for the mosaic. Resolve issue in the 
#mosaic funciton using gdal_merge to compute a reference image to mach.
#outrastnames <- "reg1_mosaic_weights.tif"
#list_param_mosaic <- list(list_r_weights,out_dir,outrastnames,file_format,NA_flag_val,out_suffix)
#r1_projected <- projectRaster(raster(list_r_weights[[1]]),r)

cmd_str <- paste("python","/usr/bin/gdal_merge.py","-o avg.tif",paste(lf_mosaic,collapse=" ")) 
system(cmd_str)

## Create raster image for original predicted images with matching resolution and extent to the mosaic (reference image)

rast_ref <- file.path(out_dir,"avg.tif") #this is a the ref
r_ref <- raster(rast_ref)
plot(r_ref)

### First match weights from linear option
lf_files <- unlist(list_linear_r_weights)

list_param_raster_match <- list(lf_files,rast_ref,file_format,out_suffix,out_dir)
names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","out_suffix","out_dir_str")

#debug(raster_match)
#r_test <- raster(raster_match(1,list_param_raster_match))

list_linear_weights_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

lf_files <- unlist(list_linear_r_weights_prod)
list_param_raster_match <- list(lf_files,rast_ref,file_format,out_suffix,out_dir)
names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","out_suffix","out_dir_str")

num_cores <-11
list_linear_weights_prod_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

### Second match wegihts from sine option

lf_files <- unlist(list_sine_r_weights)

list_param_raster_match <- list(lf_files,rast_ref,file_format,out_suffix,out_dir)
names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","out_suffix","out_dir_str")

#debug(raster_match)
#r_test <- raster(raster_match(1,list_param_raster_match))

list_sine_weights_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

lf_files <- unlist(list_sine_r_weights_prod)
list_param_raster_match <- list(lf_files,rast_ref,file_format,out_suffix,out_dir)
names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","out_suffix","out_dir_str")

num_cores <-11
list_sine_weights_prod_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

#### Third use original images
#macth file to mosaic extent using the original predictions
lf_files <- lf_mosaic
list_param_raster_match <- list(lf_files,rast_ref,file_format,out_suffix,out_dir)
names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","out_suffix","out_dir_str")

list_pred_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

#####################
###### PART 4: compute the weighted mean with the mosaic function #####

#get the list of weights into raster objects
list_args_linear_weights <- list_linear_weights_m
#list_args_weights <- (mixedsort(list.files(pattern="r_weights_m_.*.tif")))
list_args_linear_weights <- lapply(1:length(list_args_linear_weights), FUN=function(i,x){raster(x[[i]])},x=list_args_linear_weights)

#get the list of weights product into raster objects
list_args_linear_weights_prod <- list_linear_weights_prod_m
#list_args_weights_prod <- (mixedsort(list.files(pattern="r_weights_prod_m_.*.tif")))
list_args_linear_weights_prod <- lapply(1:length(list_args_linear_weights_prod), FUN=function(i,x){raster(x[[i]])},x=list_args_linear_weights_prod)

#get the list of sine weights into raster objects
list_args_sine_weights <- list_sine_weights_m
#list_args_weights <- (mixedsort(list.files(pattern="r_weights_m_.*.tif")))
list_args_sine_weights <- lapply(1:length(list_args_sine_weights), FUN=function(i,x){raster(x[[i]])},x=list_args_sine_weights)

#get the list of weights product into raster objects
list_args_sine_weights_prod <- list_sine_weights_prod_m
#list_args_weights_prod <- (mixedsort(list.files(pattern="r_weights_prod_m_.*.tif")))
list_args_sine_weights_prod <- lapply(1:length(list_args_sine_weights_prod), FUN=function(i,x){raster(x[[i]])},x=list_args_sine_weights_prod)

#get the original predicted image to raster (matched previously to the mosaic extent)
list_args_pred_m <- list_pred_m
#list_args_pred_m <- (mixedsort(list.files(pattern="^gam_CAI.*.m_mosaic_run10_1500x4500_global_analyses_03252015.tif")))
list_args_pred_m <- lapply(1:length(list_args_pred_m), FUN=function(i,x){raster(x[[i]])},x=list_args_pred_m)

list_args_linear_weights$fun <- "sum" #we want the sum to compute the weighted mean
list_args_linear_weights$na.rm <- TRUE

list_args_linear_weights_prod$fun <- "sum"
list_args_linear_weights_prod$na.rm <- TRUE

list_args_sine_weights$fun <- "sum" #we want the sum to compute the weighted mean
list_args_sine_weights$na.rm <- TRUE

list_args_sine_weights_prod$fun <- "sum"
list_args_sine_weights_prod$na.rm <- TRUE

list_args_pred_m$fun <- "mean"
list_args_pred_m$na.rm <- TRUE

#Mosaic files
r_linear_weights_sum <- do.call(mosaic,list_args_linear_weights) #weights sum image mosaiced
r_linear_prod_sum <- do.call(mosaic,list_args_linear_weights_prod) #weights sum product image mosacied

r_m_linear_weighted_mean <- r_linear_prod_sum/r_linear_weights_sum #this is the mosaic using weighted mean...

r_sine_weights_sum <- do.call(mosaic,list_args_sine_weights) #weights sum image mosaiced
r_sine_prod_sum <- do.call(mosaic,list_args_sine_weights_prod) #weights sum product image mosacied

r_m_sine_weighted_mean <- r_sine_prod_sum/r_sine_weights_sum #this is the mosaic using weighted mean...

raster_name <- file.path(out_dir,paste("r_m_linear_weighted_mean_",out_suffix,".tif",sep=""))
writeRaster(r_m_linear_weighted_mean, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  

raster_name <- file.path(out_dir,paste("r_m_sine_weighted_mean_",out_suffix,".tif",sep=""))
writeRaster(r_m_sine_weighted_mean, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  

r_m_mean <- do.call(mosaic,list_args_pred_m) #this is unweighted mean from the predicted raster

raster_name <- file.path(out_dir,paste("r_m_mean_",out_suffix,".tif",sep=""))
writeRaster(r_m_mean, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  #unweighted mean

r_diff_weighted_mean <- r_m_linear_weighted_mean - r_m_sine_weighted_mean
#r_diff_weighted_mean<-r_diff_weighted_meam

r_diff_mean_linear <- r_m_mean - r_m_linear_weighted_mean 
r_diff_mean_sine <- r_m_mean - r_m_sine_weighted_mean 

r_m_mean_terrain <- terrain(r_m_mean,opt=c("slope","aspect"),unit="degrees")
r_m_sine_weighted_mean_terrain <- terrain(r_m_sine_weighted_mean,opt=c("slope","aspect"),unit="degrees")
r_m_linear_weighted_mean_terrain <- terrain(r_m_linear_weighted_mean,opt=c("slope","aspect"),unit="degrees")

#####################
###### PART 5: Now plot of the weighted mean and unweighted mean with the mosaic function #####

res_pix <- 1200
col_mfrow <- 1 
row_mfrow <- 0.8

png(filename=paste("Figure2_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_m_mean)

dev.off()

res_pix <- 1200
col_mfrow <- 1 
row_mfrow <- 0.8

png(filename=paste("Figure2_linear_weigthed_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_m_linear_weighted_mean )

dev.off()

res_pix <- 1200
col_mfrow <- 1 
row_mfrow <- 0.8

png(filename=paste("Figure2_sine_weigthed_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_m_sine_weighted_mean)

dev.off()

res_pix <- 1200
col_mfrow <- 1 
row_mfrow <- 0.8

png(filename=paste("Figure2_diff_linear_sine_weigthed_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_diff_weighted_mean)

dev.off()

res_pix <- 1200
col_mfrow <- 1 
row_mfrow <- 0.8

png(filename=paste("Figure2_diff_mean_linear_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_diff_mean_linear)

dev.off()

res_pix <- 1200
col_mfrow <- 1 
row_mfrow <- 0.8

png(filename=paste("Figure2_diff_mean_sine_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_diff_mean_sine)

dev.off()


#### plot terrain to emphasize possible edges..
res_pix <- 1200
col_mfrow <- 1 
row_mfrow <- 0.8

png(filename=paste("Figure2_slope_mean_linear_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_m_linear_weighted_mean_terrain,y=1)

dev.off()

png(filename=paste("Figure2_aspect_mean_linear_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_m_linear_weighted_mean_terrain,y=2)

dev.off()

png(filename=paste("Figure2_slope_mean_sine_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_m_sine_weighted_mean_terrain,y=1)

dev.off()

png(filename=paste("Figure2_aspect_mean_sine_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_m_sine_weighted_mean_terrain,y=2)

dev.off()

png(filename=paste("Figure2_slope_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_m_mean_terrain,y=1)

dev.off()

png(filename=paste("Figure2_aspect_mean_for_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_m_mean_terrain,y=2)

dev.off()

##################### END OF SCRIPT ######################

#################################################
#Ok testing on fake data to experiment and check methods:

##Quick function to generate test dataset
# make_raster_from_lf <- function(i,list_lf,r_ref){
#   vect_val <- list_lf[[i]]
#   r <-  r_ref
#   values(r) <-vect_val
#   #writeRaster...
#   return(r)
# }
# 
# vect_pred1 <- c(9,4,1,3,5,9,9,9,2)
# vect_pred2 <- c(10,3,1,2,4,8,7,8,2)
# vect_pred3 <- c(10,NA,NA,3,5,9,8,9,2)
# vect_pred4 <- c(9,3,2,NA,5,8,9,9,2)
# lf_vect_pred <- list(vect_pred1,vect_pred2,vect_pred3,vect_pred4)
# 
# vect_w1 <- c(0.2,0.5,0.1,0.3,0.4,0.5,0.5,0.3,0.2)
# vect_w2 <- c(0.3,0.4,0.1,0.3,0.4,0.5,0.7,0.1,0.2)
# vect_w3 <- c(0.5,0.3,0.2,0.2,0.3,0.6,0.7,0.3,0.2)
# vect_w4 <- c(0.2,0.5,0.3,0.3,0.4,0.5,0.5,0.2,0.2)
# lf_vect_w <- list(vect_w1,vect_w2,vect_w3,vect_w4)
# df_vect_w <-do.call(cbind,lf_vect_w)
# df_vect_pred <-do.call(cbind,lf_vect_pred)
# 
# tr_ref <- raster(nrows=3,ncols=3)
# 
# r_pred_l <- lapply(1:length(lf_vect_pred),FUN=make_raster_from_lf,list_lf=lf_vect_pred,r_ref=r_ref)
# r_w_l <- lapply(1:length(lf_vect_w),FUN=make_raster_from_lf,list_lf=lf_vect_w,r_ref=r_ref)
# 
# #r_w1<- make_raster_from_lf(2,list_lf=lf_vect_w,r_ref)
# 
# list_args_pred <- r_pred_l
# list_args_pred$fun <- "sum"
# 
# list_args_w <- r_w_l
# list_args_w$fun <- prod
# 
# r_test_val <-do.call(overlay,list_args) #sum
# r_test_w <-do.call(overlay,list_args_w) #prod
# 
# #need to do sumprod
# r1<- r_w_l[[1]]*r_pred_l[[1]]
# r2<- r_w_l[[2]]*r_pred_l[[2]]
# r3<- r_w_l[[3]]*r_pred_l[[3]]
# r4<- r_w_l[[4]]*r_pred_l[[4]]
# 
# r_pred <- stack(r_pred_l)
# r_w <- stack(r_w_l)
# 
# list_args_pred <- r_pred_l
# list_args_pred$fun <- mean
# list_args_pred$na.rm <- TRUE
# #r_sum_pred <-do.call(overlay,list_args_pred) #prod
# 
# #r_sum_pred <-do.call(mean,list_args_pred) #prod
# r_sum_pred <-do.call(mosaic,list_args_pred) #prod
# 
# list_args_pred$na.rm <- FALSE
# r_sum_pred <-do.call(overlay,list_args_pred) #prod
# 
# r_sum_pred <-do.call(overlay,list_args_w) #prod
# 
# list_args_w$fun <- sum
# r_sum_w <-do.call(overlay,list_args_w) #prod
# 
# r_m_w <- ((r1+r2+r3+r4)/(r_sum_w)) #mean weiated
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
