##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Different options to explore mosaicing are tested. This script only contains functions.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 04/14/2015  
#MODIFIED ON: 04/11/2016            
#Version: 2
#PROJECT: Environmental Layers project     
#COMMENTS: first commit of function script to test mosaicing using 1500x4500km and other tiles
#TODO:
#1) Make this is a script/function callable from the shell/bash
#2) Improve performance: there will be a need to improve efficiency for the workflow.

#Error message for gdal_proximity:
#ERROR 1: Source and proximity bands are not the same size.
#gdal_proximity.py give an error when tries to replace an existent output file with different band.
#If the output already exists and was created by a different band input file, this error is displayed: ERROR 1: Source and proximity bands are not the same size.
#If you remove the existent file, it works fine.
#available:
#See below

#################################################################################################

### Loading R library and packages        
#library used in the workflow production:
library(gtools)                              # loading some useful tools 
library(mgcv)                                # GAM package by Simon Wood
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gstat)                               # Kriging and co-kriging by Pebesma et al., not on NEX
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
library(spgwr)                               # GWR method, not on NEX
library(automap)                             # Kriging automatic fitting of variogram using gstat, not on NEX
library(rgeos)                               # Geometric, topologic library of functions
library(RPostgreSQL)                         # Interface R and Postgres, not used in this script
library(gridExtra)
#Additional libraries not used in workflow
library(pgirmess)                            # Krusall Wallis test with mulitple options, Kruskalmc {pgirmess}  
library(colorRamps)
library(zoo)
library(xts)

#### FUNCTION USED IN SCRIPT

##List all the functions in this script:

#[1] "autoKrige_fun"                      
#[2] "create_accuracy_metric_raster"     
#[3] "create_accuracy_residuals_raster"   
#[4] "create_weights_fun"                
# [5] "fit_models"                         "function_mosaicing"                
# [7] "in_dir_script"                      "mosaicFiles"                       
# [9] "mosaic_m_raster_list"               "mosaic_python_merge"               
#[11] "plot_daily_mosaics"                 "plot_diff_raster"                  
#[13] "plot_mosaic"                        "plot_screen_raster_val"            
#[15] "predict_accuracy_raster_by_station" "predict_auto_krige_raster_model"   
#[17] "raster_match"                       "remove_na_spdf"                    
#[19] "select_var_stack"                   "sine_structure_fun"    

############## START FUNCTIONS DEFINITIONS ####

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

## Add numcores: done
## use mclapply
create_accuracy_metric_raster <- function(i, list_param){
  #This function generates weights from a point location on a raster layer.
  #Note that the weights are normatlized on 0-1 scale using max and min values.
  #Inputs:
  #lf: list of raster files
  #tb: data.frame table #fitting or validation table with all days
  #metric_name: accuracy metric selected to be mapped, RMSE, MAE etc.
  #pred_mod_name: model selected, such as mod1, mod_kr etc.
  #y_var_name: variable being modeled e.g."dailyTmax", dailyTmin, precip  
  #interpolation_method: names of the interpolation/modeling method
  #date_processed: day being processed , e.g. 19920101
  #num_cores : number of cores used in the parallelization
  #NA_flag_val: value used as flag in the raster 
  #file_format: e.g. tif., .rst
  #out_dir_str: output directory
  #out_suffix_str: output suffix
  #Outputs:
  #raster list of weights and product of wegihts and inuts
  #TODO: 

  # - improve efficiency
  #
  ############
  
  #Functions
  create_raster_df_centroids_fun <- function(j,list_param){
    #This function generates raster images from metrics values from a data.frame.
    #The raster layer is assigned a unique value from the pixel at every location.
    #Input Parameters:
    #  #lf: list of raster files
    #df_centroids: data.frame table #fitting or validation table with all days
    #metric_name: accuracy metric selected to be mapped, RMSE, MAE etc.
    #date_processed: day being processed , e.g. 19920101
    #num_cores : number of cores used in the parallelization
    #NA_flag_val: value used as flag in the raster 
    #file_format: e.g. tif., .rst
    #out_dir_str: output directory
    #out_suffix_str: output suffix
    #Outputs:
    
    #### PARSE arguments
    
    df_centroids <- list_param$df_centroids
    metric_name <- list_param$metric_name 
    #interpolation_method <- list_param$interpolation_method #c("gam_CAI") #PARAM3
    #date_processed <- list_param$metric_name 
    #num_cores <- list_param$num_cores
    NA_flag_val <- list_param$NA_flag_val
    file_format <- list_param$file_format
    out_dir_str <- list_param$out_dir
    out_suffix_str <- list_param$out_suffix

    ####### START SCRIPT #####    
    
    inFilename <- df_centroids$files[j]
    r1 <- raster(inFilename)
    r1[] <- df_centroids[[metric_name]][j] #improve this
    #set1f <- function(x){rep(NA, x)}
    #r_init <- init(r_in, fun=set1f)
    lf_tmp <- gsub(file_format,"",lf)
  
    extension_str <- extension(inFilename)
    raster_name_tmp <- gsub(extension_str,"",basename(inFilename))
    outFilename <- file.path(out_dir_str,paste(raster_name_tmp,"_",metric_name,"_",out_suffix,file_format,sep="")) #for use in function later...
  
    writeRaster(r1, NAflag=NA_flag_val,filename=outFilename,overwrite=TRUE)  
    #list_raster_name[[j]] <- outFilename
    return(outFilename)
  }
  
  ####### PARSE ARGUMENTS
  
  lf <- list_param$lf[[i]] #list of files to mosaic
  tb <- list_param$tb #fitting or validation table with all days
  metric_name <- list_param$metric_name #RMSE, MAE etc.
  pred_mod_name <- list_param$pred_mod_name #mod1, mod_kr etc.
  y_var_name <- list_param$y_var_name #"dailyTmax" #PARAM2
  interpolation_method <- list_param$interpolation_method #c("gam_CAI") #PARAM3
  date_processed <- list_param$days_to_process[i]
  num_cores <- list_param$num_cores #number of cores used
  NA_flag_val <- list_param$NA_flag_val
  #NAflag,file_format,out_suffix etc...
  file_format <- list_param$file_format
  out_dir_str <- list_param$out_dir
  out_suffix_str <- list_param$out_suffix
   
  ####### START SCRIPT #####
  
  lf_tmp <- gsub(file_format,"",lf)
  tx<-strsplit(as.character(lf_tmp),"_")
  #deal with the fact that we have number "1" attached to the out_suffix (centroids of tiles)
  pos_lat <- lapply(1:length(tx),function(i,x){length(x[[i]])-1},x=tx)
  pos_lon <- lapply(1:length(tx),function(i,x){length(x[[i]])},x=tx)
  lat_val <- unlist(lapply(1:length(tx),function(i,x,y){x[[i]][pos_lat[[i]]]},x=tx,y=pos_lat))
  lat <- as.character(lapply(1:length(lat_val),function(i,x){substr(x[[i]],2,nchar(x[i]))},x=lat_val)) #first number not in the coordinates
  long <- as.character(lapply(1:length(tx),function(i,x,y){x[[i]][pos_lon[[i]]]},x=tx,y=lon_lat))

  df_centroids <- data.frame(long=as.numeric(long),lat=as.numeric(lat))
  df_centroids$ID <- as.numeric(1:nrow(df_centroids))
  df_centroids$tile_coord <- paste(lat,long,sep="_")
  df_centroids$files <- lf
  df_centroids$date <- date_processed
  write.table(df_centroids,paste("df_centroids_",date_processed,"_",out_suffix,".txt",sep=""),sep=',')

  #sprintf(" %3.1f", df_centroids$lat)

  tb_date <- subset(tb,date==date_processed & pred_mod==pred_mod_name)
  tb_date$tile_coord <- as.character(tb_date$tile_coord)
  df_centroids <- merge(df_centroids,tb_date,by="tile_coord")

  #use mclapply  
  #list_raster_name <- vector("list",length=length(lf))
  list_param_raster_df_centroids <- list(df_centroids,metric_name,NA_flag_val,file_format,out_dir_str,out_suffix)
  names(list_param_raster_df_centroids) <- c("df_centroids","metric_name","NA_flag_val","file_format","out_dir","out_suffix")

  #undebug(create_raster_df_centroids_fun)
  #test_lf <- lapply(1,FUN=create_raster_df_centroids_fun,list_param=list_param_raster_df_centroids)                           
  
  list_raster_name <- mclapply(1:length(lf),FUN=create_raster_df_centroids_fun,list_param=list_param_raster_df_centroids,mc.preschedule=FALSE,mc.cores = num_cores)                           

  raster_created_obj <- list(list_raster_name,df_centroids)
  names(raster_created_obj) <- c("list_raster_name","df_centroids")
  return(raster_created_obj)
}

#### end of function

mosaic_python_merge <- function(NA_flag_val,module_path,module_name,input_file,out_mosaic_name){
  #out_mosaic_name <- r_weights_sum_raster_name <- file.path(out_dir,paste("r_weights_sum_m_",method_str,"_weighted_mean_",out_suffix,".tif",sep=""))
  cmd_str <- paste("python", file.path(module_path,module_name),
                   "--config GDAL_CACHEMAX=1500",
                   "--overwrite=TRUE",
                   paste("-o",out_mosaic_name,sep=" "),
                   paste("--optfile", input_file,sep=" "),
                   paste("-n",NA_flag_val,sep=" "))
  system(cmd_str)
  #list(out_mosaic_name,cmd_str)
  mosaic_python_merge_obj <- list(out_mosaic_name,cmd_str)
  names(mosaic_python_merge_obj) <- c("out_mosaic_name","cmd_str")
  
  return(mosaic_python_merge_obj)
}
      
      
create_weights_fun <- function(i, list_param){
  #This function generates weights from a point location on a raster layer.
  #Note that the weights are normatlized on 0-1 scale using max and min values.
  #Inputs:
  #1)lf: list of raster files
  #2)df_points: reference points from which to compute distance
  #3)r_feature: reference features as raster image from which to compute distance from
  #4)methods: options available: use_sine_weights,use_edge,use_linear_weights
  #5)NA_flag : raster flag values, e.g. -9999
  #6)file_format: raster format used, default is ".tif"
  #7)out_suffix_str: output suffix, default is NULL, it is recommended to add the variable name etc.
  #             here e.g. dailyTmax and date!!
  #8)out_dir_str: output directory, default is NULL

  #Outputs:
  #raster list of weights and product of wegihts and inuts
  #TODO: 
  # -use gdal proximity for large files and use_edge option
  # - add raster options
  # - improve efficiency
  # - change name options
  #
  ############
  
  ##### START SCRIPT #####
  
  ##### Parse out the input parameters
  
  lf <- list_param$lf
  df_points <- list_param$df_points
  r_feature <- list_param$r_feature #this should be change to a list
  padding <- TRUE #if padding true then make buffer around edges??
  method <- list_param$method #differnt methods available to create weights
  #NAflag,file_format,out_suffix etc...
  NA_flag_val <- list_param$NA_flag
  file_format <- list_param$file_format
  out_suffix_str <- list_param$out_suffix_str
  out_dir_str <- list_param$out_dir_str
    
  ##### Prepare weight layers  
  
  r_in <- raster(lf[i]) #input image
  tile_no <- i #file being processed, assuming tiles by tiles
  
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
    vx_rep <- unlist((lapply(1:n_row,FUN=function(j){rep(vx[j],time=n_col)})))  
  
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
      #r_dist <- distance(r_init)
      #out_suffix_str
      srcfile <- file.path(out_dir_str,paste("feature_target_",tile_no,"_",Sys.getpid(),out_suffix_str,file_format,sep=""))

      writeRaster(r_init,filename=srcfile,overwrite=T)
      #Sys.getpid
      dstfile <- file.path(out_dir_str,paste("feature_target_edge_distance",tile_no,"_",Sys.getpid(),out_suffix_str,file_format,sep=""))
      n_values <- "1"
        
      cmd_str <- paste("gdal_proximity.py", srcfile, dstfile,"-values",n_values,sep=" ")
      system(cmd_str)
      r_dist<- raster(dstfile)
      min_val <- cellStats(r_dist,min) 
      max_val <- cellStats(r_dist,max)
      r <- abs(r_dist - min_val)/ (max_val - min_val) #no need to inverse...
  } #too slow with R so used http://www.gdal.org/gdal_proximity.html
  
  if(method=="use_linear_weights"){
    #
    r_dist <- distance(r_init)
    min_val <- cellStats(r_dist,min) 
    max_val <- cellStats(r_dist,max)
    r <- abs(r_dist - max_val)/ (max_val - min_val)
  }
  
  extension_str <- extension(lf[i])
  raster_name_tmp <- gsub(extension_str,"",basename(lf[i]))
  raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_",method,"_weights_",out_suffix_str,file_format,sep=""))
  writeRaster(r, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  
  
  r_var_prod <- r_in*r
  raster_name_prod <- file.path(out_dir_str, paste(raster_name_tmp,"_",method,"_prod_weights_",out_suffix_str,file_format,sep=""))
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

  removeTmpFiles(h=0) #did not work if h is not set to 0
  ## end of remove section
  
  return(rast_list)
}

raster_match <- function(i,list_param){
  ### Read in parameters/arguments
  lf_files <- list_param$lf_files
  rast_ref <- list_param$rast_ref #name of reference file
  file_format <- list_param$file_format #".tif",".rst" or others
  python_bin <- list_param$python_bin
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
  python_cmd <- file.path(python_bin,"gdalwarp")

  #cmd_str <- paste("/usr/bin/gdalwarp",inFilename,raster_name,sep=" ") #this may be a problem
  cmd_str <- paste(python_cmd,inFilename,raster_name,sep=" ") #this may be a problem
  #gdalwarp -t_srs '+proj=utm +zone=11 +datum=WGS84' raw_spot.tif utm11.tif
  system(cmd_str)
  
  ##return name of file created
  return(raster_name)
}

mosaicFiles <- function(lf_mosaic,mosaic_method="unweighted",num_cores=1,r_mask_raster_name=NULL,python_bin=NULL,mosaic_python="/nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum_noDataTest.py",algorithm="R",match_extent=TRUE,df_points=NULL,NA_flag_val=-9999,file_format=".tif",out_suffix=NULL,out_dir=NULL,tmp_files=FALSE){
  #This functions mosaics tiles/files give a list of files. 
  #There are four options to mosaic:   use_sine_weights,use_edge,use_linear_weights, unweighted
  #Sine weights fits sine fuctions across rows and column producing elliptical/spherical patterns from center
  #Use edge uses the distance from the edge of the tiles/fies, higher weights towards the center
  #Linear weights use simple linear average from distance point feature (usually centroid)
  #Unweighted: average without and weigthing surface
  #In addition, option is given to the user to use R raster mosaic function or a python/gdal modified gdalmerge in mosaicing.
  #
  #INPUT Arguments: 
  #1)lf_mosaic: list of files to mosaic
  #2)mosaic_method: mosaic methods availbable:use_sine_weights,use_edge,use_linear_weights
  #3)num_cores: number of cores used in parallilization in mclapply
  #4)r_mask_raster_name: mask rference raster image
  #5)python_bin: location of python executables, defaut is NULL
  #6)mosaic_python: location/directory of python excecutable used for mosaicing with option sum/mean from Alberto Guzmann
  #7)df_points: point location used in weighting, defaul is NULL
  #8)NA_flag_val: raster flag values, e.g. -9999
  #9)file_format: raster format used, default is ".tif"
  #10)out_suffix: output suffix, default is NULL, it is recommended to add the variable name
  #             here e.g. dailyTmax and date!!
  #11)out_dir: output directory, default is NULL
  #12)algorithm: use R or python function
  #13)match extent: if TRUE match extent before mosaicing
  #14)tmp_files: if TRUE then keep temporary files
  #
  #OUTPUT:
  # Object is produced with 3 components:
  # 1) mean_mosaic: list of raster files from produced mosaic ,
  # 2) r_weights: list of raster files from weights 
  # 3) r_weights_prod: list of raster files from product weights (weights*value)
  # 4) method: weighting average used
  #

  ####################################
  ### BEGIN ####
  
  out_dir_str <- out_dir

  #if(tmp_files==T){
  out_suffix_str_tmp <- paste0(out_suffix,"_tmp")
  #}
  
  lf_r_weights <- vector("list",length=length(lf_mosaic))
  
  ###############
  ### PART 2: prepare weights using tile rasters ############
  #methods availbable:use_sine_weights,use_edge,use_linear_weights

  if(mosaic_method=="use_linear_weights"){
    method <- "use_linear_weights"
    df_points <- df_centroids
    #df_points <- NULL
    r_feature <- NULL

    #lf <- list_param$lf
    #df_points <- list_param$df_points
    #r_feature <- list_param$r_feature #this should be change to a list
    #padding <- TRUE #if padding true then make buffer around edges??
    #method <- list_param$method #differnt methods available to create weights
    #NAflag,file_format,out_suffix etc...
    #NA_flag_val <- list_param$NA_flag
    #file_format <- list_param$file_format
    #out_suffix_str <- list_param$out_suffix_str
    #out_dir_str <- list_param$out_dir_str
    list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,NA_flag_val,file_format,out_suffix_str_tmp,out_dir_str) 
    names(list_param_create_weights) <- c("lf","df_points","r_feature","method","NA_flag","file_format","out_suffix_str","out_dir_str") 
    #num_cores <- 11

    #debug(create_weights_fun)
    #weights_obj <- create_weights_fun(1,list_param=list_param_create_weights)

    #This is the function creating the weights by tile. Distance from the centroids needs to be change from distance to
    #the edges...can use rows and columsn to set edges to 1 and 0 for the others.
    linear_weights_obj_list <- mclapply(1:length(lf_mosaic),FUN=create_weights_fun,list_param=list_param_create_weights,mc.preschedule=FALSE,mc.cores = num_cores)                           

    list_linear_r_weights <- lapply(1:length(linear_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights},x=linear_weights_obj_list)
    list_linear_r_weights_prod <- lapply(1:length(linear_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights_prod},x=linear_weights_obj_list)

    list_weights <- list_linear_r_weights
    list_weights_prod <- list_linear_r_weights_prod 

  }
  if(mosaic_method=="use_sine_weights"){
    
    ### Third use sine weights
    method <- "use_sine_weights"
    #df_points <- df_centroids
    df_points <- NULL
    r_feature <- NULL

    #list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,out_dir_str) 
    #names(list_param_create_weights) <- c("lf","df_points","r_feature","method","out_dir_str") 
    list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,NA_flag_val,file_format,out_suffix_str_tmp,out_dir_str) 
    names(list_param_create_weights) <- c("lf","df_points","r_feature","method","NA_flag","file_format","out_suffix_str","out_dir_str") 

    #num_cores <- 11

    #debug(create_weights_fun)
    #weights_obj <- create_weights_fun(1,list_param=list_param_create_weights)

    #This is the function creating the weights by tile. Distance from the centroids needs to be change from distance to
    #the edges...can use rows and columsn to set edges to 1 and 0 for the others.
    sine_weights_obj_list <- mclapply(1:length(lf_mosaic),FUN=create_weights_fun,list_param=list_param_create_weights,mc.preschedule=FALSE,mc.cores = num_cores)                           

    list_sine_r_weights <- lapply(1:length(sine_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights},x=sine_weights_obj_list)
    list_sine_r_weights_prod <- lapply(1:length(sine_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights_prod},x=sine_weights_obj_list)

    list_weights <- list_sine_r_weights
    list_weights_prod <- list_sine_r_weights_prod 

  }
  
  if(mosaic_method=="use_edge_weights"){
    #this took 5 minutes for 28 tiles for reg4, South America,  4*28
    
    method <- "use_edge"
    df_points <- NULL
    r_feature <- NULL
    
    list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,NA_flag_val,file_format,out_suffix_str_tmp,out_dir_str) 
    names(list_param_create_weights) <- c("lf","df_points","r_feature","method","NA_flag","file_format","out_suffix_str","out_dir_str") 
    #list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,out_dir_str) 
    #names(list_param_create_weights) <- c("lf","df_points","r_feature","method","out_dir_str") 
    #num_cores <- 11
    #undebug(create_weights_fun)
    #weights_obj <- create_weights_fun(3,list_param=list_param_create_weights)

    #This is the function creating the weights by tile. Distance from the centroids needs to be change from distance to
    #the edges...can use rows and columsn to set edges to 1 and 0 for the others.
    use_edge_weights_obj_list <- mclapply(1:length(lf_mosaic),FUN=create_weights_fun,list_param=list_param_create_weights,mc.preschedule=FALSE,mc.cores = num_cores)                           

    #extract the list of files for weights and product weights
    list_edge_r_weights <- lapply(1:length(use_edge_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights},x=use_edge_weights_obj_list)
    list_edge_r_weights_prod <- lapply(1:length(use_edge_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights_prod},x=use_edge_weights_obj_list)
    
    #simplifly later...
    list_weights <- list_edge_r_weights
    list_weights_prod <- list_edge_r_weights_prod 
    #r_test <- raster(list_edge_r_weights[[1]])

  }
  
  ###############
  ### PART 3: prepare weightsfor mosaicing by matching extent ############

  ## Rasters tiles vary slightly in resolution, they need to be matched for the mosaic. Resolve issue in the 
  #mosaic function using gdal_merge to compute a reference image to mach.
  #This step of creating a merged raster can be avoided if a reference maks image is given
  #this needs to be changed to avoid further bugs!!!
  
  if(!is.null(r_mask_raster_name)){
    rast_ref <- r_mask_raster_name #the mask file is used as ref.
    #mask(raster(r_m_weighted_mean_raster_name),mask=r_mask_raster_name,filename=r_m_weighted_mean_mask_raster_name)
    #raster_name <- r_m_weighted_mean_mask_raster_name
  }else{
    rast_ref <- file.path(out_dir,paste("avg_",out_suffix,file_format,sep="")) #this is a the ref
    if(is.null(python_bin)){
      python_bin=""
    }
    
    python_cmd <- file.path(python_bin,"gdal_merge.py")  
    cmd_str <- paste("python",python_cmd,"-o ",rast_ref,paste(lf_mosaic,collapse=" ")) 
    system(cmd_str)
  }
  
  ## Create raster image for original predicted images with matching resolution and extent to the mosaic (reference image)

  #rast_ref <- file.path(out_dir,"avg.tif")
  r_ref <- raster(rast_ref)
  #plot(r_ref)
  
  if(mosaic_method%in%c("use_linear_weights","use_sine_weights","use_edge_weights")){
    

    #####################
    ###### PART 4: compute the weighted mean with the mosaic function #####

    if(algorithm=="python"){
      
      if(match_extent==TRUE){
        
        #If using R, it is necessary to match extent firt
        lf_files <- unlist(list_weights)

        ##Maching resolution is probably only necessary for the r mosaic function
        #Modify later to take into account option R or python...
        list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix_str_tmp,out_dir_str)
        names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")

        #undebug(raster_match)
        #r_test <- raster_match(1,list_param_raster_match)
        #r_test <- raster(raster_match(1,list_param_raster_match))

        list_weights_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

        lf_files <- unlist(list_weights_prod)
        list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix_str_tmp,out_dir_str)
        names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")

        #num_cores <-11
        list_weights_prod_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

      }else{
        list_weights_m <- list_weights
        list_weights_prod_m <- list_weights_prod 
      }
      
      #The file to do the merge is /nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum.py. Sample call below.
      #python /nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum.py --config GDAL_CACHEMAX=1500 --overwrite=TRUE -o  outputname.tif --optfile input.txt
      #lf_day_to_mosaic <- list_weights_m
      
      #pattern_str <- paste("*.","predicted_mod1",".*.",day_to_mosaic[i],".*.tif",sep="")
      #lf_day_to_mosaic <- lapply(1:length(unlist(in_dir_mosaics)),FUN=function(k){list.files(path=unlist(in_dir_mosaics)[k],pattern=pattern_str,full.names=T,recursive=T)}) 
      #lf_day_to_mosaic <- unlist(lf_day_to_mosaic)
      #write.table(lf_day_to_mosaic,file=file.path(out_dir,paste("list_to_mosaics_",day_to_mosaic[i],".txt",sep="")))
      #filename_list_mosaics <- file.path(out_dir,paste("list_to_mosaics_",day_to_mosaic[i],".txt",sep=""))

      filename_list_mosaics_weights_m <- file.path(out_dir_str,paste("list_to_mosaics_","weights_",mosaic_method,"_",out_suffix_str_tmp,".txt",sep=""))
      filename_list_mosaics_prod_weights_m <- file.path(out_dir_str,paste("list_to_mosaics_","prod_weights_",mosaic_method,"_",out_suffix_str_tmp,".txt",sep=""))
      
      #writeLines(unlist(list_weights_m),con=filename_list_mosaics_weights_m) #weights files to mosaic 
      #writeLines(unlist(list_weights_prod_m),con=filename_list_mosaics_prod_weights_m) #prod weights files to mosaic
      
      writeLines(unlist(list_weights_m),con=filename_list_mosaics_weights_m) #weights files to mosaic 
      writeLines(unlist(list_weights_prod_m),con=filename_list_mosaics_prod_weights_m) #prod weights files to mosaic

      #out_mosaic_name_weights_m <- r_weights_sum_raster_name <- file.path(out_dir,paste("r_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix,".tif",sep=""))
      #out_mosaic_name_prod_weights_m <- r_weights_sum_raster_name <- file.path(out_dir,paste("r_prod_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix,".tif",sep=""))
      out_mosaic_name_weights_m  <- file.path(out_dir_str,paste("r_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix_str_tmp,".tif",sep=""))
      out_mosaic_name_prod_weights_m <- file.path(out_dir_str,paste("r_prod_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix_str_tmp,".tif",sep=""))

      #in_file_to_mosaics <- filename_list_mosaics        
      #in_dir_mosaics <- file.path(in_dir1,region_names[i])
      #out_dir_mosaics <- "/nobackupp6/aguzman4/climateLayers/output1000x3000_km/reg5/mosaicsMean"
      #Can be changed to have mosaics in different dir..
      #out_dir_mosaics <- out_dir
      #prefix_str <- "reg4_1500x4500"
      #tile_size <- basename(dirname(in_dir[[i]]))
      #tile_size <- basename(in_dir1)

      #prefix_str <- paste(region_names[i],"_",tile_size,sep="")
      #mod_str <- "mod1" #use mod2 which corresponds to model with LST and elev
      #out_mosaic_name <- paste(region,"_mosaics_",mod_str,"_",tile_size,"_",day_to_mosaic[i],"_",out_prefix,".tif",sep="")
      
      ## Mosaic sum weights...
      #input_file <- filename_list_mosaics_weights_m
      
      module_path <- mosaic_python #this should be a parameter for the function...

      #python /nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum.py 
      #--config GDAL_CACHEMAX=1500 --overwrite=TRUE -o  outputname.tif --optfile input.txt
      #r_weights_sum_raster_name <- mosaic_python_merge(module_path=mosaic_python,
      #                                                 module_name="gdal_merge_sum.py",
      #                                                 input_file=filename_list_mosaics_weights_m,
      #                                                 out_mosaic_name=out_mosaic_name_weights_m)
      #mosaic_python_merge <- function(NA_flag_val,module_path,module_name,input_file,out_mosaic_name){
      mosaic_weights_obj <- mosaic_python_merge(NA_flag_val=NA_flag_val,
                                                module_path=mosaic_python,
                                                module_name="gdal_merge_sum.py",
                                                input_file=filename_list_mosaics_weights_m,
                                                out_mosaic_name=out_mosaic_name_weights_m)
      r_weights_sum_raster_name <- mosaic_weights_obj$out_mosaic_name
      cmd_str1 <- mosaic_weights_obj$cmd_str
      #r_prod_sum_raster_name <- mosaic_python_merge(module_path=mosaic_python,
      #                                              module_name="gdal_merge_sum.py",
      #                                              input_file=filename_list_mosaics_prod_weights_m,
      #                                              out_mosaic_name=out_mosaic_name_prod_weights_m)
      
      mosaic_prod_weights_obj <- mosaic_python_merge(NA_flag_val=NA_flag_val,
                                                     module_path=mosaic_python,
                                                     module_name="gdal_merge_sum.py",
                                                     input_file=filename_list_mosaics_prod_weights_m,
                                                     out_mosaic_name=out_mosaic_name_prod_weights_m)
      r_prod_sum_raster_name <- mosaic_prod_weights_obj$out_mosaic_name
      cmd_str2 <- mosaic_prod_weights_obj$cmd_str
      #write out python command used for mosaicing
      cmd_mosaic_logfile <- file.path(out_dir,paste("cmd_mosaic_",out_suffix,".txt",sep=""))
      writeLines(cmd_str1,con=cmd_mosaic_logfile) #weights files to mosaic 
      #writeLines(cmd_str2,con=file.path(out_dir,paste("cmd_mosaic_",out_suffix,".txt",sep=""))) #weights files to mosaic 
      cat(cmd_str2, file=cmd_mosaic_logfile, append=TRUE, sep = "\n")
    }
    
    if(algorithm=="R"){
      
      #If using R, it is necessary to match extent firt
      
      lf_files <- unlist(list_weights)

      ##Maching resolution is probably only necessary for the r mosaic function
      #MOdify later to take into account option R or python...
      list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix_str_tmp,out_dir_str)
      names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")

      #undebug(raster_match)
      #r_test <- raster_match(1,list_param_raster_match)
      #r_test <- raster(raster_match(1,list_param_raster_match))

      list_weights_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

      lf_files <- unlist(list_weights_prod)
      list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix_str_tmp,out_dir_str)
      names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")

      #num_cores <-11
      list_weights_prod_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

      
      #The file to do the merge is /nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum.py. Sample call below.
      #python /nobackupp6/aguzman4/climateLayers/sharedCode/gdal_merge_sum.py --config GDAL_CACHEMAX=1500 --overwrite=TRUE -o  outputname.tif --optfile input.txt
          ##Make this a function later
      #list_weights_m <- list(list_linear_weights_m,list_edge_weights_m,list_sine_weights_m)
      #list_weights_prod_m <- list(list_linear_weights_prod_m,list_edge_weights_prod_m,list_sine_weights_prod_m)
      #list_methods <- c("linear","edge","sine")
      list_mosaiced_files <- vector("list",length=1)

      list_args_weights <- list_weights_m
      list_args_weights_prod <- list_weights_prod_m
      method_str <- method
  
      #making a list of raster object before mosaicing
      list_args_weights <- lapply(1:length(list_args_weights), FUN=function(i,x){raster(x[[i]])},x=list_args_weights)

      #get the list of weights product into raster objects

      list_args_weights_prod <- lapply(1:length(list_args_weights_prod), FUN=function(i,x){raster(x[[i]])},x=list_args_weights_prod)
      list_args_weights_prod$fun <- "sum" #use sum while mosaicing
      list_args_weights_prod$na.rm <- TRUE #deal with NA by removal
      r_weights_sum_raster_name <- file.path(out_dir_str,paste("r_weights_sum_m_",method_str,"_weighted_mean_",out_suffix_str_tmp,".tif",sep=""))
      list_args_weights$filename <- r_weights_sum_raster_name
      list_args_weights$overwrite<- TRUE
      list_args_weights_prod$overwrite<- TRUE  #add to overwrite existing image  
    
      list_args_weights$fun <- "sum" #we want the sum to compute the weighted mean
      list_args_weights$na.rm <- TRUE
      r_prod_sum_raster_name <- file.path(out_dir_str,paste("r_prod_sum_m_",method_str,"_weighted_mean_",out_suffix_str_tmp,".tif",sep=""))
      list_args_weights_prod$filename <- r_prod_sum_raster_name

      #Mosaic files: this is where we can use Alberto Python function but modified with option for
      #sum in addition ot the current option for mean.
      #This took 23 minutes!
      r_weights_sum <- do.call(mosaic,list_args_weights) #weights sum image mosaiced
      #This took 23 minutes with the R function
      r_prod_sum <- do.call(mosaic,list_args_weights_prod) #weights sum product image mosacied

    }
    

    #r_m_weighted_mean <- r_prod_sum/r_weights_sum #this is the mosaic using weighted mean...

    r_m_weighted_mean_raster_name <- file.path(out_dir_str,paste("r_m_",mosaic_method,"_weighted_mean_",out_suffix,".tif",sep=""))

    if(is.null(python_bin)){
      python_bin=""
    }
    
    python_cmd <- file.path(python_bin,"gdal_calc.py")
    cmd_str <- paste(python_cmd, 
                     paste("-A ", r_prod_sum_raster_name,sep=""),
                     paste("-B ", r_weights_sum_raster_name,sep=""),
                     paste("--outfile=",r_m_weighted_mean_raster_name,sep=""),
                     paste("NoDataValue=",NA_flag_val,sep=""),
                     "--calc='A/B'","--overwrite",sep=" ") #division by zero is problematic...
    system(cmd_str)
    cmd_mosaic_logfile <- file.path(out_dir,paste("cmd_mosaic_",out_suffix,".txt",sep=""))
    writeLines(cmd_str1,con=cmd_mosaic_logfile) #weights files to mosaic 
    #writeLines(cmd_str2,con=file.path(out_dir,paste("cmd_mosaic_",out_suffix,".txt",sep=""))) #weights files to mosaic 
    cat(cmd_str2, file=cmd_mosaic_logfile, append=TRUE, sep = "\n")

    #cmd_str <- "/nobackupp6/aguzman4/climateLayers/sharedModules/bin/gdal_calc.py -A r_prod_weights_sum_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920101_reg4_run10_1500x4500_global_analyses_pred_1992_10052015.tif -B r_weights_sum_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19920101_reg4_run10_1500x4500_global_analyses_pred_1992_10052015.tif --outfile='test2.tif' --calc='A/B' --overwrite"
    
    #writeRaster(r_m_weighted_mean, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  
    #now use the mask
    if(!is.null(r_mask_raster_name)){
      #different extent between mask and output if match extent is false!!
      #match resolution and extent first
      
      lf_files <- c(r_m_weighted_mean_raster_name) #file(s) to be matched
      rast_ref <- r_mask_raster_name
      ##Maching resolution is probably only necessary for the r mosaic function
      #Modify later to take into account option R or python...
      list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix,out_dir)
      names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")

      #undebug(raster_match)
      r_m_weighted_mean_raster_name_matched <- raster_match(1,list_param_raster_match)

      r_m_weighted_mean_mask_raster_name <- file.path(out_dir_str,paste("r_m_",mosaic_method,"_weighted_mean_mask_",out_suffix,".tif",sep=""))
      mask(raster(r_m_weighted_mean_raster_name_matched),mask=raster(r_mask_raster_name),
           filename=r_m_weighted_mean_mask_raster_name,overwrite=TRUE)
      raster_name <- r_m_weighted_mean_mask_raster_name
    }else{
      raster_name <- r_m_weighted_mean_raster_name
    }
  }

  if(mosaic_method=="unweighted"){
    #### Fourth use original images
    #macth file to mosaic extent using the original predictions
    
    if(match_extent==TRUE){
      lf_files <- lf_mosaic
      list_param_raster_match <- list(lf_files,rast_ref,file_format,out_suffix,out_dir)
      names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","out_suffix","out_dir_str")
      list_pred_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           
    }else{
      list_pred_m <- lf_mosaic
    }
    #list_mosaiced_files <- list.files(pattern="r_m.*._weighted_mean_.*.tif")

    #names(list_mosaiced_files) <- c("edge","linear","sine")
  
    #### NOW unweighted mean mosaic

    #get the original predicted image to raster (matched previously to the mosaic extent)
    list_args_pred_m <- list_pred_m
    #list_args_pred_m <- (mixedsort(list.files(pattern="^gam_CAI.*.m_mosaic_run10_1500x4500_global_analyses_03252015.tif")))
    list_args_pred_m <- lapply(1:length(list_args_pred_m), FUN=function(i,x){raster(x[[i]])},x=list_args_pred_m)

    list_args_pred_m$fun <- "mean"
    list_args_pred_m$na.rm <- TRUE
    list_args_pred_m$overwrite<- TRUE  #add to overwrite existing image 
    #list_args_pred_m$filename <- 

    #Mosaic files using R raster mosaic 
    r_m_mean <- do.call(mosaic,list_args_pred_m) #this is unweighted mean from the predicted raster

    r_m_mean_raster_name <- file.path(out_dir,paste("r_m_mean_",out_suffix,".tif",sep=""))
    writeRaster(r_m_mean, NAflag=NA_flag_val,filename=r_m_mean_raster_name,overwrite=TRUE)  #unweighted mean

    #r_m_mean_unweighted <- paste("r_m_mean_",out_suffix,".tif",sep="")
    #list_weights <- NULL
    #list_weights_prod <- NULL

    
    if(!is.null(r_mask_raster_name)){
      #different extent between mask and output if match extent is false!!
      #match resolution and extent first
      
      lf_files <- c(r_m_mean_raster_name) #file(s) to be matched
      rast_ref <- r_mask_raster_name
      ##Maching resolution is probably only necessary for the r mosaic function
      #Modify later to take into account option R or python...
      list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix,out_dir)
      names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")

      #undebug(raster_match)
      r_m_mean_raster_name_matched <- raster_match(1,list_param_raster_match)

      r_m_mean_mask_raster_name <- file.path(out_dir,paste("r_m_",method_str,"_unweighted_mean_mask_",out_suffix,".tif",sep=""))
      mask(raster( r_m_mean_raster_name_matched),mask=raster(r_mask_raster_name),
           filename=r_m_mean_mask_raster_name,overwrite=TRUE)
      raster_name <- r_m_mean_mask_raster_name
    }else{
      raster_name <- r_m_mean_raster_name
    }

  }
  
  ########## clean up the disk/directories before ending the function ####
  
  if(tmp_files==F){ #if false...delete all files with "_tmp"
    lf_tmp <- list.files(pattern="*.*tmp*.*",path=out_dir_str,full.names=T)
    ##now delete temporary files...
    file.remove(lf_tmp)
  }
  
  #Create return object
  mosaic_obj <- list(raster_name,list_weights,list_weights_prod,mosaic_method)
  names(mosaic_obj) <- c("mean_mosaic","r_weights","r_weigths_prod","method")
  save(mosaic_obj,file=file.path(out_dir_str,paste(mosaic_method,"_","mosaic_obj_",out_suffix,".RData",sep="")))
  return(mosaic_obj)
}

plot_mosaic <- function(i,list_param){
  #Plot for mosaic test
  #Inputs:
  #method_str: method used in mosaicing
  #lf_mosaic: list of raster files from mosaicing
  #out_suffix: output suffix 
  
  method_str <- list_param$method[i]
  f_mosaic <- list_param$lf_mosaic[i]
  out_suffix_str <- list_param$out_suffix[i]
  
  r_mosaic <- raster(f_mosaic)

  r_mosaic_terrain <- terrain(r_mosaic,opt=c("slope","aspect"),unit="degrees")

  res_pix <- 1200
  col_mfrow <- 1 
  row_mfrow <- 0.8
  
  out_file1 <- paste("Figure2_mosaic_mean_",method_str,"_",out_suffix_str,".png",sep="")
  png(filename= out_file1,width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_mosaic,main=paste("mosaic mean ",method_str,sep=""))
  
  dev.off()
  
  #### plot terrain to emphasize possible edges..
  res_pix <- 1200
  col_mfrow <- 1 
  row_mfrow <- 0.8

  out_file2 <- paste("Figure2_slope_mean_",method_str,"_",out_suffix_str,".png",sep="")
  png(filename= out_file2,
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_mosaic_terrain,y=1,main=paste("slope mosaic mean ",method_str,sep=""))
  
  dev.off()

  out_file3 <- paste("Figure2_aspect_mean_",method_str,"_",out_suffix_str,".png",sep="")
  png(filename= out_file3,
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_mosaic_terrain,y=2,main=paste("aspect mean ",method_str,sep=""))
  title(paste("aspect mean ",method_str,sep=""))
  dev.off()
  
  l_out_files <- list(out_file1,out_file2,out_file3)
  return(l_out_files)
}

plot_diff_raster <- function(i,list_param){
  #Plot for mosaic differences
  #Inputs:
  #lf1: list of raster files used as reference
  #lf2: list of raster files used as second image
  #out_suffix: output suffix
  
  ### Read in parameters
  #method_str <- list_param$method[i]
  f_r1 <- list_param$lf1[i] #e.g. unweighted
  f_r2<- list_param$lf2[i] #e.g. weighted
  out_suffix_str <- list_param$out_suffix[i]
  
  ### BEGIN ####
  
  r1 <- raster(f_r1)
  r2 <- raster(f_r2)

  r_diff_raster <- r1 - r2

  res_pix <- 1200
  col_mfrow <- 1 
  row_mfrow <- 0.8
  
  out_file <- paste("Figure2_diff_raster","_",out_suffix_str,".png",sep="")
  png(filename = out_file,
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_diff_raster,main=out_suffix_str)

  dev.off()
  
  return(out_file)
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
  
  png(filename=file.path(out_dir_str,
                         paste("Figure9_clim_mosaics_day_test","_",date_proc,"_",reg_name,"_",out_suffix,".png",sep="")),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_pred,main=paste("Predicted on ",date_proc , " ", reg_name,sep=""),cex.main=1.5)
  dev.off()
  
  return(raster_name)
  
}

#Use this instead of daily mosaic plot function
#Add NAvalue flag!!
plot_screen_raster_val<-function(i,list_param){
  ##USAGE###
  #Screen raster list and produced plot as png.
  fname <-list_param$lf_raster_fname[i]
  var_name <-list_param$var_name #tmax, rmse tmax etc.  
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
  png_filename <- paste(prefix_str,"_",date_proc,"_","_",out_suffix_str,".png",sep="")
  png(filename=png_filename ,
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_pred,main=paste("Predicted on ",date_proc ," ",var_name, sep=""),cex.main=1.5)
  dev.off()
  
  return(png_filename)
}

## functions from kriging...

fit_models<-function(list_formulas,data_training){
  #This functions several models and returns model objects.
  #Arguments: - list of formulas for GAM models
  #           - fitting data in a data.frame or SpatialPointDataFrame
  #Output: list of model objects 
  list_fitted_models<-vector("list",length(list_formulas))
  for (k in 1:length(list_formulas)){
    formula<-list_formulas[[k]]
    mod<- try(gam(formula, data=data_training)) #change to any model!!
    #mod<- try(autoKrige(formula, input_data=data_s,new_data=s_sgdf,data_variogram=data_s))
    model_name<-paste("mod",k,sep="")
    assign(model_name,mod) 
    list_fitted_models[[k]]<-mod
  }
  return(list_fitted_models) 
}

select_var_stack <-function(r_stack,formula_mod,spdf=TRUE){
  ##Write function to return only the relevant layers!!
  #Note that default behaviour of the function is to remove na values in the subset 
  #of raster layers and return a spdf
  
  ### Start
  
  covar_terms<-all.vars(formula_mod) #all covariates terms...+ y_var
  if (length(covar_terms)==1){
    r_stack_covar<-subset(r_stack,1)
  } #use one layer
  if (length(covar_terms)> 1){
    r_stack_covar <-subset(r_stack,covar_terms[-1])
  }
  if (spdf==TRUE){
    s_sgdf<-as(r_stack_covar,"SpatialGridDataFrame") #Conversion to spatial grid data frame, only convert the necessary layers!!
    s_spdf<-as.data.frame(s_sgdf) #Note that this automatically removes all NA rows
    s_spdf<-na.omit(s_spdf) #removes all rows that have na...
    coords<- s_spdf[,c('s1','s2')]
    coordinates(s_spdf)<-coords
    proj4string(s_spdf)<-proj4string(s_sgdf)  #Need to assign coordinates...
    #raster_pred <- rasterize(s_spdf,r1,"pred",fun=mean)
    covar_obj<-s_spdf
  } else{
    covar_obj<-r_stack_covar
  }
  
  return(covar_obj)
}

remove_na_spdf<-function(col_names,d_spdf){
  #Purpose: remote na items from a subset of a SpatialPointsDataFrame
  x<-d_spdf
  coords <-coordinates(x)
  x$s1<-coords[,1]
  x$s2<-coords[,2]
  
  x1<-x[c(col_names,"s1","s2")]
  #x1$y_var <-data_training$y_var
  #names(x1)
  x1<-na.omit(as.data.frame(x1))
  coordinates(x1)<-x1[c("s1","s2")]
  proj4string(x1)<-proj4string(d_spdf)
  return(x1)
}

predict_auto_krige_raster_model<-function(list_formulas,r_stack,data_training,out_filename){
  #This functions performs predictions on a raster grid given input models.
  #Arguments: list of fitted models, raster stack of covariates
  #Output: spatial grid data frame of the subset of tiles
  
  list_fitted_models<-vector("list",length(list_formulas))
  list_rast_pred<-vector("list",length(list_formulas))
  #s_sgdf<-as(r_stack,"SpatialGridDataFrame") #Conversion to spatial grid data frame, only convert the necessary layers!!
  proj4string(data_training) <- projection(r_stack)
  for (k in 1:length(list_formulas)){
    formula_mod<-list_formulas[[k]]
    raster_name<-out_filename[[k]]
    #mod<- try(gam(formula, data=data_training)) #change to any model!!
    s_spdf<-select_var_stack(r_stack,formula_mod,spdf=TRUE)
    col_names<-all.vars(formula_mod)
    if (length(col_names)==1){
      data_fit <-data_training
    }else{
      data_fit <- remove_na_spdf(col_names,data_training)
    }
    #use modify version of autokrige called autokrige_fun
    mod <- try(autoKrige_fun(formula_mod, input_data=data_fit,new_data=s_spdf,data_variogram=data_fit))
    #mod <- try(autoKrige(formula_mod, input_data=data_training,new_data=s_spdf,data_variogram=data_training))
    model_name<-paste("mod",k,sep="")
    assign(model_name,mod) 
    
    if (inherits(mod,"autoKrige")) {           #change to c("gam","autoKrige")
      rpred<-mod$krige_output  #Extracting the SptialGriDataFrame from the autokrige object
      y_pred<-rpred$var1.pred                  #is the order the same?
      raster_pred <- rasterize(rpred,r_stack,"var1.pred",fun=mean)
      names(raster_pred)<-"y_pred" 
      writeRaster(raster_pred, filename=raster_name,overwrite=TRUE)  #Writing the data in a raster file format...
      #print(paste("Interpolation:","mod", j ,sep=" "))
      list_rast_pred[[k]]<-raster_name
      mod$krige_output<-NULL
      list_fitted_models[[k]]<-mod
      
    }
    if (inherits(mod,"try-error")) {
      print(paste("no autokrige model fitted:",mod,sep=" ")) #change message for any model type...
      list_fitted_models[[k]]<-mod
    }
  }
  day_prediction_obj <-list(list_fitted_models,list_rast_pred)
  names(day_prediction_obj) <-c("list_fitted_models","list_rast_pred")
  return(day_prediction_obj)
}

###MODIFIED AUTOKRIGE function from automap package
autoKrige_fun <- function (formula, input_data, new_data, data_variogram = input_data, 
    block = 0, model = c("Sph", "Exp", "Gau", "Ste"), kappa = c(0.05, 
        seq(0.2, 2, 0.1), 5, 10), fix.values = c(NA, NA, NA), 
    remove_duplicates = TRUE, verbose = FALSE, GLS.model = NA, 
    start_vals = c(NA, NA, NA), miscFitOptions = list(), ...) 
{
    if (inherits(formula, "SpatialPointsDataFrame")) {
        input_data = formula
        formula = as.formula(paste(names(input_data)[1], "~ 1"))
    }
    if (!inherits(input_data, "SpatialPointsDataFrame") | !inherits(data_variogram, 
        "SpatialPointsDataFrame")) {
        stop(paste("\nInvalid input objects: input_data or data_variogram not of class 'SpatialPointsDataFrame'.\n\tClass input_data: '", 
            class(input_data), "'", "\n\tClass data_variogram: '", 
            class(data_variogram), "'", sep = ""))
    }
    if (as.character(formula)[3] != 1 & missing(new_data)) 
        stop("If you want to use Universal Kriging, new_data needs to be specified \n  because the predictors are also required on the prediction locations.")
    if (remove_duplicates) {
        zd = zerodist(input_data)
        if (length(zd) != 0) {
            warning("Removed ", length(zd)/2, " duplicate observation(s) in input_data:", 
                immediate. = TRUE)
            print(input_data[c(zd), ])
            input_data = input_data[-zd[, 2], ]
        }
    }
    col_name = as.character(formula)[2]
    if (length(unique(input_data[[col_name]])) == 1) 
        stop(sprintf("All data in attribute '%s' is identical and equal to %s\n   Can not interpolate this data", 
            col_name, unique(input_data[[col_name]])[1]))
    if (missing(new_data)) 
        new_data = create_new_data(input_data)
    p4s_obj1 = proj4string(input_data)
    p4s_obj2 = proj4string(new_data)
    if (!all(is.na(c(p4s_obj1, p4s_obj2)))) {
        if (is.na(p4s_obj1) & !is.na(p4s_obj2)) 
            proj4string(input_data) = proj4string(new_data)
        if (!is.na(p4s_obj1) & is.na(p4s_obj2)) 
            proj4string(new_data) = proj4string(input_data)
        #if (any(!c(is.projected(input_data), is.projected(new_data)))) 
        #    stop(paste("Either input_data or new_data is in LongLat, please reproject.\n", 
        #        "  input_data: ", p4s_obj1, "\n", "  new_data:   ", 
        #        p4s_obj2, "\n"))
        if (proj4string(input_data) != proj4string(new_data)) 
            stop(paste("Projections of input_data and new_data do not match:\n", 
                "  input_data: ", p4s_obj1, "\n", "  new_data:    ", 
                p4s_obj2, "\n"))
    }
    variogram_object = autofitVariogram(formula, data_variogram, 
        model = model, kappa = kappa, fix.values = fix.values, 
        verbose = verbose, GLS.model = GLS.model, start_vals = start_vals, 
        miscFitOptions = miscFitOptions)
    krige_result = krige(formula, input_data, new_data, variogram_object$var_model, 
        block = block, ...)
    krige_result$var1.stdev = sqrt(krige_result$var1.var)
    result = list(krige_output = krige_result, exp_var = variogram_object$exp_var, 
        var_model = variogram_object$var_model, sserr = variogram_object$sserr)
    class(result) = c("autoKrige", "list")
    return(result)
}

predict_accuracy_raster_by_station <- function(var_pred,ref_rast,data_training,out_filename,out_suffix_str,out_dir,mask_file=NULL){
  #
  coord_xy<-coordinates(data_training)
  var_vals <- data_training[[var_pred]]
  fitKrig <- Krig(coord_xy,var_vals)#,theta=1e5) #use TPS or krige 
  #fitKrig <- Krig(coord_xy,var_vals,theta=1e5) #use TPS or krige 
  #mod_krtmp1<-fitKrig
  #model_name<-"mod_kr"
  #options(scipen=999)
  krig_rast <- try(interpolate(ref_rast,fitKrig)) #interpolation using function from raster package

  #Write out modeled layers
    
  if(is.null(mask_file)){
    writeRaster(krig_rast, NAflag=NA_flag_val,filename=out_filename,overwrite=TRUE)  
  }
  if(!is.null(mask_file)){
    #modify here later
    writeRaster(krig_rast, NAflag=NA_flag_val,filename=out_filename,overwrite=TRUE)  
  }

  ### Preparing return object
  kriging_obj <- list(out_filename,fitKrig)
  names(kriging_obj)<-c("raster_name","mod_obj")
  #save(kriging_obj,file= file.path(out_dir,paste("kriging_obj","_",out_suffix_str,".RData",sep="")))
  return(kriging_obj)
}

create_accuracy_residuals_raster <- function(i,list_param){

  #create_accuracy_residuals_raster <- function(i,lf_day_tiles,data_df,df_tile_processed,var_pred,list_models,date_processed,num_cores,NA_flag_val,file_format,out_dir_str,out_suffix_str){
  #This function generates surface for residuals values for a giving set of formula and data input.
  #The method used is currently kriging with two options: automap/gstat and Fields packages.
  #Inputs:
  #lf_day_tiles: list of raster files
  #df_tile_processed: processed tiles
  #data_df: data.frame table/spdf containing stations with residuals and variable
  #var_pred: variable selected to be mapped using modeling (kriging)
  #list_models: formula for the modeling (kriging) e.g. useg by autokrige
  #y_var_name: variable being modeled e.g."dailyTmax", dailyTmin, precip  
  #interpolation_method: names of the interpolation/modeling method
  #date_processed: day being processed , e.g. 19920101, can be month too!!!
  #num_cores : number of cores used in the parallelization
  #NA_flag_val: value used as flag in the raster 
  #file_format: e.g. tif., .rst
  #out_dir_str: output directory
  #out_suffix_str: output suffix
  #Outputs:
  #raster list of resdiuals surfaces and associated modeles by tiles and for one date.
  #
  #TODO: 

  # - automap/gstat for data with projection
  # - clean up
  #

  
  #### FUNCTIONS USED #####
  generate_residuals_raster <- function(j,list_param){
    ##Add documentation here later...
    
    ### PARSE ARGUMENTS ##
    
    lf <- list_param$lf
    var_pred <- list_param$var_pred
    data_df <- list_param$data_df
    df_raster_pred_tiles <- list_param$df_raster_pred_tiles
    list_formulas <- list_param$list_formulas
    use_autokrige <- list_param$use_autokrige
    NA_flag_val <- list_param$NA_flag_val
    file_format <- list_param$file_format
    out_dir_str <- list_param$out_dir_str
    out_suffix_str <- list_param$out_suffix_str

    ###### START SCRIPT
    
    #list_pred_res_obj <-vector("list",length=length(lf))
    #for(j in 1:length(lf)){

    inFilename <- lf[j]

    ref_rast <- raster(inFilename)
    #out_filename <- "test.tif"
    #create output name for predicted raster
    extension_str <- extension(inFilename)
    raster_name_tmp <- gsub(extension_str,"",basename(inFilename))
    out_filename <- file.path(out_dir_str,paste(raster_name_tmp,"_","kriged_residuals_",var_pred,"_",out_suffix_str,file_format,sep="")) #for use in function later...

    #tile_selected <- as.character(df_raster_pred_tiles$tile_id[j])
    data_df$tile_id <- as.character(data_df$tile_id)
    tile_selected <- df_raster_pred_tiles$tile_id[j]
    data_training <- subset(data_df,tile_id==tile_selected) #df_raster_pred_tiles$files
    data_training <- data_training[!is.na(data_training[[var_pred]]),] #screen NA in the independent var
    
    if(use_autokrige==TRUE){
      #this is still under development since there was an error message!!
      r_stack <- stack(inFilename)
      #debug(predict_auto_krige_raster_model) #this calls other function to clean up the inputs
      #data_training_tmp <- idw(zinc ~ 1, meuse2[!is.na(meuse2$zinc),],newdata= meuse.grid)
      
      pred_res_obj <- predict_auto_krige_raster_model(list_formulas,r_stack,data_training,out_filename)
      #mod <- try(autoKrige(formula_mod, input_data=data_fit,new_data=s_spdf,data_variogram=data_fit))
      #Error in autoKrige(formula_mod, input_data = data_fit, new_data = s_spdf,  : 
      #Either input_data or new_data is in LongLat, please reproject.
      #input_data:  +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 
      # new_data:    +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 
      #Problems using autoKrige on non projected data!!! it looks like the fields package is needed...
      #Modified autokrige function...but new error:
      #Error : dimensions do not match: locations 264 and data 125
    }
  
    if(use_autokrige==FALSE){
      #New function with the Fields package...
      #debug(predict_accuracy_raster_by_station)
      #pred_res_obj <- predict_accuracy_raster_by_station(var_pred,ref_rast,data_training,out_filename,out_suffix_str,out_dir,mask_file=NULL)
      pred_res_obj <- predict_accuracy_raster_by_station(var_pred,ref_rast,data_training,out_filename,out_suffix_str,out_dir_str,mask_file=NULL)
                                
    }  
    return(pred_res_obj)
  }

  ####### PARSE ARGUMENTS
  

  #lf <- list_param$lf[[i]] #list of files to mosaic
  lf_day_tiles <- list_param$lf_day_tiles[[i]] #list of raster files
  data_df <- list_param$data_df # data.frame table/spdf containing stations with residuals and variable
  df_tile_processed_reg <- list_param$df_tile_processed_reg #tiles processed during assessment usually by region
  var_pred <- list_param$var_pred #variable being modeled
  list_models <- list_param$list_models #formula for the modeling 
  use_autokrige <- list_param$use_autokrige #if TRUE use automap/gstat package
  y_var_name <- list_param$y_var_name #"dailyTmax" #PARAM2
  interpolation_method <- list_param$interpolation_method #c("gam_CAI") #PARAM3
  date_processed <- list_param$days_to_process[i] #can be a monthly layer
  num_cores <- list_param$num_cores #number of cores used
  NA_flag_val <- list_param$NA_flag_val
  #NAflag,file_format,out_suffix etc...
  file_format <- list_param$file_format
  out_dir_str <- list_param$out_dir_str
  out_suffix_str <- list_param$out_suffix_str
  
  ######## START SCRIPT ###############
  
  list_formulas <- lapply(list_models,as.formula,env=.GlobalEnv) #mulitple arguments passed to lapply!!

  #data_training <- data_df 
  coordinates(data_df)<-cbind(data_df$x,data_df$y)
  data_df <- subset(data_df,date==date_processed) #select the date...
  #lf_day_tiles <- lf_mosaic[[i]] #i index for time...can be monthly or daily time steps??, ok
  
  #Now match the correct tiles with data used in kriging...
  #match the correct tile!!! df_tile_processed
  #pattern_str <- as.character(unique(df_tile_processed$tile_coord))
  
  #check that all the rows are tile related (this is related to the bug of "output_test)
  df_tile_processed_reg <- df_tile_processed_reg[!is.na(df_tile_processed_reg$shp_files),]
  
  list_tile_coord <- as.character(df_tile_processed_reg$tile_coord)
  pattern_str <- glob2rx(paste("*",list_tile_coord,"*","*.tif",sep=""))
  keywords_str <- pattern_str
  tmp_str2 <-unlist(lapply(keywords_str,grep,lf_day_tiles,value=TRUE))
  list_coord_tf <- basename(dirname(dirname(tmp_str2)))
  df_raster_pred_tiles_tmp <- data.frame(files =tmp_str2, tile_coord=list_coord_tf)
   
  #df_raster_pred_tiles_tmp <- data.frame(files =tmp_str2, tile_coord=list_tile_coord)
  df_raster_pred_tiles <- merge(df_raster_pred_tiles_tmp,df_tile_processed_reg,by="tile_coord")
  df_raster_pred_tiles$path_NEX <- as.character(df_raster_pred_tiles$path_NEX)
  df_raster_pred_tiles$reg <- basename(dirname(df_raster_pred_tiles$path_NEX))
  df_raster_pred_tiles$files <- as.character(df_raster_pred_tiles$files)
  df_raster_pred_tiles$tile_id <- as.character(df_raster_pred_tiles$tile_id)
  
  #identify residuals
  
  #call kriging function
  
  ###This can be a new function here with mclapply!!!
  ## Addtional loop...
  #j <- 1 #loops across tiles from set of files/tiles

  lf <- df_raster_pred_tiles$files
  
  ##Make this loop a function later on, testing right now
  list_param_generate_residuals_raster <- list(lf,var_pred,data_df,df_raster_pred_tiles,list_formulas,use_autokrige,NA_flag_val,file_format,out_dir_str,out_suffix_str)
  names(list_param_generate_residuals_raster) <- c("lf","var_pred","data_df","df_raster_pred_tiles","list_formulas","use_autokrige","NA_flag_val","file_format","out_dir_str","out_suffix_str")

  #debug(generate_residuals_raster)
  #test_lf <- lapply(3,FUN=generate_residuals_raster,list_param=list_param_generate_residuals_raster)                           
  
  list_pred_res_obj <- mclapply(1:length(lf),FUN=generate_residuals_raster,list_param=list_param_generate_residuals_raster,mc.preschedule=FALSE,mc.cores = num_cores)                           
  ## Add to df_raster_pred_tiles
  
  
  #write output
  accuracy_residuals_obj <-list(list_pred_res_obj,data_df,df_raster_pred_tiles)
  names(accuracy_residuals_obj)<-c("list_pred_res_obj","data_df","df_raster_pred_tiles")
  save(accuracy_residuals_obj,file= file.path(out_dir_str,paste("accuracy_residuals_obj_",date_processed,"_",var_pred,
                                                            out_suffix_str,".RData",sep="")))
  
  return(accuracy_residuals_obj) 
}

##Also found in accuracy assessment script:global_run_scalingup_assessment_part1_functions_02112015.R
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

##################### END OF SCRIPT ######################
