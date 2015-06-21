##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Different options to explore mosaicing are tested. This script only contains functions.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 04/14/2015  
#MODIFIED ON: 06/21/2015            
#Version: 1
#PROJECT: Environmental Layers project     
#COMMENTS: first commit of function script to test mosaicing using 1500x4500km and other tiles
#TODO:
#1) Make this is a script/function callable from the shell/bash
#2) Improve performance: there will be a need to improve efficiency for the workflow.

#Functions: available:
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
  tile_no <- i
  
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
      srcfile <- file.path(out_dir,paste("feature_target_",tile_no,".tif",sep=""))

      writeRaster(r_init,filename=srcfile,overwrite=T)
      dstfile <- file.path(out_dir,paste("feature_target_edge_distance",tile_no,".tif",sep=""))
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

mosaicFiles <- function(lf_mosaic,mosaic_method,num_cores,python_bin=NULL,df_points=NULL,NA_flag_val=-9999,file_format=".tif",out_suffix=NULL,out_dir=NULL){
  #Add documentation!!!!
  ##
  
  ### BEGIN ####
  out_dir_str <- out_dir

  lf_r_weights <- vector("list",length=length(lf_mosaic))
  
  ###############
  ### PART 2: prepare weights using tile rasters ############
  #methods availbable:use_sine_weights,use_edge,use_linear_weights

  if(mosaic_method=="use_linear_weights"){
    method <- "use_linear_weights"
    df_points <- df_centroids
    #df_points <- NULL
    r_feature <- NULL

    list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,out_dir_str) 
    names(list_param_create_weights) <- c("lf","df_points","r_feature","method","out_dir_str") 
    #num_cores <- 11

    #debug(create_weights_fun)
    weights_obj <- create_weights_fun(1,list_param=list_param_create_weights)

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

    list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,out_dir_str) 
    names(list_param_create_weights) <- c("lf","df_points","r_feature","method","out_dir_str") 
    #num_cores <- 11

    #debug(create_weights_fun)
    weights_obj <- create_weights_fun(1,list_param=list_param_create_weights)

    #This is the function creating the weights by tile. Distance from the centroids needs to be change from distance to
    #the edges...can use rows and columsn to set edges to 1 and 0 for the others.
    sine_weights_obj_list <- mclapply(1:length(lf_mosaic),FUN=create_weights_fun,list_param=list_param_create_weights,mc.preschedule=FALSE,mc.cores = num_cores)                           

    list_sine_r_weights <- lapply(1:length(sine_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights},x=sine_weights_obj_list)
    list_sine_r_weights_prod <- lapply(1:length(sine_weights_obj_list), FUN=function(i,x){x[[i]]$r_weights_prod},x=sine_weights_obj_list)

    list_weights <- list_sine_r_weights
    list_weights_prod <- list_sine_r_weights_prod 

  }
  
  if(mosaic_method=="use_edge_weights"){
    
    method <- "use_edge"
    df_points <- NULL
    r_feature <- NULL
    list_param_create_weights <- list(lf_mosaic,df_points,r_feature,method,out_dir_str) 
    names(list_param_create_weights) <- c("lf","df_points","r_feature","method","out_dir_str") 
    #num_cores <- 11

    weights_obj <- create_weights_fun(1,list_param=list_param_create_weights)

    #This is the function creating the weights by tile. Distance from the centroids needs to be change from distance to
    #the edges...can use rows and columsn to set edges to 1 and 0 for the others.
    use_edge_weights_obj_list <- mclapply(1:length(lf_mosaic),FUN=create_weights_fun,list_param=list_param_create_weights,mc.preschedule=FALSE,mc.cores = num_cores)                           

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
  #mosaic funciton using gdal_merge to compute a reference image to mach.

  rast_ref <- file.path(out_dir,paste("avg_",out_suffix,file_format,sep="")) #this is a the ref

  cmd_str <- paste("python","/usr/bin/gdal_merge.py","-o ",rast_ref,paste(lf_mosaic,collapse=" ")) 
  system(cmd_str)

  ## Create raster image for original predicted images with matching resolution and extent to the mosaic (reference image)

  #rast_ref <- file.path(out_dir,"avg.tif")
  r_ref <- raster(rast_ref)
  #plot(r_ref)
  
  if(mosaic_method%in%c("use_linear_weights","use_sine_weights","use_edge_weights")){
    
    lf_files <- unlist(list_weights)

    list_param_raster_match <- list(lf_files,rast_ref,file_format,out_suffix,out_dir)
    names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","out_suffix","out_dir_str")

    #debug(raster_match)
    #r_test <- raster(raster_match(1,list_param_raster_match))

    list_weights_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

    lf_files <- unlist(list_weights_prod)
    list_param_raster_match <- list(lf_files,rast_ref,file_format,out_suffix,out_dir)
    names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","out_suffix","out_dir_str")

    #num_cores <-11
    list_weights_prod_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

    #####################
    ###### PART 4: compute the weighted mean with the mosaic function #####

    ##Make this a function later
    #list_weights_m <- list(list_linear_weights_m,list_edge_weights_m,list_sine_weights_m)
    #list_weights_prod_m <- list(list_linear_weights_prod_m,list_edge_weights_prod_m,list_sine_weights_prod_m)
    #list_methods <- c("linear","edge","sine")
    list_mosaiced_files <- vector("list",length=1)

    list_args_weights <- list_weights_m
    list_args_weights_prod <- list_weights_prod_m
    method_str <- method
  
    #list_args_weights <- (mixedsort(list.files(pattern="r_weights_m_.*.tif")))
    list_args_weights <- lapply(1:length(list_args_weights), FUN=function(i,x){raster(x[[i]])},x=list_args_weights)

    #get the list of weights product into raster objects
    #list_args_weights_prod <- list_weights_prod_m
    #list_args_weights_prod <- (mixedsort(list.files(pattern="r_weights_prod_m_.*.tif")))
    list_args_weights_prod <- lapply(1:length(list_args_weights_prod), FUN=function(i,x){raster(x[[i]])},x=list_args_weights_prod)

    list_args_weights_prod$fun <- "sum"
    list_args_weights_prod$na.rm <- TRUE
    
    list_args_weights$fun <- "sum" #we want the sum to compute the weighted mean
    list_args_weights$na.rm <- TRUE

    #Mosaic files
    r_weights_sum <- do.call(mosaic,list_args_weights) #weights sum image mosaiced
    r_prod_sum <- do.call(mosaic,list_args_weights_prod) #weights sum product image mosacied

    r_m_weighted_mean <- r_prod_sum/r_weights_sum #this is the mosaic using weighted mean...
    raster_name <- file.path(out_dir,paste("r_m_",method_str,"_weighted_mean_",out_suffix,".tif",sep=""))
    writeRaster(r_m_weighted_mean, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  
    
  }

  if(mosaic_method=="unweighted"){
    #### Fourth use original images
    #macth file to mosaic extent using the original predictions
    lf_files <- lf_mosaic
    list_param_raster_match <- list(lf_files,rast_ref,file_format,out_suffix,out_dir)
    names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","out_suffix","out_dir_str")

    list_pred_m <- mclapply(1:length(lf_files),FUN=raster_match,list_param=list_param_raster_match,mc.preschedule=FALSE,mc.cores = num_cores)                           

    #list_mosaiced_files <- list.files(pattern="r_m.*._weighted_mean_.*.tif")

    #names(list_mosaiced_files) <- c("edge","linear","sine")
  
    #### NOW unweighted mean mosaic

    #get the original predicted image to raster (matched previously to the mosaic extent)
    list_args_pred_m <- list_pred_m
    #list_args_pred_m <- (mixedsort(list.files(pattern="^gam_CAI.*.m_mosaic_run10_1500x4500_global_analyses_03252015.tif")))
    list_args_pred_m <- lapply(1:length(list_args_pred_m), FUN=function(i,x){raster(x[[i]])},x=list_args_pred_m)

    list_args_pred_m$fun <- "mean"
    list_args_pred_m$na.rm <- TRUE

    #Mosaic files
    r_m_mean <- do.call(mosaic,list_args_pred_m) #this is unweighted mean from the predicted raster

    raster_name <- file.path(out_dir,paste("r_m_mean_",out_suffix,".tif",sep=""))
    writeRaster(r_m_mean, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  #unweighted mean

    #r_m_mean_unweighted <- paste("r_m_mean_",out_suffix,".tif",sep="")
    list_weights <- NULL
    list_weights_prod <- NULL

  }
  
  #Create return object
  mosaic_obj <- list(raster_name,list_weights,list_weights_prod,mosaic_method)
  names(mosaic_obj) <- c("mean_mosaic","r_weights","r_weigths_prod","method")
  save(mosaic_obj,file=file.path(out_dir,paste(mosaic_method,"_","mosaic_obj_",out_suffix,".RData",sep="")))
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

##################### END OF SCRIPT ######################

