##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part2 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Accuracy methods are added in the the function scripts to evaluate results.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 03/23/2014  
#MODIFIED ON: 01/03/2016            
#Version: 5
#PROJECT: Environmental Layers project     
#COMMENTS: analyses for run 10 global analyses,all regions 1500x4500km with additional tiles to increase overlap 
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

#function_analyses_paper1 <-"contribution_of_covariates_paper_interpolation_functions_07182014.R" #first interp paper
#function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_08132014.R"

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

#turn term from list into data.frame
#name_col<-function(i,x){
#x_mat<-x[[i]]
#x_df<-as.data.frame(x_mat)
#x_df$mod_name<-rep(names(x)[i],nrow(x_df))
#x_df$term_name <-row.names(x_df)
#return(x_df)
#}
#Function to rasterize a table with coordinates and variables...,maybe add option for ref image??
rasterize_df_fun <- function(data_tb,coord_names,proj_str,out_suffix,out_dir=".",file_format=".rst",NA_flag_val=-9999,tolerance_val= 0.000120005){
  data_spdf <- data_tb
  coordinates(data_spdf) <- cbind(data_spdf[[coord_names[1]]],data_spdf[[coord_names[2]]])
  proj4string(data_spdf) <- proj_str

  data_pix <- try(as(data_spdf,"SpatialPixelsDataFrame"))
  #tolerance_val <- 0.000120005 
  #tolerance_val <- 0.000856898
  if(inherits(data_pix,"try-error")){
      data_pix <- SpatialPixelsDataFrame(data_spdf, data=data_spdf@data, tolerance=tolerance_val) 
  }
  
  #test <- as(data_spdf,"SpatialPixelsDataFrame")

  # set up an 'empty' raster, here via an extent object derived from your data
  #e <- extent(s100[,1:2])
  #e <- e + 1000 # add this as all y's are the same

  #r <- raster(e, ncol=10, nrow=2)
  # or r <- raster(xmn=, xmx=,  ...

  data_grid <- as(data_pix,"SpatialGridDataFrame") #making it a regural grid
  r_ref <- raster(data_grid) #this is the ref image
  rast_list <- vector("list",length=ncol(data_tb))
  
  for(i in 1:(ncol(data_tb))){
    field_name <- names(data_tb)[i]
    var <- as.numeric(data_spdf[[field_name]])
    data_spdf$var  <- var
    #r <-rasterize(data_spdf,r_ref,field_name)
    r <-rasterize(data_spdf,r_ref,"var",NAflag=NA_flag_val,fun=mean) #prolem with NA in NDVI!!

    data_name<-paste("r_",field_name,sep="") #can add more later...
    #raster_name<-paste(data_name,out_names[j],".tif", sep="")
    raster_name<-paste(data_name,out_suffix,file_format, sep="")
  
    writeRaster(r, NAflag=NA_flag_val,filename=file.path(out_dir,raster_name),overwrite=TRUE)
    #Writing the data in a raster file format...
    rast_list[i] <-file.path(out_dir,raster_name)
  }
  return(unlist(rast_list))
}

plot_raster_tb_diagnostic <- function(reg_layer,tb_dat,df_tile_processed,date_selected,mod_selected,var_selected,out_suffix){
  
  test <- subset(tb_dat,pred_mod==mod_selected & date==date_selected,select=c("tile_id",var_selected))

  test_data_tb <- merge(df_tile_processed,test,by="tile_id",all=T) #keep all
  test_r <- subset(test_data_tb,select=c("lat","lon","tile_id",var_selected))
  out_suffix_str <- paste(var_selected,mod_selected,date_selected,out_suffix,sep="_")
  coord_names <- c("lon","lat")
  l_rast <- rasterize_df_fun(test_r,coord_names,proj_str,out_suffix_str,out_dir=".",file_format=".tif",NA_flag_val=-9999,tolerance_val=0.000120005)
  #mod_kr_stack <- stack(mod_kr_l_rast)
  d_tb_rast <- stack(l_rast)
  names(d_tb_rast) <- names(test_r)
  #plot(d_tb_rast)
  r <- subset(d_tb_rast,"rmse")
  names(r) <- paste(mod_selected,var_selected,date_selected,sep="_")
  #plot info: with labels

  res_pix <- 1200
  col_mfrow <- 1
  row_mfrow <- 1

  png(filename=paste("Figure9_",names(r),"_map_processed_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  #plot(reg_layer)
  #p1 <- spplot(reg_layer,"ISO",colorkey=FALSE) #Use ISO instead of NAME_1 to see no color?
  title_str <- paste(names(r),"for ", region_name,sep="")

  p0 <- levelplot(r,col.regions=matlab.like(25),margin=F,main=title_str)
  p_shp <- layer(sp.polygons(reg_layer, lwd=1, col='black'))

  p <- p0 + p_shp
  print(p)

  dev.off()

}

create_raster_from_tb_diagnostic <- function(i,list_param){
  #create a raster image using tile centroids and given fields  from tb diagnostic data
  tb_dat <- list_param$tb_dat
  df_tile_processed <- list_param$df_tile_processed
  date_selected <- list_param$date_selected[i]
  mod_selected <- list_param$mod_selected
  var_selected <- list_param$var_selected
  out_suffix <- list_param$out_suffix
  
  test <- subset(tb_dat,pred_mod==mod_selected & date==date_selected,select=c("tile_id",var_selected))

  test_data_tb <- merge(df_tile_processed,test,by="tile_id",all=T) #keep all
  test_r <- subset(test_data_tb,select=c("lat","lon","tile_id",var_selected))
  out_suffix_str <- paste(var_selected,mod_selected,date_selected,out_suffix,sep="_")
  coord_names <- c("lon","lat")
  l_rast <- rasterize_df_fun(test_r,coord_names,proj_str,out_suffix_str,out_dir=".",file_format,NA_flag_val,tolerance_val=0.000120005)
  #mod_kr_stack <- stack(mod_kr_l_rast)
  #d_tb_rast <- stack(l_rast)
  #r <- subset(d_tb_rast,var_selected)
  #names(d_tb_rast) <- names(test_r)
  return(l_rast[4])
}

assign_FID_spatial_polygons_df <-function(list_spdf,ID_str=NULL){
  list_spdf_tmp <- vector("list",length(list_spdf))
  if(is.null(ID_str)){
    nf <- 0 #number of features
    #for(i in 1:length(spdf)){
    #    shp1 <- list_spdf[[i]]
    #    f <- nrow(shp1)
    #    nf <- nf + f
    #}
    #This assumes that the list has one feature per item list
    nf <- length(list_spdf)
    ID_str <- as.character(1:nf)
  }
  for(i in 1:length(list_spdf)){
    #test=spRbind(shps_tiles[[1]],shps_tiles[[2]])
    shp1 <- list_spdf[[i]]
    shp1$FID <- ID_str
    shp1<- spChFIDs(shp1, as.character(shp1$FID)) #assign ID
    list_spdf_tmp[[i]]  <-shp1
  }
  return(list_spdf_tmp)
}

combine_spatial_polygons_df_fun <- function(list_spdf_tmp,ID_str=NULL){
  if(is.null(ID_str)){
    #call function
    list_spdf_tmp <- assign_FID_spatial_polygons_df
  }
  combined_spdf <- list_spdf_tmp[[1]]
  for(i in 2:length(list_spdf_tmp)){
    combined_spdf <- rbind(combined_spdf,list_spdf_tmp[[i]])
    #sapply(slot(shps_tiles[[2]], "polygons"), function(x) slot(x, "ID"))
    #rownames(as(alaska.tract, "data.frame"))
  }
  return(combined_spdf)
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

  
##################### END OF SCRIPT ######################
