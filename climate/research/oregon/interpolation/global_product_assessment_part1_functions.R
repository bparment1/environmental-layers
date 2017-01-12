####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Assessment of product part 1 functions: mosaic and accuracy ##############################
#This script contains functions and uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Combining tables and figures for individual runs for years and tiles.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 05/24/2016  
#MODIFIED ON: 01/12/2016            
#Version: 1
#PROJECT: Environmental Layers project     
#COMMENTS: fixing bugs in extraction from raster time series and missing day functions 
#TODO:
#1) Add plot broken down by year and region 
#2) Modify code for overall assessment accross all regions and year
#3) Clean up

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rmd user:aguzman4:rwx /nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015

##COMMIT: added function for aggregate_by_id_and_coord
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
library(lubridate)
library(mosaic)

###### Function used in the script #######
  
#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

pre_process_raster_mosaic_fun <- function(i,list_param){
  
  
  ## Extract parameters
  
  lf <- list_param$lf
  python_bin <- list_param$python_bin
  infile_mask <- list_param$infile_mask
  scaling <- list_param$scaling
  mask_pred <- list_param$mask_pred
  matching <- list_param$matching
  NA_flag_val <- list_param$NA_flag_val
  out_suffix <- list_param$out_suffix
  out_dir <- list_param$out_dir
  
  raster_name_in <- lf[i]
  
  #Step 1: match extent and resolution
  if(matching==TRUE){
    lf_files <- c(raster_name_in) #match to mask
    rast_ref <- infile_mask
    ##Maching resolution is probably only necessary for the r mosaic function
    #Modify later to take into account option R or python...
    list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix,out_dir)
    names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")
    r_pred_matched <- raster_match(1,list_param_raster_match)
    raster_name_in <- c(r_pred_matched)
  }

  #Step 2: mask
  if(mask_pred==TRUE){
    r_mask <- raster(infile_mask)
    extension_str <- extension(raster_name_in)
    raster_name_tmp <- gsub(extension_str,"",basename(raster_name_in))
    raster_name <- file.path(out_dir,paste(raster_name_tmp,"_masked.tif",sep = ""))
    r_pred <- mask(raster(raster_name_in),r_mask,filename = raster_name,overwrite = TRUE)
  }
  
  NAvalue(r_pred) <- NA_flag_val
  #r_pred <- setMinMax(r_pred)

  #Step 3: remove scaling factor
  raster_name_in <- filename(r_pred)
  extension_str <- extension(raster_name_in)
  raster_name_tmp <- gsub(extension_str,"",basename(filename(r_pred)))
  raster_name_out <- file.path(out_dir,paste(raster_name_tmp,"_rescaled.tif",sep = ""))
  #r_pred <- overlay(r_pred, fun=function(x){x*1/scaling},filename=raster_name,overwrite=T)

  #raster_name_in <- "comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19990102_reg4_1999_m_gam_CAI_dailyTmax_19990102_reg4_1999_m__meeting_NASA_reg4_04292016_masked.tif"
  python_cmd <- file.path(python_bin,"gdal_calc.py")
  cmd_str3 <- paste(python_cmd, 
                     paste("-A ", raster_name_in,sep=""),
                     paste("--outfile=",raster_name_out,sep=""),
                     #paste("--type=","Int32",sep=""),
                     "--co='COMPRESS=LZW'",
                     paste("--NoDataValue=",NA_flag_val,sep=""),
                     paste("--calc='(A)*",scaling,"'",sep=""),
                     "--overwrite",sep=" ") #division by zero is problematic...
  #system(cmd_str3)
  system(cmd_str3)
  #NAvalue(r_pred) <- NA_flag_val
  #r_pred <- 
  r_pred <- setMinMax(raster(raster_name_out))
  
  return(raster_name_out)
}

plot_raster_mosaic <- function(i,list_param){
  #Function to plot mosaic for poster
  #
  l_dates <- list_param$l_dates
  r_mosaiced_scaled <- list_param$r_mosaiced_scaled
  NA_flag_val <- list_param$NA_flag_val
  out_dir <- list_param$out_dir
  out_suffix <- list_param$out_suffix
  region_name <- list_param$region_name
  variable_name <- list_param$variable_name
  zlim_val <- list_param$zlim_val

  #for (i in 1:length(nlayers(r_mosaic_scaled))){
  
  date_proc <- l_dates[i]
  r_pred <- subset(r_mosaiced_scaled,i)
  NAvalue(r_pred)<- NA_flag_val 
 
  raster_name <- filename(r_pred)
  extension_str <- extension(raster_name)
  raster_name_tmp <- gsub(extension_str,"",basename(raster_name))
  
  date_proc <- l_dates[i]
  date_val <- as.Date(strptime(date_proc,"%Y%m%d"))
  #month_name <- month.name(date_val)
  month_str <- format(date_val, "%b") ## Month, char, abbreviated
  year_str <- format(date_val, "%Y") ## Year with century
  day_str <- as.numeric(format(date_val, "%d")) ## numeric month
  date_str <- paste(month_str," ",day_str,", ",year_str,sep="")
  
  res_pix <- 1200
  #res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1
  
  
  if(is.null(zlim_val)){
    zlim_val_str <- paste(c(minValue(r_pred),maxValue(r_pred)),sep="_",collapse="_")
    #png_filename <-  file.path(out_dir,paste("Figure4_clim_mosaics_day_","_",date_proc,"_",region_name,"_zlim_",zlim_val_str,"_",out_suffix,".png",sep =""))
    raster_name_tmp
    png_filename <-  file.path(out_dir,paste("Figure_",raster_name_tmp,"_zlim_",zlim_val_str,"_",out_suffix,".png",sep =""))
    
    title_str <-  paste("Predicted ",variable_name, " on ",date_str , " ", sep = "")
  
    png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)

    plot(r_pred,main =title_str,cex.main =1.5,col=matlab.like(255),
       legend.shrink=0.8,legend.width=0.8)
       #axis.args = list(cex.axis = 1.6), #control size of legend z
       #legend.args=list(text='dNBR', side=4, line=2.5, cex=2.2))
       #legend.args=list(text='dNBR', side=4, line=2.49, cex=1.6))
    dev.off()
  }else{
    zlim_val_str <- paste(zlim_val,sep="_",collapse="_")
    #png_filename <-  file.path(out_dir,paste("Figure_mosaics_day_","_",date_proc,"_",region_name,"_",zlim_val_str,"_",out_suffix,".png",sep =""))
    png_filename <-  file.path(out_dir,paste("Figure_",raster_name_tmp,"_zlim_",zlim_val_str,"_",out_suffix,".png",sep =""))

    title_str <-  paste("Predicted ",variable_name, " on ",date_str , " ", sep = "")
    png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
    
    plot(r_pred,main =title_str,cex.main =1.5,col=matlab.like(255),zlim=zlim_val,
       legend.shrink=0.8,legend.width=0.8)
       #axis.args = list(cex.axis = 1.6), #control size of legend z
       #legend.args=list(text='dNBR', side=4, line=2.5, cex=2.2))
       #legend.args=list(text='dNBR', side=4, line=2.49, cex=1.6))
    dev.off()
  }
  return(png_filename)
}

extract_date <- function(i,x,item_no=NULL){
  y <- unlist(strsplit(x[[i]],"_"))
  if(is.null(item_no)){
    date_str <- y[length(y)-2] #count from end
  }else{
    date_str <- y[item_no]
  }
  return(date_str)
}

gClip <- function(shp, bb, keep.attribs=TRUE,outDir=NULL,outSuffix=NULL,write=FALSE){
  #Purpose: clipping SpatialPolygonsDataFrame using another SpatialPolygonsDataFrame 
  #shp: input shapefile that we would like to clip
  #bb: input shapefile used for clipping, can be and extent raster object, matrix of coordinates 
  #keep.attribs: join attributes to spatial feature
  #outDir: output directory
  #outSuffix: output suffix attached to the name of new file
  
  #Authors: Benoit Parmentier, Modified code originating at: 
  #http://robinlovelace.net/r/2014/07/29/clipping-with-r.html
  
  #Comments: 
  #- Note that some attribute should be updated: e.g. areas, length etc. of original polygons
  #- Add bbox option for spdf
  #- Send error if different carthographic projection used
  
  ### BEGIN ####
  
  #First check inputs used from clipping, output dir and suffix
  
  if(class(bb) == "matrix"){
    b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons") #if matrix of coordinates
  }
  if(class(bb)=="SpatialPolygonsDataFrame"){
    b_poly <- bb #if polygon object, keep the same
  } 
  
  if(class(bb)=="SpatialPointsDataFrame"){
    b_poly <- as(extent(bb), "SpatialPolygons") #make a Spatial Polygon from raster extent
  }
  
  if(class(bb)=="exent"){
    b_poly <- as(extent(bb), "SpatialPolygons") #make a Spatial Polygon from raster extent
  }
  rm(bb) #remove from memory in case the polygon file is a giant dataset
  
  #If no output dir present, use the current dir
  if(is.null(outDir)){
    outDir=getwd()
  }
  #if no output suffix provided, use empty string for now
  if(is.null(outSuffix)){
    outSuffix=""
  }
  
  #Second, clip using rgeos library
  new.shape <- gIntersection(shp, b_poly, byid = T)
  
  #Third, join the atrribute back to the newly created object
  if(keep.attribs){
    #create a data.frame based on the spatial polygon object
    new.attribs <- data.frame(do.call(rbind,strsplit(row.names(new.shape)," ")),stringsAsFactors = FALSE)
    #test <-over(shp,bb)
    
    #new.attrib.data <- shp[new.attribs$X1,]@data #original data slot? #not a good way to extract...
    new.attrib.data <- as.data.frame(shp[new.attribs$X1,])
    row.names(new.shape) <- row.names(new.attrib.data)
    new.shape <-SpatialPolygonsDataFrame(new.shape, new.attrib.data) #Associate Polygons and attributes

  }
  
  if(write==TRUE){

    #Writeout shapefile (default format for now)
    infile_new.shape <- paste("clipped_spdf",outSuffix,".shp",sep="")
    writeOGR(new.shape,dsn= outDir,layer= sub(".shp","",infile_new.shape), 
           driver="ESRI Shapefile",overwrite_layer="TRUE")
  }

  
  return(new.shape)
}

plot_stations_val_by_date <- function(i,list_param){
  #
  #
  #function to plot residuals by date
  #
  
  #####
  
  l_dates <- list_param$l_dates
  model_name <- list_param$model_name
  station_type_name <- list_param$station_type_name
  var_selected <- list_param$var_selected
  #list_df_points_dates <- 
  df_combined_data_points <- list_param$df_combined_data_points
  countries_shp <- list_param$countries_shp
  proj_str <- list_param$proj_str
  out_suffix <- list_param$out_suffix
  out_dir <- list_param$out_dir
  
  ### STARTFUNCTION ####
    
  date_processed <- l_dates[i]
  
  df_points <- subset(df_combined_data_points,date==date_processed)
  #df_points <- list_df_points_dates[[i]]
  #plot(df_points)
  freq_tile_id <- as.data.frame(table(df_points$tile_id))
  freq_tile_id$date <- date_processed
  names(freq_tile_id) <- c("tile_id","freq","date")
    
  #plot(df_points,add=T)
  coordinates(df_points) <- cbind(df_points$x,df_points$y)
  #proj4string(df_points) <- CRS_locs_WGS84
  proj4string(df_points) <- proj_str
  # # No spatial duplicates
  df_points_day <- remove.duplicates(df_points) #remove duplicates...
  #plot(df_points_day)
  #dim(df_points_day)
  #dim(df_points)
  #plot(df_points)
  
  ##layer to be clipped
  if(class(countries_shp)!="SpatialPolygonsDataFrame"){
    reg_layer <- readOGR(dsn=dirname(countries_shp),sub(".shp","",basename(countries_shp)))
  }else{
    reg_layer <- countries_shp
  }

  #extent_df_points_day <- extent(df_points_day)
  
  reg_layer_clipped <- gClip(shp=reg_layer, bb=df_points_day , keep.attribs=TRUE,outDir=NULL,outSuffix=NULL)
  
  #data_stations_var_pred <- aggregate(id ~ date, data = data_stations, min)
  #data_stations_var_pred <- aggregate(id ~ x + y + date + dailyTmax + mod1 + res_mod1 , data = data_stations, min)
  
  #change the formula later on to use different y_var_name (tmin and precip)
  data_stations_var_pred <- aggregate(id ~ x + y + date + dailyTmax + res_mod1 + tile_id + reg ,data = df_points, mean ) #+ mod1 + res_mod1 , data = data_stations, min)
  dim(data_stations_var_pred)
  #> dim(data_stations_var_pred)
  #[1] 11171     5

  data_stations_var_pred$date_str <- data_stations_var_pred$date
  data_stations_var_pred$date <- as.Date(strptime(data_stations_var_pred$date_str,"%Y%m%d"))
  coordinates(data_stations_var_pred) <- cbind(data_stations_var_pred$x,data_stations_var_pred$y)
  #data_stations_var_pred$constant <- c(1000,2000)
  #plot(data_stations_var_pred)
  #plot(reg_layer)
  #res_pix <- 1200
  res_pix <- 800
  col_mfrow <- 1
  row_mfrow <- 1
  
  png_filename <- paste("Figure_ac_metrics_map_stations_location_",station_type_name,"_",model_name,"_",y_var_name,"_",date_processed,
                       "_",out_suffix,".png",sep="")
  png(filename=png_filename,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  #plot(data_stations_var_pred
  #p_shp <- spplot(reg_layer_clipped,"ISO3" ,col.regions=NA, col="black") #ok problem solved!!
  #title("(a) Mean for 1 January")
  #p <- bubble(data_stations_var_pred,"constant",main=paste0("Average Residuals by validation stations ",
  #                                                      date_processed,
  #                                                      "for ",var_selected))
  #p <- spplot(data_stations_var_pred,"constant",col.regions=NA, col="black",
  #            main=paste0("Average Residuals by validation stations ",pch=3,cex=10,
  #            date_processed, "for ",var_selected))

  #p1 <- p+p_shp
  #try(print(p1)) #error raised if number of missing values below a threshold does not exist
  
  title_str <- paste0("Stations ",station_type_name," on ",date_processed, " for ",y_var_name," ", var_selected)
  plot(reg_layer_clipped,main=title_str)
  plot(data_stations_var_pred,add=T,cex=0.5,col="blue")
  legend("topleft",legend=paste("n= ", nrow(data_stations_var_pred)),bty = "n")
  
  dev.off()
  
  res_pix <- 800
  col_mfrow <- 1
  row_mfrow <- 1
  png_filename <- paste("Figure_ac_metrics_map_stations",station_type_name,"averaged","_fixed_intervals_",model_name,y_var_name,date_processed,
                       out_suffix,".png",sep="_")
  png(filename=png_filename,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
  #model_name[j]
  title_str <- paste0("Average Residuals by ",station_type_name," stations ",date_processed, " for ",y_var_name," ", var_selected)
  #p_shp <- layer(sp.polygons(reg_layer, lwd=1, col='black'))
  p_shp <- spplot(reg_layer,"ISO3" ,col.regions=NA, col="black") #ok problem solved!!
  #title("(a) Mean for 1 January")
  class_intervals <- c(-20,-10,-5,-4,-3,-2,-1,0,1,2,3,4,5,10,20)
  p <- bubble(data_stations_var_pred,zcol="res_mod1",key.entries = class_intervals , 
              fill=F, #plot circle with filling
              #col= matlab.like(length(class_intervals)),
              main=title_str)
  p1 <- p + p_shp
  try(print(p1)) #error raised if number of missing values below a threshold does not exist
  
  dev.off()

  #### Add additional plot with quantiles and min max?...
  
  
  #library(classInt)
  #breaks.qt <- classIntervals(palo_alto$PrCpInc, n = 6, style = "quantile", intervalClosure = "right")
  #spplot(palo_alto, "PrCpInc", col = "transparent", col.regions = my.palette, 
  #  at = breaks.qt$brks)

  #### histogram of values
  res_pix <- 800
  col_mfrow <- 1
  row_mfrow <- 1
  png_filename <- paste("Figure_ac_metrics_histograms_stations",station_type_name,"averaged",model_name,y_var_name,date_processed,
                       out_suffix,".png",sep="_")
  png(filename=png_filename,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  title_str <- paste0("Stations ",station_type_name," on ",date_processed, " for ",y_var_name," ", var_selected)

  h<- histogram(data_stations_var_pred$res_mod1,breaks=seq(-50,50,5),
                main=title_str)
  
  print(h)
  dev.off()
  
  ##Make data.frame with dates for later use!!
  #from libary mosaic
  df_basic_stat <- fav_stats(data_stations_var_pred$res_mod1)
  df_basic_stat$date <- date_processed
  #df_basic_stat$reg <- reg
  #quantile(data_stations_var_pred$res_mod1,c(1,5,10,90,95,99))
  df_quantile_val <- quantile(data_stations_var_pred$res_mod1,c(0.01,0.05,0.10,0.45,0.50,0.55,0.90,0.95,0.99))
  #names(df_quantile_val)
  #as.list(df_quantile_val)
  #df_test <- data.frame(names(df_quantile_val))[numeric(0), ]


  #quantile(c(1,5,10,90,95,99),data_stations_var_pred$res_mod1,)
  #rmse(data_stations_var_pred$res_mod1)
  
  plot_obj <- list(data_stations_var_pred,df_basic_stat,df_quantile_val,freq_tile_id)
  names(plot_obj) <- c("data_stations_var_pred","df_basic_stat","df_quantile_val","freq_tile_id")
   
  return(plot_obj)
}

extract_from_df <- function(x,col_selected,val_selected){
  df_tmp <- read.table(x,stringsAsFactors=F,sep=",")
  #data_subset <- subset(data_stations,col_selected==val_selected)
  data_subset <- subset(df_tmp,df_tmp$id%in%val_selected)
  return(data_subset)
}

sub_sampling_by_dist <- function(target_range_nb=c(10000,10000),dist_val=0.0,max_dist=NULL,step_val,data_in){
  #Function to select stations data that are outside a specific spatial range from each other
  #Parameters:
  #max_dist: maximum spatial distance at which to stop the pruning
  #min_dist: minimum distance to start pruning the data
  #step_val: spatial distance increment
  #Note that we are assuming that the first columns contains ID with name col of "id".
  #Note that the selection is based on unique id of original SPDF so that replicates screened.
  
  data_in$id <- as.character(data_in$id)
  data <- data_in
  
  #Now only take unique id in the shapefile!!!
  #This step is necessary to avoid the large calculation of matrix distance with replicates
  #unique(data$id)
  data <- aggregate(id ~ x + y , data=data,min)
  coordinates(data) <- cbind(data$x,data$y)
  proj4string(data) <- proj4string(data_in)
  
  target_min_nb <- target_range_nb[1]
  #target_min_nb <- target_range_day_nb[1]
  
  #station_nb <- nrow(data_in)
  station_nb <- nrow(data)
  if(is.null(max_dist)){
    while(station_nb > target_min_nb){
      data <- remove.duplicates(data, zero = dist_val) #spatially sub sample...
      dist_val <- dist_val + step_val
      station_nb <- nrow(data)
    }
    #setdiff(as.character(data$id),as.character(data_in$id))
    #ind.selected <-match(as.character(data$id),as.character(data_in$id)) #index of stations row selected
    #ind.removed  <- setdiff(1:nrow(data_in), ind.selected) #index of stations rows removed 
    id_selected <- as.character(data$id)
    id_removed <- setdiff(unique(as.character(data_in$id)),as.character(data$id))

  }
  if(!is.null(max_dist)){
    
    while(station_nb > target_min_nb & dist_val < max_dist){ 
      data <- remove.duplicates(data, zero = dist_val) #spatially sub sample...
      #id_rm <- zerodist(data, zero = dist_val, unique.ID = FALSE)
      #data_rm <- data[id_rm,]
      dist_val <- dist_val + step_val
      station_nb <- nrow(data)
    }
    #ind.selected <- match(as.character(data$id),as.character(data_in$id))
    id_selected <- as.character(data$id)
    id_removed <- setdiff(unique(as.character(data_in$id)),as.character(data$id))
  #  ind.removed  <- setdiff(1:nrow(data_in), ind.selected)
  }
    
  #data_rm <- data_in[ind.removed,]
  data_rm <- subset(data_in, id %in% id_removed)
  data_tmp <- data #store the reduced dataset with only id, for debugging purpose
  
  #data <- subset(data_in, id %in% data$id) #select based on id
  data <-subset(data_in, id %in% id_selected) #select based on id
  
  #data <- data_in[ind.selected,]
  obj_sub_sampling <- list(data,dist_val,data_rm) #data.frame selected, minimum distance, data.frame stations removed
  names(obj_sub_sampling) <- c("data","dist","data_rm")
  return(obj_sub_sampling)
}

query_for_station_lat_long <- function(point_val,df_points_spdf,step_x=0.25,step_y=0.25){
  #make this function better using gbuffer later!!!, works for now 
  #Improve: circular query + random within it
  y_val <- point_val[2]
  x_val <- point_val[1]
  
  y_val_tmp <- y_val + step_y
  if(x_val>0){
    x_val_tmp <- x_val - step_x
  }else{
    x_val_tmp <- x_val + step_x
  }


  test_df <- subset(df_points_spdf,(y_val < df_points_spdf$y & df_points_spdf$y < y_val_tmp))
  test_df2 <- subset(test_df,(x_val < test_df$x & test_df$x < x_val_tmp))
  #plot(test_df2)
 if(nrow(test_df2)>0){
   df_selected <- test_df2[1,]
 }else{
   df_selected <- NULL
 }
 
 return(df_selected)
}



finding_missing_dates <- function(date_start,date_end,list_dates){
  #this assumes daily time steps!!
  #can be improved later on
  
  #date_start <- "19840101"
  #date_end <- "19991231"
  date1 <- as.Date(strptime(date_start,"%Y%m%d"))
  date2 <- as.Date(strptime(date_end,"%Y%m%d"))
  dates_range <- seq.Date(date1, date2, by="1 day") #sequence of dates

  missing_dates <- setdiff(as.character(dates_range),as.character(list_dates))
  #df_dates_missing <- data.frame(date=missing_dates)
  #which(df_dates$date%in%missing_dates)
  #df_dates_missing$missing <- 1
  
  df_dates <- data.frame(date=as.character(dates_range),missing = 0) 

  df_dates$missing[df_dates$date %in% missing_dates] <- 1
  #a$flag[a$id %in% temp] <- 1

  missing_dates_obj <- list(missing_dates,df_dates)
  names(missing_dates_obj) <- c("missing_dates","df_dates")
  return(missing_dates_obj)
}



extract_from_time_series_raster_stack <- function(df_points,date_start,date_end,lf_raster,item_no=13,num_cores=4,pattern_str=NULL,in_dir=NULL,out_dir=".",out_suffix=""){
  #
  #This function extract value given from a raster stack time series given a spatial data.frame and a list of files
  #
  #INPUTS
  #1) df_points
  #2) date_start,num_cores=4,pattern_str=NULL,in_dir=NULL,out_dir=".",out_suffix=
  #3) date_end
  #3) lf_raster
  #4) item_no=13
  #5) num_cores=4,
  #6) pattern_str=NULL
  #7) in_dir=NULL,
  #8) out_dir="."
  #9) out_suffix=""
  #OUTPUTS
  #
  #
  
  #### Start script ####
  
  if(is.null(lf_raster)){
    #pattern_str <- ".*.tif"
    pattern_str <-"*.tif"
    lf_raster <- list.files(path=in_dir_mosaic,pattern=pattern_str,recursive=F,full.names=T)
    r_stack <- stack(lf_raster,quick=T) #this is very fast now with the quick option!
    #save(r_mosaic,file="r_mosaic.RData")
    
  }else{
    r_stack <- stack(lf_raster,quick=T) #this is very fast now with the quick option!
  }

  #df_points$files <- lf_mosaic_list
  #Use the global output??

  ##23.09 (on 05/22)
  #df_points_day_extracted <- extract(r_mosaic,data_stations,df=T)
  #df_points_day_extracted_fname <- paste0("df_points_day_extracted_fname",".txt") 
  #write.table(df_points_day_extracted,file=df_points_day_extracted_fname) #5.1 Giga
  #4.51 (on 05/23)
  #df_points_day <- data_stations_var_pred_data_s

  #15.17 (on 09/08)
  ##10.41 (on 05/22)
  #took about 7h for 5262 layers, maybe can be sped up later
  #took about 16h for 10289 layers (reg1) to extract 4800 points
  df_points_extracted <- extract(r_stack,df_points,df=T,sp=T) #attach back to the original data...

  #17.19 (on 05/23)
  #22.27 (on 09/08)
  #df_points_extracted_fname <- paste0("df_points_day_extracted_fname2",".txt")
  #17.27 (on 05/23)
  
  df_points_extracted_fname <- file.path(out_dir,paste0("df_points_extracted_",out_suffix,".txt"))
  write.table(df_points_extracted,file= df_points_extracted_fname,sep=",",row.names = F) 
  #17.19 (on 05/23)
  #browser()
  
  #### Now check for missing dates
  
  #debug(extract_date)
  #test <- extract_date(6431,lf_mosaic_list,12) #extract item number 12 from the name of files to get the data
  #list_dates_produced <- unlist(mclapply(1:length(lf_raster),FUN=extract_date,x=lf_raster,item_no=13,mc.preschedule=FALSE,mc.cores = num_cores))                         
  #list_dates_produced <-  mclapply(1:2,FUN=extract_date,x=lf_mosaic_list,item_no=13,mc.preschedule=FALSE,mc.cores = 2)      
  #takes about 1 hour for 10289 files from stack
  list_dates_produced <- unlist(mclapply(1:length(lf_raster),FUN=extract_date,x=lf_raster,item_no=item_no,
                                         mc.preschedule=FALSE,mc.cores = num_cores))                         

  list_dates_produced_date_val <- as.Date(strptime(list_dates_produced,"%Y%m%d"))
  month_str <- format(list_dates_produced_date_val, "%b") ## Month, char, abbreviated
  year_str <- format(list_dates_produced_date_val, "%Y") ## Year with century
  day_str <- as.numeric(format(list_dates_produced_date_val, "%d")) ## numeric month

  df_raster <- data.frame(lf=basename(lf_raster),
                          date=list_dates_produced_date_val,
                          month_str=month_str,
                          year=year_str,
                          day=day_str,
                          dir=dirname(lf_raster))

  df_raster_fname <- file.path(out_dir,paste0("df_raster_",out_suffix,".txt"))
  write.table(df_raster,file= df_raster_fname,sep=",",row.names = F) 

  missing_dates_obj <- finding_missing_dates(date_start,date_end,list_dates_produced_date_val)
  
  df_time_series <- missing_dates_obj$df_dates
  df_time_series$date <- as.character(df_time_series$date)  
  df_raster$date <- as.character(df_raster$date)
  
  df_time_series <- merge(df_time_series,df_raster,by="date",all=T) #outer join to keep missing dates
  
  df_time_series_fname <- file.path(out_dir,paste0("df_time_series_",out_suffix,".txt")) #add the name of var later (tmax)
  write.table(df_time_series,file= df_time_series_fname,sep=",",row.names = F) 
  
  extract_obj <- list(df_points_extracted_fname,df_raster_fname,df_time_series_fname)
  names(extract_obj) <- c("df_points_extracted_fname","df_raster_fname","df_time_series_fname")
  
  #extract_obj_fname <- file.path(out_path,paste("raster_extract_obj_",interpolation_method,"_", y_var_name,out_prefix_str[i],".RData",sep=""))
  extract_obj_fname <- file.path(out_dir,paste("raster_extract_obj_",out_suffix,".RData",sep=""))
  save( extract_obj,file= extract_obj_fname)

  return(extract_obj)
}


combine_measurements_and_predictions_df <- function(i,df_raster,df_time_series,df_ts_pix,data_var,list_selected_ID,r_ts_name,var_name,var_pred,out_dir=".",out_suffix="",plot_fig=T){
  
  # Input arguments:
  # i : selected station
  # df_ts_pix_data : data extracted from raster layer
  # data_var : data with station measurements (tmin,tmax or precip)
  # list_selected_ID : list of selected station
  # plot_fig : if T, figures are plotted
  # Output
  #
  
  ##### START FUNCTION ############
  
  #get the relevant station
  id_name <- list_selected_ID[i] # e.g. WS037.00,1238099999
  #id_selected <- df_ts_pix[[var_ID]]==id_name
  id_selected <- df_ts_pix[["id"]]== id_name
  
  ### Not get the data from the time series
  data_pixel <- df_ts_pix[id_selected,] #this should be a unique row!!!
  #data_pixel <- data_pixel[1,]
  data_pixel <- as.data.frame(data_pixel)
  
  ##Transpose data to have rows as date and one unique column
  pix_ts <- t(as.data.frame(subset(data_pixel,select=r_ts_name))) #can subset to range later
  #pix_ts <- subset(as.data.frame(pix_ts),select=r_ts_name)
  pix_ts <- (as.data.frame(pix_ts))
  names(pix_ts) <- paste(var_pred,"_mosaic",sep="")
  #add scaling option
  #!is.null(scaling)
  ## Process the measurements data (with tmax/tmin/precip)
  
  #there are several measurements per day for some stations !!!
  #id_name <- data_pixel[[var_ID]]
  
  #df_tmp  <-data_var[data_var$LOCATION_ID==id_name,]
  df_tmp <- subset(data_var,data_var$id==id_name)
  #if(da)
  #aggregate(df_tmp
  #if(nrow(df_tmp)>1){
  #  
  #  formula_str <- paste(var_name," ~ ","TRIP_START_DATE_f",sep="")
  #  #var_pix <- aggregate(COL_SCORE ~ TRIP_START_DATE_f, data = df_tmp, mean) #aggregate by date
  #  var_pix <- try(aggregate(as.formula(formula_str), data = df_tmp, FUN=mean)) #aggregate by date
  #  #length(unique(test$TRIP_START_DATE_f))
  #  #var_pix_ts <- t(as.data.frame(subset(data_pixel,select=var_name)))
  #  #pix <- t(data_pixel[1,24:388])#can subset to range later
  #}else{
  #  var_pix <- as.data.frame(df_tmp) #select only dates and var_name!!!
  #}
  #var_pix <- subset(as.data.frame(data_id_selected,c(var_name,"TRIP_START_DATE_f")])) #,select=var_name)
  var_pix <- as.data.frame(df_tmp) #select only dates and var_name!!!
  var_pix$date_str <- as.character(var_pix$date)
  #match from 20011231 to 2001-12-31 to date format
  var_pix$date <- as.character(as.Date(var_pix$date_str,"%Y%m%d")) #format back to the relevant date format for files
  
  #dates_val <- df_time_series$date
  dates_val <- df_raster$date
  pix_ts$date <- dates_val 
  #pix_ts <- merge(df_raster,pix_ts,by="date")
  
  pix_ts$lf <- df_raster$lf
  #pix_ts$
  pix_ts <- merge(df_time_series,pix_ts,by="date",all=T)
  
  #check for duplicates in extracted values (this can happen if there is a test layer or repetition
  if(nrow(pix_ts)!=length(unique(pix_ts$date))){
    var_pred_tmp <- paste0(var_pred,"_mosaic")
    md <- melt(pix_ts, id=(c("date")),measure.vars=c(var_pred_tmp, "missing")) #c("x","y","dailyTmax","mod1","res_mod1"))
    #formula_str <- "id + date ~ x + y + dailyTmax + mod1 + res_mod1"
    pix_ts <- cast(md, date ~ variable, fun.aggregate = mean, 
    na.rm = TRUE)

  }
  
  #if(nrow(var_pix)!=length(unique(var_pix$date))){
  #
  #  md <- melt(var_pix, id=(c("date")),measure.vars=c(var_pred, "missing")) #c("x","y","dailyTmax","mod1","res_mod1"))
  #  #formula_str <- "id + date ~ x + y + dailyTmax + mod1 + res_mod1"
  #  test <- cast(md, date ~ variable, fun.aggregate = mean, 
  #  na.rm = TRUE)
  #
  #  
  #}
  df_pix_ts <- merge(pix_ts,var_pix,by="date",all=T)
  #Create time series object from extract pixel time series

  df_pix_ts$month_str <- format(as.Date(df_pix_ts$date), "%b") ## Month, char, abbreviated
  df_pix_ts$year <- format(as.Date(df_pix_ts$date), "%Y") ## Year with century
  df_pix_ts$day <- as.numeric(format(as.Date(df_pix_ts$date), "%d")) ## numeric month
  
  #compute residuals from mosaics
  df_pix_ts[[paste0("res_",var_pred_tmp)]] <- df_pix_ts[[var_pred_tmp]] - df_pix_ts[[y_var_name]]
  
  id_name <- list_selected_ID[i]
  df_pix_ts_filename <- file.path(out_dir,paste0("df_pix_ts_",id_name,y_var_name,out_suffix,".txt"))
  write.table(df_pix_ts,df_pix_ts_filename,sep=",")


  nb_zero <- sum( df_pix_ts[[var_pred_tmp]]==0) #relevant for precip
  #nb_NA <- sum(is.na(df2$COL_SCORE))
  nb_NA <- sum(is.na( df_pix_ts[[var_pred_tmp]])) #for ID 394 DMR it is 361 missing values for 2012!!
  ##Add quantile, and range info later on...
  
  ## Cumulated precip and lag?
  #Keep number of  0 for every year for rainfall
  #summarize by month
  #Kepp number of NA for scores... 
  #Summarize by season...
  ## Threshold?
  station_summary_obj <- list(nb_zero,nb_NA, df_pix_ts,df_pix_ts_filename )
  names(station_summary_obj) <- c("nb_zero_precip","nb_NA_var","df_pix_ts","df_pix_ts_filename")
  return(station_summary_obj)
}

aggregate_by_id_and_coord <- function(i,list_df_data,list_out_suffix,out_dir){
  #This functions aggrate iput data.frame based on the ID of the station and coordinates x,y
  #
  #
  
  df_points_data <- list_df_data[[i]]
  out_suffix_str <- list_out_suffix[i]
  ##Begin
  if(class(df_points_data)!="data.frame"){
    df_points_data <- read.table(df_points_data,stringsAsFactors=F,sep=",")
  }
  #df_points_data <- (list_df_v_stations[[1]])
  #test3 <- aggregate(id  ~x + y,data=test,FUN=mean)
  df_station_data <- aggregate(id  ~ x + y,data=df_points_data,FUN=mean)
  out_filename <- file.path(out_dir,paste0("df_station_data_",out_suffix_str,".txt"))
  write.table(df_station_data,out_filename,sep=",")
  return(df_station_data)
}

############################ END OF SCRIPT ##################################