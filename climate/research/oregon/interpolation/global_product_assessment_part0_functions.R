####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part 0 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This script checks the number of predictions by tiles and years.
#with the goal of predicting potential gaps or missing predictions in fugure mosaics by region.
#The general logic is to check the number of overlap by shapefile polyon tiles
#along with the predicitons for every day of the year (*.tif)
#Summary tables and data are also produced in the script.
#
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/31/2016  
#MODIFIED ON: 11/04/2016            
#Version: 1
#PROJECT: Environmental Layers project     
#COMMENTS: removing unused functions and clean up for part0 global prodduct assessment part0 
#TODO:#PROJECT: Environmental Layers project     
#COMMENTS:
#TODO:
#1) Add plot broken down by year and region 
#2) Modify code for overall assessment accross all regions and year
#3) Clean up

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rmd user:aguzman4:rwx /nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015

##COMMIT: adding necessary function for assessment from other script

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
#library(mosaic)

###### Function used in the script #######
  
#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
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
check_missing <- function(lf, pattern_str=NULL,in_dir=".",date_start="1984101",date_end="20141231",item_no=13,out_suffix="",num_cores=1,out_dir="."){
  #Function to check for missing files such as mosaics or predictions for tiles etc.
  #The function assumes the name of the files contain "_".
  #INPUTS:
  #1) lf
  #2) pattern_str
  #3) in_dir
  #4) date_start
  #5) date_end
  #6) item_no
  #7) num_cores
  #8) out_suffix
  #9) out_dir
  #OUTPUTS
  #
  #
  
  ##### Start script #####
  
  out_dir <- in_dir
  
  list_dates_produced <- unlist(mclapply(1:length(lf),
                                         FUN = extract_date,
                                         x = lf,
                                         item_no = item_no,
                                         mc.preschedule = FALSE,
                                         mc.cores = num_cores))
  
  list_dates_produced_date_val <- as.Date(strptime(list_dates_produced, "%Y%m%d"))
  month_str <- format(list_dates_produced_date_val, "%b") ## Month, char, abbreviated
  year_str <- format(list_dates_produced_date_val, "%Y") ## Year with century
  day_str <- as.numeric(format(list_dates_produced_date_val, "%d")) ## numeric month
  df_files <- data.frame(lf = basename(lf),
                         date = list_dates_produced_date_val,
                         month_str = month_str,
                         year = year_str,
                         day = day_str,
                         dir = dirname(lf))
  
  df_files_fname <- file.path(out_dir, paste0("df_files_", out_suffix, ".txt"))
  write.table(df_files,file = df_files_fname,sep = ",",row.names = F)
  
  #undebug(finding_missing_dates )
  missing_dates_obj <- finding_missing_dates (date_start,date_end,list_dates_produced_date_val)
  
  df_time_series <- missing_dates_obj$df_dates
  df_time_series$date <- as.character(df_time_series$date)  
  df_files$date <- as.character(df_files$date)
  
  df_time_series <- merge(df_time_series,df_files,by="date",all=T) #outer join to keep missing dates
  
  df_time_series$month_str <- format(as.Date(df_time_series$date), "%b") ## Month, char, abbreviated
  df_time_series$year_str <- format(as.Date(df_time_series$date), "%Y") ## Year with century
  df_time_series$day <- as.numeric(format(as.Date(df_time_series$date), "%d")) ## numeric month
  
  df_time_series_fname <- file.path(out_dir,paste0("df_time_series_",out_suffix,".txt")) #add the name of var later (tmax)
  write.table(df_time_series,file= df_time_series_fname,sep=",",row.names = F) 
  
  df_time_series_obj <- list(df_time_series_fname,df_time_series_fname,df_time_series)
  names(df_time_series_obj) <- c("df_time_series_fname","df_time_series_fname","df_time_series")
  
  ## report in text file missing by year and list of dates missing in separate textfile!!
  return(df_time_series_obj)
}

centroids_shp_fun <- function(i,list_shp_reg_files){
  
  #
  shp_filename <- list_shp_reg_files[[i]]
  layer_name <- sub(".shp","",basename(shp_filename))
  path_to_shp <- dirname(shp_filename)
  shp1 <- try(readOGR(path_to_shp, layer_name)) #use try to resolve error below
  #shp_61.0_-160.0
  #Geographical CRS given to non-conformant data: -186.331747678
    
  #shp1<-readOGR(dirname(list_shp_reg_files[[i]]),sub(".shp","",basename(list_shp_reg_files[[i]])))
  if (!inherits(shp1,"try-error")) {
    pt <- gCentroid(shp1)
    #centroids_pts[[i]] <- pt
  }else{
    pt <- shp1
    #centroids_pts[[i]] <- pt
  }
    
  #shps_tiles[[i]] <- shp1
  #centroids_pts[[i]] <- centroids
    
  shp_obj <- list(shp1,pt)
  names(shp_obj) <- c("spdf","centroid")
  return(shp_obj)
}

rasterize_tile_day <- function(i,list_spdf,df_missing,list_r_ref,col_name,date_val){
  #
  #
  tile_spdf <- list_spdf[[i]]
  tile_coord <- names(list_spdf)[i]
  r_ref <- list_r_ref[[i]]
    
  df_tmp <- subset(df_missing,date==date_val,select=tile_coord)
  #for each row (date)
  val <- df_tmp[[tile_coord]]
  if(val==1){
    val<-0 #missing then not predicted
  }else{
    val<-1
  }
  
  tile_spdf$predicted <- val
  tile_spdf$tile_coord <- tile_coord
  tile_spdf$overlap <- 1
    
  #r <- rasterize(tile_spdf,r_ref,"predicted")
  #r <- rasterize(tile_spdf,r_ref,col_name)
  #r <- raster(r_ref,crs=projection(r_ref)) #new layer without values
  r <- raster(r_ref) #new layer without values
  
  if(col_name=="overlap"){
    set1f <- function(x){rep(1, x)}
   	r <- init(r, fun=set1f, overwrite=TRUE)
  }
  if(col_name=="predicted"){
    set1f <- function(x){rep(val, x)}
   	r <- init(r, fun=set1f, overwrite=TRUE)
  }
    
  return(r)
}

############################ END OF SCRIPT ##################################