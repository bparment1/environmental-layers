####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Assessment of product part 1 functions: mosaic and accuracy ##############################
#This script contains functions and uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Combining tables and figures for individual runs for years and tiles.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 05/24/2016  
#MODIFIED ON: 01/27/2017            
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

##COMMIT: clearning up functions and splitting with part3 product assessment
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
library(lubridate)
library(mosaic)

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

#Functions used in the validation metric script
calc_val_metrics<-function(x,y){
  #This functions calculates accurayc metrics on given two vectors.
  #Arguments: list of fitted models, raster stack of covariates
  #Output: spatial grid data frame of the subset of tiles
  #s_sgdf<-as(r_stack,"SpatialGridDataFrame") #Conversion to spatial grid data frame
  
  residuals<-x-y
  mae<-mean(abs(residuals),na.rm=T)
  rmse<-sqrt(mean((residuals)^2,na.rm=T))
  me<-mean(residuals,na.rm=T)
  r<-cor(x,y,use="complete")
  m50<-median(residuals,na.rm=T)
  metrics_dat<-as.data.frame(cbind(mae,rmse,me,r,m50))
  names(metrics_dat)<-c("mae","rmse","me","r","m50")
  metrics_obj<-list(metrics_dat,as.data.frame(residuals))
  names(metrics_obj)<-c("metrics_dat","residuals")
  return(metrics_obj)
}

combine_measurements_and_predictions_df <- function(i,df_raster,df_time_series, df_points_extracted,data_var,list_selected_ID,r_ts_name,var_name,var_pred,num_cores=1,scaling=NULL,out_dir=".",out_suffix="",plot_fig=T){
  
  # Input arguments:
  #1) i : selected station
  #2) df_raster:
  #3) df_time_series
  #4) df_points_extracted : data extracted from raster layer
  #5) data_var : data with station measurements (tmin,tmax or precip), this can be a list of files
  #6) list_selected_ID : list of selected station
  #7) r_ts_name: names of the raster layers extracted and relevant to combine
  #8) var_name: variable predicted (e.g. dailyTmax)
  #9) var_pred: model predicted (e.g. "mod1)
  #10) sclaing: if NUll no scaing, the extracted value is multiplied by scaling
  #11) num_cores: 1, number of cores used in reading data
  #10) out_dir: output directory
  #11) out_suffix_str: output suffix
  #12) plot_fig : if T, figures are plotted
  # Output
  #
  
  ##### START FUNCTION ############
  
  #######
  ## STEP 1: Select stations by ID 
  
  id_name <- list_selected_ID[i] # e.g. WS037.00,1238099999
  #id_selected <- df_ts_pix[[var_ID]]==id_name
  id_selected <- df_points_extracted[["id"]]== id_name #select specific station from data extracted from raster stack
  #id_selected <- df_ts_pix[["id"]]== id_name #select specific station from data extracted from raster stack
  
  #######
  ## STEP 2: Select stations with data extracted from raster stack 
  ## Now get time series
  data_pixel <- df_points_extracted[id_selected,] #this should be a unique row!!!, e.g. 1x10295
  data_pixel <- as.data.frame(data_pixel)
  
  ##Transpose data to have rows as date and one unique column
  pix_ts <- t(as.data.frame(subset(data_pixel,select=r_ts_name))) #subset to relevant layers, e.g. climate time series
  #pix_ts <- subset(as.data.frame(pix_ts),select=r_ts_name)
  pix_ts <- (as.data.frame(pix_ts)) # transform the matrix transposed into data.frame
  
  names(pix_ts) <- paste(var_pred,"_mosaic",sep="")
  #add scaling option
  if(!is.null(scaling)){
    pix_ts[paste(var_pred,"_mosaic",sep="")] <- pix_ts[paste(var_pred,"_mosaic",sep="")]*scaling
  }
  #!is.null(scaling)
  
  ########
  ## STEP 3: Process the measurements data (with tmax/tmin/precip)
  ## 
  
  if(class(data_var)!="data.frame"){
    #
    browser()
    lf_data_var_subset <- mclapply(data_var,
                           FUN=extract_from_df,
                           col_selected="id",
                           val_selected=id_name,
                           mc.preschedule=FALSE,
                           mc.cores = num_cores)   
    #took less than 2 minutes to extract one station data for 31 years (one year by file)
    df_tmp <- do.call(rbind,lf_data_var_subset)
  }else{
    df_tmp <- subset(data_var,data_var$id==id_name)
  }
  #there are several measurements per day for some stations !!!

  var_pix <- as.data.frame(df_tmp) #select only dates and var_name!!!
  var_pix$date_str <- as.character(var_pix$date)
  #match from 20011231 to 2001-12-31 to date format
  var_pix$date_str <- as.character(as.Date(var_pix$date_str,"%Y%m%d")) #format back to the relevant date format for files
  
  dates_val <- df_raster$date
  pix_ts$date <- dates_val 
  #pix_ts <- merge(df_raster,pix_ts,by="date")
  
  pix_ts$lf <- df_raster$lf
  #Combine data with list of range and missing: this is fast
  pix_ts <- merge(df_time_series,pix_ts,by="date",all=T)
  pix_ts$id_val <- id_name
      
  #var_pred_tmp <- paste0(var_pred,"_mosaic")
  #check for duplicates in extracted values (this can happen if there is a test layer or repetition
  if(nrow(pix_ts)!=length(unique(pix_ts$date))){

    md <- melt(pix_ts, id=(c("date")),measure.vars=c(var_pred_tmp,"missing")) #c("x","y","dailyTmax","mod1","res_mod1"))
    #formula_str <- "id + date ~ x + y + dailyTmax + mod1 + res_mod1"
    pix_ts <- cast(md, date ~ variable, fun.aggregate = mean, na.rm = TRUE)
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
  
  ## Combined extracted values (pix_ts) with data observation (var_pix)
  df_pix_ts <- merge(pix_ts,var_pix,by="date",all=T)
  #Create time series object from extract pixel time series

  df_pix_ts$month_str <- format(as.Date(df_pix_ts$date), "%b") ## Month, char, abbreviated
  df_pix_ts$year <- format(as.Date(df_pix_ts$date), "%Y") ## Year with century
  df_pix_ts$day <- as.numeric(format(as.Date(df_pix_ts$date), "%d")) ## numeric month
  
  #compute residuals from mosaics
  var_pred_tmp <- paste0(var_pred,"_mosaic")
  df_pix_ts[[paste0("res_",var_pred_tmp)]] <- df_pix_ts[[var_pred_tmp]] - df_pix_ts[[var_name]]
  #browser()
  #id_name <- list_selected_ID[i]
  df_pix_ts_filename <- file.path(out_dir,paste0("df_pix_ts_id_",id_name,"_",var_name,"_",out_suffix_str,".txt"))
  write.table(df_pix_ts,df_pix_ts_filename,sep=",")
  
  #Compute accuracy metrics
  metrics_obj <- calc_val_metrics(x=df_pix_ts[[var_pred_tmp]],y=df_pix_ts[[var_name]])
  
  nb_zero <- sum( df_pix_ts[[var_pred_tmp]]==0,na.rm=T) #relevant for precip
  nb_NA_pred <- sum(is.na( df_pix_ts[[var_pred_tmp]])) #for ID 394 DMR it is 361 missing values for 2012!!
  nb_NA_obs <- sum(is.na( df_pix_ts[[var_name]])) #for ID 394 DMR it is 361 missing values for 2012!!
  n <- length(df_pix_ts[[var_name]])
  
  number_data_df <- data.frame(id=id_name,n_zero=nb_zero,n_NA_pred=nb_NA_pred,n_NA_obs=nb_NA_obs,n=n)
  metric_dat <- metrics_obj$metrics_dat
  
  metric_stat_df <- cbind(number_data_df,metric_dat)
  station_summary_obj <- list(metric_stat_df,df_pix_ts,df_pix_ts_filename )
  names(station_summary_obj) <- c("metric_stat_df","df_pix_ts","df_pix_ts_filename")
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


### combine training and testing by year
#also count testing and training used by day 
combine_and_aggregate_df_data_fun <- function(i,list_data_df_training,list_data_df_testing,selected_var=NULL,fun_selected_var="mean",list_out_suffix=NULL,out_dir="."){
  #i,list_data_df_training,list_data_df_testing,selected_var=NULL,fun_selected_var="mean",out_suffix="",out_dir="."
  #
  
  if(is.null(list_out_suffix)){
    out_suffix_str <- ""
  }else{
    out_suffix_str <- list_out_suffix[i]
  }
  
  data_s_df  <- read.table(list_data_df_training[i],header=T,stringsAsFactors=F,sep=",")
  data_v_df <- read.table(list_data_df_testing[i],header=T,stringsAsFactors=F,sep=",")
  
  data_s_df$training <- 1
  data_v_df$testing <- 1
  
  ## use merge function
  #df_combined_data <- do.call(rbind,list(data_df1,data_df2)) #reading only the years related to the the dates e.g. 1999
  data_stations <- rbind.fill(data_v_df, data_s_df) #should work?
  ### Write out combined training and testing data
  filename_data_stations_combined_v_s <- file.path(out_dir,paste0("data_stations_combined_v_s_",out_suffix_str,".txt"))
  write.table(data_stations,file=filename_data_stations_combined_v_s,sep=",")
    
  ##Add tile id here...and check if data stations was training or testing.

  #01/23 21:58  
  md <- melt(data_stations, id=(c("id", "date")),measure.vars=c("x","y","dailyTmax","mod1","res_mod1"))
  data_stations_var_pred <- cast(md, id + date ~ variable, fun.aggregate = mean, na.rm = TRUE)
  #01/23 22:08
  
  #01/23 22:10
  md2 <- melt(data_stations, id=(c("id", "date")),measure.vars=c("testing","training"))
  data_stations_var_pred_tmp2 <- cast(md2, id + date ~ variable, fun.aggregate = sum,na.rm = TRUE)
  #01/23 22:14  
  
  ### Now combine both
  #write.table(data_stations_var_pred,
  #            file=file.path(out_dir,paste0("data_stations_var_pred_tmp_",out_suffix,".txt",
  #                                                                 sep=",")))
  
  #test <- merge(data_stations_var_pred2,data_stations_var_pred,by=c(("id")))  
  #test <- data_stations_var_pred  

  data_stations_var_pred$testing <- data_stations_var_pred_tmp2$testing
  data_stations_var_pred$training <- data_stations_var_pred_tmp2$training
  #dim(data_stations_var_pred)
    
  #An inner join of df1 and df2:
  #Return only the rows in which the left table have matching keys in the right table.
   
  if(!is.null(selected_var)){
    
    #8:42
    #test2<- merge(test,data_stations[,c("id",selected_var)],by=c("id"),all=F)
    md3 <- melt(data_stations, id=(c("id", "date")),measure.vars=selected_var)
    data_stations_var_pred_tmp3 <- cast(md3, id + date ~ variable, fun.aggregate = min,na.rm = TRUE)
    #10:14
    #Error in Summary.factor(integer(0), na.rm = FALSE) : 
    #‘min’ not meaningful for factors
    data_stations_var_pred <- cbind(data_stations_var_pred, data_stations_var_pred_tmp3)
    dim(data_stations_var_pred)
  }
             
  data_stations_var_pred$date_str <- data_stations_var_pred$date
  data_stations_var_pred$date <- as.Date(strptime(data_stations_var_pred$date_str,"%Y%m%d"))
  #dim(data_stations_var_pred)
  #[1] 462885     10
  
  #> length(unique(data_stations_var_pred$id))
  #[1] 1464
  #length(unique(data_s_df$id))
  #[1] 1464
  #length(unique(data_v_df$id))
  #[1] 1458
  #interesect(unique(data_v_df$id),unique(data_s_df$id))

  #data_s_df
  
  #data_stations_var_pred$year <- as.Date(strptime(data_stations_var_pred$date_str,"%Y%m%d"))

  ### Write out data combined
  filename_data_stations_var_pred <- file.path(out_dir,paste0("data_stations_var_pred_",out_suffix_str,".txt"))
  write.table(data_stations_var_pred,file=filename_data_stations_var_pred,sep=",")
  ### Write out combined training and testing data
  #filename_data_stations_training_testing <- file.path(out_dir,paste0("data_stations_training_testing_",out_suffix,".txt"))
  #write.table(data_stations_training_testing,file=filename_data_stations_training_testing,sep=",")
  
  #Prepare return object
  combine_data_obj <- list(filename_data_stations_var_pred,filename_data_stations_combined_v_s)
  names(combine_data_obj) <- c("data_stations_var_pred","data_stations_combined_v_s")
  return(combine_data_obj)
}


############################ END OF SCRIPT ##################################