####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part 0 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This script checks the number of predictions by tiles and years.
#with the goal of predicting potential gaps or missing predictions in fugure mosaics by region.
#The general logic is to check the number of overlap by shapefile polyon tiles
#along with the predictions for every day of the year (*.tif)
#Summary tables and data are also produced in the script.
#
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/31/2016  
#MODIFIED ON: 04/28/2017            
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

##COMMIT: modifying output figures and files 

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
  df_files <- data.frame(lf =lf,
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

rasterize_tile_day <- function(i,list_spdf,df_missing,list_r_ref,col_name,date_val,out_dir=".",out_suffix=""){
  #
  #This function creates a raster image from tiles and missing information data.frame.
  #
  #INPUTS:
  #1) i : counter to consider tile being processed
  #2) list_spdf: list of spatial polygon data.frame to convert in raster
  #3) df_missing: data.frame with information used in rasterization
  #4) list_r_ref: reference raster to use
  #5) col_name: name of column containing value used in the rasteriation, the column is stored in the df_missing data.frame
  #6) date_val:
  #7) out_dir:
  #8) out_suffix
  #OUTPUTS
  #1) raster_name: output raster_name 
  #
  
  #### Start script ####
  
  tile_spdf <- list_spdf[[i]]
  tile_coord <- names(list_spdf)[i]
  r_ref <- list_r_ref[[i]]
    
  tile_id <- paste0("tile_",i)
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
  tile_spdf$tile_id <- tile_id
    
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
    
  #out_dir <- "." #can set this up later
  out_suffix_str <- paste0(col_name,out_suffix) # can set this parameter later
  raster_name <- file.path(out_dir,paste("r_",tile_id,"_",tile_coord,"_",out_suffix_str,".tif",sep=""))
  #raster_name <- 
  writeRaster(r, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  #unweighted mean

  return(raster_name)
}

generate_raster_number_of_prediction_by_day <- function(i,list_param){
  #
  ##This function generates raster of missing pixels and number of predictions for days with missing tiles for a given region.
  #
  #INPUTS
  #1) list_tiles_predicted_masked 
  #2) df_missing_tiles_day     
  #3) r_overlap_m 
  #4) num_cores
  #5) region_name 
  #6) NA_flag_val 
  #7)  scaling <- list_param$scaling
  #8) python_bin <- list_param$python_bin
  #9) data_type <- list_param$data_type
  #10) plotting figures: if True plot png for missing day predictions
  #11) out_suffix 
  #12) out_dir 
  #OUTPUTS
    
    
  ###### Start script #####
    
  #### read in paramters    
  list_tiles_predicted_masked <- list_param$list_tiles_predicted_masked
  df_missing_tiles_day <- list_param$df_missing_tiles_day    
  r_overlap_m <- list_param$r_overlap_m
  item_no <- list_param$item_no
  num_cores <- list_param$num_cores # 6 #PARAM 14
  region_name <- list_param$region_name #<- "world" #PARAM 15
  NA_flag_val <-list_param$NA_flag_val
  scaling <- list_param$scaling
  python_bin <- list_param$python_bin
  data_type <- list_param$data_type
  plotting_figures <- list_param$plotting_figures
  out_suffix  <- list_param$out_suffix
  out_dir  <- list_param$out_dir

    
  ### To add or explore later...could have differences between predictions and rmse
  layers_option <- c("var_pred") #arg 17 ,param 17, options are:#res_training, res_testing,ac_training, ac_testing, var_pred
  NA_value <- NA_flag_val 
  #metric_name <- "rmse" #to be added to the code later...
  #data_type <- "Int16" #, param 19, use int32 for output layers mosaiced

  if(data_type=="Int16"){
    data_type_str <- "INT2S"
  }

  ##### Select relevant day and create stack of missing tiles
    
  missing_tiles <- df_missing_tiles_day[i,]
  date_str <- missing_tiles$date
  selected_col <- names(list_tiles_predicted_masked)
  missing_tiles_subset <- subset(missing_tiles,select=selected_col)
  selected_missing <- missing_tiles_subset==1

  list_missing_tiles_raster <- list_tiles_predicted_masked[selected_missing]
  r_tiles_s <- stack(list_missing_tiles_raster)
    
  ##### Sum missing tiles in the stack and generate number of predictions by pixels
  ## This stores files in the temp dir
  raster_name_data_sum <- file.path(out_dir,paste("r_data_sum","_",region_name,"_masked_",date_str,"_tmp",file_format,sep=""))
  r_data_sum <- stackApply(r_tiles_s, 1:nlayers(r_tiles_s), fun = sum,filename=raster_name_data_sum,overwrite=TRUE)
     
  ### then substract missing tiles...
  raster_name_number_prediction <- file.path(out_dir,paste("r_day_number_of_prediction_sum_day_mosaiced","_",region_name,"_masked_",date_str,file_format,sep=""))

  #r_day_predicted <- r_overlap_m - r_data_sum
  #r_day_predicted <- overlay(r_overlap_m, r_data_sum, fun=function(x,y) x - y, filename=raster_name_number_prediction, overwrite=TRUE)

  #raster_name_day_predicted <- file.path(out_dir,paste("r_day_predicted","_",region_name,"_masked_",date_str,file_format,sep=""))
  ## do this in gdalcalc or overlay function to go faster?
  python_cmd <- file.path(python_bin,"gdal_calc.py")
  cmd_str1 <- paste(python_cmd, 
                     paste("-A ", filename(r_overlap_m),sep=""),
                     paste("-B ", raster_name_data_sum,sep=""),
                     paste("--outfile=", raster_name_number_prediction,sep=""),
                     paste("--type=",data_type,sep=""),
                     "--co='COMPRESS=LZW'",
                     paste("--NoDataValue=",NA_flag_val,sep=""),
                     paste("--calc='(A-B)*",scaling,"'",sep=""),
                     "--overwrite",sep=" ") #division by zero is problematic...
  system(cmd_str1)
  
  #raster_name_number_prediction <- file.path(out_dir,paste("r_day_number_of_prediction_sum_day_mosaiced","_",region_name,"_masked_",date_str,file_format,sep=""))
  #writeRaster(r_day_predicted, NAflag=NA_flag_val,filename=raster_name_number_prediction,overwrite=TRUE)  
  r_day_predicted <- raster(raster_name_number_prediction)
  ### Uses too much resources, change here...
  #r_table <- ratify(r_day_predicted) # build the Raster Attibute table
  #rat <- levels(r_table)[[1]]#get the values of the unique cell frot the attribute table
  tb_freq <- as.data.frame(freq(r_day_predicted))
  
  #rat$legend <- paste0("tile_",1:26)
  #tb_freq <- as.data.frame(freq(r_table))
  #rat$legend <- tb_freq$value
  #levels(r_table) <- rat

  if(plotting_figures==TRUE){
    
    res_pix <- 800
    #res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
  
    png_filename_number_of_predictions <-  file.path(out_dir,paste("Figure_number_of_predictions_by_pixel_",date_str,"_",region_name,"_",out_suffix,".png",sep =""))
    title_str <-  paste("Number of predicted pixels for ",variable_name," on ",date_str, sep = "")
  
    png(filename=png_filename_number_of_predictions,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
    #my_col=c('blue','red','green')
    my_col <- rainbow(length(tb_freq$value))
    plot(r_day_predicted,col=my_col,legend=F,box=F,axes=F,main=title_str)
    legend(x='topright', legend =tb_freq$value,fill = my_col,cex=0.8)
    dev.off()
  }else{
    png_filename_number_of_predictions <- NULL
  }

  ### Day missing reclass above
  ## change here
  #r_missing_day <- r_day_predicted == 0
  raster_name_missing <- file.path(out_dir,paste("r_missing_day_mosaiced","_",region_name,"_masked_",date_str,file_format,sep=""))
  
  ## do this in gdalcalc or overlay function to go faster?
  python_cmd <- file.path(python_bin,"gdal_calc.py")
  cmd_str2 <- paste(python_cmd, 
                     paste("-A ", raster_name_number_prediction,sep=""),
                     paste("--outfile=",raster_name_missing,sep=""),
                     paste("--type=",data_type,sep=""),
                     "--co='COMPRESS=LZW'",
                     paste("--NoDataValue=",NA_flag_val,sep=""),
                     paste("--calc='(A<1)*",scaling,"'",sep=""),
                     "--overwrite",sep=" ") #division by zero is problematic...
  system(cmd_str2)
  
  r_missing_day <- raster(raster_name_missing)  
  
  if(plotting_figures==TRUE){
    
    res_pix <- 800
    #res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
  
    png_filename_missing_predictions <-  file.path(out_dir,paste("Figure_missing_predictions_by_pixel_",date_str,"_",region_name,"_",out_suffix,".png",sep =""))
    title_str <-  paste("Number of predicted pixels for ",variable_name," on ",date_str, sep = "")
  
    png(filename=png_filename_missing_predictions,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
    #my_col=c('blue','red','green')
    my_col <- c("black","red")
    plot(r_missing_day,col=my_col,legend=F,box=F,axes=F,main=title_str)
    legend(x='topright', legend =c("prediced","missing"),fill = my_col,cex=0.8)
  
    dev.off()
    
  }else{
    png_filename_missing_predictions <- NULL
  }

  ### writeout data
  #extension_str <- extension(lf_files)
  #raster_name_tmp <- gsub(extension_str,"",basename(lf_files))
  #out_suffix_str <- paste0(region_name,"_",out_suffix)
  #raster_name_missing <- file.path(out_dir,paste("r_missing_day_mosaiced","_",region_name,"_masked_",date_str,file_format,sep=""))
  #writeRaster(r_missing_day, NAflag=NA_flag_val,filename=raster_name_missing,overwrite=TRUE)  
    

  ### generate return object
  obj_number_day_predicted <- list(raster_name_number_prediction,raster_name_missing,tb_freq,
                                   png_filename_number_of_predictions,png_filename_missing_predictions)
  names(obj_number_day_predicted) <- c("raster_name_number_prediction","raster_name_missing","tb_freq",
                                       "png_filename_number_of_predictions","png_missing_predictions")
    
  return(obj_number_day_predicted)
}

predictions_tiles_missing_fun <- function(list_param){
  #Add documentation

  ##############################
  #### Parameters and constants  
  

  in_dir1 <- list_param$in_dir1 
  region_name <- list_param$region_name #e.g. c("reg23","reg4") #run only for one region
  y_var_name <- list_param$y_var_name # e.g. dailyTmax" #PARAM3
  interpolation_method <- list_param$interpolation_method #c("gam_CAI") #PARAM4
  out_suffix <- list_param$out_suffix #output suffix e.g."run_global_analyses_pred_12282015" #PARAM5
  out_dir <- list_param$out_dir #<- "/nobackupp8/bparmen1/" #PARAM6
  create_out_dir_param <-list_param$create_out_dir_param #if TRUE output dir created #PARAM7
  proj_str <- list_param$proj_str # CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84, #PARAM8
  list_year_predicted <- list_param$list_year_predicted # 1984:2004
  file_format <- list_param$file_format #<- ".tif" #format for mosaiced files #PARAM10
  NA_flag_val <- list_param$NA_flag_val #<- -9999  #No data value, #PARAM11
  num_cores <- list_param$num_cores #<- 6 #number of cores used #PARAM13
  plotting_figures <- list_param$plotting_figures #if true run generate png for missing date #PARAm 14
  
  ##for plotting assessment function
  
  item_no <- list_param$item_no  #PARAM15
  day_to_mosaic_range <- list_param$day_to_mosaic_range #PARAM16
  countries_shp <- list_param$countries_shp #PARAM17
  plotting_figures <- list_param$plotting_figures #PARAM18
  #threshold_missing_day <- list_param$threshold_missing_day #PARAM20
  pred_mod_name <- list_param$pred_mod_name #PARAM21
  metric_name <- list_param$metric_name #PARAM22
  
  year_predicted <- list_param$year_predicted #selected year #PARAM 23
  data_type <- list_param$data_type #PARAM 24
  scaling <- list_param$scaling #PARAM 25
  tmp_files <- list_param$tmp_files #PARAM 26
  
  ## generate raster
  raster_overlap <- list_param$raster_overlap #PARAM 27
  raster_pred <-list_param$raster_pred #PARAM 28
  
  ########################## START SCRIPT #########################################
  
  #browser()
  #system("ls /nobackup/bparmen1"
  #out_dir <- in_dir #use directory of out_dir
  if(create_out_dir_param==TRUE){
    out_dir <- create_dir_fun(out_dir,out_suffix)
    setwd(out_dir)
  }else{
    setwd(out_dir) #use previoulsy defined directory
  }
  
  setwd(out_dir)
  
  if(is.null(scaling)){
    scaling <- 1
  }
  #valid_range <- values_range #if NULL don't screen values!!
  #valid_range <- c(-100,100) #pass this as parameter!! (in the next update)
  if(data_type=="Int16"){
    data_type_str <- "INT2S"
  }

  rasterOptions(tmpdir=out_dir)#trying to resolve long time when use from command line
  
  #list_outfiles <- vector("list", length=35) #collect names of output files, this should be dynamic?
  #list_outfiles_names <- vector("list", length=35) #collect names of output files
  #browser(()
  
  in_dir1_reg <- file.path(in_dir1,region_name)
  list_outfiles <- vector("list", length=14) #collect names of output files
  in_dir_list <- list.dirs(path=in_dir1_reg,recursive=FALSE) #get the list regions processed for this run

  #in_dir_list_all  <- unlist(lapply(in_dir_list,function(x){list.dirs(path=x,recursive=F)}))
  in_dir_list_all <- in_dir_list
  in_dir_subset <- in_dir_list_all[grep("subset",basename(in_dir_list_all),invert=FALSE)] #select directory with shapefiles...
  in_dir_shp <- file.path(in_dir_subset,"shapefiles")
  
  #select only directories used for predictions
  #nested structure, we need to go to higher level to obtain the tiles...
  in_dir_reg <- in_dir_list[grep(".*._.*.",basename(in_dir_list),invert=FALSE)] #select directory with shapefiles...
  #in_dir_reg <- in_dir_list[grep("july_tiffs",basename(in_dir_reg),invert=TRUE)] #select directory with shapefiles...
  in_dir_list <- in_dir_reg
  
  in_dir_list <- in_dir_list[grep("bak",basename(basename(in_dir_list)),invert=TRUE)] #the first one is the in_dir1
  #list of shapefiles used to define tiles
  in_dir_shp_list <- list.files(in_dir_shp,".shp",full.names=T)
  
  ## load problematic tiles or additional runs
  #modify later...

  ##raster_prediction object : contains testing and training stations with RMSE and model object
  in_dir_list_tmp <- file.path(in_dir_list,year_predicted)
  list_raster_obj_files <- mclapply(in_dir_list_tmp,
                                    FUN=function(x){list.files(path=x,pattern="^raster_prediction_obj.*.RData",full.names=T)},
                                    mc.preschedule=FALSE,mc.cores = num_cores)
  
  #list_raster_obj_files <- try(lapply(in_dir_list_tmp,FUN=function(x){list.files(path=x,pattern="^raster_prediction_obj.*.RData",full.names=T)}))
  #Add stop message here...if no raster object in any tiles then break from the function
  
  list_names_tile_coord <- lapply(list_raster_obj_files,FUN=function(x){basename(dirname(x))})
  list_names_tile_id <- paste("tile",1:length(list_raster_obj_files),sep="_")
  names(list_raster_obj_files)<- list_names_tile_id
  
  #pred_mod_name <- "mod1"
  list_lf_raster_tif_tiles <- mclapply(in_dir_list_tmp,
                                    FUN=function(x){list.files(path=x,pattern=paste0("gam_CAI_",y_var_name,"_predicted_",pred_mod_name,".*.tif"),full.names=T)},
                                    mc.preschedule=FALSE,mc.cores = num_cores)
  list_names_tile_coord <- lapply(list_lf_raster_tif_tiles,FUN=function(x){basename(dirname(dirname(x)))})
  list_names_tile_id <- paste("tile",1:length(list_lf_raster_tif_tiles),sep="_")
  names(list_lf_raster_tif_tiles)<- list_names_tile_id
  
  #one level up
  #lf_covar_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar_obj.*.RData",full.names=T)})
  #lf_covar_tif <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar.*.tif",full.names=T)})
  
  #lf_sub_sampling_obj_files <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern=paste("^sub_sampling_obj_",interpolation_method,".*.RData",sep=""),full.names=T)})
  #lf_sub_sampling_obj_daily_files <- lapply(in_dir_list_tmp,FUN=function(x){list.files(path=x,pattern="^sub_sampling_obj_daily.*.RData",full.names=T)})
  year_processed <- year_predicted
  if(is.null(day_to_mosaic_range)){
  #  start_date <- #first date
     start_date <- paste0(year_processed,"0101") #change this later!!
     end_date <-   paste0(year_processed,"1231") #change this later!!
     day_to_mosaic <- seq(as.Date(strptime(start_date,"%Y%m%d")), as.Date(strptime(end_date,"%Y%m%d")), 'day')
     day_to_mosaic <- format(day_to_mosaic,"%Y%m%d") #format back to the relevant date format for files
  }else{
    start_date <- day_to_mosaic_range[1]
    end_date <- day_to_mosaic_range[2]
    day_to_mosaic <- seq(as.Date(strptime(start_date,"%Y%m%d")), as.Date(strptime(end_date,"%Y%m%d")), 'day')
    day_to_mosaic <- format(day_to_mosaic,"%Y%m%d") #format back to the relevant date format for files
  }
  
  in_dir_tiles_tmp <- in_dir1 #
  #in_dir_tiles_tmp <- in_dir_reg
  
  ### Do this by tile!!!
  
  #gam_CAI_dailyTmax_predicted_mod1_0_1_20001231_30_1-39.7_165.1.tif
  
  #undebug(check_missing)

  test_missing <- try(lapply(1:length(list_lf_raster_tif_tiles),function(i){check_missing(lf=list_lf_raster_tif_tiles[[i]], 
                                                                                          pattern_str=NULL,
                                                                                          in_dir=in_dir,
                                                                                          date_start=start_date,
                                                                                          date_end=end_date,
                                                                                          item_no=item_no, #9 for predicted tiles
                                                                                          out_suffix=out_suffix,
                                                                                          num_cores=num_cores,
                                                                                          out_dir=out_dir)}))

 
  #test_missing <- try(lapply(1:1,function(i){check_missing(lf=list_lf_raster_tif_tiles[[i]], 
  #                                                                    pattern_str=NULL,
  #                                                                    in_dir=out_dir,
  #                                                                    date_start=start_date,
  #                                                                    date_end=end_date,
  #                                                                    item_no=item_no, #9 for predicted tiles
  #                                                                    out_suffix=out_suffix,
  #                                                                    num_cores=num_cores,
  #                                                                    out_dir=".")}))
  
  #browser()
  
  #df_time_series <- test_missing[[1]]$df_time_series
  #head(df_time_series)

  #table(df_time_series$missing)
  #table(df_time_series$year)
  #browser()
  #####
  #Now combine df_time_series in one table
  
  #dim(test_missing[[1]]$df_time_series)
  list_lf <- lapply(1:length(test_missing),FUN=function(i){df_time_series <- as.character(test_missing[[i]]$df_time_series$lf)})
  df_lf_tiles_time_series <- as.data.frame(do.call(cbind,list_lf))
  #http://stackoverflow.com/questions/26220913/replace-na-with-na
  #Use dfr[dfr=="<NA>"]=NA where dfr is your dataframe.
  names(df_lf_tiles_time_series) <- unlist(basename(in_dir_reg))
  filename_df_lf_tiles <- file.path(out_dir,paste0("df_files_by_tiles_predicted_tif_",region_name,"_",pred_mod_name,"_",out_suffix,".txt"))
  write.table(df_lf_tiles_time_series,file=filename_df_lf_tiles)

  ###Now combined missing in one table?
  
  list_missing <- lapply(1:length(test_missing),FUN=function(i){df_time_series <- test_missing[[i]]$df_time_series$missing})
  
  #browser()
  df_missing <- as.data.frame(do.call(cbind,list_missing))
  names(df_missing) <- unlist(basename(in_dir_reg))
  df_missing$tot_missing <- rowSums (df_missing, na.rm = FALSE, dims = 1)
  df_missing$reg <- region_name
  df_missing$date <- day_to_mosaic

  #This contains in rows date and tiles columns
  filename_df_missing <- file.path(out_dir,paste0("df_missing_by_dates_tiles_predicted_tif_",region_name,"_",pred_mod_name,"_",out_suffix,".txt"))
  write.table(df_missing,file=filename_df_missing)
  
  ###### Assessment ####
  #### Generate summary from missing table
  #1) Select dates with missing tiles
  #2) Plot number of days missing in maps over 365 days
  #3) Plot number of days missing in histograms/barplots
  #4) Keep raster of number pix predictions of overlap
  ### do sum across tiles to find number of missing per tiles and map it
  
  df_missing_tiles_sp <- t(df_missing[,1:length(test_missing)]) #transpose to get df with lines centroids of tiles
  df_missing_tiles_sp <- as.data.frame(df_missing_tiles_sp)
  names(df_missing_tiles_sp) <- df_missing$date
  df_missing_tiles_sp$tot_missing <- rowSums(df_missing_tiles_sp) #total missing over a year by tile
  df_missing_tiles_sp$tot_pred <- length(df_missing$date) - df_missing_tiles_sp$tot_missing
  
  list_xy <- strsplit(unlist(basename(in_dir_reg)),"_")
  #list_xy <- lapply(centroids_pts,function(x){coordinates(x)})
  coord_xy <- do.call(rbind,list_xy)
  y_val <- as.numeric(coord_xy[,1]) #lat, long! need to reverse
  x_val <- as.numeric(coord_xy[,2])
  coordinates(df_missing_tiles_sp) <- as.matrix(cbind(x_val,y_val))
  
  #This contains in rows date and tiles columns
  filename_df_missing_tiles_sp <- file.path(out_dir,paste0("df_missing_by_centroids_tiles_and_dates",region_name,"_",pred_mod_name,"_",out_suffix,".txt"))
  write.table(as.data.frame(df_missing_tiles_sp),file=filename_df_missing_tiles_sp )
 
  ### Now generate plots if missing tiles for specific dates in the region
  df_missing_tiles_day <- subset(df_missing,tot_missing > 0)
  path_to_shp <- dirname(countries_shp)
  layer_name <- sub(".shp","",basename(countries_shp))
  reg_layer <- readOGR(path_to_shp, layer_name) #outlines of the region
  r_mask <- raster(infile_mask)
  
  if(nrow(df_missing_tiles_day)>0){

    res_pix <- 800
    #res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
    png_filename_histogram <-  file.path(out_dir,paste("Figure_histogram_",region_missing_tiles_name,"_",out_suffix,".png",sep =""))
    
    png(filename=png_filename_maximum_overlap,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
    hist(df_missing$tot_missing,
         ylab="frequency of missing",
         xlab="tiles",
         main="Number of missing predictions over a year by tile")
    dev.off()
    
    #check for this in every output, if not present then there are no missing tiles over the full year for 
    #the specific region
    write.table(df_missing_tiles_day,file=paste0("df_missing_tiles_day_mosaic_",out_suffix,".txt"))
    
    #do spplot after that on tot sum
    
    png(filename=paste("Figure_total_missing_days_map_centroids_tile_",pred_mod_name,"_",
                       "_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    #spplot for reg_layer not working on NEX again
    #p_shp <- spplot(reg_layer,"ISO3" ,col.regions=NA, col="black") #ok problem solved!!
    p_r <-levelplot(r_mask,colorkey=F) #no key legend
    p <- bubble(df_missing_tiles_sp,"tot_missing",main=paste0("Missing per tile and by ",pred_mod_name,
                                                             " for ",y_var_name))
    #p1 <- p+p_shp
    p_c <- p + p_r + p #set the legend first by using p first
    
    try(print(p_c)) #error raised if number of missing values below a threshold does not exist
    dev.off()
    
  }else{
    #do spplot after that on tot sum
    res_pix <- 800
    #res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
    
    png(filename=paste("Figure_total_predicted_days_predicted_map_centroids_tile_",pred_mod_name,"_",
                       "_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    #bubble(pp, "att",
    #       panel=function(...) {
    #         sp.polygons(SpP, fill="blue")
    #         sp:::panel.bubble(...)
    #       }) 
    
    #p_shp <- spplot(reg_layer,"ISO3" ,sp.layout=list("sp.polygons", reg_layer, fill="white"))
    #p_shp spplot not working on NEX anymore, use raster as background
    p_r <-levelplot(r_mask,colorkey=F) #no key legend
    p <- bubble(df_missing_tiles_sp,"tot_pred",main=paste0("Prediction per tile and by ",pred_mod_name,
                                                          " for ", y_var_name))
    p_c <- p + p_r + p #set the legend first by using p first
    #p1 <- p+p_shp
    try(print(p_c)) #error raised if number of missing values below a threshold does not exist
    dev.off()
  }
  
  ########################
  #### Step 2: Examine tiles layout for the region 
  #browser()
  
  #collect info: read in all shapefiles
  #obj_centroids_shp <- centroids_shp_fun(1,list_shp_reg_files=in_dir_shp_list)
                                         
  obj_centroids_shp <- mclapply(1:length(in_dir_shp_list),
                                FUN=centroids_shp_fun,
                                list_shp_reg_files=in_dir_shp_list,
                                mc.preschedule=FALSE,
                                mc.cores = num_cores)

  centroids_pts <- lapply(obj_centroids_shp, FUN=function(x){x$centroid})
  shps_tiles <-   lapply(obj_centroids_shp, FUN=function(x){x$spdf})

  #remove try-error polygons...we loose three tiles because they extend beyond -180 deg
  tmp <- shps_tiles
  shps_tiles <- remove_errors_list(shps_tiles) #[[!inherits(shps_tiles,"try-error")]]
  #shps_tiles <- tmp
  length(tmp)-length(shps_tiles) #number of tiles with error message
  
  tmp_pts <- centroids_pts 
  centroids_pts <- remove_errors_list(centroids_pts) #[[!inherits(shps_tiles,"try-error")]]
  #centroids_pts <- tmp_pts 
  

  
  ### Plot locations of tiles after?
  #plot(r)
  #plot(shps_tiles[[1]],add=T,border="blue",usePolypath = FALSE) #added usePolypath following error on brige and NEX
  #browser()
  
  ### preparing inputs for raster_overlap production
  names(shps_tiles) <- basename(unlist(in_dir_reg))
  r_ref <- raster(list_lf_raster_tif_tiles[[1]][1])
  list_r_ref <- lapply(1:length(in_dir_reg), function(i){raster(list_lf_raster_tif_tiles[[i]][1])})
  tile_spdf <- shps_tiles[[1]]
  tile_coord <- basename(in_dir_reg[1])
  date_val <- df_missing$date[1]
  
  browser()

  if(raster_overlap==TRUE){
    
    ### use rasterize
    #spdf_tiles <- do.call(bind, shps_tiles) #bind all tiles together in one shapefile
    #Error in (function (classes, fdef, mtable)  : 
    #unable to find an inherited method for function 'bind' for signature '"missing", "missing"'
  
    #undebug(rasterize_tile_day)
    #list_predicted <- rasterize_tile_day(1,
    #        list_spdf=shps_tiles,
    #         df_missing=df_missing,
    #         list_r_ref=list_r_ref,
    #         col_name="overlap",
    #         date_val=df_missing$date[1])
    #list_predicted <- mclapply(1:6,
    #         FUN=rasterize_tile_day,
    #         list_spdf=shps_tiles,
    #         df_missing=df_missing,
    #         list_r_ref=list_r_ref,
    #         col_name = "overlap",
    #         date_val=df_missing$date[1],
    #          mc.preschedule=FALSE,
    #         mc.cores = num_cores)
  
    list_predicted <- mclapply(1:length(shps_tiles),
           FUN=rasterize_tile_day,
           list_spdf=shps_tiles,
           df_missing=df_missing,
           list_r_ref=list_r_ref,
           col_name = "overlap",
           date_val=df_missing$date[1],
           out_dir = out_dir,
           #out_suffix = "",
           out_suffix = "_tmp",
            mc.preschedule=FALSE,
           mc.cores = num_cores)

    ##check that everything is correct:
    #plot(r_mask)
    #plot(raster(list_predicted[[1]]),add=T)
    #plot(spdf_tiles_test,add=T,border="green",usePolypath = FALSE) #added usePolypath following error on brige and NEX

    ### Make a list of file
    out_suffix_str_tmp <- paste0(region_name,"_",out_suffix,"_tmp")
    out_dir_str <- out_dir
    filename_list_predicted <- file.path(out_dir_str,paste("list_to_mosaics_",out_suffix_str_tmp,".txt",sep=""))
    writeLines(unlist(list_predicted),con=filename_list_predicted) #weights files to mosaic 
    
    #writeLines(unlist(list_weights_m),con=filename_list_mosaics_weights_m) #weights files to mosaic 
    #writeLines(unlist(list_weights_prod_m),con=filename_list_mosaics_prod_weights_m) #prod weights files to mosaic
      
    #browser()
  
    #out_mosaic_name_weights_m <- r_weights_sum_raster_name <- file.path(out_dir,paste("r_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix,".tif",sep=""))
    #out_mosaic_name_prod_weights_m <- r_weights_sum_raster_name <- file.path(out_dir,paste("r_prod_weights_sum_m_",mosaic_method,"_weighted_mean_",out_suffix,".tif",sep=""))
    out_mosaic_name_predicted_m  <- file.path(out_dir_str,paste("r_overlap_sum_m_",out_suffix_str_tmp,"_tmp",".tif",sep=""))
    rast_ref_name <- infile_mask
    mosaic_python <- "/nobackupp6/aguzman4/climateLayers/sharedCode/"
    rast_ref_name <- infile_mask
    #python /nobackupp6/aguzman4/climateLayers/sharedCode//gdal_merge_sum.py --config GDAL_CACHEMAX=1500 --overwrite=TRUE -o /nobackupp8/bparmen1/climateLayers/out
    mosaic_overlap_tiles_obj <- mosaic_python_merge(NA_flag_val=NA_flag_val,
                                                module_path=mosaic_python,
                                                module_name="gdal_merge_sum.py",
                                                input_file=filename_list_predicted,
                                                out_mosaic_name=out_mosaic_name_predicted_m,
                                                raster_ref_name = rast_ref_name) ##if NA, not take into account
    r_overlap_raster_name <- mosaic_overlap_tiles_obj$out_mosaic_name
    cmd_str1 <-   mosaic_overlap_tiles_obj$cmd_str

    r_overlap <- raster(r_overlap_raster_name)
    r_mask <- raster(infile_mask)
    out_mosaic_name_overlap_masked  <- file.path(out_dir_str,paste("r_overlap_sum_masked_",region_name,"_",out_suffix,".tif",sep=""))

    r_overlap_m <- mask(r_overlap,
                  mask=r_mask,
                  filename=out_mosaic_name_overlap_masked,
                  datatype=data_type_str,
                  #datatype=data_type,
                  options=c("COMPRESS=LZW"),#compress tif
                  overwrite=TRUE,
                  NAflag=NA_flag_val)

    #r_overlap_m <- mask(r_overlap,r_mask,filename=out_mosaic_name_overlap_masked,overwrite=T)
    #plot(r_overlap_m)
    #plot(spdf_tiles_test,add=T,border="green",usePolypath = FALSE) #added usePolypath following error on brige and NEX
  
    #r_table <- ratify(r_overlap_m) # build the Raster Attibute table
    #rat <- levels(r_table)[[1]]#get the values of the unique cell frot the attribute table
    #rat$legend <- paste0("tile_",1:26)
    tb_freq_overlap <- as.data.frame(freq(r_overlap_m))
    write.table(tb_freq_overlap,file=paste0("tb_freq_overlap_",out_suffix,".txt"))
    
    res_pix <- 800
    #res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
  
    png_filename_maximum_overlap <-  file.path(out_dir,paste("Figure_maximum_overlap_",region_name,"_",out_suffix,".png",sep =""))
    title_str <-  paste("Maximum overlap: Number of predicted pixels for ",variable_name, sep = "")
  
    png(filename=png_filename_maximum_overlap,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
    #my_col=c('blue','red','green')
    my_col <- rainbow(length(tb_freq_overlap$value))
    plot(r_overlap_m,col=my_col,legend=F,box=F,axes=F,main=title_str)
    legend(x='topright', legend =tb_freq_overlap$value,fill = my_col,cex=0.8)
    dev.off()
    
    #browser()
    
  }else{ #if raster_overalp==FALSE
    out_mosaic_name_overlap_masked <- NULL
    tb_freq_overlap <- NULL
    png_filename_maximum_overlap <- NULL
  }

  ########################
  #### Step 3: combine overlap information and number of predictions by day
  ##Now loop through every day if missing then generate are raster showing map of number of prediction
  
  #df_missing_tiles_day <- subset(df_missing,tot_missing > 0) #already above
  browser()
  
  if(raster_pred==TRUE & (nrow(df_missing_tiles_day)>0)){
    
    #browser()
    
    ### now assign id and match extent for tiles
    
    lf_files <- unlist(list_predicted)
    rast_ref_name <- infile_mask
    rast_ref <- rast_ref_name
    
    ##Maching resolution is probably only necessary for the r mosaic function
    #Modify later to take into account option R or python...
    list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix_str_tmp,out_dir_str)
    names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")
    

    #r_test <- raster(raster_match(1,list_param_raster_match))
    list_tiles_predicted_m <- unlist(mclapply(1:length(lf_files),
                                              FUN=raster_match,list_param=list_param_raster_match,
                                              mc.preschedule=FALSE,mc.cores = num_cores))                           
    
    extension_str <- extension(lf_files)
    raster_name_tmp <- gsub(extension_str,"",basename(lf_files))
    out_suffix_str <- paste0(region_name,"_",out_suffix)
    raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_","masked_",out_suffix_str,file_format,sep=""))
    
    
    #writeRaster(r, NAflag=NA_flag_val,filename=raster_name,overwrite=TRUE)  
    
    #r_stack <- stack(list_tiles_predicted_m)
    list_mask_out_file_name <- raster_name
    list_tiles_predicted_masked <- unlist(mclapply(1:length(list_tiles_predicted_m),
                                                   FUN=function(i){mask(raster(list_tiles_predicted_m[i]),
                                                                        r_mask,filename=list_mask_out_file_name[i],
                                                                        overwrite=T,
                                                                        datatype=data_type_str,                  
                                                                        options=c("COMPRESS=LZW"))},
                                                   mc.preschedule=FALSE,
                                                   mc.cores = num_cores))                         
    #r_stack_masked <- mask(r, m2) #, maskvalue=TRUE)
    
    #Debugged on 12/16
    #r_tiles_stack <- stack(list_tiles_predicted_masked)
    #names(r_tiles_stack) <- basename(in_dir_reg) #this does not work, X. is added to the name, use list instead
  
    #names(list_tiles_predicted_masked) <- basename(in_dir_reg)
    #df_missing_tiles_day <- subset(df_missing,tot_missing > 0)
    #r_tiles_s <- r_tiles_stack
    #names_tiles <- basename(in_dir_reg)
  

    list_param_generate_raster_number_pred <- list(list_tiles_predicted_masked,df_missing_tiles_day,r_overlap_m,
                                                 num_cores,region_name,data_type,scaling,python_bin,
                                                 plotting_figures,
                                                 NA_flag_val,out_suffix,out_dir)
  
    names(list_param_generate_raster_number_pred) <- c("list_tiles_predicted_masked","df_missing_tiles_day","r_overlap_m",
                                                     "num_cores","region_name","data_type","scaling","python_bin",
                                                     "plotting_figures",
                                                      "NA_flag_val","out_suffix","out_dir")
  
    #function_product_assessment_part0_functions <- "global_product_assessment_part0_functions_11152016b.R"
    #source(file.path(script_path,function_product_assessment_part0_functions)) #source all functions used in this script 

    #undebug(generate_raster_number_of_prediction_by_day)
    #4.51pm
    #browser()
    #5.10pm
    #test_number_pix_predictions <- generate_raster_number_of_prediction_by_day(1,list_param=list_param_generate_raster_number_pred)
    obj_number_pix_predictions <- mclapply(1:nrow(df_missing_tiles_day),
                                        FUN=generate_raster_number_of_prediction_by_day,
                                        list_param=list_param_generate_raster_number_pred,
                                        mc.preschedule=FALSE,
                                        mc.cores = num_cores)
    
  }else{
    obj_number_pix_predictions <- NULL
  }
    
  browser()
  #Delete temporary files : Fix this part later...
  #rasterOptions(), find where tmp dir are stored
  
  if(tmp_files==F){ #if false...delete all files with "_tmp"
      lf_tmp <- list.files(path=out_dir,pattern=".*._tmp.*")
      #lf_tmp <- unlist(lf_accuracy_training_raster)
      ##now delete temporary files...
      file.remove(lf_tmp)
  }
  
  predictions_tiles_missing_obj <- list(df_lf_tiles_time_series,df_missing_tiles_day,out_mosaic_name_overlap_masked,
                                        tb_freq_overlap,png_filename_maximum_overlap,obj_number_pix_predictions)
  names(predictions_tiles_missing_obj) <- c("df_lf_tiles_time_series","df_missing_tiles_day","raster_name_overlap",
                                            "tb_freq_overlap","png_maximum_overlap","obj_number_pix_predictions")
  
  
  predictions_tiles_missing_obj_filename <- file.path(out_dir,paste("obj_predictions_tiles_missing_fun_",interpolation_method,y_var_name,region_name,out_suffix,".RData",sep=""))
  save(predictions_tiles_missing_obj,file=predictions_tiles_missing_obj_filename)

  return(predictions_tiles_missing_obj)
}


############################ END OF SCRIPT ##################################