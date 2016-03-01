##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part5 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Combining tables and figures for individual runs for years and tiles.
#This script complements part1 and part2 of the accuracy assessment and group tables and outputs 
#from run of accuracy assessement generated earlier.
#It focuses on additional analyses, figures, tables of accuracy values. In particular, extreme values by tiles.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 02/25/2016  
#MODIFIED ON: 02/28/2016            
#Version: 5
#PROJECT: Environmental Layers project     
#COMMENTS: Initial commit, script based on part 3 of assessment, will be modified further for overall assessment 
#TODO:
#1) 

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rmd user:aguzman4:rwx /nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015

#################################################################################################

#### FUNCTION USED IN SCRIPT

#function_analyses_paper1 <-"contribution_of_covariates_paper_interpolation_functions_07182014.R" #first interp paper
#function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_08132014.R"
#function_global_run_assessment_part2 <- "global_run_scalingup_assessment_part2_functions_0923015.R"

############################################
#### Parameters and constants  

#on ATLAS
#in_dir1 <- "/data/project/layers/commons/NEX_data/test_run1_03232014/output" #On Atlas
#parent output dir : contains subset of the data produced on NEX
#in_dir1 <- "/data/project/layers/commons/NEX_data/output_run6_global_analyses_09162014/output20Deg2"
# parent output dir for the curent script analyes
#out_dir <-"/data/project/layers/commons/NEX_data/output_run3_global_analyses_06192014/" #On NCEAS Atlas
# input dir containing shapefiles defining tiles
#in_dir_shp <- "/data/project/layers/commons/NEX_data/output_run5_global_analyses_08252014/output/subset/shapefiles"

#On NEX
#contains all data from the run by Alberto
#in_dir1 <- " /nobackupp6/aguzman4/climateLayers/out_15x45/" #On NEX
#parent output dir for the current script analyes
#out_dir <- "/nobackup/bparmen1/" #on NEX
#in_dir_shp <- "/nobackupp4/aguzman4/climateLayers/output4/subset/shapefiles/"

#in_dir <- "/data/project/layers/commons/NEX_data/reg4_assessment"
#list_in_dir_run <-
#in_dir_list <-  c("output_run_global_analyses_pred_2009_reg4","output_run_global_analyses_pred_2010_reg4",
#                  "output_run_global_analyses_pred_2011_reg4","output_run_global_analyses_pred_2012_reg4",
#                  "output_run_global_analyses_pred_2013_reg4","output_run_global_analyses_pred_2014_reg4")
#in_dir_list_filename <- "/data/project/layers/commons/NEX_data/reg4_assessment/stage6_reg4_in_dir_list_02072016.txt"
#in_dir <- "" #PARAM 0
#y_var_name <- "dailyTmax" #PARAM1
#interpolation_method <- c("gam_CAI") #PARAM2
#out_suffix <- "global_analyses_overall_assessment_reg4_02072016"
#out_suffix <- "output_run10_1000x3000_global_analyses_02102015"
#out_suffix <- "run10_1500x4500_global_analyses_pred_1992_10052015" #PARAM3
#out_dir <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1992_10052015" #PARAM4
#create_out_dir_param <- TRUE #PARAM 5
#mosaic_plot <- FALSE #PARAM6
#if daily mosaics NULL then mosaicas all days of the year
#day_to_mosaic <- c("19920101","19920102","19920103") #PARAM7
#CRS_WGS84 <-    CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84 #CONSTANT1
#CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
#proj_str<- CRS_WGS84 #PARAM 8 #check this parameter
#file_format <- ".rst" #PARAM 9
#NA_value <- -9999 #PARAM10
#NA_flag_val <- NA_value
#multiple_region <- TRUE #PARAM 12
#region_name <- "world" #PARAM 13
#countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp" #PARAM 13, copy this on NEX too
#plot_region <- TRUE
#num_cores <- 6 #PARAM 14
#region_name <- c("reg4") #reference region to merge if necessary, if world all the regions are together #PARAM 16
#use previous files produced in step 1a and stored in a data.frame
#df_assessment_files <- "df_assessment_files_reg4_1984_run_global_analyses_pred_12282015.txt" #PARAM 17
#threshold_missing_day <- c(367,365,300,200) #PARM18

#list_param_run_assessment_plottingin_dir <- list(in_dir,y_var_name, interpolation_method, out_suffix, 
#                      out_dir, create_out_dir_param, mosaic_plot, proj_str, file_format, NA_value,
#                      multiple_region, countries_shp, plot_region, num_cores, 
#                      region_name, df_assessment_files, threshold_missing_day) 

#names(list_param_run_assessment_plottingin_dir) <- c("in_dir","y_var_name","interpolation_method","out_suffix", 
#                      "out_dir","create_out_dir_param","mosaic_plot","proj_str","file_format","NA_value",
#                      "multiple_region","countries_shp","plot_region","num_cores", 
#                      "region_name","df_assessment_files","threshold_missing_day") 

#run_assessment_plotting_prediction_fun(list_param_run_assessment_plottingin_dir) 

run_assessment_combined_region_plotting_prediction_fun <-function(list_param_run_assessment_plotting){
  
  ####
  #1) in_dir: input directory containing data tables and shapefiles for plotting #PARAM 0
  #2) y_var_name : variables being predicted e.g. dailyTmax #PARAM1
  #3) interpolation_method: method used #c("gam_CAI") #PARAM2
  #4) out_suffix: output suffix #PARAM3
  #5) out_dir  #
  #6) create_out_dir_param # FALSE #PARAM 5
  #7) mosaic_plot  #FALSE #PARAM6
  #8) proj_str # projection/coordinates system e.g. CRS_WGS84 #PARAM 8 #check this parameter
  #9) file_format #".rst" #PARAM 9
  #10) NA_value #-9999 #PARAM10
  #11) multiple_region  # <- TRUE #PARAM 12
  #12) countries_shp  #<- "world" #PARAM 13
  #13) plot_region  #<- TRUE
  #14) num_cores <- number of cores used # 6 #PARAM 14
  #15) region_name  #<- c("reg4"), world if full assessment #reference region to merge if necessary #PARAM 16
  #16) df_assessment_files  #PARAM 16
  #17) threshold_missing_day  #PARM18
  #18) year_predicted
  
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
  
  ####### Function used in the script #######
  
  #script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts"
  #function_assessment_part2_functions <- "global_run_scalingup_assessment_part2_functions_01032016.R"
  #source(file.path(script_path,function_assessment_part2_functions)) #source all functions used in this script 

  ####### PARSE INPUT ARGUMENTS/PARAMETERS #####
  in_dir_list_filename <- list_param_run_assessment_plotting$in_dir_list_filename #PARAM 0
  in_dir <- list_param_run_assessment_plotting$in_dir #PARAM 1
  y_var_name <- list_param_run_assessment_plotting$y_var_name #PARAM2
  interpolation_method <- list_param_run_assessment_plotting$interpolation_method #c("gam_CAI") #PARAM3
  out_suffix <- list_param_run_assessment_plotting$out_suffix #PARAM4
  out_dir <- list_param_run_assessment_plotting$out_dir # PARAM5
  create_out_dir_param <- list_param_run_assessment_plotting$create_out_dir_param # FALSE #PARAM 6
  mosaic_plot <- list_param_run_assessment_plotting$mosaic_plot #FALSE #PARAM7
  proj_str<- list_param_run_assessment_plotting$proj_str #CRS_WGS84 #PARAM 8 #check this parameter
  file_format <- list_param_run_assessment_plotting$file_format #".rst" #PARAM 9
  NA_flag_val <- list_param_run_assessment_plotting$NA_flag_val #-9999 #PARAM10
  multiple_region <- list_param_run_assessment_plotting$multiple_region # <- TRUE #PARAM 11
  countries_shp <- list_param_run_assessment_plotting$countries_shp #<- "world" #PARAM 12
  plot_region <- list_param_run_assessment_plotting$plot_region # PARAM13 
  num_cores <- list_param_run_assessment_plotting$num_cores # 6 #PARAM 14
  region_name <- list_param_run_assessment_plotting$region_name #<- "world" #PARAM 15
  #df_assessment_files_name <- list_param_run_assessment_plotting$df_assessment_files_name #PARAM 16
  threshold_missing_day <- list_param_run_assessment_plotting$threshold_missing_day #PARM17
  year_predicted <- list_param_run_assessment_plotting$year_predicted
 
  NA_value <- NA_flag_val 
  metric_name <- "rmse" #to be added to the code later...
  
  ##################### START SCRIPT #################
  
  ####### PART 1: Read in data ########
  out_dir <- in_dir
  if(create_out_dir_param==TRUE){
    out_dir <- create_dir_fun(out_dir,out_suffix)
    setwd(out_dir)
  }else{
    setwd(out_dir) #use previoulsy defined directory
  }

  setwd(out_dir)
  
  list_outfiles <- vector("list", length=35) #collect names of output files, this should be dynamic?
  list_outfiles_names <- vector("list", length=35) #collect names of output files
  counter_fig <- 0 #index of figure to collect outputs
  
  #i <- year_predicted
  ###Table 1: Average accuracy metrics
  ###Table 2: daily accuracy metrics for all tiles

  in_dir_list <- as.list(read.table(in_dir_list_filename,stringsAsFactors=F)[,1])
  
  ##Read in data list from in_dir_list
  list_tb_fname <- list.files(path=file.path(in_dir,in_dir_list),"tb_diagnostic_v_NA_.*.txt",full.names=T)
  list_df_fname <- list.files(path=file.path(in_dir,in_dir_list),"df_tile_processed_.*..txt",full.names=T)
  list_summary_metrics_v_fname <- list.files(path=file.path(in_dir,in_dir_list),"summary_metrics_v2_NA_.*.txt",full.names=T)
  list_tb_s_fname <- list.files(path=file.path(in_dir,in_dir_list),"tb_diagnostic_s_NA.*.txt",full.names=T)
  list_tb_month_s_fname <- list.files(path=file.path(in_dir,in_dir_list),"tb_month_diagnostic_s.*.txt",full.names=T)
  list_data_month_s_fname <- list.files(path=file.path(in_dir,in_dir_list),"data_month_s.*.txt",full.names=T)
  list_data_s_fname <- list.files(path=file.path(in_dir,in_dir_list),"data_day_s.*.txt",full.names=T)
  list_data_v_fname <- list.files(path=file.path(in_dir,in_dir_list),"data_day_v.*.txt",full.names=T)
  list_pred_data_month_info_fname <- list.files(path=file.path(in_dir,in_dir_list),"pred_data_month_info.*.txt",full.names=T)
  list_pred_data_day_info_fname <- list.files(path=file.path(in_dir,in_dir_list),"pred_data_day_info.*.txt",full.names=T)
  
  #need to fix this !! has all of the files in one list (for a region)
  #list_shp <- list.files(path=file.path(in_dir,file.path(in_dir_list,"shapefiles")),"*.shp",full.names=T)

  ## Step 2: only read what is necessary at this stage...
  list_tb <- lapply(list_tb_fname,function(x){read.table(x,stringsAsFactors=F,sep=",")})
  tb <- do.call(rbind,list_tb)
  list_tb_s <- lapply(list_tb_s_fname,function(x){read.table(x,stringsAsFactors=F,sep=",")})
  tb_s <- do.call(rbind,list_tb_s)
  
  list_df_tile_processed <- lapply(list_df_fname,function(x){read.table(x,stringsAsFactors=F,sep=",")})
  df_tile_processed <- do.call(rbind,list_df_tile_processed)  
  list_summary_metrics_v <- lapply(list_summary_metrics_v_fname,function(x){read.table(x,stringsAsFactors=F,sep=",")})
  summary_metrics_v <- do.call(rbind,list_summary_metrics_v)  

  list_tb_month_s <- lapply(list_tb_month_s_fname,function(x){read.table(x,stringsAsFactors=F,sep=",")})
  tb_month_s <- do.call(rbind,list_tb_month_s)  
  
  list_tb_data_s <- lapply(list_data_s_fname,function(x){read.table(x,stringsAsFactors=F,sep=",")})
  tb_data_s <- do.call(rbind,list_tb_data_s) 

  list_tb_data_v <- lapply(list_data_v_fname,function(x){read.table(x,stringsAsFactors=F,sep=",")})
  tb_data_v <- do.call(rbind,list_tb_data_v) 

  ##Stop added
  ##Screen for non shapefiles tiles due to dir
  df_tile_processed <- df_tile_processed[!is.na(df_tile_processed$shp_files),] 
  
  #add column for tile size later on!!!
  
  tb$pred_mod <- as.character(tb$pred_mod)
  summary_metrics_v$pred_mod <- as.character(summary_metrics_v$pred_mod)
  summary_metrics_v$tile_id <- as.character(summary_metrics_v$tile_id)
  df_tile_processed$tile_id <- as.character(df_tile_processed$tile_id)
  
  tb_month_s$pred_mod <- as.character(tb_month_s$pred_mod)
  tb_month_s$tile_id<- as.character(tb_month_s$tile_id)
  tb_s$pred_mod <- as.character(tb_s$pred_mod)
  tb_s$tile_id <- as.character(tb_s$tile_id)
  
  #multiple regions? #this needs to be included in the previous script!!!
  #if(multiple_region==TRUE){
  df_tile_processed$reg <- as.character(df_tile_processed$reg)
  tb <- merge(tb,df_tile_processed,by="tile_id")
  tb_s <- merge(tb_s,df_tile_processed,by="tile_id")
  tb_month_s<- merge(tb_month_s,df_tile_processed,by="tile_id")
  summary_metrics_v <- merge(summary_metrics_v,df_tile_processed,by="tile_id")
  #test <- merge(summary_metrics_v,df_tile_processed,by="tile_id",all=F)
  #duplicate columns...need to be cleaned up later
  try(tb$year_predicted <- tb$year_predicted.x)
  try(tb$reg <- tb$reg.x)
  try(summary_metrics_v$year_predicted <- summary_metrics_v$year_predicted.x)
  try(summary_metrics_v$reg <- summary_metrics_v$reg.x)  
  try(summary_metrics_v$lat <- summary_metrics_v$lat.x)
  try(summary_metrics_v$lon <- summary_metrics_v$lon.x)

  ############ PART 2: PRODUCE FIGURES ################
  
  
  tb_subset <- subset(tb,pred_mod=="mod1")
  
  #sqrt(var(tb_subset$rmse))
  
  std_dev_val <- sqrt(var(tb_subset[[metric_name]])) #e.g. rmse
  mean_val <- mean(tb_subset[[metric_name]])
  median_val <- median(tb_subset[[metric_name]])
  max_val <- max(tb_subset[[metric_name]])
  min_val <- min(tb_subset[[metric_name]])
  n_val <- length(tb_subset[[metric_name]])
    
  #mode(tb_subset[[metric_name]])
  
  stat_name <- c("std_dev","mean","median","max","min","n")
  stat_val <- c(std_dev_val,mean_val,median_val,max_val,min_val,n_val)
  df_stat <- data.frame(stat_name=stat_name,stat_val=stat_val)
  
  model_name <- c("mod1","mod_kr")
  threshold_val <- c(5,10,20,50)
  
  ###########################
  ### Figure 1: plot location of the study area with tiles processed
  
  #df_tiled_processed <- na.omit(df_tile_processed) #remove other list of folders irrelevant
  #list_shp_reg_files <- df_tiled_processed$shp_files
  
  #list_shp_reg_files<- as.character(df_tile_processed$shp_files) #this could be the solution!!
  #Use melt!! quick solution
  list_shp_reg_files <- as.character(basename(unique(df_tile_processed$shp_files))) #this could be the solution!!
  list_tile_id <- as.character((unique(df_tile_processed$tile_id))) #this is in the order of appearance
  #list_tile_lat <- as.character((unique(df_tile_processed$lat))) #this is in the order of appearance
  #list_tile_lon <- as.character((unique(df_tile_processed$lon))) #this is in the order of appearance

  #df_tile_processed$shp_files2 <- basename(df_tile_processed$shp_files)
  df_tiles_reg <- data.frame(shp_files=(list_shp_reg_files),tile_id=list_tile_id)
  df_tiles_reg$tile_id <- as.character(df_tiles_reg$tile_id)
  #,lat=list_tile_lat,lon=list_tile_lon)
  #cast(df_tile_processed, shp_files ~ tile_id+lat+lon, mean, value = 'income')

  ### Potential function starts here:
  #function(in_dir,out_dir,list_shp_reg_files,title_str,region_name,num_cores,out_suffix,out_suffix)
  
  ### First get background map to display where study area is located
  #can make this more general later on..should have this already in a local directory on Atlas or NEX!!!!
  
  #http://www.diva-gis.org/Data
  #countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp"
  path_to_shp <- dirname(countries_shp)
  layer_name <- sub(".shp","",basename(countries_shp))
  reg_layer <- readOGR(path_to_shp, layer_name)
  #proj4string(reg_layer) <- CRS_locs_WGS84
  #reg_shp<-readOGR(dirname(list_shp_reg_files[[i]]),sub(".shp","",basename(list_shp_reg_files[[i]])))
  
  centroids_pts <- vector("list",length(list_shp_reg_files))
  shps_tiles <- vector("list",length(list_shp_reg_files))
  #collect info: read in all shapfiles
  #This is slow...make a function and use mclapply??
  #/data/project/layers/commons/NEX_data/output_run6_global_analyses_09162014/shapefiles
  
  in_dir_shp <- file.path(in_dir,in_dir_list[[1]],"shapefiles") #this should be set as a input parameter!!!
  for(i in 1:length(list_shp_reg_files)){
    #path_to_shp <- dirname(list_shp_reg_files[[i]])
    path_to_shp <- in_dir_shp
    layer_name <- sub(".shp","",basename(list_shp_reg_files[[i]]))
    shp1 <- try(readOGR(path_to_shp, layer_name)) #use try to resolve error below
    #shp_61.0_-160.0
    #Geographical CRS given to non-conformant data: -186.331747678
    
    #shp1<-readOGR(dirname(list_shp_reg_files[[i]]),sub(".shp","",basename(list_shp_reg_files[[i]])))
    if (!inherits(shp1,"try-error")) {
      pt <- gCentroid(shp1)
      centroids_pts[[i]] <- pt
    }else{
      pt <- shp1
      centroids_pts[[i]] <- pt
    }
    shps_tiles[[i]] <- shp1
    #centroids_pts[[i]] <- centroids
  }
  
  #fun <- function(i,list_shp_files)
  #coord_names <- c("lon","lat")
  #l_ras#t <- rasterize_df_fun(test,coord_names,proj_str,out_suffix=out_suffix,out_dir=".",file_format,NA_flag_val,tolerance_val=0.000120005)
  
  #remove try-error polygons...we loose three tiles because they extend beyond -180 deg
  tmp <- shps_tiles
  shps_tiles <- remove_errors_list(shps_tiles) #[[!inherits(shps_tiles,"try-error")]]
  #shps_tiles <- tmp
  length(tmp)-length(shps_tiles) #number of tiles with error message
  
  tmp_pts <- centroids_pts 
  centroids_pts <- remove_errors_list(centroids_pts) #[[!inherits(shps_tiles,"try-error")]]
  #centroids_pts <- tmp_pts 
  
  df_pts <- as.data.frame(do.call(rbind,tmp_pts))
  #df_pts <- cbind(df_pts,df_tiles_reg)
  df_tiles_reg <- cbind(df_pts,df_tiles_reg) #(shp_files=(list_shp_reg_files),tile_id=list_tile_id)
  df_tiles_reg$id <- as.numeric(unlist(lapply(strsplit(df_tiles_reg$tile_id,"_"),FUN=function(x){x[2]})))
  
  coordinates(df_tiles_reg)<- cbind(df_tiles_reg$x,df_tiles_reg$y)
  
  #plot info: with labels
  res_pix <-1200
  col_mfrow <- 1 
  row_mfrow <- 1
  
  png(filename=paste("Figure1a_tile_processed_region_",region_name,"_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  plot(reg_layer)
  #Add polygon tiles...
  for(i in 1:length(shps_tiles)){
    shp1 <- shps_tiles[[i]]
    pt <- centroids_pts[[i]]
    #if(!inherits(shp1,"try-error")){
    #  plot(shp1,add=T,border="blue")
    #  #plot(pt,add=T,cex=2,pch=5)
    #  label_id <- df_tile_processed$tile_id[i]
    #  text(coordinates(pt)[1],coordinates(pt)[2],labels=i,cex=1.3,font=2,col=c("red"))
    #}
    #to be able to run on NEX set font and usePolypath, maybe add option NEX?
    if(!inherits(shp1,"try-error")){
      plot(shp1,add=T,border="blue",usePolypath = FALSE) #added usePolypath following error on brige and NEX
      #plot(pt,add=T,cex=2,pch=5)
      label_id <- df_tile_processed$tile_id[i]
      text(coordinates(pt)[1],coordinates(pt)[2],labels=i,cex=1.3,font=2,col=c("red"),family="HersheySerif")
    }

  }
  #title(paste("Tiles ", tile_size,region_name,sep=""))
  
  dev.off()
  
  res_pix <-1200
  col_mfrow <- 1 
  row_mfrow <- 1
  
  png(filename=paste("Figure1b_tile_processed_centroids_region_",region_name,"_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  #plot(reg_layer)
  #Add polygon tiles...
  #title(paste("Tiles ", tile_size,region_name,sep=""))
  #plot(df_tiles_reg,add=T,pch=2)
  #label_id <- df_tiles_reg$id
  #text(coordinates(df_tiles_reg)[1],coordinates(df_tiles_reg)[2],labels=i,cex=1.3,font=2,col=c("red"),family="HersheySerif")

  p_shp <- spplot(reg_layer,"ISO3" ,col.regions=NA, col="black") #ok problem solved!!
  #title("(a) Mean for 1 January")
  df_tiles_reg$lab <- 1
  sl1 <- list('sp.points',df_tiles_reg, pch=19, cex=.8, col='red')
  sl2 <- list('sp.pointLabel', df_tiles_reg, label=list_id,
            cex=2.4, font=2,col='red',col.regions="red",
            fontfamily='Palatino') #Add labels at centroids
 
  p <- spplot(df_tiles_reg,"id",main=paste("Tile id processed",sep=""),sp.layout=list(sl1, sl2))
  #spplot(meuse.grid["dist"], col.regions=myCols, sp.layout=list(sl1, sl2)
  #p <- spplot(df_tiles_reg,"lab",main=paste("Tile id processed",sep=""))
  p1 <- p+p_shp
  try(print(p1)) #error raised if number of missing values below a threshold does not exist
  #ltext(coordinates(df_tiles_reg)[,1],coordinates(df_tiles_reg)[,2],labels=list_tile_id)

  #list_id <- df_tiles_reg$id
  dev.off()

  list_outfiles[[counter_fig+1]] <- paste("Figure1_tile_processed_region_",region_name,"_",out_suffix,".png",sep="")
  counter_fig <- counter_fig+1
  #this will be changed to be added to data.frame on the fly
  r1 <-c("figure_1","Tiles processed for the region",NA,NA,region_name,year_predicted,list_outfiles[[1]]) 


  ######################
  ### Figure 2: Number of predictions: daily and monthly
  
  ## Figure 2a
 
  #Plot location of extremes and select them for further analyses?


  j<-1 #for model name 1,mod1
  #model_name <- c("mod1","mod_kr")
  #threshold_val <- c(5,10,20,50)
  
  for(i in 1:length(threshold_val)){
    
    tb_tmp <- try(subset(tb_subset,tb_subset[[metric_name]]>threshold_val[i]))
    hist(tb_tmp$rmse)
    fig_filename1 <-paste("Figure2a_barplot_extremes_val_centroids_tile_",model_name[j],"_",threshold_val[i],
                       "_",out_suffix,".png",sep="")
    fig_filename2 <- paste("Figure2b_ac_metrics_extremes_map_centroids_tile_",model_name[j],"_",threshold_val[i],
                       "_",out_suffix,".png",sep="")
    list_outfiles[[counter_fig+i]] <- fig_filename1
    list_outfiles[[counter_fig+i]] <- fig_filename2
        
    if(nrow(tb_subset)>0){
      
      df_extremes <- as.data.frame(table(tb_tmp$tile_id))
      names(df_extremes)<- c("tile_id","freq_extremes")
      tb_sorted <- merge(tb_tmp,df_extremes,"tile_id")
      tb_sorted <- arrange(tb_sorted,desc(freq_extremes)) #[,c("pred_mod","rmse","mae","tile_id")]
      coordinates(tb_sorted) <- c("lon","lat")
      df_extremes <- arrange(df_extremes,desc(freq_extremes))
    
      fig_filename <- paste("Figure2a_barplot_extremes_val_centroids_tile_",model_name[j],"_",threshold_val[i],
                       "_",out_suffix,".png",sep="")
      
      #res_pix <- 1200
      res_pix <- 960
      col_mfrow <- 1
      row_mfrow <- 1
      #only mod1 right now
      png(filename=fig_filename,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)

      #barplot(df_ac_mod$rmse, names.arg=x)
      barplot(df_extremes$freq_extremes,main=paste("Extremes threshold val for ",threshold_val[i],"deg C",sep=""),
            ylab="frequency",xlab="tile_id",names.arg=df_extremes$tile_id,las=2)    
      #barplot(table(tb_subset$tile_id),main=paste("Extremes threshold val for ",threshold_val[i],sep=""),
       #        ylab="frequency",xlab="tile_id",las=2)
      dev.off()

      #test<-(subset(tb,tb$tile_id==unique(tb_subset$tile_id)))
      #df_ac_mod <- arrange(as.data.frame(ac_mod),desc(rmse))[,c("pred_mod","rmse","mae","tile_id")]
    
      #plot top three, then all,and histogram...make this a function...
      #list_df_ac_mod[[i]] <- arrange(as.data.frame(ac_mod),desc(rmse))[,c("rmse","mae","tile_id")]

      fig_filename <- paste("Figure2b_ac_metrics_extremes_map_centroids_tile_",model_name[j],"_",threshold_val[i],
                       "_",out_suffix,".png",sep="")

      #res_pix <- 1200
      res_pix <- 960
      col_mfrow <- 1
      row_mfrow <- 1
      #only mod1 right now
      png(filename=fig_filename,
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
      model_name[j]
    
      #p_shp <- layer(sp.polygons(reg_layer, lwd=1, col='black'))
      p_shp <- spplot(reg_layer,"ISO3" ,col.regions=NA, col="black") #ok problem solved!!
      #title("(a) Mean for 1 January")
      p <- bubble(tb_sorted,"freq_extremes",main=paste("Extremes per tile and by ",model_name[j]," for ",
                                                                threshold_val[i]))
      p1 <- p+p_shp
      try(print(p1)) #error raised if number of missing values below a threshold does not exist
      dev.off()

    }

    #list_outfiles[[counter_fig+i]] <- fig_filename
  }
  counter_fig <- counter_fig+length(threshold_missing_day) #currently 4 days...

  r18 <-c("figure_7","Number of missing days threshold1 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[18]])
  r19 <-c("figure_7","Number of missing days threshold2 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[19]])  
  r20 <-c("figure_7","Number of missing days threshold3 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[20]])
  r21 <-c("figure_7","Number of missing days threshold4 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[21]])  

  ##### Figure 3 ###
  
  test<- subset(tb_subset,tb_subset$tile_id=="tile_14")
  test2 <- aggregate(rmse~date,test,min)
  idx <- test$date #transform this format...
  d_z <- zoo(test,idx) #make a time series ...
  #add horizontal line...

  plot(test$rmse,type="b")
  unique(test$year_predicted)
   unique(tb$year_predicted)
   
  ######################################################
  ##### Prepare objet to return ####

  assessment_obj <- list(list_df_assessment_files, df_assessment_figures_files)
  names(assessment_obj) <- c("df_assessment_files", "df_assessment_figures_files")
  ## Prepare list of files to return...
  return(assessment_obj)
 
}
  
##################### END OF SCRIPT ######################
