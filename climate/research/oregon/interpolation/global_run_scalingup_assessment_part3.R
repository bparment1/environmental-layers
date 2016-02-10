##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part3 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Combining tables and figures for individual runs for years and tiles.
#This script complements part1 and part2 of the accuracy assessment and group tables and outputs 
#from run of accuracy assessement generated earlier.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 03/23/2014  
#MODIFIED ON: 02/09/2016            
#Version: 5
#PROJECT: Environmental Layers project     
#COMMENTS: Initial commit, script based on part 2 of assessment, will be modified further for overall assessment 
#TODO:
#1) Add plot broken down by year and region 
#2) Modify code for overall assessment accross all regions and year
#3) Clean up

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
  
  list_outfiles <- vector("list", length=31) #collect names of output files, this should be dynamic?
  list_outfiles_names <- vector("list", length=31) #collect names of output files
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
  
  ###########################
  ### Figure 1: plot location of the study area with tiles processed
  
  #df_tiled_processed <- na.omit(df_tile_processed) #remove other list of folders irrelevant
  #list_shp_reg_files <- df_tiled_processed$shp_files
  
  #list_shp_reg_files<- as.character(df_tile_processed$shp_files) #this could be the solution!!
  list_shp_reg_files <- as.character(basename(unique(df_tile_processed$shp_files))) #this could be the solution!!
  #the shapefiles must be copied in the proper folder!!!
  #list_shp_reg_files<- file.path(in_dir,in_dir_list[1],"shapefile",list_shp_reg_files)
  #list_shp_reg_files <- file.path("/data/project/layers/commons/NEX_data/",out_dir,
  #          as.character(df_tile_processed$tile_coord),"shapefiles",basename(list_shp_reg_files))
  #list_shp_reg_files <- file.path("/data/project/layers/commons/NEX_data/",out_dir,
                                  #"shapefiles",basename(list_shp_reg_files))
  
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
  
  #plot info: with labels
  res_pix <-1200
  col_mfrow <- 1 
  row_mfrow <- 1
  
  png(filename=paste("Figure1_tile_processed_region_",region_name,"_",out_suffix,".png",sep=""),
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
  
  list_outfiles[[counter_fig+1]] <- paste("Figure1_tile_processed_region_",region_name,"_",out_suffix,".png",sep="")
  counter_fig <- counter_fig+1
  #this will be changed to be added to data.frame on the fly
  r1 <-c("figure_1","Tiles processed for the region",NA,NA,region_name,year_predicted,list_outfiles[[1]]) 

  ###############
  ### Figure 2: boxplot of average accuracy by model and by tiles
  
  tb$tile_id <- factor(tb$tile_id, levels=unique(tb$tile_id))
  model_name <- as.character(unique(tb$pred_mod))
  
  ## Figure 2a
  for(i in  1:length(model_name)){
    
    res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
    fig_filename <-  paste("Figure2a_boxplot_with_oultiers_by_tiles_",model_name[i],"_",out_suffix,".png",sep="")
    png(filename=paste("Figure2a_boxplot_with_oultiers_by_tiles_",model_name[i],"_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    boxplot(rmse~tile_id,data=subset(tb,tb$pred_mod==model_name[i]))
    title(paste("RMSE per ",model_name[i]))
    
    dev.off()
    list_outfiles[[counter_fig+i]] <- fig_filename
  }
  counter_fig <- counter_fig + length(model_name)
  
  r2 <-c("figure_2a","Boxplot of accuracy with outliers by tiles","mod1",metric_name,region_name,year_predicted,list_outfiles[[2]]) 
  r3 <-c("figure_2a","boxplot of accuracy with outliers by tiles","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[3]])
  
  ## Figure 2b
  #with ylim and removing trailing...
  for(i in  1:length(model_name)){ #there are two models!!
    
    res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
    fig_filename <- paste("Figure2b_boxplot_scaling_by_tiles","_",model_name[i],"_",out_suffix,".png",sep="")
    png(filename=paste("Figure2b_boxplot_scaling_by_tiles","_",model_name[i],"_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    model_name <- unique(tb$pred_mod)
    boxplot(rmse~tile_id,data=subset(tb,tb$pred_mod==model_name[i])
            ,ylim=c(0,4),outline=FALSE)
    title(paste("RMSE per ",model_name[i]))
    dev.off()
    #we already stored one figure
    list_outfiles[[counter_fig+i]] <- fig_filename
  }
  counter_fig <- counter_fig + length(model_name)
  #bwplot(rmse~tile_id, data=subset(tb,tb$pred_mod=="mod1"))
  r4 <-c("figure_2b","Boxplot of accuracy with scaling by tiles","mod1",metric_name,region_name,year_predicted,list_outfiles[[4]])  
  r5 <-c("figure_2b","Boxplot of accuracy with scaling by tiles","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[5]])  
 
  ###############
  ### Figure 3: boxplot of average RMSE by model acrosss all tiles
  
  #Ok fixed..now selection of model but should also offer an option for using both models!! so make this a function!!
  for(i in  1:length(model_name)){ #there are two models!!
    ## Figure 3a
    res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
    
    png(filename=paste("Figure3a_boxplot_overall_region_with_oultiers_",model_name[i],"_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    #boxplot(rmse~pred_mod,data=tb)#,names=tb$pred_mod)
    boxplot(rmse~pred_mod,data=subset(tb,tb$pred_mod==model_name[i]))#,names=tb$pred_mod)
    title(paste("RMSE with outliers for all tiles: ",model_name[i],sep=""))
    dev.off()
    list_outfiles[[counter_fig+1]] <- paste("Figure3a_boxplot_overall_region_with_oultiers_",model_name[i],"_",out_suffix,".png",sep="")
    
    ## Figure 3b
    png(filename=paste("Figure3b_boxplot_overall_region_scaling_",model_name[i],"_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    #boxplot(rmse~pred_mod,data=tb,ylim=c(0,5),outline=FALSE)#,names=tb$pred_mod)
    boxplot(rmse~pred_mod,data=subset(tb,tb$pred_mod==model_name[i]),ylim=c(0,5),outline=FALSE)#,names=tb$pred_mod)
    #title("RMSE per model over all tiles")
    title(paste("RMSE with scaling for all tiles: ",model_name[i],sep=""))
    dev.off()
    list_outfiles[[counter_fig+2]] <- paste("Figure3b_boxplot_overall_region_scaling_",model_name[i],"_",out_suffix,".png",sep="")
  }
  counter_fig <- counter_fig + length(model_name)
  r6 <-c("figure_3a","Boxplot overall accuracy with outliers","mod1",metric_name,region_name,year_predicted,list_outfiles[[6]])  
  r7 <-c("figure_3b","Boxplot overall accuracy with scaling","mod1",metric_name,region_name,year_predicted,list_outfiles[[7]])  
  r8 <-c("figure_3a","Boxplot overall accuracy with outliers","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[8]])
  r9 <-c("figure_3b","Boxplot overall accuracy with scaling","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[9]])  

  ################ 
  ### Figure 4: plot predicted tiff for specific date per model
  ## Replace by break out by season?
  
  #y_var_name <-"dailyTmax"
  #index <-244 #index corresponding to Sept 1
  
  # if (mosaic_plot==TRUE){
  #   index  <- 1 #index corresponding to Jan 1
  #   date_selected <- "20100901"
  #   name_method_var <- paste(interpolation_method,"_",y_var_name,"_",sep="")
  # 
  #   pattern_str <- paste("mosaiced","_",name_method_var,"predicted",".*.",date_selected,".*.tif",sep="")
  #   lf_pred_list <- list.files(pattern=pattern_str)
  # 
  #   for(i in 1:length(lf_pred_list)){
  #     
  #   
  #     r_pred <- raster(lf_pred_list[i])
  #   
  #     res_pix <- 480
  #     col_mfrow <- 1
  #     row_mfrow <- 1
  #   
  #     png(filename=paste("Figure4_models_predicted_surfaces_",model_name[i],"_",name_method_var,"_",data_selected,"_",out_suffix,".png",sep=""),
  #        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  #   
  #     plot(r_pred)
  #     title(paste("Mosaiced",model_name[i],name_method_var,date_selected,sep=" "))
  #     dev.off()
  #   }
  #   #Plot Delta and clim...
  # 
  #    ## plotting of delta and clim for later scripts...
  # 
  # }
  
  
  ######################
  ### Figure 5: plot accuracy ranked 
  
  #Turn summary table to a point shp
  
  list_df_ac_mod <- vector("list",length=length(model_name))
  for (i in 1:length(model_name)){
    
    ac_mod <- summary_metrics_v[summary_metrics_v$pred_mod==model_name[i],]
    ### Ranking by tile...
    df_ac_mod <- arrange(as.data.frame(ac_mod),desc(rmse))[,c("pred_mod","rmse","mae","tile_id")]
    list_df_ac_mod[[i]] <- arrange(as.data.frame(ac_mod),desc(rmse))[,c("rmse","mae","tile_id")]
    
    res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
    fig_filename <- paste("Figure5_ac_metrics_ranked_",model_name[i],"_",out_suffix,".png",sep="")

    png(filename=paste("Figure5_ac_metrics_ranked_",model_name[i],"_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    x<- as.character(df_ac_mod$tile_id)
    barplot(df_ac_mod$rmse, names.arg=x)
    #plot(ac_mod1,cex=sqrt(ac_mod1$rmse),pch=1,add=T)
    #plot(ac_mod1,cex=(ac_mod1$rmse1)*2,pch=1,add=T)
    title(paste("RMSE ranked by tile for ",model_name[i],sep=""))
    
    dev.off()
    list_outfiles[[counter_fig+i]] <- fig_filename
  }
  
  counter_fig <- counter_fig + length(model_name)
  
  r10 <-c("figure_5","Barplot of accuracy metrics ranked by tiles","mod1",metric_name,region_name,year_predicted,list_outfiles[[8]])
  r11 <-c("figure_5","Barplot of accuracy metrics ranked by tiles","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[9]])  

  ######################
  ### Figure 6: plot map of average RMSE per tile at centroids
  
  ### Without 
  
  #list_df_ac_mod <- vector("list",length=length(lf_pred_list))
  list_df_ac_mod <- vector("list",length=length(model_name))
  
  for (i in 1:length(model_name)){
    
    ac_mod <- summary_metrics_v[summary_metrics_v$pred_mod==model_name[i],]
    #r_pred <- raster(lf_list[i])
    
    res_pix <- 1200
    #res_pix <- 480
    
    col_mfrow <- 1
    row_mfrow <- 1
    fig_filename <- paste("Figure6_ac_metrics_map_centroids_tile_",model_name[i],"_",out_suffix,".png",sep="")
    png(filename=paste("Figure6_ac_metrics_map_centroids_tile_",model_name[i],"_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
    coordinates(ac_mod) <- ac_mod[,c("lon","lat")] 
    #coordinates(ac_mod) <- ac_mod[,c("lon.x","lat.x")] #solve this later
    #p_shp <- layer(sp.polygons(reg_layer, lwd=1, col='black'))
    #p_shp <- spplot(reg_layer,border="black",col="transparent")
    p_shp <- spplot(reg_layer,"ISO3" ,col.regions=NA, col="black") #ok problem solved!!

    #title("(a) Mean for 1 January")
    p <- bubble(ac_mod,"rmse",main=paste("Average RMSE per tile and by ",model_name[i]))
    p1 <- p+p_shp
    print(p1)
    #plot(ac_mod1,cex=(ac_mod1$rmse1)*2,pch=1,add=T)
    #title(paste("Averrage RMSE per tile and by ",model_name[i]))
    
    dev.off()
    
    ### Ranking by tile...
    #df_ac_mod <- 
    list_df_ac_mod[[i]] <- arrange(as.data.frame(ac_mod),desc(rmse))[,c("rmse","mae","tile_id")]
    list_outfiles[[counter_fig+i]] <- fig_filename
  }
  counter_fig <- counter_fig+length(model_name)

  r12 <-c("figure_6","Average accuracy metrics map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[8]])
  r13 <-c("figure_6","Average accuracy metrics map at centroids","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[9]])  
  
  
  ######################
  ### Figure 7: Number of predictions: daily and monthly
  
  ## Figure 7a
 
  ## Number of tiles with information:
  sum(df_tile_processed$metrics_v) #26,number of tiles with raster object
  length(df_tile_processed$metrics_v) #26,number of tiles in the region
  sum(df_tile_processed$metrics_v)/length(df_tile_processed$metrics_v) #80 of tiles with info
  
  #coordinates
  try(coordinates(summary_metrics_v) <- c("lon","lat"))
  #try(coordinates(summary_metrics_v) <- c("lon.y","lat.y"))
  
  #threshold_missing_day <- c(367,365,300,200)
  
  nb<-nrow(subset(summary_metrics_v,model_name=="mod1"))  
  sum(subset(summary_metrics_v,model_name=="mod1")$n_missing)/nb #33/35
  
  ## Make this a figure...
  
  #plot(summary_metrics_v)
  #Make this a function later so that we can explore many thresholds...
  #Problem here
  #Browse[3]> c
   #Error in grid.Call.graphics(L_setviewport, pvp, TRUE) : 
  #non-finite location and/or size for viewport

  j<-1 #for model name 1,mod1
  for(i in 1:length(threshold_missing_day)){
    
    #summary_metrics_v$n_missing <- summary_metrics_v$n == 365
    #summary_metrics_v$n_missing <- summary_metrics_v$n < 365
    summary_metrics_v$n_missing <- as.numeric(summary_metrics_v$n < threshold_missing_day[i])
    summary_metrics_v_subset <- subset(summary_metrics_v,model_name=="mod1")
    
    fig_filename <- paste("Figure7a_ac_metrics_map_centroids_tile_",model_name[j],"_","missing_day_",threshold_missing_day[i],
                       "_",out_suffix,".png",sep="")

    if(sum(summary_metrics_v_subset$n_missing) > 0){#then there are centroids to plot!!!
      
      #res_pix <- 1200
      res_pix <- 960
      col_mfrow <- 1
      row_mfrow <- 1
      #only mod1 right now
      png(filename=paste("Figure7a_ac_metrics_map_centroids_tile_",model_name[j],"_","missing_day_",threshold_missing_day[i],
                       "_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
      model_name[j]
    
      #p_shp <- layer(sp.polygons(reg_layer, lwd=1, col='black'))
      p_shp <- spplot(reg_layer,"ISO3" ,col.regions=NA, col="black") #ok problem solved!!
      #title("(a) Mean for 1 January")
      p <- bubble(summary_metrics_v_subset,"n_missing",main=paste("Missing per tile and by ",model_name[j]," for ",
                                                                threshold_missing_day[i]))
      p1 <- p+p_shp
      try(print(p1)) #error raised if number of missing values below a threshold does not exist
      dev.off()

    } 
     
    list_outfiles[[counter_fig+i]] <- fig_filename
  }
  counter_fig <- counter_fig+length(threshold_missing_day) #currently 4 days...

  r14 <-c("figure_7","Number of missing days threshold1 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[8]])
  r15 <-c("figure_7","Number of missing days threshold2 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[9]])  
  r16 <-c("figure_7","Number of missing days threshold3 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[8]])
  r17 <-c("figure_7","Number of missing days threshold4 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[9]])  

  ### Potential
  png(filename=paste("Figure7b_number_daily_predictions_per_models","_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  xyplot(n~pred_mod | tile_id,data=subset(as.data.frame(summary_metrics_v),
                                          pred_mod!="mod_kr"),type="h")
  dev.off()
  
  list_outfiles[[counter_fig+1]] <- paste("Figure7b_number_daily_predictions_per_models","_",out_suffix,".png",sep="")
  counter_fig <- counter_fig + 1
  r18 <-c("figure_7b","Number of daily predictions per_models","mod1",metric_name,region_name,year_predicted,list_outfiles[[9]])  
  
  table(tb$pred_mod)
  table(tb$index_d)
  #table(subset(tb,pred_mod!="mod_kr"))
  table(subset(tb,pred_mod=="mod1")$index_d)
  #aggregate()
  tb$predicted <- 1
  test <- aggregate(predicted~pred_mod+tile_id,data=tb,sum)
  #xyplot(predicted~pred_mod | tile_id,data=subset(as.data.frame(test),
  #                                                pred_mod!="mod_kr"),type="h")
  
  as.character(unique(test$tile_id)) #141 tiles
  
  dim(subset(test,test$predicted==365 & test$pred_mod=="mod1"))
  #histogram(subset(test, test$pred_mod=="mod1")$predicted)
  unique(subset(test, test$pred_mod=="mod1")$predicted)
  table((subset(test, test$pred_mod=="mod1")$predicted))
  
  #LST_avgm_min <- aggregate(LST~month,data=data_month_all,min)
  png(filename=paste("Figure7c_histogram_number_daily_predictions_per_models","_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  histogram(test$predicted~test$tile_id)
  dev.off()
  
  list_outfiles[[counter_fig+1]] <- paste("Figure7c_histogram_number_daily_predictions_per_models","_",out_suffix,".png",sep="")
  counter_fig <- counter_fig + 1
  r19 <-c("figure_7c","Histogram number daily predictions per models","mod1",metric_name,region_name,year_predicted,list_outfiles[[9]])  

  
  #table(tb)
  ## Figure 7b
  #png(filename=paste("Figure7b_number_daily_predictions_per_models","_",out_suffix,".png",sep=""),
  #    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  #xyplot(n~month | tile_id + pred_mod,data=subset(as.data.frame(tb_month_s),
  #                                           pred_mod!="mod_kr"),type="h")
  #xyplot(n~month | tile_id,data=subset(as.data.frame(tb_month_s),
  #                                           pred_mod="mod_1"),type="h")
  #test=subset(as.data.frame(tb_month_s),pred_mod="mod_1")
  #table(tb_month_s$month)
  #dev.off()
  #
  
  ##########################################################
  ##### Figure 8: Breaking down accuracy by regions!! #####
  
  #summary_metrics_v <- merge(summary_metrics_v,df_tile_processed,by="tile_id")
  
  ##################
  ##First plot with all models together
  
  ## Figure 8a
  res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1
  
  fig_filename <- paste("Figure8a_boxplot_overall_separated_by_region_with_oultiers_",out_suffix,".png",sep="")
  png(filename=fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  p<- bwplot(rmse~pred_mod | reg, data=tb,
             main="RMSE per model and region over all tiles with outliers")
  print(p)
  dev.off()
  
  list_outfiles[[counter_fig+1]] <- fig_filename
  counter_fig <- counter_fig + 1
  
  ## Figure 8b
  fig_filename <- paste("Figure8b_boxplot_overall_separated_by_region_scaling_",out_suffix,".png",sep="")
  png(filename=fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  #boxplot(rmse~pred_mod,data=tb,ylim=c(0,5),outline=FALSE)#,names=tb$pred_mod)
  #title("RMSE per model over all tiles")
  p<- bwplot(rmse~pred_mod | reg, data=tb,ylim=c(0,5),
             main="RMSE per model and region over all tiles with scaling")
  print(p)
  dev.off()
  
  list_outfiles[[counter_fig+1]] <- fig_filename
  counter_fig <- counter_fig + 1
  
  
  r20 <-c("figure 8a","Boxplot overall accuracy by model separated by region with outliers",NA,metric_name,region_name,year_predicted,list_outfiles[[20]])  
  r21 <-c("figure 8b","Boxplot overall accuracy by model separated by region with scaling",NA,metric_name,region_name,year_predicted,list_outfiles[[21]])  

  #######
  ##Second, plot for each model separately
  
  for(i in 1:length(model_name)){
    
    tb_subset <- subset(tb,pred_mod==model_name[i])#mod1 is i=1, mod_kr is last
    ## Figure 8c
  
    res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
  
    fig_filename <- paste("Figure8c_boxplot_overall_accuracy_separated_by_region_with_outliers_",model_name[i],"_",out_suffix,".png",sep="")
    png(filename=fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
    p<- bwplot(rmse~pred_mod | reg, data=tb_subset,
             main="RMSE per model and region over all tiles with outliers")
    print(p)
    dev.off()
  
    list_outfiles[[counter_fig+1]] <- fig_filename
    counter_fig <- counter_fig + 1
  
    ## Figure 8d
    fig_filename <- paste("Figure8d_boxplot_overall_accuracy_separated_by_region_scaling_",model_name[i],"_",out_suffix,".png",sep="")
    png(filename=fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
    #boxplot(rmse~pred_mod,data=tb,ylim=c(0,5),outline=FALSE)#,names=tb$pred_mod)
    #title("RMSE per model over all tiles")
    p<- bwplot(rmse~pred_mod | reg, data=tb_subset,ylim=c(0,5),
             main="RMSE per model and region over all tiles with scaling")
    print(p)
    dev.off()
  
    list_outfiles[[counter_fig+1]] <- fig_filename
    counter_fig <- counter_fig + 1

  }
  
  r22 <-c("figure 8c","Boxplot overall accuracy separated by region with outliers","mod1",metric_name,region_name,year_predicted,list_outfiles[[22]])  
  r23 <-c("figure 8d","Boxplot overall accuracy separated by region with scaling","mod1",metric_name,region_name,year_predicted,list_outfiles[[23]])  
  r24 <-c("figure 8c","Boxplot overall accuracy separated by region with outliers","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[24]])  
  r25 <-c("figure 8d","Boxplot overall accuracy separated by region with scaling","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[25]])  

  #####################################################
  #### Figure 9: plotting boxplot by year and regions ###########
  
  ## Figure 9a
  res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1

  png(filename=paste("Figure9a_boxplot_overall_separated_by_year_and_model_with_oultiers_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  #This will need to be changed, the layout is difficult at this stage 
  #p<- bwplot(rmse~pred_mod + reg + year_predicted, data=tb,
  #           main="RMSE per model and region over all tiles")
  p<- bwplot(rmse~pred_mod | year_predicted, data=tb,
             main="RMSE per model and region over all tiles")
  print(p)
  dev.off()
  
  ## Figure 9b
  png(filename=paste("Figure9b_boxplot_overall_separated_by_year_and_model_scaling_",out_suffix,".png",sep=""),
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  #boxplot(rmse~pred_mod,data=tb,ylim=c(0,5),outline=FALSE)#,names=tb$pred_mod)
  #title("RMSE per model over all tiles")
  p<- bwplot(rmse~pred_mod | year_predicted, data=tb,ylim=c(0,5),
             main="RMSE per model and region over all tiles")
  print(p)
  dev.off()

  for(i in 1:length(model_name)){
    
    tb_subset <- subset(tb,pred_mod==model_name[i])#mod1 is i=1, mod_kr is last
    ## Figure 9c
    fig_filename <- paste("Figure9c_boxplot_overall_accuracy_separated_by_year_with_outliers_",model_name[i],"_",out_suffix,".png",sep="")
    png(filename=fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
    boxplot(rmse~year_predicted,data=tb_subset,ylab=metric_name,xlab="year predicted")#,names=tb$pred_mod)
    title(paste("Overall accuracy for ", model_name[i], " by year for ",region_name,sep=""))
    
    #p<- bwplot(rmse~year_predicted | reg , data=tb_subset,ylim=c(0,5),
             #main="RMSE per model and region over all tiles")
    #print(p)
    dev.off()
    
    list_outfiles[[counter_fig+1]] <- fig_filename
    counter_fig <- counter_fig + 1

    fig_filename <- paste("Figure9d_boxplot_overall_separated_by_year_scaling_",model_name[i],"_",out_suffix,".png",sep="")
    png(filename=fig_filename,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
    boxplot(rmse~year_predicted,data=tb_subset,ylim=c(0,5),outline=FALSE,ylab=metric_name,xlab="year predicted")#,names=tb$pred_mod)
    title(paste("Overall accuracy for ", model_name[i], " by year for ",region_name,sep=""))
    #p<- bwplot(rmse~year_predicted | reg , data=tb_subset,ylim=c(0,5),
             #main="RMSE per model and region over all tiles")
    #print(p)
    dev.off()

    list_outfiles[[counter_fig+1]] <- fig_filename
    counter_fig <- counter_fig + 1
  }
  
  r26 <-c("figure 9a","Boxplot overall accuracy separated_by year and model with oultiers",NA,metric_name,region_name,year_predicted,list_outfiles[[22]])  
  r27 <-c("figure 9b","Boxplot overall accuracy separated_by year and model with scaling",NA,metric_name,region_name,year_predicted,list_outfiles[[23]])  
  r28 <-c("figure 9c","Boxplot overall accuracy separated by year with outliers","mod1",metric_name,region_name,year_predicted,list_outfiles[[24]])  
  r29 <-c("figure 9d","Boxplot overall accuracy separated by year with scaling","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[25]])  
  r30 <-c("figure 9c","Boxplot overall accuracy separated by year with outliers","mod1",metric_name,region_name,year_predicted,list_outfiles[[24]])  
  r31 <-c("figure 9d","Boxplot overall accuracy separated by year with scaling","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[25]])  

  ##############################################################
  ############## Prepare object to return
  ############## Collect information from assessment ##########
  
  # This is hard coded and can be improved later on for flexibility. It works for now...                                                                 
  #This data.frame contains all the files from the assessment

  #Should have this at the location of the figures!!! will be done later?    
  #r1 <-c("figure_1","Tiles processed for the region",NA,NA,region_name,year_predicted,list_outfiles[[1]])
  #r2 <-c("figure_2a","Boxplot of accuracy with outliers by tiles","mod1",metric_name,region_name,year_predicted,list_outfiles[[2]]) 
  #r3 <-c("figure_2a","boxplot of accuracy with outliers by tiles","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[3]])
  #r4 <-c("figure_2b","Boxplot of accuracy with scaling by tiles","mod1",metric_name,region_name,year_predicted,list_outfiles[[4]])  
  #r5 <-c("figure_2b","Boxplot of accuracy with scaling by tiles","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[5]])  
  #r6 <-c("figure_3a","Boxplot overall accuracy with outliers","mod1",metric_name,region_name,year_predicted,list_outfiles[[6]])  
  #r7 <-c("figure_3b","Boxplot overall accuracy with scaling","mod1",metric_name,region_name,year_predicted,list_outfiles[[7]])  
  #r8 <-c("figure_3a","Boxplot overall accuracy with outliers","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[8]])
  #r9 <-c("figure_3b","Boxplot overall accuracy with scaling","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[9]])  
  #r10 <-c("figure_5","Barplot of accuracy metrics ranked by tiles","mod1",metric_name,region_name,year_predicted,list_outfiles[[10]])
  #r11 <-c("figure_5","Barplot of accuracy metrics ranked by tiles","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[11]])  
  #r12 <-c("figure_6","Average accuracy metrics map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[12]])
  #r13 <-c("figure_6","Average accuracy metrics map at centroids","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[13]])  
  #r14 <-c("figure_7","Number of missing days threshold1 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[14]])
  #r15 <-c("figure_7","Number of missing days threshold2 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[15]])  
  #r16 <-c("figure_7","Number of missing days threshold3 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[16]])
  #r17 <-c("figure_7","Number of missing days threshold4 map at centroids","mod1",metric_name,region_name,year_predicted,list_outfiles[[17]])  
  #r18 <-c("figure_7b","Number of daily predictions per_models","mod1",metric_name,region_name,year_predicted,list_outfiles[[18]])  
  #r19 <-c("figure_7c","Histogram number daily predictions per models","mod1",metric_name,region_name,year_predicted,list_outfiles[[19]])  
  #r20 <-c("figure 8a","Boxplot overall accuracy by model separated by region with outliers",NA,metric_name,region_name,year_predicted,list_outfiles[[20]])  
  #r21 <-c("figure 8b","Boxplot overall accuracy by model separated by region with scaling",NA,metric_name,region_name,year_predicted,list_outfiles[[21]])  
  #r22 <-c("figure 8c","Boxplot overall accuracy separated by region with outliers","mod1",metric_name,region_name,year_predicted,list_outfiles[[22]])  
  #r23 <-c("figure 8d","Boxplot overall accuracy separated by region with scaling","mod1",metric_name,region_name,year_predicted,list_outfiles[[23]])  
  #r24 <-c("figure 8c","Boxplot overall accuracy separated by region with outliers","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[24]])  
  #r25 <-c("figure 8d","Boxplot overall accuracy separated by region with scaling","mod_kr",metric_name,region_name,year_predicted,list_outfiles[[25]])  

  #Assemble all the figures description and information in a data.frame for later use
  list_rows <-list(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,
                   r25,r26,r27,r28,r29,r30,r31)
  df_assessment_figures_files <- as.data.frame(do.call(rbind,list_rows))
  names(df_assessment_figures_files) <- c("figure_no","comment","model_name","reg","metric_name","year_predicted","filename")
  
  ###Prepare files for copying back?
  df_assessment_figures_files_names <- file.path(out_dir,paste("df_assessment_figures_files_",region_name,"_",year_predicted,"_",out_suffix,".txt",sep=""))
  write.table(df_assessment_figures_files,
              file=df_assessment_figures_files_names ,sep=",")

  #df_assessment_figures_files_names
  
  ######################################################
  ##### Prepare objet to return ####

  assessment_obj <- list(df_assessment_files, df_assessment_figures_files)
  names(assessment_obj) <- c("df_assessment_files", "df_assessment_figures_files")
  ## Prepare list of files to return...
  return(assessment_obj)
 
}
  
##################### END OF SCRIPT ######################
