####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up and mosaicing on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Different options to explore mosaicing are tested.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 04/14/2015  
#MODIFIED ON: 02/13/2017            
#Version: 6
#PROJECT: Environmental Layers project     
#COMMENTS: analyses run for reg4 1991 for test of mosaicing using 1500x4500km and other tiles
#TODO:
#1) Make this is a script/function callable from the shell/bash
#2) clean up temporary files, it builds currently on the disk
#3) fix output folder for some of output files: create a mosaic output folder if doesn't exist?
#4) create a helper function for inputs/arguments to automate (optparse pacakge)...?? 
    #Could also be in the assessment stage

### Before running, the gdal modules and other environment parameters need to be set if on NEX-NASA.
### This can be done by running the following commands:
#
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#MODULEPATH=$MODULEPATH:/nex/modules/files
#module load pythonkits/gdal_1.10.0_python_2.7.3_nex
#
#setfacl -Rmd user:aguzman4:rwx /nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015
#
#reg1   : North America
#reg23  : Europe + Asia
#reg4   : South America
#reg5   : Africa
#reg6   : Oceania+ South East Asia
#

### COMMIT: moving function to list mosaics by regions and additional changes

#################################################################################################


#### FUNCTION USED IN SCRIPT


############################################
#### Parameters and constants  


########################## START SCRIPT ##############################

run_mosaicing_prediction_fun <-function(i,list_param_run_mosaicing_prediction){

  ##This is a general function to mosaic predicted tiles and accuracy layers with 34 input parameters.
  #
  #
  #Input Parameters:
  #1) in_dir: input directory, parent directory containing predictions for all regions #PARAM1
  #2) y_var_name: climate variable predicted (dailyTmax, dailyTmin, dailyPrecip) # #PARAM2
  #3) interpolation_method: interpolation method being used, gam_CAI currently ##PARAM3
  #4) region_name: region_name e.g. reg4 South America #PARAM4
  #5) mosaicing_method: options include "unweighted","use_edge_weights" #PARAM5
  #6) out_suffix : output suffix #PARAM 6
  #7) out_suffix_str: additional output suffix with region name #PARAM 7
  #8) metric_name: metric or columns to use for additional mosaicing: "rmse" #RMSE, MAE etc. #PARAM 8
  #9) pred_mod_name : model name used e.g. "mod1" #PARAM 9
  #10) var_pred : variable for use in residuals mapping (e.g. "res_mod1") #PARAM 10
  #11) out_dir: output directory #PARAM 11
  #12) create_out_dir_param: if TRUE then create a new dir #PARAM 12
  #13) day_to_mosaic_range: start and end date for daily mosaics, if NULL then mosaic all days of the year #PARAM 12
  #14) year_predicted: year of the prediction being mosaiced (process is done by year)
  #15) proj_str :porjection used by tiles e.g. CRS_WGS84 #PARAM 13
  #16) file_format: output file format used for raster eg ".tif" #PARAM 14
  #17) NA_value: NA value used e.g. -9999 #PARAM 15
  #18) num_cores: number of cores #PARAM 16                 
  #19) use_autokrige: use_autokrige if FALSE use kriging from Fields package #PARAM 18
  #20) infile_mask: input file mask used for the region under process #PARAM 19
  #21) tb_accuracy_name: daily accuracy from testing/validation stations by tiles #PARAM 20
  #22) data_month_s_name: training stations for climatology time steps  #PARAM 21
  #23) data_day_v_name:  testing stations for daily predictions combined #PARAM 22
  #24) data_day_s_name: training stations for daily predictions cominbed #PARAM 23
  #25) df_tile_processed_name: processed tiles from the accuracy assessment ##PARAM 24
  #26) mosaic_python: python script used in the mosoicing (gdalmerge script from Alberto Guzmann) #PARAM 25
  #27) python_bin: directory for general python "/usr/bin" #PARAM 26
  #28) algorithm: python or R, if R use mosaic function for R, if python use modified gdal merge, PARAM 27
  #29) match_extent : if "FALSE" try without matching geographic extent #PARAM 28 
  #30) list_models : if NULL use y~1 formula #PARAM 29
  #31) layers_option: mosaic to create as a layer from var_pred (e.g. TMax), res_training, res_testing, ac_testing
  #32) tmp_files: if TRUE keep temporary files generated during mosaicing
  #33) data_type: if NULL, use Float32, other possibilities are gdal based in the final output
  #34) scaling: scaling factor to multiply the original variable before conversation to int
  #35) values_range: valid range for predicted values and mosaic e.g. -100,100
  #36) infile_reg_mosaics: input mosaic files
  
  ###OUTPUT
  # 
  #
  
  ###Loading R library and packages     
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
  #library(spgwr)                               # GWR method
  library(automap)                             # Kriging automatic fitting of variogram using gstat
  library(rgeos)                               # Geometric, topologic library of functions
  #RPostgreSQL                                 # Interface R and Postgres, not used in this script
  library(gridExtra)
  #Additional libraries not used in workflow
  library(pgirmess)                            # Krusall Wallis test with mulitple options, Kruskalmc {pgirmess}  
  library(colorRamps)
  library(zoo)
  library(xts)
  
  #Data is on ATLAS or NASA NEX
  #PARAM 1
  in_dir <- list_param_run_mosaicing_prediction$in_dir #PARAM1
  y_var_name <- list_param_run_mosaicing_prediction$y_var_name # #PARAM2
  interpolation_method <- list_param_run_mosaicing_prediction$interpolation_method ##PARAM3
  region_name <- list_param_run_mosaicing_prediction$region_name #"reg4" #PARAM 4 #reg4 South America, Africa reg5,Europe reg2, North America reg1, Asia reg3
  mosaicing_method <- list_param_run_mosaicing_prediction$mosaicing_method # c("unweighted","use_edge_weights") #PARAM5
  out_suffix <- list_param_run_mosaicing_prediction$out_suffix #paste(region_name,"_","run10_1500x4500_global_analyses_pred_1992_12072015",sep="") #PARAM 6
  out_suffix_str <- list_param_run_mosaicing_prediction$out_suffix_str #"run10_1500x4500_global_analyses_pred_1992_12072015" #PARAM 7
  metric_name <- list_param_run_mosaicing_prediction$metric_name # "rmse" #RMSE, MAE etc. #PARAM 8
  pred_mod_name <- list_param_run_mosaicing_prediction$pred_mod_name #"mod1" #PARAM 9
  var_pred <- list_param_run_mosaicing_prediction$var_pred # "res_mod1" #used in residuals mapping #PARAM 10
  out_dir <- list_param_run_mosaicing_prediction$out_dir #PARAM 11
  create_out_dir_param <- list_param_run_mosaicing_prediction$create_out_dir_param # FALSE #PARAM 12
  
  #if daily mosaics NULL then mosaicas all days of the year #PARAM 13
  day_to_mosaic_range <- list_param_run_mosaicing_prediction$day_to_mosaic_range # c("19920101","19920102","19920103") #,"19920104","19920105") #PARAM14, two dates note in /tiles for now on NEX
  year_processed <- list_param_run_mosaicing_prediction$year_predicted #PARAM 15
  #CRS_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84 #CONSTANT1
  #CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
  proj_str <- list_param_run_mosaicing_prediction$proj_str# CRS_WGS84 #PARAM 16 #check this parameter
  
  file_format <- list_param_run_mosaicing_prediction$file_format # ".tif" #PARAM 17
  NA_value <- list_param_run_mosaicing_prediction$NA_value # -9999 #PARAM 18
  
  num_cores <- list_param_run_mosaicing_prediction$num_cores  #6 #PARAM 19                 
  #region_names <- list_param_run_mosaicing_prediction$region_names # c("reg23","reg4") #selected region names, ##PARAM 18 
  use_autokrige <- list_param_run_mosaicing_prediction$use_autokrige # F #PARAM 20
  
  ###Separate folder for masks by regions, should be listed as just the dir!!... #PARAM 21
  #infile_mask <- "/nobackupp8/bparmen1/regions_input_files/r_mask_reg4.tif"
  infile_mask <- list_param_run_mosaicing_prediction$infile_mask # input mask used in defining the region
  
  #in_dir can be on NEX or Atlas
  
  ##skip this for now
  df_assessment_files_name <- list_param_run_mosaicing_prediction$df_assessment_files_name # data.frame with all files used in assessmnet, PARAM 21

  #python script and gdal on NEX NASA:
  #mosaic_python <- "/nobackupp6/aguzman4/climateLayers/sharedCode/"
  #python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin"
  #python script and gdal on Atlas NCEAS
  mosaic_python <- list_param_run_mosaicing_prediction$mosaic_python # "/data/project/layers/commons/NEX_data/sharedCode" #PARAM 22
  python_bin <- list_param_run_mosaicing_prediction$python_bin # "/usr/bin" #PARAM 23
  
  algorithm <- list_param_run_mosaicing_prediction$algorithm #"python" #PARAM 24 #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
  #algorithm <- "R" #if R use mosaic function for R, if python use modified gdalmerge script from Alberto Guzmann
  
  match_extent <- list_param_run_mosaicing_prediction$match_extent #"FALSE" #PARAM 25 #try without matching!!!
  
  #for residuals...
  list_models <- list_param_run_mosaicing_prediction$list_models #  NULL #PARAM 26
  #list_models <- paste(var_pred,"~","1",sep=" ") #if null then this is the default...
  layers_option <- list_param_run_mosaicing_prediction$layers_option #PARAM 27
  tmp_files <- list_param_run_mosaicing_prediction$tmp_files  #PARAM 28
  data_type <- list_param_run_mosaicing_prediction$data_type #PARAM 29
  scaling <- list_param_run_mosaicing_prediction$scaling 
  values_range <- list_param_run_mosaicing_prediction$values_range
  infile_reg_mosaics <- list_param_run_mosaicing_prediction$infile_reg_mosaics
  
  #################################################################
  ####### PART 1: Read in data and process data ########
  ########################################################
  
  #browser()
  #out_dir <- in_dir #PARAM 11
  #in_dir_tiles <- file.path(in_dir,"tiles") #this is valid both for Atlas and NEX
  NA_flag_val <- NA_value #PARAM 16
  
  #in_dir <- file.path(in_dir,region_name)
  #out_dir <- in_dir
  if(create_out_dir_param==TRUE){
    out_dir_tmp <- file.path(out_dir,"mosaic")
    #   #create if does not exists
    if(!file.exists(out_dir_tmp)){
      dir.create(out_dir_tmp)
    }
    out_dir <- create_dir_fun(out_dir_tmp,out_suffix)
    setwd(out_dir)
  }else{
    setwd(out_dir) #use previoulsy defined directory
  }
  
  setwd(out_dir)
  
  #browser()
  
  ### Read in assessment and accuracy files if not null (not world mosaic)
  if(!is.null(df_assessment_files_name)){

    df_assessment_files <- read.table(df_assessment_files_name,stringsAsFactors=F,sep=",")
    #browser()
    tb_v_accuracy_name <- df_assessment_files$files[2] 
    tb_s_accuracy_name <- df_assessment_files$files[4] 
    tb_s_month_accuracy_name <- df_assessment_files$files[3] 
    data_month_s_name <- df_assessment_files$files[5] 
    data_day_s_name <- df_assessment_files$files[6] 
    data_day_v_name <- df_assessment_files$files[7] 

    ##data_month_v_name <- file.path(in_dir,basename(df_assessment_files$files[8])) 
    pred_data_month_info_name <- df_assessment_files$files[10]
    pred_data_day_info_name <- df_assessment_files$files[11]
    df_tile_processed_name <- df_assessment_files$files[12]
    # accuracy table by tiles
    tb <- read.table(tb_v_accuracy_name,sep=",")
    tb_s <- read.table(tb_s_accuracy_name,sep=",")
    data_month_s <- read.table(data_month_s_name,sep=",") # textfiles of stations by month
    data_day_s <- read.table(data_day_s_name,sep=",") #daily testing/validation stations by dates and tiles
    data_day_v <- read.table(data_day_v_name,sep=",") #daily training stations by dates and tiles
    df_tile_processed <- read.table(df_tile_processed_name,sep=",")
    
  }

  ##Read additional data from table assessment, add later
  #pred_data_day_info_1999_reg4_1999.txt
  #pred_data_day_info_name <- df_assessment_files$files[11]
  
  #this part needs to be improve make this a function and use multicore to loop through files...
  #give a range of dates to run...
  #browser()
  
  if(is.null(day_to_mosaic_range)){
  #  start_date <- #first date
     start_date <- paste0(year_processed,"0101") #change this later!!
     end_date <-   paste0(year_processed,"0101") #change this later!!
     day_to_mosaic <- seq(as.Date(strptime(start_date,"%Y%m%d")), as.Date(strptime(end_date,"%Y%m%d")), 'day')
     day_to_mosaic <- format(day_to_mosaic,"%Y%m%d") #format back to the relevant date format for files
  }else{
    start_date <- day_to_mosaic_range[1]
    end_date <- day_to_mosaic_range[2]
    day_to_mosaic <- seq(as.Date(strptime(start_date,"%Y%m%d")), as.Date(strptime(end_date,"%Y%m%d")), 'day')
    day_to_mosaic <- format(day_to_mosaic,"%Y%m%d") #format back to the relevant date format for files
  }
  
  if(is.null(infile_reg_mosaics)){
    
    in_dir_tiles_tmp <- file.path(in_dir, region_name)
    lf_mosaic <- mclapply(1:length(day_to_mosaic),FUN=function(i){
    searchStr = paste(in_dir_tiles_tmp,"/*/",year_processed,"/gam_CAI_dailyTmax_predicted_",pred_mod_name,"*",day_to_mosaic[i],"*.tif",sep="")
    Sys.glob(searchStr)},mc.preschedule=FALSE,mc.cores = num_cores)
    
  }else{
    
    tb_reg_mosaic_input <- read.table(infile_reg_mosaics,sep=",")
    in_dir_mosaic_reg_list <- as.character(tb_reg_mosaic_input[,1])
    in_dir_tiles_tmp <- in_dir_mosaic_reg_list
    ### need to modify this

    #debug(get_mosaic_files_fun)
    #test_df <- get_mosaic_files_fun(1,day_to_mosaic_range,in_dir_tiles_tmp,year_processed)
    
    lf_mosaic <- mclapply(1:length(day_to_mosaic),
                          FUN=get_mosaic_files_fun,
                          day_to_mosaic=day_to_mosaic,
                          in_dir_tiles_tmp=in_dir_tiles_tmp,
                          year_processed=year_processed,
                          mc.preschedule=FALSE,
                          mc.cores = num_cores)
  }
  #browser()
  
  #########################################################################
  ##################### PART 2: produce the mosaic ##################
  ######################################################################
  
  #This is is assuming a list of file for a region!! 
  #this is where the main function for mosaicing region starts!!
  #use reg4 to test the code for now, redo later for any regions!!!

  region_selected <- region_name

  #There a 28 files for reg4, South America
    
  ######################################################
  #### PART 3: GENERATE MOSAIC FOR LIST OF FILES #####
  #################################
  #### Mosaic tiles for the variable predicted and accuracy metrics, residuals surfaces or other options
    
  browser()

  #methods availbable:use_sine_weights,use_edge,use_linear_weights
  #only use edge method for now
  #loop to dates..., make this a function...
  #This is a loop but uses multicores when calling the mosaic function
  list_mosaic_obj <- vector("list",length=length(day_to_mosaic))
  
  for(i in 1:length(day_to_mosaic)){
    #
    
    if(layers_option=="var_pred"){
      
      #mosaic_method <- "use_edge_weights" #this is distance from edge
      mosaic_method <- mosaicing_method
      out_suffix_tmp <- paste(interpolation_method,y_var_name,day_to_mosaic[i],out_suffix,sep="_")
      #debug(mosaicFiles)
      #can also loop through methods!!!
      #python_bin <- "/usr/bin/" #python gdal bin, on Atlas NCEAS
      #python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules/bin" #on NEX
      #gdal_merge_sum_noDataTest.py
      
      mosaic_obj <- mosaicFiles(lf_mosaic[[i]],
                                mosaic_method="use_edge_weights",
                                num_cores=num_cores,
                                r_mask_raster_name=infile_mask,
                                python_bin=python_bin,
                                mosaic_python=mosaic_python,
                                algorithm=algorithm,
                                match_extent=match_extent,
                                df_points=NULL,
                                NA_flag=NA_flag_val,
                                file_format=file_format,
                                out_suffix=out_suffix_tmp,
                                out_dir=out_dir,
                                tmp_files=tmp_files,
                                data_type=data_type,
                                scaling=scaling,
                                values_range=values_range)
      #runs in 15-16 minutes for 3 dates and mosaicing of 28 tiles...
      list_mosaic_obj[[i]] <- mosaic_obj
    }
    #browser()
    
    if(layers_option=="ac_testing"){
    
      #try running the function now:
      lf <- lf_mosaic[i]
      df <- tb #for ac_testing
      days_to_process <- day_to_mosaic[i]
  
      #browser()

      lf_accuracy_testing_raster <- generate_ac_assessment_layers_by_tile(lf,layers_option,df,df_tile_processed,metric_name,
                                                    var_pred,list_models,use_autokrige,pred_mod_name,
                                                    y_var_name,interpolation_method,region_selected,
                                                    days_to_process,num_cores,NA_flag_val,file_format,
                                                    out_dir,out_suffix)   #### create a function to generate accuracy layers by tiles
      #if (inherits(mod,"try-error"))
      ## Now accuracy based on center of centroids
      mosaic_method <- "use_edge_weights" #this is distance from edge
      #Adding metric name in the name...
      #out_suffix_tmp <- paste(interpolation_method,metric_name,day_to_mosaic[i],out_suffix,sep="_")
      out_suffix_tmp <- paste(interpolation_method,metric_name,layers_option,day_to_mosaic[i],out_suffix,sep="_")
      
      #Can modify here....add creation of data for the specific date here rather than beforehand!!!!
      
      ## To be inserted...
      
      #undebug(mosaicFiles)
      #can also loop through methods!!!
      mosaic_obj <- mosaicFiles(unlist(lf_accuracy_testing_raster),
                                mosaic_method="use_edge_weights",
                                num_cores=num_cores,
                                r_mask_raster_name=infile_mask,
                                python_bin=python_bin,
                                mosaic_python=mosaic_python,
                                algorithm=algorithm,
                                match_extent = match_extent,
                                df_points=NULL,
                                NA_flag=NA_flag_val,
                                file_format=file_format,
                                out_suffix=out_suffix_tmp,
                                out_dir=out_dir,
                                tmp_files=tmp_files,
                                data_type=data_type,
                                scaling=scaling,
                                values_range=values_range)
      ##Took 12-13 minutes for 28 tiles and one date...!!! 
      #browser()
      
      if(tmp_files==F){ #if false...delete all files with "_tmp"
        lf_tmp <- unlist(lf_accuracy_testing_raster)
        ##now delete temporary files...
        file.remove(lf_tmp)
      }
      list_mosaic_obj[[i]] <- mosaic_obj
    }
    
    if(layers_option=="ac_training"){
      ## Now accuracy based on center of centroids
      
      #try running the function now:
      lf <- lf_mosaic[i]
      df <- tb_s #for ac_testing
      days_to_process <- day_to_mosaic[i]
  
      #browser()
      #debug(generate_ac_assessment_layers_by_tile)
      lf_accuracy_training_raster<- generate_ac_assessment_layers_by_tile(lf,layers_option,df,df_tile_processed,metric_name,
                                                    var_pred,list_models,use_autokrige,pred_mod_name,
                                                    y_var_name,interpolation_method,region_selected,
                                                    days_to_process,num_cores,NA_flag_val,file_format,
                                                    out_dir,out_suffix)   #### create a function to generate accuracy layers by tiles
      
      mosaic_method <- "use_edge_weights" #this is distance from edge
      #Adding metric name in the name...
      out_suffix_tmp <- paste(interpolation_method,metric_name,layers_option,day_to_mosaic[i],out_suffix,sep="_")
      #out_suffix_tmp <- paste(interpolation_method,"kriged_residuals","data_day_v",day_to_mosaic[i],out_suffix,sep="_")

      #Can modify here....add creation of data for the specific date here rather than beforehand!!!!
      
      ## To be inserted...
      
      #undebug(mosaicFiles)
      #can also loop through methods!!!
      #browser()
      mosaic_obj <- mosaicFiles(unlist(lf_accuracy_training_raster),
                                mosaic_method="use_edge_weights",
                                num_cores=num_cores,
                                r_mask_raster_name=infile_mask,
                                python_bin=python_bin,
                                mosaic_python=mosaic_python,
                                algorithm=algorithm,
                                match_extent = match_extent,
                                df_points=NULL,
                                NA_flag=NA_flag_val,
                                file_format=file_format,
                                out_suffix=out_suffix_tmp,
                                out_dir=out_dir,
                                tmp_files=tmp_files,
                                data_type=data_type,
                                scaling=scaling,
                                values_range=values_range)
      ##Took 13-14 minutes for 28 tiles and one date...!!! 
      if(tmp_files==F){ #if false...delete all files with "_tmp"
        lf_tmp <- unlist(lf_accuracy_training_raster)
        ##now delete temporary files...
        file.remove(lf_tmp)
      }
      list_mosaic_obj[[i]] <- mosaic_obj
    }

    #list_mosaic_obj[[i]] <- list(prediction=mosaic_edge_obj_prediction,accuracy=mosaic_edge_obj_accuracy)

    ### produce residuals mosaics
    if(layers_option=="res_testing"){
      
      lf <- lf_mosaic[i]
      df <- data_day_v #for ac_testing
      days_to_process <- day_to_mosaic[i]
  
      #browser()
      #print("browser")
      #debug(generate_ac_assessment_layers_by_tile)
      lf_accuracy_residuals_testing_raster <- generate_ac_assessment_layers_by_tile(lf,layers_option,df,df_tile_processed,metric_name,
                                                    var_pred,list_models,use_autokrige,pred_mod_name,
                                                    y_var_name,interpolation_method,region_selected,
                                                    days_to_process,num_cores,NA_flag_val,file_format,
                                                    out_dir,out_suffix)   #### create a function to generate accuracy layers by tiles

      #for now add data_day_s in the name!!
      mosaic_method <- "use_edge_weights" #this is distance from edge
      out_suffix_tmp <- paste(interpolation_method,"kriged_",layers_option,day_to_mosaic[i],out_suffix,sep="_")
      #lf_tmp<-list.files(pattern="*kriged_residuals.*.tif",full.names=T)
      lf_tmp <- unlist(lf_accuracy_residuals_testing_raster)
      #lf_accuracy_residuals_raster[[i]]
      #debug(mosaicFiles)
      mosaic_obj <- mosaicFiles(lf_tmp,
                                mosaic_method="use_edge_weights",
                                num_cores=num_cores,
                                r_mask_raster_name=infile_mask,
                                python_bin=python_bin,
                                mosaic_python=mosaic_python,
                                algorithm=algorithm,
                                match_extent=match_extent,
                                df_points=NULL,
                                NA_flag=NA_flag_val,
                                file_format=file_format,
                                out_suffix=out_suffix_tmp,
                                out_dir=out_dir,
                                tmp_files=tmp_files,
                                data_type=data_type,
                                scaling=scaling,
                                values_range=values_range)
      #Took 11 to 19 minues for one day and 28 tiles in region 4
      if(tmp_files==F){ #if false...delete all files with "_tmp"
        #lf_tmp <- unlist(lf_accuracy_residuals_testing_raster)
        ##now delete temporary files...
        file.remove(lf_tmp)
      }
      list_mosaic_obj[[i]] <- mosaic_obj
    }      
    
    ### produce residuals mosaics
    if(layers_option=="res_training"){
      
      lf <- lf_mosaic[i]
      df <- data_day_s #for ac_training
      days_to_process <- day_to_mosaic[i]
  
      #browser()
      #debug(generate_ac_assessment_layers_by_tile)
      lf_accuracy_residuals_training_raster <- generate_ac_assessment_layers_by_tile(lf,layers_option,df,df_tile_processed,metric_name,
                                                    var_pred,list_models,use_autokrige,pred_mod_name,
                                                    y_var_name,interpolation_method,region_selected,
                                                    days_to_process,num_cores,NA_flag_val,file_format,
                                                    out_dir,out_suffix)   #### create a function to generate accuracy layers by tiles

      #for now add data_day_s in the name!!
      mosaic_method <- "use_edge_weights" #this is distance from edge
      out_suffix_tmp <- paste(interpolation_method,"kriged_",layers_option,day_to_mosaic[i],out_suffix,sep="_")
      #lf_tmp<-list.files(pattern="*kriged_residuals.*.tif",full.names=T)
      lf_tmp <- unlist(lf_accuracy_residuals_training_raster)
      #lf_accuracy_residuals_raster[[i]]
      #debug(mosaicFiles)
      mosaic_obj <- mosaicFiles(lf_tmp,
                                mosaic_method="use_edge_weights",
                                num_cores=num_cores,
                                r_mask_raster_name=infile_mask,
                                python_bin=python_bin,
                                mosaic_python=mosaic_python,
                                algorithm=algorithm,
                                match_extent=match_extent,
                                df_points=NULL,
                                NA_flag=NA_flag_val,
                                file_format=file_format,
                                out_suffix=out_suffix_tmp,
                                out_dir=out_dir,
                                tmp_files=tmp_files,
                                data_type=data_type,
                                scaling=scaling,
                                values_range=values_range)
      list_mosaic_obj[[i]] <- mosaic_obj
      
      #Took 11 to 12 minues for one day and 28 tiles in region 4
      if(tmp_files==F){ #if false...delete all files with "_tmp"
        lf_tmp <- unlist(lf_accuracy_residuals_training_raster)
        ##now delete temporary files...
        file.remove(lf_tmp)
      }

    }
    
    ##End of mosaicing function for region predictions
  }## end of day_to_mosaic loop
  
  ##Create return object
  #list of mosaiced files: get the list of files now to include in the output object!!
  mosaicing_prediction_obj <- list(list_mosaic_obj,layers_option) #debugged
  names(mosaicing_prediction_obj) <- c("list_mosaic_obj","layers_option")
  fname_mosaicing_prediction_obj <- file.path(out_dir,paste("mosaicing_prediction_obj_",out_suffix_str,".RData",sep=""))
  save(mosaicing_prediction_obj,file= fname_mosaicing_prediction_obj)

  return(mosaicing_prediction_obj)
}

###############

##################### END OF SCRIPT ######################

