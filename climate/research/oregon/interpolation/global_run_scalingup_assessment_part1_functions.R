####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for assessment of scaling up on NEX: part 1 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#The purpose is to create as set of functions to diagnose and assess quickly a set of predictd tiles.
#Part 1 create summary tables and inputs for figure in part 2 and part 3.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 03/23/2014  
#MODIFIED ON: 02/05/2015            
#Version: 4
#PROJECT: Environmental Layers project  
#TO DO:
# - generate delta and clim mosaic
# - clean up

#First source file:
#source /nobackupp4/aguzman4/climateLayers/sharedModules/etc/environ.sh
#MODULEPATH=$MODULEPATH:/nex/modules/files
#module load /nex/modules/files/pythonkits/gdal_1.10.0_python_2.7.3_nex
# These are the names and number for the current subset regions used for global runs:
#reg1 - North America (NAM)
#reg2 - Western Europe (WE)
#reg3 - Eastern Europe to East Asia (EE_EA)
#reg4 - South America (SAM)
#reg5 - Africa (AF)
#reg6 - South East Asia and Australia (SEA_AUS)

#################################################################################################

### Loading R library and packages        
#library used in the workflow production:
library(gtools)                              # loading some useful tools 
library(mgcv)                                # GAM package by Simon Wood
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
#library(gstat)                               # Kriging and co-kriging by Pebesma et al. #not on NEX NASA
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
#library(automap)                             # Kriging automatic fitting of variogram using gstat #not on NEX NASA
#library(rgeos)                               # Geometric, topologic library of functions # not on NEX NASA
#RPostgreSQL                                 # Interface R and Postgres, not used in this script
library(gridExtra)
#Additional libraries not used in workflow
library(pgirmess)                            # Krusall Wallis test with mulitple options, Kruskalmc {pgirmess}  
#library(colorRamps)
  
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
  #create if does not exists: create the output dir as defined 
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  return(out_dir)
}

extract_list_from_list_obj<-function(obj_list,list_name){
  library(plyr)
  
  #Create a list of an object from a given list of object using a name prodived as input
  
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp <- obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]] <- as.data.frame(tmp) #deal with spdf...
  }
  return(list_tmp) #this is  a data.frame
}

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

## Function to mosaic modis or other raster images

mosaic_m_raster_list<-function(j,list_param){
  #This functions returns a subset of tiles from the modis grid.
  #Arguments: modies grid tile,list of tiles
  #Output: spatial grid data frame of the subset of tiles
  #Note that rasters are assumed to be in the same projection system!!
  
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
  }
  input.rasters <- lapply(m_list, raster)
  mosaiced_rast<-input.rasters[[1]]
  
  for (k in 2:length(input.rasters)){
    mosaiced_rast<-mosaic(mosaiced_rast,input.rasters[[k]], fun=mean)
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

### Function:
pred_data_info_fun <- function(k,list_data,pred_mod,sampling_dat_info){
    #Summarizing input info from sampling and df used in training/testing
      
  data <- list_data[[k]]
  sampling_dat <- sampling_dat_info[[k]]
  if(data!="try-error"){
    n <- nrow(data)
    n_mod <- vector("numeric",length(pred_mod))
    for(j in 1:length(pred_mod)){
      n_mod[j] <- sum(!is.na(data[[pred_mod[j]]]))
    }
    n <- rep(n,length(pred_mod))
    sampling_dat <- sampling_dat[rep(seq_len(nrow(sampling_dat)), each=length(pred_mod)),]
    row.names(sampling_dat) <- NULL
    df_n <- data.frame(n,n_mod,pred_mod)
    df_n <- cbind(df_n,sampling_dat)
  }else{        
    n <- rep(NA,length(pred_mod))
    n_mod <- vector("numeric",length(pred_mod))
    n_mod <- rep(NA,length(pred_mod))
    df_n <- data.frame(n,n_mod,pred_mod)
    sampling_dat <- sampling_dat[rep(seq_len(nrow(sampling_dat)), each=length(pred_mod)),]
    row.names(sampling_dat) <- NULL
    df_n <- data.frame(n,n_mod,pred_mod)
    df_n <- cbind(df_n,sampling_dat)
  
  }  
  return(df_n)
}

extract_daily_training_testing_info <- function(i,list_param){
  #This function extracts training and testing information from the raster object produced for each tile
  #This is looping through tiles...
  ##Functions used
  #Defined outside this script:
  #pred_data_info_fun
  
  ##### Parse parameters and arguments ####
  
  raster_obj_file <- list_param$list_raster_obj_files[i]
  use_month <- list_param$use_month
  use_day <- list_param$use_day
  tile_id <- list_param$list_names_tile_id[i]
  
  ### Start script ##
  
  raster_obj <- load_obj(unlist(raster_obj_file)) #may not need unlist
  nb_models <- length((raster_obj$clim_method_mod_obj[[1]]$formulas))
  pred_mod <- paste("mod",c(1:nb_models,"_kr"),sep="")
  #we are assuming no monthly hold out...
  #we are assuming only one specific daily prop?
  nb_models <- length(pred_mod)
  #names(raster_obj$method_mod_obj[[1]])
  var_interp <- unique(raster_obj$tb_diagnostic_s$var_interp)
  method_interp <- unique(raster_obj$tb_diagnostic_s$method_interp)
  
  if(use_day==TRUE){
    
    list_data_day_v <- try(extract_list_from_list_obj(raster_obj$validation_mod_obj,"data_v"))
    list_data_day_s <- try(extract_list_from_list_obj(raster_obj$validation_mod_obj,"data_s"))
    sampling_dat_day <- extract_list_from_list_obj(raster_obj$method_mod_obj,"daily_dev_sampling_dat")
    #debug(pred_data_info_fun)
    #list_pred_data_day_s_info <- pred_data_info_fun(1,list_data=list_data_day_s,pred_mod=pred_mod,sampling_dat_info=sampling_dat_day)
    list_pred_data_day_s_info <- lapply(1:length(sampling_dat_day),FUN=pred_data_info_fun,
           list_data=list_data_day_s,pred_mod=pred_mod,sampling_dat_info=sampling_dat_day)
    list_pred_data_day_v_info <- lapply(1:length(sampling_dat_day),FUN=pred_data_info_fun,
           list_data=list_data_day_v,pred_mod=pred_mod,sampling_dat_info=sampling_dat_day)
    pred_data_day_s_info <- do.call(rbind,list_pred_data_day_s_info)
    pred_data_day_v_info <- do.call(rbind,list_pred_data_day_v_info)
    pred_data_day_s_info$training <- rep(1,nrow(pred_data_day_s_info)) 
    pred_data_day_v_info$training <- rep(0,nrow(pred_data_day_v_info)) 
    pred_data_day_info <-rbind(pred_data_day_v_info,pred_data_day_s_info)
    
    pred_data_day_info$method_interp <- rep(method_interp,nrow(pred_data_day_info)) 
    pred_data_day_info$var_interp <- rep(var_interp,nrow(pred_data_day_info)) 
    pred_data_day_info$tile_id <- rep(tile_id,nrow(pred_data_day_info))                                       
  }
  
  if(use_month==TRUE){
    
    list_data_month_s <- try(extract_list_from_list_obj(raster_obj$validation_mod_month_obj,"data_s"))
    list_data_month_v <- try(extract_list_from_list_obj(raster_obj$validation_mod_month_obj,"data_v"))
    sampling_dat_month <- extract_list_from_list_obj(raster_obj$clim_method_mod_obj,"sampling_month_dat")
    list_pred_data_month_s_info <- lapply(1:length(sampling_dat_month),FUN=pred_data_info_fun,
           list_data=list_data_month_s,pred_mod=pred_mod,sampling_dat_info=sampling_dat_month)
    list_pred_data_month_v_info <- mclapply(1:length(sampling_dat_month),FUN=pred_data_info_fun,
           list_data=list_data_month_v,pred_mod=pred_mod,sampling_dat_info=sampling_dat_month,mc.preschedule=FALSE,mc.cores = 1)
    #Note for list_pred_data_month_v_info it will be try-error when there is no holdout.
    #combine training and testing later? also combined with accuracy
    pred_data_month_s_info <- do.call(rbind,list_pred_data_month_s_info)
    #Adding a column for training: if data item used as training then it is 1
    pred_data_month_s_info$training <- try(rep(1,nrow(pred_data_month_s_info)))

    #Remove try-error??
    list_pred_data_month_v_info_tmp <- remove_from_list_fun(list_pred_data_month_v_info)
    list_pred_data_month_v_info <- list_pred_data_month_v_info_tmp$list
    if (length(list_pred_data_month_v_info)>0){
      pred_data_month_v_info <- do.call(rbind,list_pred_data_month_v_info)
      pred_data_month_v_info$training <- try(rep(0,nrow(pred_data_month_v_info)))
      pred_data_month_info <- rbind(pred_data_month_v_info,pred_data_month_s_info)
    }else{
      pred_data_month_info <- pred_data_month_s_info
    }
    
    pred_data_month_info$method_interp <- rep(method_interp,nrow(pred_data_month_info)) 
    pred_data_month_info$var_interp <- rep(var_interp,nrow(pred_data_month_info)) 
    pred_data_month_info$tile_id <- rep(tile_id,nrow(pred_data_month_info)) 
  }    
    
  if(use_month==FALSE){
    pred_data_month_info <- NULL
  }
  if(use_day==FALSE){
    pred_data_day_info <- NULL
  }  
  #prepare object to return
  pred_data_info_obj <- list(pred_data_month_info,pred_data_day_info)
  names(pred_data_info_obj) <- c("pred_data_month_info","pred_data_day_info")
  #could add data.frame data_s and data_v later...
  ###
  return(pred_data_info_obj)
}

remove_from_list_fun <- function(l_x,condition_class ="try-error"){
  index <- vector("list",length(l_x))
  for (i in 1:length(l_x)){
    if (inherits(l_x[[i]],condition_class)){
      index[[i]] <- FALSE #remove from list
    }else{
      index[[i]] <- TRUE
    }
  }
  l_x<-l_x[unlist(index)] #remove from list all elements using subset
  
  obj <- list(l_x,index)
  names(obj) <- c("list","valid")
  return(obj)
}

##Function to list predicted tif
list_tif_fun <- function(i,in_dir_list,pattern_str){
  #list.files(path=x,pattern=".*predicted_mod1_0_1.*20100101.*.tif",full.names=T)})
  pat_str<- pattern_str[i]
  list_tif_files_dates <-lapply(in_dir_list,
         FUN=function(x,pat_str){list.files(path=x,pattern=pat_str,full.names=T)},pat_str=pat_str)
  return(list_tif_files_dates)
} 

calculate_summary_from_tb_diagnostic <-function(tb_diagnostic,metric_names){#,out_prefix,out_path){
  #now boxplots and mean per models
  library(gdata) #Nesssary to use cbindX
  
  ### Start script
  y_var_name<-unique(tb_diagnostic$var_interp) #extract the name of interpolated variable: dailyTmax, dailyTmin
  
  mod_names<-sort(unique(tb_diagnostic$pred_mod)) #models that have accuracy metrics
  t<-melt(tb_diagnostic,
          #measure=mod_var, 
          id=c("date","pred_mod","prop"),
          na.rm=F)
  t$value<-as.numeric(t$value) #problem with char!!!
  avg_tb<-cast(t,pred_mod~variable,mean)
  avg_tb$var_interp<-rep(y_var_name,times=nrow(avg_tb))
  median_tb<-cast(t,pred_mod~variable,median)
  
  #avg_tb<-cast(t,pred_mod~variable,mean)
  tb<-tb_diagnostic
 
  #mod_names<-sort(unique(tb$pred_mod)) #kept for clarity
  tb_mod_list<-lapply(mod_names, function(k) subset(tb, pred_mod==k)) #this creates a list of 5 based on models names
  names(tb_mod_list)<-mod_names
  #mod_metrics<-do.call(cbind,tb_mod_list)
  #debug here
  if(length(tb_mod_list)>1){
    mod_metrics<-do.call(cbindX,tb_mod_list) #column bind the list??
  }else{
    mod_metrics<-tb_mod_list[[1]]
  }
  
  test_names<-lapply(1:length(mod_names),function(k) paste(names(tb_mod_list[[1]]),mod_names[k],sep="_"))
  #test names are used when plotting the boxplot for the different models
  names(mod_metrics)<-unlist(test_names)
  rows_total<-lapply(tb_mod_list,nrow)
  avg_tb$n<-rows_total #total number of predictions on which the mean is based
  median_tb$n<-rows_total
  summary_obj<-list(avg_tb,median_tb)
  names(summary_obj)<-c("avg","median")
  return(summary_obj)  
}

#boxplot_month_from_tb(tb_diagnostic,metric_names,out_prefix,out_path)
## Function to display metrics by months/seasons
calculate_summary_from_tb_month_diagnostic <-function(tb_diagnostic,metric_names){ #,out_prefix,out_path){
  
  #Generate boxplot per month for models and accuracy metrics
  #Input parameters:
  #1) df: data frame containing accurayc metrics (RMSE etc.) per day)
  #2) metric_names: metrics used for validation
  #3) out_prefix
  #
  
  #################
  ## BEGIN
  y_var_name<-unique(tb_diagnostic$var_interp) #extract the name of interpolated variable: dailyTmax, dailyTmin  
  date_f<-strptime(tb_diagnostic$date, "%Y%m%d")   # interpolation date being processed
  tb_diagnostic$month<-strftime(date_f, "%m")          # current month of the date being processed
  mod_names<-sort(unique(tb_diagnostic$pred_mod)) #models that have accuracy metrics
  tb_mod_list<-lapply(mod_names, function(k) subset(tb_diagnostic, pred_mod==k)) #this creates a list of 5 based on models names
  names(tb_mod_list)<-mod_names
  t<-melt(tb_diagnostic,
          #measure=mod_var, 
          id=c("date","pred_mod","prop","month"),
          na.rm=F)
  t$value<-as.numeric(t$value) #problem with char!!!
  tb_mod_m_avg <-cast(t,pred_mod+month~variable,mean) #monthly mean for every model
  tb_mod_m_avg$var_interp<-rep(y_var_name,times=nrow(tb_mod_m_avg))
  
  tb_mod_m_sd <-cast(t,pred_mod+month~variable,sd)   #monthly sd for every model  
  tb_mod_m_list <-lapply(mod_names, function(k) subset(tb_mod_m_avg, pred_mod==k)) #this creates a list of 5 based on models names
  
  summary_month_obj <-c(tb_mod_m_list,tb_mod_m_avg,tb_mod_m_sd)
  names(summary_month_obj)<-c("tb_list","metric_month_avg","metric_month_sd")
  return(summary_month_obj)  
}

create_raster_prediction_obj<- function(in_dir_list,interpolation_method, y_var_name,out_prefix,out_path_list=NULL){
  
  
  #gather all necessary info:
  lf_validation_mod_month_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="gam_CAI_validation_mod_month_obj_dailyTmax.*.RData",full.names=T)})
  lf_validation_mod_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="gam_CAI_validation_mod_obj_dailyTmax.*.RData",full.names=T)})
  lf_clim_method_mod_obj <- lapply(in_dir_list,FUN=function(x){mixedsort(list.files(path=x,pattern="clim_obj_CAI_month_.*._TMAX_0_1_.*.RData",full.names=T))})
  lf_method_mod_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="method_mod_obj_gam_CAI_dailyTmax.*.RData",full.names=T)})

  tb_month_diagnostic_v_list <- mclapply(lf_validation_mod_month_obj,FUN=function(x){try( x<- load_obj(x)); try(extract_from_list_obj(x,"metrics_v"))},mc.preschedule=FALSE,mc.cores = 6)                           
  tb_month_diagnostic_s_list <- mclapply(lf_validation_mod_month_obj,FUN=function(x){try( x<- load_obj(x)); try(extract_from_list_obj(x,"metrics_s"))},mc.preschedule=FALSE,mc.cores = 6)                           
  tb_diagnostic_v_list <- mclapply(lf_validation_mod_obj,FUN=function(x){try( x<- load_obj(x)); try(extract_from_list_obj(x,"metrics_v"))},mc.preschedule=FALSE,mc.cores = 6)                           
  tb_diagnostic_s_list <- mclapply(lf_validation_mod_obj,FUN=function(x){try( x<- load_obj(x)); try(extract_from_list_obj(x,"metrics_s"))},mc.preschedule=FALSE,mc.cores = 6)                           

  #rownames(tb_month_diagnostic_v)<-NULL #remove row names
  #tb_month_diagnostic_v$method_interp <- interpolation_method

  #Call functions to create plots of metrics for validation dataset
  metric_names<-c("rmse","mae","me","r","m50")
  summary_metrics_v_list <- lapply(tb_diagnostic_v_list,FUN=function(x){try(calculate_summary_from_tb_diagnostic(x,metric_names))})
  summary_month_metrics_v_list <- lapply(tb_diagnostic_v_list,FUN=function(x){try(calculate_summary_from_tb_month_diagnostic(x,metric_names))}) 
  #list_summalist_summary_metrics_vry_month_metrics_v <- calculate_summary_from_tb_month_diagnostic(list_tb_diagnostic_v[[1]],metric_names)

  #if (interpolation_method %in% c("gam_CAI","kriging_CAI","gwr_CAI","gam_fusion","kriging_fusion","gwr_fusion")){
  lf_raster_obj <- vector("list",length=length(in_dir_list))
  for (i in 1:length(in_dir_list)){
    
    clim_method_mod_obj <- try(lapply(lf_clim_method_mod_obj[[i]],FUN=try(load_obj)))
    method_mod_obj <- try(load_obj(lf_method_mod_obj[[i]]))
    validation_mod_month_obj <- try(load_obj(lf_validation_mod_month_obj[[i]]))   
    validation_mod_obj <- try(load_obj(lf_validation_mod_obj[[i]]))

    tb_month_diagnostic_v <- try(tb_month_diagnostic_v_list[[i]])
    tb_month_diagnostic_s <- try(tb_month_diagnostic_s_list[[i]])
    tb_diagnostic_v <- try(tb_diagnostic_v_list[[i]])
    tb_diagnostic_s <- try(tb_diagnostic_s_list[[i]])

    summary_metrics_v <- summary_metrics_v_list[[i]]
    summary_month_metrics_v <- summary_month_metrics_v_list[[i]] 
    
    raster_prediction_obj <- list(clim_method_mod_obj,method_mod_obj,validation_mod_obj,validation_mod_month_obj, tb_diagnostic_v,
                                tb_diagnostic_s,tb_month_diagnostic_v,tb_month_diagnostic_s,summary_metrics_v,summary_month_metrics_v)
    names(raster_prediction_obj) <-c ("clim_method_mod_obj","method_mod_obj","validation_mod_obj","validation_mod_month_obj","tb_diagnostic_v",
                                "tb_diagnostic_s","tb_month_diagnostic_v","tb_month_diagnostic_s","summary_metrics_v","summary_month_metrics_v") 
    if(is.null(out_path_list)){
      out_path <- in_dir_list[[i]]
    }  
    if(!is.null(out_path_list)){
      out_path <- out_path_list[[i]]
    }  
    lf_raster_obj[[i]] <- file.path(out_path,paste("raster_prediction_obj_",interpolation_method,"_", y_var_name,out_prefix_str[i],".RData",sep=""))  

    save(raster_prediction_obj,file= file.path(out_path,paste("raster_prediction_obj_",interpolation_method,"_", y_var_name,out_prefix_str[i],".RData",sep="")))
  }
  return(lf_raster_obj)
}

##################### END OF SCRIPT ######################


