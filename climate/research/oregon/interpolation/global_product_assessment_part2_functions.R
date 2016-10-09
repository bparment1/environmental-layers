####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Assessment of product part 2 functions: mosaic and accuracy ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This part 2 of the assessment focuses on graphics to explore the spatial patterns of raster times series as figures and movie.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 05/24/2016  
#MODIFIED ON: 10/09/2016            
#Version: 2
#PROJECT: Environmental Layers project     
#COMMENTS: Initial commit, script based on part NASA biodiversity conference 
#TODO:
#1) Add plot broken down by year and region 
#2) Modify code for overall assessment accross all regions and year
#3) Clean up

#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rmd user:aguzman4:rwx /nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015

##COMMIT: adding function to generate animation

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
  if(class(date_proc)!="Date"){
    date_val <- as.Date(strptime(date_proc,"%Y%m%d"))
    #month_name <- month.name(date_val)
  }else{
    date_val <- date_proc
  }
  
  month_str <- format(date_val, "%b") ## Month, char, abbreviated
  year_str <- format(date_val, "%Y") ## Year with century
  day_str <- as.numeric(format(date_val, "%d")) ## numeric month
  date_str <- paste(month_str," ",day_str,", ",year_str,sep="")
  
  res_pix <- 1200
  #res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1
  
  
  if(is.null(zlim_val)){
    
    if(is.na(minValue(r_pred))){
      r_pred <- setMinMax(r_pred)
    }
    
    zlim_val_str <- paste(c(minValue(r_pred),maxValue(r_pred)),sep="_",collapse="_")
    #png_filename <-  file.path(out_dir,paste("Figure4_clim_mosaics_day_","_",date_proc,"_",region_name,"_zlim_",zlim_val_str,"_",out_suffix,".png",sep =""))
    #raster_name_tmp
    png_filename <-  file.path(out_dir,paste("Figure_",raster_name_tmp,"_zlim_",zlim_val_str,"_",out_suffix,".png",sep =""))
    
    title_str <-  paste("Predicted ",variable_name, " on ",date_str , " ", sep = "")
    #browser()
    
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
  return(missing_dates)
}

#create animation from figures:
generate_animation_from_figures_fun <- function(filenames_figures,frame_speed=60,format_file=".gif",out_suffix="",out_dir=".",out_filename_figure_animation=NULL){

  if(is.null(out_filename_figure_animation)){
    #out_filename_figure_movies <- file.path(out_dir,paste("mosaic_movie_",out_suffix_movie,".gif",sep=""))
    out_filename_figure_animation <- file.path(out_dir,paste("animation_frame_",frame_speed,"_",out_suffix,format_file,sep=""))
  }
  
  if(class(filenames_figures)=="list"){
    #filename_figures_mosaic <- file.path(out_dir,"mosaic_plot_fig.txt")
    out_filenames_figures <- file.path(out_dir,paste("list_figures_animation_",out_suffix,".txt",sep=""))
    write.table(unlist(filenames_figures),out_filenames_figures,row.names = F,col.names = F,quote = F)
    filenames_figures <- out_filenames_figures
  }
  
  #now generate movie with imageMagick

  #-delay 20
  #delay_option <- 60
  delay_option <- frame_speed
  
  cmd_str <- paste("convert",
                   paste("-delay",delay_option),
                   paste0("@",filenames_figures),
                   out_filename_figure_animation)
  #convert @myimages.txt mymovie.gif
  #save cmd_str in text file!!!
  
  system(cmd_str)
  
  return(out_filename_figure_animation)

}

############################ END OF SCRIPT ##################################