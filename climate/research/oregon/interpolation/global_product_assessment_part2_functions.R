####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Assessment of product part 2 functions: mosaic and accuracy ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This part 2 of the product assessment focuses on graphics to explore the spatial patterns of raster times series as figures.
#The file contains functions to genrate figures and animation (movie).
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/03/2016  
#MODIFIED ON: 12/13/2017            
#Version: 2
#PROJECT: Environmental Layers project     
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

##COMMIT: modifying function to check missing files and dates for predictions and others

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
  #Function to plot raster image
  #
  #INPUTS
  #1) l_dates
  #2) r_stack
  #3) NA_flag_val
  #4) out_dir,
  #5) out_suffix_str
  #6) region_name
  #7) variable_name
  #8) zlim_val
  #
  
  ############# Start script #########
  
  l_dates <- list_param$l_dates
  r_mosaiced_scaled <- list_param$r_mosaiced_scaled
  NA_flag_val <- list_param$NA_flag_val
  out_dir <- list_param$out_dir
  out_suffix <- list_param$out_suffix
  region_name <- list_param$region_name
  variable_name <- list_param$variable_name
  zlim_val <- list_param$zlim_val
  stat_opt <- list_param$stat_opt #if TRUE computer stats for the image
  
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
    
    min_max_val <- c(minValue(r_pred),maxValue(r_pred))
    min_max_df <- data.frame(min=min_max_val[1],max=min_max_val[2])
    
    zlim_val_str <- paste(min_max_val,sep="_",collapse="_")
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
    min_max_val <- NULL
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
  
  ### Option to compute stats for the raster image
  if(stat_opt==TRUE){
    if(is.null(min_max_val)){
      min_val <- cellStats(r_pred,stat="min",na.rm=T)
      max_val <- cellStats(r_pred,stat="max",na.rm=T)
    }else{
      min_val <- min_max_val[1]
      max_val <- min_max_val[2]
    }

    NA_no <- freq(r_pred,value=NA)
    n_cell <- ncell(r_pred)
    perc_NA <- 100*(NA_no/n_cell)
    if(dataType(r_pred)=="INT2S" | dataType(r_pred)=="INT4S"){
      #integer overflow for sum
      #filename(r_pred)
      #r_pred_tmp <- r_pred
      raster_name_tmp <- paste(raster_name_tmp,"_tmp",file_format,sep="")
      writeRaster(r_pred_tmp,raster_name_tmp,datatype="FLT8S",overwrite=T) #this very slow
      r_pred_tmp <-raster(raster_name_tmp)
      #dataType(r_pred_tmp)
      mean_val <- cellStats(r_pred_tmp,stat="mean",na.rm=T)
      sd_val <- cellStats(r_pred_tmp,stat="sd",na.rm=T)
      file.remove(raster_name_tmp)
    }else{
      mean_val <- cellStats(r_pred,stat="mean",na.rm=T)
      sd_val <- cellStats(r_pred,stat="sd",na.rm=T)
    }
    stat_df <- data.frame(min=min_val,
                          max=max_val,
                          sd=sd_val,
                          mean=mean_val,
                          NA_no=NA_no,
                          n_cell=n_cell,
                          perc_NA)
  }else{
    stat_df <- NULL
  }
  
  fig_obj <- c(png_filename,min_max_df,stat_df)
  names(fig_obj) <- c("png_filename","min_max_df","stat_df")
  
  #return(png_filename)
  return(fig_obj)
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
                                         x = basename(lf),
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


#create animation from figures:
generate_animation_from_figures_fun <- function(filenames_figures,frame_speed=50,format_file=".gif",out_suffix="",out_dir=".",out_filename_figure_animation=NULL){
  #This function generates an animation given a list of files or textfile. The default format is .gif.
  #The function requires ImageMagick to produce animation.
  #INPUTS:
  #1) filenames_figures: list of files as "list" or "character, or file name of a text file containing the list of figures.
  #2) frame_speed: delay option in constructing the animation, default is 50,
                  #the unit is 1/100th of second so 50 is 2 frame per second
  #3) format_file=".gif"
  #4) out_suffix=""
  #5) out_dir=".",
  #6) out_filename_figure_animation=NUL
  #OUTPUTS:
  #

  if(is.null(out_filename_figure_animation)){
    #out_filename_figure_movies <- file.path(out_dir,paste("mosaic_movie_",out_suffix_movie,".gif",sep=""))
    out_filename_figure_animation <- file.path(out_dir,paste("animation_frame_",frame_speed,"_",out_suffix,format_file,sep=""))
  }
  
  if(class(filenames_figures)=="list" | (length(filenames_figures)>1)){ #if length of object greater than 1 then assume a list of files
    #filename_figures_mosaic <- file.path(out_dir,"mosaic_plot_fig.txt")
    out_filenames_figures <- file.path(out_dir,paste("list_figures_animation_",out_suffix,".txt",sep=""))
    write.table(unlist(filenames_figures),out_filenames_figures,row.names = F,col.names = F,quote = F)
    filenames_figures <- out_filenames_figures
  }
  
  #now generate movie with imageMagick

  if(file_format==".gif"){
    #-delay 20
    #delay_option <- 60
    delay_option <- frame_speed
    
    cmd_str <- paste("convert",
                     paste("-delay",delay_option),
                     paste0("@",filenames_figures),
                     out_filename_figure_animation)
    #convert @myimages.txt mymovie.gif
    #save cmd_str in text file!!!
    
  }
  
  #if format is .mp4
  if(file_format==".mp4"){
    #ffmpeg -f image2 -pattern_type glob -i '*.png' out.mp4
    #ffmpeg -f image2 -r 1 -pattern_type glob -i '*.png' out.mp4
    
    #-r1 one second by frame
    #rate of one per second:
    #ffmpeg -f image2 -r 1 -pattern_type glob -i '*.png' out.mp4
    #crf: used for compression count rate factor is between 18 to 24, the lowest is the highest quality
    #ffmpeg -f image2 -r 1 -vcodec libx264 -crf 24 -pattern_type glob -i '*.png' out.mp4
    frame_rate <- delay_option/100 # convert frame rate in second
    cmd_str <- paste("ffmpeg",
                       paste("-f","image2",sep=" "),
                       paste("-r",frame_rate,sep=" "),
                       paste("-pattern_type glob -i","'*.png'",sep=" "),
                       #paste("-pattern_type glob -i ",filenames_figures),
                       out_filename_figure_animation,sep=" ")
    
  }
  system(cmd_str)
  
  return(out_filename_figure_animation)

}


plot_and_animate_raster_time_series <- function(lf_raster, item_no,region_name,var_name,metric_name,NA_flag_val,filenames_figures=NULL,frame_speed=60,animation_format=".gif",zlim_val=NULL,plot_figure=T,generate_animation=T,stat_opt,num_cores=2,out_suffix="",out_dir="."){
  #Function to generate figures and animation for a list of time series raster files.
  #The list assumes that dates can be extracted from the filename (as in the NEX worklow)
  #
  #INPUTS
  #1) lf_raster
  #2) item_no
  #3) region_name,
  #4) var_name
  #5) metric_name
  #6) NAvalue
  #7) filenames_figures
  #8) frame_speed
  #9) animation_format
  #10) zlim_val
  #11) plot_figure
  #12) generate_animation
  #13) stat_opt: if TRUE compute general stat for raster (mi, max, mean etc.)
  #14) num_cores
  #15) out_suffix
  #16) out_dir
  #OUTPUTS
  #An object as list called figure_animation_obj with:
  #1) filenames_figures_mosaic
  #2) out_filename_figure_animation
  
  ####### Start script #############
  
  
  #lf_mosaic_list <- lf_raster
  variable_name <- var_name
  
  if(!is.null(plot_figure)){
    #item_no <- 13
    list_dates_produced <- unlist(mclapply(1:length(lf_raster),
                                           FUN = extract_date,
                                           x = basename(lf_raster),
                                           item_no = item_no,
                                           mc.preschedule = FALSE,
                                           mc.cores = num_cores))
    
    list_dates_produced_date_val <- as.Date(strptime(list_dates_produced, "%Y%m%d"))
    month_str <- format(list_dates_produced_date_val, "%b") ## Month, char, abbreviated
    year_str <- format(list_dates_produced_date_val, "%Y") ## Year with century
    day_str <- as.numeric(format(list_dates_produced_date_val, "%d")) ## numeric month
    df_raster <- data.frame(lf = basename(lf_raster),
                            date = list_dates_produced_date_val,
                            month_str = month_str,
                            year = year_str,
                            day = day_str,
                            dir = dirname(lf_raster))
    
    df_raster_fname <- file.path(out_dir, paste0("df_raster_", out_suffix, ".txt"))
    write.table(df_raster,file = df_raster_fname,sep = ",",row.names = F)
    
    ############### PART5: Make raster stack and display maps #############
    #### Extract corresponding raster for given dates and plot
    
    r_stack <- stack(lf_raster,quick=T)
    l_dates <- list_dates_produced_date_val #[1:11]
    
    #undebug(plot_raster_mosaic)
    #zlim_val <- zlim_val
    
    ### Now run for the full time series
    #13.26 Western time: start
    #l_dates <- list_dates_produced_date_val
    #r_stack_subset <- r_stack
    #zlim_val <- NULL
    out_suffix_str <- paste0(var_name,"_",metric_name,"_",out_suffix)
    list_param_plot_raster_mosaic <- list(l_dates,r_stack,
                                          NA_flag_val,out_dir,
                                          out_suffix_str,region_name,
                                          variable_name,zlim_val,stat_opt)
    names(list_param_plot_raster_mosaic) <- c("l_dates","r_mosaiced_scaled",
                                              "NA_flag_val_mosaic","out_dir",
                                              "out_suffix", "region_name",
                                              "variable_name","zlim_val","stat_opt")
    #debug(plot_raster_mosaic)
    browser()
    lf_mosaic_plot_fig <- plot_raster_mosaic(1,
                                   list_param = list_param_plot_raster_mosaic)
    
    ##start at 12.29
    ##finished at 15.23 (for reg 6 with 2,991 figures)
    lf_mosaic_plot_fig <- mclapply(1:length(l_dates),
                                   FUN = plot_raster_mosaic,
                                   list_param = list_param_plot_raster_mosaic,
                                   mc.preschedule = FALSE,
                                   mc.cores = num_cores)
    
    if (is.null(zlim_val)) {
      out_suffix_movie <- paste("min_max_", out_suffix_str, sep = "")
    } else{
      zlim_val_str <- paste(zlim_val, sep = "_", collapse = "_")
      out_suffix_movie <- paste(zlim_val_str, "_", out_suffix, sep = "")
    }
    filenames_figures_mosaic <- paste0("list_figures_animation_", out_suffix_movie, ".txt")
    
    write.table(unlist(lf_mosaic_plot_fig),filenames_figures_mosaic,row.names = F,col.names = F,quote = F)
    
  }
  
  ### Part 2 generate movie
  
  if(generate_animation==TRUE){
    
    out_suffix_str <- paste0(var_name,"_",metric_name,"_",out_suffix)
    
    if (is.null(zlim_val)) {
      out_suffix_movie <- paste("min_max_", out_suffix, sep = "")
    } else{
      zlim_val_str <- paste(zlim_val, sep = "_", collapse = "_")
      out_suffix_movie <- paste(zlim_val_str, "_", out_suffix, sep = "")
    }
    
    #already provided as a parameter
    #filenames_figures_mosaic <- paste0("list_figures_animation_", out_suffix_movie, ".txt")
    #write.table(unlist(lf_mosaic_plot_fig),filenames_figures_mosaic,row.names = F,col.names = F,quote = F)
    
    #now generate movie with imageMagick
    #frame_speed <- 60
    #animation_format <- ".gif"
    #out_suffix_str <- out_suffix
    #started
    #debug(generate_animation_from_figures_fun)
    #lf_mosaic_plot_fig <- read.table(filenames_figure,sep=",")
    
    #out_filename_figure_animation <- generate_animation_from_figures_fun(filenames_figures = unlist(lf_mosaic_plot_fig[1:11]),
    #                                    frame_speed = frame_speed,
    #                                    format_file = animation_format,
    #                                    out_suffix = out_suffix_str,
    #                                    out_dir = out_dir,
    #                                    out_filename_figure_animation = "test2_reg6_animation.gif")
    
    #started 17.36 Western time on Oct 10 and 18.18
    #        15.58                 oct 11 for 16.38 for reg6 pred (about 2991)
    #out_filename_figure_animation
    out_filename_figure_animation <- generate_animation_from_figures_fun(filenames_figures = filenames_figures_mosaic,
                                                                         frame_speed = frame_speed,
                                                                         format_file = animation_format,
                                                                         out_suffix = out_suffix_movie,
                                                                         out_dir = out_dir,
                                                                         out_filename_figure_animation = NULL)
  }
  
  ## prepare object to return
  
  figure_animation_obj <- list(filenames_figures_mosaic,out_filename_figure_animation)
  names(figure_animation_obj) <- c("filenames_figures_mosaic","out_filename_figure_animation")
  return(figure_animation_obj)
  
}

############################ END OF SCRIPT ##################################