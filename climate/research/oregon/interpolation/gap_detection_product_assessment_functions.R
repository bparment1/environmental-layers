####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part 0 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#This script provides function for an overall assessment of gaps by regions.
#It uses global product assessment of the number of predictions by tiles and years.
#The general logic is to check the number of overlap by shapefile polyon tiles
#along with the predictions for every day of the year (*.tif)
#Summary tables and data are also produced in the script for each region over 31 years.
#
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/31/2016  
#MODIFIED ON: 06/19/2017            
#Version: 1
#PROJECT: Environmental Layers project     
#COMMENTS: removing unused functions and clean up for part0 global product assessment part0 
#TODO:#PROJECT: Environmental Layers project     
#COMMENTS:
#TODO:
#
#First source these files:
#Resolved call issues from R.
#source /nobackupp6/aguzman4/climateLayers/sharedModules2/etc/environ.sh 
#
#setfacl -Rmd user:aguzman4:rwx /nobackupp8/bparmen1/output_run10_1500x4500_global_analyses_pred_1992_10052015

##COMMIT: testing command line

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

plot_raster_poly_overlap <- function(shps_tiles,list_lf_raster_tif_tiles,df_missing,num_cores=1,
                                     mosaic_python_script,data_tye_str,
                                     region_name="",out_suffix="",out_dir="."){
  
  ###This functions generate a mosaic of overlap for tiles of SpatialPolygonsDataFrame.
  ##INPUTS
  ##
  ##TO DO: modify to drop df_missing
  
  ########################### Beging script ##################
  
  ### if list of raster exist check that none are missing
  if(!is.null(list_lf_raster_tif_tiles)){
    list_r_ref <- mclapply(1:length(list_lf_raster_tif_tiles), 
                           function(i,x){try(raster(x[[i]][1]))},
                           x = list_lf_raster_tif_tiles,
                           mc.preschedule=FALSE,
                           mc.cores = num_cores)
    
    #find try-error: missing raster
    list_r_ref_error <- as.numeric(unlist(lapply(list_r_ref,function(x){class(x)=="try-error"})))
    
    ## Select tiles without raster, generate raster and from r_mask
    
    ref_test<- mclapply(1:length(shps_tiles),
                        FUN=generate_raster_tile_ref,
                        shps_tiles = shps_tiles,
                        r_mask=r_mask,
                        list_r_ref_error=list_r_ref_error,
                        mc.preschedule=FALSE,
                        mc.cores = num_cores)
    
    #now fill in list_r_ref with ref_test
    
    list_r_ref <- lapply(1:length(list_r_ref),
                         FUN = replace_raster_ref,
                         list_r_ref=list_r_ref,
                         ref_missing = ref_test)
    
  }else{
    #generate raster from polygon tiles for al the list
    ref_test<- mclapply(1:length(shps_tiles),
                        FUN=generate_raster_tile_ref,
                        shps_tiles = shps_tiles,
                        r_mask=r_mask,
                        list_r_ref_error=list_r_ref_error,
                        mc.preschedule=FALSE,
                        mc.cores = num_cores)
  }
  
  tile_spdf <- shps_tiles[[1]]
  #tile_coord <- basename(in_dir_reg[1])
  date_val <- df_missing$date[1] # need to get rid of this parameters later to make this a general function
  
  #browser()
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
                             df_missing=df_missing, #modify later to get rid of this input
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
  #mosaic_python <- "/nobackupp6/aguzman4/climateLayers/sharedCode/"
  mosaic_python <- dirname(mosaic_python_script)
  rast_ref_name <- infile_mask
  #python /nobackupp6/aguzman4/climateLayers/sharedCode//gdal_merge_sum.py --config GDAL_CACHEMAX=1500 --overwrite=TRUE -o /nobackupp8/bparmen1/climateLayers/out
  mosaic_overlap_tiles_obj <- mosaic_python_merge(NA_flag_val=NA_flag_val,
                                                  module_path=mosaic_python,
                                                  #module_name="gdal_merge_sum.py",
                                                  module_name = basename(mosaic_python_script),
                                                  input_file=filename_list_predicted,
                                                  out_mosaic_name=out_mosaic_name_predicted_m,
                                                  raster_ref_name = rast_ref_name) ##if NA, not take into account
  r_overlap_raster_name <- mosaic_overlap_tiles_obj$out_mosaic_name
  cmd_str1 <-   mosaic_overlap_tiles_obj$cmd_str
  
  #browser()
  
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
  
  ####
  
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
  #### prepare object to return
  
  return()
}



gap_tiles_assessment_fun <- function(in_dir,y_var_name,region_name,num_cores,NA_flag_val,data_type_str,
                                     shps_tiles,list_lf_raster_tif_tiles,infile_mask,countries_shp,
                                     moscaic_python_script,out_dir,out_suffix){
  #
  #This function assesses missing tiles over 31 years of predictions using output from product assessment.
  #It is to be run region by region.
  #
  #INPUTS
  #1) in_dir
  #2) y_var_name
  #3) region_name
  #4) num_cores
  #5) NA_flag_val
  #6) data_type_str
  #7) shps_tile
  #8) list_lf_raster_tif_tiles
  #9) infile_mask
  #10) countries_shp,
  #11) moscaic_python_script
  #12) out_dir
  #13) out_suffix
  #OUTPUTS
  #
  #

  ###################### Begin script ##########################
  
  ##Select relevant folder/dir by region given input dir
  in_dir_list_tmp <- list.files(pattern =paste0(".*.",region_name,".*."),full.names=T,path=in_dir)
  ##collect data.frame with missing information from year assessments
  list_lf_df_missing_tiles <- unlist(mclapply(in_dir_list_tmp,
                                              FUN=function(x){list.files(path=x,
                                                                         pattern=paste0("df_missing_by_dates_tiles_predicted_tif_",".*.",region_name,".*.txt"),full.names=T)},
                                              mc.preschedule=FALSE,mc.cores = num_cores))
  ## read in data.frame with missing information
  list_df_missing_tiles <- mclapply(list_lf_df_missing_tiles,
                                    FUN=function(x){read.table(x,sep=",",stringsAsFactors = F,check.names = F)},
                                    mc.preschedule=FALSE,mc.cores = num_cores)
  ## Combine all annual data.frame over the years
  df_missing_tiles_reg <- do.call(rbind,list_df_missing_tiles)

  ###### Assessment ####
  #### Generate summary from missing table
  #1) Select dates with missing tiles
  #2) Plot number of days missing in maps over 365 days
  #3) Plot number of days missing in histograms/barplots
  #4) Keep raster of number pix predictions of overlap
  ### do sum across tiles to find number of missing per tiles and map it
  
  n_tiles <- length(df_missing_tiles_reg) - 3 #number of tiles in this region
  df_missing_tiles_reg_sp <- t(df_missing_tiles_reg[,1:n_tiles]) #transpose to get df with lines centroids of tiles
  df_missing_tiles_reg_sp <- as.data.frame(df_missing_tiles_reg_sp)
  names(df_missing_tiles_reg_sp) <- df_missing_tiles_reg$date
  df_missing_tiles_reg_sp$tot_missing <- rowSums(df_missing_tiles_reg_sp) #total missing over a year by tile
  df_missing_tiles_reg_sp$tot_pred <- length(df_missing_tiles_reg$date) - df_missing_tiles_reg_sp$tot_missing
  
  tiles_names <- names(df_missing_tiles_reg)[1:n_tiles]
  list_xy <- strsplit(unlist(basename(tiles_names)),"_")

  coord_xy <- do.call(rbind,list_xy)
  y_val <- as.numeric(coord_xy[,1]) #lat, long! need to reverse
  x_val <- as.numeric(coord_xy[,2])
  coordinates(df_missing_tiles_reg_sp) <- as.matrix(cbind(x_val,y_val))
  
  #This contains in rows date and tiles columns
  filename_df_missing_tiles_sp <- file.path(out_dir,paste0("df_missing_by_centroids_tiles_and_dates",region_name,"_",pred_mod_name,"_",out_suffix,".txt"))
  write.table(as.data.frame(df_missing_tiles_reg_sp),file=filename_df_missing_tiles_sp,sep=",")

  #browser()
  
  ### Now generate plots if missing tiles for specific dates in the region
  df_missing_tiles_day <- subset(df_missing_tiles_reg,tot_missing > 0)
  
  path_to_shp <- dirname(countries_shp)
  layer_name <- sub(".shp","",basename(countries_shp))
  reg_layer <- readOGR(path_to_shp, layer_name) #outlines of the region
  r_mask <- raster(infile_mask)
  
  ##### write basic stats in textfile?
  #sum(df_missing_tiles_reg)
  range(df_missing_tiles_reg$tot_missing)
  #sink()
  missing_val <- table(df_missing_tiles_reg$tot_missing) #save this info!!
   
  if(nrow(df_missing_tiles_day)>0){
    
    res_pix <- 800
    #res_pix <- 480
    col_mfrow <- 1
    row_mfrow <- 1
    png_filename_histogram <-  file.path(out_dir,paste("Figure_barplot_",region_name,"region_missing_tiles","_",out_suffix,".png",sep =""))
    
    png(filename=png_filename_histogram,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
    barplot(missing_val,
            ylab="frequency of missing",
            xlab="number of missing day predictions",
            main=paste0("Number of missing predictions over 31 years for ",region_name))
    dev.off()
    #check for this in every output, if not present then there are no missing tiles over the full year for 
    #the specific region
    write.table(df_missing_tiles_day,file=paste0("df_missing_tiles_day_mosaic_",out_suffix,".txt"))
    
    ### Now output the table:
    out_filename <- paste0("missing_val_table_barplot_day_tiles_",out_suffix,".txt")
    write.table(missing_val,file= file.path(out_dir,out_filename))
    #do spplot after that on tot sum
    
    png(filename=paste("Figure_total_missing_days_map_centroids_tile_",pred_mod_name,"_",
                       "_",out_suffix,".png",sep=""),
        width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    #spplot for reg_layer not working on NEX again
    #p_shp <- spplot(reg_layer,"ISO3" ,col.regions=NA, col="black") #ok problem solved!!
    p_r <-levelplot(r_mask,colorkey=F) #no key legend
    p <- bubble(df_missing_tiles_reg_sp,"tot_missing",main=paste0("Missing per tile and by ",pred_mod_name,
                                                              " for ",y_var_name))
    #p1 <- p+p_shp
    p_c <- p + p_r + p #set the legend first by using p first
    
    try(print(p_c)) #error raised if number of missing values below a threshold does not exist
    dev.off()
    
  }
  
  #do spplot after that on tot sum
  res_pix <- 800
  #res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1
    
  png(filename=paste("Figure_total_predicted_days_predicted_map_centroids_tile_",pred_mod_name,"_",
                       "_",out_suffix,".png",sep=""),
  width=col_mfrow*res_pix,height=row_mfrow*res_pix)
    
  p_r <-levelplot(r_mask,colorkey=F) #no key legend
  p <- bubble(df_missing_tiles_sp,"tot_pred",main=paste0("Prediction per tile and by ",pred_mod_name,
                                                           " for ", y_var_name))
  p_c <- p + p_r + p #set the legend first by using p first
  #p1 <- p+p_shp
  try(print(p_c)) #error raised if number of missing values below a threshold does not exist
  dev.off()
  
  #### add more here
  if(raster_overlap==TRUE){
    
    ### preparing inputs for raster_overlap production
    names(shps_tiles) <- tiles_names
    #r_ref <- raster(list_lf_raster_tif_tiles[[15]][1])
    #list_r_ref <- lapply(1:length(in_dir_reg), function(i){raster(list_lf_raster_tif_tiles[[i]][1])})
    
    ####  

    plot_raster_poly_overlap(shps_tiles,list_lf_raster_tif_tiles,df_missing,num_cores=1,
                                       mosaic_python_script,data_tye_str,
                                       region_name="",out_suffix="",out_dir=".")
    
  }else{ #if raster_overalp==FALSE
    out_mosaic_name_overlap_masked <- NULL
    tb_freq_overlap <- NULL
    png_filename_maximum_overlap <- NULL
  }
  
  ##### Prepare object to return
  return()
}


############################# END OF SCRIPT ###################################