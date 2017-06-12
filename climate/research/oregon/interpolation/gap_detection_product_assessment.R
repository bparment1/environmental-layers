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
#MODIFIED ON: 06/07/2017            
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




var <- "TMIN" # variable being interpolated #PARAM 1, arg 1
in_dir <- "/nobackupp6/aguzman4/climateLayers/tMinOut/reg1/assessment" #PARAM2
region_name <- c("reg1") #PARAM 3, arg 3
out_suffix <- "predictions_gaps_tiles_assessment_reg1_1993" #PARAM 4
#out_suffix_str <- region_name #PARAM 4, CONST 3
out_dir <- "/nobackupp8/bparmen1/climateLayers/tMinOut/reg1/assessment"
create_out_dir_param <- TRUE #PARAM 12, arg 6
year_predicted <- c(1993) #PARAM 7, arg7
num_cores <- 6 #number of cores used # PARAM 8, arg 8
max_mem <- 1e+07 #PARAM 9
item_no <- 9 #PARAM 10, arg 10
metric_name <- "rmse" # "mae", "r" for MAE, R etc.; can also be ns or nv? #PARAM 11, arg 11
day_start <- "19930101" #PARAM 12, arg 12
day_end <- "19931231" #PARAM 13, arg 13
infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg1.tif" #PARAM 14, arg 14
in_dir1 <- "/nobackupp6/aguzman4/climateLayers/tMinOut" # PARAM 15 On NEX
layers_option <- c("var_pred") #PARAM 16, arg 16
tmp_files <- FALSE #PARAM 17, arg 17
plotting_figures <- TRUE #PARAm 18, arg 18
raster_overlap <- FALSE # PARAM 19, if TRUE, raster overlap is generated
raster_pred <- FALSE # PARAM 20, if TRUE, raster prediction is generated


#################### Begin Script ######
in_dir <- "/nobackupp6/aguzman4/climateLayers/tMinOut/testGaps"
# ./*reg1*/df_missing_by_dates_tiles_predicted_tif_reg1_mod1_predictions_gaps_tiles_assessment_reg1_*.txt 
num_cores <- 6
pattern_str <- "df_missing_by_dates_tiles_predicted_tif_reg1_mod1_predictions_gaps_tiles_assessment_reg1_.*.txt"

function(in_dir,region_name,num_cores,out_dir, out_suffix){
  #d
  #d
  
}
## select relevant files for region
in_dir_list_tmp <- list.files(pattern =paste0(".*.",region_name,".*."),full.names=T,path=in_dir)

list_lf_df_missing_tiles <- unlist(mclapply(in_dir_list_tmp,
                                     FUN=function(x){list.files(path=x,
                                                                pattern=paste0("df_missing_by_dates_tiles_predicted_tif_",".*.",region_name,".*.txt"),full.names=T)},
                                     mc.preschedule=FALSE,mc.cores = num_cores))
list_df_missing_tiles <- mclapply(list_lf_df_missing_tiles,
                                     FUN=function(x){read.table(x,sep=",",stringsAsFactors = F,check.names = F)},
                                     mc.preschedule=FALSE,mc.cores = num_cores)

df_missing_tiles_reg <- do.call(rbind,list_df_missing_tiles)
sum(df_missing_tiles_reg)

range(df_missing_tiles_reg$tot_missing)
table(df_missing_tiles_reg$tot_missing)
histogram(df_missing_tiles_reg$tot_missing)


###### Assessment ####
#### Generate summary from missing table
#1) Select dates with missing tiles
#2) Plot number of days missing in maps over 365 days
#3) Plot number of days missing in histograms/barplots
#4) Keep raster of number pix predictions of overlap
### do sum across tiles to find number of missing per tiles and map it

n_tiles <- length(df_missing_tiles_reg) - 3
df_missing_tiles_reg_sp <- t(df_missing_tiles_reg[,1:n_tiles]) #transpose to get df with lines centroids of tiles
df_missing_tiles_reg_sp <- as.data.frame(df_missing_tiles_reg_sp)
names(df_missing_tiles_reg_sp) <- df_missing_tiles_reg$date
df_missing_tiles_reg_sp$tot_missing <- rowSums(df_missing_tiles_reg_sp) #total missing over a year by tile
df_missing_tiles_reg_sp$tot_pred <- length(df_missing_tiles_reg$date) - df_missing_tiles_reg_sp$tot_missing

tiles_names <- names(df_missing_tiles_reg)[1:n_tiles]
list_xy <- strsplit(unlist(basename(tiles_names)),"_")
sub("X","",list_xy)

#list_xy <- lapply(centroids_pts,function(x){coordinates(x)})
coord_xy <- do.call(rbind,list_xy)
y_val <- as.numeric(coord_xy[,1]) #lat, long! need to reverse
x_val <- as.numeric(coord_xy[,2])
coordinates(df_missing_tiles_sp) <- as.matrix(cbind(x_val,y_val))

#This contains in rows date and tiles columns
filename_df_missing_tiles_sp <- file.path(out_dir,paste0("df_missing_by_centroids_tiles_and_dates",region_name,"_",pred_mod_name,"_",out_suffix,".txt"))
write.table(as.data.frame(df_missing_tiles_sp),file=filename_df_missing_tiles_sp,sep=",")

#browser()

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
  png_filename_histogram <-  file.path(out_dir,paste("Figure_histogram_","region_missing_tiles","_",out_suffix,".png",sep =""))
  
  png(filename=png_filename_histogram,width = col_mfrow * res_pix,height = row_mfrow * res_pix)
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

## combine a make count+ plot




