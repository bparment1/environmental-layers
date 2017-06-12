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

## selecte relevant files for region
in_dir_list_tmp <- list.files(pattern =paste0(".*.",region_name,".*."),full.names=T,path=in_dir)

list_lf_df_missing_tiles <- unlist(mclapply(in_dir_list_tmp,
                                     FUN=function(x){list.files(path=x,
                                                                pattern=paste0("df_missing_by_dates_tiles_predicted_tif_",".*.",region_name,".*.txt"),full.names=T)},
                                     mc.preschedule=FALSE,mc.cores = num_cores))
list_df_missing_tiles <- mclapply(list_lf_df_missing_tiles,
                                     FUN=function(x){read.table(x,sep=",",stringsAsFactors = F)},
                                     mc.preschedule=FALSE,mc.cores = num_cores)

df_missing_tiles_reg <- do.call(rbind,list_df_missing_tiles)
sum(df_missing_tiles_reg)

range(df_missing_tiles_reg$tot_missing)
table(df_missing_tiles_reg$tot_missing)
histogram(df_missing_tiles_reg$tot_missing)



## combine a make count+ plot




