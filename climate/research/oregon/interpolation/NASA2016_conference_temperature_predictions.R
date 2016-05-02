####################################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  NASA 2016 Meeting: biodiversity and ecological forecasting ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Combining tables and figures for individual runs for years and tiles.
#Figures and data for the AAG conference are also produced.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 05/01/2016  
#MODIFIED ON: 05/02/2016            
#Version: 1
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

###### Function used in the script #######
  
#script_path <- "/nobackupp8/bparmen1/env_layers_scripts" #path to script
script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts" #path to script

#Mosaic related
#script_path <- "/home/parmentier/Data/IPLANT_project/env_layers_scripts"
function_mosaicing_functions <- "global_run_scalingup_mosaicing_function_04232016.R" #PARAM12
function_mosaicing <-"global_run_scalingup_mosaicing_05012016.R"
source(file.path(script_path,function_mosaicing)) #source all functions used in this script 
source(file.path(script_path,function_mosaicing_functions)) #source all functions used in this script 

#Assessment
function_assessment_part1_functions <- "global_run_scalingup_assessment_part1_functions_02112015.R" #PARAM12
function_assessment_part1a <-"global_run_scalingup_assessment_part1a_01042016.R"
function_assessment_part2 <- "global_run_scalingup_assessment_part2_02092016.R"
function_assessment_part2_functions <- "global_run_scalingup_assessment_part2_functions_01032016.R"
function_assessment_part3 <- "global_run_scalingup_assessment_part3_04292016b.R"
source(file.path(script_path,function_assessment_part1_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part1a)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part2_functions)) #source all functions used in this script 
source(file.path(script_path,function_assessment_part3)) #source all functions used in this script 

### Parameters, constants and arguments ###

CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #constant 1

var<-"TMAX" # variable being interpolated #param 1, arg 1

##Add for precip later...
if (var == "TMAX") {
  y_var_name <- "dailyTmax"
  y_var_month <- "TMax"
}
if (var == "TMIN") {
  y_var_name <- "dailyTmin"
  y_var_month <- "TMin"
}


##Add for precip later...
if (var == "TMAX") {
  variable_name <- "maximum temperature"
}
if (var == "TMIN") {
  variable_name <- "minimum temperature"
}

#interpolation_method<-c("gam_fusion") #other otpions to be added later
interpolation_method<-c("gam_CAI") #param 2
CRS_interp <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #param 3
#CRS_interp <-"+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs";

out_region_name<-""
list_models<-c("y_var ~ s(lat,lon,k=5) + s(elev_s,k=3) + s(LST,k=3)") #param 4

#reg1 (North Am), reg2(Europe),reg3(Asia), reg4 (South Am), reg5 (Africa), reg6 (Australia-Asia)
#master directory containing the definition of tile size and tiles predicted
in_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/assessment"
in_dir_mosaic <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics"

region_name <- c("reg4") #param 6, arg 3

create_out_dir_param <- TRUE #param 9, arg 6
out_suffix <- "_meeting_NASA_reg4_04292016"

out_dir <- "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/assessment"

create_out_dir_param <- TRUE #param 9, arg 

#run_figure_by_year <- TRUE # param 10, arg 7
list_year_predicted <- "1984,2014"

file_format <- ".tif" #format for mosaiced files # param 11
NA_flag_val <- -32768  #No data value, # param 12
#-32768
#num_cores <- 6 #number of cores used # param 13, arg 8
plotting_figures <- TRUE #running part2 of assessment to generate figures... # param 14
#num_cores <- args[8] #number of cores used # param 13, arg 8
num_cores <- 11 #number of cores used # param 13, arg 8
#python_bin <- "/nobackupp6/aguzman4/climateLayers/sharedModules2/bin" #PARAM 30
python_bin <- "/usr/bin" #PARAM 30

day_start <- "19990101" #PARAM 12 arg 12
day_end <- "19990103" #PARAM 13 arg 13

#infile_mask <- "/nobackupp8/bparmen1/NEX_data/regions_input_files/r_mask_LST_reg4.tif"
infile_mask <- "/data/project/layers/commons/NEX_data/regions_input_files/r_mask_LST_reg4.tif"

#run_figure_by_year <- TRUE # param 10, arg 7
list_year_predicted <- "1984,2014"
scaling <- 100

##################### START SCRIPT #################

####### PART 1: Read in data ########
out_dir <- in_dir
if (create_out_dir_param == TRUE) {
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

setwd(out_dir)


###########  ####################

#start_date <- day_to_mosaic_range[1]
#end_date <- day_to_mosaic_range[2]
start_date <- day_start #PARAM 12 arg 12
end_date <- day_end #PARAM 13 arg 13

date_to_plot <- seq(as.Date(strptime(start_date,"%Y%m%d")), as.Date(strptime(end_date,"%Y%m%d")), 'day')
l_dates <- format(date_to_plot,"%Y%m%d") #format back to the relevant date format for files

raster_name_lf <- c("/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19990101_reg4_1999_m_gam_CAI_dailyTmax_19990101_reg4_1999.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19990102_reg4_1999_m_gam_CAI_dailyTmax_19990102_reg4_1999.tif",
                    "/data/project/layers/commons/NEX_data/climateLayers/out/reg4/mosaic/int_mosaics/comp_r_m_use_edge_weights_weighted_mean_gam_CAI_dailyTmax_19990103_reg4_1999_m_gam_CAI_dailyTmax_19990103_reg4_1999.tif")

pred_temp_s <- raster(raster_name_lf[1])

lf_files <- c(raster_name_in) #match to mask
rast_ref <- infile_mask
##Maching resolution is probably only necessary for the r mosaic function
#Modify later to take into account option R or python...
list_param_raster_match <- list(lf_files,rast_ref,file_format,python_bin,out_suffix,out_dir)
names(list_param_raster_match) <- c("lf_files","rast_ref","file_format","python_bin","out_suffix","out_dir_str")
r_pred_matched <- raster_match(1,list_param_raster_match)
raster_name_in <- c(r_pred_matched)

if(mask_pred==TRUE){
  
  r_mask <- raster(infile_mask)
  extension_str <- extension(raster_name_in)
  raster_name_tmp <- gsub(extension_str,"",basename(raster_name_in))
  raster_name <- file.path(out_dir,paste(raster_name_tmp,"_masked.tif",sep = ""))
  r_pred <- mask(raster(raster_name_in),r_mask,filename = raster_name,overwrite = TRUE)
}

extension_str <- extension(filename(r_pred))
raster_name_tmp <- gsub(extension_str,"",basename(filename(r_pred)))
raster_name <- file.path(out_dir,paste(raster_name_tmp,"_rescaled.tif",sep = ""))

r_pred <- overlay(r_pred, fun=function(x){x*1/scaling},filename=raster_name)
NAvalue(r_pred) <- NA_flag_val
r_pred <- setMinMax(r_pred)

month_name <- month.name()
l_dates <- as.Date(strptime(date_proc,"%Y%m%d"))

#s.range <- c(min(minValue(pred_temp_s)), max(maxValue(pred_temp_s)))
#s.range <- s.range+c(5,-5)
#col.breaks <- pretty(s.range, n=200)
#lab.breaks <- pretty(s.range, n=100)
temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
max_val<-s.range[2]
min_val <-s.range[1]
#max_val<- -10
min_val <- 0

#layout_m<-c(1,3) #one row two columns
date_proc <- l_dates[i]

#png(paste("Figure7a__spatial_pattern_tmax_prediction_levelplot_",date_selected,out_prefix,".png", sep=""),
#    height=480*layout_m[1],width=480*layout_m[2])
#plot(r_pred,col=temp.colors(255),zlim=c(-3500,4500))
plot(r_pred,col=matlab.like(255),zlim=c(-35,45))
#paste(raster_name[1:7],collapse="_")
#add filename option later

res_pix <- 1200
#res_pix <- 480

col_mfrow <- 1
row_mfrow <- 1
date_proc <- 
png_filename <-  file.path(out_dir_str,paste("Figure4_clim_mosaics_day_","_",date_proc,"_",reg_name,"_",out_suffix,".png",sep =""))

png(filename=png_filename,width = col_mfrow * res_pix,height = row_mfrow * res_pix)

plot(r_pred,main = paste("Predicted ",variable_name, " on ",date_proc , " ", ,sep = ""),cex.main =1.5)

dev.off()

#col.regions=temp.colors(25))
#dev.off()



############################ END OF SCRIPT ##################################