  ##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part2 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Accuracy methods are added in the the function scripts to evaluate results.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 03/23/2014  
#MODIFIED ON: 09/21/2015            
#Version: 4
#PROJECT: Environmental Layers project     
#COMMENTS: analyses for run 10 global analyses,all regions 1500x4500km with additional tiles to increase overlap 
#TODO:
#1) Split functions and master script
#2) Make this is a script/function callable from the shell/bash
#3) Check image format for tif

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

#### FUNCTION USED IN SCRIPT

function_analyses_paper1 <-"contribution_of_covariates_paper_interpolation_functions_07182014.R" #first interp paper
function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_08132014.R"

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
  #create if does not exists
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  return(out_dir)
}

 #Remove models that were not fitted from the list
#All modesl that are "try-error" are removed
remove_errors_list<-function(list_items){
  
  #This function removes "error" items in a list
  list_tmp<-list_items
    if(is.null(names(list_tmp))){
    names(list_tmp) <- paste("l",1:length(list_tmp),sep="_")
    names(list_items) <- paste("l",1:length(list_tmp),sep="_")
  }

  for(i in 1:length(list_items)){
    if(inherits(list_items[[i]],"try-error")){
      list_tmp[[i]]<-0
    }else{
    list_tmp[[i]]<-1
   }
  }
  cnames<-names(list_tmp[list_tmp>0])
  x <- list_items[match(cnames,names(list_items))]
  return(x)
}

#turn term from list into data.frame
#name_col<-function(i,x){
#x_mat<-x[[i]]
#x_df<-as.data.frame(x_mat)
#x_df$mod_name<-rep(names(x)[i],nrow(x_df))
#x_df$term_name <-row.names(x_df)
#return(x_df)
#}
#Function to rasterize a table with coordinates and variables...,maybe add option for ref image??
rasterize_df_fun <- function(data_tb,coord_names,proj_str,out_suffix,out_dir=".",file_format=".rst",NA_flag_val=-9999,tolerance_val= 0.000120005){
  data_spdf <- data_tb
  coordinates(data_spdf) <- cbind(data_spdf[[coord_names[1]]],data_spdf[[coord_names[2]]])
  proj4string(data_spdf) <- proj_str

  data_pix <- try(as(data_spdf,"SpatialPixelsDataFrame"))
  #tolerance_val <- 0.000120005 
  #tolerance_val <- 0.000856898
  if(inherits(data_pix,"try-error")){
      data_pix <- SpatialPixelsDataFrame(data_spdf, data=data_spdf@data, tolerance=tolerance_val) 
  }
  
  #test <- as(data_spdf,"SpatialPixelsDataFrame")

  # set up an 'empty' raster, here via an extent object derived from your data
  #e <- extent(s100[,1:2])
  #e <- e + 1000 # add this as all y's are the same

  #r <- raster(e, ncol=10, nrow=2)
  # or r <- raster(xmn=, xmx=,  ...

  data_grid <- as(data_pix,"SpatialGridDataFrame") #making it a regural grid
  r_ref <- raster(data_grid) #this is the ref image
  rast_list <- vector("list",length=ncol(data_tb))
  
  for(i in 1:(ncol(data_tb))){
    field_name <- names(data_tb)[i]
    var <- as.numeric(data_spdf[[field_name]])
    data_spdf$var  <- var
    #r <-rasterize(data_spdf,r_ref,field_name)
    r <-rasterize(data_spdf,r_ref,"var",NAflag=NA_flag_val,fun=mean) #prolem with NA in NDVI!!

    data_name<-paste("r_",field_name,sep="") #can add more later...
    #raster_name<-paste(data_name,out_names[j],".tif", sep="")
    raster_name<-paste(data_name,out_suffix,file_format, sep="")
  
    writeRaster(r, NAflag=NA_flag_val,filename=file.path(out_dir,raster_name),overwrite=TRUE)
    #Writing the data in a raster file format...
    rast_list[i] <-file.path(out_dir,raster_name)
  }
  return(unlist(rast_list))
}

plot_raster_tb_diagnostic <- function(reg_layer,tb_dat,df_tile_processed,date_selected,mod_selected,var_selected,out_suffix){
  
  test <- subset(tb_dat,pred_mod==mod_selected & date==date_selected,select=c("tile_id",var_selected))

  test_data_tb <- merge(df_tile_processed,test,by="tile_id",all=T) #keep all
  test_r <- subset(test_data_tb,select=c("lat","lon","tile_id",var_selected))
  out_suffix_str <- paste(var_selected,mod_selected,date_selected,out_suffix,sep="_")
  coord_names <- c("lon","lat")
  l_rast <- rasterize_df_fun(test_r,coord_names,proj_str,out_suffix_str,out_dir=".",file_format=".tif",NA_flag_val=-9999,tolerance_val=0.000120005)
  #mod_kr_stack <- stack(mod_kr_l_rast)
  d_tb_rast <- stack(l_rast)
  names(d_tb_rast) <- names(test_r)
  #plot(d_tb_rast)
  r <- subset(d_tb_rast,"rmse")
  names(r) <- paste(mod_selected,var_selected,date_selected,sep="_")
  #plot info: with labels

  res_pix <- 1200
  col_mfrow <- 1
  row_mfrow <- 1

  png(filename=paste("Figure9_",names(r),"_map_processed_region_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  #plot(reg_layer)
  #p1 <- spplot(reg_layer,"ISO",colorkey=FALSE) #Use ISO instead of NAME_1 to see no color?
  title_str <- paste(names(r),"for ", region_name,sep="")

  p0 <- levelplot(r,col.regions=matlab.like(25),margin=F,main=title_str)
  p_shp <- layer(sp.polygons(reg_layer, lwd=1, col='black'))

  p <- p0 + p_shp
  print(p)

  dev.off()

}

create_raster_from_tb_diagnostic <- function(i,list_param){
  #create a raster image using tile centroids and given fields  from tb diagnostic data
  tb_dat <- list_param$tb_dat
  df_tile_processed <- list_param$df_tile_processed
  date_selected <- list_param$date_selected[i]
  mod_selected <- list_param$mod_selected
  var_selected <- list_param$var_selected
  out_suffix <- list_param$out_suffix
  
  test <- subset(tb_dat,pred_mod==mod_selected & date==date_selected,select=c("tile_id",var_selected))

  test_data_tb <- merge(df_tile_processed,test,by="tile_id",all=T) #keep all
  test_r <- subset(test_data_tb,select=c("lat","lon","tile_id",var_selected))
  out_suffix_str <- paste(var_selected,mod_selected,date_selected,out_suffix,sep="_")
  coord_names <- c("lon","lat")
  l_rast <- rasterize_df_fun(test_r,coord_names,proj_str,out_suffix_str,out_dir=".",file_format,NA_flag_val,tolerance_val=0.000120005)
  #mod_kr_stack <- stack(mod_kr_l_rast)
  #d_tb_rast <- stack(l_rast)
  #r <- subset(d_tb_rast,var_selected)
  #names(d_tb_rast) <- names(test_r)
  return(l_rast[4])
}

assign_FID_spatial_polygons_df <-function(list_spdf,ID_str=NULL){
  list_spdf_tmp <- vector("list",length(list_spdf))
  if(is.null(ID_str)){
    nf <- 0 #number of features
    #for(i in 1:length(spdf)){
    #    shp1 <- list_spdf[[i]]
    #    f <- nrow(shp1)
    #    nf <- nf + f
    #}
    #This assumes that the list has one feature per item list
    nf <- length(list_spdf)
    ID_str <- as.character(1:nf)
  }
  for(i in 1:length(list_spdf)){
    #test=spRbind(shps_tiles[[1]],shps_tiles[[2]])
    shp1 <- list_spdf[[i]]
    shp1$FID <- ID_str
    shp1<- spChFIDs(shp1, as.character(shp1$FID)) #assign ID
    list_spdf_tmp[[i]]  <-shp1
  }
  return(list_spdf_tmp)
}

combine_spatial_polygons_df_fun <- function(list_spdf_tmp,ID_str=NULL){
  if(is.null(ID_str)){
    #call function
    list_spdf_tmp <- assign_FID_spatial_polygons_df
  }
  combined_spdf <- list_spdf_tmp[[1]]
  for(i in 2:length(list_spdf_tmp)){
    combined_spdf <- rbind(combined_spdf,list_spdf_tmp[[i]])
    #sapply(slot(shps_tiles[[2]], "polygons"), function(x) slot(x, "ID"))
    #rownames(as(alaska.tract, "data.frame"))
  }
  return(combined_spdf)
}

plot_daily_mosaics <- function(i,list_param_plot_daily_mosaics){
  #Purpose:
  #This functions mask mosaics files for a default range (-100,100 deg).
  #It produces a masked tif in a given dataType format (FLT4S)
  #It creates a figure of mosaiced region being interpolated.
  #Author: Benoit Parmentier
  #Parameters:
  #lf_m: list of files 
  #reg_name:region name with tile size included
  #To do:
  #Add filenames
  #Add range
  #Add output dir
  #Add dataType_val option
  
  ##### BEGIN ########
  
  #Parse the list of parameters
  lf_m <- list_param_plot_daily_mosaics$lf_m
  reg_name <- list_param_plot_daily_mosaics$reg_name
  out_dir_str <- list_param_plot_daily_mosaics$out_dir_str
  out_suffix <- list_param_plot_daily_mosaics$out_suffix
  l_dates <- list_param_plot_daily_mosaics$l_dates


  #list_param_plot_daily_mosaics <- list(lf_m=lf_m,reg_name=reg_name,out_dir_str=out_dir_str)

  
  #rast_list <- vector("list",length=length(lf_m))
  r_test<- raster(lf_m[i])

  m <- c(-Inf, -100, NA,  
         -100, 100, 1, 
         100, Inf, NA) #can change the thresholds later
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  rc <- reclassify(r_test, rclmat,filename=paste("rc_tmp_",i,".tif",sep=""),dataType="FLT4S",overwrite=T)
  file_name <- unlist(strsplit(basename(lf_m[i]),"_"))
  
  #date_proc <- file_name[7] #specific tot he current naming of files
  date_proc <- l_dates[i]
  #paste(raster_name[1:7],collapse="_")
  #add filename option later
  extension_str <- extension(filename(r_test))
  raster_name_tmp <- gsub(extension_str,"",basename(filename(r_test)))
  raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_masked.tif",sep=""))
  r_pred <- mask(r_test,rc,filename=raster_name,overwrite=TRUE)
  
  res_pix <- 1200
  #res_pix <- 480

  col_mfrow <- 1
  row_mfrow <- 1
  
  png(filename=paste("Figure9_clim_mosaics_day_test","_",date_proc,"_",reg_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_pred,main=paste("Predicted on ",date_proc , " ", reg_name,sep=""),cex.main=1.5)
  dev.off()
  
  return(raster_name)
  
}

plot_screen_raster_val<-function(i,list_param){
  ##USAGE###
  #Screen raster list and produced plot as png.
  fname <-list_param$lf_raster_fname[i]
  screenRast <- list_param$screenRast
  l_dates<- list_param$l_dates
  out_dir_str <- list_param$out_dir_str
  prefix_str <-list_param$prefix_str
  out_suffix_str <- list_param$out_suffix_str
  
  ### START SCRIPT ####
  date_proc <- l_dates[i]
  
  if(screenRast==TRUE){
    r_test <- raster(fname)

    m <- c(-Inf, -100, NA,  
         -100, 100, 1, 
         100, Inf, NA) #can change the thresholds later
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    rc <- reclassify(r_test, rclmat,filename=paste("rc_tmp_",i,".tif",sep=""),dataType="FLT4S",overwrite=T)
    #file_name <- unlist(strsplit(basename(lf_m[i]),"_"))
    extension_str <- extension(filename(r_test))
    raster_name_tmp <- gsub(extension_str,"",basename(filename(r_test)))
    raster_name <- file.path(out_dir_str,paste(raster_name_tmp,"_masked.tif",sep=""))
  
    r_pred <- mask(r_test,rc,filename=raster_name,overwrite=TRUE)
  }else{
    r_pred <- raster(fname)
  }
  
  #date_proc <- file_name[7] #specific tot he current naming of files
  #date_proc<- "2010010101"

  #paste(raster_name[1:7],collapse="_")
  #add filename option later

  res_pix <- 960
  #res_pix <- 480

  col_mfrow <- 1
  row_mfrow <- 1
  
#  png(filename=paste("Figure10_clim_world_mosaics_day_","_",date_proc,"_",tile_size,"_",out_suffix,".png",sep=""),
#    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  png(filename=paste(prefix_str,"_",date_proc,"_",tile_size,"_",out_suffix_str,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_pred,main=paste("Predicted on ",date_proc , " ", tile_size,sep=""),cex.main=1.5)
  dev.off()

}

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
#in_dir1 <- "/nobackupp4/aguzman4/climateLayers/output4" #On NEX
#parent output dir for the current script analyes
#out_dir <- "/nobackup/bparmen1/" #on NEX
#in_dir_shp <- "/nobackupp4/aguzman4/climateLayers/output4/subset/shapefiles/"

y_var_name <- "dailyTmax" #PARAM1
interpolation_method <- c("gam_CAI") #PARAM2
#out_suffix<-"run10_global_analyses_01282015"
#out_suffix <- "output_run10_1000x3000_global_analyses_02102015"
out_suffix <- "run10_1500x4500_global_analyses_pred_1982_09152015" #PARAM3
out_dir <- "/data/project/layers/commons/NEX_data/output_run10_1500x4500_global_analyses_pred_1982_09152015" #PARAM4
create_out_dir_param <- FALSE #PARAM 5

mosaic_plot <- FALSE #PARAM6

#if daily mosaics NULL then mosaicas all days of the year

day_to_mosaic <- c("19820101","19820102","19820103","19820104","19820105",
                   "19820106","19820107","19820108","19820109","19820110",
                   "1982011")

CRS_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84 #CONSTANT1
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

proj_str<- CRS_WGS84 #PARAM 8 #check this parameter
file_format <- ".rst" #PARAM 9
NA_value <- -9999 #PARAM10
NA_flag_val <- NA_value
                                   
tile_size <- "1500x4500" #PARAM 11
multiple_region <- TRUE #PARAM 12

region_name <- "world" #PARAM 13
plot_region <- TRUE
num_cores <- 10 #PARAM 14
reg_modified <- TRUE
region <- c("reg4") #reference region to merge if necessary #PARAM 16

########################## START SCRIPT ##############################


####### PART 1: Read in data ########

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

setwd(out_dir)

###Table 1: Average accuracy metrics
###Table 2: daily accuracy metrics for all tiles

summary_metrics_v <- read.table(file=file.path(out_dir,paste("summary_metrics_v2_NA_",out_suffix,".txt",sep="")),sep=",")
#fname <- file.path(out_dir,paste("summary_metrics_v2_NA_",out_suffix,".txt",sep=""))
tb <- read.table(file=file.path(out_dir,paste("tb_diagnostic_v_NA","_",out_suffix,".txt",sep="")),sep=",")
#tb_diagnostic_s_NA_run10_global_analyses_11302014.txt
tb_s <- read.table(file=file.path(out_dir,paste("tb_diagnostic_s_NA","_",out_suffix,".txt",sep="")),sep=",")

tb_month_s <- read.table(file=file.path(out_dir,paste("tb_month_diagnostic_s_NA","_",out_suffix,".txt",sep="")),sep=",")
pred_data_month_info <- read.table(file=file.path(out_dir,paste("pred_data_month_info_",out_suffix,".txt",sep="")),sep=",")
pred_data_day_info <- read.table(file=file.path(out_dir,paste("pred_data_day_info_",out_suffix,".txt",sep="")),sep=",")
df_tile_processed <- read.table(file=file.path(out_dir,paste("df_tile_processed_",out_suffix,".txt",sep="")),sep=",")

#add column for tile size later on!!!

tb$pred_mod <- as.character(tb$pred_mod)
summary_metrics_v$pred_mod <- as.character(summary_metrics_v$pred_mod)
summary_metrics_v$tile_id <- as.character(summary_metrics_v$tile_id)
df_tile_processed$tile_id <- as.character(df_tile_processed$tile_id)

tb_month_s$pred_mod <- as.character(tb_month_s$pred_mod)
tb_month_s$tile_id<- as.character(tb_month_s$tile_id)
tb_s$pred_mod <- as.character(tb_s$pred_mod)
tb_s$tile_id <- as.character(tb_s$tile_id)


#multiple regions?
if(multiple_region==TRUE){
  df_tile_processed$reg <- basename(dirname(as.character(df_tile_processed$path_NEX)))
  
  tb <- merge(tb,df_tile_processed,by="tile_id")
  tb_s <- merge(tb_s,df_tile_processed,by="tile_id")
  tb_month_s<- merge(tb_month_s,df_tile_processed,by="tile_id")
  summary_metrics_v <- merge(summary_metrics_v,df_tile_processed,by="tile_id")

}

tb_all <- tb

summary_metrics_v_all <- summary_metrics_v 
#deal with additional tiles...
# if(reg_modified==T){
#   
#   summary_metrics_v_tmp <- summary_metrics_v
#   #summary_metrics_v_tmp$reg[summary_metrics_v_tmp$reg=="reg_1b"] <- "reg1"
#   #summary_metrics_v_tmp$reg[summary_metrics_v_tmp$reg=="reg_1c"] <- "reg1"
#   #summary_metrics_v_tmp$reg[summary_metrics_v_tmp$reg=="reg_3b"] <- "reg3"
#   summary_metrics_v_tmp$reg[summary_metrics_v_tmp$reg=="reg5b"] <- "reg5"
# 
#   summary_metrics_v_tmp$reg_all <- summary_metrics_v$reg
#   ###
#   summary_metrics_v<- summary_metrics_v_tmp
#   
#   ###
#   tb_tmp <- tb
#   #tb_tmp$reg[tb_tmp$reg=="reg_1b"] <- "reg1"
#   #tb_tmp$reg[tb_tmp$reg=="reg_1c"] <- "reg1"
#   #tb_tmp$reg[tb_tmp$reg=="reg_3b"] <- "reg3"
#   tb_tmp$reg[tb_tmp$reg=="reg5b"] <- "reg5"
# 
#   ###
#   tb <- tb_tmp
# }

table(summary_metrics_v_all$reg)
table(summary_metrics_v$reg)
table(tb_all$reg)
table(tb$reg)

############ PART 2: PRODUCE FIGURES ################

###########################
### Figure 1: plot location of the study area with tiles processed

#df_tiled_processed <- na.omit(df_tile_processed) #remove other list of folders irrelevant
#list_shp_reg_files <- df_tiled_processed$shp_files
list_shp_reg_files<- as.character(df_tile_processed$shp_files)
#list_shp_reg_files <- file.path("/data/project/layers/commons/NEX_data/",out_dir,
#          as.character(df_tile_processed$tile_coord),"shapefiles",basename(list_shp_reg_files))
list_shp_reg_files <- file.path("/data/project/layers/commons/NEX_data/",out_dir,
          "shapefiles",basename(list_shp_reg_files))

#table(summary_metrics_v$reg)
#table(summary_metrics_v$reg)

### Potential function starts here:
#function(in_dir,out_dir,list_shp_reg_files,title_str,region_name,num_cores,out_suffix,out_suffix)

### First get background map to display where study area is located
#can make this more general later on..should have this already in a local directory on Atlas or NEX!!!!

if (region_name=="USA"){
  usa_map <- getData('GADM', country='USA', level=1) #Get US map
  #usa_map <- getData('GADM', country=region_name,level=1) #Get US map, this is not working right now
  usa_map <- usa_map[usa_map$NAME_1!="Alaska",] #remove Alaska
  reg_layer <- usa_map[usa_map$NAME_1!="Hawaii",] #remove Hawai 
}

if (region_name=="world"){
  #http://www.diva-gis.org/Data
  countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp"
  path_to_shp <- dirname(countries_shp)
  layer_name <- sub(".shp","",basename(countries_shp))
  reg_layer <- readOGR(path_to_shp, layer_name)
  #proj4string(reg_layer) <- CRS_locs_WGS84
  #reg_shp<-readOGR(dirname(list_shp_reg_files[[i]]),sub(".shp","",basename(list_shp_reg_files[[i]])))
}

centroids_pts <- vector("list",length(list_shp_reg_files))
shps_tiles <- vector("list",length(list_shp_reg_files))
#collect info: read in all shapfiles
#This is slow...make a function and use mclapply??
#/data/project/layers/commons/NEX_data/output_run6_global_analyses_09162014/shapefiles
  
for(i in 1:length(list_shp_reg_files)){
  #path_to_shp <- dirname(list_shp_reg_files[[i]])
  path_to_shp <- file.path(out_dir,"/shapefiles")
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
  if(!inherits(shp1,"try-error")){
    plot(shp1,add=T,border="blue")
    #plot(pt,add=T,cex=2,pch=5)
    label_id <- df_tile_processed$tile_id[i]
    text(coordinates(pt)[1],coordinates(pt)[2],labels=i,cex=1.3,font=2,col=c("red"))
  }
}
#title(paste("Tiles ", tile_size,region_name,sep=""))

dev.off()
      
#unique(summaty_metrics$tile_id)
#text(lat-shp,)
#union(list_shp_reg_files[[1]],list_shp_reg_files[[2]])

###############
### Figure 2: boxplot of average accuracy by model and by tiles


tb$tile_id <- factor(tb$tile_id, levels=unique(tb$tile_id))

model_name <- as.character(unique(tb$pred_mod))


## Figure 2a

for(i in  1:length(model_name)){
  
  res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1

  png(filename=paste("Figure2a_boxplot_with_oultiers_by_tiles_",model_name[i],"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  boxplot(rmse~tile_id,data=subset(tb,tb$pred_mod==model_name[i]))
  title(paste("RMSE per ",model_name[i]))
  
  dev.off()
}

## Figure 2b
#wtih ylim and removing trailing...
for(i in  1:length(model_name)){

  res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1
  png(filename=paste("Figure2b_boxplot_scaling_by_tiles","_",model_name[i],"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  model_name <- unique(tb$pred_mod)
  boxplot(rmse~tile_id,data=subset(tb,tb$pred_mod==model_name[i])
          ,ylim=c(0,4),outline=FALSE)
  title(paste("RMSE per ",model_name[i]))
  dev.off()
}
#bwplot(rmse~tile_id, data=subset(tb,tb$pred_mod=="mod1"))

###############
### Figure 3: boxplot of average RMSE by model acrosss all tiles

## Figure 3a
res_pix <- 480
col_mfrow <- 1
row_mfrow <- 1

png(filename=paste("Figure3a_boxplot_overall_region_with_oultiers_",model_name[i],"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

boxplot(rmse~pred_mod,data=tb)#,names=tb$pred_mod)
title("RMSE per model over all tiles")

dev.off()

## Figure 3b
png(filename=paste("Figure3b_boxplot_overall_region_scaling_",model_name[i],"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

boxplot(rmse~pred_mod,data=tb,ylim=c(0,5),outline=FALSE)#,names=tb$pred_mod)
title("RMSE per model over all tiles")

dev.off()


################ 
### Figure 4: plot predicted tiff for specific date per model

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

  png(filename=paste("Figure5_ac_metrics_ranked_",model_name[i],"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  x<- as.character(df_ac_mod$tile_id)
  barplot(df_ac_mod$rmse, names.arg=x)
  #plot(ac_mod1,cex=sqrt(ac_mod1$rmse),pch=1,add=T)
  #plot(ac_mod1,cex=(ac_mod1$rmse1)*2,pch=1,add=T)
  title(paste("RMSE ranked by tile for ",model_name[i],sep=""))

  dev.off()
  
}

######################
### Figure 6: plot map of average RMSE per tile at centroids

#Turn summary table to a point shp

# coordinates(summary_metrics_v) <- cbind(summary_metrics_v$lon,summary_metrics_v$lat)
# proj4string(summary_metrics_v) <- CRS_WGS84
# #lf_list <- lf_pred_list
# list_df_ac_mod <- vector("list",length=length(lf_pred_list))
# for (i in 1:length(lf_list)){
#   
#   ac_mod <- summary_metrics_v[summary_metrics_v$pred_mod==model_name[i],]
#   r_pred <- raster(lf_list[i])
#   
#   res_pix <- 480
#   col_mfrow <- 1
#   row_mfrow <- 1
# 
#   png(filename=paste("Figure6_ac_metrics_map_centroids_tile_",model_name[i],"_",out_suffix,".png",sep=""),
#     width=col_mfrow*res_pix,height=row_mfrow*res_pix)
# 
#   plot(r_pred)  
#   
#   #plot(ac_mod1,cex=sqrt(ac_mod1$rmse),pch=1,add=T)
#   plot(ac_mod,cex=(ac_mod$rmse^2)/10,pch=1,add=T)
#   #plot(ac_mod1,cex=(ac_mod1$rmse1)*2,pch=1,add=T)
#   title(paste("Averrage RMSE per tile and by ",model_name[i]))
# 
#   dev.off()
#   
#   ### Ranking by tile...
#   #df_ac_mod <- 
#   list_df_ac_mod[[i]] <- arrange(as.data.frame(ac_mod),desc(rmse))[,c("rmse","mae","tile_id")]
# }
  
#quick kriging...
#autokrige(rmse~1,r2,)


### Without 

#list_df_ac_mod <- vector("list",length=length(lf_pred_list))
list_df_ac_mod <- vector("list",length=2)

for (i in 1:length(model_name)){
  
  ac_mod <- summary_metrics_v[summary_metrics_v$pred_mod==model_name[i],]
  #r_pred <- raster(lf_list[i])
  
  res_pix <- 1200
  #res_pix <- 480

  col_mfrow <- 1
  row_mfrow <- 1

  png(filename=paste("Figure6_ac_metrics_map_centroids_tile_",model_name[i],"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  #plot(r_pred)  
  #plot(reg_layer)
  #plot(ac_mod1,cex=sqrt(ac_mod1$rmse),pch=1,add=T)
  #plot(ac_mod,cex=(ac_mod$rmse^2)/10,pch=1,col="red",add=T)

  #coordinates(ac_mod) <- ac_mod[,c("lon","lat")] 
  coordinates(ac_mod) <- ac_mod[,c("lon.x","lat.x")] #solve this later
  p_shp <- layer(sp.polygons(reg_layer, lwd=1, col='black'))
  #title("(a) Mean for 1 January")
  p <- bubble(ac_mod,"rmse",main=paste("Averrage RMSE per tile and by ",model_name[i]))
  p1 <- p+p_shp
  print(p1)
  #plot(ac_mod1,cex=(ac_mod1$rmse1)*2,pch=1,add=T)
  #title(paste("Averrage RMSE per tile and by ",model_name[i]))

  dev.off()
  
  ### Ranking by tile...
  #df_ac_mod <- 
  list_df_ac_mod[[i]] <- arrange(as.data.frame(ac_mod),desc(rmse))[,c("rmse","mae","tile_id")]
}

## Number of tiles with information:
sum(df_tile_processed$metrics_v) #26,number of tiles with raster object
length(df_tile_processed$metrics_v) #26,number of tiles in the region
sum(df_tile_processed$metrics_v)/length(df_tile_processed$metrics_v) #80 of tiles with info

#coordinates
#coordinates(summary_metrics_v) <- c("lon","lat")
coordinates(summary_metrics_v) <- c("lon.y","lat.y")

threshold_missing_day <- c(367,365,300,200)

nb<-nrow(subset(summary_metrics_v,model_name=="mod1"))  
sum(subset(summary_metrics_v,model_name=="mod1")$n_missing)/nb #33/35

## Make this a figure...

#plot(summary_metrics_v)
#Make this a function later so that we can explore many thresholds...

j<-1 #for model name 1
for(i in 1:length(threshold_missing_day)){
  
  #summary_metrics_v$n_missing <- summary_metrics_v$n == 365
  #summary_metrics_v$n_missing <- summary_metrics_v$n < 365
  summary_metrics_v$n_missing <- summary_metrics_v$n < threshold_missing_day[i]
  summary_metrics_v_subset <- subset(summary_metrics_v,model_name=="mod1")

  #res_pix <- 1200
  res_pix <- 960

  col_mfrow <- 1
  row_mfrow <- 1

  png(filename=paste("Figure7a_ac_metrics_map_centroids_tile_",model_name[j],"_","missing_day_",threshold_missing_day[i],
                     "_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  model_name[j]
  
  p_shp <- layer(sp.polygons(reg_layer, lwd=1, col='black'))
  #title("(a) Mean for 1 January")
  p <- bubble(summary_metrics_v_subset,"n_missing",main=paste("Missing per tile and by ",model_name[j]," for ",
                                                              threshold_missing_day[i]))
  p1 <- p+p_shp
  print(p1)
  #plot(ac_mod1,cex=(ac_mod1$rmse1)*2,pch=1,add=T)
  #title(paste("Averrage RMSE per tile and by ",model_name[i]))

  dev.off()
}

######################
### Figure 7: Number of predictions: daily and monthly

#xyplot(rmse~pred_mod | tile_id,data=subset(as.data.frame(summary_metrics_v),
#                                           pred_mod!="mod_kr"),type="b")

#xyplot(n~pred_mod | tile_id,data=subset(as.data.frame(summary_metrics_v),
#                                           pred_mod!="mod_kr"),type="h")

#cor

# 
## Figure 7a
png(filename=paste("Figure7a_number_daily_predictions_per_models","_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

xyplot(n~pred_mod | tile_id,data=subset(as.data.frame(summary_metrics_v),
                                           pred_mod!="mod_kr"),type="h")
dev.off()

table(tb$pred_mod)
table(tb$index_d)
table(subset(tb,pred_mod!="mod_kr"))
table(subset(tb,pred_mod=="mod1")$index_d)
#aggregate()
tb$predicted <- 1
test <- aggregate(predicted~pred_mod+tile_id,data=tb,sum)
xyplot(predicted~pred_mod | tile_id,data=subset(as.data.frame(test),
                                           pred_mod!="mod_kr"),type="h")

test

as.character(unique(test$tile_id)) #141 tiles

dim(subset(test,test$predicted==365 & test$pred_mod=="mod1"))
histogram(subset(test, test$pred_mod=="mod1")$predicted)
unique(subset(test, test$pred_mod=="mod1")$predicted)
table((subset(test, test$pred_mod=="mod1")$predicted))

#LST_avgm_min <- aggregate(LST~month,data=data_month_all,min)
histogram(test$predicted~test$tile_id)
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

## Figure 8a
res_pix <- 480
col_mfrow <- 1
row_mfrow <- 1

png(filename=paste("Figure8a_boxplot_overall_separated_by_region_with_oultiers_",model_name[i],"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

p<- bwplot(rmse~pred_mod | reg, data=tb,
           main="RMSE per model and region over all tiles")
print(p)
dev.off()

## Figure 8b
png(filename=paste("Figure8b_boxplot_overall_separated_by_region_scaling_",model_name[i],"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

boxplot(rmse~pred_mod,data=tb,ylim=c(0,5),outline=FALSE)#,names=tb$pred_mod)
title("RMSE per model over all tiles")
p<- bwplot(rmse~pred_mod | reg, data=tb,ylim=c(0,5),
           main="RMSE per model and region over all tiles")
print(p)
dev.off()

#####################################################
#### Figure 9: plotting mosaics of regions ###########
## plot mosaics...

#First collect file names


#names(lf_mosaics_reg) <- l_reg_name

#This part should be automated...
#plot Australia
#lf_m <- lf_mosaics_reg[[2]]
#out_dir_str <- out_dir
#reg_name <- "reg6_1000x3000"
#lapply()
#list_param_plot_daily_mosaics <- list(lf_m=lf_m,reg_name=reg_name,out_dir_str=out_dir_str,out_suffix=out_suffix)
#list_param_plot_daily_mosaics <- list(lf_m=lf_m,reg_name=reg_name,out_dir_str=out_dir_str,out_suffix=out_suffix,l_dates=day_to_mosaic)

#lf_m_mask_reg4_1500x4500 <- mclapply(1:2,FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics,mc.preschedule=FALSE,mc.cores = 6)
#debug(plot_daily_mosaics)
#lf_m_mask_reg6_1000x3000 <- plot_daily_mosaics(1,list_param=list_param_plot_daily_mosaics)

#lf_m_mask_reg6_1000x3000 <- mclapply(1:length(lf_m),FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics,mc.preschedule=FALSE,mc.cores = 10)

plot_screen_raster_val<-function(i,list_param){
  ##USAGE###
  #Screen raster list and produced plot as png.
  fname <-list_param$lf_raster_fname[i]
  screenRast <- list_param$screenRast
  l_dates<- list_param$l_dates
  out_dir_str <- list_param$out_dir_str
  prefix_str <-list_param$prefix_str
  out_suffix_str <- list_param$out_suffix_str
  
if(is.null(day_to_mosaic)){
  
  #idx <- seq(as.Date('2010-01-01'), as.Date('2010-12-31'), 'day')
  #idx <- seq(as.Date('20100101'), as.Date('20101231'), 'day')
  #date_l <- strptime(idx[1], "%Y%m%d") # interpolation date being processed
  #dates_l <- format(idx, "%Y%m%d") # interpolation date being processed
  day_to_mosaic <- dates_predicted #should be 365 days...
  #l_dates <- day_to_mosaic
}
  
if(plot_region==TRUE){
  
  #get the files
  l_reg_name <- unique(df_tile_processed$reg)
  l_reg_name <- c("reg4") #use this for the time
  #lf_mosaics_reg5 <- mixedsort(list.files(path="/data/project/layers/commons/NEX_data/output_run10_global_analyses_11302014/mosaics/reg5",
  #           pattern="CAI_TMAX_clim_month_.*_mod1_all.tif", full.names=T))
  lf_mosaics_reg <- vector("list",length=length(l_reg_name))
  for(i in 1:length(l_reg_name)){
    lf_mosaics_reg[[i]] <- try(mixedsort(
    list.files(
    path=file.path(out_dir,"mosaics"),
    #pattern="reg6_.*._CAI_TMAX_clim_month_.*._mod1_all_mean.tif",
    pattern=paste(l_reg_name[i],".*._CAI_TMAX_clim_month_.*._mod1_all_mean.tif",sep=""), 
    full.names=T))
    )
  }
  
  #now mask and potl
  lf_mosaics_mask_reg <- vector("list",length=length(l_reg_name))
  for(i in 1:length(l_reg_name)){
    
    #
    lf_m <- lf_mosaics_reg[[i]]
    out_dir_str <- out_dir
    reg_name <- paste(l_reg_name[i],"_",tile_size,sep="") #make this automatic
    #lapply()
    #list_param_plot_daily_mosaics <- list(lf_m=lf_m,reg_name=reg_name,out_dir_str=out_dir_str,out_suffix=out_suffix,l_dates=day_to_mosaic)
    #lf_m_mask_reg4_1500x4500 <- mclapply(1:2,FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics,mc.preschedule=FALSE,mc.cores = 6)
    #lf_mosaics_mask_reg[[i]] <- lapply(1:1,FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics)
    lf_raster_fname <- lf_m
    prefix_str <- "Figure10_clim_reg4_mosaics_day_"
    l_dates <-day_to_mosaic
    screenRast=FALSE
    list_param_plot_screen_raster <- list(lf_raster_fname,screenRast,l_dates,out_dir,prefix_str,out_suffix)
    names(list_param_plot_screen_raster) <- c("lf_raster_fname","screenRast","l_dates","out_dir_str","prefix_str","out_suffix_str")

    #undebug(plot_screen_raster_val)

    #world_m_list1<- plot_screen_raster_val(1,list_param_plot_screen_raster)
    #world_m_list <- mclapply(11:30, list_param=list_param_plot_screen_raster, plot_screen_raster_val,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
    lf_mosaics_mask_reg[[i]] <- mclapply(1:length(l_dates), list_param=list_param_plot_screen_raster, plot_screen_raster_val,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement


    #lf_mosaics_mask_reg[[i]] <- mclapply(1:length(lf_m),FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics,mc.preschedule=FALSE,mc.cores = num_cores)
  }
}

################## WORLD MOSAICS NEEDS MAJOR CLEAN UP OF CODE HERE
##make functions!!
###Combine mosaics with modified code from Alberto

#use list from above!!

# test_list <-list.files(path=file.path(out_dir,"mosaics"),    
#            pattern=paste("^world_mosaics.*.tif$",sep=""), 
# )
# #world_mosaics_mod1_output1500x4500_km_20101105_run10_1500x4500_global_analyses_03112015.tif
# 
# #test_list<-lapply(1:30,FUN=function(i){lapply(1:x[[i]]},x=lf_mosaics_mask_reg)
# #test_list<-unlist(test_list)
# #mosaic_list_mean <- vector("list",length=1)
# mosaic_list_mean <- test_list 
# out_rastnames <- "world_test_mosaic_20100101"
# out_path <- out_dir
# 
# list_param_mosaic<-list(mosaic_list_mean,out_path,out_rastnames,file_format,NA_flag_val,out_suffix)
# names(list_param_mosaic)<-c("mosaic_list","out_path","out_rastnames","file_format","NA_flag_val","out_suffix")
# #mean_m_list <-mclapply(1:12, list_param=list_param_mosaic, mosaic_m_raster_list,mc.preschedule=FALSE,mc.cores = 6) #This is the end bracket from mclapply(...) statement
# 
# lf <- mosaic_m_raster_list(1,list_param_mosaic)
# 
# debug(mosaic_m_raster_list)
# mosaic_list<-list_param$mosaic_list
# out_path<-list_param$out_path
# out_names<-list_param$out_rastnames
# file_format <- list_param$file_format
# NA_flag_val <- list_param$NA_flag_val
# out_suffix <- list_param$out_suffix

##Now mosaic for mean: should reorder files!!
#out_rastnames_mean<-paste("_",lst_pat,"_","mean",out_suffix,sep="")
#list_date_names<-c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
#mosaic_list<-split(tmp_str1,list_date_names)
#new_list<-vector("list",length(mosaic_list))
# for (i in 1:length(list_date_names)){
#    j<-grep(list_date_names[i],mosaic_list,value=FALSE)
#    names(mosaic_list)[j]<-list_date_names[i]
#    new_list[i]<-mosaic_list[j]
#  }
#  mosaic_list_mean <-new_list #list ready for mosaicing
#  out_rastnames_mean<-paste(list_date_names,out_rastnames_mean,sep="")
  
  ### Now mosaic tiles...Note that function could be improved to use less memory


################## PLOTTING WORLD MOSAICS ################

#lf_world_pred <-list.files(path=file.path(out_dir,"mosaics"),    
#           pattern=paste("^world_mosaics.*.tif$",sep=""),full.names=T) 

lf_world_pred <-list.files(path=file.path(out_dir,"mosaics"),    
           pattern=paste("^reg4.*.",".tif$",sep=""),full.names=T) 
l_reg_name <- unique(df_tile_processed$reg)
lf_world_pred <-list.files(path=file.path(out_dir,l_reg_name[[i]]),    
           pattern=paste(".tif$",sep=""),full.names=T) 

#mosaic_list_mean <- test_list 
#out_rastnames <- "world_test_mosaic_20100101"
#out_path <- out_dir

#lf_world_pred <- list.files(pattern="world.*2010090.*.tif$")
#lf_raster_fname <- list.files(pattern="world.*2010*.*02162015.tif$",full.names=T)
lf_raster_fname <- lf_world_pred
prefix_str <- "Figure10_clim_world_mosaics_day_"
l_dates <-day_to_mosaic
screenRast=TRUE
list_param_plot_screen_raster <- list(lf_raster_fname,screenRast,l_dates,out_dir,prefix_str,out_suffix)
names(list_param_plot_screen_raster) <- c("lf_raster_fname","screenRast","l_dates","out_dir_str","prefix_str","out_suffix_str")

#undebug(plot_screen_raster_val)

#world_m_list1<- plot_screen_raster_val(1,list_param_plot_screen_raster)
#world_m_list <- mclapply(11:30, list_param=list_param_plot_screen_raster, plot_screen_raster_val,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
world_m_list <- mclapply(1:length(l_dates), list_param=list_param_plot_screen_raster, plot_screen_raster_val,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement

s_pred <- stack(lf_raster_fname)

res_pix <- 1500
col_mfrow <- 3 
row_mfrow <- 2

png(filename=paste("Figure10_levelplot_combined_",region_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

levelplot(s_pred,layers=1:6,col.regions=rev(terrain.colors(255)),cex=4)

dev.off()

#   blues<- designer.colors(6, c( "blue", "white") )
# reds <- designer.colors(6, c( "white","red")  )
# colorTable<- c( blues[-6], reds)
# breaks with a gap of 10 to 17 assigned the white color
# brks<- c(seq( 1, 10,,6), seq( 17, 25,,6)) 
# image.plot( x,y,z,breaks=brks, col=colorTable)
#

#lf_world_mask_reg <- vector("list",length=length(lf_world_pred))
#for(i in 1:length(lf_world_pred)){
  #
#  lf_m <- lf_mosaics_reg[i]
#  out_dir_str <- out_dir
  #reg_name <- paste(l_reg_name[i],"_",tile_size,sep="") #make this automatic
  #lapply()
#  list_param_plot_daily_mosaics <- list(lf_m=lf_m,reg_name=tile_size,out_dir_str=out_dir_str,out_suffix=out_suffix,l_dates=day_to_mosaic)
  #lf_m_mask_reg4_1500x4500 <- mclapply(1:2,FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics,mc.preschedule=FALSE,mc.cores = 6)

  #lf_world_mask_reg[[i]] <- mclapply(1:length(lf_m),FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics,mc.preschedule=FALSE,mc.cores = 10)
}

############# NEW MASK AND DATA
## Plot areas and day predicted as first check
  
l_reg_name <- unique(df_tile_processed$reg)
#(subset(df_tile_processed$reg == l_reg_name[i],date)

for(i in 1:length(l_reg_name)){
  lf_world_pred<-list.files(path=file.path(out_dir,l_reg_name[[i]]),    
           pattern=paste(".tif$",sep=""),full.names=T) 

  #mosaic_list_mean <- test_list 
  #out_rastnames <- "world_test_mosaic_20100101"
  #out_path <- out_dir

  #lf_world_pred <- list.files(pattern="world.*2010090.*.tif$")
  #lf_raster_fname <- list.files(pattern="world.*2010*.*02162015.tif$",full.names=T)
  lf_raster_fname <- lf_world_pred
  prefix_str <- paste("Figure10_",l_reg_name[i],sep="")

  l_dates <- basename(lf_raster_fname)
  tmp_name <- gsub(".tif","",l_dates)
  tmp_name <- gsub("gam_CAI_dailyTmax_predicted_mod1_0_1_","",tmp_name)
  #l_dates <- tmp_name
  l_dates <- paste(l_reg_name[i],"_",tmp_name,sep="")

  screenRast=TRUE
  list_param_plot_screen_raster <- list(lf_raster_fname,screenRast,l_dates,out_dir,prefix_str,out_suffix)
  names(list_param_plot_screen_raster) <- c("lf_raster_fname","screenRast","l_dates","out_dir_str","prefix_str","out_suffix_str")

  #undebug(plot_screen_raster_val)

  #world_m_list1<- plot_screen_raster_val(1,list_param_plot_screen_raster)
  #world_m_list <- mclapply(11:30, list_param=list_param_plot_screen_raster, plot_screen_raster_val,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
  world_m_list <- mclapply(1:length(l_dates), list_param=list_param_plot_screen_raster, plot_screen_raster_val,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement

  #s_pred <- stack(lf_raster_fname)

  #res_pix <- 1500
  #col_mfrow <- 3 
  #row_mfrow <- 2

  #png(filename=paste("Figure10_levelplot_combined_",region_name,"_",out_suffix,".png",sep=""),
  #  width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  #levelplot(s_pred,layers=1:6,col.regions=rev(terrain.colors(255)),cex=4)

  #dev.off()
}


  
##################### END OF SCRIPT ######################
