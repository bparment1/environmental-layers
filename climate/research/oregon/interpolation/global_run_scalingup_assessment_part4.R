##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part2 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Accuracy methods are added in the the function scripts to evaluate results.
#Analyses, focusing on comparison across tile size and resolution
#AUTHOR: Benoit Parmentier 
#CREATED ON: 12/31/2014  
#MODIFIED ON: 01/14/2015            
#Version: 3
#PROJECT: Environmental Layers project     
#COMMENTS: analyses for run 10 global analyses, Europe 1000x300km
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

#Not working yet...prolem with resampling!!
generate_difference_images_fun(i,list_param){

  #Purpose:
  #This functions produce difference image between raster series.
  #It creates a figure of image differencing.
  #Author: Benoit Parmentier
  #Parameters:
  #lf_1: list of files for the first raster series/stack.
  #reg_name:region name with tile size included
  #To do:
  #Add filenames
  #Add range
  #Add output dir
  #Add dataType_val option
  
  ##### BEGIN ########
  
  #Parse the list of parameters
  lf1 <- list_param$lf1
  lf2 <- list_param$lf2
  lf_out <- list_param$lf_out
  
  reg_name_str <- list_param$reg_name_str
  out_dir_str <- list_param$out_dir_str
  out_suffix <- list_param$out_suffix
  
  #rast_list <- vector("list",length=length(lf_m))
  r1 <- raster(lf1[i])
  r2 <- raster(lf2[i])

  #paste(raster_name[1:7],collapse="_")
  #add filename option later
  #extension_str <- extension(filename(r1))
  #raster_name_tmp <- gsub(extension_str,"",basename(filename(r_test)))
  raster_name <- file.path(out_dir_str,lf_out[i])

  #Find common extent (union)
  
  #if(extent(r1),extent(r2) #if extent is different then...
  r1_out <- merge(r1,extent(r2))
  r1_out <- extend(r1,extent(r2))
  
  r2_out <- extend(r2,extent(r1_out))

  #Extend first image if necssary?
  r2_out <- projectRaster(r2_out,to=r1_out,alignOnly=T)
  projectRaster(r1_out)
  r_test <- resample(r2,r1_out)
  #if res different then
  
  #reproject to match res...use the resolution of the first image...
  
  #
  
  
  r_out <- overlay(r1,r2, fun=function(x,y){x-y},filename=raster_name,
                 dataType="FLT4S",overwrite=T)
  
  res_pix <- 1200
  #res_pix <- 480

  col_mfrow <- 1
  row_mfrow <- 1
  
  
  png(filename=paste("Figure_image_diff","_",date_proc,"_",reg_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_out,main=paste("climatology month ",date_proc , " ", reg_name,sep=""),cex.main=1.5)
  dev.off()
  
  return(raster_name)

}

#Plot image difference for smaller zoom windows
plot_diff_w_fun <- function(i,list_param){
  #plot difference and windows...
  #Assumes a stack of three images:
  #First iamge: difference
  #Second: first ref image e.g. mosaic of predictions at tile size of 1000x3000km
  #Third: second ref image e.g. mosaic of predicions at tile size of 1500x3500km
  
  ############
  
  r_stack <- list_param$r_stack
  l_windows <- list_param$l_windows
  date_proc <- list_param$date_proc
  reg_name <- list_param$reg_name
  out_suffix_str <- list_param$out_suffix
  
  #r_w_d2 <- list_param$r_w_d2
  
  r_w <- raster(l_windows[[i]])
  r_diff <- subset(r_stack,1)
  r_w_diff <- crop(r_diff,r_w)
  r_w_d1 <- crop(subset(r_stack,2),r_w)
  r_w_d2 <- crop(subset(r_stack,3),r_w)
  r_s <- stack(r_w_diff,r_w_d1,r_w_d2)
  names(r_s) <- names(r_stack)
  
  p0 <- levelplot(r_s,col.regions=matlab.like(100),margin=FALSE,layer=1)
  p1 <- levelplot(r_s,col.regions=matlab.like(100),margin=FALSE,layer=2)
  p2 <- levelplot(r_s,col.regions=matlab.like(100),margin=FALSE,layer=3)

  layout_m <- c(1,3) # works if set to this?? ok set the resolution...
  png(paste("Figure_diff_","_comparison_levelplot_","window_",i,"_",date_proc,"_",reg_name,out_suffix_str,".png", sep=""),
      height=480*layout_m[1],width=480*layout_m[2])
  grid.arrange(p0,p1,p2,ncol=3)
  dev.off() 

  layout_m <- c(1,3) # works if set to this?? ok set the resolution...
  png(paste("Figure_diff_","_comparison_plot_","window_",i,"_",date_proc,"_",reg_name,out_suffix_str,".png", sep=""),
      height=480*layout_m[1],width=480*layout_m[2])
  plot(r_s)
  dev.off() 

  return(r_s)
}

############################################
#### Parameters and constants  

out_dir <-"/data/project/layers/commons/NEX_data/output_run10_global_analyses_12232014/mosaics_differences"

y_var_name <- "dailyTmax"
interpolation_method <- c("gam_CAI")
out_prefix<-"run10_global_analyses_12232014"
mosaic_plot <- FALSE

#CRS_locs_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
CRS_WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

proj_str<- CRS_WGS84
file_format <- ".rst"
NA_value <- -9999
NA_flag_val <- NA_value
out_suffix <- "_12232014" 

#out_dir <-paste(out_dir,"_",out_prefix,sep="")
create_out_dir_param <- FALSE

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_prefix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

setwd(out_dir)
                                   
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
CRS_WGS84<-c("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

region_names <- c("reg5") #Africa

l_windows <- c("reg5_1500x4500_window_1.rst",
               "reg5_1500x4500_window_2.rst",
               "reg5_1500x4500_window_3.rst",
               "reg5_1500x4500_window_4.rst")

########################## START SCRIPT ##############################


############### PART I: CREATE MOSAICS MASKED WITH FIGURES ##############

### for reg5_1500x4500: use "mod1 in this case

out_prefix_str <- "reg5_1500x4500"
lf_mosaics_reg5 <- mixedsort(
           list.files(path=
                        "/data/project/layers/commons/NEX_data/output_run10_global_analyses_12152014/mosaics/reg5_1500x4500",
           pattern=".*._mod1_all_mean.tif$",full.names=T)
           )
lf_m <- lf_mosaics_reg5
out_dir_str <- out_dir
reg_name <- "reg5_1500x4500"
#lapply()
list_param_plot_daily_mosaics <- list(lf_m=lf_m,reg_name=reg_name,out_dir_str=out_dir_str,out_suffix=out_prefix)
#lf_m_mask_reg4_1500x4500 <- mclapply(1:2,FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics,mc.preschedule=FALSE,mc.cores = 6)

lf_m_mask_reg5_1500x4500 <- mclapply(1:length(lf_m),FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics,mc.preschedule=FALSE,mc.cores = 6)

### for reg5_1500x4500: use "mod1 in this case, this is Africa

out_prefix_str <- "reg5_1000x3000"
lf_mosaics_reg5 <- mixedsort(
           list.files(path=
                        "/data/project/layers/commons/NEX_data/output_run10_global_analyses_12152014/mosaics/reg5_1000x3000",
           pattern=".*._mod1_all_mean.tif$",full.names=T)
           )
lf_m <- lf_mosaics_reg5
out_dir_str <- out_dir
reg_name <- "reg5_1000x3000"
#lapply()
list_param_plot_daily_mosaics <- list(lf_m=lf_m,reg_name=reg_name,out_dir_str=out_dir_str,out_suffix=out_prefix)
#lf_m_mask_reg4_1500x4500 <- mclapply(1:2,FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics,mc.preschedule=FALSE,mc.cores = 6)

lf_m_mask_reg5_1000x3000 <- mclapply(1:length(lf_m),FUN=plot_daily_mosaics,list_param=list_param_plot_daily_mosaics,mc.preschedule=FALSE,mc.cores = 6)


############### PART II: IMAGE DIFFERENCE FOR MOSAICS MASKED WITH FIGURES ##############

#Uing temporary solution with IDRISI!

rt2 <-raster("reg5_1500x4500_CAI_TMAX_clim_month_20100101_mod1_all_mean_20100101_masked_rec_out.rst")
rt1 <-raster("reg5_1000x300_CAI_TMAX_clim_month_20100101_mod1_all_mean_20100101_masked_rec.rst")
r_test <- rt1 - rt2

#> r_test
#class       : RasterLayer 
#dimensions  : 9439, 9108, 85970412  (nrow, ncol, ncell)
#resolution  : 0.008333433, 0.008334769  (x, y)
#extent      : -21.41903, 54.48188, -38.67189, 40  (xmin, xmax, ymin, ymax)
#coord. ref. : +proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs 
#data source : /tmp/R_raster_tmp/parmentier/raster_tmp_2015-01-05_110752_51464.grd 
#names       : layer 
#values      : -11.3716, 8.497058  (min, max)
plot(r_test)

writeRaster(r_test,file.path(out_dir,"image_diff_reg5_1000x3000_by_1500x4500_20100101.rst"))

res_pix <- 1200
#res_pix <- 480

col_mfrow <- 1
row_mfrow <- 1

reg_name <- "reg5_1000x3000_by_1500x4500"
date_proc <- "20100101"
png(filename=paste("Figure_image_diff","_",date_proc,"_",reg_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_test,main=paste("Image differencing ",date_proc, " ", reg_name,sep=""),cex.main=1.5)
dev.off()

r_stack <- stack(r_test,rt1,rt2)
names(r_stack) <- c("r_diff","reg5_1000x3000","reg5_1500x4500")
p <- levelplot(r_stack,col.regions=matlab.like(200),layers=2:3)

res_pix <- 600

col_mfrow <- 2
row_mfrow <- 1
png(filename=paste("Figure_image_diff_raster_stack","_",date_proc,"_",reg_name,"_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
print(p)
dev.off()

histogram(r_stack)

############### PART III: CREATE IMAGE ZOOM WINDOW AND FIGURES ##############

l_windows <- c("reg5_1500x4500_window_1.rst",
               "reg5_1500x4500_window_2.rst",
               "reg5_1500x4500_window_3.rst",
               "reg5_1500x4500_window_4.rst")

r_w <- lapply(1:length(l_windows),FUN=function(i){raster(file.path(out_dir,l_windows[i]))})
#r_w_diff <- lapply(1:length(l_windows),FUN=function(i){crop(r_diff,r_w[[i]])})
#r_w_1000x3000 <- lapply(1:length(l_windows),FUN=function(i){crop(subset(r_stack,1),r_w[[i]])})
#r_w_1500x4500 <- lapply(1:length(l_windows),FUN=function(i){crop(subset(r_stack,1),r_w[[i]])})

plot(r_w_diff[[1]])

#p1<-levelplot(r_w_diff[[1]],col.regions=matlab.like(100),margin=FALSE)
#print(p1)

r_w_rec <- lapply(1:length(r_w),FUN=function(i){values(x[[i]]) <- rep(1,ncell(x[[i]]))})
#spp_w <- polygonFromExtent(exten(r_w[[1]]), sp=TRUE)
pol <- rasterToPolygons(r_w_rec[[1]])

res_pix <- 1200
#res_pix <- 480

col_mfrow <- 1
row_mfrow <- 1

plot(r_stack,y=1)
plot(pol,add=T)
plot(extent(r_w_diff),add=T)


reg_name <- "reg5_1000x3000_by_1500x4500"
date_proc <- "20100101"

list_param_plot_diff <- list(r_stack,l_windows,date_proc,reg_name,out_suffix)
names(list_param_plot_diff) <- c("r_stack","l_windows","date_proc","reg_name","out_suffix")

#debug(plot_diff_w_fun)

plot_diff_w_fun(1,list_param=list_param_plot_diff)

l_r_s <- lapply(1:length(l_windows),FUN=plot_diff_w_fun,list_param=list_param_plot_diff)


##########################################################
##################### END OF SCRIPT ######################
