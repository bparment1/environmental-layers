####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for assessment of scaling up on NEX ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Accuracy methods are added in the the function scripts to evaluate results.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 03/23/2014  
#MODIFIED ON: 05/01/2014            
#Version: 2
#PROJECT: Environmental Layers project                                     
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

#### FUNCTION USED IN SCRIPT

function_analyses_paper1 <-"contribution_of_covariates_paper_interpolation_functions_10222013.R" #first interp paper
function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_03182014.R"

load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

extract_list_from_list_obj<-function(obj_list,list_name){
  #Create a list of an object from a given list of object using a name prodived as input
  
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp<-obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]]<-tmp
  }
  return(list_tmp) #this is  a data.frame
}

#This extract a data.frame object from raster prediction obj and combine them in one data.frame 
extract_from_list_obj<-function(obj_list,list_name){
  #extract object from list of list. This useful for raster_prediction_obj
  library(plyr)
  
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp<-obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]]<- as.data.frame(tmp) #if spdf
  }
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


##############################
#### Parameters and constants  

#in_dir1 <- "/data/project/layers/commons/NEX_data/test_run1_03232014/output" #On Atlas

in_dir1 <- "/nobackupp4/aguzman4/climateLayers/output" #On NEX
in_dir_list <- list.dirs(path=in_dir1) #get the list of directories with resutls by 10x10 degree tiles
#in_dir_list <- as.list(in_dir_list[-1])
in_dir_list <- in_dir_list[grep("output",basename(in_dir_list),invert=TRUE)] #the first one is the in_dir1
in_dir_shp <- in_dir_list[grep("shapefiles",basename(in_dir_list),invert=FALSE)] #select directory with shapefiles...

in_dir_list <- in_dir_list[grep("shapefiles",basename(in_dir_list),invert=TRUE)] 
#the first one is the in_dir1
# the last directory contains shapefiles 
y_var_name <- "dailyTmax"
interpolation_method <- c("gam_fusion")
out_prefix<-"run1_NA_analyses_03232013"
#out_dir<-"/data/project/layers/commons/NEX_data/"
out_dir <- "/nobackup/bparmen1/"
out_dir <-paste(out_dir,"_",out_prefix,sep="")

#system("ls /nobackup/bparmen1")

if (!file.exists(out_dir)){
  dir.create(out_dir)
  #} else{
  #  out_path <-paste(out_path..)
}

setwd(out_dir)
                                   
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84

##raster_prediction object : contains testing and training stations with RMSE and model object

list_raster_obj_files <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="^raster_prediction_obj.*.RData",full.names=T)})
basename(dirname(list_raster_obj_files[[1]]))
list_names_tile_coord <- lapply(list_raster_obj_files,FUN=function(x){basename(dirname(x))})
list_names_tile_id <- paste("tile",1:length(list_raster_obj_files),sep="_")
#names(list_raster_obj_files)<- 
                                                                                                             
names(list_raster_obj_files)<- list_names_tile_coord


###################### PART I: Generate tables to collect information over all tiles in North America ##########
#Table 1: Average accuracy metrics
#Table 2: daily accuracy metrics for all tiles

##Quick exploration of raster object
robj1 <- load_obj(list_raster_obj_files[[1]])
names(robj1)
names(robj1$method_mod_obj[[1]]) #for January 1, 2010
names(robj1$method_mod_obj[[1]]$dailyTmax) #for January

names(robj1$clim_method_mod_obj[[1]]$data_month) #for January
names(robj1$validation_mod_month_obj[[1]]$data_s) #for January with predictions

data_month_list <- lapply(list_raster_obj_files,FUN=function(x){x<-load_obj(x);x[["summary_metrics_v"]]$avg})                           

robj1$tb_diagnostic_v[1:10,] #first 10 rows of accuarcy metrics per day and model (for specific tile)
robj1$summary_metrics_v #first 10 rows of accuarcy metrics per day and model (for specific tile)

#summary_metrics_v_list <- lapply(list_raster_obj_files,FUN=function(x){x<-load_obj(x);x[["summary_metrics_v"]]$avg$rmse})                           
summary_metrics_v_list <- lapply(list_raster_obj_files,FUN=function(x){x<-load_obj(x);x[["summary_metrics_v"]]$avg})                           
#summary_metrics_v_NA <- do.call(rbind,summary_metrics_v_list) #create a df for NA tiles with all accuracy metrics
names(summary_metrics_v_list) <- list_names_tile_coord
summary_metrics_v_NA <- do.call(rbind.fill,summary_metrics_v_list) #create a df for NA tiles with all accuracy metrics
list_names_tile_id <- paste("tile",1:length(list_raster_obj_files),sep="_")
tile_coord <- lapply(1:length(summary_metrics_v_list),FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=summary_metrics_v_list)
tile_id <- lapply(1:length(summary_metrics_v_list),
                     FUN=function(i,x,y){rep(y[i],nrow(x[[i]]))},x=summary_metrics_v_list,y=list_names_tile_id)

summary_metrics_v_NA$tile_id <-unlist(tile_id)
summary_metrics_v_NA$tile_coord <-unlist(tile_coord)

summary_metrics_v_NA$n <- as.integer(summary_metrics_v_NA$n)
write.table(as.data.frame(summary_metrics_v_NA),
            file=file.path(out_dir,paste("summary_metrics_v2_NA_",out_prefix,".txt",sep="")),sep=",")
#Function to collect all the tables from tiles into a table

tb_diagnostic_v_list <- lapply(list_raster_obj_files,FUN=function(x){x<-load_obj(x);x[["tb_diagnostic_v"]]})                           
names(tb_diagnostic_v_list)

tb_diagnostic_v_NA <- do.call(rbind.fill,tb_diagnostic_v_list) #create a df for NA tiles with all accuracy metrics

tile_coord <- lapply(1:length(tb_diagnostic_v_list),
                     FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=tb_diagnostic_v_list)
tile_id <- lapply(1:length(tb_diagnostic_v_list),
                     FUN=function(i,x,y){rep(y[i],nrow(x[[i]]))},x=tb_diagnostic_v_list,y=list_names_tile_id)

tb_diagnostic_v_NA$tile_id <- unlist(tile_id) #adding identifier for tile
tb_diagnostic_v_NA$tile_coord <- unlist(tile_coord) #adding identifier for tile

write.table((tb_diagnostic_v_NA),
            file=file.path(out_dir,paste("tb_diagnostic_v2_NA","_",out_prefix,".txt",sep="")),sep=",")

#load data_month for specific tiles
data_month <- extract_from_list_obj(robj1$clim_method_mod_obj,"data_month")

names(data_month) #this contains LST means (mm_1, mm_2 etc.) as well as TMax and other info

#problem with tile 12...the raster ojbect has missing sub object
#data_month_list <- lapply(1:length(list_raster_obj_files),x=list_raster_obj_files,
#                          FUN=function(i,x){x<-load_obj(x[[i]]);
#                                            extract_from_list_obj(x$validation_mod_month_obj,"data_s")})                           

data_month_list <- lapply(1:length(list_raster_obj_files),x=list_raster_obj_files,
                          FUN=function(i,x){x<-load_obj(x[[i]]);
                                            extract_from_list_obj(x$clim_method_mod_obj,"data_month")})                           

names(data_month_list) <- paste("tile","_",1:length(data_month_list),sep="")


#names(data_month_list) <- basename(in_dir_list) #use folder id instead

tile_id <- lapply(1:length(data_month_list),
                  FUN=function(i,x){rep(names(x)[i],nrow(x[[i]]))},x=data_month_list)
data_month_NAM <- do.call(rbind.fill,data_month_list) #combined data_month for "NAM" North America
data_month_NAM$tile_id <- unlist(tile_id)

#plot(mm_01 ~ elev_s,data=data_month_NAM) #Jan across all tiles
#plot(mm_06 ~ elev_s,data=data_month_NAM) #June across all tiles
#plot(TMax ~ mm_01,data=data_month_NAM) #monthly tmax averages and LST across all tiles

#copy back to atlas
system("scp -p ./*.txt ./*.tif parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/_run1_NA_analyses_03232013")

####### PART 2 CREATE MOSAIC OF PREDICTIONS PER DAY ###

#list.files(path=in_dir_list[[1]],pattern=".*predicted.*20100101.*.tif")

#some files are
#"/nobackupp4/aguzman4/climateLayers/output/30.0_-115.0/gam_fusion_dailyTmax_predicted_mod_kr_0_1_20100901_30_130.0_-115.0.tif"
#"/nobackupp4/aguzman4/climateLayers/output/30.0_-115.0/gam_fusion_dailyTmax_predicted_mod_kr_20100901_30_130.0_-115.0.tif"
#list_tif_files_dates <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern=".*month.*.tif",full.names=T)})
#list_tif_files_dates <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern=".*predicted.*20100101.*.tif",full.names=T)})

#".*predicted_mod1_0_1.*20100101.*.tif"

tb <- tb_diagnostic_v_NA
#get specific dates from tb
dates_l <- unique(tb$date)

list_tif_fun <- function(i,in_dir_list,pattern_str){
  #list.files(path=x,pattern=".*predicted_mod1_0_1.*20100101.*.tif",full.names=T)})
  pat_str<- pattern_str[i]
  list_tif_files_dates <-lapply(in_dir_list,
         FUN=function(x,pat_str){list.files(path=x,pattern=pat_str,full.names=T)},pat_str=pat_str)
  return(list_tif_files_dates)
} 

#First mosaic mod1

## make this a function? report on number of tiles used for mosaic...

#inputs
l_pattern_models <- lapply(c(".*predicted_mod1_0_1.*",".*predicted_mod2_0_1.*",".*predicted_mod3_0_1.*"),
                           FUN=function(x){paste(x,dates_l,".*.tif",sep="")})
out_prefix_s <- paste(c("gam_fusion_dailyTmax_"),c("predicted_mod1_0_01","predicted_mod2_0_01","predicted_mod3_0_01"),sep="")
dates_l #list of predicted dates
                      
for (i in 1:lenth(l_pattern_models)){
  
  l_pattern_mod <- l_pattern_models[[i]]
  out_prefix_s <-    
    
  l_out_rastnames_var <- paste("gam_fusion_dailyTmax_","predicted_mod1_0_01_",dates_l,sep="")

  #list_tif_files_dates <- list_tif_fun(1,in_dir_list,l_pattern_mod)

  ##List of predicted tif ...
  list_tif_files_dates <-lapply(1:length(l_pattern_mod1),FUN=list_tif_fun, 
                              in_dir_list=in_dir_list,pattern_str=l_pattern_mod)
  #save(list_tif_files_dates, file=paste("list_tif_files_dates","_",out_prefix,".RData",sep=""))


  mosaic_list_var <- list_tif_files_dates
  out_rastnames_var <- l_out_rastnames_var  

  file_format <- ".tif"
  NA_flag_val <- -9999

  j<-1
  list_param_mosaic<-list(j,mosaic_list_var,out_rastnames_var,out_dir,file_format,NA_flag_val)
  names(list_param_mosaic)<-c("j","mosaic_list","out_rastnames","out_path","file_format","NA_flag_val")
  list_var_mosaiced <- mclapply(1:365,FUN=mosaic_m_raster_list,list_param=list_param_mosaic,mc.preschedule=FALSE,mc.cores = 2)

}





### Now find out how many files were predicted

l_pattern_mod1<-paste(".*predicted_mod1_0_1.*",dates_l,".*.tif",sep="")

l_f_t12 <- list.files(path=in_dir_list[12],".*predicted_mod1_0_1.*")


l_f_bytiles<-lapply(in_dir_list,function(x){list.files(path=x,pattern=".*predicted_mod1_0_1.*")})


unlist(lapply(l_f_bytiles,length))

### Create a combined shape file for all region?

#get centroid
#plot the average RMSE at the centroid??
#quick kriging for every tile?
                    
### Create a combined boxplot for every tile (can also do that in pannel)
### Create a quick mosaic (mean method...)
#mean predicitons
#mean of kriging error?
                    
tb <- read.table(file.path(out_dir,"tb_diagnostic_v2_NA_run1_NA_analyses_03232013.txt"),sep=",")
summary_metrics_v <- read.table(file.path(out_dir,"summary_metrics_v2_NA_run1_NA_analyses_03232013.txt"),sep=",")

strsplit(unique(tb$tile_id),"_")
tb$tile_id <- factor(tb$tile_id, levels=unique(tb$tile_id))
#tb$tile_id <- as.character(tb$tile_id)
boxplot(rmse~tile_id,data=subset(tb,tb$pred_mod=="mod1"),
        names=1:24)
title("RMSE per tile")
#bwplot(rmse~as.factor(tile_id), data=subset(tb,tb$pred_mod=="mod1"))

#Boxplot of RMSE by model
boxplot(rmse~pred_mod,data=tb)#,names=tb$pred_mod)
title("RMSE per model over 24 tiles")

#Turn summary table to a point shp

tx<-strsplit(as.character(summary_metrics_v$tile_coord),"_")
lat<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][1]},x=tx))
long<- as.numeric(lapply(1:length(tx),function(i,x){x[[i]][2]},x=tx))
summary_metrics_v$lat <- lat
summary_metrics_v$lon <- long

coordinates(summary_metrics_v) <- cbind(long,lat) 
 
ac_mod1 <- summary_metrics_v[summary_metrics_v$pred_mod=="mod1",]
  
plot(r2)  
#plot(ac_mod1,cex=sqrt(ac_mod1$rmse),pch=1,add=T)
plot(ac_mod1,cex=(ac_mod1$rmse^2)/10,pch=1,add=T)
plot(ac_mod1,cex=(ac_mod1$rmse1)*2,pch=1,add=T)

df <-arrange(as.data.frame(ac_mod1),desc(rmse))[,c("rmse","mae","tile_id")]
#View(df)
#quick kriging...
autokrige(rmse~1,r2,)
### COMBINED SHAPEFILES AND EXAMING CENTROID

##Get the list of shapefiles
in_dir_list_NEX <- read.table("/data/project/layers/commons/NEX_data/test_run1_03232014/output/in_dir_list.txt",sep=" ")

list.files(path=in_dir_shp,pattern=paste("*.",as.character(in_dir_list_NEX[1,1]),".*.shp",sep=""))
list_shp_reg_files <- list.files(path=in_dir_shp,pattern="*.shp")
list_shp_reg_files2 <- list_shp_reg_files[grep("shapefiles",basename(in_dir_list),invert=FALSE)] #select directory with shapefiles...

tx<-strsplit(as.character(in_dir_list_NEX[,1]),"_")
lat_shp<- lapply(1:length(tx),function(i,x){x[[i]][1]},x=tx)
long_shp<- lapply(1:length(tx),function(i,x){x[[i]][2]},x=tx)

summary_metrics
#print.numeric<-function(x, digits = 2) formatC(x, digits = digits, format = "f")
#print.numeric(long_shp[[1]])

pattern_str<-lapply(1:length(long_shp),function(i,y,x){paste(y[[i]],"0_",x[[i]],"0",sep="")},y=lat_shp,x=long_shp)
#Select the 24 tiles that were used in the predictions based on lat, long
list_shp_reg_files <- lapply(pattern_str,FUN=function(x){list.files(path=in_dir_shp,
                                                                        pattern=paste("*.",x,".*.shp$",sep=""),full.names=T)})
###
#OK now get the shapefiles merged
#ghcn_dat <- readOGR(dsn=dirname(met_stations_obj$monthly_covar_ghcn_data),
#        sub(".shp","",basename(met_stations_obj$monthly_covar_ghcn_data)))

#plot(r44)

#usa_map <- getData('GADM', country='USA', level=1) #Get US map
usa_map <- getData('GADM', country='USA',level=1) #Get US map, this is not working right now
usa_map_2 <- usa_map[usa_map$NAME_1!="Alaska",] #remove Alaska
usa_map_2 <- usa_map_2[usa_map_2$NAME_1!="Hawaii",] #remove Hawai 

centroids_pts <- vector("list",length(list_shp_reg_files))
shps_tiles <- vector("list",length(list_shp_reg_files))
#collect info
for(i in 1:length(list_shp_reg_files)){
  shp1<-readOGR(dirname(list_shp_reg_files[[i]]),sub(".shp","",basename(list_shp_reg_files[[i]])))
  pt <- gCentroid(shp1)
  centroids_pts[[i]] <-pt
  shps_tiles[[i]] <- shp1
}
#plot info: with labels
plot(usa_map_2)
for(i in 1:length(list_shp_reg_files)){
  shp1 <- shps_tiles[[i]]
  pt <- centroids_pts[[i]]
  plot(shp1,add=T,border="blue")
  #plot(pt,add=T,cex=2,pch=5)
  text(coordinates(pt)[1],coordinates(pt)[2],labels=i,cex=1,col=c("red"))
}
title("Tiles location 10x10 degrees for NAM")

## Now plot RMSE: with labels
plot(usa_map_2)
for(i in 1:length(list_shp_reg_files)){
  shp1 <- shps_tiles[[i]]
  pt <- centroids_pts[[i]]
  plot(shp1,add=T,border="blue")
  #plot(pt,add=T,cex=2,pch=5)
  text(coordinates(pt)[1],coordinates(pt)[2],labels=i,cex=1,col=c("red"))
}
title("Tiles location 10x10 degrees for NAM")


#unique(summaty_metrics$tile_id)
#text(lat-shp,)
#union(list_shp_reg_files[[1]],list_shp_reg_files[[2]])

l_m_tif <- list.files(pattern="*mosaiced.*.tif")
r1 <- raster("mosaiced_gam_fusion_dailyTmax_predicted_mod1_0_01_20100101.tif")
r2 <- raster("mosaiced_gam_fusion_dailyTmax_predicted_mod1_0_01_20100901.tif")

pred_s <- stack(l_m_tif[2:5])
levelplot(pred_s)
