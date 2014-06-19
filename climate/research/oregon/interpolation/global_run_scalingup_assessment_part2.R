##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part2 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#Accuracy methods are added in the the function scripts to evaluate results.
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 03/23/2014  
#MODIFIED ON: 06/19/2014            
#Version: 3
#PROJECT: Environmental Layers project     
#COMMENTS: analyses for run 3 global using 2 specific tiles
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
function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_05052014.R"

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


##############################
#### Parameters and constants  

#on ATLAS
#in_dir1 <- "/data/project/layers/commons/NEX_data/test_run1_03232014/output" #On Atlas
#parent output dir : contains subset of the data produced on NEX
in_dir1 <- "/data/project/layers/commons/NEX_data/output_run3_global_analyses_06192014/output"
# parent output dir for the curent script analyes
#out_dir <-"/data/project/layers/commons/NEX_data/output_run3_global_analyses_06192014/" #On NCEAS Atlas
out_dir <-"/data/project/layers/commons/NEX_data/output_run3_global_analyses_06192014/"
# input dir containing shapefiles defining tiles
in_dir_shp <- "/data/project/layers/commons/NEX_data/output_run3_global_analyses_06192014/output/subset/shapefiles"

#On NEX
#contains all data from the run by Alberto
#in_dir1 <- "/nobackupp4/aguzman4/climateLayers/output4" #On NEX
#parent output dir for the current script analyes
#out_dir <- "/nobackup/bparmen1/" #on NEX
#in_dir_shp <- "/nobackupp4/aguzman4/climateLayers/output4/subset/shapefiles/"

df_tile_processed <- read.table(file=file.path(out_dir,paste("df_tile_processed_",out_prefix,".txt",sep="")),sep=",")
in_dir_list <- file.path(in_dir1,read.table(file.path(in_dir1,"processed.txt"))$V1)
y_var_name <- "dailyTmax"
interpolation_method <- c("gam_CAI")
out_prefix<-"run3_global_analyses_06192014"

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

region_name <- "world"

###Table 1: Average accuracy metrics
###Table 2: daily accuracy metrics for all tiles
#lf_tables <- list.files(out_dir,)
summary_metrics_v <- read.table(file=file.path(out_dir,paste("summary_metrics_v2_NA_",out_prefix,".txt",sep="")),sep=",")
tb <- read.table(file=file.path(out_dir,paste("tb_diagnostic_v_NA","_",out_prefix,".txt",sep="")),sep=",")
df_tile_processed <- read.table(file=file.path(out_dir,paste("df_tile_processed_",out_prefix,".txt",sep="")),sep=",")

########################## START SCRIPT ##############################

###############
### Figure 1: plot location of the study area with tiles processed
list_shp_world <- list.files(in_dir_shp,".shp")
l_shp <- unlist(lapply(1:length(list_shp_world),FUN=function(i){paste(strsplit(list_shp_world[i],"_")[[1]][2:3],collapse="_")}))
list_shp_reg_files <- as.character(df_tile_processed$tile_coord)
#matching_index <- match(l_shp,list_shp_reg_files)
matching_index <- match(list_shp_reg_files,l_shp)
df_tile_processed$shp_files <-list_shp_world[matching_index]

df_tiled_processed <- na.omit(df_tile_processed) #remove other list of folders irrelevant
list_shp_reg_files <- df_tiled_processed$shp_files
list_shp_reg_files<- file.path(in_dir_shp,list_shp_reg_files)

### First get background map to display where study area is located
#can make this more general later on..       
if region_name=="USA"{
  usa_map <- getData('GADM', country='USA', level=1) #Get US map
  #usa_map <- getData('GADM', country=region_name,level=1) #Get US map, this is not working right now
  #usa_map_2 <- usa_map[usa_map$NAME_1!="Alaska",] #remove Alaska
  reg_layer <- usa_map[usa_map$NAME_1!="Hawaii",] #remove Hawai 
}

if region_name=="world"{
  #http://www.diva-gis.org/Data
  countries_shp <-"/data/project/layers/commons/NEX_data/countries.shp"
  path_to_shp <- dirname(countries_shp)
  layer_name <- sub(".shp","",basename(countries_shp))
  reg_layer <- readOGR(path_to_shp, layer_name)
  proj4string(reg_layer) <- CRS_locs_WGS84
  #reg_shp<-readOGR(dirname(list_shp_reg_files[[i]]),sub(".shp","",basename(list_shp_reg_files[[i]])))

}

centroids_pts <- vector("list",length(list_shp_reg_files))
shps_tiles <- vector("list",length(list_shp_reg_files))
#collect info: read in all shapfiles
#This is slow...make a function and use mclapply??
for(i in 1:length(list_shp_reg_files)){
  path_to_shp <- dirname(list_shp_reg_files[[i]])
  layer_name <- sub(".shp","",basename(list_shp_reg_files[[i]]))
  shp1 <- readOGR(path_to_shp, layer_name)
  #shp1<-readOGR(dirname(list_shp_reg_files[[i]]),sub(".shp","",basename(list_shp_reg_files[[i]])))
  pt <- gCentroid(shp1)
  centroids_pts[[i]] <-pt
  shps_tiles[[i]] <- shp1
}

#plot info: with labels
res_pix <- 480
col_mfrow <- 1
row_mfrow <- 1

png(filename=paste("Figure1_tile_processed_region_",region_name,"_",out_prefix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(reg_layer)
#Add polygon tiles...
for(i in 1:length(list_shp_reg_files)){
  shp1 <- shps_tiles[[i]]
  pt <- centroids_pts[[i]]
  plot(shp1,add=T,border="blue")
  #plot(pt,add=T,cex=2,pch=5)
  label_id <- df_tile_processed$tile_id[i]
  text(coordinates(pt)[1],coordinates(pt)[2],labels=i,cex=1,col=c("red"))
}
title(paste("Tiles location 10x10 degrees for ", region_name,sep=""))

dev.off()
      
#unique(summaty_metrics$tile_id)
#text(lat-shp,)
#union(list_shp_reg_files[[1]],list_shp_reg_files[[2]])

###############
### Figure 2: boxplot of average accuracy by model and by tiles

tb$tile_id <- factor(tb$tile_id, levels=unique(tb$tile_id))
model_name <- unique(tb$pred_mod)

## Figure 2a

for(i in  1:length(model_name)){
  
  res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1

  png(filename=paste("Figure2a_boxplot_with_oultiers_by_tiles_",model_name[i],"_",out_prefix,".png",sep=""),
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
  png(filename=paste("Figure2b_boxplot_scaling_by_tiles","_",model_name[i],"_",out_prefix,".png",sep=""),
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

png(filename=paste("Figure3a_boxplot_overall_region_with_oultiers_",model_name[i],"_",out_prefix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

boxplot(rmse~pred_mod,data=tb)#,names=tb$pred_mod)
title("RMSE per model over all tiles")

dev.off()

## Figure 3b
png(filename=paste("Figure3b_boxplot_overall_region_scaling_",model_name[i],"_",out_prefix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

boxplot(rmse~pred_mod,data=tb,ylim=c(0,5),outline=FALSE)#,names=tb$pred_mod)
title("RMSE per model over all tiles")

dev.off()

################ 
### Figure 4: plot predicted tiff for specific date per model

#y_var_name <-"dailyTmax"
#index <-244 #index corresponding to Sept 1
index  <- 1 #index corresponding to Jan 1
date_selected <- "20100901"
name_method_var <- paste(interpolation_method,"_",y_var_name,"_",sep="")

pattern_str <- paste("mosaiced","_",name_method_var,"predicted",".*.",date_selected,".*.tif",sep="")
lf_pred_list <- list.files(pattern=pattern_str)

for(i in 1:length(lf_pred_list)){
  
  r_pred <- raster(lf_pred_list[i])
  
  res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1
  
  png(filename=paste("Figure4_models_predicted_surfaces_",model_name[i],"_",name_method_var,"_",data_selected,"_",out_prefix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
  plot(r_pred)
  title(paste("Mosaiced",model_name[i],name_method_var,date_selected,sep=" "))
  dev.off()
  
}

####### Figure 5...
### Adding tiles do a plot of mod1 with tiles

r_pred <- raster(lf_pred_list[i])
  
res_pix <- 480
col_mfrow <- 1
row_mfrow <- 1
  
png(filename=paste("Figure5_tiles_with_models_predicted_surfaces_",model_name[i],"_",name_method_var,out_prefix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  
plot(r_pred)
title(paste("Mosaiced",model_name[i],name_method_var,date_selected,"with tiles",sep=" "))

#Add polygon tiles...
for(i in 1:length(list_shp_reg_files)){
  shp1 <- shps_tiles[[i]]
  pt <- centroids_pts[[i]]
  plot(shp1,add=T,border="blue")
  #plot(pt,add=T,cex=2,pch=5)
  label_id <- df_tile_processed$tile_id[i]
  text(coordinates(pt)[1],coordinates(pt)[2],labels=i,cex=1,col=c("red"))
}
#title(paste("Prediction with tiles location 10x10 degrees for ", region_name,sep=""))
dev.off()

### 

#### Now combined plot...
#Use the function to match extent...

#pred_s <- stack(lf_list) #problem different extent!!
#methods_names <-c("gam","kriging","gwr")
#methods_names <- interpolation_method

#names_layers<-methods_names[1]
#names_layers <-c("mod1 = lat*long + elev","mod2 = lat*long + elev + LST",
#                 "mod3 = lat*long + elev + LST*FOREST")#, "mod_kr")
#nb_fig<- c("4a","4b")
#list_plots_spt <- vector("list",length=length(lf))
#png(filename=paste("Figure4_models_predicted_surfaces_",date_selected,"_",out_prefix,".png",sep=""),
#    height=480*layout_m[1],width=480*layout_m[2])
#  max_val <- 40
#  min_val <- -40
#  layout_m <- c(1,3) #one row two columns
#  no_brks <- length(seq(min_val,max_val,by=0.1))-1
#  #temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
#  #temp.colors <- colorRampPalette(c('blue', 'lightgoldenrodyellow', 'red'))
#  #temp.colors <- matlab.like(no_brks)
#  temp.colors <- colorRampPalette(c('blue', 'khaki', 'red'))
#  
#  p <- levelplot(pred_s,main=methods_names[i], 
#                 ylab=NULL,xlab=NULL,
#          par.settings = list(axis.text = list(font = 2, cex = 2),layout=layout_m,
#                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),par.strip.text=list(font=2,cex=2),
#          names.attr=names_layers,
#                 col.regions=temp.colors(no_brks),at=seq(min_val,max_val,by=0.1))
##col.regions=temp.colors(25))
#print(p)
#dev.off()

######################
### Figure 6: plot map of average RMSE per tile at centroids

#Turn summary table to a point shp

coordinates(summary_metrics_v) <- cbind(summary_metrics_v$lon,summary_metrics_v$lat)
proj4string(summary_metrics_v) <- CRS_WGS84
lf_list <- lf_pred_list
list_df_ac_mod <- vector("list",length=length(lf_pred_list))
for (i in 1:length(lf_list)){
  
  ac_mod <- summary_metrics_v[summary_metrics_v$pred_mod==model_name[i],]
  r_pred <- raster(lf_list[i])
  
  res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1

  png(filename=paste("Figure6_ac_metrics_map_centroids_tile_",model_name[i],"_",out_prefix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

  plot(r_pred)  

  #plot(ac_mod1,cex=sqrt(ac_mod1$rmse),pch=1,add=T)
  plot(ac_mod,cex=(ac_mod$rmse^2)/10,pch=1,add=T)
  #plot(ac_mod1,cex=(ac_mod1$rmse1)*2,pch=1,add=T)
  title(paste("Averrage RMSE per tile and by ",model_name[i]))

  dev.off()
  
  ### Ranking by tile...
  #df_ac_mod <- 
  list_df_ac_mod[[i]] <- arrange(as.data.frame(ac_mod),desc(rmse))[,c("rmse","mae","tile_id")]
}
  
#quick kriging...
#autokrige(rmse~1,r2,)

######################
### Figure 7: Delta and clim...

## plotting of delta and clim for later scripts...

######################
### Figure 8: Histograms...


##################### END OF SCRIPT ######################
