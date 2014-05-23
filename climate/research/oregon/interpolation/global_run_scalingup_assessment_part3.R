##############################  INTERPOLATION OF TEMPERATURES  #######################################
#######################  Script for assessment of scaling up on NEX : part3 ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#In part 3, models and results are assessed on a tile basis using modifications of method assessment 
#Analyses, figures, tables and data are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 05/21/2014  
#MODIFIED ON: 05/21/2014            
#Version: 1
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

function_analyses_paper1 <- "contribution_of_covariates_paper_interpolation_functions_05212014.R" #first interp paper
function_analyses_paper2 <- "multi_timescales_paper_interpolation_functions_05052014.R"
function_assessment_by_tile <- "results_interpolation_date_output_analyses_05212014.R"
#source(file.path(script_path,"results_interpolation_date_output_analyses_08052013.R"))

script_path<-"/home/parmentier/Data/IPLANT_project/env_layers_scripts/" #path to script
source(file.path(script_path,function_analyses_paper1)) #source all functions used in this script 1.
source(file.path(script_path,function_analyses_paper2)) #source all functions used in this script 2.
source(file.path(script_path,function_assessment_by_tile)) #source all functions used in this script 2.

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

Atlas_server <- TRUE #data on NCEAS Atlas
NEX_sever <- TRUE #data on NEX NASA

#on ATLAS
#if(Atlas_server==TRUE){
#
#}

#in_dir1 <- "/data/project/layers/commons/NEX_data/test_run1_03232014/output" #On Atlas
#parent output dir : contains subset of the data produced on NEX
in_dir1 <- "/data/project/layers/commons/NEX_data/output_run2_05122014/output/"
# parent output dir for the curent script analyes
out_dir <-"/data/project/layers/commons/NEX_data/" #On NCEAS Atlas
# input dir containing shapefiles defining tiles
in_dir_shp <- "/data/project/layers/commons/NEX_data/output_run2_05122014/output/subset/shapefiles"

#On NEX
#contains all data from the run by Alberto
#in_dir1 <- "/nobackupp4/aguzman4/climateLayers/output4" #On NEX
#parent output dir for the current script analyes
#out_dir <- "/nobackup/bparmen1/" #on NEX
#in_dir_shp <- "/nobackupp4/aguzman4/climateLayers/output4/subset/shapefiles/"

y_var_name <- "dailyTmax"
interpolation_method <- c("gam_CAI")
out_prefix<-"run2_global_analyses_05122014"

#out_dir <-paste(out_dir,"_",out_prefix,sep="")
create_out_dir_param <- TRUE

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_prefix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}
setwd(out_dir)

df_tile_processed <- read.table(file=file.path(out_dir,paste("df_tile_processed_",out_prefix,".txt",sep="")),sep=",")
#in_dir_list <- file.path(in_dir1,read.table(file.path(in_dir1,"processed.txt"))$V1)
                              
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
region_name <- "USA"

###Table 1: Average accuracy metrics
###Table 2: daily accuracy metrics for all tiles

summary_metrics_v <- read.table(file=file.path(out_dir,paste("summary_metrics_v2_NA_",out_prefix,".txt",sep="")),sep=",")
tb <- read.table(file=file.path(out_dir,paste("tb_diagnostic_v_NA","_",out_prefix,".txt",sep="")),sep=",")
#df_tile_processed <- read.table(file=file.path(out_dir,paste("df_tile_processed_",out_prefix,".txt",sep="")),sep=",")

########################## START SCRIPT ##############################

#Now add things here...
#
selected_tiles <- c("45.0_-120.0","35.0_-115.0")
##raster_prediction object : contains testing and training stations with RMSE and model object
in_dir_list <- list.files(path=in_dir1,full.names=T)
in_dir_list <- in_dir_list[grep("subset",basename(basename(in_dir_list)),invert=TRUE)] #the first one is the in_dir1
in_dir_list <- in_dir_list[match(selected_tiles,basename(basename(in_dir_list)))] #the first one is the in_dir1

list_raster_obj_files <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="^raster_prediction_obj.*.RData",full.names=T)})
#list_names_tile_coord <- lapply(list_raster_obj_files,FUN=function(x){basename(dirname(x))})
#list_names_tile_id <- paste("tile",1:length(list_raster_obj_files),sep="_")
#names(list_raster_obj_files)<- list_names_tile_id

lf_covar_obj <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar_obj.*.RData",full.names=T)})
lf_covar_tif <- lapply(in_dir_list,FUN=function(x){list.files(path=x,pattern="covar.*.tif",full.names=T)})

list_shp_reg_files <- file.path(in_dir_shp,df_tile_processed$shp_files)
list_shp_reg_files <- file.path(in_dir_shp,df_tile_processed$shp_files)

###############
### Figure 1: plot location of the study area with tiles processed

### Figures diagnostic tile:
#Use stage 5 modified/updated code

##Quick exploration of raster object

date_selected_results <- c("20100101") 
raster_prediction_obj <- list_raster_obj_files[[1]]
in_path_tile <- in_dir_list[[1]] #Oregon tile
#in_path_tile <- NULL # set to NULL if the script is run on the NEX node as part of job
covar_obj <- lf_covar_obj[[1]] 

var <- "TMAX"
list_param_results_analyses<-list(out_dir,in_path,script_path,raster_prediction_obj,interpolation_method,
                                  covar_obj,date_selected_results,var,out_prefix)
names(list_param_results_analyses)<-c("out_path","in_path_tile","script_path","raster_prediction_obj","interpolation_method",
                     "covar_obj","date_selected_results","var","out_prefix")
list_param <- list_param_results_analyses

#Run modified code from stage 5...
#plots_assessment_by_date<-function(j,list_param){
#Use lapply or mclapply
#debug(plots_assessment_by_date)
summary_v_day <- plots_assessment_by_date(1,list_param_results_analyses)
#Call as function...

#Boxplots...etc...

#This is a repeat from earlier code.

##Create data.frame with validation and fit metrics for a full year/full numbe of runs

#Call functions to create plots of metrics for validation dataset
tile_selected <- 6
tb_diagnostic_v <- subset(tb,tile_id==6)
metric_names<-c("rmse","mae","me","r","m50")

summary_metrics_v<- boxplot_from_tb(tb_diagnostic_v,metric_names,out_prefix,out_path) #if adding for fit need to change outprefix

names(summary_metrics_v)<-c("avg","median")

summary_month_metrics_v<- boxplot_month_from_tb(tb_diagnostic_v,metric_names,out_prefix,out_path)

#Call functions to create plots of metrics for validation dataset

metric_names<-c("rmse","mae","me","r","m50")

summary_metrics_v<- boxplot_from_tb(tb_diagnostic_v,metric_names,out_prefix,out_path) #if adding for fit need to change outprefix

names(summary_metrics_v)<-c("avg","median")

summary_month_metrics_v<- boxplot_month_from_tb(tb_diagnostic_v,metric_names,out_prefix,out_path)



#### Function to plot boxplot from data.frame table of accuracy metrics


### need to improve these
boxplot_from_tb <-function(tb_diagnostic,metric_names,out_prefix,out_path){
  #now boxplots and mean per models
  library(gdata) #Nesssary to use cbindX
  
  ### Start script
  y_var_name<-unique(tb_diagnostic$var_interp) #extract the name of interpolated variable: dailyTmax, dailyTmin
  
  mod_names<-sort(unique(tb_diagnostic$pred_mod)) #models that have accuracy metrics
  t<-melt(tb_diagnostic,
          #measure=mod_var, 
          id=c("date","pred_mod","prop"),
          na.rm=F)
  t$value<-as.numeric(t$value) #problem with char!!!
  avg_tb<-cast(t,pred_mod~variable,mean)
  avg_tb$var_interp<-rep(y_var_name,times=nrow(avg_tb))
  median_tb<-cast(t,pred_mod~variable,median)
  
  #avg_tb<-cast(t,pred_mod~variable,mean)
  tb<-tb_diagnostic
 
  #mod_names<-sort(unique(tb$pred_mod)) #kept for clarity
  tb_mod_list<-lapply(mod_names, function(k) subset(tb, pred_mod==k)) #this creates a list of 5 based on models names
  names(tb_mod_list)<-mod_names
  #mod_metrics<-do.call(cbind,tb_mod_list)
  #debug here
  if(length(tb_mod_list)>1){
    mod_metrics<-do.call(cbindX,tb_mod_list) #column bind the list??
  }else{
    mod_metrics<-tb_mod_list[[1]]
  }
  
  test_names<-lapply(1:length(mod_names),function(k) paste(names(tb_mod_list[[1]]),mod_names[k],sep="_"))
  #test names are used when plotting the boxplot for the different models
  names(mod_metrics)<-unlist(test_names)
  rows_total<-lapply(tb_mod_list,nrow)
  for (j in 1:length(metric_names)){
    metric_ac<-metric_names[j]
    mod_pat<-glob2rx(paste(metric_ac,"_*",sep=""))   
    mod_var<-grep(mod_pat,names(mod_metrics),value=TRUE) # using grep with "value" extracts the matching names     
    #browser()
    test<-mod_metrics[mod_var]
    png(file.path(out_path,paste("boxplot_metric_",metric_ac, out_prefix,".png", sep="")))
    #boxplot(test,outline=FALSE,horizontal=FALSE,cex=0.5,
    #        ylab=paste(metric_ac,"in degree C",sep=" "))
    
    boxplot(test,outline=FALSE,horizontal=FALSE,cex=0.5,
              ylab=paste(metric_ac,"in degree C",sep=" "),axisnames=FALSE,axes=FALSE)
    axis(1, labels = FALSE)
    ## Create some text labels
    labels <- labels<- names(test)
    ## Plot x axis labels at default tick marks
    text(1:ncol(test), par("usr")[3] - 0.25, srt = 45, adj = 1,
         labels = labels, xpd = TRUE)
    axis(2)
    box()
    #legend("bottomleft",legend=paste(names(rows_total),":",rows_total,sep=""),cex=0.7,bty="n")
    #title(as.character(t(paste(t(names(rows_total)),":",rows_total,sep=""))),cex=0.8)
    title(paste(metric_ac,"for",y_var_name,sep=" "),cex=0.8)
    dev.off()
  }
  
  avg_tb$n<-rows_total #total number of predictions on which the mean is based
  median_tb$n<-rows_total
  summary_obj<-list(avg_tb,median_tb)
  names(summary_obj)<-c("avg","median")
  return(summary_obj)  
}
#boxplot_month_from_tb(tb_diagnostic,metric_names,out_prefix,out_path)
## Function to display metrics by months/seasons
boxplot_month_from_tb <-function(tb_diagnostic,metric_names,out_prefix,out_path){
  
  #Generate boxplot per month for models and accuracy metrics
  #Input parameters:
  #1) df: data frame containing accurayc metrics (RMSE etc.) per day)
  #2) metric_names: metrics used for validation
  #3) out_prefix
  #
  
  #################
  ## BEGIN
  y_var_name<-unique(tb_diagnostic$var_interp) #extract the name of interpolated variable: dailyTmax, dailyTmin  
  date_f<-strptime(tb_diagnostic$date, "%Y%m%d")   # interpolation date being processed
  tb_diagnostic$month<-strftime(date_f, "%m")          # current month of the date being processed
  mod_names<-sort(unique(tb_diagnostic$pred_mod)) #models that have accuracy metrics
  tb_mod_list<-lapply(mod_names, function(k) subset(tb_diagnostic, pred_mod==k)) #this creates a list of 5 based on models names
  names(tb_mod_list)<-mod_names
  t<-melt(tb_diagnostic,
          #measure=mod_var, 
          id=c("date","pred_mod","prop","month"),
          na.rm=F)
  t$value<-as.numeric(t$value) #problem with char!!!
  tb_mod_m_avg <-cast(t,pred_mod+month~variable,mean) #monthly mean for every model
  tb_mod_m_avg$var_interp<-rep(y_var_name,times=nrow(tb_mod_m_avg))
  
  tb_mod_m_sd <-cast(t,pred_mod+month~variable,sd)   #monthly sd for every model
  
  tb_mod_m_list <-lapply(mod_names, function(k) subset(tb_mod_m_avg, pred_mod==k)) #this creates a list of 5 based on models names
  
  for (k in 1:length(mod_names)){
    mod_metrics <-tb_mod_list[[k]]
    current_mod_name<- mod_names[k]
    for (j in 1:length(metric_names)){    
      metric_ac<-metric_names[j]
      col_selected<-c(metric_ac,"month")
      test<-mod_metrics[col_selected]
      png(file.path(out_path,paste("boxplot_metric_",metric_ac,"_",current_mod_name,"_by_month_",out_prefix,".png", sep="")))
      boxplot(test[[metric_ac]]~test[[c("month")]],outline=FALSE,horizontal=FALSE,cex=0.5,
              ylab=paste(metric_ac,"in degree C",sep=" "),,axisnames=FALSE,axes=FALSE)
      #boxplot(test[[metric_ac]]~test[[c("month")]],outline=FALSE,horizontal=FALSE,cex=0.5,
      #        ylab=paste(metric_ac,"in degree C",sep=" "))
      axis(1, labels = FALSE)
      ## Create some text labels
      labels <- month.abb # abbreviated names for each month
      ## Plot x axis labels at default tick marks
      text(1:length(labels), par("usr")[3] - 0.25, srt = 45, adj = 1,
           labels = labels, xpd = TRUE)
      axis(2)
      box()
      #legend("bottomleft",legend=paste(names(rows_total),":",rows_total,sep=""),cex=0.7,bty="n")
      title(paste(metric_ac,"for",current_mod_name,"by month",sep=" "))
      dev.off()
    }  
    
  }
  summary_month_obj <-c(tb_mod_m_list,tb_mod_m_avg,tb_mod_m_sd)
  names(summary_month_obj)<-c("tb_list","metric_month_avg","metric_month_sd")
  return(summary_month_obj)  
}


##################### END OF SCRIPT ######################
