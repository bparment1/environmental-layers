####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for experimentation of GAM model fitting ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#The purpose is to test different parameters for the GAM function to allow . 
#AUTHOR: Benoit Parmentier 
#CREATED ON: 07/30/2014  
#MODIFIED ON: 07/31/2014            
#Version: 1
#PROJECT: Environmental Layers project  
#TO DO:
# - experiment with gamma
# - experiment with k
#
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

function_analyses_paper1 <-"contribution_of_covariates_paper_interpolation_functions_05212014.R" #first interp paper
function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_05052014.R"

load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

predict_raster_model<-function(in_models,r_stack,out_filename){

  #This functions performs predictions on a raster grid given input models.
  #Arguments: list of fitted models, raster stack of covariates
  #Output: spatial grid data frame of the subset of tiles

  list_rast_pred<-vector("list",length(in_models))

  for (i in 1:length(in_models)){

    mod <-in_models[[i]] #accessing GAM model ojbect "j"
    raster_name<-out_filename[[i]]
    if (inherits(mod,"gam")) {           #change to c("gam","autoKrige")
      raster_pred<- predict(object=r_stack,model=mod,na.rm=FALSE) #Using the coeff to predict new values.
      names(raster_pred)<-"y_pred"  
      writeRaster(raster_pred, filename=raster_name,overwrite=TRUE)  #Writing the data in a raster file format...(IDRISI)
      #print(paste("Interpolation:","mod", j ,sep=" "))
      list_rast_pred[[i]]<-raster_name
    }
  }
  if (inherits(mod,"try-error")) {
    print(paste("no gam model fitted:",mod[1],sep=" ")) #change message for any model type...
  }
  return(list_rast_pred)
}

fit_models<-function(list_formulas,data_training){

  #This functions several models and returns model objects.
  #Arguments: - list of formulas for GAM models
  #           - fitting data in a data.frame or SpatialPointDataFrame

  #Output: list of model objects 

  list_fitted_models<-vector("list",length(list_formulas))
  for (k in 1:length(list_formulas)){
    formula<-list_formulas[[k]]
    mod<- try(gam(formula, data=data_training)) #change to any model!!
    #mod<- try(autoKrige(formula, input_data=data_s,new_data=s_sgdf,data_variogram=data_s))

    model_name<-paste("mod",k,sep="")
    assign(model_name,mod) 
    list_fitted_models[[k]]<-mod
  }
  return(list_fitted_models) 
}


create_dir_fun <- function(out_dir,out_suffix){
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    out_dir <- file.path(out_dir,out_name)
  }
  #create if does not exists: create the output dir as defined 
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  return(out_dir)
}

extract_list_from_list_obj<-function(obj_list,list_name){
  library(plyr)
  
  #Create a list of an object from a given list of object using a name prodived as input
  
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp <- obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]] <- as.data.frame(tmp) #deal with spdf...
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

#scp -rp raster_prediction_obj_gam_CAI_dailyTmax15.0_20.0.RData parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/output_regions/15.0_20.0"
#scp -rp reg*.tif parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/output_regions/15.0_20.0"

in_dir <- "/data/project/layers/commons/NEX_data/output_regions/15.0_20.0"
raster_obj_infile <- "raster_prediction_obj_gam_CAI_dailyTmax15.0_20.0.RData"

setwd(in_dir)

raster_obj<- load_obj(raster_obj_infile)

#raster_obj <- load_obj(unlist(raster_obj_file)) #may not need unlist
nb_models <- length((raster_obj$clim_method_mod_obj[[1]]$formulas))
l_formulas <- (raster_obj$clim_method_mod_obj[[1]]$formulas)

#y_var ~ s(lat, lon) + s(elev_s)
#y_var ~ s(lat, lon) + s(elev_s) + s(LST)
#y_var ~ s(lat, lon) + s(elev_s) + s(N_w, E_w) + s(LST) + ti(LST,LC1) + s(LC1)

pred_mod <- paste("mod",c(1:nb_models,"_kr"),sep="")
#we are assuming no monthly hold out...
#we are assuming only one specific daily prop?
nb_models <- length(pred_mod)

j <- 7
clim_method_mod_obj <- raster_obj$clim_method_mod_obj[[j]]
#this is made of "clim",data_month, data_month_v , sampling_month_dat, mod and formulas
clim_method_mod_obj$clim #file predicted

l_mod <- clim_method_mod_obj$mod #file predicted

reg_rast <- stack(list.files(pattern="*.tif"))
plot(reg_rast,y=15)

names(clim_method_mod_obj)
clim_method_mod_obj$data_month
clim_method_mod_obj$data_month_v

## check for residual pattern, removeable by increasing `k'
## typically `k', below, chould be substantially larger than 
## the original, `k' but certainly less than n/2.
vis.gam(mod1)
vis.gam(mod1,view=c("lat","lon"),theta= 35) # plot against by variable
#http://stats.stackexchange.com/questions/12223/how-to-tune-smoothing-in-mgcv-gam-model

fit_models<-function(list_formulas,data_training){

  #This functions several models and returns model objects.
  #Arguments: - list of formulas for GAM models
  #           - fitting data in a data.frame or SpatialPointDataFrame

  #Output: list of model objects 

  list_fitted_models<-vector("list",length(list_formulas))
  for (k in 1:length(list_formulas)){
    formula<-list_formulas[[k]]
    mod<- try(gam(formula, data=data_training)) #change to any model!!
    mod<- try(gam(formula, data=data_training,k=5)) #change to any model!!
    mod<- try(gam(y_var ~ s(lat, lon,k=14) + s(elev_s) + s(LST), data=data_training)) #change to any model!!

    #mod<- try(autoKrige(formula, input_data=data_s,new_data=s_sgdf,data_variogram=data_s))

    model_name<-paste("mod",k,sep="")
    assign(model_name,mod) 
    list_fitted_models[[k]]<-mod
  }
  return(list_fitted_models) 
}

mod1<- try(gam(y_var ~ s(lat, lon,k=22) + s(elev_s) + s(LST), data=data_training)) #change to any model!!

> sink("test.txt")
> (gam.check(mod1))
> getwd(0)
Error in getwd(0) : unused argument (0)
> getwd()
> sink()
> getwd(0)
Error in getwd(0) : unused argument (0)
> getwd()
[1] "/data/project/layers/commons/NEX_data/output_regions/15.0_20.0"
> system("more test.txt")

format_formula_k_fun <- function(formula,l_k){
  #l_k: list of k parameter for GAM method with the value of k and the name of the variable
  #split the formula by term 
  #add k parameter to formula
  names(l_k) <- attr(terms(formula),"term.labels")
  var_names <- all.vars(formula)
  term_names <- attr(terms(formula),"term.labels")
  term_names_tmp <- strsplit(term_names,")")
  term_names_tmp <- paste(term_names_tmp,", k=",l_k,")",sep="")
  term_names_tmp <- paste(term_names_tmp,collapse= " + ")
  formula_k <- as.formula(paste(var_names[1] ,"~",term_names_tmp))
  #nb_covar <- length(var_names) - 1
  #paste(var_names[2:(nb_covar+1)]),l_k)
  return(formula_k)
  #mod <- try(gam(y_var ~ s(lat, lon,k=22) + s(elev_s,k=10) + s(LST,k=10), data=data_training)) #change to any model!!
  #mod <- try(gam(formula_k, data=data_training)) #change to any model!!
}

test_k_gam <- function(formula,l_k){

  l_k_tmp <- as.numeric(l_k)/2
  #l_k_tmp <- l_k
  #maybe go in  a loop to  fit? could this in a loop for each term then select the term with highest k and k_index >1 ??
  #list_mod <- vector("list",0)
  #list_l_k <- vector("list",0)
  list_mod <- list()
  list_l_k <- list()

  j<- 1
  formula_k<-format_formula_k_fun(formula,l_k_tmp)
  mod <- try(gam(formula_k, data=data_training)) #change to any model!! 
  while(!inherits(mod,"try-error")){
    formula_k<-format_formula_k_fun(formula,l_k_tmp)
    mod <- try(gam(formula_k, data=data_training)) #change to any model!! 
    list_l_k[[j]] <- l_k_tmp
    list_mod[[j]] <- mod
    l_k_tmp <- l_k_tmp + rep(1,length(l_k_tmp))
    j <- j+1
    limit_l_k <- (l_k_tmp > l_k) #now check that c(30,10,10 is not gone above!!)
    for (i in 1:length(limit_l_k)){
      if (l_k_tmp[i] >= l_k[i]){
        l_k_tmp[i] <- l_k[i]
      }
    }
  }
  #store mod and l_k in a object
  l_k_obj <- list(list_mod,l_k_tmp)
  names(l_k_obj) <- c("list_mod","l_k_tmp")
  
  return(l_k_obj)
}

#d
#Now get k_index for each model and store it in a table.
#function gam.check with
#select the  right mode based on k-index, edf or other criteria?

l_k <- vector("list",3)
l_k <- list(25,10,10)
l_k <- c(30,10,10)


#1. Fit with standard settings:
#(Start at k=10 or k=30 if 2 interactive term.)
#start with the highest k first...
#If  does not fit then

#mod <- try(gam(formula_k, data=data_training)) #change to any model!!

#         K=15 for interactive term and k=5 unique term:
#         - do gam.check:
#                    If k-index < 1
#                     Then increase k by 1 for interactive term and  fit again
#          If k-index >1 or k=10 for unique term and k=30  
#Break…


#Now function to run mod depending on conditions