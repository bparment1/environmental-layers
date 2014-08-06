####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for experimentation of GAM model fitting ##############################
#This script uses the worklfow code applied to the globe. Results currently reside on NEX/PLEIADES NASA.
#The purpose is to test different parameters for the GAM function to allow .
## - experiment with gamma
# - experiment with k
# A table of diagnostic is created.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 07/30/2014  
#MODIFIED ON: 08/06/2014            
#Version: 1
#PROJECT: Environmental Layers project  
#TO DO:
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


### Function to fit using  the training  data  from the  workflow
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
    #This is where to change the function??
    model_name<-paste("mod",k,sep="")
    assign(model_name,mod) 
    list_fitted_models[[k]]<-mod
  }
  return(list_fitted_models) 
}

## Maybe write a new function to experiment gam fitting based on the above results??

##Function to predict using a mod object from the workflow...
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

## Create a new directory
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

#Extract info from object
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

### New functions to set the k dimension in GAM

##Format formula using prescribed k values for gam
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

## Fit model using training  data, formula and k dimensions with increment
test_k_gam <- function(formula,l_k,data_training){

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
    #if opt_increase{}
    l_k_tmp <- l_k_tmp + rep(1,length(l_k_tmp)) #can modify this option to go -1 instead of plus 1?
    #if opt_decrease{}
    l_k_tmp <- l_k_tmp + rep(1,length(l_k_tmp)) #can modify this option to go -1 instead of plus 1?
    
    j <- j+1
    limit_l_k <- (l_k_tmp > l_k) #now check that c(30,10,10 is not gone above!!)
    for (i in 1:length(limit_l_k)){
      if (l_k_tmp[i] >= l_k[i]){
        l_k_tmp[i] <- l_k[i]
      }
    }
  }
  #store mod and l_k in a object
  l_k_obj <- list(list_mod,list_l_k)
  names(l_k_obj) <- c("list_mod","list_l_k")
  
  return(l_k_obj)
}

## generate table of edf and k_index based on fitted models  
create_gam_check_table <-function(l_k_obj){
  list_mod <- l_k_obj$list_mod
  list_l_k <- l_k_obj$list_l_k

  #remove try-error model!!
  list_mod<- list_mod[unlist(lapply(list_mod,FUN=function(x){!inherits(x,"try-error")}))]
  list_df_k <- vector("list",length(list_mod))
  
  for(i in 1:length(list_mod)){
    mod <- list_mod[[i]]
    l_k <- list_l_k[[i]]
    sink("test.txt")
    gam.check(mod)
    sink()
    f <- readLines("test.txt") #read all lines
    f_tmp <- strsplit(f," ")
    names_term <- attr(terms(formula),"term.labels")
    match(f_tmp,names_term[1])
    grepl(glob2rx(paste(names_term[1],"*",sep="")),f_tmp)
    selected_lines_s <- grep("s\\(",f_tmp)
    selected_lines_ti <- grep("ti\\(",f_tmp)
    selected_lines_te <- grep("te\\(",f_tmp)    
    selected_lines <- sort(c(selected_lines_s,selected_lines_te,selected_lines_ti))
    #selected_lines <- c(min(selected_lines)-1,selected_lines)
    
    rows <- f[selected_lines]
    list_df <- vector("list",length(rows))
    for (j in  1:length(rows)){
      tmp <- unlist(strsplit(rows[j]," ",fixed=T))
      index <- unlist(lapply(tmp,FUN=function(x){x!=""}))
      tmp<- tmp[index]
      list_df[[j]] <- data.frame(term=tmp[1],kp=tmp[2],edf=tmp[3],k_index=tmp[4],p_value=tmp[5])
    }
    df_k <-do.call(rbind,list_df)
    df_k$k  <- l_k 
    df_k$fit_no <- paste("t_",i,sep="")
    list_df_k[[i]] <- df_k
  }
  return(list_df_k)
}

extract_fitting_mod_gam_stat <- function(mod){
  #Note that this assumes that we are using a gam mod object
  gcv_val <- mod$gcv.ubre
  aic_val <- AIC(mod)
  rmse_val <- sqrt(mean((residuals(mod))^2))
  mae_val <- mean(abs(residuals(mod)))
  bias_val <- mean(residuals(mod))
  n_val <- length(residuals(mod))
  #Now create a data.frame
  df_val <- data.frame(gcv=gcv_val,
                   aic=aic_val,
                   rmse=rmse_val,
                   mae=mae_val,
                   bias=bias_val,
                   n=n_val)  
  return(df_val)
} 

## Run whole process of diagnostic creation for given formula (without k), training data and model name
fit_gam_model_with_diagnostics <- function(l_k,data_training,formula,names_mod){
  #Input parameters:
  #l_k: vector of k value for each smooth term
  #data_training: data.frame containing  data for fitting of  gam
  #names_mod: character (string) with name of the model fitted.
  
  library(mgcv) #to avoid error below?? Maybe an indepent graphic devices needs to be opened??
  #Error in get(name, envir = asNamespace(pkg), inherits = FALSE) : 
  #object 'rversion' not found
  #Graphics error: Error in get(name, envir = asNamespace(pkg), inherits = FALSE) : 
  #object 'rversion' not found
  
  #STEP 1: fit with a range for k values
  
  #Fit models using given k values, formula and training dataset   
  #This is done starting at l_k/2 e.g. c(30,10,10) becomes c(15,5,5) where k is assigned to each term in the order given
  l_k_obj <- test_k_gam(formula,l_k,data_training)
  
  #STEP 2: get  k-index diagnostic for ech model
  
  #Now get k_index for each model and store it in a table.
  #function gam.check is sued
  
  list_df_k <- create_gam_check_table(l_k_obj)
  
  #STEP 3: produce additional metrics for diagnostis of fit
  #this includes the calculation of RMSE, MAE, bias, gcv,aic for fitted model
  list_mod <- l_k_obj$list_mod
  #remove try-error model!!
  list_mod<- list_mod[unlist(lapply(list_mod,FUN=function(x){!inherits(x,"try-error")}))]

  list_mod_gam_stat <- lapply(list_mod,FUN=extract_fitting_mod_gam_stat)
  
  #STEP 4: combine tables of diagnostics, format and add additional information  

  #combine tables...
  l_df_diagnostics <- lapply(1:length(list_df_k),FUN=function(i,df_1,df_2){cbind(df_1[[i]],df_2[[i]])},df_1=list_df_k,df_2=list_mod_gam_stat)
  df_diagnostics <- do.call(rbind,l_df_diagnostics)
  df_diagnostics$date <- rep(unique(data_training$date),nrow(df_diagnostics))
  df_diagnostics$month <- rep(unique(data_training$month))
  df_diagnostics$pred_mod <-names_mod #defined earlier... (input argutment )
  
  #Return the diagnostic table and model object
  names(list_mod) <- paste("t_",1:length(list_mod),sep="")
  diagnostics_obj <- list(df_diagnostics,list_mod)
  names(diagnostics_obj) <- c("df_diagnostics","list_mod")
  
  return(diagnostics_obj)
}

##############################s
#### Parameters and constants  

#scp -rp raster_prediction_obj_gam_CAI_dailyTmax15.0_20.0.RData parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/output_regions/15.0_20.0"
#scp -rp reg*.tif parmentier@atlas.nceas.ucsb.edu:/data/project/layers/commons/NEX_data/output_regions/15.0_20.0"

## laod data from Northen Africa
in_dir <- "/data/project/layers/commons/NEX_data/output_regions/15.0_20.0"
raster_obj_infile <- "raster_prediction_obj_gam_CAI_dailyTmax15.0_20.0.RData"

setwd(in_dir)

########################## START SCRIPT ##############################

########################################
#### PART I: Explore fitting of GAM with k and gamma parameters ####

raster_obj<- load_obj(raster_obj_infile)

#raster_obj <- load_obj(unlist(raster_obj_file)) #may not need unlist
nb_models <- length((raster_obj$clim_method_mod_obj[[1]]$formulas))
list_formulas <- (raster_obj$clim_method_mod_obj[[1]]$formulas)

#Models used:
#y_var ~ s(lat, lon) + s(elev_s)
#y_var ~ s(lat, lon) + s(elev_s) + s(LST)
#y_var ~ s(lat, lon) + s(elev_s) + s(N_w, E_w) + s(LST) + ti(LST,LC1) + s(LC1)

pred_mod <- paste("mod",c(1:nb_models,"_kr"),sep="")
#we are assuming no monthly hold out...
#we are assuming only one specific daily prop?
nb_models <- length(pred_mod)

#Select one month to play around:
j <- 7 # July
clim_method_mod_obj <- raster_obj$clim_method_mod_obj[[j]]
#this is made of "clim",data_month, data_month_v , sampling_month_dat, mod and formulas
clim_method_mod_obj$clim #file predicted

l_mod <- clim_method_mod_obj$mod #file predicted

#Quick look at the study area: Equatorial to Northern Africa
reg_rast <- stack(list.files(pattern="*.tif"))
plot(reg_rast,y=15)

names(clim_method_mod_obj)
clim_method_mod_obj$data_month
clim_method_mod_obj$data_month_v


k <- 2 #select model 2 with LST
formula <-list_formulas[[k]]
mod<- try(gam(formula, data=data_training)) #does not fit!! as expected

### TRY with different gamma

mod_g1<- try(gam(y_var ~ s(lat, lon) + s(elev_s),gamma=1.4, data=data_training)) #change to any model!!
gam.check(mod_g1)
#               k'    edf k-index p-value
#s(lat,lon) 29.000  9.470   0.967    0.32
#s(elev_s)   9.000  2.020   0.732    0.02

mod_g2<- try(gam(y_var ~ s(lat, lon) + s(elev_s),gamma=10, data=data_training)) #change to any model!!
gam.check(mod_g2) #increase in k-index!!

#              k'   edf k-index p-value
#s(lat,lon) 29.00 29.00    1.50    1.00
#s(elev_s)   9.00  9.00    1.11    0.74

mod<- try(gam(y_var ~ s(lat, lon) + s(elev_s) + s(LST),gamma=1.4, data=data_training)) #does not fit!
mod<- try(gam(y_var ~ s(lat, lon) + s(elev_s) + s(LST),gamma=10, data=data_training)) #does not fitl!!

### TRY with different k

mod_t1<- try(gam(y_var ~ s(lat, lon,k=14) + s(elev_s) + s(LST), data=data_training)) #change to any model!!
gam.check(mod_t1)
mod_t2<- try(gam(y_var ~ s(lat, lon,k=5) + s(elev_s) + s(LST), data=data_training)) #change to any model!!
gam.check(mod_t2) #in this case k=5 is too small for the interactive  term as k-index is less than 1

## check for residual pattern, removeable by increasing `k'
## typically `k', below, chould be substantially larger than 
## the original, `k' but certainly less than n/2.
vis.gam(mod_t1)
vis.gam(mod_t1,view=c("lat","lon"),theta= 35) # plot against by variable
#http://stats.stackexchange.com/questions/12223/how-to-tune-smoothing-in-mgcv-gam-model

mod_t3 <- try(gam(y_var ~ s(lat, lon,k=22) + s(elev_s) + s(LST), data=data_training)) #change to any model!!
gam.check(mod_t3)

#Explore mod object
mod_t1$gcv.ubre
mod_t1$aic
mod_t1$edf

########################################
#### PART II: Use the new functions to explore fitting with k-dimension in GAM  ####

j <- 7 # July
clim_method_mod_obj <- raster_obj$clim_method_mod_obj[[j]]
#this is made of "clim",data_month, data_month_v , sampling_month_dat, mod and formulas

l_mod <- clim_method_mod_obj$mod #file predicted
nb_models <- length((raster_obj$clim_method_mod_obj[[1]]$formulas))
list_formulas <- (raster_obj$clim_method_mod_obj[[1]]$formulas)
pred_mod <- paste("mod",c(1:nb_models,"_kr"),sep="")

#model_name<-paste("mod",k,sep="")
#we are assuming no monthly hold out...
#we are assuming only one specific daily prop?
nb_models <- length(pred_mod)

data_training <- clim_method_mod_obj$data_month

#list_fitted_models<-vector("list",length(list_formulas))
k <-2 #select model 3 with LST
names_mod <- paste("mod_",k,sep="")
model_name<-paste("mod",k,sep="")

formula <-list_formulas[[k]]

l_k <- c(30,10,10) #default values
#Create 
l_k_obj <- test_k_gam(formula,l_k,data_training)
#d
#Now get k_index for each model and store it in a table.
#function gam.check with

list_df_k <- create_gam_check_table(l_k_obj)

list_mod_gam_stat <- lapply(list_mod,FUN=extract_fitting_mod_gam_stat)
df_mod_stat <- list_mod_gam_stat[[7]]
#combine tables...
l_df_diagnostics<- lapply(1:length(list_df_k),FUN=function(i,df_1,df_2){cbind(df_1[[i]],df_2[[2]])},df_1=list_df_k,df_2=list_mod_gam_stat)

df_diagnostics <- do.call(rbind,l_df_diagnostics)

df_diagnostics$date <- rep(unique(data_training$date),nrow(df_diagnostics))
df_diagnostics$month <- rep(unique(data_training$month))

df_diagnostics$pred_mod <-names_mod #defined earlier... 

#select the  right mode based on k-index, edf or other criteria?

#choice of the lowest k for fit and k_index > 1
#then use model 1

#choice of the highest k for fit and k_index > 1
#then use model 7

########################################
#### PART III: Use the general functions to explore fitting with k-dimension in GAM  ####

#This can be done accross different months ...
#The general function provides a quick call to all function and formatting of diagnostic table
#first make a function for one model (could be month)

j <- 7 # July
clim_method_mod_obj <- raster_obj$clim_method_mod_obj[[j]]
#this is made of "clim",data_month, data_month_v , sampling_month_dat, mod and formulas
clim_method_mod_obj$clim #file predicted

l_mod <- clim_method_mod_obj$mod #file predicted

names(clim_method_mod_obj)
data_training <- clim_method_mod_obj$data_month
clim_method_mod_obj$data_month_v

list_fitted_models<-vector("list",length(list_formulas))
k <-3 #select model 3 with LST
formula <-list_formulas[[k]]

j <- 7 # July
clim_method_mod_obj <- raster_obj$clim_method_mod_obj[[j]]
#this is made of "clim",data_month, data_month_v , sampling_month_dat, mod and formulas

l_mod <- clim_method_mod_obj$mod #file predicted
nb_models <- length((raster_obj$clim_method_mod_obj[[1]]$formulas))
list_formulas <- (raster_obj$clim_method_mod_obj[[1]]$formulas)
pred_mod <- paste("mod",c(1:nb_models,"_kr"),sep="")
#we are assuming no monthly hold out...
#we are assuming only one specific daily prop?
nb_models <- length(pred_mod)

data_training <- clim_method_mod_obj$data_month

#list_fitted_models<-vector("list",length(list_formulas))
k <-2 #select model 2 with LST
#names_mod <- paste("mod_",k,sep="")
model_name<-paste("mod",k,sep="")

formula <-list_formulas[[k]]

l_k <- c(30,10,10) #default values

#test_df <- fit_gam_model_with_diagnostics(l_k,data_training,formula,model_name)
test_obj <- fit_gam_model_with_diagnostics(l_k,data_training,formula,model_name)

#Now do this over 12 months??

################## END OF SCRIPT ###############