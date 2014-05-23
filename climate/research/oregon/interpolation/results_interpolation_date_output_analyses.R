##################    Validation and analyses of results  #######################################
############################ Covariate production for a given tile/region ##########################################
#This script examines inputs and outputs from the interpolation step.                             
#Part 1: Script produces plots for every selected date
#Part 2: Examine 
#AUTHOR: Benoit Parmentier                                                                       
#DATE: 08/05/2013                                                                                 
#DATE MODIFIED: 05/21/2014                                                                                 

#PROJECT: NCEAS INPLANT: Environment and Organisms --TASK#???--   

##################################################################################################

## Function(s) used in script

load_obj <- function(f) 
{
  env <- new.env()
  nm <- load(f, env)[1]  
  env[[nm]]
}

#Also used in validation script...place new files for common functions used in interpolation

extract_from_list_obj<-function(obj_list,list_name){
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp<-obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]]<-tmp
  }
  tb_list_tmp<-do.call(rbind,list_tmp) #long rownames
  return(tb_list_tmp) #this is  a data.frame
}

plots_assessment_by_date<-function(j,list_param){
  ###Function to assess results from interpolation predictions
  #AUTHOR: Benoit Parmentier                                                                       
  #DATE: 08/05/2013                                                                                 
  #PROJECT: NCEAS INPLANT: Environment and Organisms --TASK#363--   
  
  #1) in_path_tile: location of files, if NULL then code is run on NEX node or as stage 5
  #2) out_path
  #3) script_path
  #4) raster_prediction_obj
  #5) interpolation_method
  #7) covar_obj: covariates object contains file name and names of covariates
  #8) date_selected_results
  #9) var
  #10) out_prefix
  
  ###Loading R library and packages   
  library(RPostgreSQL)
  library(sp)                                             # Spatial pacakge with class definition by Bivand et al.
  library(spdep)                                          # Spatial pacakge with methods and spatial stat. by Bivand et al.
  library(rgdal)                                          # GDAL wrapper for R, spatial utilities
  library(raster)
  library(gtools)
  library(rasterVis)
  library(graphics)
  library(grid)
  library(lattice)
  
  ### BEGIN SCRIPT
  #Parse input parameters
  #date_selected_results <- c("x")
  date_selected<-list_param$date_selected_results #dates for plot creation
  var<-list_param$var #variable being interpolated
  out_path <- list_param$out_path
  interpolation_method <- list_param$interpolation_method

  in_path_tile <- list_param$in_path_tile
  
  if(!is.null(in_path_tile)){
    covar_obj <- load_obj(list_param$covar_obj)
    infile_covariates <- file.path(in_path_tile,basename(covar_obj$infile_covariates))
    covar_names <- covar_obj$covar_names
  }else{ #we are on the node or running as stage 5
    infile_covariates <- list_param$covar_obj$infile_covariates #already loaded in memory
    covar_names<-list_param$covar_obj$covar_names
  }
  
  #if raster_obj has not been loaded in memory then we have
  #the name of the RData object for a specific tile
  raster_prediction_obj<-list_param$raster_prediction_obj
  if(class(raster_prediction_obj)=="character"){
    raster_prediction_obj <- load_obj(raster_prediction_obj)
  }

  method_mod_obj <- raster_prediction_obj$method_mod_obj
  validation_mod_obj <-raster_prediction_obj$validation_mod_obj
  
  if(interpolation_method %in% c("gam_CAI","kriging_CAI","gwr_CAI","gam_fusion","kriging_fusion","gwr_fusion")){
    multi_timescale <- TRUE
  }
  #This should not be set here...? master script
  if (var=="TMAX"){
    y_var_name<-"dailyTmax"
    y_var_month<-"TMax"
  }
  if (var=="TMIN"){
    y_var_name<-"dailyTmin"
    y_var_month <-"TMin"
  }
  #add precip option later...
  
  ## Read covariate brick...
  s_raster <- brick(infile_covariates) #stack produced for specific tile
  names(s_raster)<- covar_names               #Assigning names to the raster layers: making sure it is included in the extraction
  
  ## Prepare study area  mask: based on LC12 (water)
  
  LC_mask <- subset(s_raster,"LC12")
  LC_mask[LC_mask==100]<- NA
  LC_mask <- LC_mask < 100
  LC_mask_rec<-LC_mask
  LC_mask_rec[is.na(LC_mask_rec)]<- 0
    
  #determine index position matching date selected
  i_dates<-vector("list",length(date_selected))
  for (j in 1:length(date_selected)){
    for (i in 1:length(method_mod_obj)){
      if(method_mod_obj[[i]]$sampling_dat$date==date_selected[j]){  
        i_dates[[j]]<-i
      }
    }
  }
  #Examine the first select date add loop or function later
  #j=1
  date <- strptime(date_selected[j], "%Y%m%d")   # interpolation date being processed
  month <- strftime(date, "%m")          # current month of the date being processed
  
  #Get raster stack of interpolated surfaces
  index <- i_dates[[j]]
  ##The path of production is not the same if input_path_tile is not NULL
  if(!is.null(in_path_tile)){
    #infile_covariates <- file.path(in_path_tile,basename(list_param$covar_obj$infile_covariates))
    pred_temp <- basename(as.character(method_mod_obj[[index]][[y_var_name]])) #list of daily prediction files with path included
    pred_temp <- file.path(in_path_tile,pred_temp)
  }else{
    pred_temp <- as.character(method_mod_obj[[index]][[y_var_name]]) #list of daily prediction files with path included
  }

  rast_pred_temp_s <- stack(pred_temp) #stack of temperature predictions from models (daily)
  rast_pred_temp <- mask(rast_pred_temp_s,LC_mask,file=file.path(out_path,"test.tif"),overwrite=TRUE)
  
  #Get validation metrics, daily spdf training and testing stations, monthly spdf station input
  sampling_dat<-method_mod_obj[[index]]$sampling_dat
  metrics_v<-validation_mod_obj[[index]]$metrics_v
  metrics_s<-validation_mod_obj[[index]]$metrics_s
  data_v<-validation_mod_obj[[index]]$data_v #testing with residuals
  data_s<-validation_mod_obj[[index]]$data_s
  #no formula if multi-timescale method
  if(multi_timescale==TRUE){
    formulas<- raster_prediction_obj$clim_method_mod_obj[[as.integer(month)]]$formulas #models ran
  }else{
    formulas<- method_mod_obj[[index]]$formulas #models ran
  }
  
  #Adding layer LST to the raster stack of covariates
  #The names of covariates can be changed...
  
  LST_month<-paste("mm_",month,sep="") # name of LST month to be matched
  pos<-match("LST",names(s_raster)) #Find the position of the layer with name "LST", if not present pos=NA
  s_raster<-dropLayer(s_raster,pos)      # If it exists drop layer
  pos<-match(LST_month,names(s_raster)) #Find column with the current month for instance mm12
  r1<-raster(s_raster,layer=pos)             #Select layer from stack
  names(r1)<-"LST"
  #Get mask image!!
  
  date_proc<-strptime(sampling_dat$date, "%Y%m%d")   # interpolation date being processed
  mo<-as.integer(strftime(date_proc, "%m"))          # current month of the date being processed
  day<-as.integer(strftime(date_proc, "%d"))
  year<-as.integer(strftime(date_proc, "%Y"))
  datelabel=format(ISOdate(year,mo,day),"%b %d, %Y")
  
  ## Figure 1: LST_TMax_scatterplot 
  
  rmse<-metrics_v$rmse[nrow(metrics_v)]
  rmse_f<-metrics_s$rmse[nrow(metrics_s)]  
  
  #Set as constant in master script ?? : c("gam_CAI","kriging_CAI","gwr_CAI","gam_fusion","kriging_fusion","gwr_fusion")
  #if (interpolation_method %in% c("gam_CAI","kriging_CAI","gwr_CAI","gam_fusion","kriging_fusion","gwr_fusion")){
  if (multi_timescale==TRUE){
    #if multi-time scale method then the raster prediction object contains a "clim_method_mod_obj"
    clim_method_mod_obj <- raster_prediction_obj$clim_method_mod_obj
    data_month <-clim_method_mod_obj[[mo]]$data_month
    
    png(file.path(out_path,paste("LST_",y_var_month,"_scatterplot_",sampling_dat$date,"_",sampling_dat$prop,"_",sampling_dat$run_samp,
                                 out_prefix,".png", sep="")))
    plot(data_month[[y_var_month]],data_month$LST,xlab=paste("Station mo ",y_var_month,sep=""),ylab=paste("LST month ",mo," ",sep=""))
    title(paste("LST vs ", y_var_month,"for month ",mo,sep=" "))
    abline(0,1)
    nb_point<-paste("n=",length(data_month[[y_var_month]]),sep="")
    LSTD_bias <- data_month$TMax - data_month$LST #in case it is a CAI method, calculate bias
    mean_bias<-paste("Mean LST bias= ",format(mean(LSTD_bias,na.rm=TRUE),digits=3),sep="")
    #Add the number of data points on the plot
    legend("topleft",legend=c(mean_bias,nb_point),bty="n")
    dev.off()
    
    ## Figure 2: Daily_tmax_monthly_TMax_scatterplot, modify for TMin!!
    #This is not stored in data_s$TMax?
    #png(file.path(out_path,paste("Month_day_scatterplot_",y_var_name,"_",y_var_month,"_",sampling_dat$date,"_",sampling_dat$prop,"_",sampling_dat$run_samp,
    #                             out_prefix,".png", sep="")))
    #plot(data_s[[y_var_name]]~data_s[[y_var_month]],xlab=paste("Month") ,ylab=paste("Daily for",datelabel),main="across stations in VE")
    #nb_point<-paste("ns=",length(data_s[[y_var_month]]),sep="")
    #nb_point2<-paste("ns_obs=",length(data_s[[y_var_month]])-sum(is.na(data_s[[y_var_name]])),sep="")
    #nb_point3<-paste("n_month=",length(data_month[[y_var_month]]),sep="")
    #Add the number of data points on the plot
    #legend("topleft",legend=c(nb_point,nb_point2,nb_point3),bty="n",cex=0.8)
    #dev.off()
    
    ## Figure 3: monthly stations used
    
    png(file.path(out_path,paste("Monthly_data_study_area_", y_var_name,
                                 out_prefix,".png", sep="")))
    plot(raster(rast_pred_temp,layer=1))
    plot(data_month,col="black",cex=1.2,pch=4,add=TRUE)
    title("Monthly ghcn station in tile for January")
    dev.off()
  
  } #End of if multi_timescale=TRUE
  
  ## Figure 4: Predicted_tmax_versus_observed_scatterplot 
  
  names_mod <- names(method_mod_obj[[index]][[y_var_name]]) #names of models to plot
  #model_name<-"mod_kr" #can be looped through models later on...
  
  for (k in 1:length(names_mod)){
    model_name <- names_mod[k]
    png(file.path(out_path,paste("Predicted_versus_observed_scatterplot_",y_var_name,"_",model_name,"_",sampling_dat$date,"_",sampling_dat$prop,"_",
                                 sampling_dat$run_samp,out_prefix,".png", sep="")))
    y_range<-range(c(data_s[[model_name]],data_v[[model_name]]),na.rm=T)
    x_range<-range(c(data_s[[y_var_name]],data_v[[y_var_name]]),na.rm=T)
    col_t<- c("black","red")
    pch_t<- c(1,2)
    plot(data_s[[model_name]],data_s[[y_var_name]], 
         xlab=paste("Actual daily for",datelabel),ylab="Pred daily", 
         ylim=y_range,xlim=x_range,col=col_t[1],pch=pch_t[1])
    points(data_v[[model_name]],data_v[[y_var_name]],col=col_t[2],pch=pch_t[2])
    grid(lwd=0.5, col="black")
    abline(0,1)
    legend("topleft",legend=c("training","testing"),pch=pch_t,col=col_t,bty="n",cex=0.8)
    title(paste("Predicted_versus_observed_",y_var_name,"_",model_name,"_",datelabel,sep=" "))
    nb_point1<-paste("ns_obs=",length(data_s[[y_var_name]])-sum(is.na(data_s[[model_name]])),sep="")
    nb_point2<-paste("nv_obs=",length(data_v[[y_var_name]])-sum(is.na(data_v[[model_name]])),sep="")
    
    rmse_str1<-paste("RMSE= ",format(rmse,digits=3),sep="")
    rmse_str2<-paste("RMSE_f= ",format(rmse_f,digits=3),sep="")
    
    #Add the number of data points on the plot
    legend("bottomright",legend=c(nb_point1,nb_point2,rmse_str1,rmse_str2),bty="n",cex=0.8)
    dev.off()
  }
  
  ## Figure 5a: prediction raster images
  png(file.path(out_path,paste("Raster_prediction_levelplot_",y_var_name,"_",sampling_dat$date,"_",sampling_dat$prop,"_",sampling_dat$run_samp,
            out_prefix,".png", sep="")))
  #paste(metrics_v$pred_mod,format(metrics_v$rmse,digits=3),sep=":")
  names(rast_pred_temp)<-paste(metrics_v$pred_mod,format(metrics_v$rmse,digits=3),sep=":")
  #plot(rast_pred_temp)
  levelplot(rast_pred_temp) #not working...takes too long to plot?
  dev.off()
  
  ## Figure 5b: prediction raster images
  png(file.path(out_path,paste("Raster_prediction_plot_",y_var_name,"_",sampling_dat$date,"_",sampling_dat$prop,"_",sampling_dat$run_samp,
            out_prefix,".png", sep="")))
  #paste(metrics_v$pred_mod,format(metrics_v$rmse,digits=3),sep=":")
  names(rast_pred_temp)<-paste(metrics_v$pred_mod,format(metrics_v$rmse,digits=3),sep=":")
  plot(rast_pred_temp)
  dev.off()
  
  ## Figure 6: training and testing daily stations used
  png(file.path(out_path,paste("Training_testing_stations_map_",y_var_name,"_",sampling_dat$date,"_",sampling_dat$prop,"_",sampling_dat$run_samp,
            out_prefix,".png", sep="")))
  plot(raster(rast_pred_temp,layer=1))
  plot(data_s,col="black",cex=1.2,pch=2,add=TRUE)
  plot(data_v,col="red",cex=1.2,pch=1,add=TRUE)
  legend("topleft",legend=c("training stations", "testing stations"), 
         cex=1, col=c("black","red"),
         pch=c(2,1),bty="n")
  title(paste("Daily stations ", datelabel,sep=""))
  nb_point1<-paste("ns_obs=",nrow(data_s),sep="")
  nb_point2<-paste("nv_obs=",nrow(data_v),sep="")
  legend("bottomright",legend=c(nb_point1,nb_point2),bty="n",cex=0.8)

  dev.off()
  
  ## Figure 7: delta surface and bias
  
  if (interpolation_method%in% c("gam_fusion","kriging_fusion","gwr_fusion")){
    png(file.path(out_path,paste("Bias_delta_surface_",y_var_name,"_",sampling_dat$date[i],"_",sampling_dat$prop[i],
              "_",sampling_dat$run_samp[i],out_prefix,".png", sep="")))
      ##The path of production is not the same if input_path_tile is not NULL
    if(!is.null(in_path_tile)){
      #infile_covariates <- file.path(in_path_tile,basename(list_param$covar_obj$infile_covariates))
      bias_lf <- basename(as.character(clim_method_mod_obj[[index]]$bias)) #list of daily prediction files with path included
      bias_lf <- file.path(in_path_tile,bias_lf)
      delta_lf <- basename(unlist(method_mod_obj[[index]]$delta))
      delta_lf <- file.path(in_path,delta_lf)
    }else{
      bias_lf <- clim_method_mod_obj[[index]]$bias #list of daily prediction files with path included
      delta_lf <- method_mod_obj[[index]]$delta
    }

    bias_rast<-stack(bias_lf)
    delta_rast<-raster(delta_lf) #only one delta image!!!
    names(delta_rast)<-"delta"
    rast_temp_date<-stack(bias_rast,delta_rast)
    layers_names <- names(rast_temp_date)
    rast_temp_date<-mask(rast_temp_date,LC_mask,file=file.path(out_path,"test.tif"),overwrite=TRUE)
    #bias_d_rast<-raster("fusion_bias_LST_20100103_30_1_10d_GAM_fus5_all_lstd_02082013.rst")
    names(rast_temp_date) <-layers_names
    plot(rast_temp_date)
    dev.off()
  }
  #if CAI
  if (interpolation_method %in% c("gam_CAI","kriging_CAI","gwr_CAI","gam_fusion","kriging_fusion","gwr_fusion")){
    png(file.path(out_path,paste("clim_surface_",y_var_name,"_",sampling_dat$date[i],"_",sampling_dat$prop[i],
              "_",sampling_dat$run_samp[i],out_prefix,".png", sep="")))
      ##The path of production is not the same if input_path_tile is not NULL
    if(!is.null(in_path_tile)){
      #infile_covariates <- file.path(in_path_tile,basename(list_param$covar_obj$infile_covariates))
      clim_lf <- basename(as.character(clim_method_mod_obj[[mo]]$clim)) #list of daily prediction files with path included
      clim_lf <- file.path(in_path_tile,clim_lf)
      delta_lf <- basename(unlist(method_mod_obj[[index]]$delta))
      delta_lf <- file.path(in_path,delta_lf)
    }else{
      clim_lf <- clim_method_mod_obj[[index]]$clim #list of monthly prediction files with path included
      delta_lf <- method_mod_obj[[index]]$delta
    }    
    clim_rast<-stack(clim_lf)
    delta_rast<-stack(delta_lf) #this is a stack now... delta images!!!
     
    names(delta_rast)<- paste(names_mod,"_delta",sep="")
    names(clim_rast) <- paste(names_mod,"_month",mo,sep="")
    #layers_names <- c(names(clim_rast),"delta")
    #rast_temp_date<-stack(clim_rast,delta_rast)
    #rast_temp_date<-mask(rast_temp_date,LC_mask,file=file.path(out_path,"test.tif"),overwrite=TRUE) #loosing names here
    #bias_d_rast<-raster("fusion_bias_LST_20100103_30_1_10d_GAM_fus5_all_lstd_02082013.rst")
    #names(rast_temp_date) <-layers_names
    #plot(rast_temp_date)
    plot(clim_rast)
    #title("Climatology for month ", mo, sep="")

    dev.off()
    
    png(file.path(out_path,paste("delta_surface_",y_var_name,"_",sampling_dat$date[i],"_",sampling_dat$prop[i],
              "_",sampling_dat$run_samp[i],out_prefix,".png", sep="")))
    plot(delta_rast) 
    dev.off()
  }
  
  ### Figure 9: map of residuals...

  for (k in 1:length(names_mod)){
    model_name <- names_mod[k]
    #fig_name <- file.path(out_path,paste("Figure_residuals_map_",y_var_name,"_",model_name,"_",sampling_dat$date,"_",sampling_dat$prop,"_",
    #                             sampling_dat$run_samp,out_prefix,".png", sep=""))
    
    png(file.path(out_path,paste("Figure_residuals_map_",y_var_name,"_",model_name,"_",sampling_dat$date,"_",sampling_dat$prop,"_",
                                 sampling_dat$run_samp,out_prefix,".png", sep="")))
    res_model_name <- paste("res",model_name,sep="_")
    elev <- subset(s_raster,"elev_s")
    p1 <- levelplot(elev,scales = list(draw = FALSE), colorkey = FALSE,col.regions=rev(terrain.colors(255)))
    #add legend..
    cx <- ((data_v[[res_model_name]])^2)/10
    p2 <- spplot(data_v,res_model_name, cex=1 * cx,main=paste("Residuals for ",res_model_name," ",datelabel,sep=""))
    p3 <-p2+p1+p2 #to force legend...
    #p2
    print(p3)
    dev.off()
  }
  
  #Figure 9: histogram for all images/residuals...
  ## Add later...? distance to closest fitting station?
  
  #tb_diagnostic_v <- raster_prediction_obj$tb_diagnostic_v
  #raster_prediction_obj$summary_metrics_v
  #raster_prediction_obj$summary_month_metrics_v$metric_month_avg
  #raster_prediction_obj$summary_month_metrics_v$metric_month_sd
  
  #Write out accuracy information:
  
  #add sd later...
  #write.table(tb_diagnostic_v,)
  #write.table(tb_diagnostic_v, file= file.path(out_path,interpolation_method,"_tb_diagnostic_v",out_prefix,".txt",sep=""), sep=",",overwrite=FALSE)
  #write.table(tb, file= paste(path,"/","results2_gwr_Assessment_measure_all",out_prefix,".txt",sep=""), sep=",")
  
  #histogram(rast_pred_temp)
  list_output_analyses<-list(metrics_s,metrics_v)
  names(list_output_analyses) <- c("metrics_s", "metrics_v")
  return(list_output_analyses)
  
}

