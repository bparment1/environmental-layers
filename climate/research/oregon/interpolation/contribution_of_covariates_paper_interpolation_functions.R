####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for manuscript analyses,tables and figures #######################################
#This script reads information concerning the Oregon case study to adapt data for the revised 
# interpolation code.
#Figures and data for the contribution of covariate paper are also produced.
#AUTHOR: Benoit Parmentier                                                                      #
#DATE: 07/18/2014            
#Version: 2
#PROJECT: Environmental Layers project                                       #
#################################################################################################

###Loading R library and packages                                                      
library(gtools)                              # loading some useful tools 
library(mgcv)                                # GAM package by Simon Wood
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gstat)                               # Kriging and co-kriging by Pebesma et al.
library(fields)                              # NCAR Spatial Interpolation methods such as kriging, splines
library(raster)                              # Hijmans et al. package for raster processing
library(gdata)                               # various tools with xls reading
library(rasterVis)
library(parallel)
library(maptools)
library(maps)
library(reshape)
library(plotrix)
library(plyr)

#### FUNCTION USED IN SCRIPT

function_analyses_paper <-"contribution_of_covariates_paper_interpolation_functions_05252014.R"

load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

create_dir_fun <- function(out_dir,out_suffix){
  if(!is.null(out_suffix)){
    out_dir <- paste(out_dir,"_",out_suffix,sep="")
    #out_dir <- file.path(out_dir,out_name)
  }
  #create if does not exists
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  return(out_dir)
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
  #the output is a data.frame
  library(plyr)
  
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp<-obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]]<- as.data.frame(tmp) #deal with spdf cases
  }
  
  tb_list_tmp<-do.call(rbind.fill,list_tmp) #long rownames
  #tb_list_tmp<-do.call(rbind,list_tmp) #long rownames
  
  return(tb_list_tmp) #this is  a data.frame
}

add_rownames_list_df <-function(obj_list){
  #extract object from list of list. This useful for raster_prediction_obj
  #rownames are loaded in a column called rownames
  #the output is a data.frame
  library(plyr)
  
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    tmp <- obj_list[[i]] #double bracket to return data.frame
    rownames_char <- names(obj_list)[i]   
    rownames_col <- rep(rownames_char,nrow(tmp))
    tmp$rownames <- rownames_col
    list_tmp[[i]]<- as.data.frame(tmp) #deal with spdf cases
  }
  #tb_list_tmp<-do.call(rbind.fill,list_tmp) #long rownames
  #tb_list_tmp<-do.call(rbind,list_tmp) #long rownames
  return(list_tmp) #this is  a data.frame
}

#df must have a specific format: "%Y%m%d
add_month_tag<-function(tb,date_col){
  date<-strptime(tb[[date_col]], "%Y%m%d")   # interpolation date being processed
  month<-strftime(date, "%m")          # current month of the date being processed
}
  
calc_stat_from_raster_prediction_obj <-function(raster_prediction_obj,stat,training=FALSE){
  #Calculate statistics from validation and training out of raster_prediction_obj
  #If training is TRUE, then using training dataset
  
  #extract relevant information
  if(training==TRUE){
    tb <- extract_from_list_obj(raster_prediction_obj$validation_mod_obj,"metrics_s")
    rownames(tb)<-NULL #remove row names
  }else{
    tb <- raster_prediction_obj$tb_diagnostic_v
  }
  
  #Now summarize
  
  t<-melt(tb,
          measure=c("mae","rmse","r","me","m50"), 
          id=c("pred_mod"),
          na.rm=T)
  
  stat_tb<-cast(t,pred_mod~variable,stat)
  return(stat_tb)
}

## Concatenate datal.frames with numbers...
table_combined_symbol <-function(dfrm1,dfrm2,symbol_char){
  #concatenate element in a table format combining elements from 2 tables/data.frames
  #For instance an element can become: 2.3±0.32
  #Note that it is assumed that dfrm1 and dfrm2 have the same size/dimension!!
  
  rows<-nrow(dfrm1)
  cols<-ncol(dfrm1)
  #symbol_char <- "±" #set in arguments
  
  table_combined<-dfrm1
  
  for (i in 1:rows){
    for (j in 1:cols){
      table_combined[i,j]<-paste(dfrm1[i,j],symbol_char,dfrm2[i,j],sep="") #on ± windows xp char "\261", what is it on ubuntu?
    }
  }
  return(table_combined)
}

## Produce data.frame with residuals for models and distance to closest fitting station
calc_dist_ref_data_point <- function(i,list_param){
  #This function creates a list of data.frame containing the distance to teh closest
  # reference point (e.g. fitting station) for a give data frame. 
  #Inputs:
  #data_s: given data.frame from wich distance is computed
  #data_v: reference data.frame, destination, often the fitting points used in analyses
  #i: index variable to operate on list
  #names_var: 
  #Outputs:
  #list_dstspat_er
  
  #Parsing input arguments
  data_s<-list_param$data_s[[i]]
  data_v<-list_param$data_v[[i]]
  
  names_var<-list_param$names_var
  
  ######
  
  names_var<-intersect(names_var,names(data_v)) #there may be missing columns
  #use columns that exists
  
  d_s_v<-matrix(0,nrow(data_v),nrow(data_s))
  for(k in 1:nrow(data_s)){
    pt<-data_s[k,]
    d_pt<-(spDistsN1(data_v,pt,longlat=FALSE))/1000  #Distance to station k in km
    d_s_v[,k]<-d_pt
  }
  
  #Create data.frame with position, ID, dst and residuals...
  
  pos<-vector("numeric",nrow(data_v))
  y<-vector("numeric",nrow(data_v))
  dst<-vector("numeric",nrow(data_v))
  
  for (k in 1:nrow(data_v)){
    pos[k]<-match(min(d_s_v[k,]),d_s_v[k,])
    dst[k]<-min(d_s_v[k,]) 
  }
  
  dstspat_er<-as.data.frame(cbind(v_id=as.vector(data_v$id),s_id=as.vector(data_s$id[pos]),
                                  pos=pos, lat=data_v$lat, lon=data_v$lon, x=data_v$x,y=data_v$y,
                                  dst=dst,
                                  as.data.frame(data_v[,names_var])))
  
  return(dstspat_er)  
}  

### Main function to call to obtain distance to closest fitting stations for valiation dataset
distance_to_closest_fitting_station<-function(raster_prediction_obj,names_mod,dist_classes=c(0)){
  #This function computes the distance between training and testing points and returns and data frame
  #with distance,residuals, ID of training and testing points
  #Input: raster_prediction_obj, names of residuals models to return, distance classes
  #Output: data frame
  
  #Functions used in the script
  
  mae_fun<-function(x){mean(abs(x))} #Mean Absolute Error give a residuals vector
  sd_abs_fun<-function(x){sd(abs(x))} #sd Absolute Error give a residuals vector
  rmse_fun<-function(x){sqrt(mean(x^2))} #Mean Absolute Error give a residuals vector
  
  #calc_dist_ref_data_point : defined outside this function
  
  ##BEGIN
  
  ##### PART I: generate data.frame with residuals in term of distance to closest fitting station
  
  #return list of training and testing data frames
  list_data_s <- extract_list_from_list_obj(raster_prediction_obj$validation_mod_obj,"data_s") #training
  list_data_v <- extract_list_from_list_obj(raster_prediction_obj$validation_mod_obj,"data_v") #testing (validation)
  
  i<-1
  names_var<-c(names_mod,"dates")
  list_param_dst<-list(i,list_data_s,list_data_v,names_mod)
  names(list_param_dst) <- c("index","data_s","data_v","names_var")
  
  #call function "calc_dist_ref_data_point" over 365 dates
  #note that this function depends on other functions !!! see script
  
  list_dstspat_er <-lapply(1:length(list_data_v),FUN=calc_dist_ref_data_point,list_param=list_param_dst)
  #now assemble in one data.frame
  dstspat_dat<-do.call(rbind.fill,list_dstspat_er)

  ########### PART II: generate distance classes and summary statistics
  
  if (length(dist_classes)==1){
    range_val<-range(dstspat_dat$dst)
    max_val<-round_any(range_val[2],10, f=ceiling) #round max value to nearest 10 (from plyr package)
    min_val<-0
    limit_val<-seq(min_val,max_val, by=10)
  }else{
    limit_val<-dist_classes
  }
 
  dstspat_dat$dst_cat1 <- cut(dstspat_dat$dst,include.lowest=TRUE,breaks=limit_val)
  
  names_var <- intersect(names_mod,names(dstspat_dat)) #some of the models may not have been predictd so subselect
  t<-melt(dstspat_dat,
          measure=names_var, 
          id=c("dst_cat1"),
          na.rm=T)
  
  mae_tb <-cast(t,dst_cat1~variable,mae_fun)
  rmse_tb <-cast(t,dst_cat1~variable,rmse_fun)
  sd_abs_tb<-cast(t,dst_cat1~variable,sd_abs_fun)
  
  avg_tb<-cast(t,dst_cat1~variable,mean)
  sd_tb<-cast(t,dst_cat1~variable,sd)
  n_tb<-cast(t,dst_cat1~variable,length)
  #n_NA<-cast(t,dst_cat1~variable,is.na)
  
  #### prepare returning object
  dstspat_obj<-list(dstspat_dat,mae_tb,rmse_tb,sd_abs_tb,avg_tb,sd_tb,n_tb)
  names(dstspat_obj) <-c("dstspat_dat","mae_tb","rmse_tb","sd_abs_tb","avg_tb","sd_tb","n_tb")
  
  return(dstspat_obj)
  
}

# create plot of accuracy in term of distance to closest fitting station
plot_dst_spat_fun<-function(stat_tb,names_var,cat_val){
  
  range_y<-range(as.vector(unlist(stat_tb[,names_var])),na.rm=T) #flatten data.frame
  col_t<-rainbow(length(names_var))
  pch_t<- 1:length(names_var)
  plot(stat_tb[,names_var[1]], ylim=range_y,pch=pch_t[1],col=col_t[1],type="b",
       yla="MAE (in degree C)",xlab="",xaxt="n")
  #points((stat_tb[,names_var[k]], ylim=range_y,pch=pch_t[1],col=col_t[1]),type="p")
  for (k in 2:length(names_var)){
    lines(stat_tb[,names_var[k]], ylim=range_y,pch=pch_t[k],col=col_t[k],type="b",
          xlab="",axes=F)
    #points((stat_tb[,names_var[k]], ylim=range_y,pch=pch_t[k],col=col_t[k]),type="p")
  }
  legend("topleft",legend=names_var, 
         cex=1.2, pch=pch_t,col=col_t,lty=1,bty="n")
  axis(1,at=1:length(stat_tb[,1]),labels=stat_tb[,1])
}

plot_dst_MAE <-function(list_param){
  #
  #list_dist_obj: list of dist object 
  #col_t: list of color for each 
  #pch_t: symbol for line
  #legend_text: text for line and symbol
  #mod_name: selected models
  #
  ## BEGIN ##
  
  list_dist_obj<-list_param$list_dist_obj
  col_t<-list_param$col_t 
  pch_t<- list_param$pch_t 
  legend_text <- list_param$legend_text
  list_mod_name<-list_param$mod_name
  metric_name <-list_param$metric_name
  title_plot <- list_param$title_plot
  y_lab_text <- list_param$y_lab_text
  add_CI <- list_param$add_CI  
  
  for (i in 1:length(list_dist_obj)){
    
    l<-list_dist_obj[[i]]
    metric_tb<-l[[metric_name]]
    n_tb<-l$n_tb
    sd_abs_tb<-l$sd_abs_tb
    mod_name<-list_mod_name[i]
    #xlab_text<-"distance to fitting station (km) "
    
    n <- unlist(n_tb[,c(mod_name)])
    y <- unlist(metric_tb[,c(mod_name)])
    
    #x<- 1:length(y)
    x <- x_tick_labels
    y_sd <- unlist(sd_abs_tb[,c(mod_name)])
    
    ciw <-y_sd
    #ciw2   <- qt(0.975, n) * y_sd2 / sqrt(n)
    ciw   <- qt(0.975, n) * y_sd / sqrt(n)
    
    #if(i==1){
    #  plotCI(y=y, x=x, uiw=ciw, col=col_t[i], barcol="blue", lwd=1,
    #         ylab=y_lab_text, xlab="")
    #  lines(y~x, pch=pch_t[i],col=col_t[i],type="b")
    #  
    #}else{
    #  lines(y~x, pch=pch_t[i],col=col_t[i],type="b")
    #}
    
    if(add_CI[i]==TRUE){
      
      if (i==1){
        plotCI(y=y, x=x, uiw=ciw, col=col_t[i], barcol="blue", lwd=1,pch=pch_t[i],
               ylab=y_lab_text, xlab="") 
        lines(y~x, pch=pch_t[i],col=col_t[i],type="b")
      }else{
        lines(y~x, pch=pch_t[i],col=col_t[i],type="b")
        plotCI(y=y, x=x, uiw=ciw, col=col_t[i], barcol="blue", lwd=1,pch=pch_t[i],
               ylab=y_lab_text, xlab="",add=TRUE) 
      }
 
    }else{
      if (i==1){
        plot(y~x, pch=pch_t[i],col=col_t[i],type="b")
      }else{
        lines(y~x, pch=pch_t[i],col=col_t[i],type="b")
      }
    }
    
    
  }
  legend("topleft",legend=legend_text, 
         cex=1.2, pch=pch_t,col=col_t,lty=1,bty="n")
  #axis(1,at=1:length(stat_tb[,1]),labels=stat_tb[,1])
  title(title_plot)
  #paste(" Comparison of MAE in ",mod_name,sep="")
  
}

calc_stat_prop_tb <-function(names_mod,raster_prediction_obj,testing=TRUE){
  
  #add for testing??
  if (testing==TRUE){
    tb <-raster_prediction_obj$tb_diagnostic_v #use testing accuracy information
  }else{
    tb <-raster_prediction_obj$tb_diagnostic_s #use training accuracy information
  }
  
  t<-melt(subset(tb,pred_mod==names_mod),
          measure=c("mae","rmse","r","me","m50"), 
          id=c("pred_mod","prop"),
          na.rm=T)
  
  avg_tb<-cast(t,pred_mod+prop~variable,mean)
  sd_tb<-cast(t,pred_mod+prop~variable,sd)
  n_tb<-cast(t,pred_mod+prop~variable,length)
  #n_NA<-cast(t,dst_cat1~variable,is.na)
  
  #### prepare returning object
  prop_obj<-list(tb,avg_tb,sd_tb,n_tb)
  names(prop_obj) <-c("tb","avg_tb","sd_tb","n_tb")
  
  return(prop_obj)
}

#ploting 
plot_prop_metrics <-function(list_param){
  #
  #list_dist_obj: list of dist object 
  #col_t: list of color for each 
  #pch_t: symbol for line
  #legend_text: text for line and symbol
  #mod_name: selected models
  #
  ## BEGIN ##
  
  list_obj<-list_param$list_prop_obj
  col_t <-list_param$col_t 
  pch_t <- list_param$pch_t 
  legend_text <- list_param$legend_text
  list_mod_name<-list_param$mod_name
  metric_name<-list_param$metric_name
  add_CI <- list_param$add_CI
  CI_bar <- list_param$CI_bar
  
  for (i in 1:length(list_obj)){
    
    l<-list_obj[[i]]
    mod_name<-list_mod_name[i]
    avg_tb<-subset(l$avg_tb,pred_mod==mod_name,select=metric_name) #selecte relevant accuarcy metric
    n_tb<-subset(l$n_tb,pred_mod==mod_name,select=metric_name) 
    sd_tb<-subset(l$sd_tb,pred_mod==mod_name,select=metric_name) #l$sd_abs_tb[,metric_name]
    
    xlab_text<-""
    
    no <- unlist(as.data.frame(n_tb))
    y <- unlist(as.data.frame(avg_tb))
    
    x<- l$avg_tb$prop
    y_sd <- unlist(as.data.frame(sd_tb)) #sd_tb
    
    #ciw2   <- qt(0.975, n) * y_sd2 / sqrt(n)
    if (CI_bar[i]==TRUE){
      ciw   <- qt(0.975, no) * y_sd / sqrt(no)
    }else{
      ciw <-y_sd
    }
    
    #plotCI(y=y, x=x, uiw=ciw, col="red", main=paste(" MAE for ",mod_name,sep=""), barcol="blue", lwd=1,
    #       ylab="RMSE (C)", xlab=xlab_text)
    
    #ciw   <- qt(0.975, no) * y_sd / sqrt(no)
    
    #if(i==1){
    #  plotCI(y=y, x=x, uiw=ciw, col=col_t[i], main=paste(" Comparison of ",metric_name," in ",mod_name,sep=""), barcol="blue", lwd=1,
    #         ylab="RMSE (C)", xlab=xlab_text)
    #  lines(y~x, col=col_t[i])
      
    #}else{
    #  lines(y~x, col=col_t[i])
    #}
    
    if(add_CI[i]==TRUE){
      
      if (i==1){
        plotCI(y=y, x=x, uiw=ciw, col=col_t[i], barcol="blue", lwd=1,pch=pch_t[i],
               ylab="", xlab=xlab_text) 
        lines(y~x, pch=pch_t[i],col=col_t[i],type="b")
      }else{
        lines(y~x, pch=pch_t[i],col=col_t[i],type="b")
        plotCI(y=y, x=x, uiw=ciw, col=col_t[i], barcol="blue", lwd=1,pch=pch_t[i],
               ylab="", xlab=xlab_text,add=TRUE) 
      }
      
    }else{
      if (i==1){
        plot(y~x, pch=pch_t[i],col=col_t[i],type="b")
      }else{
        lines(y~x, pch=pch_t[i],col=col_t[i],type="b")
      }
    }
    
    
  }
  legend("topleft",legend=legend_text, 
         cex=1.2, pch=pch_t,col=col_t,lty=1,bty="n")
  #axis(1,at=1:length(stat_tb[,1]),labels=stat_tb[,1])
  
}

plot_dst_spat_fun<-function(stat_tb,names_var,cat_val){
  
  range_y<-range(as.vector(unlist(stat_tb[,names_var])),na.rm=T) #flatten data.frame
  col_t<-rainbow(length(names_var))
  pch_t<- 1:length(names_var)
  plot(stat_tb[,names_var[1]], ylim=range_y,pch=pch_t[1],col=col_t[1],type="b",
       yla="MAE (in degree C)",xlab="",xaxt="n")
  #points((stat_tb[,names_var[k]], ylim=range_y,pch=pch_t[1],col=col_t[1]),type="p")
  for (k in 2:length(names_var)){
    lines(stat_tb[,names_var[k]], ylim=range_y,pch=pch_t[k],col=col_t[k],type="b",
          xlab="",axes=F)
    #points((stat_tb[,names_var[k]], ylim=range_y,pch=pch_t[k],col=col_t[k]),type="p")
  }
  legend("topleft",legend=names_var, 
         cex=1.2, pch=pch_t,col=col_t,lty=1,bty="n")
  axis(1,at=1:length(stat_tb[,1]),labels=stat_tb[,1])
}

#Calculate the difference between training and testing in two different data.frames. Columns to substract are provided.
diff_df<-function(tb_s,tb_v,list_metric_names){
  tb_diff<-vector("list", length(list_metric_names))
  for (i in 1:length(list_metric_names)){
    metric_name<-list_metric_names[i]
    tb_diff[[i]] <-tb_s[,c(metric_name)] - tb_v[,c(metric_name)]
  }
  names(tb_diff)<-list_metric_names
  tb_diff<-as.data.frame(do.call(cbind,tb_diff))
  return(tb_diff)
}

create_s_and_p_table_term_models <-function(i,list_myModels){
  #Purpose:
  #Function to extract smooth term  table,parameter table and AIC from a list of models.
  #Originally created to processed GAM models run over a full year.
  #Inputs: 
  #1) list_myModels: list of fitted GAM models
  #2) i: index of list to run with lapply or mcapply
  #Outputs: list of 
  #1)s.table.term
  #2)p.table.term
  #3)AIC list
  #Authors: Benoit Parmentier
  #date: 08/142013
  
  ### Functions used in the scritp:
  
  #Remove models that were not fitted from the list
  #All modesl that are "try-error" are removed
  remove_errors_list<-function(list_items){
    #This function removes "error" items in a list
    list_tmp<-list_items
    for(i in 1:length(list_items)){
      
      if(inherits(list_items[[i]],"try-error")){
        list_tmp[[i]]<-0
      }else{
        list_tmp[[i]]<-1
      }
    }
    cnames<-names(list_tmp[list_tmp>0])
    x<-list_items[match(cnames,names(list_items))]
    return(x)
  }
  
  #turn term from list into data.frame
  name_col<-function(i,x){
    x_mat<-x[[i]]
    x_df<-as.data.frame(x_mat)
    x_df$mod_name<-rep(names(x)[i],nrow(x_df))
    x_df$term_name <-row.names(x_df)
    return(x_df)
  }
  
  ##BEGIN ###
  
  myModels <- list_myModels[[i]]
  myModels<-remove_errors_list(myModels)
  #could add AIC, GCV to the table as well as ME, RMSE...+dates...
  
  summary_list <- lapply(myModels, summary)
  s.table_list <- lapply(summary_list, `[[`, 's.table')
  p.table_list <- lapply(summary_list, `[[`, 'p.table')
  AIC_list <- lapply(myModels, AIC)
  
  #now put in one table
  
  s.table_list2<-lapply(1:length(myModels),name_col,s.table_list)
  p.table_list2<-lapply(1:length(myModels),name_col,p.table_list)
  s.table_term <-do.call(rbind,s.table_list2)
  p.table_term <-do.call(rbind,p.table_list2)
  
  #Now get AIC
  AIC_list2<-lapply(1:length(myModels),name_col,AIC_list)
  AIC_models <- do.call(rbind,AIC_list2)
  names(AIC_models)[1]<-"AIC"
  
  #Set up return object
  
  s_p_table_term_obj<-list(s.table_term,p.table_term,AIC_models)
  names(s_p_table_term_obj) <-c("s.table_term","p.table_term","AIC_models")
  return(s_p_table_term_obj)
  
}

convert_spdf_to_df_from_list <-function(obj_list,list_name){
  #extract object from list of list. This useful for raster_prediction_obj
  #output: data.frame formed by rbinding sp data.frame in the list
  library(plyr)
  
  list_tmp<-vector("list",length(obj_list))
  for (i in 1:length(obj_list)){
    #tmp<-obj_list[[i]][[list_name]] #double bracket to return data.frame
    list_tmp[[i]]<-as.data.frame(obj_list[[i]])
  }
  tb_list_tmp<-do.call(rbind.fill,list_tmp) #long rownames
  #tb_list_tmp<-do.call(rbind,list_tmp) #long rownames
  
  return(tb_list_tmp) #this is  a data.frame
}

calc_stat_month_from_raster_prediction_obj <-function(raster_prediction_obj,stat){
  #Calculate monthly averages from tb
  #Function
  #tb must have a month column month
  add_month_tag<-function(tb){
    date<-strptime(tb$date, "%Y%m%d")   # interpolation date being processed
    month<-strftime(date, "%m")          # current month of the date being processed
  }
  
  ## Beging ##
  tb<-raster_prediction_obj$tb_diagnostic_v
  tb$month<-add_month_tag(tb)
  
  t<-melt(tb,
          measure=c("mae","rmse","r","me","m50"), 
          id=c("pred_mod","month"),
          na.rm=T)
  
  stat_tb<-cast(t,pred_mod+month~variable,stat)
  return(stat_tb)
}

### generate filter for Moran's I function in raster package
autocor_filter_fun <-function(no_lag=1,f_type="queen"){
  if(f_type=="queen"){
    no_rows <- 2*no_lag +1
    border_row <-rep(1,no_rows)
    other_row <- c(1,rep(0,no_rows-2),1)
    other_rows <- rep(other_row,no_rows-2)
    mat_data<- c(border_row,other_rows,border_row)
    autocor_filter<-matrix(mat_data,nrow=no_rows)
  }
  #if(f_type=="rook){} #add later
  return(autocor_filter)
}

#Now run Moran's I for raster image given a list of  filters for different lags and raster stack
moran_multiple_fun<-function(i,list_param){
  #Parameters:
  #list_filters: list of filters with different lags in the image
  #r_stack: stack of raster image, only the selected layer is used...
  list_filters <-list_param$list_filters
  r <- subset(list_param$r_stack,i)
  moran_list <- lapply(list_filters,FUN=Moran,x=r)
  moran_v <-as.data.frame(unlist(moran_list))
  names(moran_v)<-names(r)
  return(moran_v)
}


################### END OF SCRIPT ###################


