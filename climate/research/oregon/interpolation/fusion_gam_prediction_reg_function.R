runGAMFusion <- function(i) {            # loop over dates
  
  date<-strptime(dates[i], "%Y%m%d")   # interpolation date being processed
  month<-strftime(date, "%m")          # current month of the date being processed
  LST_month<-paste("mm_",month,sep="") # name of LST month to be matched
  
  #Adding layer LST to the raster stack
  
  pos<-match(LST_month,layerNames(s_raster)) #Find column with the current month for instance mm12
  r1<-raster(s_raster,layer=pos)             #Select layer from stack
  layerNames(r1)<-"LST"
  s_raster<-addLayer(s_raster,r1)            #Adding current month
  
  ###Regression part 1: Creating a validation dataset by creating training and testing datasets
  
  mod_LST <-ghcn.subsets[[i]][,match(LST_month, names(ghcn.subsets[[i]]))]  #Match interpolation date and monthly LST average
  ghcn.subsets[[i]] = transform(ghcn.subsets[[i]],LST = mod_LST)            #Add the variable LST to the subset dataset
  #n<-nrow(ghcn.subsets[[i]])
  #ns<-n-round(n*prop)   #Create a sample from the data frame with 70% of the rows
  #nv<-n-ns              #create a sample for validation with prop of the rows
  #ind.training <- sample(nrow(ghcn.subsets[[i]]), size=ns, replace=FALSE) #This selects the index position for 70% of the rows taken randomly
  ind.training<-sampling[[i]]
  ind.testing <- setdiff(1:nrow(ghcn.subsets[[i]]), ind.training)
  data_s <- ghcn.subsets[[i]][ind.training, ]   #Training dataset currently used in the modeling
  data_v <- ghcn.subsets[[i]][ind.testing, ]    #Testing/validation dataset using input sampling
  
  ns<-nrow(data_s)
  nv<-nrow(data_v)
  #i=1
  date_proc<-dates[i]
  date_proc<-strptime(dates[i], "%Y%m%d")   # interpolation date being processed
  mo<-as.integer(strftime(date_proc, "%m"))          # current month of the date being processed
  day<-as.integer(strftime(date_proc, "%d"))
  year<-as.integer(strftime(date_proc, "%Y"))

  datelabel=format(ISOdate(year,mo,day),"%b %d, %Y")
  
  ###########
  #  STEP 1 - LST 10 year monthly averages
  ###########

  themolst<-raster(molst,mo) #current month being processed saved in a raster image
  plot(themolst)
  
  ###########
  # STEP 2 - Weather station means across same days: Monthly mean calculation
  ###########
  
  modst=dst[dst$month==mo,] #Subsetting dataset for the relevant month of the date being processed
  
  ##########
  # STEP 3 - get LST at stations
  ##########
  
  sta_lola=modst[,c("lon","lat")] #Extracting locations of stations for the current month..
  
  proj_str="+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs";
  lookup<-function(r,lat,lon) {
    xy<-project(cbind(lon,lat),proj_str);
    cidx<-cellFromXY(r,xy);
    return(r[cidx])
  }
  sta_tmax_from_lst=lookup(themolst,sta_lola$lat,sta_lola$lon) #Extracted values of LST for the stations
  
  #########
  # STEP 4 - bias at stations     
  #########
  
  sta_bias=sta_tmax_from_lst-modst$TMax; #That is the difference between the monthly LST mean and monthly station mean
  #Added by Benoit
  modst$LSTD_bias<-sta_bias  #Adding bias to data frame modst containning the monthly average for 10 years
  
  bias_xy=project(as.matrix(sta_lola),proj_str)
  png(paste("LST_TMax_scatterplot_",dates[i],out_prefix,".png", sep=""))
  plot(modst$TMax,sta_tmax_from_lst,xlab="Station mo Tmax",ylab="LST mo Tmax",main=paste("LST vs TMax for",datelabel,sep=" "))
  abline(0,1)
  dev.off()
  
  #added by Benoit 
  #x<-ghcn.subsets[[i]]  #Holds both training and testing for instance 161 rows for Jan 1
  x<-data_v
  d<-data_s
  
  pos<-match("value",names(d)) #Find column with name "value"
  #names(d)[pos]<-c("dailyTmax")
  names(d)[pos]<-y_var_name
  names(x)[pos]<-y_var_name
  #names(x)[pos]<-c("dailyTmax")
  d$dailyTmax=(as.numeric(d$dailyTmax))/10 #stored as 1/10 degree C to allow integer storage
  x$dailyTmax=(as.numeric(x$dailyTmax))/10 #stored as 1/10 degree C to allow integer storage
  pos<-match("station",names(d)) #Find column with name "value"
  names(d)[pos]<-c("id")
  names(x)[pos]<-c("id")
  names(modst)[1]<-c("id")       #modst contains the average tmax per month for every stations...
  dmoday=merge(modst,d,by="id")  #LOOSING DATA HERE!!! from 113 t0 103
  xmoday=merge(modst,x,by="id")  #LOOSING DATA HERE!!! from 48 t0 43
  names(dmoday)[4]<-c("lat")
  names(dmoday)[5]<-c("lon")     #dmoday contains all the the information: BIAS, monn
  names(xmoday)[4]<-c("lat")
  names(xmoday)[5]<-c("lon")     #dmoday contains all the the information: BIAS, monn
  
  data_v<-xmoday
  ###
  
  #dmoday contains the daily tmax values for training with TMax being the monthly station tmax mean
  #xmoday contains the daily tmax values for validation with TMax being the monthly station tmax mean
  
  # windows()
  #png(paste("LST_TMax_scatterplot_",dates[i],out_prefix,".png", sep=""))
  png(paste("Daily_tmax_monthly_TMax_scatterplot_",dates[i],out_prefix,".png", sep=""))
  plot(dailyTmax~TMax,data=dmoday,xlab="Mo Tmax",ylab=paste("Daily for",datelabel),main="across stations in OR")
  #savePlot(paste("Daily_tmax_monthly_TMax_scatterplot_",dates[i],out_prefix,".png", sep=""), type="png")
  #png(paste("LST_TMax_scatterplot_",dates[i],out_prefix,".png", sep=""))
  dev.off()
  
  ########
  # STEP 5 - interpolate bias
  ########
  
  # ?? include covariates like elev, distance to coast, cloud frequency, tree height
  #library(fields)
  #windows()
  quilt.plot(sta_lola,sta_bias,main="Bias at stations",asp=1)
  US(add=T,col="magenta",lwd=2)
  #fitbias<-Tps(bias_xy,sta_bias) #use TPS or krige
  
  #Adding options to use only training stations: 07/11/2012
  bias_xy=project(as.matrix(sta_lola),proj_str)
  #bias_xy2=project(as.matrix(c(dmoday$lon,dmoday$lat),proj_str)
  if(bias_val==1){
    sta_bias<-dmoday$LSTD_bias
    bias_xy<-cbind(dmoday$x_OR83M,dmoday$y_OR83M)
  }
  
  fitbias<-Krig(bias_xy,sta_bias,theta=1e5) #use TPS or krige 
  #The output is a krig object using fields
  mod9a<-fitbias
  # Creating plot of bias surface and saving it
  #X11()
  png(paste("Bias_surface_LST_TMax_",dates[i],out_prefix,".png", sep="")) #Create file to write a plot
  datelabel2=format(ISOdate(year,mo,day),"%B ") #added by Benoit, label
  surface(fitbias,col=rev(terrain.colors(100)),asp=1,main=paste("Interpolated bias for",datelabel2,sep=" ")) #Plot to file
  #savePlot(paste("Bias_surface_LST_TMax_",dates[i],out_prefix,".png", sep=""), type="png")
  dev.off()  #Release the hold to the file
  
  #US(add=T,col="magenta",lwd=2)
  
  ##########
  # STEP 7 - interpolate delta across space
  ##########
  
  daily_sta_lola=dmoday[,c("lon","lat")] #could be same as before but why assume merge does this - assume not
  daily_sta_xy=project(as.matrix(daily_sta_lola),proj_str)
  daily_delta=dmoday$dailyTmax-dmoday$TMax
  #windows()
  quilt.plot(daily_sta_lola,daily_delta,asp=1,main="Station delta for Jan 15")
  US(add=T,col="magenta",lwd=2)
  #fitdelta<-Tps(daily_sta_xy,daily_delta) #use TPS or krige
  fitdelta<-Krig(daily_sta_xy,daily_delta,theta=1e5) #use TPS or krige
  #Kriging using fields package
  mod9b<-fitdelta
  # Creating plot of bias surface and saving it
  #X11()
  png(paste("Delta_surface_LST_TMax_",dates[i],out_prefix,".png", sep=""))
  surface(fitdelta,col=rev(terrain.colors(100)),asp=1,main=paste("Interpolated delta for",datelabel,sep=" "))
  #savePlot(paste("Delta_surface_LST_TMax_",dates[i],out_prefix,".png", sep=""), type="png")
  dev.off()
  #US(add=T,col="magenta",lwd=2)
  #
  
  #### Added by Benoit on 06/19
  data_s<-dmoday #put the 
  data_s$daily_delta<-daily_delta
  
  #data_s$y_var<-daily_delta  #y_var is the variable currently being modeled, may be better with BIAS!!
  #data_s$y_var<-data_s$LSTD_bias
  #### Added by Benoit ends
  
  #########
  # STEP 8 - assemble final answer - T=LST+Bias(interpolated)+delta(interpolated)
  #########

  bias_rast=interpolate(themolst,fitbias) #interpolation using function from raster package
  #themolst is raster layer, fitbias is "Krig" object from bias surface
  #plot(bias_rast,main="Raster bias") #This not displaying...
  
  #Saving kriged surface in raster images
  data_name<-paste("bias_LST_",dates[[i]],sep="")
  raster_name<-paste("fusion_",data_name,out_prefix,".rst", sep="")
  writeRaster(bias_rast, filename=raster_name,overwrite=TRUE)  #Writing the data in a raster file format...(IDRISI)
  
  daily_delta_rast=interpolate(themolst,fitdelta) #Interpolation of the bias surface...
  
  #plot(daily_delta_rast,main="Raster Daily Delta")
  
  #Saving kriged surface in raster images
  data_name<-paste("daily_delta_",dates[[i]],sep="")
  raster_name<-paste("fusion_",data_name,out_prefix,".rst", sep="")
  writeRaster(daily_delta_rast, filename=raster_name,overwrite=TRUE)  #Writing the data in a raster file format...(IDRISI)
  
  tmax_predicted=themolst+daily_delta_rast-bias_rast #Final surface  as a raster layer...
  #tmax_predicted=themolst+daily_delta_rast+bias_rast #Added by Benoit, why is it -bias_rast
  #plot(tmax_predicted,main="Predicted daily")
  
  #Saving kriged surface in raster images
  data_name<-paste("tmax_predicted_",dates[[i]],sep="")
  raster_name<-paste("fusion_",data_name,out_prefix,".rst", sep="")
  writeRaster(tmax_predicted, filename=raster_name,overwrite=TRUE)  #Writing the data in a raster file format...(IDRISI)
  
  ########
  # check: assessment of results: validation
  ########
  RMSE<-function(x,y) {return(mean((x-y)^2)^0.5)}
  MAE_fun<-function(x,y) {return(mean(abs(x-y)))}
  #ME_fun<-function(x,y){return(mean(abs(y)))}
  #FIT ASSESSMENT
  sta_pred_data_s=lookup(tmax_predicted,data_s$lat,data_s$lon)
  rmse_fit=RMSE(sta_pred_data_s,data_s$dailyTmax)
  mae_fit=MAE_fun(sta_pred_data_s,data_s$dailyTmax)
    
  sta_pred=lookup(tmax_predicted,data_v$lat,data_v$lon)
  #sta_pred=lookup(tmax_predicted,daily_sta_lola$lat,daily_sta_lola$lon)
  #rmse=RMSE(sta_pred,dmoday$dailyTmax)
  #pos<-match("value",names(data_v)) #Find column with name "value"
  #names(data_v)[pos]<-c("dailyTmax")
  tmax<-data_v$dailyTmax
  #data_v$dailyTmax<-tmax
  rmse=RMSE(sta_pred,tmax)
  mae<-MAE_fun(sta_pred,tmax)
  r2<-cor(sta_pred,tmax)^2              #R2, coef. of var
  me<-mean(sta_pred-tmax)
 
  #plot(sta_pred~dmoday$dailyTmax,xlab=paste("Actual daily for",datelabel),ylab="Pred daily",main=paste("RMSE=",rmse))
  
  png(paste("Predicted_tmax_versus_observed_scatterplot_",dates[i],out_prefix,".png", sep=""))
  plot(sta_pred~tmax,xlab=paste("Actual daily for",datelabel),ylab="Pred daily",main=paste("RMSE=",rmse))
  abline(0,1)
  #savePlot(paste("Predicted_tmax_versus_observed_scatterplot_",dates[i],out_prefix,".png", sep=""), type="png")
  dev.off()
  #resid=sta_pred-dmoday$dailyTmax
  resid=sta_pred-tmax
  quilt.plot(daily_sta_lola,resid)
  
  ### END OF BRIAN's code
  
  ### Added by benoit
  
  ###BEFORE GAM prediction the data object must be transformed to SDF
  
  coords<- data_v[,c('x_OR83M','y_OR83M')]
  coordinates(data_v)<-coords
  proj4string(data_v)<-CRS  #Need to assign coordinates...
  coords<- data_s[,c('x_OR83M','y_OR83M')]
  coordinates(data_s)<-coords
  proj4string(data_s)<-CRS  #Need to assign coordinates..
  
  ns<-nrow(data_s) #This is added to because some loss of data might have happened because of the averaging...
  nv<-nrow(data_v)
  
  ###GAM PREDICTION
  
  #data_s$y_var<-data_s$dailyTmax  #This shoudl be changed for any variable!!!
  #data_v$y_var<-data_v$dailyTmax
  data_v$y_var<-data_v[[y_var_name]]
  data_s$y_var<-data_s[[y_var_name]]
  
  #Model and response variable can be changed without affecting the script
  
  formula1 <- as.formula("y_var ~ s(lat) + s(lon) + s(ELEV_SRTM)", env=.GlobalEnv)
  formula2 <- as.formula("y_var~ s(lat,lon)+ s(ELEV_SRTM)", env=.GlobalEnv)
  formula3 <- as.formula("y_var~ s(lat) + s (lon) + s (ELEV_SRTM) +  s (Northness)+ s (Eastness) + s(DISTOC)", env=.GlobalEnv)
  formula4 <- as.formula("y_var~ s(lat) + s (lon) + s(ELEV_SRTM) + s(Northness) + s (Eastness) + s(DISTOC) + s(LST)", env=.GlobalEnv)
  formula5 <- as.formula("y_var~ s(lat,lon) +s(ELEV_SRTM) + s(Northness,Eastness) + s(DISTOC) + s(LST)", env=.GlobalEnv)
  formula6 <- as.formula("y_var~ s(lat,lon) +s(ELEV_SRTM) + s(Northness,Eastness) + s(DISTOC) + s(LST)+s(LC1)", env=.GlobalEnv)
  formula7 <- as.formula("y_var~ s(lat,lon) +s(ELEV_SRTM) + s(Northness,Eastness) + s(DISTOC) + s(LST)+s(LC3)", env=.GlobalEnv)
  formula8 <- as.formula("y_var~ s(lat,lon) +s(ELEV_SRTM) + s(Northness,Eastness) + s(DISTOC) + s(LST) + s(LC1,LC3)", env=.GlobalEnv)
  
  mod1<- try(gam(formula1, data=data_s))
  mod2<- try(gam(formula2, data=data_s)) #modified nesting....from 3 to 2
  mod3<- try(gam(formula3, data=data_s))
  mod4<- try(gam(formula4, data=data_s))
  mod5<- try(gam(formula5, data=data_s))
  mod6<- try(gam(formula6, data=data_s))
  mod7<- try(gam(formula7, data=data_s))
  mod8<- try(gam(formula8, data=data_s))
  
#   mod1<- try(gam(formula1, data=data_s))
#   mod2<- try(gam(formula2, data=data_s)) #modified nesting....from 3 to 2
#   mod3<- try(gam(y_var~ s(lat) + s (lon) + s (ELEV_SRTM) +  s (Northness)+ s (Eastness) + s(DISTOC), data=data_s))
#   mod4<- try(gam(y_var~ s(lat) + s (lon) + s(ELEV_SRTM) + s(Northness) + s (Eastness) + s(DISTOC) + s(LST), data=data_s))
#   mod5<- try(gam(y_var~ s(lat,lon) +s(ELEV_SRTM) + s(Northness,Eastness) + s(DISTOC) + s(LST), data=data_s))
#   mod6<- try(gam(y_var~ s(lat,lon) +s(ELEV_SRTM) + s(Northness,Eastness) + s(DISTOC) + s(LST)+s(LC1), data=data_s))
#   mod7<- try(gam(y_var~ s(lat,lon) +s(ELEV_SRTM) + s(Northness,Eastness) + s(DISTOC) + s(LST)+s(LC3), data=data_s))
#   mod8<- try(gam(y_var~ s(lat,lon) +s(ELEV_SRTM) + s(Northness,Eastness) + s(DISTOC) + s(LST) + s(LC1,LC3), data=data_s))
#   
  #Added
  #tmax_predicted=themolst+daily_delta_rast-bias_rast #Final surface?? but daily_rst
  
  ### Added by benoit
  #Store results using TPS
  j=nmodels+1
  results_RMSE[1]<- dates[i]    #storing the interpolation dates in the first column
  results_RMSE[2]<- ns          #number of stations used in the training stage
  results_RMSE[3]<- "RMSE"

  results_RMSE[j+3]<- rmse  #Storing RMSE for the model j
  
  results_RMSE_f[1]<- dates[i]    #storing the interpolation dates in the first column
  results_RMSE_f[2]<- ns          #number of stations used in the training stage
  results_RMSE_f[3]<- "RMSE_f"
  results_RMSE_f[j+3]<- rmse_fit  #Storing RMSE for the model j
  
  results_MAE_f[1]<- dates[i]    #storing the interpolation dates in the first column
  results_MAE_f[2]<- ns          #number of stations used in the training stage
  results_MAE_f[3]<- "RMSE_f"
  results_MAE_f[j+3]<- mae_fit  #Storing RMSE for the model j

  results_MAE[1]<- dates[i]    #storing the interpolation dates in the first column
  results_MAE[2]<- ns          #number of stations used in the training stage
  results_MAE[3]<- "MAE"
  results_MAE[j+3]<- mae  #Storing RMSE for the model j

  results_ME[1]<- dates[i]    #storing the interpolation dates in the first column
  results_ME[2]<- ns          #number of stations used in the training stage
  results_ME[3]<- "ME"
  results_ME[j+3]<- me  #Storing RMSE for the model j
  
  results_R2[1]<- dates[i]    #storing the interpolation dates in the first column
  results_R2[2]<- ns          #number of stations used in the training stage
  results_R2[3]<- "R2"
  results_R2[j+3]<- r2  #Storing RMSE for the model j
  
  #ns<-nrow(data_s) #This is added to because some loss of data might have happened because of the averaging...
  #nv<-nrow(data_v)
  
  
  for (j in 1:nmodels){
    
    ##Model assessment: specific diagnostic/metrics for GAM
    
    name<-paste("mod",j,sep="")  #modj is the name of The "j" model (mod1 if j=1) 
    mod<-get(name)               #accessing GAM model ojbect "j"
    
    #If mod "j" is not a model object
    if (inherits(mod,"try-error")) {
      results_AIC[1]<- dates[i]  #storing the interpolation dates in the first column
      results_AIC[2]<- ns        #number of stations used in the training stage
      results_AIC[3]<- "AIC"
      results_AIC[j+3]<- NA
      
      results_GCV[1]<- dates[i]  #storing the interpolation dates in the first column
      results_GCV[2]<- ns        #number of stations used in the training 
      results_GCV[3]<- "GCV"
      results_GCV[j+3]<- NA
      
      results_DEV[1]<- dates[i]  #storing the interpolation dates in the first column
      results_DEV[2]<- ns        #number of stations used in the training stage
      results_DEV[3]<- "DEV"
      results_DEV[j+3]<- NA
      
      results_RMSE_f[1]<- dates[i]  #storing the interpolation dates in the first column
      results_RMSE_f[2]<- ns        #number of stations used in the training stage
      results_RMSE_f[3]<- "RSME_f"
      results_RMSE_f[j+3]<- NA
      
      results_MAE_f[1]<- dates[i]  #storing the interpolation dates in the first column
      results_MAE_f[2]<- ns        #number of stations used in the training stage
      results_MAE_f[3]<- "MAE_f"
      results_MAE_f[j+3]<-NA
      
      results_RMSE[1]<- dates[i]    #storing the interpolation dates in the first column
      results_RMSE[2]<- ns          #number of stations used in the training stage
      results_RMSE[3]<- "RMSE"
      results_RMSE[j+3]<- NA  #Storing RMSE for the model j
      results_MAE[1]<- dates[i]     #storing the interpolation dates in the first column
      results_MAE[2]<- ns           #number of stations used in the training stage
      results_MAE[3]<- "MAE"
      results_MAE[j+3]<- NA    #Storing MAE for the model j
      results_ME[1]<- dates[i]      #storing the interpolation dates in the first column
      results_ME[2]<- ns            #number of stations used in the training stage
      results_ME[3]<- "ME"
      results_ME[j+3]<- NA      #Storing ME for the model j
      results_R2[1]<- dates[i]      #storing the interpolation dates in the first column
      results_R2[2]<- ns            #number of stations used in the training stage
      results_R2[3]<- "R2"
      results_R2[j+3]<- NA      #Storing R2 for the model j
      
    }
    
    #If mod is a modelobject
    
    #If mod "j" is not a model object
    if (inherits(mod,"gam")) {
      
      results_AIC[1]<- dates[i]  #storing the interpolation dates in the first column
      results_AIC[2]<- ns        #number of stations used in the training stage
      results_AIC[3]<- "AIC"
      results_AIC[j+3]<- AIC (mod)
      
      results_GCV[1]<- dates[i]  #storing the interpolation dates in the first column
      results_GCV[2]<- ns        #number of stations used in the training 
      results_GCV[3]<- "GCV"
      results_GCV[j+3]<- mod$gcv.ubre
      
      results_DEV[1]<- dates[i]  #storing the interpolation dates in the first column
      results_DEV[2]<- ns        #number of stations used in the training stage
      results_DEV[3]<- "DEV"
      results_DEV[j+3]<- mod$deviance
      
      y_var_fit= mod$fit
    
      results_RMSE_f[1]<- dates[i]  #storing the interpolation dates in the first column
      results_RMSE_f[2]<- ns        #number of stations used in the training stage
      results_RMSE_f[3]<- "RSME_f"
      #results_RMSE_f[j+3]<- sqrt(sum((y_var_fit-data_s$y_var)^2)/ns)
      results_RMSE_f[j+3]<-sqrt(mean(mod$residuals^2,na.rm=TRUE))
      
      results_MAE_f[1]<- dates[i]  #storing the interpolation dates in the first column
      results_MAE_f[2]<- ns        #number of stations used in the training stage
      results_MAE_f[3]<- "MAE_f"
      #results_MAE_f[j+3]<-sum(abs(y_var_fit-data_s$y_var))/ns
      results_MAE_f[j+3]<-mean(abs(mod$residuals),na.rm=TRUE)
      
      ##Model assessment: general diagnostic/metrics
      ##validation: using the testing data
      if (predval==1) {
      
        ##Model assessment: specific diagnostic/metrics for GAM
        
        name<-paste("mod",j,sep="")  #modj is the name of The "j" model (mod1 if j=1) 
        mod<-get(name)               #accessing GAM model ojbect "j"
        
        s_sgdf<-as(s_raster,"SpatialGridDataFrame") #Conversion to spatial grid data frame
        
        rpred<- predict(mod, newdata=s_sgdf, se.fit = TRUE) #Using the coeff to predict new values.
        y_pred<-rpred$fit
        raster_pred<-r1
        layerNames(raster_pred)<-"y_pred"
        values(raster_pred)<-as.numeric(y_pred)
        data_name<-paste("predicted_mod",j,"_",dates[[i]],sep="")
        raster_name<-paste("GAM_",data_name,out_prefix,".rst", sep="")
        writeRaster(raster_pred, filename=raster_name,overwrite=TRUE)  #Writing the data in a raster file format...(IDRISI)
        #writeRaster(r2, filename=raster_name,overwrite=TRUE)  #Writing the data in a raster file format...(IDRISI)
        
        pred_sgdf<-as(raster_pred,"SpatialGridDataFrame") #Conversion to spatial grid data frame
        #rpred_val_s <- overlay(raster_pred,data_s)             #This overlays the kriged surface tmax and the location of weather stations
        
        rpred_val_s <- overlay(pred_sgdf,data_s)             #This overlays the kriged surface tmax and the location of weather stations
        rpred_val_v <- overlay(pred_sgdf,data_v)             #This overlays the kriged surface tmax and the location of weather stations
        
        pred_mod<-paste("pred_mod",j,sep="")
        #Adding the results back into the original dataframes.
        data_s[[pred_mod]]<-rpred_val_s$y_pred
        data_v[[pred_mod]]<-rpred_val_v$y_pred  
        
        #Model assessment: RMSE and then krig the residuals....!
        
        res_mod_s<- data_s$y_var - data_s[[pred_mod]]           #Residuals from kriging training
        res_mod_v<- data_v$y_var - data_v[[pred_mod]]           #Residuals from kriging validation
        
      }
      
      if (predval==0) {
      
        y_mod<- predict(mod, newdata=data_v, se.fit = TRUE) #Using the coeff to predict new values.
        
        pred_mod<-paste("pred_mod",j,sep="")
        #Adding the results back into the original dataframes.
        data_s[[pred_mod]]<-as.numeric(mod$fit)
        data_v[[pred_mod]]<-as.numeric(y_mod$fit)
        
        #Model assessment: RMSE and then krig the residuals....!
        
        res_mod_s<- data_s$y_var - data_s[[pred_mod]]           #Residuals from kriging training
        res_mod_v<- data_v$y_var - data_v[[pred_mod]]           #Residuals from kriging validation
      }

      ####ADDED ON JULY 20th
      res_mod<-res_mod_v
      
      #RMSE_mod <- sqrt(sum(res_mod^2)/nv)                 #RMSE FOR REGRESSION STEP 1: GAM  
      RMSE_mod<- sqrt(mean(res_mod^2,na.rm=TRUE))
      #MAE_mod<- sum(abs(res_mod),na.rm=TRUE)/(nv-sum(is.na(res_mod)))        #MAE from kriged surface validation
      MAE_mod<- mean(abs(res_mod), na.rm=TRUE)
      #ME_mod<- sum(res_mod,na.rm=TRUE)/(nv-sum(is.na(res_mod)))                    #ME, Mean Error or bias FOR REGRESSION STEP 1: GAM
      ME_mod<- mean(res_mod,na.rm=TRUE)                            #ME, Mean Error or bias FOR REGRESSION STEP 1: GAM
      #R2_mod<- cor(data_v$y_var,data_v[[pred_mod]])^2              #R2, coef. of var FOR REGRESSION STEP 1: GAM
      R2_mod<- cor(data_v$y_var,data_v[[pred_mod]], use="complete")^2
      results_RMSE[1]<- dates[i]    #storing the interpolation dates in the first column
      results_RMSE[2]<- ns          #number of stations used in the training stage
      results_RMSE[3]<- "RMSE"
      results_RMSE[j+3]<- RMSE_mod  #Storing RMSE for the model j
      results_MAE[1]<- dates[i]     #storing the interpolation dates in the first column
      results_MAE[2]<- ns           #number of stations used in the training stage
      results_MAE[3]<- "MAE"
      results_MAE[j+3]<- MAE_mod    #Storing MAE for the model j
      results_ME[1]<- dates[i]      #storing the interpolation dates in the first column
      results_ME[2]<- ns            #number of stations used in the training stage
      results_ME[3]<- "ME"
      results_ME[j+3]<- ME_mod      #Storing ME for the model j
      results_R2[1]<- dates[i]      #storing the interpolation dates in the first column
      results_R2[2]<- ns            #number of stations used in the training stage
      results_R2[3]<- "R2"
      results_R2[j+3]<- R2_mod      #Storing R2 for the model j
      
      #Saving residuals and prediction in the dataframes: tmax predicted from GAM

      name2<-paste("res_mod",j,sep="")
      data_v[[name2]]<-as.numeric(res_mod_v)
      data_s[[name2]]<-as.numeric(res_mod_s)
      #end of loop calculating RMSE
    }
  }
  
  #if (i==length(dates)){
  
  #Specific diagnostic measures related to the testing datasets

  results_table_RMSE<-as.data.frame(results_RMSE)
  results_table_MAE<-as.data.frame(results_MAE)
  results_table_ME<-as.data.frame(results_ME)
  results_table_R2<-as.data.frame(results_R2)
  results_table_RMSE_f<-as.data.frame(results_RMSE_f)
  results_table_MAE_f<-as.data.frame(results_MAE_f)
  
  results_table_AIC<-as.data.frame(results_AIC)
  results_table_GCV<-as.data.frame(results_GCV)
  results_table_DEV<-as.data.frame(results_DEV)
  
  tb_metrics1<-rbind(results_table_RMSE,results_table_MAE, results_table_ME, results_table_R2,results_table_RMSE_f,results_table_MAE_f)   #
  tb_metrics2<-rbind(results_table_AIC,results_table_GCV, results_table_DEV)
  cname<-c("dates","ns","metric","mod1", "mod2","mod3", "mod4", "mod5", "mod6", "mod7","mod8","mod9")
  colnames(tb_metrics1)<-cname
  cname<-c("dates","ns","metric","mod1", "mod2","mod3", "mod4", "mod5", "mod6", "mod7","mod8")
  colnames(tb_metrics2)<-cname
  #colnames(results_table_RMSE)<-cname
  #colnames(results_table_RMSE_f)<-cname
  #tb_diagnostic1<-results_table_RMSE      #measures of validation
  #tb_diagnostic2<-results_table_RMSE_f    #measures of fit
  
  #write.table(tb_diagnostic1, file= paste(path,"/","results_fusion_Assessment_measure1",out_prefix,".txt",sep=""), sep=",")
  
  #}  
  print(paste(dates[i],"processed"))
  # end of the for loop1
  mod_obj<-list(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9a,mod9b)
  names(mod_obj)<-c("mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8","mod9a","mod9b")
  #results_list<-list(data_s,data_v,tb_metrics1,tb_metrics2)
  results_list<-list(data_s,data_v,tb_metrics1,tb_metrics2,mod_obj)
  names(results_list)<-c("data_s","data_v","tb_metrics1","tb_metrics2","mod_obj")
  save(results_list,file= paste(path,"/","results_list_metrics_objects_",dates[i],out_prefix,".RData",sep=""))
  return(results_list)
  #return(tb_diagnostic1)
}