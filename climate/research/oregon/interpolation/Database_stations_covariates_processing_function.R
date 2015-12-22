##################    Data preparation for interpolation  #######################################
############################ Extraction of station data##########################################


database_covariates_preparation<-function(list_param_prep){
  #This function performs queries on the Postgres ghcnd database for stations matching the             
  #interpolation area. It requires 11 inputs:                                           
  # 1)  db.name :  Postgres database name containing the meteorological stations
  # 2)  var: the variable of interest - "TMAX","TMIN" or "PRCP" 
  # 3)  range_years: range of records used in the daily interpolation, note that upper bound year is not included               
  # 4)  range_years_clim: range of records used in the monthly climatology interpolation, note that upper bound is not included
  # 5)  infile_reg_outline: region outline as a shape file - used in the interpolation  stage too                              
  # 6)  infile_ghncd_data: ghcnd stations locations as a textfile name with lat-long fields                                                                                   
  # 7)  infile_covarariates: tif file of raser covariates for the interpolation area: it should have a local projection                                                                                           
  # 8)  CRS_locs_WGS84: longlat EPSG 4326 used as coordinates reference system (proj4)for stations locations
  # 9)  in_path: input path for covariates data and other files, this is also the output?
  # 10) out_path: output path created in master script
  # 11) covar_names: names of covariates used for the interpolation --may be removed later? (should be stored in the brick)
  # 12) qc_flags_stations: flag used to screen out at this stage only two values...
  # 13) out_prefix: output suffix added to output names--it is the same in the interpolation script
  #
  # 
  # 14)
  # 15)..
  #The output is a list of four shapefile names produced by the function:
  #1) loc_stations: outfile1, locations of stations as shapefile in EPSG 4326 
  #2) loc_stations_ghcn: outfile2, ghcn daily data for the year range of interpolation (locally projected) 
  #3) daily_query_ghcn_data: outfile3, ghcn daily data from daily query before application of quality flag
  #4) daily_covar_ghcn_data: outfile4, ghcn daily data with covariates for the year range of interpolation (locally projected)
  #5) monthly_query_ghcn_data: outfile5, ghcn daily data from monthly query before application of quality flag
  #6) monthly_covar_ghcn_data: outfile6, ghcn monthly averaged data with covariates for the year range of interpolation (locally projected)
  
  #AUTHOR: Benoit Parmentier                                                                       
  #CREATED ON: 02/22/2013
  #MODIFIED ON: 12/21/2015
  #PROJECT: NCEAS INPLANT: Environment and Organisms --TASK#363--     
  #Comments:
  #Modifications to take into account the existence of clim query, to avoid querying many times the database at the climatology step.
  #TODO
  #-Add buffer option: solved by overlapping tiles
  #-Add screening for value predicted: var
  # 
  ##################################################################################################
  
  ###Loading R library and packages: should it be read in before???   
  
  library(RPostgreSQL)
  library(sp)                                           # Spatial pacakge with class definition by Bivand et al.
  library(spdep)                                          # Spatial pacakge with methods and spatial stat. by Bivand et al.
  library(rgdal)                                          # GDAL wrapper for R, spatial utilities
  library(rgeos)
  library(rgdal)
  library(raster)
  library(rasterVis)
  library(maps)
  library(maptools)
  
  ### Functions used in the script
  
  format_s <-function(s_ID){
    #Format station ID in a vector format/tuple that is used in a psql query.
    # Argument 1: vector of station ID
    # Return: character of station ID
    tx2<-s_ID
    tx2<-as.character(tx2)
    stat_list<-tx2
    temp<-shQuote(stat_list)
    t<-paste(temp, collapse= " ")
    t1<-gsub(" ", ",",t)
    sf_ID<-paste("(",t1,")",sep="") #vector containing the station ID to query
    return(sf_ID)
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
  
  #parsing input arguments
  
  db.name <- list_param_prep$db.name             #name of the Postgres database
  var <- list_param_prep$var                     #name of the variables to keep: TMIN, TMAX or PRCP
  year_start <-list_param_prep$range_years[1] #"2010"               #starting year for the query (included)
  year_end <-list_param_prep$range_years[2] #"2011"                 #end year for the query (excluded)
  year_start_clim <-list_param_prep$range_years_clim[1] #right bound not included in the range!! starting year for monthly query to calculate clime
  year_end_clim <-list_param_prep$range_years_clim[2] #right bound not included in the range!! starting year for monthly query to calculate clime
  infile_reg_outline<- list_param_prep$infile_reg_outline  #This is the shape file of outline of the study area                                                      #It is an input/output of the covariate script
  infile_ghncd_data<-list_param_prep$infile_ghncd_data      #"/home/layers/data/climate/ghcn/v2.92-upd-2012052822/ghcnd-stations.txt"                              #This is the textfile of station locations from GHCND
  infile_covariates<-list_param_prep$infile_covariates #"covariates__venezuela_region__VE_01292013.tif" #this is an output from covariate script
  CRS_locs_WGS84<-list_param_prep$CRS_locs_WGS84 #Station coords WGS84: same as earlier
  in_path <- list_param_prep$in_path #CRS_locs_WGS84"/home/parmentier/Data/IPLANT_project/Venezuela_interpolation/Venezuela_01142013/input_data/"
  out_path <- list_param_prep$out_path #CRS_locs_WGS84"/home/parmentier/Data/IPLANT_project/Venezuela_interpolation/Venezuela_01142013/input_data/"
  covar_names<-list_param_prep$covar_names # names should be written in the tif file!!!
  qc_flags_stations <- list_param_prep$qc_flags_stations #flags allowed for the query from the GHCND??
  out_prefix<-list_param_prep$out_prefix #"_365d_GAM_fus5_all_lstd_03012013"                #User defined output prefix
  
  ## New parameters added for  sub samplineg in areas with important density of meteo stations
  sub_sampling <- list_param_prep$sub_sampling  #if TRUE then monthly stations data are resampled
  sub_sample_rnd <- list_param_prep$sub_sample_rnd #if  TRUE use random sampling  in addition to spatial  sub-sampling
  target_range_nb <- list_param_prep$target_range_nb # number of stations desired as min and max, convergence to  min  for  now
  
  #needs to be added in master script!!	
  sub_sampling_day <- list_param_prep$sub_sampling_day
  target_range_daily_nb <- list_param_prep$target_range_daily_nb #desired number range of daily stations

  dist_range <- list_param_prep$dist_range #distance range  for pruning,  usually (0,5) in km or 0,0.009*5 for  degreee
  step_dist <- list_param_prep$step_dist #stepping distance used in pruning  spatially, use 1km or 0.009 for degree data

  ## working directory is the same for input and output for this function  
  #setwd(in_path) 
  ##### Create additional output dir for clim datasets
  
  out_path_clim <- file.path(out_path,paste("clim_",year_start_clim,"_",year_end_clim,sep=""))
  create_dir_fun(out_dir=out_path_clim,out_suffix=NULL) #create dir if it does not exists
  setwd(out_path) #thisvaries for clim datasets in the code
  
  out_path <- file.path(out_path,year_start) #if predictions for year 1992, then a 1992 dir is created
  create_dir_fun(out_dir=out_path,out_suffix=NULL) #create dir if it does not exists
  setwd(out_path) #thisvaries for clim datasets in the code
  
  ##### STEP 1: Select station in the study area
    
  filename<-sub(".shp","",fixed=TRUE,infile_reg_outline)             #Removing the extension from file.
  interp_area <- readOGR(dsn=dirname(filename),basename(filename))
  CRS_interp<-proj4string(interp_area)         #Storing the coordinate information: geographic coordinates longlat WGS84
  
  #Read in GHCND database station locations
  dat_stat <- read.fwf(infile_ghncd_data, 
                       widths = c(11,9,10,7,3,31,4,4,6),fill=TRUE)

  colnames(dat_stat)<-c("STAT_ID","latitude","longitude","elev","state","name","GSNF","HCNF","WMOID")
  coords<- dat_stat[,c('longitude','latitude')]
  coordinates(dat_stat)<-coords
  proj4string(dat_stat)<-CRS_locs_WGS84 #this is the WGS84 projection
  #proj4string(dat_stat)<-CRS_interp
  interp_area_WGS84 <-spTransform(interp_area,CRS_locs_WGS84)         # Project from WGS84 to new coord. system
  # Spatial query to find relevant stations

  inside <- !is.na(over(dat_stat, as(interp_area_WGS84, "SpatialPolygons")))  #Finding stations contained in the current interpolation area
  stat_reg<-dat_stat[inside,]              #Selecting stations contained in the current interpolation area
  
  ####
  ##TODO: Add buffer option? Not needed right now as tiles overlap
  ####
  
  #### STEP 2: Connecting to the database and query for relevant data 
  
  ##This is the query for daily time step predictions
  drv <- dbDriver("PostgreSQL")
  db <- dbConnect(drv, dbname=db.name)

  time1<-proc.time()    #Start stop watch
  list_s<-format_s(stat_reg$STAT_ID)
  data2<-dbGetQuery(db, paste("SELECT *
                              FROM ghcn
                              WHERE element=",shQuote(var),
                              "AND year>=",year_start,
                              "AND year<",year_end,
                              "AND station IN ",list_s,";",sep=""))  #Selecting station using a SQL query
  time_duration<-proc.time()-time1             #Time for the query may be long given the size of the database
  time_minutes<-time_duration[3]/60
  
  dbDisconnect(db)
 
  data_table<-merge(data2,as.data.frame(stat_reg), by.x = "station", by.y = "STAT_ID")

  #Transform the subset data frame in a spatial data frame and reproject
  data_reg<-data_table                               #Make a copy of the data frame
  coords<- data_reg[c('longitude','latitude')]              #Define coordinates in a data frame: clean up here!!
  coordinates(data_reg)<-coords                      #Assign coordinates to the data frame
  proj4string(data_reg)<-CRS_locs_WGS84                #Assign coordinates reference system in PROJ4 format
  #data_reg<-spTransform(data_reg,CRS(CRS_interp))     #Project from WGS84 to new coord. system
  
  data_d <-data_reg  #data_d: daily data containing the query without screening
  #data_reg <-subset(data_d,mflag=="0" | mflag=="S") #should be input arguments!!
  #Transform the query to be depending on the number of flags
  
  data_reg <-subset(data_d, mflag %in% qc_flags_stations) #screening using flags
  #data_reg2 <-subset(data_d,mflag==qc_flags_stations[1] | mflag==qc_flags_stations[2]) #screening using flags

  ##################################################################
  ### STEP 3: Save results and outuput in textfile and a shape file
  #browser()
  #Save shape files of the locations of meteorological stations in the study area
  #outfile1<-file.path(in_path,paste("stations","_",out_prefix,".shp",sep=""))
  outfile1<-file.path(out_path,paste("stations","_",out_prefix,".shp",sep=""))
  writeOGR(stat_reg,dsn= dirname(outfile1),layer= sub(".shp","",basename(outfile1)), driver="ESRI Shapefile",overwrite_layer=TRUE)
  #writeOGR(dst,dsn= ".",layer= sub(".shp","",outfile4), driver="ESRI Shapefile",overwrite_layer=TRUE)
  #Also save as rds
  outfile1_rds<-sub(".shp",".rds",outfile1)
  saveRDS(stat_reg,outfile1)

  #Also save in clim data
  outfile1_clim <-file.path(out_path_clim,paste("stations","_",out_prefix,".shp",sep=""))
  writeOGR(stat_reg,dsn= dirname(outfile1),layer= sub(".shp","",basename(outfile1)), driver="ESRI Shapefile",overwrite_layer=TRUE)
  #writeOGR(dst,dsn= ".",layer= sub(".shp","",outfile4), driver="ESRI Shapefile",overwrite_layer=TRUE)
  #Also save as rds
  outfile1_rds<-sub(".shp",".rds",outfile1)
  saveRDS(stat_reg,outfile1)
  
  ##This goes in the daily output dir (e.g. for year 1992)
  #/nobackupp6/aguzman4/climateLayers/out/reg4/20_-60/1992/
  
  outfile2<-file.path(out_path,paste("ghcn_data_query_",var,"_",year_start,"_",year_end,out_prefix,".shp",sep=""))         #Name of the file
  #writeOGR(data_proj, paste(outfile, "shp", sep="."), outfile, driver ="ESRI Shapefile") #Note that the layer name is the file name without extension
  writeOGR(data_d,dsn= dirname(outfile2),layer= sub(".shp","",basename(outfile2)), driver="ESRI Shapefile",overwrite_layer=TRUE)
  outfile2_rds<-sub(".shp",".rds",outfile2)
  saveRDS(data_d,outfile2_rds) 

  outfile3<-file.path(out_path,paste("ghcn_data_",var,"_",year_start,"_",year_end,out_prefix,".shp",sep=""))         #Name of the file
  #writeOGR(data_proj, paste(outfile, "shp", sep="."), outfile, driver ="ESRI Shapefile") #Note that the layer name is the file name without extension
  writeOGR(data_reg,dsn= dirname(outfile3),layer= sub(".shp","",basename(outfile3)), driver="ESRI Shapefile",overwrite_layer=TRUE)
  outfile3_rds<-sub(".shp",".rds",outfile3)
  saveRDS(data_reg,outfile3_rds)

  ###################################################################
  ### STEP 4: Extract values at stations from covariates stack of raster images
  #Eventually this step may be skipped if the covariates information is stored in the database...
  #s_raster<-stack(file.path(in_path,infile_covariates))                   #read in the data stack
  
  s_raster<-brick(infile_covariates)                   #read in the data stack
  names(s_raster)<-covar_names               #Assigning names to the raster layers: making sure it is included in the extraction
  stat_val<- extract(s_raster, data_reg)        #Extracting values from the raster stack for every point location in coords data frame.
  #stat_val_test<- extract(s_raster, data_reg,def=TRUE)
  
  #create a shape file and data_frame with names
  data_RST<-as.data.frame(stat_val)                                            #This creates a data frame with the values extracted
  data_RST_SDF<-cbind(data_reg,data_RST)

  coordinates(data_RST_SDF)<-coordinates(data_reg) #Transforming data_RST_SDF into a spatial point dataframe
  CRS_reg<-proj4string(data_reg)
  proj4string(data_RST_SDF)<-CRS_reg  #Need to assign coordinates...

  #Creating a date column
  date1<-ISOdate(data_RST_SDF$year,data_RST_SDF$month,data_RST_SDF$day) #Creating a date object from 3 separate column
  date2<-gsub("-","",as.character(as.Date(date1)))
  data_RST_SDF$date<-date2                                              #Date format (year,month,day) is the following: "20100627"
 
  #This allows to change only one name of the data.frame
  pos<-match("value",names(data_RST_SDF)) #Find column with name "value"
  if (var=="TMAX"){
    #names(data_RST_SDF)[pos]<-c("TMax")
    data_RST_SDF$value<-data_RST_SDF$value/10                #TMax is the average max temp for monthy data
  }
  
  if (var=="TMIN"){
    #names(data_RST_SDF)[pos]<-c("TMin")
    data_RST_SDF$value<-data_RST_SDF$value/10                #TMax is the average max temp for monthy data
  }
  
   ## Adding subsampling for daily stations...
  
    #This must be set up in master script
  #target_max_nb <- 100,000 #this is not actually used yet in the current implementation,can be set to very high value...
  #target_min_nb <- 600 #this is the target number of stations we would like for daily and 1000x3000 tiles   
                        #to be set by Alberto...THIS differentent than monthly number!!!
  
  ##max_dist <- 1000 # the maximum distance used for pruning ie removes stations that are closer than 1000m, this in degree...? 
  #max_dist <- 0.009*5 #5km in degree
  #min_dist <- 0    #minimum distance to start with
  #step_dist <- 0.009 #iteration step to remove the stations

  #test5 <- sub_sampling_by_dist_nb_stat(target_range_nb=target_range_nb,dist_range=dist_range,step_dist=step_dist,data_in=data_month,sampling=T,combined=F)
  
  #Daily range set at the begining...line 199!
  
  #target_range_daily_nb <- list_param$target_range_daily_nb #desired number range of daily stations
  
  if(sub_sampling_day==TRUE){
    print("daily subsample")
    sub_sampling_obj <- sub_sampling_by_dist_nb_stat(target_range_nb=target_range_nb,dist_range=dist_range,step_dist=step_dist,data_in=data_RST_SDF,sampling=T,combined=F)
    data_RST_SDF <- sub_sampling_obj$data #get sub-sampled data...for monhtly stations
    #save the information for later use (validation at monthly step!!)
    save(sub_sampling_obj,file= file.path(out_path,paste("sub_sampling_obj_","daily_",interpolation_method,"_", out_prefix,".RData",sep="")))
  }

  #Make sure this is still a shapefile...!! This might need to be uncommented...

  #coordinates(data_RST_SDF)<-cbind(data_RST_SDF$x,data_RST_SDF$y) #Transforming data_RST_SDF into a spatial point dataframe
  #CRS_reg<-proj4string(data_RST_SDF)
  #proj4string(data_RST_SDF)<-CRS_reg  #Need to assign coordinates...

  #write out a new shapefile (including .prj component)
  outfile4<-file.path(out_path,paste("daily_covariates_ghcn_data_",var,"_",range_years[1],"_",range_years[2],out_prefix,".shp",sep=""))         #Name of the file
  writeOGR(data_RST_SDF,dsn= dirname(outfile4),layer= sub(".shp","",basename(outfile4)), driver="ESRI Shapefile",overwrite_layer=TRUE)
  outfile4_rds<-sub(".shp",".rds",outfile4)
  saveRDS(data_RST_SDF,outfile4_rds)  

  ###############################################################
  ######## STEP 5: Preparing monthly averages from the ProstGres database
  
  ## Only query for monthly climatology data if data does not exist in the directory
  
  #This file contains the input data for climatology predictions
  outfile6 <- file.path(out_path_clim,paste("monthly_covariates_ghcn_data_",var,"_",year_start_clim,"_",year_end_clim,out_prefix,".shp",sep=""))  #Name of the file

  #Could also be replaced by run_clim_query parameter.
  if(!(file.exists(outfile6))){
    
    drv <- dbDriver("PostgreSQL")
    db <- dbConnect(drv, dbname=db.name)
  
    #year_start_clim: set at the start of the script
    time1<-proc.time()    #Start stop watch
    list_s<-format_s(stat_reg$STAT_ID)
 
    data_m <- dbGetQuery(db, paste("SELECT *
                               FROM ghcn
                               WHERE element=",shQuote(var),
                               "AND year>=",year_start_clim,
                               "AND station IN ",list_s,";",sep=""))  #Selecting station using a SQL query
    time_duration<-proc.time()-time1             #Time for the query may be long given the size of the database
    time_minutes<-time_duration[3]/60
    print(time_minutes)
    dbDisconnect(db)
  
    #Clean out this section!!
    date1<-ISOdate(data_m$year,data_m$month,data_m$day) #Creating a date object from 3 separate column
    date2<-as.POSIXlt(as.Date(date1))
    data_m$date<-date2
    #Save the query data here...
    data_m<-merge(data_m, stat_reg, by.x="station", by.y="STAT_ID")   #Inner join all columns are retained
    #Extracting covariates from stack for the monthly dataset...
    coords<- data_m[c('longitude','latitude')]              #Define coordinates in a data frame
    coordinates(data_m)<-coords                      #Assign coordinates to the data frame
    proj4string(data_m)<-CRS_locs_WGS84                  #Assign coordinates reference system in PROJ4 format
    data_m<-spTransform(data_m,CRS(CRS_interp))     #Project from WGS84 to new coord. system
  
    outfile5<-file.path(out_path_clim,paste("monthly_query_ghcn_data_",var,"_",year_start_clim,"_",year_end_clim,out_prefix,".shp",sep=""))  #Name of the file
  
    #outfile5<-file.path(out_path,paste("monthly_query_ghcn_data_",var,"_",year_start_clim,"_",year_end_clim,out_prefix,".shp",sep=""))  #Name of the file
    writeOGR(data_m,dsn= dirname(outfile5),layer= sub(".shp","",basename(outfile5)), driver="ESRI Shapefile",overwrite_layer=TRUE)
  
    outfile5_rds<-sub(".shp",".rds",outfile5)
    saveRDS(data_m,outfile5_rds) 

    #In Venezuela and other regions where there are not many stations...mflag==S should be added..see Durenne etal.2010.

    #d<-subset(data_m,mflag==qc_flags_stations[1] | mflag==qc_flags_stations[2])
    d<-subset(data_m,mflag %in% qc_flags_stations) 
  
    #Add screening here ...May need some screeing??? i.e. range of temp and elevation...
  
    d1<-aggregate(value~station+month, data=d, mean)  #Calculate monthly mean for every station in OR
    #d2<-aggregate(value~station+month, data=d, length)  #Calculate monthly mean for every station in OR
    is_not_na_fun<-function(x) sum(!is.na(x)) #count the number of available observation
    d2<-aggregate(value~station+month, data=d, is_not_na_fun)  #Calculate monthly mean for every station in OR
    #names(d2)[names(d2)=="value"] <-"nobs_station"
    d1$nobs_station<-d2$value
    dst<-merge(d1, stat_reg, by.x="station", by.y="STAT_ID")   #Inner join all columns are retained
  
    #This allows to change only one name of the data.frame
    pos<-match("value",names(dst)) #Find column with name "value"
    if (var=="TMAX"){
      names(dst)[pos]<-c("TMax")           #renaming the column "value" extracted from the Postgres database
      dst$TMax<-dst$TMax/10                #TMax is the average max temp for monthy data
    }
  
    if (var=="TMIN"){
      names(dst)[pos]<-c("TMin")
      dst$TMin<-dst$TMin/10                #TMin is the average min temp for monthy data
    }
  
    #Extracting covariates from stack for the monthly dataset...
    #names(dst)[5:6] <-c('latitude','longitude')
    coords<- dst[c('longitude','latitude')]              #Define coordinates in a data frame

    coordinates(dst)<-coords                      #Assign coordinates to the data frame
    proj4string(dst)<-CRS_locs_WGS84                  #Assign coordinates reference system in PROJ4 format
    dst_month<-spTransform(dst,CRS(CRS_interp))     #Project from WGS84 to new coord. system
  
    stations_val<-extract(s_raster,dst_month,df=TRUE)  #extraction of the information at station location in a data frame
    #dst_extract<-spCbind(dst_month,stations_val) #this is in sinusoidal from the raster stack
    dst_extract<-cbind(dst_month,stations_val) #this is in sinusoidal from the raster stack
    dst<-dst_extract #problem!!! two column named elev!!! use elev_s??
  
    #browser()
    coords<- dst[c('x','y')]              #Define coordinates in a data frame, this is the local x,y
    index<-!is.na(coords$x)               #remove if NA, may need to revisit at a later stage
    dst<-dst[index,]
    coords<- dst[c('x','y')]              #Define coordinates in a data frame, this is the local x,y
    coordinates(dst)<-coords                    #Assign coordinates to the data frame
    proj4string(dst)<-CRS_interp        #Assign coordinates reference system in PROJ4 format

    ##Added 01-07-2015
    dst$id <- dst$station #the id field is needed for possible subsampling
  
    ### ADD SCREENING HERE BEFORE WRITING OUT DATA
    #Covariates ok since screening done in covariate script
    #screening on var i.e. value, TMIN, TMAX...
  
    #### Adding  subsampling for regions  that  have  too  many stations...
  
    #This must be set up in master script
    #target_max_nb <- 100,000 #this is not actually used yet in the current implementation,can be set to very high value...
    #target_min_nb <- 8,000 #this is the target number of stations we would like: to be set by Alberto...
    ##max_dist <- 1000 # the maximum distance used for pruning ie removes stations that are closer than 1000m, this in degree...? 
    #max_dist <- 0.009*5 #5km in degree
    #min_dist <- 0    #minimum distance to start with
    #step_dist <- 0.009 #iteration step to remove the stations

    #test5 <- sub_sampling_by_dist_nb_stat(target_range_nb=target_range_nb,dist_range=dist_range,step_dist=step_dist,data_in=data_month,sampling=T,combined=F)
    #note that  this is for monthly stations.

    if(sub_sampling==TRUE){
      print("monthly subsample")
      sub_sampling_obj <- sub_sampling_by_dist_nb_stat(target_range_nb=target_range_nb,dist_range=dist_range,step_dist=step_dist,data_in=dst,sampling=T,combined=F)
      dst <- sub_sampling_obj$data #get sub-sampled data...for monhtly stations
      #save the information for later use (validation at monthly step!!)
      #save(sub_sampling_obj,file= file.path(out_path,paste("sub_sampling_obj_",interpolation_method,"_", out_prefix,".RData",sep="")))
      save(sub_sampling_obj,file= file.path(out_path_clim,paste("sub_sampling_obj_",interpolation_method,"_", out_prefix,".RData",sep="")))
    }
 
    #Make sure this is still a shapefile...!! This might need to be uncommented...
    #coordinates(dst)<-cbind(dst$x,dst$y) #Transforming data_RST_SDF into a spatial point dataframe
    #CRS_reg<-proj4string(data_reg)
    #proj4string(dst)<-CRS_reg  #Need to assign coordinates...

    ####
    #write out a new shapefile (including .prj component)
    dst$OID<-1:nrow(dst) #need a unique ID?
    #outfile6<-file.path(out_path,paste("monthly_covariates_ghcn_data_",var,"_",year_start_clim,"_",year_end_clim,out_prefix,".shp",sep=""))  #Name of the file
    outfile6<-file.path(out_path_clim,paste("monthly_covariates_ghcn_data_",var,"_",year_start_clim,"_",year_end_clim,out_prefix,".shp",sep=""))  #Name of the file

    writeOGR(dst,dsn= dirname(outfile6),layer= sub(".shp","",basename(outfile6)), driver="ESRI Shapefile",overwrite_layer=TRUE)
  
    outfile6_rds<-sub(".shp",".rds",outfile6)
    saveRDS(dst,outfile6_rds)

  }
    
  ###If clim exists, just create name for outfile5 to object used later...
  #outfile6 <- file.path(out_path_clim,paste("monthly_covariates_ghcn_data_",var,"_",year_start_clim,"_",year_end_clim,out_prefix,".shp",sep=""))  #Name of the file
  #Could also be replaced by run_clim_query parameter.
  if(file.exists(outfile6)){
    outfile5<-file.path(out_path_clim,paste("monthly_query_ghcn_data_",var,"_",year_start_clim,"_",year_end_clim,out_prefix,".shp",sep=""))  #Name of the file
  }

  outfiles_obj<-list(outfile1,outfile2,outfile3,outfile4,outfile5,outfile6)
  names(outfiles_obj)<- c("loc_stations","loc_stations_ghcn","daily_query_ghcn_data","daily_covar_ghcn_data","monthly_query_ghcn_data","monthly_covar_ghcn_data")
  save(outfiles_obj,file= file.path(out_path,paste("met_stations_outfiles_obj_",interpolation_method,"_", out_prefix,".RData",sep="")))
  
  return(outfiles_obj)
  
  #END OF FUNCTION # 
}

########## END OF SCRIPT ##########
