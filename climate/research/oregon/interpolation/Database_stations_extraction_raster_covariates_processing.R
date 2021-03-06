##################    Data preparation for interpolation   #######################################
############################ Extraction of station data ##########################################


database_covaratiates_preparation<-function(list_param_prep){
  #This function performs queries on the Postgres ghcnd database for stations matching the             
  #interpolation area. It requires 11 inputs:                                           
  # 1) db.name :  Postgres database name containing the meteorological stations
  # 2) var: the variable of interest - "TMAX","TMIN" or "PRCP" 
  # 3) range_years: range of records used in the daily interpolation, note that upper bound year is not included               
  # 4) range_years_clim: range of records used in the monthly climatology interpolation, note that upper bound is not included
  # 5) infile1: region outline as a shape file - used in the interpolation  stage too                              
  # 6) infile2: ghcnd stations locations as a textfile name with lat-long fields                                                                                   
  # 7) infile_covarariates: tif file of raser covariates for the interpolation area: it should have a local projection                                                                                           
  # 8) CRS_locs_WGS84: longlat EPSG 4326 used as coordinates reference system (proj4)for stations locations
  # 9) in_path: input path for covariates data and other files, this is also the output?
  # 10) covar_names: names of covariates used for the interpolation --may be removed later? (should be stored in the brick)
  # 11) out_prefix: output suffix added to output names--it is the same in the interpolation script
  #
  #The output is a list of four shapefile names produced by the function:
  #1) loc_stations: locations of stations as shapefile in EPSG 4326
  #2) loc_stations_ghcn: ghcn daily data for the year range of interpolation (locally projected)
  #3) daily_covar_ghcn_data: ghcn daily data with covariates for the year range of interpolation (locally projected)
  #4) monthly_covar_ghcn_data: ghcn daily data with covariates for the year range of interpolation (locally projected)
  
  #AUTHOR: Benoit Parmentier                                                                       
  #DATE: 03/01/2013                                                                                 
  #PROJECT: NCEAS INPLANT: Environment and Organisms --TASK#363--     
  #Comments and TODO
  #-Add buffer option...
  #-Add output path argument option
  #-Add qc flag options
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
  
  #parsing input arguments
  
  db.name <- list_param_prep$db.name             #name of the Postgres database
  var <- list_param_prep$var                     #name of the variables to keep: TMIN, TMAX or PRCP
  year_start <-list_param_prep$range_years[1] #"2010"               #starting year for the query (included)
  year_end <-list_param_prep$range_years[2] #"2011"                 #end year for the query (excluded)
  year_start_clim <-list_param_prep$range_years_clim[1] #right bound not included in the range!! starting year for monthly query to calculate clime
  infile1<- list_param_prep$infile1  #This is the shape file of outline of the study area                                                      #It is an input/output of the covariate script
  infile2<-list_param_prep$infile2      #"/home/layers/data/climate/ghcn/v2.92-upd-2012052822/ghcnd-stations.txt"                              #This is the textfile of station locations from GHCND
  infile3<-list_param_prep$infile_covariates #"covariates__venezuela_region__VE_01292013.tif" #this is an output from covariate script
  CRS_locs_WGS84<-list_param_prep$CRS_locs_WGS84 #Station coords WGS84: same as earlier
  in_path <- list_param_prep$in_path #CRS_locs_WGS84"/home/parmentier/Data/IPLANT_project/Venezuela_interpolation/Venezuela_01142013/input_data/"
  out_prefix<-list_param_prep$out_prefix #"_365d_GAM_fus5_all_lstd_03012013"                #User defined output prefix
  #qc_flags<-list_param_prep$qc_flags    flags allowed for the query from the GHCND??
  covar_names<-list_param_prep$covar_names # names should be written in the tif file!!!
  
  ## working directory is the same for input and output for this function  
  setwd(in_path) 
  
  ##### STEP 1: Select station in the study area
  
  filename<-sub(".shp","",infile1)             #Removing the extension from file.
  interp_area <- readOGR(".",filename)
  CRS_interp<-proj4string(interp_area)         #Storing the coordinate information: geographic coordinates longlat WGS84
  
  #Read in GHCND database station locations
  dat_stat <- read.fwf(infile2, 
                       widths = c(11,9,10,7,3,31,4,4,6),fill=TRUE)
  colnames(dat_stat)<-c("STAT_ID","lat","lon","elev","state","name","GSNF","HCNF","WMOID")
  coords<- dat_stat[,c('lon','lat')]
  coordinates(dat_stat)<-coords
  proj4string(dat_stat)<-CRS_locs_WGS84 #this is the WGS84 projection
  #proj4string(dat_stat)<-CRS_interp
  dat_stat2<-spTransform(dat_stat,CRS(CRS_interp))         # Project from WGS84 to new coord. system
  
  # Spatial query to find relevant stations
  inside <- !is.na(over(dat_stat2, as(interp_area, "SpatialPolygons")))  #Finding stations contained in the current interpolation area
  stat_reg<-dat_stat2[inside,]              #Selecting stations contained in the current interpolation area
  
  ####
  ##TODO: Add buffer option? 
  ####
  
  #### STEP 2: Connecting to the database and query for relevant data 
  
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
  ###
  #Add month query and averages here...
  ###
  
  #data2 contains only 46 stations for Venezueal area??
  data_table<-merge(data2,as.data.frame(stat_reg), by.x = "station", by.y = "STAT_ID")
  
  #Transform the subset data frame in a spatial data frame and reproject
  data_reg<-data_table                               #Make a copy of the data frame
  coords<- data_reg[c('lon','lat')]              #Define coordinates in a data frame: clean up here!!
  #Wrong label...it is in fact projected...
  coordinates(data_reg)<-coords                      #Assign coordinates to the data frame
  #proj4string(data3)<-locs_coord                  #Assign coordinates reference system in PROJ4 format
  proj4string(data_reg)<-CRS_locs_WGS84                #Assign coordinates reference system in PROJ4 format
  data_reg<-spTransform(data_reg,CRS(CRS_interp))     #Project from WGS84 to new coord. system
  
  ##################################################################
  ### STEP 3: Save results and outuput in textfile and a shape file
  #browser()
  #Save shape files of the locations of meteorological stations in the study area
  outfile1<-file.path(in_path,paste("stations","_",out_prefix,".shp",sep=""))
  writeOGR(stat_reg,dsn= dirname(outfile1),layer= sub(".shp","",basename(outfile1)), driver="ESRI Shapefile",overwrite_layer=TRUE)
  #writeOGR(dst,dsn= ".",layer= sub(".shp","",outfile4), driver="ESRI Shapefile",overwrite_layer=TRUE)
  
  outfile2<-file.path(in_path,paste("ghcn_data_",var,"_",year_start_clim,"_",year_end,out_prefix,".shp",sep=""))         #Name of the file
  #writeOGR(data_proj, paste(outfile, "shp", sep="."), outfile, driver ="ESRI Shapefile") #Note that the layer name is the file name without extension
  writeOGR(data_reg,dsn= dirname(outfile2),layer= sub(".shp","",basename(outfile2)), driver="ESRI Shapefile",overwrite_layer=TRUE)
  
  ###################################################################
  ### STEP 4: Extract values at stations from covariates stack of raster images
  #Eventually this step may be skipped if the covariates information is stored in the database...
  
  s_raster<-stack(infile3)                   #read in the data stack
  names(s_raster)<-covar_names               #Assigning names to the raster layers: making sure it is included in the extraction
  stat_val<- extract(s_raster, data_reg)        #Extracting values from the raster stack for every point location in coords data frame.
  
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
  
  #write out a new shapefile (including .prj component)
  outfile3<-file.path(in_path,paste("daily_covariates_ghcn_data_",var,"_",range_years[1],"_",range_years[2],out_prefix,".shp",sep=""))         #Name of the file
  writeOGR(data_RST_SDF,dsn= dirname(outfile3),layer= sub(".shp","",basename(outfile3)), driver="ESRI Shapefile",overwrite_layer=TRUE)
  
  ###############################################################
  ######## STEP 5: Preparing monthly averages from the ProstGres database
  
  drv <- dbDriver("PostgreSQL")
  db <- dbConnect(drv, dbname=db.name)
  
  #year_start_clim: set at the start of the script
  year_end<-2011
  time1<-proc.time()    #Start stop watch
  list_s<-format_s(stat_reg$STAT_ID)
  data_m<-dbGetQuery(db, paste("SELECT *
                               FROM ghcn
                               WHERE element=",shQuote(var),
                               "AND year>=",year_start_clim,
                               "AND year<",year_end,
                               "AND station IN ",list_s,";",sep=""))  #Selecting station using a SQL query
  time_duration<-proc.time()-time1             #Time for the query may be long given the size of the database
  time_minutes<-time_duration[3]/60
  dbDisconnect(db)
  #Clean out this section!!
  date1<-ISOdate(data_m$year,data_m$month,data_m$day) #Creating a date object from 3 separate column
  date2<-as.POSIXlt(as.Date(date1))
  data_m$date<-date2
  #In Venezuela and other regions where there are not many stations...mflag==S should be added..see Durenne etal.2010.
  #d<-subset(data_m,year>=2000 & mflag=="0" ) #Selecting dataset 2000-2010 with good quality: 193 stations
  d<-subset(data_m,mflag=="0" | mflag=="S") #should be input arguments!!
  #May need some screeing??? i.e. range of temp and elevation...
  d1<-aggregate(value~station+month, data=d, mean)  #Calculate monthly mean for every station in OR
  id<-as.data.frame(unique(d1$station))     #Unique station in OR for year 2000-2010: 193 but 7 loss of monthly avg    
  
  dst<-merge(d1, stat_reg, by.x="station", by.y="STAT_ID")   #Inner join all columns are retained
  
  #This allows to change only one name of the data.frame
  pos<-match("value",names(dst)) #Find column with name "value"
  if (var=="TMAX"){
    names(dst)[pos]<-c("TMax")
    dst$TMax<-dst$TMax/10                #TMax is the average max temp for monthy data
  }
  
  #Extracting covariates from stack for the monthly dataset...
  coords<- dst[c('lon','lat')]              #Define coordinates in a data frame
  coordinates(dst)<-coords                      #Assign coordinates to the data frame
  proj4string(dst)<-CRS_locs_WGS84                  #Assign coordinates reference system in PROJ4 format
  dst_month<-spTransform(dst,CRS(CRS_interp))     #Project from WGS84 to new coord. system
  
  stations_val<-extract(s_raster,dst_month)  #extraction of the infomration at station location
  stations_val<-as.data.frame(stations_val)
  dst_extract<-cbind(dst_month,stations_val) #this is in sinusoidal from the raster stack
  dst<-dst_extract
  
  coords<- dst[c('x','y')]              #Define coordinates in a data frame, this is the local x,y
  coordinates(dst)<-coords                    #Assign coordinates to the data frame
  proj4string(dst)<-projection(s_raster)        #Assign coordinates reference system in PROJ4 format
  
  ####
  #write out a new shapefile (including .prj component)
  dst$OID<-1:nrow(dst) #need a unique ID?
  outfile4<-file.path(in_path,paste("monthly_covariates_ghcn_data_",var,"_",year_start_clim,"_",year_end,out_prefix,".shp",sep=""))  #Name of the file
  writeOGR(dst,dsn= dirname(outfile4),layer= sub(".shp","",basename(outfile4)), driver="ESRI Shapefile",overwrite_layer=TRUE)
  
  ### list of outputs return
  
  outfiles_obj<-list(outfile1,outfile2,outfile3,outfile4)
  names(outfiles_obj)<- c("loc_stations","loc_stations_ghcn","daily_covar_ghcn_data","monthly_covar_ghcn_data")
  return(outfiles_obj)
  
  #END OF FUNCTION # 
}

##### END OF SCRIPT ##########
