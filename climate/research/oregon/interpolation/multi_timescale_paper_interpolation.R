####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for manuscript analyses,tables and figures: MULTI TIME SCALE  ##############################
#This script uses the worklfow code applied to the Oregon case study. Single and Multi timescale procedures are compared.
#Single and Multitime scale (CAI and FSS)procedures are tested with different covariate combinations
#using GAM,GWR, Kriging. Results conurrently reside on NCEAS ATLAS server.
#Accuracy methods are added in the the function scripts to evaluate results.
#Analyses, figures, tables and data for the  paper are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/31/2013  
#MODIFIED ON: 05/05/2014            
#Version: 5
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
function_analyses_paper2 <-"multi_timescales_paper_interpolation_functions_05052014.R"

##############################
#### Parameters and constants  

script_path<-"/home/parmentier/Data/IPLANT_project/env_layers_scripts/" #path to script
source(file.path(script_path,function_analyses_paper1)) #source all functions used in this script 1.
source(file.path(script_path,function_analyses_paper2)) #source all functions used in this script 2.

#Singlte time scale -direct methods: gam, kriging, gwr (no monthly holdout but 30%daily hold out)
in_dir1 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gam_daily_lst_comb5_11012013"
in_dir2 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_kriging_daily_lst_comb5_11022013"
in_dir3a <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gwr_daily_lst_comb5p1_3_11062013"
in_dir3b <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gwr_daily_lst_comb5p4_7_11292013"

#Multi time scale - CAI: gam, kriging, gwr
in_dir4 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gam_cai_lst_comb5_11032013"
in_dir5 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_kriging_cai_lst_comb5_11032013"
in_dir6 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gwr_cai_lst_comb5_11042013"
#Multi time scale - FSS: gam, kriging, gwr
in_dir7 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gam_fss_lst_comb5_11062013"
in_dir8 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_kriging_fss_lst_comb5_11052013"
in_dir9 <- "/data/project/layers/commons/Oregon_interpolation/output_data_365d_gwr_fss_lst_comb5_11052013"
###

### Multi time scale: Results from monthly hold out 0-70%
in_dir10 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gam_fss_lst_mults_0_70_comb5_11082013"
in_dir11 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_kriging_fss_lst_mults_0_70_comb5_11132013"
in_dir12 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gwr_fss_lst_mults_0_70_comb5_11162013"
in_dir13 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gam_cai_lst_mults_0_70_comb5_11192013"
in_dir14 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_kriging_cai_lst_mults_0_70_comb5_11272013"
in_dir15 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gwr_cai_lst_mults_0_70_comb5_11222013"

##raster_prediction object for comb5
#direct methods
raster_obj_file_1 <- "raster_prediction_obj_gam_daily_dailyTmax_365d_gam_daily_lst_comb5_11012013.RData" 
raster_obj_file_2 <- "raster_prediction_obj_kriging_daily_dailyTmax_365d_kriging_daily_lst_comb5_11022013.RData"
raster_obj_file_3a <- "raster_prediction_obj_gwr_daily_dailyTmax_365d_gwr_daily_lst_comb5p1_3_11062013.RData"
raster_obj_file_3b <- "raster_prediction_obj_gwr_daily_dailyTmax_365d_gwr_daily_lst_comb5p4_7_11292013.RData"
#CAI
raster_obj_file_4 <- "raster_prediction_obj_gam_CAI_dailyTmax_365d_gam_cai_lst_comb5_11032013.RData"
raster_obj_file_5 <- "raster_prediction_obj_kriging_CAI_dailyTmax_365d_kriging_cai_lst_comb5_11032013.RData"
raster_obj_file_6 <- "raster_prediction_obj_gwr_CAI_dailyTmax_365d_gwr_cai_lst_comb5_11042013.RData"
#FSS
raster_obj_file_7 <- "raster_prediction_obj_gam_fusion_dailyTmax_365d_gam_fss_lst_comb5_11062013.RData"
raster_obj_file_8 <- "raster_prediction_obj_kriging_fusion_dailyTmax_365d_kriging_fss_lst_comb5_11052013.RData"
raster_obj_file_9 <- "raster_prediction_obj_gwr_fusion_dailyTmax_365d_gwr_fss_lst_comb5_11052013.RData"

## Results from monthly holdout 0 to 70%

raster_obj_file_10 <- "raster_prediction_obj_gam_fusion_dailyTmax_365d_gam_fss_lst_mults_0_70_comb5_11082013.RData"
raster_obj_file_11 <- "raster_prediction_obj_kriging_fusion_dailyTmax_365d_kriging_fss_lst_mults_0_70_comb5_11132013.RData"
raster_obj_file_12 <- "raster_prediction_obj_gwr_fusion_dailyTmax_365d_gwr_fss_lst_mults_0_70_comb5_11162013.RData"
raster_obj_file_13 <- "raster_prediction_obj_gam_CAI_dailyTmax_365d_gam_cai_lst_mults_0_70_comb5_11192013.RData"
raster_obj_file_14 <- "raster_prediction_obj_kriging_CAI_dailyTmax_365d_kriging_cai_lst_mults_0_70_comb5_11272013.RData"
raster_obj_file_15 <- "raster_prediction_obj_gwr_CAI_dailyTmax_365d_gwr_cai_lst_mults_0_70_comb5_11222013.RData"

#List of object with 0% monthly hold out
list_raster_obj_files  <- list(file.path(in_dir1,raster_obj_file_1),file.path(in_dir2,raster_obj_file_2),
                               file.path(in_dir3a,raster_obj_file_3a),file.path(in_dir3b,raster_obj_file_3b),
                               file.path(in_dir4,raster_obj_file_4),file.path(in_dir5,raster_obj_file_5),
                               file.path(in_dir6,raster_obj_file_6),file.path(in_dir7,raster_obj_file_7),
                               file.path(in_dir8,raster_obj_file_8),file.path(in_dir9,raster_obj_file_9))
 
names(list_raster_obj_files)<- c("gam_daily","kriging_daily","gwr_daily_a","gwr_daily_b",
                                 "gam_CAI","kriging_CAI","gwr_CAI",
                                 "gam_fss","kriging_fss","gwr_fss")

y_var_name <- "dailyTmax"
out_prefix<-"analyses_05062014"
out_dir<-"/home/parmentier/Data/IPLANT_project/paper_multitime_scale__analyses_tables"
create_out_dir_param = TRUE

#Create output directory

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_prefix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

infile_reg_outline <- "/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/OR83M_state_outline.shp"  #input region outline defined by polygon: Oregon
met_stations_outfiles_obj_file<-"/data/project/layers/commons/data_workflow/output_data_365d_gam_fus_lst_test_run_07172013/met_stations_outfiles_obj_gam_fusion__365d_gam_fus_lst_test_run_07172013.RData"
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
ref_rast_name<- "/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/mean_day244_rescaled.rst"                     #This is the shape file of outline of the study area. #local raster name defining resolution, exent, local projection--. set on the fly??
infile_reg_outline <- "/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/OR83M_state_outline.shp"  #input region outline defined by polygon: Oregon
ref_rast_name <-"/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/mean_day244_rescaled.rst"  #local raster name defining resolution, exent: oregon
ref_rast_d001 <-"/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/mean_day001_rescaled.rst"
transect_list <-c("/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/transect_OR_1.shp",
                  "/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/transect_OR_2.shp",
                  "/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/transect_OR_3.shp")

#method_interpolation <- "gam_daily"
#covar_obj_file_1<- list.files(path=in_dir1,pattern="covar_obj.*.RData")
covar_obj_file_1 <- "covar_obj__365d_gam_daily_lst_comb5_11012013.RData"
met_obj_file_1 <- "met_stations_outfiles_obj_gam_daily__365d_gam_daily_lst_comb5_11012013.RData"

#Load objects containing training, testing, models objects 
met_stations_obj <- load_obj(file.path(in_dir1,met_obj_file_1))
covar_obj <-load_obj(file.path(in_dir1,covar_obj_file_1)) #Reading covariates object for GAM daily method
infile_covariates <- covar_obj$infile_covariates
infile_reg_outline <- covar_obj$infile_reg_outline
covar_names<- covar_obj$covar_names
#####
s_raster <- brick(infile_covariates)
names(s_raster)<-covar_names

###################################################################
############### BEGIN SCRIPT ##################################

###################### PART I: Generate tables for paper ##########
#Table 1: Covariates used (not produced in R)
#Table 2: Combinations of models tested (not produced in R)
#Table 3. Average RMSE for all models across year 2010 wiht 30% daily hold out
#Table 4. RMSE and monthly holdout using linear model with cat var
#Table 5: Average monhtly RMSE for CAI model 7, climatology surfaces,  CAI and single time scale methods 
#Table 6: Monthly correlations TMax-elev and TMax-LST  

##### Table 3: Contribution of covariates using validation accuracy metrics
## This is a table of accuracy  

summary_metrics_v_list <- lapply(list_raster_obj_files,FUN=function(x){x<-load_obj(x);x[["summary_metrics_v"]]$avg$rmse})                           
names(summary_metrics_v_list)


#for gam models
table_gam<- summary_metrics_v_list[grep("gam",names(summary_metrics_v_list))]
table_gam <-do.call(cbind,table_gam)
table_gam<-table_gam[1:7,]                           
#for kriging models
table_kriging <- summary_metrics_v_list[grep("kriging",names(summary_metrics_v_list))]
table_kriging <- do.call(cbind,table_kriging)
table_kriging <- table_kriging[1:7,]                         

#for gwr models
table_gwr <- summary_metrics_v_list[grep("gwr",names(summary_metrics_v_list))] #select models related to gwr
table_gwr$gwr_daily <- c(table_gwr[[1]],table_gwr[[2]]) #combine both gwr object results
table_gwr <- table_gwr[c(5,3,4)] #reorder list and drop first two objects
table_gwr <- do.call(cbind,table_gwr)
table_gwr <- table_gwr[1:7,]                         

table3_paper <- as.data.frame(do.call(rbind,list(table_gam,table_kriging,table_gwr)))    
table3_paper <- round(table3_paper,digit=3) #roundto three digits teh differences
table3_paper$Methods <- c(rep("gam",7),
                          rep("kriging",7),
                          rep("gwr",7))    
#Drop fss
#grep("gam_fss",names(table3_paper),invert=T,value=F)
table3_paper <- subset(table3_paper,
                       select=grep("gam_fss",names(table3_paper),invert=T,value=T))

list_tb <-lapply(list_raster_oj_files,FUN=function(x){x<-load_obj(x);x[["tb_diagnostic_v"]]})                           

stat<-"sd"
training <- "FALSE"

list_param_calc_stat <- list(list_tb,stat,training)
names(list_param_calc_stat) <- c("list_tb","stat","training")

#Calculate standard deviation for each metric
#undebug(calc_stat_from_tb_list)
#sd1 <- calc_stat_from_tb_list(1,list_param=list_param_calc_stat) # see function script

list_sd <- lapply(1:length(list_tb),FUN=calc_stat_from_tb_list,list_param=list_param_calc_stat)
names(list_sd) <- names(list_tb) #now select all the daily in one column and CAI in antoehr
#table_sd <- do.call(rbind,list_sd)
tab_sd <- do.call(rbind,list_sd)
tab_sd <- subset(tab_sd,pred_mod!="mod_kr")

table_sd <- as.data.frame(tab_sd[grep("daily",rownames(tab_sd),value=F),c("rmse")])
names(table_sd) <- c("rmse_s")
rmse_m <- (tab_sd[grep("CAI",rownames(tab_sd),value=F),c("rmse")])
table_sd$rmse_m <- rmse_m

#Combined sd in one table for mod1 (baseline 2)
#table_sd <- do.call(rbind,list(sd1[1,],sd3[1,],sd4[1,])) #table containing the sd for the three mdethods for baseline 2
table_sd <- round(table_sd,digit=3) #round to three digits the differences

#combined tables with accuracy metrics and their standard deviations
table_ac <- table3_paper[,1:2] #only use the rmse values
table_ac_combined <-table_combined_symbol(table_ac,table_sd,"±")
#lapply(lapply(table_ac,FUN=paste,table_sd,FUN=paste,sep="±"),FUN=paste)

table_ac_combined$Methods <- table3_paper$Methods
table3_paper <- table_ac_combined


#Check input covariates and model formula:
#list_formulas <-raster_prediction_obj_2$method_mod_obj[[1]]$formulas #formulas for models run comb5
list_formulas <-unlist(lapply(list_raster_obj_files[[1]],FUN=function(x){x<-load_obj(x);x$method_mod_obj[[1]]$formulas}))
#strsplit(list_formulas,"~")
                             
table3_paper$Formulas<-rep(list_formulas,3)                             
table3_paper<-table3_paper[(c(4,3,1,2))]  #reordering columns                           

#Testing siginificance between models                    
#mod_compk1 <-kruskal.test(tb1$rmse ~ as.factor(tb1$pred_mod)) #Kruskal Wallis test
#mod_compk2 <-kruskal.test(tb2$rmse ~ as.factor(tb2$pred_mod))
#print(mod_compk1) #not significant

###Table 3, writeout
names_table_col <-c("Models","Method","Single time scale","Multiple time scale CAI") # 

names(table3_paper)<- names_table_col
print(table3_paper) #show resulting table

#Now write out table 3

file_name<-paste("table3_multi_timescale_paper","_",out_prefix,".txt",sep="")
write.table(table3_paper,file=file_name,sep=",")

################ Table 4 ################
#### Test for monthly hold out  

#names(tb_v_list) <- paste("tb_v_",names(tb_v_list),sep="")
tb_v <- do.call(rbind.fill,tb_v_list)
#Create a table...
list_interp_method <- c("gam_CAI","kriging_CAI","gwr_CAI")
#loop through and extract p val and estimates...make this a function!!! to get a table
#tb_mv contains everything!!
list_param_extract_term <- list(list_interp_method,tb_v)
names(list_param_extract_term) <- c("list_interp_method","tb_mv")

#debug(extract_table_term_factor)
#test1 <- extract_table_term_factor(1,list_param=list_param_extract_term)

list_table_method_holdout_obj <-lapply(1:length(list_interp_method),FUN=extract_table_term_factor,list_param=list_param_extract_term)
list_table_method_holdout<- lapply(list_table_method_holdout_obj,
       FUN=function(x){x$table_method})
table4 <- do.call(rbind,list_table_method_holdout)
table4_paper <- table4[,c(9,1,2,3,4,5,6,7,8)]

file_name<-paste("table4_multi_timescale_paper","_",out_prefix,".txt",sep="")
write.table(table4_paper,file=file_name,sep=",")

################ Table 5 ################
#### Monthly RMSE for GAM CAI Climatology, GAM CAI daily and Single  time scale
## This is a table of accuracy  

#raster_prediction_obj_4 <-load_obj(file.path(in_dir4,raster_obj_file_4)) #comb5 gam_dir

tb_m_gam_CAI <- load_obj(file.path(in_dir13,raster_obj_file_13))$tb_month_diagnostic_v #getting objet
tb_gam_CAI <- load_obj(file.path(in_dir13,raster_obj_file_13))$tb_diagnostic_v #getting objet
tb_gam_daily <- load_obj(file.path(in_dir1,raster_obj_file_1))$tb_diagnostic_v #getting objet

names_mod <- paste("mod",1:7,sep="")
list_stat_tb_gam_CAI <-calc_stat_by_month_tb(names_mod,tb_gam_CAI,month_holdout=T)
list_stat_tb_gam_daily <-calc_stat_by_month_tb(names_mod,tb_gam_daily,month_holdout=F)

list_stat_tb_gam_CAI$avg
list_stat_tb_gam_CAI$sd

clim_rmse <- subset(tb_m_gam_CAI,prop==30 & pred_mod== "mod7",select="rmse")
CAI_rmse <- subset(list_stat_tb_gam_CAI$avg,prop_month==30 & pred_mod== "mod7",select="rmse")
daily_rmse <- subset(list_stat_tb_gam_daily$avg, pred_mod== "mod7",select="rmse")
#now get sd: no sd for clim
#clim_rmse_sd <- subset(tb_m_gam_CAI,prop==30 & pred_mod== "mod7",select="rmse")
#CAI_rmse_sd <- subset(list_stat_tb_gam_CAI$avg,prop_month==30 & pred_mod== "mod7",select="rmse")
#daily_rmse_sd <- subset(list_stat_tb_gam_daily$avg, pred_mod== "mod7",select="rmse")

table5 <- data.frame(clim=clim_rmse,CAI_month=CAI_rmse,daily_month=daily_rmse)
table5 <- round(table5,digit=3) #roundto three digits teh differences

table5$month <- month.abb
table5 <- table5[,c(4,1,2,3)]

###Table 5, writeout
names_table_col <-c("Month","Long term Climatology RMSE",
                    "Average Monthly RMSE from Daily pred (CAI)"
                    ,"Average Monthly RMSE from Daily pred (Single time scale)") # 

names(table5)<- names_table_col
print(table5) #show resulting table

#Now write out table 5

file_name<-paste("table5_multi_timescale_paper","_",out_prefix,".txt",sep="")
write.table(table5,file=file_name,sep=",")

################ Table 6 ################

### NOW EXMAINE LST-TMAx and elev cor
         
#Use data from accuracy with no monthly hold out
raster_obj <- load_obj(list_raster_obj_files$gam_CAI)
data_month <- (raster_obj$method_mod_obj[[1]]$data_month_s)
dates<-unique(raster_obj$tb_diagnostic_v$date)
#Extract monthly data frame used in fitting 
list_data_month_s <-extract_list_from_list_obj(raster_obj$validation_mod_month_obj,"data_s") 
year_nbs <- sapply(list_data_month_s,FUN=length) #number of observations per month
#Convert sp data.frame and combined them in one unique df, see function define earlier
data_month_all <-convert_spdf_to_df_from_list(list_data_month_s) #long rownames

names_tmp<-c("mm_01","mm_02","mm_03","mm_04","mm_05","mm_06","mm_07","mm_08",
             "mm_09","mm_10","mm_11","mm_12")#,"LST","elev_s")

list_cor_x <- vector("list",length=12)
for(i in 1:12){
  data_month<-list_data_month_s[[i]]
  data_x <- data_month[,c("TMax","elev_s",names_tmp[i])] #,"LC1","LC6","LC9")]
  data_x <- as.data.frame(data_x)[,1:3]
  cor_x<- cor(data_x,use="complete.obs")
  list_cor_x[[i]] <- cor_x
}

cor_tab<-lapply(list_cor_x,function(x){x[1,2:3]})
#cor_tab<-lapply(list_cor_x,function(x){x[1,2:4]})
         
cor_tab<-do.call(rbind,cor_tab)
table6 <- as.data.frame(cor_tab)
table6 <- round(table6,digit=3) #roundto three digits teh differences

table6$month <- month.abb
table6 <- table6[,c(3,1,2)]

###Table 6, writeout
names_table_col <-c("Month","Elevation-TMax",
                    "LST-TMax")
                
names(table6)<- names_table_col
print(table6) #show resulting table

#Now write out table 6

file_name<-paste("table6_correlation_multi_timescale_paper","_",out_prefix,".txt",sep="")
write.table(table6,file=file_name,sep=",")

#######################################################
####### PART 2: generate figures for paper #############

#figure 1: Study area Oregon state with GHCND stations
#Figure 2. Accuracy and monthly hold out for CAI procedures and GWR, Kriging and GAM methods.
#Figure 3. Overtraining tendency, difference between training and testing
#Figure 4: Spatial pattern of prediction for one day (3 maps)
#Figure 5: Spatial lag profiles for predicted surfaces 
#Figure 6: Granularity, Moran's I and standard deviation averaged monthly 
#Figure 7: Daily deviation and long term surfaces
#Figure 8: LST and air temperature, elevation and land cover across OR

#ANNEX : (end of script)

# Annex 1: LST averages daily and monthly for January 1 and January
# Annex 2: monthly holdout model
# Annex 3: Transect location (transect map)
# Annex 4: Transect profiles for one day 

################################################
######### Figure 1: Oregon study area, add labeling names to  Willamette Valley an Mountain Ranges?
#3 parameters from output
#infile_monthly<-list_outfiles$monthly_covar_ghcn_data #outile4 from database_covar script
#infile_daily<-list_outfiles$daily_covar_ghcn_data  #outfile3 from database_covar script
#infile_locs<- list_outfiles$loc_stations_ghcn #outfile2? from database covar script

ghcn_dat <- readOGR(dsn=dirname(met_stations_obj$monthly_covar_ghcn_data),
        sub(".shp","",basename(met_stations_obj$monthly_covar_ghcn_data)))
ghcn_dat_WGS84 <-spTransform(ghcn_dat,CRS_locs_WGS84)         # Project from WGS84 to new coord. system

interp_area <- readOGR(dsn=dirname(infile_reg_outline),sub(".shp","",basename(infile_reg_outline)))
interp_area_WGS84 <-spTransform(interp_area,CRS_locs_WGS84)         # Project from WGS84 to new coord. system

#usa_map <- getData('GADM', country='USA', level=1) #Get US map
usa_map <- getData('GADM', country='USA',level=1) #Get US map, this is not working right now

usa_map_2 <- usa_map[usa_map$NAME_1!="Alaska",] #remove Alaska
usa_map_2 <- usa_map_2[usa_map_2$NAME_1!="Hawaii",] #remove Hawai 
usa_map_OR <- usa_map_2[usa_map_2$NAME_1=="Oregon",] #get OR

elev <- subset(s_raster,"elev_s")
elev_WGS84 <- projectRaster(from=elev,crs=CRS_locs_WGS84,method="ngb")

#set up the output file to plot
res_pix<-960
col_mfrow<-1
row_mfrow<-1
png(filename=paste("Figure1_paper_study_area_",out_prefix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
par(mfrow=c(1,1))

#plot(elev_WGS84,cex.axis=2.1,cex.z=2.1)
plot(elev_WGS84,cex.axis=1.4,ylab="Latitude (degree)",xlab="Longitude (degree)",
     cex.lab=1.4,
     axis.args=list(at=seq(0,3500,500),
               labels=seq(0, 3500, 500), 
               cex.axis=1.4))
plot(interp_area_WGS84,add=T)
plot(ghcn_dat_WGS84,add=T)
title_text <-c("Elevation (m) and meteorological"," stations in Oregon")
legend("topleft",legend=title_text,cex=2.1,bty="n")

#this work on non sp plot too
arrow_position<-c(-118,41)
SpatialPolygonsRescale(layout.north.arrow(), offset = arrow_position, 
                       scale = 0.5, fill=c("transparent","black"),plot.grid=FALSE)
#note that scale in SpatialPolygonRescale sets the size of the north arrow!!

#Add region label
par(mar = c(0,0,0,0)) # remove margin
#opar <- par(fig=c(0.9,0.95,0.8, 0.85), new=TRUE)
#opar <- par(fig=c(0.85,0.95,0.8, 0.9), new=TRUE)
#opar <- par(fig=c(0.8,0.95,0.75, 0.9), new=TRUE)
#opar <- par(fig=c(0.75,0.95,0.7, 0.9), new=TRUE)
opar <- par(fig=c(0.65,0.9,0.8, 0.92), new=TRUE)

plot(usa_map_2,border="black") #border and lwd are options of graphics package polygon object
plot(usa_map_OR,col="dark grey",add=T)
box()
dev.off()

################################################
######### Figure 3. RMSE multi-timescale mulitple hold out for FSS and CAI

raster_prediction_obj_9 <-load_obj(file.path(in_dir9,raster_obj_file_9)) #comb5 gwr_fss

list_raster_obj_files_holdout  <- list(file.path(in_dir10,raster_obj_file_10),file.path(in_dir11,raster_obj_file_11),
                                  file.path(in_dir12,raster_obj_file_12),file.path(in_dir13,raster_obj_file_13),
                                  file.path(in_dir14,raster_obj_file_14),file.path(in_dir15,raster_obj_file_15))
 
names(list_raster_obj_files_holdout)<- c("gam_fss","kriging_fss","gwr_fss",
                                         "gam_CAI","kriging_CAI","gwr_CAI")
tb_v_list<-lapply(list_raster_obj_files_holdout,FUN=function(x){x<-load_obj(x);x[["tb_diagnostic_v"]]})                           
tb_s_list<-lapply(list_raster_obj_files_holdout,FUN=function(x){x<-load_obj(x);x[["tb_diagnostic_s"]]})                           
#tb_s_list<-mclapply(list_raster_obj_files_holdout,FUN=function(x){x<-load_obj(x);x[["tb_diagnostic_s"]]},mc.preschedule=FALSE,mc.cores = 6)                           
#tb_v_list<-mclapply(list_raster_obj_files_holdout,FUN=function(x){x<-load_obj(x);x[["tb_diagnostic_v"]]},mc.preschedule=FALSE,mc.cores = 6)                           
names(tb_s_list) <- paste("tb_s_",names(tb_s_list),sep="")
names(tb_v_list) <- paste("tb_v_",names(tb_v_list),sep="")

#Extract object element...
tb_mv_list<-lapply(list_raster_obj_files_holdout,FUN=function(x){x<-load_obj(x);x[["tb_month_diagnostic_v"]]})                           
tb_ms_list<-lapply(list_raster_obj_files_holdout,FUN=function(x){x<-load_obj(x);x[["tb_month_diagnostic_s"]]})                           
#tb_mv_list<-mclapply(list_raster_obj_files_holdout,FUN=function(x){x<-load_obj(x);x[["tb_month_diagnostic_v"]]},mc.preschedule=FALSE,mc.cores = 6)                           
#tb_ms_list<-mclapply(list_raster_obj_files_holdout,FUN=function(x){x<-load_obj(x);x[["tb_month_diagnostic_s"]]},mc.preschedule=FALSE,mc.cores = 6)                           
names(tb_ms_list) <- paste("tb_ms_",names(tb_ms_list),sep="") #monthly training accuracy
names(tb_mv_list) <- paste("tb_mv_",names(tb_mv_list),sep="") #monthly testing accuracy

list_tb <- c(tb_s_list,tb_v_list,tb_ms_list,tb_mv_list) #combined in one list
ac_metric <- "rmse"
#plot_names <- c("RMSE for gam_FSS","RMSE for kriging_FSS","RMSE for gwr_FSS",
#               "RMSE for gam_CAI","RMSE for kriging_CAI","RMSE for gwr_CAI")
names_mod <- paste("mod",1:7,sep="")
plot_names <- names(list_tb) #this is the default names for the plots
#now replace names for relevant figure used later on.
plot_names[7:12] <- c("RMSE for gam_FSS","RMSE for kriging_FSS","RMSE for gwr_FSS",
               "RMSE for gam_CAI","RMSE for kriging_CAI","RMSE for gwr_CAI")
#Quick function to explore accuracy make this a function to create solo figure...and run only a subset...

#list_plots <- plot_accuracy_by_holdout_fun(list_tb,ac_metric,plot_names,names_mod)
tb_obj <- plot_accuracy_by_holdout_fun(list_tb,ac_metric,plot_names,names_mod)

selected_df <- 7:12
names(tb_obj$list_avg_tb[selected_df]) #names of method  and validation set selected
avg_tb_ac <- do.call(rbind.fill,tb_obj$list_avg_tb[selected_df])

## We use only CAI results

avg_tb_ac <- subset(avg_tb_ac, method_interp %in% c("gam_CAI","kriging_CAI","gwr_CAI"))

#avg_tb_ac <- avg_tb_ac[, avg_tb_ac$method_interp %in% c("gam_CAI","kriging_CAI","gwr_CAI")]

layout_m <-c(1.25,2.5) #one row two columns

png(paste("Figure2_paper_accuracy_",ac_metric,"_prop_month","_",out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

p <- xyplot(rmse~prop_month|method_interp, # |set up pannels using method_interp
            group=pred_mod,data=avg_tb_ac, #group by model (covariates)
            main="",
            type="b",as.table=TRUE,
            index.cond=list(c(1,3,2)),    #this provides the order of the panels)
            pch=1:length(avg_tb_ac$pred_mod),
            par.settings=list(superpose.symbol = list(
            pch=1:length(avg_tb_ac$pred_mod),cex=2.3),
            axis.text = list(font = 2, cex = 2),
            par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
            par.strip.text=list(font=2,cex=2),
            #auto.key=list(columns=1,space="right",title="Model",cex=1.5),
            auto.key=list(columns=5,space="top",cex=2.5,font=2),
            #auto.key=list(columns=5),
            xlab=list(label="Monthly Hold out proportion",cex=2.5,font=2),
            ylab=list(label="RMSE (°C)",cex=2.5, font=2))

print(p) #need to use print to plot...

dev.off()

       
#Correlation tmax elev and LST...

################################################
######### Figure 3. RMSE multi-timescale mulitple hold out Overtraining tendency

#x<-tb_ms_list[[1]]
tb_ms_list_diff <- lapply(tb_ms_list,FUN=function(x){x[x$prop!=0,]})                           
tb_mv_list_diff <- lapply(tb_mv_list,FUN=function(x){x[x$prop!=0,]})                           

tb_s_list_diff <- lapply(tb_s_list,FUN=function(x){x[x$prop_month!=0,]})                           
tb_v_list_diff <- lapply(tb_v_list,FUN=function(x){x[x$prop_month!=0,]}) 

#For paper...
#Combine figures... tb_v for GWR, Kriging and GAM for both FSS and CAI

metric_names <- c("mae","rmse","me","r")
list_metric_names <- vector("list", length=6) #list(metric_names)
list_metric_names[[1]] <- metric_names
list_metric_names <-lapply(1:6,FUN=function(i,list_metric_names,metric_names){list_metric_names[[i]]<-metric_names},
       list_metric_names,metric_names)

##### Calculate differences between training and testing for all multi-timescale mehtods
list_diff <-lapply(1:6,FUN=list_diff_df_fun,list_tb_s=tb_s_list, list_tb_v=tb_v_list,list_of_list_metric_names=list_metric_names)  
names(list_diff) <- c("gam_fss","kriging_fss","gwr_fss",
                                         "gam_CAI","kriging_CAI","gwr_CAI")

#Check results
#metric_names <- c("mae","rmse","me","r")
#diff_kriging_CAI <- diff_df(list_tb[["tb_s_kriging_CAI"]],list_tb[["tb_v_kriging_CAI"]],metric_names)
#diff_gwr_CAI <- diff_df(list_tb[["tb_s_gwr_CAI"]],list_tb[["tb_v_gwr_CAI"]],metric_names)
#head(diff_gwr_CAI)
#head(list_diff[[6]]

list_m_diff <-lapply(1:6,FUN=list_diff_df_fun,list_tb_s=tb_ms_list_diff, list_tb_v=tb_mv_list_diff,list_of_list_metric_names=list_metric_names)  
names(list_m_diff) <- c("gam_fss","kriging_fss","gwr_fss",
                                         "gam_CAI","kriging_CAI","gwr_CAI")

## Now create boxplots...
layout_m<-c(1,2) #one row two columns

png(paste("Figure3_paper_boxplot_overtraining_",out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
#boxplot(diff_kriging_CAI$rmse,diff_gam_CAI$rmse,diff_gwr_CAI$rmse,names=c("kriging_CAI","gam_CAI","gwr_CAI"),
#        main="Difference between training and testing daily rmse")
par(mfrow=layout_m)

#monthly CAI
boxplot(list_m_diff$kriging_CAI$rmse,list_m_diff$gam_CAI$rmse,list_m_diff$gwr_CAI$rmse,
        ylim=c(-6,1), outline=TRUE,
        names=c("kriging_CAI","gam_CAI","gwr_CAI"),ylab="RMSE (°C)",xlab="Interpolation Method",
        main="Difference between training and testing for monhtly rmse")
legend("topleft",legend=c("a"),bty="n") #bty="n", don't put box around legend

#monthly fss
#boxplot(list_m_diff$kriging_fss$rmse,list_m_diff$gam_fss$rmse,list_diff$gwr_fss$rmse,
#        ylim=c(-9,1),outline=TRUE,
#        names=c("kriging_FSS","gam_FSS","gwr_FSS"),ylab="RMSE (°C)",xlab="Interpolation Method",
#        main="Difference between training and testing for monhtly rmse")
#legend("topleft",legend=c("c"),bty="n")

#daily CAIdras
boxplot(list_diff$kriging_CAI$rmse,list_diff$gam_CAI$rmse,list_diff$gwr_CAI$rmse,
        ylim=c(-6,1.5),outline=TRUE,
        names=c("kriging_CAI","gam_CAI","gwr_CAI"),ylab="RMSE (°C)",xlab="Interpolation Method",
        main="Difference between training and testing daily rmse")
legend("topleft",legend=c("b"),bty="n")


#daily fss
#boxplot(list_diff$kriging_fss$rmse,list_diff$gam_fss$rmse,list_diff$gwr_fss$rmse,
#        ylim=c(-12,1.5),outline=TRUE,
#        names=c("kriging_FSS","gam_FSS","gwr_FSS"),ylab="RMSE (°C)",xlab="Interpolation Method",
#        main="Difference between training and testing daily rmse")
#legend("topleft",legend=c("d"),bty="n")

dev.off()


################################################
######### Figure 4: Spatial pattern of prediction for one day (maps)

y_var_name <-"dailyTmax"
index<-244 #index corresponding to Sept 1

lf_list<-lapply(list_raster_obj_files[c("gam_daily","gam_CAI","gam_fss")],
                               FUN=function(x){x<-load_obj(x);x$method_mod_obj[[index]][[y_var_name]]})                           

date_selected <- "20109101"
#methods_names <-c("gam","kriging","gwr")
methods_names <-c("gam_daily","gam_CAI","gam_FSS")

names_layers<-methods_names[1:2]
#lf <- (list(lf1,lf4[1:7],lf7[1:7]))
#lf<-list(lf_list[[1]],lf_list[[2]][1:7],lf_list[[3]][1:7])
lf<-list(lf_list[[1]],lf_list[[2]][1:7])

names_layers <-c("mod1 = lat*long","mod2 = lat*long + LST","mod3 = lat*long + elev","mod4 = lat*long + elev + N_w*E_w",
                 "mod5 = lat*long + elev + DISTOC","mod6 = lat*long + elev + LST","mod7 = lat*long + elev + LST*FOREST")
nb_fig<- c("4a","4b")
list_plots_spt <- vector("list",length=length(lf))
for (i in 1:length(lf)){
  pred_temp_s <-stack(lf[[i]])

  #lf2 <- raster_prediction_obj_2$method_mod_obj[[index]][[y_var_name]]
  #lf2 #contains the models for gam

  #s.range <- c(min(minValue(pred_temp_s)), max(maxValue(pred_temp_s)))
  #s.range <- s.range+c(5,-5)
  #col.breaks <- pretty(s.range, n=200)
  #lab.breaks <- pretty(s.range, n=100)
  
  #max_val<-s.range[2]
  #min_val <-s.range[1]
  max_val <- 40
  min_val <- -10
  layout_m <- c(2,4) #one row two columns
  no_brks <- length(seq(min_val,max_val,by=0.1))-1
  #temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
  #temp.colors <- colorRampPalette(c('blue', 'lightgoldenrodyellow', 'red'))
  #temp.colors <- matlab.like(no_brks)
  temp.colors <- colorRampPalette(c('blue', 'khaki', 'red'))
  
  png(paste("Figure_",nb_fig[i],"_spatial_pattern_tmax_prediction_models_gam_levelplot_",date_selected,out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

  p <- levelplot(pred_temp_s,main=methods_names[i], 
                 ylab=NULL,xlab=NULL,
          par.settings = list(axis.text = list(font = 2, cex = 2),layout=layout_m,
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),par.strip.text=list(font=2,cex=2),
          names.attr=names_layers,
                 col.regions=temp.colors(no_brks),at=seq(min_val,max_val,by=0.1))
  #col.regions=temp.colors(25))
  print(p)
  dev.off()
  list_plots_spt[[i]] <- p
}

layout_m <- c(2,4) # works if set to this?? ok set the resolution...
#layout_m<-c(2*3,4) # works if set to this?? ok set the resolution...

png(paste("Figure4_paper_","_spatial_pattern_tmax_prediction_models_gam_levelplot_",date_selected,out_prefix,".png", sep=""),
    height=960*layout_m[1],width=960*layout_m[2])
    #height=480*6,width=480*4)

p1 <- list_plots_spt[[1]]
p2 <- list_plots_spt[[2]]
#p3 <- list_plots_spt[[3]]

#grid.arrange(p1,p2,p3,ncol=1)
grid.arrange(p1,p2,ncol=1)

dev.off()

################################################
#### Figure 5: Spatial lag profiles  
#This figure is generated to show the spatial Moran'I for 10 spatial 
#for Jan 1 and Sept 1 in 2010 for all models (1 to 7) and methods

index <- 1 #index corresponding to Jan 1 #For now create Moran's I for only one date...
lf_moran_list_date1 <-lapply(list_raster_obj_files[c("gam_daily","gam_CAI","gam_fss")],
                               FUN=function(x){x<-load_obj(x);x$method_mod_obj[[index]][[y_var_name]]})                           
index <- 244 #index corresponding to Sept 1 #For now create Moran's I for only one date...
lf_moran_list_date2 <-lapply(list_raster_obj_files[c("gam_daily","gam_CAI","gam_fss")],
                               FUN=function(x){x<-load_obj(x);x$method_mod_obj[[index]][[y_var_name]]})                           

#date_selected <- "20100101"
#methods_names <-c("gam","kriging","gwr")
methods_names <-c("gam_daily","gam_CAI","gam_FSS") #drop gam FSS later

names_layers<-methods_names
#Subset images to eliminate mod_kr
lf1 <- list(lf_moran_list_date1[[1]],lf_moran_list_date1[[2]][1:7]) #,lf_moran_list_date1[[3]][1:7])
lf2 <- list(lf_moran_list_date2[[1]],lf_moran_list_date2[[2]][1:7]) #,lf_moran_list_date2[[3]][1:7])
names(lf1)<-c("gam_daily","gam_CAI") #,"gam_FSS")
names(lf2)<-c("gam_daily","gam_CAI") #,"gam_FSS")

### Now extract Moran's I for a range of lag using a list of images

#set maximum lag range
nb_lag <-10
#Provide list of raster images:
list_lf <- list(lf1,lf2)

list_moran_df1 <- calculate_moranI_profile(list_lf[[1]],nb_lag) #for January 1
list_moran_df2 <- calculate_moranI_profile(list_lf[[2]],nb_lag) #for September 1
names(list_moran_df1)<-c("gam_daily","gam_CAI")#,"gam_FSS")
names(list_moran_df2)<-c("gam_daily","gam_CAI")#,"gam_FSS")

#Run accross two dates...
#list_moran lapply(list_lf,FUN=calculate_moranI_profile,nb_lag=nb_lag)

### Prepare to plot lag Moran's I profiles

#generate automatic title for exploration if necessary!!
list_title_plot<- list(c("Spatial lag profile on January 1, 2010"),
                  c("Spatial lag profile on September 1, 2010"))
#name used in the panel!!!
names_panel_plot <-c("mod1 = lat*long","mod2 = lat*long + LST","mod3 = lat*long + elev","mod4 = lat*long + elev + N_w*E_w",
                 "mod5 = lat*long + elev + DISTOC","mod6 = lat*long + elev + LST","mod7 = lat*long + elev + LST*FOR")
layout_m<-c(2,4) # works if set to this?? ok set the resolution...
list_moran_df <- list(list_moran_df1,list_moran_df2)
list_param_plot_moranI_profile_fun <- list(list_moran_df,list_title_plot,names_panel_plot,layout_m)
names(list_param_plot_moranI_profile_fun) <- c("list_moran_df","list_title_plot","names_panel_plot","layout_m")

#debug(plot_moranI_profile_fun)
#p<- plot_moranI_profile_fun(1,list_param=list_param_plot_moranI_profile_fun)
  
list_moranI_plots <- lapply(1:2,FUN=plot_moranI_profile_fun,list_param=list_param_plot_moranI_profile_fun)


#This function uses list moran_df object from calculate_moranI_profile function!!

#layout_m<-c(2,4) # works if set to this?? ok set the resolution...
#layout_m<-c(2*3,4) # works if set to this?? ok set the resolution...

png(paste("Figure5_paper_spatial_correlogram_prediction_models_levelplot_",out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
    #height=3*480*layout_m[1],width=2*480*layout_m[2])
    #height=480*6,width=480*4)
#png(paste("Figure11_paper_spatial_correlogram_prediction_models_levelplot_",out_prefix,".png", sep=""),
#    height=480,width=480)
    #height=480*6,width=480*4)

#p1 <- list_moranI_plots[[1]]
p2 <- list_moranI_plots[[2]]
#grid.arrange(p2,ncol=1)
print(p2)
#grid.arrange(p1,p2,ncol=1)
dev.off()

##################################################
####### Figure 6: Granularity -MoranI and std averaged monthly over 2010  

#get the  list of all prediction for a specific method
#list_lf <-extract_list_from_list_obj(load_obj(list_raster_obj_files$gam_CAI)$method_mod_obj,"dailyTmax") #getting objet
list_lf_gam_CAI <-extract_list_from_list_obj(load_obj(list_raster_obj_files[["gam_CAI"]])$method_mod_obj,"dailyTmax") #getting objet
#list_lf_gam_fss <-extract_list_from_list_obj(load_obj(list_raster_obj_files[["gam_fss"]])$method_mod_obj,"dailyTmax") #getting objet
list_lf_gam_daily <-extract_list_from_list_obj(load_obj(list_raster_obj_files[["gam_daily"]])$method_mod_obj,"dailyTmax") #getting objet

nb_lag <- 10
list_filters<-lapply(1:nb_lag,FUN=autocor_filter_fun,f_type="queen") #generate lag 10 filters
list_param_stat_moran_CAI <- list(filter=list_filters[[10]],lf_list=list_lf_gam_CAI)
list_param_stat_moran_gam_daily <- list(filter=list_filters[[10]],lf_list=list_lf_gam_daily)
#tt <- stat_moran_std_raster_fun(1,list_param=list_param_stat_moran)
#tt <- mclapply(1:11, list_param=list_param_stat_moran_CAI, FUN=stat_moran_std_raster_fun,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
tt_gam_daily <- mclapply(1:365, list_param=list_param_stat_moran_gam_daily, FUN=stat_moran_std_raster_fun,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
save(tt_gam_daily,file=file.path(out_dir,"moran_std_tt_gam_daily.RData"))

#Takes 1 hour to get the average moran's I for the whole year so load moran_std_tt_fss2.RData
#list_param_stat_moran <- list(filter=list_filters[[10]],lf_list=list_lf_gam_fss)
#tt_fss <- mclapply(1:365, list_param=list_param_stat_moran, FUN=stat_moran_std_raster_fun,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
#tt_gam_daily <- mclapply(1:365, list_param=list_param_stat_moran_dam_daily, FUN=stat_moran_std_raster_fun,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement

#save(tt_fss,file=file.path(out_dir,"moran_std_tt_fss.RData"))
#list_param_stat_moran <- list(filter=list_filters[[10]],lf_list=list_lf_gam_CAI)
tt_gam_CAI <- mclapply(1:365, list_param=list_param_stat_moran_CAI, FUN=stat_moran_std_raster_fun,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
save(tt_gam_CAI,file=file.path(out_dir,"moran_std_tt_gam_CAI.RData"))

tx1 <- do.call(rbind,tt_gam_daily)
tx2 <- do.call(rbind,tt_gam_CAI)
#tt_fss2 <- load_obj("moran_std_tt_fss2.RData")
#tx<- do.call(rbind,tt_fss2)
list_tb_std_moran <- list(tx1,tx2)
list_plot_std_moranI <- vector("list",length=2)

for(j in 1:length(list_tb_std_moran)){
  tb_std_moran <- list_tb_std_moran[[j]]
  dates <- unique(raster_obj$tb_diagnostic_v$date)
  dates_proc <- strptime(dates, "%Y%m%d")   # interpolation date being processed
  dates_proc <- as.data.frame(dates_proc)
  dates_proc$index <- 1:365
  tx <- merge(tb_std_moran,dates_proc,by="index")
  tx$month <- as.integer(strftime(tx$dates_proc, "%m"))          # current month of the date being processed
  names(tx) <- c("index","moranI","std","pred_mod","date","month")
  t <- melt(tx,
            #measure=mod_var, 
            id=c("pred_mod","month"),
            na.rm=F)
  #t$value<-as.numeric(t$value) #problem with char!!!
  tb_mod_m_avg2 <-cast(t,pred_mod+month~variable,mean) #monthly mean for every model
  
  #mo<-as.integer(strftime(date_proc, "%m"))          # current month of the date being processed
  #day<-as.integer(strftime(date_proc, "%d"))
  #year<-as.integer(strftime(date_proc, "%Y"))
  # end of pasted
  
  tb_mod_m_avg <- subset(tb_mod_m_avg2,pred_mod!=c("mod_kr"))
  
  p_moranI <- xyplot(moranI~month, # |set up pannels using method_interp
                     group=pred_mod,data=tb_mod_m_avg, #group by model (covariates)
                     main="Moran's I in daily predictions average by model and month",
                     type="b",as.table=TRUE,
                     #index.cond=list(c(1,3,2)),    #this provides the order of the panels)
                     pch=1:length(tb_mod_m_avg$pred_mod),
                     par.settings=list(superpose.symbol = list(
                       pch=1:length(tb_mod_m_avg$pred_mod))),
                     auto.key=list(columns=1,space="right",title="Model",cex=1),
                     #auto.key=list(columns=5),
                     xlab="Month",
                     ylab="Moran I")
  #print(p)
  
  p_std <- xyplot(std~month, # |set up pannels using method_interp
                  group=pred_mod,data=tb_mod_m_avg, #group by model (covariates)
                  main="Standard devation in daily predictions average by model and month",
                  type="b",as.table=TRUE,
                  #index.cond=list(c(1,3,2)),    #this provides the order of the panels)
                  pch=1:length(tb_mod_m_avg$pred_mod),
                  par.settings=list(superpose.symbol = list(
                    pch=1:length(tb_mod_m_avg$pred_mod))),
                  auto.key=list(columns=1,space="right",title="Model",cex=1),
                  #auto.key=list(columns=5),
                  xlab="Month",
                  ylab="Standard Deviation")
  list_plot_std_moranI[[j]] <- list(p_moranI,p_std)  
  
}

#xyplot(moranI~month,data=tb_mod_m_avg,groups=pred_mod,type="b")
#title("Spatial variation in FSS surfaces")
#layout_m <- c(1,2)
layout_m <- c(2,2)
png(paste("Figure6_paper_moranI_std__predictions_models_levelplot_",out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
    #height=3*480*layout_m[1],width=2*480*layout_m[2])
    #height=480*6,width=480*4)
#png(paste("Figure11_paper_spatial_correlogram_prediction_models_levelplot_",out_prefix,".png", sep=""),
#    height=480,width=480)
    #height=480*6,width=480*4)
p3 <- list_plot_std_moranI[[1]][[1]] #moran I for gam daily
p4 <- list_plot_std_moranI[[1]][[2]] #std for gam daily
p1 <- list_plot_std_moranI[[2]][[1]] #moran I for gam CAI
p2 <- list_plot_std_moranI[[2]][[2]] #std for gam CAI

p1$main <- "Moran I in daily predictions averaged by month"
p2$main <- "Standard devation in daily predictions averaged by month"
p3$main <- "Moran I devation in daily predictions average by month"
p4$main <- "Standard devation in daily predictions average by month"

p1$legend$right$args$title <- "gam_CAI"
p2$legend$right$args$title <- "gam_CAI"
p3$legend$right$args$title <- "gam_s"
p4$legend$right$args$title <- "gam_s"
#p1 <- list_moranI_plots[[1]]
#p2 <- list_moranI_plots[[2]]
#grid.arrange(p_moranI,p_std,ncol=2)

grid.arrange(p1,p2,p3,p4,ncol=2)

dev.off()

################################################
#Figure 7: Monthly climatology, Daily deviation 

#The list of files for 365 days in 2010...
#y_var_name <- dailyTmax
lf_delta_gam_CAI <-extract_list_from_list_obj(load_obj(list_raster_obj_files[["gam_CAI"]])$method_mod_obj,"delta") #getting objet
lf_daily_gam_CAI <-extract_list_from_list_obj(load_obj(list_raster_obj_files[["gam_CAI"]])$method_mod_obj,y_var_name) #getting objet
lf_clim_gam_CAI <-extract_list_from_list_obj(load_obj(list_raster_obj_files[["gam_CAI"]])$method_mod_obj,"clim") #getting objet

index <- 289
temp_day_s <- stack(lf_clim_gam_CAI[[index]]$mod7,lf_delta_gam_CAI[[index]]$mod7,lf_daily_gam_CAI[[index]]$mod7)
names(temp_day_s) <- c("Monhtly Climatology","Daily Devation","Daily Prediction")


no_brks=25
temp.colors <- colorRampPalette(c('blue', 'khaki', 'red'))

names_layers <- c("Monthly Climatology","Daily deviation","Daily prediction")
fig_nb <- c("7a","7b","7c")
list_p <- vector("list",length=length(names_layers))

for(i in 1:length(names_layers)){
  p <- levelplot(temp_day_s,layer=i, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),par.strip.text=list(font=2,cex=2),
                 main=names_layers[i],col.regions=temp.colors(no_brks))
  
  png(paste("Figure",fig_nb[i],"_paper_climatology_daily_deviation_daily_prediction_model7_levelplot_",out_prefix,".png", sep=""),
    height=480*1,width=480*1)
  print(p) #to plot in a loop!!  
  dev.off()
  
  list_p[[i]] <- p

}

layout_m <- c(1.4,3.2)
png(paste("Figure7_paper_climatology_daily_deviation_daily_prediction_model7_levelplot_",out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
    #height=3*480*layout_m[1],width=2*480*layout_m[2])
    #height=480*6,width=480*4)

#,legend.shrink=2)
grid.arrange(list_p[[1]],list_p[[2]] ,list_p[[3]] ,ncol=3)

dev.off()

######################################
####### Figure 8 : Disussion section

############################
#Get Tmax and LST averages

names_tmp<-c("mm_01","mm_02","mm_03","mm_04","mm_05","mm_06","mm_07","mm_08",
             "mm_09","mm_10","mm_11","mm_12")
LST_s<-subset(s_raster,names_tmp)
names_tmp<-c("nobs_01","nobs_02","nobs_03","nobs_04","nobs_05","nobs_06","nobs_07","nobs_08",
             "nobs_09","nobs_10","nobs_11","nobs_12")
LST_nobs<-subset(s_raster,names_tmp)
plot(LST_s)
plot(LST_nobs)
y_range_nobs <- c(0,385)
y_range_avg <- c(-15,60)
var_name <- "LST"

#It works for any variable with stack of monthly values...!!
LST_stat <- plot_mean_nobs_r_stack_by_month(var_s=LST_s,var_nobs=LST_nobs,y_range_nobs,y_range_avg,var_name,out_prefix)

## open device for plotting
layout_m<-c(2,2) # works if set to this?? ok set the resolution...

png(paste("Figure8_paper_TMax_and_LST_relationship_",out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
par(mfrow=layout_m)    

#######################################
####### Figure 8a: LST and monthly air temperature  

#Use data from accuracy with no monthly hold out
raster_obj <- load_obj(list_raster_obj_files$gam_CAI)
data_month <- (raster_obj$method_mod_obj[[1]]$data_month_s)
dates<-unique(raster_obj$tb_diagnostic_v$date)
#Extract monthly data frame used in fitting 
list_data_month_s <-extract_list_from_list_obj(raster_obj$validation_mod_month_obj,"data_s") 
year_nbs <- sapply(list_data_month_s,FUN=length) #number of observations per month
#Convert sp data.frame and combined them in one unique df, see function define earlier
data_month_all <-convert_spdf_to_df_from_list(list_data_month_s) #long rownames
#LSTD_bias_avgm<-aggregate(LSTD_bias~month,data=data_month_all,mean)
#LSTD_bias_sdm<-aggregate(LSTD_bias~month,data=data_month_all,sd)
LST_avgm <- aggregate(LST~month,data=data_month_all,mean)
LST_avgm_min <- aggregate(LST~month,data=data_month_all,min)
LST_avgm_max <- aggregate(LST~month,data=data_month_all,max)
TMax_avgm <- aggregate(TMax~month,data=data_month_all,mean)

#plot(LST_avgm,type="b")
#lines(TMax_avgm,col="blue",add=T)
#plot(data_month_all$month,data_month_all$LST, xlab="month",ylab="LST at station")                  

statistics_LST_s<- LST_stat$avg #extracted earlier!!!

##Now plot Figure 8a
plot(TMax_avgm,type="b",ylim=c(-10,50),col="red",xlab="month",
     ylab="tmax (degree C)",cex=2,cex.lab=1.5) #,cex.axis=1.8)
lines(1:12,LST_avgm$LST,type="b",col="green",lty="dashed",pch=2,cex=2)
#lines(1:12,statistics_LST_s$mean,type="b",col="darkgreen")
lines(1:12,LST_avgm_min$LST,type="b",col="black",lty="dashed",pch=3,cex=2)
lines(1:12,LST_avgm_max$LST,type="b",col="black",lty="dashed",pch=4,cex=2)
#text(1:12,LST_avgm$LST,month.abb,cex=0.8,pos=2)
legend("topleft",legend=c("TMax","LST","LST min","LST max"), cex=1.5, 
       col=c("red","green","black","black"),
              lty=c("solid","dashed","dashed","dashed"),pch=1:4,bty="n")
#legend("topright",legend=c("TMax","LST","LST min","LST max"), cex=1.5, 
#       col=c("red","green","black","black"),
#              lty=c("solid","dashed","dashed","dashed"),pch=1:4,bty="n")
title("Monthly average LST and tmax at stations in Oregon",cex=2.4)           

##############################
####### Figure 8b: spatial variation -MoranI and  std  
#LST for area with Forest 50> and grass >50

#Account for forest
data_month_all_forest<-data_month_all[data_month_all$LC1>=50,]
data_month_all_grass<-data_month_all[data_month_all$LC6>=50,]
#data_month_all_grass<-data_month_all[data_month_all$LC6>=10,]
data_month_all_urban<-data_month_all[data_month_all$LC9>=50,]

LST_avgm_forest <-aggregate(LST~month,data=data_month_all_forest,mean)
LST_avgm_grass <-aggregate(LST~month,data=data_month_all_grass,mean)
LST_avgm_urban <-aggregate(LST~month,data=data_month_all_urban,mean)

##Now plot Figure 8b
plot(TMax_avgm,type="b",ylim=c(-7,50),col="red",xlab="month",
     ylab="tmax (degree C)",cex=2,cex.lab=1.5)#,cex.axis=1.8)
lines(1:12,LST_avgm$LST,type="b",col="green",lty="dashed",pch=2,cex=2)
lines(LST_avgm_forest,type="b",col="black",lty="dotted",pch=3,cex=2)
lines(LST_avgm_grass,type="b",col="black",lty="dotdash",pch=4,cex=2)
#lines(LST_avgm_urban,type="b",col="black",lty="longdash",pch=4)
legend("topleft",legend=c("TMax","LST","LST_forest","LST_grass"), cex=1.5, 
       col=c("red","green","black","black","black"),
       lty=1:4,pch=1:4,bty="n")
title("Monthly average LST and land cover types for stations ",cex=2.4)
       
########################################
####### Figure 8c: correlation between elevation and LST

LST_s <- subset(s_raster,c("mm_01","mm_02","mm_03","mm_04","mm_05","mm_06","mm_07","mm_08","mm_09","mm_10","mm_11","mm_12"))
elev <- subset(s_raster,"elev_s")
LC1 <- subset(s_raster,"LC1")

#### Now elev

r_tmp1 <-stack(LST_s,elev)
names(r_tmp1) <- c("mm_01","mm_02","mm_03","mm_04","mm_05","mm_06","mm_07","mm_08","mm_09","mm_10","mm_11","mm_12","elev")
x_cor_tmp1 <-layerStats(r_tmp1,stat="pearson",na.rm=T)

r_tmp2 <-stack(LST_s,LC1)
names(r_tmp2) <- c("mm_01","mm_02","mm_03","mm_04","mm_05","mm_06","mm_07","mm_08","mm_09","mm_10","mm_11","mm_12","LC1")
x_cor_tmp2 <-layerStats(r_tmp2,stat="pearson",na.rm=T)

##Now plot Figure 8d
plot(1:12,x_cor_tmp1[[1]][1:12,13],ylim=c(-0.9,0.9),
     type="b",lty="solid",pch=1,
     ylab="Pearson corelation",xlab="month",
     cex=2,cex.lab=1.5)#,cex.axis=1.8)
lines(1:12,x_cor_tmp2[[1]][1:12,13],type="b",lty="dashed",pch=2,cex=2)

legend("topleft",legend=c("Correlation LST-Elevation","Correlation LST-Forest"), 
       cex=1.5, col=c("black","black"),
       lty=1:2,pch=1:2,bty="n")
title("Correlation between LST-elevation and LST-Forest",cex=2.4)
       
########################################
####### Figure 8d: LST, TMax and RMSE for GAM_CAI (monthly averages)

##Now plot Figure 8d
plot(TMax_avgm,type="b",ylim=c(-7,40),col="red",xlab="month",
     ylab="tmax (degree C)",cex=2,cex.lab=1.5)#,cex.axis=1.8)
lines(1:12,LST_avgm$LST,type="b",col="green",lty="dashed",pch=2,cex=2)
par(new=TRUE)              # key: ask for new plot without erasing old
plot(1:12,table5[,2],type="b",pch=3,col="black",xlab="",ylab="",lty="dotted",
     ,ylim=c(0.5,2),axes=F) #plotting fusion profile
     #axis(4,xlab="",ylab="elevation(m)")  
     axis(4,cex=2)

#lines(LST_avgm_grass,type="b",col="black",lty="dotdash",pch=4)
#lines(LST_avgm_urban,type="b",col="black",lty="longdash",pch=4)
legend("topleft",legend=c("TMax","LST","RMSE Clim"), cex=1.5, 
       col=c("red","green","black"),
       lty=1:3,pch=1:3,bty="n")
title("RMSE for GAM_CAI Climatology in relation to LST and TMax",
      cex=2.4)

#closing file for figure 8
dev.off()

################################################
############ GENERATE ANNEX FIGURES 

##########################
#######Figure Annex 1: LST averaging: daily mean compared to monthly mean

lst_md <- raster(ref_rast_name)
projection(lst_md)<-projection(s_raster)
lst_mm_09<-subset(s_raster,"mm_09")

lst_md<-raster(ref_rast_d001)
lst_md<- lst_md - 273.16
lst_mm_01<-subset(s_raster,"mm_01")

png(filename=paste("Figure_annex1_paper_Comparison_daily_monthly_mean_lst",out_prefix,".png",sep=""),width=960,height=480)
par(mfrow=c(1,2))
plot(lst_md)
plot(interp_area,add=TRUE)
title("Mean for January 1")
plot(lst_mm_01)
plot(interp_area,add=TRUE)
title("Mean for month of January")
dev.off()

############################################
######Figure annex 2: Map of transects

## Transects image location in OR  
layout_m<-c(1,1)
png(paste("Figure_annex2_paper_elevation_transect_paths_",out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
vx_text <- c(395770.1,395770.1,395770.1) #location of labels
vy_text  <-c(473000,297045.9,112611.5)

plot(elev)
for(i in 1:length(transect_list)){
  filename<-sub(".shp","",transect_list[i])             #Removing the extension from file.
  transect<-readOGR(dirname(filename), basename(filename))                 #reading shapefile 
  #transect2 <- elide(transect, shift=c(0, -24000))
  #plot(transect2,add=TRUE)
  #writeOGR(transect2,dsn=".",layer="transect_OR_1_12282013",driver="ESRI Shapefile")
  plot(transect,add=TRUE)
}
text(x=vx_text,y=vy_text,labels=c("Transect 1","Transect 2","Transect 3"))
title("Transect Oregon")
dev.off()

##############################
######Figure annex 3: TRANSECT PROFILES

nb_transect <- 3
list_transect2<-vector("list",nb_transect)
list_transect3<-vector("list",nb_transect)
list_transect4<-vector("list",nb_transect)

#names_layers <-c("mod1 = lat*long","mod2 = lat*long + LST","mod3 = lat*long + elev","mod4 = lat*long + N_w*E_w",
#                 "mod5 = lat*long + elev + DISTOC","mod6 = lat*long + elev + LST","mod7 = lat*long + elev + LST*FOREST")
y_var_name <-"dailyTmax"
index<-244 #index corresponding to Sept 1

lf_list<-lapply(list_raster_obj_files[c("gam_daily","gam_CAI","gam_fss")],
                               FUN=function(x){x<-load_obj(x);x$method_mod_obj[[index]][[y_var_name]]})                           

date_selected <- "20109101"
#methods_names <-c("gam","kriging","gwr")
methods_names <-c("gam_daily","gam_CAI","gam_FSS")

names_layers<-methods_names[1:2]
#lf <- (list(lf1,lf4[1:7],lf7[1:7]))
#lf<-list(lf_list[[1]],lf_list[[2]][1:7],lf_list[[3]][1:7])
lf<-list(lf_list[[1]],lf_list[[2]][1:7])

rast_pred<-stack(lf[[2]]) #GAM_CAI
rast_pred_selected2<-subset(rast_pred,c(1,2,6)) #3 is referring to FSS, plot it first because it has the
rast_pred_selected3<-subset(rast_pred,c(1,3,6)) #3 is referring to FSS, plot it first because it has the
                                          # the largest range.
rast_pred2 <- stack(rast_pred_selected2,subset(s_raster,"elev_s"))
rast_pred3 <- stack(rast_pred_selected3,subset(s_raster,"elev_s"))

#layers_names<-layerNames(rast_pred2)<-c("lat*lon","lat*lon + elev + LST","elev")
layers_names2 <- names(rast_pred2)<-c("mod1","mod2","mod6","elev")
layers_names3 <- names(rast_pred3)<-c("mod1","mod3","mod6","elev")
#pos<-c(1,2) # postions in the layer prection
#transect_list
list_transect2[[1]]<-c(transect_list[1],paste("Figure_annex3a_paper_tmax_elevation_transect1_OR_",date_selected,
                                           paste("mod1_2_6",collapse="_"),out_prefix,sep="_"))
list_transect2[[2]]<-c(transect_list[2],paste("Figure_annex3b_tmax_elevation_transect2_OR_",date_selected,
                                           paste("mod1_2_6",collapse="_"),out_prefix,sep="_"))
list_transect2[[3]]<-c(transect_list[3],paste("Figure_annex3c_paper_tmax_elevation_transect3_OR_",date_selected,
                                           paste("mod1_2_6",collapse="_"),out_prefix,sep="_"))

list_transect3[[1]]<-c(transect_list[1],paste("Figure_annex3a_tmax_elevation_transect1_OR_",date_selected,
                                           paste("mod1_3_6",collapse="_"),out_prefix,sep="_"))
list_transect3[[2]]<-c(transect_list[2],paste("Figure_annex3b_tmax_elevation_transect2_OR_",date_selected,
                                           paste("mod1_3_6",collapse="_"),out_prefix,sep="_"))
list_transect3[[3]]<-c(transect_list[3],paste("Figure_annex3c_tmax_elevation_transect3_OR_",date_selected,
                                           paste("mod1_3_6",collapse="_"),out_prefix,sep="_"))

names(list_transect2)<-c("Oregon Transect 1","Oregon Transect 2","Oregon Transect 3")
names(list_transect3)<-c("Oregon Transect 1","Oregon Transect 2","Oregon Transect 3")

names(rast_pred2)<-layers_names2


names(rast_pred3)<-layers_names3

title_plot2<-paste(names(list_transect2),"on",date_selected,sep=" ")
#title_plot2<-paste(names(list_transect2),"on",date_selected,sep=" ")

#title_plot2<-paste(rep("Oregon transect on ",3), date_selected,sep="")
#title_plot3<-paste(names(list_transect3),date_selected,sep=" ")
#title_plot3<-paste(rep("Oregon transect on ",3), date_selected,sep="")

#r_stack<-rast_pred
#m_layers_sc<-c(3) #elevation in the third layer
m_layers_sc<-c(4) #elevation in the third layer

#title_plot2
#rast_pred2
#debug(plot_transect_m2)
trans_data2 <- plot_transect_m2(list_transect2,rast_pred2,title_plot2,disp=FALSE,m_layers_sc)
trans_data3 <- plot_transect_m2(list_transect3,rast_pred3,title_plot2,disp=FALSE,m_layers_sc)

         
         
#### DATA MONTH
         

###################### END OF SCRIPT #######################

# #LAND COVER INFORMATION

# LC1: Evergreen/deciduous needleleaf trees
# LC2: Evergreen broadleaf trees
# LC3: Deciduous broadleaf trees
# LC4: Mixed/other trees
# LC5: Shrubs
# LC6: Herbaceous vegetation
# LC7: Cultivated and managed vegetation
# LC8: Regularly flooded shrub/herbaceous vegetation
# LC9: Urban/built-up
# LC10: Snow/ice
# LC11: Barren lands/sparse vegetation
# LC12: Open water
#1,5,79,11
###


