####################################  INTERPOLATION OF TEMPERATURES  #######################################
############################  Script for manuscript analyses,tables and figures #######################################
#This script uses the worklfow code applied to the Oregon case study. Daily methods (GAM,GWR, Kriging) are tested with
#different covariates using two baselines. Accuracy methods are added in the the function script to evaluate results.
#Figures, tables and data for the contribution of covariate paper are also produced in the script.
#AUTHOR: Benoit Parmentier                                                                      
#MODIFIED ON: 07/09/2014            
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
library(colorRamps)
#Additional libraries not used in workflow
library(pgirmess)                            # Krusall Wallis test with mulitple options, Kruskalmc {pgirmess}  
library(ncf)
#### FUNCTION USED IN SCRIPT

function_analyses_paper <-"contribution_of_covariates_paper_interpolation_functions_05212014.R"

##############################
#### Parameters and constants  

script_path<-"/home/parmentier/Data/IPLANT_project/env_layers_scripts/" #path to script
source(file.path(script_path,function_analyses_paper)) #source all functions used in this script.

in_dir1 <-"/home/parmentier/Data/IPLANT_project/Oregon_interpolation/Oregon_03142013/output_data_365d_gam_day_lst_comb3_08132013"
in_dir2 <-"/home/parmentier/Data/IPLANT_project/Oregon_interpolation/Oregon_03142013/output_data_365d_gam_day_lst_comb4_08152013"
#kriging results:
in_dir3 <-"/home/parmentier/Data/IPLANT_project/Oregon_interpolation/Oregon_03142013/output_data_365d_kriging_day_lst_comb3_07112013"
#gwr results:
in_dir4 <-"/home/parmentier/Data/IPLANT_project/Oregon_interpolation/Oregon_03142013/output_data_365d_gwr_day_lst_comb3_part1_07122013"
#multisampling results (gam)
#in_dir5<- "/home/parmentier/Data/IPLANT_project/Oregon_interpolation/Oregon_03142013/output_data_365d_gam_daily_mults10_lst_comb3_08082013"
#in_dir6<- "/home/parmentier/Data/IPLANT_project/Oregon_interpolation/Oregon_03142013/output_data_365d_kriging_daily_mults10_lst_comb3_08062013"
#in_dir7<- "/home/parmentier/Data/IPLANT_project/Oregon_interpolation/Oregon_03142013/output_data_365d_gwr_daily_mults10_lst_comb3_08072013"
#Hold-out every two days over 365 days
in_dir5 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gam_daily_mults1_lst_comb3_10122013"
in_dir6 <- "/data/project/layers/commons/Oregon_interpolation/output_data_365d_kriging_daily_mults1_lst_comb3_10112013"
in_dir7 <-"/data/project/layers/commons/Oregon_interpolation/output_data_365d_gwr_daily_mults1_lst_comb3_10132013"


infile_reg_outline <- "/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/OR83M_state_outline.shp"  #input region outline defined by polygon: Oregon
met_stations_outfiles_obj_file<-"/data/project/layers/commons/data_workflow/output_data_365d_gam_fus_lst_test_run_07172013/met_stations_outfiles_obj_gam_fusion__365d_gam_fus_lst_test_run_07172013.RData"
CRS_locs_WGS84<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0") #Station coords WGS84
y_var_name <- "dailyTmax"
out_prefix<-"_07032014"
out_dir<-"/home/parmentier/Data/IPLANT_project/paper_contribution_covar_analyses_tables_fig"
setwd(out_dir)

#out_dir<-"/home/parmentier/Data/IPLANT_project/paper_multitime_scale__analyses_tables"
create_out_dir_param = TRUE

#Create output directory

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_prefix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

#method_interpolation <- "gam_daily"
covar_obj_file_1 <- "covar_obj__365d_gam_day_lst_comb3_08132013.RData"
met_obj_file_1 <- "met_stations_outfiles_obj_gam_daily__365d_gam_day_lst_comb3_08132013.RData"
#met_stations_outfiles_obj_gam_daily__365d_gam_day_lst_comb3_08132013.RData

#raster_prediciton object for baseline 1 () s(lat,lon) + s(elev)) and baseline 2 (slat,lon))
raster_obj_file_1 <- "raster_prediction_obj_gam_daily_dailyTmax_365d_gam_day_lst_comb3_08132013.RData" 
raster_obj_file_2 <- "raster_prediction_obj_gam_daily_dailyTmax_365d_gam_day_lst_comb4_08152013.RData"

raster_obj_file_3 <- "raster_prediction_obj_kriging_daily_dailyTmax_365d_kriging_day_lst_comb3_07112013.RData"
raster_obj_file_4 <- "raster_prediction_obj_gwr_daily_dailyTmax_365d_gwr_day_lst_comb3_part1_07122013.RData"
#multisampling using baseline lat,lon + elev
#raster_obj_file_5 <- "raster_prediction_obj_gam_daily_dailyTmax_365d_gam_daily_mults10_lst_comb3_08082013.RData"
#raster_obj_file_6 <- "raster_prediction_obj_kriging_daily_dailyTmax_365d_kriging_daily_mults10_lst_comb3_08062013.RData"
#raster_obj_file_7 <- "raster_prediction_obj_gwr_daily_dailyTmax_365d_gwr_daily_mults10_lst_comb3_08072013.RData"

raster_obj_file_5 <- "raster_prediction_obj_gam_daily_dailyTmax_365d_gam_daily_mults1_lst_comb3_10122013.RData"
raster_obj_file_6 <- "raster_prediction_obj_kriging_daily_dailyTmax_365d_kriging_daily_mults1_lst_comb3_10112013.RData"
raster_obj_file_7 <- "raster_prediction_obj_gwr_daily_dailyTmax_365d_gwr_daily_mults1_lst_comb3_10132013.RData"

#Load objects containing training, testing, models objects 

met_stations_obj <- load_obj(met_stations_outfiles_obj_file)
covar_obj <-load_obj(file.path(in_dir1,covar_obj_file_1)) #Reading covariates object for GAM daily method
infile_covariates <- covar_obj$infile_covariates
infile_reg_outline <- covar_obj$infile_reg_outline
covar_names<- covar_obj$covar_names
#####
s_raster <- brick(infile_covariates)
names(s_raster)<-covar_names
ref_rast_name <-"/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/mean_day244_rescaled.rst"  #local raster name defining resolution, exent: oregon
ref_rast_d001 <-"/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/mean_day001_rescaled.rst"

raster_prediction_obj_1 <-load_obj(file.path(in_dir1,raster_obj_file_1)) #comb3 (baseline 2)
raster_prediction_obj_2 <-load_obj(file.path(in_dir2,raster_obj_file_2)) #comb4 (baseline 1)
raster_prediction_obj_3 <-load_obj(file.path(in_dir3,raster_obj_file_3)) #comb3/mod1 baseline 2, kriging
raster_prediction_obj_4 <-load_obj(file.path(in_dir4,raster_obj_file_4)) #comb3/mod1 baseline 2, gwr
raster_prediction_obj_5 <-load_obj(file.path(in_dir5,raster_obj_file_5)) #gam daily multisampling 10 to 70%
raster_prediction_obj_6 <-load_obj(file.path(in_dir6,raster_obj_file_6)) #kriging daily multisampling 10 to 70% 
raster_prediction_obj_7 <-load_obj(file.path(in_dir7,raster_obj_file_7)) #gwr daily multisampling 10 to 70%


############### BEGIN SCRIPT #################

############
##### 1) Generate: Table 3. Contribution of covariates using validation accuracy metrics
## This is a table of accuracy compared to baseline by difference 

#Check input covariates and model formula:
raster_prediction_obj_1$method_mod_obj[[1]]$formulas #models run for baseline 2
raster_prediction_obj_2$method_mod_obj[[1]]$formulas #models run for baseline 1

summary_metrics_v1<-raster_prediction_obj_1$summary_metrics_v
summary_metrics_v2<-raster_prediction_obj_2$summary_metrics_v
tb1 <-raster_prediction_obj_1$tb_diagnostic_v #365 days accuracy table for baseline 2
tb2 <-raster_prediction_obj_2$tb_diagnostic_v #365 days accuracy table for baseline 1

table_data1 <-summary_metrics_v1$avg[,c("mae","rmse","me","r")] #select relevant columns from data.frame
table_data2 <-summary_metrics_v2$avg[,c("mae","rmse","me","r")]

###Table 3a, baseline 1: s(lat,lon) 
#Change here !!! need  to reorder rows based on  mod first
model_col<-c("Baseline1","Elevation","Northing","Easting","LST","DISTOC","Forest","CANHEIGHT","LST*Forest","LST*CANHEIGHT") # 
names_table_col<-c("ΔMAE","ΔRMSE","ΔME","Δr","Model")
df3a<- as.data.frame(sapply(table_data2,FUN=function(x) x-x[1]))
df3a<- round(df3a,digit=3) #roundto three digits teh differences
df3a$Model <-model_col
names(df3a)<- names_table_col
print(df3a) #show resulting table

###Table 3b, baseline 2: s(lat,lon) + s(elev)

model_col<-c("Baseline2","Northness","Eastness","LST","DISTOC","Forest","CANHEIGHT","LST*Forest","LST*CANHEIGHT")

df3b <- as.data.frame(sapply(table_data1,FUN=function(x) x-x[1])) #difference between baseline (line 1) and other models
df3b <- round(df3b,digit=3) #roundto three digits the differences
df3b$Model <- model_col
names(df3b)<- names_table_col
print(df3b) #Part b of table 3

sd2_v <-calc_stat_from_raster_prediction_obj(raster_prediction_obj_2,"sd",training=FALSE)

#Testing siginificance between models

mod_compk1 <-kruskal.test(tb1$rmse ~ as.factor(tb1$pred_mod)) #Kruskal Wallis test
mod_compk2 <-kruskal.test(tb2$rmse ~ as.factor(tb2$pred_mod))
print(mod_compk1) #not significant
print(mod_compk2) #not significant

#Multiple Kruskal Wallis
mod_compk1 <-kruskalmc(tb1$rmse ~ as.factor(tb1$pred_mod))
mod_compk2 <-kruskalmc(tb2$rmse ~ as.factor(tb2$pred_mod))

print(mod_compk1)
print(mod_compk2)

#Now write out table 3

file_name<-paste("table3a_baseline1_paper","_",out_prefix,".txt",sep="")
write.table(df3a,file=file_name,sep=",")

file_name<-paste("table3b_baseline2_paper","_",out_prefix,".txt",sep="")
write.table(df3b,file=file_name,sep=",")

############
##### 2) Generate: Table 4. Comparison between interpolation methods using validation accuracy metrics
## This is a table of accuracy metric for the optimal model (baseline2) as identified in the previous step 

##Table 4: Interpolation methods comparison

#get sd for kriging, gam and gwr
#tb1 <-raster_prediction_obj_1$tb_diagnostic_v  HGAM baseline 1, loaded
#tb2 <-raster_prediction_obj_2$tb_diagnostic_v  #GAM baseline 2, loaded
tb3 <-raster_prediction_obj_3$tb_diagnostic_v  #Kriging methods
tb4 <-raster_prediction_obj_4$tb_diagnostic_v  #GWR methods

names_mod<-c("mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8","mod9","mod10")

#Calculate standard deviation for each metric
sd1 <- calc_stat_from_raster_prediction_obj(raster_prediction_obj_1,"sd") # see function script
sd2 <- calc_stat_from_raster_prediction_obj(raster_prediction_obj_2,"sd") # standard deviation for baseline 2
sd3 <- calc_stat_from_raster_prediction_obj(raster_prediction_obj_3,"sd") # kriging
sd4 <- calc_stat_from_raster_prediction_obj(raster_prediction_obj_4,"sd") #gwr

#Combined sd in one table for mod1 (baseline 2)
table_sd <- do.call(rbind,list(sd1[1,],sd3[1,],sd4[1,])) #table containing the sd for the three mdethods for baseline 2
table_sd <- round(table_sd[,-1],digit=3) #round to three digits the differences
table_sd <- table_sd[,c("mae","rmse","me","r")] #column 5 contains m50, remove it

summary_metrics_v3<-raster_prediction_obj_3$summary_metrics_v  #Kriging methods
summary_metrics_v4<-raster_prediction_obj_4$summary_metrics_v  # GWR method

table_data3 <-summary_metrics_v3$avg[1,c("mae","rmse","me","r")] #first line mod1 (baseline)
table_data4 <-summary_metrics_v4$avg[1,c("mae","rmse","me","r")] #first line mod1 (baseline)
table_data1 <- table_data1[1,]

table_ac <-do.call(rbind, list(table_data1,table_data3,table_data4))
table_ac <- round(table_ac,digit=3) #roundto three digits teh differences

#combined tables with accuracy metrics and their standard deviations
table4_paper <-table_combined_symbol(table_ac,table_sd,"±")
#lapply(lapply(table_ac,FUN=paste,table_sd,FUN=paste,sep="±"),FUN=paste)

model_col<-c("GAM","Kriging","GWR")
names_table_col<-c("MAE","RMSE","ME","R","Model")

table4_paper$Model <-model_col
names(table4_paper)<- names_table_col

file_name<-paste("table4_compariaons_interpolation_methods_avg_paper","_",out_prefix,".txt",sep="")
write.table(as.data.frame(table4_paper),file=file_name,sep=",")

####################################
####### Now create figures #############

#figure 1: study area
#figure 2: methodological worklfow (generated outside R)
#figure 3: LST daily and monthly climatology
#figure 4: MAE/RMSE and distance to closest fitting station.
#Figure 5. RMSE and MAE, mulitisampling and hold out for FSS and GAM.
#Figure 6. Overtraining tendency
#Figure 7a: Spatial pattern of prediction for one day
#figure 7b: Spatial autocorrelation profile : Moran's vs lag distance
#Figure 8: Figure 8. (a) Monthly MAE averages for the three interpolation methods: GAM, GWR and Kriging.(b) Monthly MAE boxplot for GAM.
#Figure 9: difference image

### Figure 1: Oregon study area
#3 parameters from output
#infile_monthly<-list_outfiles$monthly_covar_ghcn_data #outile4 from database_covar script
#infile_daily<-list_outfiles$daily_covar_ghcn_data  #outfile3 from database_covar script
#infile_locs<- list_outfiles$loc_stations_ghcn #outfile2? from database covar script
#
#
#ghcn_day_dat <- readOGR(dsn=dirname(met_stations_obj$daily_covar_ghcn_data),
#                    sub(".shp","",basename(met_stations_obj$daily_covar_ghcn_data)))
#ghcn_dayq_dat <- readOGR(dsn=dirname(met_stations_obj$daily_query_ghcn_data),
#                        sub(".shp","",basename(met_stations_obj$daily_query_ghcn_data)))


ghcn_dat <- readOGR(dsn=dirname(met_stations_obj$monthly_covar_ghcn_data),
        sub(".shp","",basename(met_stations_obj$monthly_covar_ghcn_data)))
ghcn_dat_WGS84 <-spTransform(ghcn_dat,CRS_locs_WGS84)         # Project from WGS84 to new coord. system


interp_area <- readOGR(dsn=dirname(infile_reg_outline),sub(".shp","",basename(infile_reg_outline)))
interp_area_WGS84 <-spTransform(interp_area,CRS_locs_WGS84)         # Project from WGS84 to new coord. system

usa_map <- getData('GADM', country='USA', level=1) #Get US map
usa_map_2 <- usa_map[usa_map$NAME_1!="Alaska",] #remove Alaska
usa_map_2 <- usa_map_2[usa_map_2$NAME_1!="Hawaii",] #remove Hawai 
usa_map_OR <- usa_map_2[usa_map_2$NAME_1=="Oregon",] #get OR

elev <- subset(s_raster,"elev_s")
elev_WGS84 <- projectRaster(from=elev,crs=CRS_locs_WGS84,method="ngb")

#set up the output file to plot
res_pix<-960
col_mfrow<-1
row_mfrow<-1
png(filename=paste("Figure1_contribution_covariates_study_area_",out_prefix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
par(mfrow=c(1,1))

#plot(elev_WGS84,cex.axis=2.1,cex.z=2.1)
plot(elev_WGS84,cex.axis=1.4,axis.args=list(at=seq(0,3500,500),
               labels=seq(0, 3500, 500), 
               cex.axis=1.4))
plot(interp_area_WGS84,add=T)
plot(ghcn_dat_WGS84,add=T)
title_text <-c("Elevation (m) and meteorological"," stations in Oregon")
legend("topleft",legend=title_text,cex=2.1,bty="n")

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

#########################
### Figure 2:  Method comparison workflow 

# Workflow figure is not generated in R

##########################
### Figure 3: LST averaging: daily mean compared to monthly mean

#interp_area <- readOGR(dsn=dirname(infile_reg_outline),sub(".shp","",basename(infile_reg_outline)))

lst_md <- raster(ref_rast_name)
projection(lst_md)<-projection(s_raster)
lst_mm_09<-subset(s_raster,"mm_09")

lst_md<-raster(ref_rast_d001)
lst_md<- lst_md - 273.16
lst_mm_01<-subset(s_raster,"mm_01")

#no_brks <- length(seq(min_val,max_val,by=0.1))-1
#temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
#temp.colors <- colorRampPalette(c('blue', 'lightgoldenrodyellow', 'red'))
#temp.colors <- matlab.like(no_brks)
temp.colors <- colorRampPalette(c('blue', 'khaki', 'red'))

png(filename=paste("Figure_3_paper_Comparison_daily_monthly_mean_lst",out_prefix,".png",sep=""),width=960,height=480)
par(mfrow=c(1,2))
plot(lst_md,col=temp.colors(25),axes=F) #use axes=F to remove lat and lon or coordinates
plot(interp_area,add=TRUE)
title("Mean for January 1")
plot(lst_mm_01,col=temp.colors(25),axes=F)
plot(interp_area,add=TRUE)
title("Mean for month of January")
dev.off()

path_in_tmp <- "/home/parmentier/Data/IPLANT_project/region_outlines_ref_files"
interp_area_tmp <- readOGR(dsn=path_in_tmp,"OR83M_state_outline")
proj4string(interp_area_tmp ) <- proj4string(interp_area)
projection(lst_md) <- proj4string(interp_area)
r_region_ref <- rasterize(x=interp_area_tmp,y=lst_md,field="State_ID_s",fun=max)

################################################
################### Figure 4. MAE/RMSE and distance to closest fitting station.

#Analysis accuracy in term of distance to closest station
#Assign model's names

names_mod <- paste("res_mod",1:10,sep="")
names(raster_prediction_obj_1$validation_mod_obj[[1]])
limit_val<-seq(0,150, by=10)

#Call function to extract residuals in term of distance to closest fitting station and summary statistics
l1 <- distance_to_closest_fitting_station(raster_prediction_obj_1,names_mod,dist_classes=limit_val) #GAM method
l3 <- distance_to_closest_fitting_station(raster_prediction_obj_3,names_mod,dist_classes=limit_val) #Kriging method
l4 <- distance_to_closest_fitting_station(raster_prediction_obj_4,names_mod,dist_classes=limit_val) #GWR method

l1$mae_tb #contains

list_dist_obj <-list(l1,l3,l4)
head(list_dist_obj[[1]]$dstspat_dat)
#should add method_interp in data.frame

#Prepare parameters/arguments for plotting

col_t <- c("red","blue","black")
pch_t <- 1:length(col_t)
legend_text <- c("GAM","Kriging","GWR")
mod_name <- c("res_mod1","res_mod1","res_mod1")#selected models
x_tick_labels <- limit_val<-seq(5,125, by=10)
metric_name <-"rmse_tb"
title_plot <- "RMSE and distance to closest station for baseline 2"
#y_lab_text <- "RMSE (C)"
#y_lab_text <- expression(paste("RMSE", ("^o","C")"), sep="")) 
y_lab_text <- "RMSE (°C)"
#quick test
add_CI <- c(TRUE,TRUE,TRUE)

list_param_plot<-list(list_dist_obj,col_t,pch_t,legend_text,mod_name,x_tick_labels,metric_name,title_plot,y_lab_text,add_CI)
names(list_param_plot)<-c("list_dist_obj","col_t","pch_t","legend_text","mod_name","x_tick_labels","metric_name","title_plot","y_lab_text","add_CI")
plot_dst_MAE(list_param_plot)

metric_name <-"mae_tb"
title_plot <- "MAE and distance to closest fitting station"
y_lab_text <- "MAE (°C)"
add_CI <- c(TRUE,TRUE,TRUE)
#Now set up plotting device
res_pix<-480
col_mfrow<-2
row_mfrow<-1
png_file_name<- paste("Figure_3_accuracy_and_distance_to_closest_fitting_station_",out_prefix,".png", sep="")
png(filename=file.path(out_dir,png_file_name),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
par(mfrow=c(row_mfrow,col_mfrow))

#Figure 3a
list_param_plot<-list(list_dist_obj,col_t,pch_t,legend_text,mod_name,x_tick_labels,metric_name,
                      title_plot,y_lab_text,add_CI)
names(list_param_plot)<-c("list_dist_obj","col_t","pch_t","legend_text","mod_name","x_tick_labels","metric_name",
                          "title_plot","y_lab_text","add_CI")
#undebug(plot_dst_MAE)
plot_dst_MAE(list_param_plot)
title(xlab="Distance to closest fitting station (km)")

#Figure 3b: histogram
barplot(l1$n_tb$res_mod1,names.arg=limit_val,
        ylab="Number of observations",
        xlab="Distance to closest fitting station (km)")
title("Number of observation in term of distance bins")
box()
dev.off()

####################################################
#########Figure 4. RMSE and MAE, mulitisampling and hold out for single time scale methods.

#Using baseline 2: lat,lon and elev

#Use run of 7 hold out proportions, 10 to 70% with 10 random samples and 12 dates...
#Use gam_day method
#Use comb3 i.e. using baseline s(lat,lon)+s(elev)

#names_mod<-c("mod1","mod2","mod3","mod4","mod5","mod6","mod7","mod8","mod9")
names_mod<-c("mod1")

#debug(calc_stat_prop_tb)
prop_obj_gam<-calc_stat_prop_tb(names_mod,raster_prediction_obj_5)
prop_obj_kriging<-calc_stat_prop_tb(names_mod,raster_prediction_obj_6)
prop_obj_gwr<-calc_stat_prop_tb(names_mod,raster_prediction_obj_7)

list_prop_obj<-list(prop_obj_gam,prop_obj_kriging,prop_obj_gwr)
#Testing siginificance between models

mod_compk1 <-kruskal.test(prop_obj_gwr$tb$rmse ~ as.factor(prop_obj_gwr$tb$prop)) #Kruskal Wallis test
mod_compk2 <-kruskal.test(tb2$rmse ~ as.factor(tb2$pred_mod))
print(mod_compk1) #not significant
print(mod_compk2) #not significant

#Multiple Kruskal Wallis
mod_compk1 <-kruskalmc(prop_obj_gwr$tb$rmse ~ as.factor(prop_obj_gwr$tb$prop))
mod_compk2 <-kruskalmc(tb2$rmse ~ as.factor(tb2$pred_mod))

print(mod_compk1)
print(mod_compk2)

prop_tb <- rbind(prop_obj_gam$tb,prop_obj_kriging$tb,prop_obj_gwr$tb)

mod_compk1 <-kruskal.test(prop_tb$rmse ~ as.factor(prop_tb$method_interp)+as.factor(prop_tb$prop)) #Kruskal Wallis test
mod_prop <-lm(prop_tb$rmse ~ as.factor(prop_tb$method_interp)+as.factor(prop_tb$prop)) #This is significant!!

## plot setting for figure 4

col_t<-c("red","blue","black")
pch_t<- 1:length(col_t)
legend_text <- c("GAM","Kriging","GWR")
mod_name<-c("mod1","mod1","mod1")#selected models
add_CI <- c(TRUE,TRUE,TRUE)
CI_bar <- c(TRUE,TRUE,TRUE)
#add_CI <- c(TRUE,FALSE,FALSE)

##### plot figure 4 for paper
####

res_pix<-480
col_mfrow<-2
row_mfrow<-1
png_file_name<- paste("Figure_4_proportion_of_holdout_and_accuracy_",out_prefix,".png", sep="")
png(filename=file.path(out_dir,png_file_name),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
par(mfrow=c(row_mfrow,col_mfrow))
#par(mfrow=c(1,1))
metric_name<-"mae"
list_param_plot<-list(list_prop_obj,col_t,pch_t,legend_text,mod_name,metric_name,add_CI,CI_bar)
names(list_param_plot)<-c("list_prop_obj","col_t","pch_t","legend_text","mod_name","metric_name","add_CI","CI_bar")

#debug(plot_prop_metrics)
plot_prop_metrics(list_param_plot)
title(main="MAE for hold out and methods",
      xlab="Hold out validation/testing proportion",
      ylab="MAE (°C)")

#now figure 4b
metric_name<-"rmse"
list_param_plot<-list(list_prop_obj,col_t,pch_t,legend_text,mod_name,metric_name,add_CI,CI_bar)
names(list_param_plot)<-c("list_prop_obj","col_t","pch_t","legend_text","mod_name","metric_name","add_CI","CI_bar")
plot_prop_metrics(list_param_plot)
title(main="RMSE for hold out and methods",
      xlab="Hold out validation/testing proportion",
      ylab="RMSE (°C)")

dev.off()

####################################################
#########Figure 5. Overtraining tendency

#read in relevant data:
## Calculate average difference for RMSE for all three methods
#read in relevant data:

tb1_s<-extract_from_list_obj(raster_prediction_obj_1$validation_mod_obj,"metrics_s")
rownames(tb1_s)<-NULL #remove row names
tb1_s$method_interp <- "gam_daily" #add type of interpolation...out_prefix too??

tb3_s<-extract_from_list_obj(raster_prediction_obj_3$validation_mod_obj,"metrics_s")
rownames(tb1_s)<-NULL #remove row names
tb3_s$method_interp <- "kriging_daily" #add type of interpolation...out_prefix too??

tb4_s<-extract_from_list_obj(raster_prediction_obj_4$validation_mod_obj,"metrics_s")
rownames(tb4_s)<-NULL #remove row names
tb4_s$method_interp <- "gwr_daily" #add type of interpolation...out_prefix too??

tb5_s<-extract_from_list_obj(raster_prediction_obj_5$validation_mod_obj,"metrics_s")
rownames(tb5_s)<-NULL #remove row names
tb5_s$method_interp <- "gam_daily" #add type of interpolation...out_prefix too??

tb6_s<-extract_from_list_obj(raster_prediction_obj_6$validation_mod_obj,"metrics_s")
rownames(tb6_s)<-NULL #remove row names
tb6_s$method_interp <- "kriging_daily" #add type of interpolation...out_prefix too??

tb7_s<-extract_from_list_obj(raster_prediction_obj_7$validation_mod_obj,"metrics_s")
rownames(tb7_s)<-NULL #remove row names
tb7_s$method_interp <- "gwr_daily" #add type of interpolation...out_prefix too??

#tb1_s <-raster_prediction_obj_1$tb_diagnostic_s  #gam dailycontains the accuracy metrics for each run...
#tb3_s <-raster_prediction_obj_3$tb_diagnostic_s  #Kriging daily methods
#tb4_s <-raster_prediction_obj_4$tb_diagnostic_s  #gwr daily methods

tb1 <-raster_prediction_obj_1$tb_diagnostic_v  #gam dailycontains the accuracy metrics for each run...
tb3 <-raster_prediction_obj_3$tb_diagnostic_v  #Kriging daily methods
tb4 <-raster_prediction_obj_4$tb_diagnostic_v  #gwr daily methods

tb5 <-raster_prediction_obj_5$tb_diagnostic_v  #gam dailycontains the accuracy metrics for each run...
tb6 <-raster_prediction_obj_6$tb_diagnostic_v  #Kriging daily methods
tb7 <-raster_prediction_obj_7$tb_diagnostic_v  #gwr daily methods

#Calculate difference in MAE and RMSE for training and testing data: call diff_df function
diff_tb1 <-diff_df(tb1_s[tb1_s$pred_mod=="mod1",],tb1[tb1$pred_mod=="mod1",],c("mae","rmse")) #gam select differences for mod1
diff_tb3 <-diff_df(tb3_s[tb3_s$pred_mod=="mod1",],tb3[tb3$pred_mod=="mod1",],c("mae","rmse")) #kriging
diff_tb4 <-diff_df(tb4_s[tb4_s$pred_mod=="mod1",],tb4[tb4$pred_mod=="mod1",],c("mae","rmse")) #gwr

diff_tb5 <-diff_df(tb5_s[tb5_s$pred_mod=="mod1",],tb5[tb5$pred_mod=="mod1",],c("mae","rmse")) #gam select differences for mod1
diff_tb6 <-diff_df(tb6_s[tb6_s$pred_mod=="mod1",],tb6[tb6$pred_mod=="mod1",],c("mae","rmse")) #kriging
diff_tb7 <-diff_df(tb7_s[tb7_s$pred_mod=="mod1",],tb7[tb7$pred_mod=="mod1",],c("mae","rmse")) #gwr


diff_mae_data <-data.frame(gam=diff_tb1$mae,kriging=diff_tb3$mae,gwr=diff_tb4$mae)
diff_rmse_data <-data.frame(gam=diff_tb1$rmse,kriging=diff_tb3$rmse,gwr=diff_tb4$rmse)


diff_mae_data_mult  <-data.frame(gam=diff_tb5$mae,kriging=diff_tb6$mae,gwr=diff_tb7$mae)
diff_rmse_data_mult <-data.frame(gam=diff_tb5$rmse,kriging=diff_tb6$rmse,gwr=diff_tb7$rmse)

diff_mae_data_mult$prop <- tb5$prop

boxplot(diff_mae_data_mult$gam~diff_mae_data_mult$prop)
boxplot(diff_mae_data_mult$kriging~diff_mae_data_mult$prop)

plot(diff_mae_data_mult$gam~diff_mae_data_mult$prop,type="p")



#Test the plot
boxplot(diff_mae_data)
boxplot(diff_rmse_data) #plot differences in training and testing accuracies for three methods
title(main="Training and testing RMSE for hold out and interpolation methods",
      xlab="Interpolation method",
      ylab="RMSE (°C)")

boxplot(diff_mae_data_mult)
boxplot(diff_rmse_data_mult) #plot differences in training and testing accuracies for three methods


boxplot(diff_mae_data_mult)

tb5 <-raster_prediction_obj_5$tb_diagnostic_v  #gam dailycontains the accuracy metrics for each run...
tb6 <-raster_prediction_obj_6$tb_diagnostic_v  #Kriging daily methods
tb7 <-raster_prediction_obj_7$tb_diagnostic_v  #gwr daily methods

prop_obj_gam_s<-calc_stat_prop_tb(names_mod,raster_prediction_obj_5,testing=FALSE) #training accuracy with hold out proportion
prop_obj_kriging_s<-calc_stat_prop_tb(names_mod,raster_prediction_obj_6,testing=FALSE)
prop_obj_gwr_s<-calc_stat_prop_tb(names_mod,raster_prediction_obj_7,testing=FALSE)

plot(prop_obj_gam_s$avg_tb$rmse ~ prop_obj_gam_s$avg_tb$prop, type="b",)

###############
#### plot figure 5
#set up the output file to plot
res_pix<-480
col_mfrow<-2
row_mfrow<-1
png_file_name<- paste("Figure_5_overtraining_tendency_and_holdout_proportion_",out_prefix,".png", sep="")
png(filename=file.path(out_dir,png_file_name),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
par(mfrow=c(row_mfrow,col_mfrow))

col_t<-c("red","blue","black")
pch_t<- 1:length(col_t)
##Make this a function???
y_range<-range(prop_obj_kriging$avg_tb$rmse,prop_obj_kriging_s$avg_tb$rmse)
#y_range<-range(prop_obj_gam$avg_tb$rmse,prop_obj_gam_s$avg_tb$rmse)
#plot(prop_obj_gam$avg_tb$rmse ~ prop_obj_gam$avg_tb$prop, ylab="",xlab="",type="b",col=c("red"),pch=pch_t[1],ylim=y_range,lty=2)
#lines(prop_obj_gam_s$avg_tb$rmse ~ prop_obj_gam_s$avg_tb$prop, ylab="",xlab="",type="b",pch=pch_t[1],ylim=y_range,col=c("red"))
#lines(prop_obj_gwr_s$avg_tb$rmse ~ prop_obj_gwr_s$avg_tb$prop, ylab="",xlab="",type="b",ylim=y_range,pch=pch_t[3],col=c("black"))
#lines(prop_obj_gwr$avg_tb$rmse ~ prop_obj_gam$avg_tb$prop, ylab="",xlab="",type="b",ylim=y_range,pch=pch_t[3],,col=c("black"),lty=2)
#lines(prop_obj_kriging$avg_tb$rmse ~ prop_obj_kriging$avg_tb$prop,ylab="",xlab="", type="b",ylim=y_range,pch=pch_t[2],,col=c("blue"),lty=2)
#lines(prop_obj_kriging_s$avg_tb$rmse ~ prop_obj_kriging_s$avg_tb$prop,ylab="",xlab="",type="b",ylim=y_range,pch=pch_t[2],col=c("blue"))

plot_ac_holdout_prop<- function(l_prop,l_col_t,l_pch_t,add_CI,y_range){
  
  for(i in 1:length(l_prop)){
    if(i==1){
      plot(l_prop[[i]]$avg_tb$rmse ~ l_prop[[i]]$avg_tb$prop,ylab="",xlab="",type="b",pch=l_pch_t[i],ylim=y_range,col=l_col_t[i],lty=l_lty_t[i])
    }else{
      lines(l_prop[[i]]$avg_tb$rmse ~ l_prop[[i]]$avg_tb$prop,ylab="",xlab="",type="b",pch=l_pch_t[i],ylim=y_range,col=l_col_t[i],lty=l_lty_t[i])
    }
    if(add_CI[i]==TRUE){
      ciw   <- qt(0.975, l_prop[[i]]$n_tb$rmse) * l_prop[[i]]$sd_tb$rmse / sqrt(l_prop[[i]]$n_tb$rmse)
      plotCI(y=l_prop[[i]]$avg_tb$rmse, x=l_prop[[i]]$avg_tb$prop, uiw=ciw, col=l_col_t[i], barcol="blue", lwd=1,pch=l_pch_t[i],
             add=TRUE) 
    }
  }
}

l_prop<- list(prop_obj_gam,prop_obj_gam_s,prop_obj_gwr_s,prop_obj_gwr,prop_obj_kriging,prop_obj_kriging_s)

l_col_t <- c("red","red","black","black","blue","blue")
l_pch_t <- c(1,1,3,3,2,2)
l_lty_t <- c(2,1,1,2,2,1)
add_CI <- c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
add_CI <- c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
y_range<-c(0.5,3)

plot_ac_holdout_prop(l_prop,l_col_t,l_pch_t,add_CI,y_range)

legend_text <- c("GAM","Kriging","GWR")
#legend_text <- c("GAM","Kriging","GWR","training","testing")

#legend("topleft",legend=legend_text, 
#       cex=0.9, pch=c(pch_t,NA,NA),col=c(col_t,"white","white"),lty=c(NA,NA,NA,1,2),bty="n")
legend("topleft",legend=legend_text, 
       cex=0.9, pch=c(pch_t),col=c(col_t),lty=c(1,1,1),bty="n")
legend_text_data <-c("training","testing")
legend("top",legend=legend_text_data, 
       cex=0.9, lty=c(1,2),bty="n")

title(main="Training and testing RMSE for hold out and methods",
      xlab="Hold out validation/testing proportion",
      ylab="RMSE (°C)")


boxplot(diff_mae_data_mult[-4]) #plot differences in training and testing accuracies for three methods

title(main="Difference between training and testing MAE",
      xlab="Interpolation method",
      ylab="MAE (°C)")

dev.off()

############### STUDY TIME AND accuracy
#########Figure 6: Monthly RMSE averages for the three interpolation methods: GAM, GWR and Kriging.

mae_tmp<- data.frame(gam=tb1[tb1$pred_mod=="mod1",c("mae")],
                     kriging=tb3[tb3$pred_mod=="mod1",c("mae")],
                     gwr=tb4[tb4$pred_mod=="mod1",c("mae")])

plot(mae_tmp$gam,col=c("red"),type="b",pch=1)
lines(mae_tmp$kriging,col=c("blue"),type="b",pch=2)
lines(mae_tmp$gwr,col=c("black"),type="b",pch=2)
legend("topleft",legend=legend_text, 
       cex=1.2, pch=pch_t,col=col_t,lty=1,bty="n")

max(mae_tmp$gam)

x2<-tb1[tb1$pred_mod=="mod1",c("mae","date")]
arrange(x2,desc(mae))

#kriging=tb3[tb3$pred_mod=="mod1",c("mae")],
#                     gwr=tb4[tb4$pred_mod=="mod1",c("mae")])

##### MONTHLY AVERAGES

tb1_month<-raster_prediction_obj_1$summary_month_metrics_v[[1]] #note that this is for model1
tb3_month<-raster_prediction_obj_3$summary_month_metrics_v[[1]]
tb4_month<- raster_prediction_obj_4$summary_month_metrics_v[[1]]

y_range<-range(tb1_month$mae,tb3_month$mae,tb4_month$mae)
plot(1:12,tb1_month$mae,col=c("red"),type="b",ylim=y_range)
lines(1:12,tb3_month$mae,col=c("blue"),type="b")
lines(1:12,tb4_month$mae,col=c("black"),type="b")

add_month_tag<-function(tb){
  date<-strptime(tb$date, "%Y%m%d")   # interpolation date being processed
  month<-strftime(date, "%m")          # current month of the date being processed
}
tb1$month<-add_month_tag(tb1)
tb3$month<-add_month_tag(tb3)
tb4$month<-add_month_tag(tb4)

metric_name<-"mae"
month_data_list<-list(gam=tb1[tb1$pred_mod=="mod1",c(metric_name,"month")],
                      kriging=tb3[tb3$pred_mod=="mod1",c(metric_name,"month")],
                      gwr=tb4[tb4$pred_mod=="mod1",c(metric_name,"month")])
y_range<-range(unlist(month_data_list))


#boxplot(test[[metric_ac]]~test[[c("month")]],outline=FALSE,horizontal=FALSE,cex=0.5,
#        ylab=paste(metric_ac,"in degree C",sep=" "),,axisnames=FALSE,axes=FALSE)
#boxplot(test[[metric_ac]]~test[[c("month")]],outline=FALSE,horizontal=FALSE,cex=0.5,
#        ylab=paste(metric_ac,"in degree C",sep=" "))
#axis(1, labels = FALSE)
## Create some text labels
#labels <- month.abb # abbreviated names for each month
## Plot x axis labels at default tick marks
#text(1:length(labels), par("usr")[3] - 0.25, srt = 45, adj = 1,
#     labels = labels, xpd = TRUE)
#axis(2)
#box()

#Now plot figure 6
res_pix<-480
col_mfrow<-2
#row_mfrow<-2
row_mfrow<-1

png_file_name<- paste("Figure_6_monthly_accuracy_",out_prefix,".png", sep="")
png(filename=file.path(out_dir,png_file_name),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
par(mfrow=c(row_mfrow,col_mfrow))

col_t<-c("red","blue","black")
pch_t<- 1:length(col_t)
legend_text <- c("GAM","Kriging","GWR")

y_range<-range(tb1_month$mae,tb3_month$mae,tb4_month$mae)
xlab_tick <- unique(tb1$month)
xlab_text <-"Month"
ylab_text <- "MAE (°C)"
plot(1:12,tb1_month$mae,col=c("red"),type="b",ylim=y_range,xlab=xlab_text,ylab=ylab_text,xaxt="n")
lines(1:12,tb3_month$mae,col=c("blue"),type="b")
lines(1:12,tb4_month$mae,col=c("black"),type="b")
axis(1,at=1:12,labels=xlab_tick)
title(main="Monthly average MAE")
legend("topleft",legend=legend_text, 
       cex=0.9, pch=c(pch_t),col=c(col_t),lty=c(1,1,1),bty="n")

#Second plot
ylab_text<-"MAE (°C)"
xlab_text<-"Month"
#y_range<-range(month_data_list$gam$mae,month_data_list$kriging$mae,month_data_list$gwr$mae)
#y_range<-range(month_data_list$gam$mae)
boxplot(mae~month,data=month_data_list$gam,main="Monthly MAE boxplot", xlab=xlab_text,ylab=ylab_text,outline=FALSE)
#boxplot(mae~month,data=month_data_list$kriging,ylim=y_range,main="Kriging",ylab=ylab_text,outline=FALSE)
#boxplot(mae~month,data=month_data_list$gwr,ylim=y_range,main="GWR",ylab=ylab_text,outline=FALSE)

dev.off()

#Now generate table 5

test<-boxplot(mae~month,data=month_data_list$gam,main="GAM",ylab=ylab_text,outline=FALSE)

length(tb1_month$mae)
names(tb1_month)

#Calculate standard deviation for each metric
sd1 <- calc_stat_month_from_raster_prediction_obj(raster_prediction_obj_1,"sd") # see function script
sd2 <- calc_stat_month_from_raster_prediction_obj(raster_prediction_obj_2,"sd") # standard deviation for baseline 2
sd3 <- calc_stat_month_from_raster_prediction_obj(raster_prediction_obj_3,"sd") # kriging
sd4 <- calc_stat_month_from_raster_prediction_obj(raster_prediction_obj_4,"sd") #gwr

avg_v1 <- calc_stat_month_from_raster_prediction_obj(raster_prediction_obj_1,"mean") # see function script
avg_v2 <- calc_stat_month_from_raster_prediction_obj(raster_prediction_obj_2,"mean") # standard deviation for baseline 2
avg_v3 <- calc_stat_month_from_raster_prediction_obj(raster_prediction_obj_3,"mean") # kriging
avg_v4 <- calc_stat_month_from_raster_prediction_obj(raster_prediction_obj_4,"mean") #gwr

#Combined sd in one table for mod1 (baseline 2)
table_sd <- do.call(cbind,list(gam=sd1[sd1$pred_mod=="mod1",c("rmse")],
                               kriging=sd3[sd3$pred_mod=="mod1",c("rmse")],
                               gwr=sd4[sd4$pred_mod=="mod1",c("rmse")])) #table containing the sd for the three mdethods for baseline 2
table_sd <- as.data.frame(round(table_sd,digit=3)) #round to three digits the differences

#Combined mean in one table for mod1 (baseline 2)
table_avg <- do.call(cbind,list(gam=avg_v1[avg_v1$pred_mod=="mod1",c("rmse")],
                               kriging=avg_v3[sd3$pred_mod=="mod1",c("rmse")],
                               gwr=avg_v4[avg_v4$pred_mod=="mod1",c("rmse")])) #table containing the sd for the three mdethods for baseline 2
table_avg <- as.data.frame(round(table_avg,digit=3)) #round to three digits the differences

#combined tables with accuracy metrics and their standard deviations
table5_paper <-table_combined_symbol(table_avg,table_sd,"±")
table5_paper$month <- month.abb

file_name<-paste("table5_comparisons_monthly_averages_interpolation_methods_paper","_",out_prefix,".txt",sep="")
write.table(as.data.frame(table5_paper),file=file_name,sep=",")

####### FIGURE 7: Spatial pattern ######

y_var_name <-"dailyTmax"
index<-244 #index corresponding to Sept 1

lf1 <- raster_prediction_obj_1$method_mod_obj[[index]][[y_var_name]] #select relevant raster images for the given dates
lf3 <- raster_prediction_obj_3$method_mod_obj[[index]][[y_var_name]]
lf4 <- raster_prediction_obj_4$method_mod_obj[[index]][[y_var_name]]

date_selected <- "20109101"
methods_names <-c("gam","kriging","gwr")
names_layers<-methods_names
lf <-list(lf1$mod1,lf3$mod1,lf4$mod1)
#lf <-lf[[1]]

pred_temp_s <-stack(lf)
#predictions<-mask(predictions,mask_rast)
names(pred_temp_s)<-names_layers
s.range <- c(min(minValue(pred_temp_s)), max(maxValue(pred_temp_s)))
#s.range <- s.range+c(5,-5)
col.breaks <- pretty(s.range, n=200)
lab.breaks <- pretty(s.range, n=100)
temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
max_val<-s.range[2]
min_val <-s.range[1]
#max_val<- -10
min_val <- 0
layout_m<-c(1,3) #one row two columns

#png(paste("Figure7__spatial_pattern_tmax_prediction_levelplot_",date_selected,out_prefix,".png", sep=""),
#    height=480*layout_m[1],width=480*layout_m[2])

p<-levelplot(pred_temp_s,main="Interpolated Surfaces Method Comparison", ylab=NULL,xlab=NULL,
          par.settings = list(axis.text = list(font = 2, cex = 1.3),layout=layout_m,
                              par.main.text=list(font=2,cex=2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
          names.attr=names_layers,col.regions=temp.colors,at=seq(max_val,min_val,by=0.01))
print(p)
#col.regions=temp.colors(25))
#dev.off()

## FIGURE COMPARISON OF  MODELS COVARRIATES: Figure 7...

lf2 <- raster_prediction_obj_2$method_mod_obj[[index]][[y_var_name]]
lf2 #contains the models for gam

pred_temp_s <-stack(lf2)
date_selected <- "20109101"
#names_layers <-c("mod1=s(lat,long)+s(elev)","mod4=s(lat,long)+s(LST)","diff=mod1-mod4")
names_layers <-c("mod1 = s(lat,long)","mod2 = s(lat,long)+s(elev)","mod3 = s(lat,long)+s(N_w)","mod4 = s(lat,long)+s(E_w)",
                 "mod5 = s(lat,long)+s(LST)","mod6 = s(lat,long)+s(DISTOC)","mod7 = s(lat,long)+s(FOREST)",
                 "mod8 = s(lat,long)+s(CANHGHT)","mod9 = s(lat,long)+s(LST,FOREST)","mod10 = s(lat,long)+s(LST,CANHGHT)")

#names_layers<-names(pred_temp_s)
#names(pred_temp_s)<-names_layers

s.range <- c(min(minValue(pred_temp_s)), max(maxValue(pred_temp_s)))
#s.range <- s.range+c(5,-5)
col.breaks <- pretty(s.range, n=200)
lab.breaks <- pretty(s.range, n=100)
temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
max_val<-s.range[2]
min_val <-s.range[1]
#max_val<- -10
#min_val <- 0
layout_m<-c(4,3) #one row two columns

png(paste("Figure_7_spatial_pattern_tmax_prediction_models_baseline1_gam_levelplot_",date_selected,out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

p<- levelplot(pred_temp_s,main="Interpolated Surfaces Model Comparison baseline 1", ylab=NULL,xlab=NULL,
          par.settings = list(axis.text = list(font = 2, cex = 1.3),layout=layout_m,
                              par.main.text=list(font=2,cex=2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
          names.attr=names_layers,col.regions=temp.colors,at=seq(max_val,min_val,by=0.01))
#col.regions=temp.colors(25))
print(p)
dev.off()

########################################################
#### Examining spatial correlograms for each model...use spedep but can be problematic when many datapoints
#The follwing method using Moran's I raster command works for now but does not provide CI/std dev...

r_stack <-pred_temp_s

r<- subset(r_stack,"mod1")
Moran(r) #with lag 1 and default rooks lag correlation

#generate filters for 10 lags: quick solution

list_filters<-lapply(1:10,FUN=autocor_filter_fun,f_type="queen") #generate 10 filters
#moran_list <- lapply(list_filters,FUN=Moran,x=r)

list_param_moran <- list(list_filters=list_filters,r_stack=r_stack)
#moran_r <-moran_multiple_fun(1,list_param=list_param_moran)
nlayers(r_stack) 
moran_I_df <-mclapply(1:nlayers(r_stack), list_param=list_param_moran, FUN=moran_multiple_fun,mc.preschedule=FALSE,mc.cores = 10) #This is the end bracket from mclapply(...) statement

moran_df <- do.call(cbind,moran_I_df) #bind Moran's I value 10*nlayers data.frame
moran_df$lag <-1:nrow(moran_df)

#prepare to automate the plotting of   all columns
mydata<-moran_df
dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
dd$lag <- mydata$lag 

layout_m<-c(4,3) #one row two columns

png(paste("Figure_7b_spatial_correlogram_tmax_prediction_models_baseline1_gam_levelplot_",date_selected,out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

p<-xyplot(data ~ lag | which, dd,type="b",
          par.settings = list(axis.text = list(font = 2, cex = 1.3),layout=layout_m,
                              par.main.text=list(font=2,cex=2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
          strip=strip.custom(factor.levels=names_layers),
          xlab=list(label="Spatial lag neighbor", cex=2,font=2),
          ylab=list(label="Moran's I", cex=2, font=2))
print(p)

dev.off()

### median for first neighbour distances in data
ghcn_day_dat <- readOGR(dsn=dirname(met_stations_obj$daily_covar_ghcn_data),
                        sub(".shp","",basename(met_stations_obj$daily_covar_ghcn_data)))
loc_spdf <- readOGR(dsn=dirname(met_stations_obj$loc_stations),
                        sub(".shp","",basename(met_stations_obj$loc_stations)))
loc_2010_spdf<-subset(loc_spdf,loc_spdf$STAT_ID%in%ghcn_day_dat$station) #stations in 2010
loc_2010_spdf <- spTransform(loc_2010_spdf, CRS=CRS(proj4string(ghcn_day_dat)))
loc_knn1 <- knearneigh(coordinates(loc_2010_spdf), k=1) #lag1
class(loc_knn1) #knn object
loc_nb1 <- knn2nb(loc_knn1) #nb object
dsts_nb1 <-unlist(nbdists(loc_nb1,coordinates(loc_2010_spdf))) #vector with first nearest neighbour distance for each station

median(dsts_nb1) #median first nearest neighbour distance: 20023.03 m
mean(dsts_nb1) #mean first nearest neighbour distance : 22480.91 meters
sd(dsts_nb1) #
hist(dsts_nb1)

################ #FIGURE 8: difference image

lf1 <- raster_prediction_obj_1$method_mod_obj[[index]][[y_var_name]]
lf1 #contains the models for gam

pred_temp_s <-stack(lf1$mod1,lf1$mod4)
date_selected <- "20109101"
#names_layers <-c("mod1=s(lat,long)+s(elev)","mod4=s(lat,long)+s(LST)","diff=mod1-mod4")
names_layers <-c("mod1 = s(lat,long)+s(elev)","mod4 = s(lat,long)+s(LST)")
names(pred_temp_s)<-names_layers
s.range <- c(min(minValue(pred_temp_s)), max(maxValue(pred_temp_s)))
#s.range <- s.range+c(5,-5)
col.breaks <- pretty(s.range, n=200)
lab.breaks <- pretty(s.range, n=100)
temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
max_val<-s.range[2]
min_val <-s.range[1]
#max_val<- -10
min_val <- 0
layout_m<-c(1,2) #one row two columns

png(paste("Figure_8a_spatial_pattern_tmax_prediction_models_gam_levelplot_",date_selected,out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

p<-levelplot(pred_temp_s,main="Interpolated Surfaces Model Comparison", ylab=NULL,xlab=NULL,
          par.settings = list(axis.text = list(font = 2, cex = 1.3),layout=layout_m,
                              par.main.text=list(font=2,cex=2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
          names.attr=names_layers,col.regions=temp.colors,at=seq(max_val,min_val,by=0.01))
#col.regions=temp.colors(25))
print(p)
#col.regions=temp.colors(25))
dev.off()


diff<-raster(lf1$mod1)-raster(lf1$mod4)
names_layers <- c("difference=mod1-mod4")
names(diff) <- names_layers


png(paste("Figure_8b_spatial_pattern_tmax_prediction_models_gam_levelplot_",date_selected,out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

plot(diff,col=temp.colors(100),main=names_layers)
#levelplot(diff,main="Interpolated Surfaces Model Comparison", ylab=NULL,xlab=NULL,
#          par.settings = list(axis.text = list(font = 2, cex = 1.3),layout=c(1,1),
#                              par.main.text=list(font=2,cex=2),strip.background=list(col="white")),par.strip.text=list(font=2,cex=1.5),
#          names.attr=names_layers,col.regions=temp.colors)
dev.off

###############################
########## Prepare table 6
# Now get p values and other things...
###baseline 2: s(lat,lon) + s(elev)
      
l_obj<-vector("list",length=2)
l_obj[[1]]<-raster_prediction_obj_1
l_obj[[2]]<-raster_prediction_obj_2
l_table <- vector("list",length=length(l_obj))   
for (k in 1:length(l_obj)){
  raster_prediction_obj<- l_obj[[k]]
  list_myModels <- extract_list_from_list_obj(raster_prediction_obj$method_mod_obj,"mod")
    
  list_models_info <-lapply(1:length(list_myModels),FUN=create_s_and_p_table_term_models,list_myModels)
  dates<-(extract_from_list_obj(raster_prediction_obj_1$method_mod_obj,"sampling_dat"))$date #get vector of dates
  names(list_models_info)<-dates #adding names to the list
    
  #Prepare and process p. value information regarding models: count number of times values were above a threshold.
  s.table_term_tb <-extract_from_list_obj(list_models_info,"s.table_term")
    
  threshold_val<-c(0.01,0.05,0.1)
  s.table_term_tb$p_val_rec1 <- s.table_term_tb[["p-value"]] < threshold_val[1]
  s.table_term_tb$p_val_rec2 <- s.table_term_tb[["p-value"]] < threshold_val[2]
  s.table_term_tb$p_val_rec3 <- s.table_term_tb[["p-value"]] < threshold_val[3]
  
  names_var <- c("p_val_rec1","p_val_rec2","p_val_rec3")
  t2<-melt(s.table_term_tb,
           measure=names_var, 
           id=c("mod_name","term_name"),
           na.rm=T)
  
  summary_s.table_term2 <- cast(t2,term_name+mod_name~variable,sum)
  
  #Now add AIC
  AIC_models_tb <-extract_from_list_obj(list_models_info,"AIC_models")
  AIC_models_tb
  names_var <- c("AIC")
  #id_var <- 
  t3<-melt(AIC_models_tb,
           measure=names_var, 
           id=c("mod_name","term_name"),
           na.rm=T)
  
  summary_AIC <- cast(t3,term_name+mod_name~variable,median)
  summary_AIC$AIC <- round(summary_AIC[,c("AIC")],digit=2) #roundto three digits teh differences
  
  #Now combine tables and drop duplicate columns the combined table can be modified for the paper...
  avg_s <- calc_stat_from_raster_prediction_obj(raster_prediction_obj,"mean",training=TRUE)
  avg_v <- calc_stat_from_raster_prediction_obj(raster_prediction_obj,"mean",training=FALSE)
  avg <- cbind(avg_s[,c("pred_mod","mae")],avg_v[,c("mae")])
  names(avg)<-c("pred_mod","mae_s","mae_v")
  avg[,c("mae_s","mae_v")] <- round(avg[,c("mae_s","mae_v")],digit=2)
  table <- merge(summary_AIC,avg,by.x="mod_name",by.y="pred_mod")
  
  tables_AIC_ac_p_val <-list(table,summary_s.table_term2)
  names(tables_AIC_ac_p_val) <-c("table","s.table_p_val_term")
  l_table[[k]] <- tables_AIC_ac_p_val
}

## Now prepare table
s.table_p_val_term <- l_table[[1]]$s.table_p_val_term[-c(10:26),]
s.table_p_val_term <- s.table_p_val_term[order(s.table_p_val_term$mod_name),]
#summary_s.table_term2 <- summary_s.table_term2[-c(8,10),]
table <- l_table[[1]]$table
      
table6a <- merge(s.table_p_val_term,table,by="mod_name")  
table6a <- table6a[,-match(c("term_name.y"),names(table6a))]
#model_col<-c("Baseline2","Northness","Eastness","LST","DISTOC","Forest","CANHEIGHT","LST*Forest") # removed ,"LST*CANHEIGHT")
#names_table_col<-c("DiffMAE","DiffRMSE","DiffME","Diffr","Model")
       
s.table_p_val_term <- l_table[[2]]$s.table_p_val_term[-c(11:19),]
s.table_p_val_term <- s.table_p_val_term[order(s.table_p_val_term$mod_name),]
#summary_s.table_term2 <- summary_s.table_term2[-c(8,10),]
tableb <- l_table[[2]]$table
      
table6b <- merge(s.table_p_val_term,tableb,by="mod_name")  
table6b <- table6b[,-match(c("term_name.y"),names(table6b))]
#model_col<-c("Baseline2","Northness","Eastness","LST","DISTOC","Forest","CANHEIGHT","LST*Forest") # removed ,"LST*CANHEIGHT")
#names_table_col<-c("DiffMAE","DiffRMSE","DiffME","Diffr","Model")
      
#table6b[order(table6b$mod_name),]
      
#Now write out table...

file_name<-paste("table6a_paper","_",out_prefix,".txt",sep="")
write.table(table6a,file=file_name,sep=",")

file_name<-paste("table6b_paper","_",out_prefix,".txt",sep="")
write.table(table6b,file=file_name,sep=",")

########################
### Prepare table 7: correlation matrix between covariates      

names(s_raster)

list_formulas<-raster_prediction_obj_2$method_mod_obj[[1]]$formulas
list_formulas <- lapply(list_formulas,FUN=as.formula)
covar_names_model <- unique(unlist(lapply(list_formulas,FUN=all.vars)))[-1]  
covar_names_model <- c(covar_names_model[-6],c("mm_01","mm_02","mm_03","mm_04","mm_05","mm_06",
                                               "mm_07","mm_08","mm_09","mm_10","mm_11","mm_12"))
covar_raster <- subset(s_raster,covar_names_model)

names(covar_raster) <- covar_names_model
corr_layers_covar <-layerStats(covar_raster,"pearson",na.rm=TRUE)
corr_mat <-round(corr_layers_covar[[1]], digit=2)
      
file_name<-paste("table7_paper","_",out_prefix,".txt",sep="")
write.table(corr_mat,file=file_name,sep=",")
      
#met_obj <-load_obj(file.path(in_dir1,met_obj_file_1))
#stat_loc<-readOGR(dsn=in_dir1,layer=sub(".shp","",basename(met_obj$loc_stations)))
#dim(stat_loc)      

############################
#### Figure 9: variograms for 2 dates  and model 1 ###
      
#Make a function to extract and plot variogram models...one per day and model
#plot.gstatVariogram
#Plot a sample variogram, and possibly a fitted model
#model 1 lat,lon and elev
layout_m <- c(1,2) # works if set to this?? ok set the resolution...
      
date_selected <- c("20100101","20100901")
png(paste("Figure9_paper_","_variogram_",date_selected[1],"_",date_selected[2],"_",out_prefix,".png", sep=""),
    height=960*layout_m[1],width=960*layout_m[2])
    #height=480*6,width=480*4)

#p3 <- list_plots_spt[[3]]
p1<-plot(raster_prediction_obj_3$method_mod_obj[[1]]$mod[[1]]$exp_var,raster_prediction_obj_3$method_mod_obj[[1]]$mod[[1]]$var_model,
         ylim=c(0,9),
         ylab=list(label="Semivariance",cex=1.5),
         xlab=list(label="Distance (meter)",cex=1.5),
         main=list(label="Mod1 January 1, 2010",cex=1.8))
#plot(p1)
p241<-plot(raster_prediction_obj_3$method_mod_obj[[241]]$mod[[1]]$exp_var,raster_prediction_obj_3$method_mod_obj[[241]]$mod[[1]]$var_model,
         ylim=c(0,9),
         ylab=list(label="Semivariance",cex=1.5),
         xlab=list(label="Distance (meter)",cex=1.5),
         main=list(label="Mod1 September 1, 2010",cex=1.8))
#plot(p241)

#grid.arrange(p1,p2,p3,ncol=1)
grid.arrange(p1,p241,ncol=2)

dev.off()
#Combine both plot?     + plot info on sill, nugget and range? and most frequent model selected

list_kg_var_model <- lapply(1:365,function(i){raster_prediction_obj_3$method_mod_obj[[i]]$mod[[1]]$var_model})
list_kg_var_model <- lapply(1:365,function(i){raster_prediction_obj_3$method_mod_obj[[i]]$mod[[1]]$var_model})

list_kg_var_model
   
test <- do.call(rbind,list_kg_var_model)
tb_variogram <- subset(test,model!="Nug")      
tt2 <- subset(test,model=="Nug")      
tb_variogram["Nug"] <- (tt2$psill)

add_month_tag<-function(tb){
  date<-strptime(tb$date, "%Y%m%d")   # interpolation date being processed
  month<-strftime(date, "%m")          # current month of the date being processed
}

dates<-(extract_from_list_obj(raster_prediction_obj_1$method_mod_obj,"sampling_dat"))$date #get vector of dates
tb_variogram["date"] <- dates
tb_variogram$month <- add_month_tag(tb_variogram)
#add dates and month??


#layout_m <- c(2,2) # works if set to this?? ok set the resolution...
res_pix<-480
col_mfrow<-2
#row_mfrow<-2
row_mfrow<-2

png_file_name<- paste("Figure_14_paper_variogram_statistics_",out_prefix,".png", sep="")
png(filename=file.path(out_dir,png_file_name),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
par(mfrow=c(row_mfrow,col_mfrow))
     
#date_selected <- c("20100101","20100901")
#png(paste("Figure14_paper_","_variogram_statistics","_",out_prefix,".png", sep=""),
#    height=960*layout_m[1],width=960*layout_m[2])
    #height=480*6,width=480*4)
#mfrow(2,2)
p_hist<-histogram(tb_variogram$model,col="grey",
          main=list(label="Variogram model type",cex=1.8),
          ylab=list(label="Percent of total",cex=1.5),
          xlab=list(label="Variogram model",cex=1.5))      
print(p_hist)
#hist(tt1$range)
#hist(tt1$psill)
#hist(tt1$Nug)
boxplot(tb_variogram$range~tb_variogram$month,outline=F)
boxplot(tb_variogram$Nug~tb_variogram$month,outline=F)
boxplot(tb_variogram$psill~tb_variogram$month,outline=F)
#boxplot(tb_variogram$model~tb_variogram$month,outline=F)

dev.off()

###########################################
### Figure 10: map of residuals...for a specific date...
index <- 244
names_mod <- names(raster_prediction_obj_2$method_mod_obj[[index]][[y_var_name]]) #names of models to plot
#in_dir2 <-"/home/parmentier/Data/IPLANT_project/Oregon_interpolation/Oregon_03142013/output_data_365d_gam_day_lst_comb4_08152013"
#date_selected <-  c("20100101")
date_selected <-  c("20100901")

data_v <- raster_prediction_obj_2$validation_mod_obj[[index]]$data_v #testing with residuals
#data_s<-validation_mod_obj[[index]]$data_s
date_proc<-strptime(date_selected, "%Y%m%d")   # interpolation date being processed
mo<-as.integer(strftime(date_proc, "%m"))          # current month of the date being processed
day<-as.integer(strftime(date_proc, "%d"))
year<-as.integer(strftime(date_proc, "%Y"))
datelabel=format(ISOdate(year,mo,day),"%b %d, %Y")

#height=480*6,width=480*4)
list_p <- vector("list", length(names_mod))
for (k in 1:length(names_mod)){
  model_name <- names_mod[k]
  #fig_name <- file.path(out_path,paste("Figure_residuals_map_",y_var_name,"_",model_name,"_",sampling_dat$date,"_",sampling_dat$prop,"_",
  #                             sampling_dat$run_samp,out_prefix,".png", sep=""))
    
  png(file.path(out_dir,paste("Figure_residuals_map_",y_var_name,"_",model_name,"_",date_selected,out_prefix,".png", sep="")))
  res_model_name <- paste("res",model_name,sep="_")
  elev <- subset(s_raster,"elev_s")
  #p1 <- levelplot(elev,scales = list(draw = FALSE), colorkey = FALSE,col.regions=rev(terrain.colors(255)),contour=T)
  #add legend..par.settings = GrTheme
  cx <- ((data_v[[res_model_name]])*2)
  p1 <- levelplot(elev,scales = list(draw = FALSE), colorkey = FALSE,par.settings = GrTheme)

  #p2 <- spplot(data_v,res_model_name, cex=1 * cx,main=paste("Residuals for ",res_model_name," ",datelabel,sep=""),
  #             col.regions=matlab.like(25))
  p2 <- bubble(data_v,res_model_name, main=paste("Residuals for ",res_model_name," ",datelabel,sep=""),
               col=matlab.like(25))  
  p3 <- p2 + p1 + p2 #to force legend...
  #p1 <- spplot(interp_area)
  #p3 <- p1+p2 #to force legend...
  #p + layer(sp.lines(mapaSHP, lwd=0.8, col='darkgray'))

  #p2
  print(p3)

  dev.off()
  list_p[[k]] <- p3
}

layout_m<-c(1,3) # works if set to this?? ok set the resolution...
png(paste("Figure13_paper_","_residuals_",date_selected,"_",out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

grid.arrange(list_p_mae[[1]],list_p_mae[[2]],list_p_mae[[3]],ncol=3)
      
dev.off()   

######## NOW GET A ACCURACY BY STATIONS
l_formulas<-(extract_from_list_obj(raster_prediction_obj_2$method_mod_obj,"formulas")) #get vector of dates
#                            y_var ~ s(lat,lon)
#                y_var ~ s(lat,lon) + s(elev_s)
#                   y_var ~ s(lat,lon) + s(N_w)
#                   y_var ~ s(lat,lon) + s(E_w)
#                   y_var ~ s(lat,lon) + s(LST)
#                y_var ~ s(lat,lon) + s(DISTOC)
#                   y_var ~ s(lat,lon) + s(LC1)
#               y_var ~ s(lat,lon) + s(CANHGHT)
#     y_var ~ s(lat,lon) + s(LST) + ti(LST,LC1)
# y_var ~ s(lat,lon) + s(LST) + ti(LST,CANHGHT)
list_data_s <-extract_list_from_list_obj(raster_prediction_obj_2$validation_mod_obj,"data_s")
list_data_v <-extract_list_from_list_obj(raster_prediction_obj_2$validation_mod_obj,"data_v")
names_mod <- names(raster_prediction_obj_2$method_mod_obj[[1]][[y_var_name]]) #names of models to plot

#number of observations per day
year_nbv <- sapply(list_data_v,FUN=length)
year_nbs <- sapply(list_data_s,FUN=length)
nb_df <- data.frame(nv=year_nbv,ns=year_nbs)
nb_df$n_tot <- year_nbv + year_nbs
range(nb_df$n_tot)

data_v_test <- list_data_v[[1]]

#Convert sp data.frame and combined them in one unique df, see function define earlier
data_v_combined <-convert_spdf_to_df_from_list(list_data_v) #long rownames
names_var<-c("res_mod1","res_mod2","res_mod3","res_mod4","res_mod5","res_mod6","res_mod7","res_mod8","res_mod9","res_mod10")

t<-melt(data_v_combined,
        measure=names_var, 
        id=c("id"),
        na.rm=T)

#hist(data_v_combined)
names(data_v_combined)

mae_fun<-function(x){mean(abs(x))} #Mean Absolute Error give a residuals vector
sd_abs_fun<-function(x){sd(abs(x))} #sd Absolute Error give a residuals vector

mae_tb<-cast(t,id~variable,mae_fun) #join to station location...

met_obj <-load_obj(file.path(in_dir1,met_obj_file_1))
stat_loc<-readOGR(dsn=in_dir1,layer=sub(".shp","",basename(met_obj$loc_stations)))

data_v_mae <-merge(mae_tb,stat_loc,by.x=c("id"),by.y=c("STAT_ID"))
hist(data_v_mae$res_mod1)
mean(data_v_mae$res_mod1)

coords<- data_v_mae[c('longitude','latitude')]              #Define coordinates in a data frame
CRS_interp<-proj4string(data_v_test)
coordinates(data_v_mae)<-coords                      #Assign coordinates to the data frame
proj4string(data_v_mae)<- proj4string(stat_loc)                #Assign coordinates reference system in PROJ4 format
data_v_mae<-spTransform(data_v_mae,CRS(CRS_interp))     #Project from WGS84 to new coord. system
elev <- subset(s_raster,"elev_s")
list_p_mae <- vector("list", length(names_mod))
    
for (k in 1:length(names_mod)){
  
  model_name <- names_mod[k]
  res_model_name <- paste("res",model_name,sep="_")
  #res_model_name <- "res_mod1"
  p1 <- levelplot(elev,scales = list(draw = FALSE), colorkey = FALSE,par.settings = GrTheme)
  tt <- na.omit(data_v_mae)
  tt <- data_v_mae[data_v_mae[res_model_name],]
  p2 <- bubble(data_v_mae,res_model_name, main=paste("MAE per station for ",model_name,sep=""),
               col=matlab.like(25),na.rm=TRUE)
  p3 <- p2 + p1 + p2 #to force legend...
  print(p3)
  #p <- bubble(data_v_mae,"res_mod1",maxsize=4,col=c("red"),fill=FALSE)
  #p<-bubble(data_v_mae,"res_mod1",maxsize=4,col=c("red"),fill=FALSE,key.entries=c(1,1.5,2,2.5,3,3.5,4,4.5))
  list_p_mae[[k]] <- p3
}

layout_m<-c(1,3) # works if set to this?? ok set the resolution...
png(paste("Figure13_paper_","_residuals_",date_selected,"_",out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

grid.arrange(list_p_mae[[1]],list_p_mae[[2]],list_p_mae[[5]],ncol=3)
      
dev.off()   

infile_reg_outline <- "/data/project/layers/commons/data_workflow/inputs/region_outlines_ref_files/OR83M_state_outline.shp"  #input region outline defined by polygon: Oregon
reg_outline <- readOGR(dsn=dirname(infile_reg_outline),layer=sub(".shp","",basename(infile_reg_outline)))

p + layer(sp.polygons(reg_outline,lwd=0.9,col="darkgray"))

p4 <- bubble(data_v_mae,"res_mod4",maxsize=4,col=c("red"),fill=FALSE)
p4 + layer(sp.polygons(reg_outline,lwd=0.9,col="darkgray"))

col_t <- colorRampPalette(c('blue', 'white', 'red'))

p_elev <-levelplot(subset(s_raster,"elev_s"),margin=FALSE)
p4 <-bubble(data_v_mae[data_v_mae$res_mod4>2.134,],"res_mod4",maxsize=4,col=c("blue"),fill=FALSE)
p_elev + p4 + layer(sp.polygons(reg_outline,lwd=0.9,col="green"))
title("mod4")

p_elev <-levelplot(subset(s_raster,"elev_s"))
p1 <-bubble(data_v_mae[data_v_mae$res_mod1>2.109,],"res_mod1",maxsize=4,col=c("blue"),fill=FALSE)
p_elev + p1 + layer(sp.polygons(reg_outline,lwd=0.9,col="green"))
#bubble(data_v_mae,"res_mod1")
#p<-spplot(data_v_mae,"res_mod1",maxsize=4,col=c("red"))
#p
#stations that are outliers in one but not the other...
id_setdiff<-setdiff(data_v_mae[data_v_mae$res_mod1>2.109,]$id,data_v_mae[data_v_mae$res_mod4>2.134,]$id)

data_id_setdiff <- data_v_mae[data_v_mae$id %in% id_setdiff,]

p_elev +layer(sp.polygons(reg_outline,lwd=0.9,col="green")) + layer(sp.points(data_id_setdiff,pch=4,cex=2,col="pink"))


#raster_prediciton object for baseline 1 () s(lat,lon) + s(elev)) and baseline 2 (slat,lon))      
list_data_v <- lapply(1:365,function(i){raster_prediction_obj_2$validation_mod_obj[[i]]$data_v} #testing with residuals
l_formulas<-(extract_from_list_obj(raster_prediction_obj_2$method_mod_obj,"formulas")) #get vector of dates

test <- do.call(rbind,list_data_v)        
data_v_ag <-test                
plot(res_mod1~elev,data=test)                    
plot(res_mod2~elev,data=test)                    
cor(test$res_mod1,test$elev)
cor(test$res_mod2,test$elev,use="complete.obs") #decrease in corellation when using elev
cor(test$res_mod1,test$LST,,use="complete.obs")
cor(test$res_mod5,test$LST,use="complete.obs") #decrease in correlation when using LST

plot(res_mod1~LST,data=test)                    
plot(res_mod5~LST,data=test)                    
                     
brks<-c(0,500,1000,1500,2000,2500,4000)
lab_brks<-1:6
brks<-c(0,500,1000,1500,2000)
lab_brks<-1:4
elev_rcstat<-cut(data_v_ag$elev,breaks=brks,labels=lab_brks,right=F)
boxplot(data_v_ag$res_mod2~elev_rcstat,ylim=c(-15,15))           
boxplot(data_v_ag$res_mod1~elev_rcstat,ylim=c(-15,15))                      
boxplot(data_v_ag$res_mod5~elev_rcstat,ylim=c(-15,15))                      

#brks<-c(0,20,40,60,80,100)
#lab_brks<-1:5
#rcstat<-cut(data_v_ag$LC1,breaks=brks,labels=lab_brks,right=F)
#plot(data_v_ag$res_mod5~rcstat,ylim=c(-5,5))                      
#plot(data_v_ag$res_mod5~rcstat,ylim=c(-5,5))
                      
brks<-c(-20,-10,0,10,20,30,40)
lab_brks<-1:6
rcstat<-cut(data_v_ag$LST,breaks=brks,labels=lab_brks,right=F)
boxplot(data_v_ag$res_mod4 ~ rcstat,ylim=c(-15,15))           
                    
y_range<-range(c(diff_fc))
x_range<-range(c(elev_rcstat))
plot(elev_rcstat,res_mod, ylab="diff_cf", xlab="ELEV_SRTM (m) ",
       ylim=y_range, xlim=x_range)
text(elev_rcstat,diff_cf,labels=data_vf$idx,pos=3)
grid(lwd=0.5,col="black")
title(paste("Testing stations residuals fusion vs Elevation",date_selected,sep=" "))

lat<- aggregate(lat~id,data=data_v_ag,FUN=mean)
lon<- aggregate(lon~id,data=data_v_ag,FUN=mean)                     

data_v_mae <- data.frame(lon=lon$lon,lat=lat$lat)
coordinates(data_v_mae)<- cbind(lon$lon,lat$lat)

mae_fun<-function(x){mean(abs(x))}                      
test<- aggregate(res_mod1~id,data=data_v_ag,FUN=mae_fun)
data_v_mae$mod1 <- test$res_mod1  
test<- aggregate(res_mod2~id,data=data_v_ag,FUN=mae_fun)
data_v_mae$mod2 <- test$res_mod2   
p<- bubble(data_v_mae,"mod1")
print(p)                    
p1 <- levelplot(elev,scales = list(draw = FALSE), colorkey = FALSE,par.settings = GrTheme)
hist(data_v_mae$mod1)
p<- bubble(data_v_mae,"mod1")
spplot()
                     
p_elev <- levelplot(elev)
pc <-p+p1 + p
print(pc)

t<-melt(data_v_ag,
        measure=names_var, 
        id=c("res_mod1","res_mod2,"id"),
        na.rm=T)

n_tb<-cast(t,res_mod1_rc1+id~variable,mae_fun)
n_tb_tot <-cast(t,id~variable,length) #number of times the stations was used for validation

###################### END OF SCRIPT #######################


