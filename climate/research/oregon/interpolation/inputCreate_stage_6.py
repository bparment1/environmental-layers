import sys, os
from optparse import OptionParser
import glob
import shutil
import math

from osgeo import ogr

#Creates a shell script that can be used to run stage 2 and 3 of the climate layer scripts. It sets up all the 
#input arguments for master_script based on provided arguments.
#Inputs:
#    0) inLatLonList: A text file that includes lat_lon of the tiles you want to process,one per line in format 10.0_-120.0 
#    1) outputFolder:  The path to the dir where you want to place the output files (e.g. /nobackupp4/climateLayers/output20Deg/reg4/
#    2) outDirForShFiles: The directory where you want to place the shell script serial files.
#    3) filePerDir: How many sh files to put in each subdirectory, this is useful if you want to run a subset.
#    4) lowStatCount: What's the minimum station count you want to use
#    5) higStatCount: The max station count you want to use
#    6) scriptFn: The filename of the master script you want to use to run stage4. 
#              For example 'Rscript /climateCode/environmental-layers/climate/research/world/interpolation/master_script_09012014_normal.R'
#    7) prefix: e.g. "r4" prefix for the job 
#    8) year: year predicted, this can be changed into a list
#
#Output:
#    Places number of files into directory.
#/nobackupp6/aguzman4/climateLayers/out/
#
#
#
#

def main():
  #usage = "usage: %prog [options] <inLatLonList> <outputFolder> <outDirForShFiles> <filesPerDir> <lowFeatureCount> <highFeatureCount> <scriptFn> <prefix> <year>\n"+\
  #          "Prepare input for stage 4 "
  #var <- args[1] # variable being interpolated #param 1, arg 1
  #in_dir1 <- args[2] #param 5, arg 2
  #region_name <- args[3] #param 6, arg 3
  #out_prefix <- args[4] #param 7, arg 4
  #out_dir <- args[5] #param 8, arg 5
  #create_out_dir_param <- args[6] #param 9, arg 6
  #list_year_predicted <- args[7] # param 10, arg 7, replaced by yearInt here
  #num_cores <- args[8] #number of cores used # param 13, arg 8
  #max_mem <- args[9] #param 21   

  #python /nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/runSetup/inputCreate_stage_4.py 
  #/nobackupp6/aguzman4/climateLayers/inputLayers/regionCSV/list_reg4.csv 
  #/nobackupp6/aguzman4/climateLayers/out/reg4/ 
  #/nobackupp6/aguzman4/climateLayers/out/reg4/serialStage6/2010/ 
  #10 
  #0 
  #250000 
  #/nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/master_script_stage_4.R 
  #r4 
  #2010

  inFn = "/nobackupp6/aguzman4/climateLayers/inputLayers/regionCSV/list_reg4.csv"
  outputFolder = "/nobackupp6/aguzman4/climateLayers/out/reg4/" 
  outDir = "/nobackupp6/aguzman4/climateLayers/out/reg4/serialStage6/2010/"
  fPer = 10 #subdivide the job for NEX
  lowTresh = 0
  hiTresh = int(250000)
  scriptFn = "/nobackupp8/bparmen1/env_layers_scripts/master_script_stage_6_01182016.R"
  #scriptFn = "/nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/master_script_stage_4.R" 
  prefix = r4
  year = 2010 #give a range of date here
  var = "TMAX" # variable being interpolated #param 1, arg 1
  in_dir1 = "/nobackupp6/aguzman4/climateLayers/out/" #param 5, arg 2
  region_name = "reg4" #param 6, arg 3
  out_prefix = "run_global_analyses_pred_12282015" #param 7, arg 4
  out_dir = "/nobackupp8/bparmen1/" #param 8, arg 5
  create_out_dir_param = "TRUE" #param 9, arg 6
  list_year_predicted = 2010 # param 10, arg 7
  num_cores = 6 #number of cores used # param 13, arg 8
  max_mem = 1e+07 #param 21, arg 9
  
  usage = "usage: %prog [options] <inLatLonList> <outputFolder> <outDirForShFiles> <filesPerDir> <lowFeatureCount> <highFeatureCount> <scriptFn> <prefix> "+\         
          " <var> <in_dir1> <region_name> <out_prefix> <out_dir> <create_out_dir_param> <yearInt> <num_cores> <max_mem> "+\
            "Prepare input for stage 4 "
  parser = OptionParser(usage=usage)
  opts,args=parser.parse_args()
  if len(args)<2:
     sys.exit(parser.print_help())

  inFn = args[0]
  outputFolder = args[1]
  outDir = args[2]
  fPer = int(args[3]) #subdivide the job for NEX
  lowTresh = int(args[4])
  hiTresh = int(args[5])
  scriptFn = args[6]
  prefix = args[7]
  year = args[8] #give a range of date here

  inPtr = open(inFn,"r") 
  if inPtr is None:
    print "Error: %s doesnt' exist" % inFn
  lines = inPtr.readlines() #tile directory one per line, 28 for region 4
  inPtr.close()

  found = 0
  notFound = 0
  
  if os.path.exists(outDir) is False:
    print "Warning: %s doesn't exists, creating" % outDir
    os.mkdir(outDir)
    
  wrap = "%s/wrapper.sh" % outDir
  if os.path.exists(wrap) is False:
    print "Warning: %s doesn't exists, creating" % wrap
    wrapFn = os.path.abspath(os.path.dirname(sys.argv[0]))+"/wrapper.sh"
    shutil.copyfile(wrapFn,wrap)
    os.chmod(wrap,0750) #set permission...

  qs = "%s/serial_TEMPLATE.qsub" % outDir
  if os.path.exists(qs) is False:
    print "Warning: %s doesn't exists, creating" % qs
    qsFn = os.path.abspath(os.path.dirname(sys.argv[0]))+"/serial_TEMPLATE.qsub"
    shutil.copyfile(qsFn,qs)
    os.chmod(wrap,0750)

  #rankCount = 0
  #subDirNum = (len(lines)/fPer) + 2


  #for i in range(0,subDirNum):
  #  subDirName = "%s/dirSub_%d" % (outDir,i)
  #  if os.path.exists(subDirName) is False:
  #    print "Creating %s" % subDirName
  #    os.mkdir(subDirName)
  subDirName = "%s/dirSub_*" % outDir
  dirs = glob.glob(subDirName)
  for d in dirs:
    print "Removing %s" % d
    shutil.rmtree(d)

  
  dirCounter = 0
  
  shpFiles = []
  
  #Open the shapefiles to look at actual number of stations
  for l in lines:
    l = l.rstrip().split(",")
    print l[0]
    #41.3_-84.9/daily_covariates_ghcn_data_TMAX_2010_201141.3_-84.9.shp'
    searchStr = "%s/%s/%s/daily_covariates_ghcn_data_*%s.shp" % (outputFolder,l[0],year,l[0])

    shpFn = glob.glob(searchStr)
    if len(shpFn) > 0:
      fCount = featureCount(shpFn[0])
      
      if fCount is not None:
        shpFiles.append([shpFn[0],fCount,l[0]])
    else:
       #notFound += 1
       print "Error finding %s" % searchStr
  
  sortedShpFiles = sorted(shpFiles,key=lambda x: x[1]) #sort by numer of stations...
  #print sortedShpFiles

  stCFn = "%s/stationCounts.txt" % outDir
  fPtrSt = open(stCFn,"w+")
  if fPtrSt is None:
    print "Error opening %s" % stCFn
    sys.exit(0)

  for s in sortedShpFiles:
   outSt = "%s,%d\n" % (s[2],s[1])
   fPtrSt.write(outSt)
  fPtrSt.close()


  stDim = len(sortedShpFiles)-1
   
  stMin = int(sortedShpFiles[0][1])
  stMax = int(sortedShpFiles[stDim][1])
  
  print "Min %d, Max %d " % (stMin,stMax)
 
  bins = []

  #Every 15k about 2 mins for daily
  minNeed = 3
  minAdd = 2
  tStep = 30000
  tStart = 15000
  tEnd = 300000
  for s in range(tStart,tEnd,tStep):
    print "%d-%d range needs %d minutes" % (s-tStep,s,minNeed)
    #continue

    bin1 =  [i for i in sortedShpFiles if ((i[1] < s) and (i[1] > (s-tStep)))]
    #print bin1
    rows = len(bin1)
    if rows > 0:
      #cols = len(bin1[0])
      #print rows,cols
      for b in bin1:
         b.append(minNeed)
      bins.append(bin1)
    else:
      print "No values for %d-%d range" % (s-tStep,s) 
    
    minNeed = minNeed + minAdd  
  
  dirCount = -1
  count = 0
  for binSingle in bins:
    rankCount = 0

    subDirName = "%s/dirSub_%d" % (outDir,dirCounter)
    if os.path.exists(subDirName) is True:
      if os.listdir(subDirName):
         dirCounter += 1
 

    subDirName = "%s/dirSub_%d" % (outDir,dirCounter)
    if os.path.exists(subDirName) is False:
      print "Creating %s" % subDirName
      os.mkdir(subDirName)

    print count,subDirName
    if len(binSingle) > 0:
      b0 = float(binSingle[0][3])
      #(minutesPerTile*(365/cores)+(climatologyStep) + assesment)/minutes
      hours = math.ceil(((b0 * (365.0/10.0))+60.0+30.0)/60.0)
      fPtrMinu = open(subDirName+"/numMinutes.txt","w+")
      if fPtrMinu is None:
        print "Error opening "+ subDirName+"/numMinutes.txt"
      else:
        print "Creating "+ subDirName+"/numMinutes.txt"
        fPtrMinu.write(str(hours)+'\n')
        fPtrMinu.close() 

    count += 1
    
    ##Input files that contain Rscript command
    #binSingle: this contains a list of tiles to run ranked by time 
    #run by year 
    #for b in binSingle:
    for year in yearInt:
      #ll = b[2]#"%s_%s" % (baseSpl[1],baseSpl[2]) 
      #check for existence of files from stage2
      #metTest = "%s/%s/%s/met_stations_outfiles_obj_gam_CAI_%s.RData" % (outputFolder,ll,year,ll)
 
      #if os.path.exists(metTest) is False:
      #   print "No met object %s" % metTest
      #   continue     
      #prediction object check
      #methodTest =  "%s/%s/%s/method_mod_obj_gam_CAI_dailyTmax%s.RData" % (outputFolder,ll,year,ll)
      #if os.path.exists(methodTest) is True:
      #   print "Method object exists"
      #   continue 
    
      #yearInt = int(year)+1 #interval of year e.g. 201-2015
      #var <- args[1] # variable being interpolated #param 1, arg 1
      #in_dir1 <- args[2] #param 5, arg 2
      #region_name <- args[3] #param 6, arg 3
      #out_prefix <- args[4] #param 7, arg 4
      #out_dir <- args[5] #param 8, arg 5
      #create_out_dir_param <- args[6] #param 9, arg 6
      #list_year_predicted <- args[7] # param 10, arg 7
      #num_cores <- args[8] #number of cores used # param 13, arg 8
      #max_mem <- args[9] #param 21
      list_year_predicted = list(year) #list with one element
      #"Rscript %s %s %s %s %s %s %s %s" % (var,in_dir1,region_name,out_prefix,out_dir,create_out_dir_param,list_year_predicted,num_cores,max_mem)
      #var = "TMAX" # variable being interpolated #param 1, arg 1
      #in_dir1 = "/nobackupp6/aguzman4/climateLayers/out/" #param 5, arg 2
      #region_name = "reg4" #param 6, arg 3
      #out_prefix = "run_global_analyses_pred_12282015" #param 7, arg 4
      #out_dir = "/nobackupp8/bparmen1/" #param 8, arg 5
      #create_out_dir_param = "TRUE" #param 9, arg 6
      #list_year_predicted = c(2010) # param 10, arg 7
      #num_cores = 6 #number of cores used # param 13, arg 8
      #max_mem = 1e+07 #param 21, arg 9
         
      #scriptFn= /nobackupp8/bparmen1/env_layers_scripts/master_script_stage_6_01182016.R   
      outStr = "Rscript %s %s %s %s %s %s %s %s" % (scriptFn,var,in_dir1,region_name,out_prefix,out_dir,create_out_dir_param,list_year_predicted,num_cores,max_mem)

      #outStr = "Rscript %s %s wgs84Grid %s %s %s %s/subset/mean_LST_%s_jan_6_wgs84.tif FALSE %s/%s/covar_obj_%s.RData %s/%s/%s/met_stations_outfiles_obj_gam_CAI_%s.RData 10 4800  %s %s > %s/outLogs/%s_stage4_%s.log 2>  %s/outLogs/%s_stage4_err_%s.log" % (scriptFn,ll,ll,outputFolder,b[0],outputFolder,ll,outputFolder,ll,ll,outputFolder,ll,year,ll,year,yearInt,outputFolder,ll,year,outputFolder,ll,year)
      #print outStr

      outFnSt = "%s/dirSub_%d/input_%d.in" % (outDir,dirCounter,rankCount)
      outPtr = open(outFnSt,"w+")
      if outPtr is None:
         print "Error: Opening %s" % outFnSt
      else:
         outPtr.write("#!/bin/bash"+"\n")
         outPtr.write(outStr+"\n")
         outPtr.close()
      
         if rankCount >= fPer:
           dirCounter += 1
           rankCount = 0

           subDirName = "%s/dirSub_%d" % (outDir,dirCounter)
           if os.path.exists(subDirName) is False:
             print "Creating %s" % subDirName
             os.mkdir(subDirName)

           #if len(binSingle) > 0:
           b0 = float(binSingle[0][3])
           hours = math.ceil(((b0 * (365.0/10.0))+60.0+30.0)/60.0)
           fPtrMinu = open(subDirName+"/numMinutes.txt","w+")
           if fPtrMinu is None:
             print "Error opening "+ subDirName+"/numMinutes.txt"
           else:
             print "Creating "+ subDirName+"/numMinutes.txt"
             fPtrMinu.write(str(hours)+'\n')
             fPtrMinu.close() 
         else:
           rankCount += 1

      found += 1
      outPtr.close()
      #found += 1
  prefix3 = prefix + "_" + year[2:4]
  postQsubCreator(outDir,prefix3)

  print "Found %d" % found
  print "Not found %d" % notFound
  print "End" 
   
  sys.exit(0)
  
def postQsubCreator(outDir,prefix):
  subDirName = "%s/dirSub_*" % outDir
  dirs = glob.glob(subDirName)
  
  if len(dirs) > 0:
   for d in dirs:
     fn = "%s/numMinutes.txt" % d
     if os.path.exists(fn): 
       fPtr = open(fn)
       numHrs = float(fPtr.readlines()[0])
       outDirTmp = os.path.dirname(fn)
       tmpRpl = "%sdirSub_" % outDir
       rank = int(outDirTmp.replace(tmpRpl,''))
       cores = len(glob.glob(outDirTmp+"/input*"))
       prefix2 = prefix + "_"+str(rank)
       print outDir,rank,cores,numHrs,prefix2           

       createQsub(outDir,rank,cores,numHrs,prefix2)
       fPtr.close()
     else:
       print "File %s doesn't exists" % fn
       continue

def featureCount(fn):
  driver = ogr.GetDriverByName('ESRI Shapefile')
  ds = driver.Open(fn, 0) # 0 means read-only. 1 means writeable.

  # Check to see if shapefile is found.
  if ds is None:
    print 'Could not open %s' % (fn)
    return None
  else:
    print 'Opened %s' % (fn)
    layer = ds.GetLayer()
    featureCount = layer.GetFeatureCount()
    print "Number of features in %s: %d" % (os.path.basename(fn),featureCount)

  return featureCount

def createQsub(outDir,rank,cores,time,prefix):
  fn = "%s/qsub_%d.qsub" % (outDir,rank)
  if time <= 2:
    qType = "devel"
  elif time <= 8:
    qType = "normal"
  elif time > 8:
    qType = "long" 
  else:
    return None
   
  fPtr = open(fn,"w+")
  if fPtr is None:
    print "Error opening %s " % fn
    return None
  
  fPtr.write("#PBS -S /bin/bash\n")
  out1 = "#PBS -l select=%d:ncpus=14:model=ivy\n" % cores
  fPtr.write(out1)
  out1 = "#PBS -l walltime=%d:00:00\n" % time
  fPtr.write(out1)
  fPtr.write("#PBS -j n\n")
  fPtr.write("#PBS -m be\n")
  out1 = "#PBS -N %s\n" % prefix # 
  fPtr.write(out1)
  fPtr.write("#PBS -V\n")
  out1 = "#PBS -q %s\n" % qType
  fPtr.write(out1)
  
  #ELayers Project ID
  fPtr.write("#PBS -W group_list=s1557\n\n")

  out1 = "CORES=%d\n" % cores
  fPtr.write(out1)
  out1 = "SUBDIRNUM=%d\n" % rank
  fPtr.write(out1)
  out1 = "HDIR=%s\n\n" % outDir
  fPtr.write(out1)

  fPtr.write("source /u/aguzman4/sharedModules/etc/environ.sh\n\n")

  fPtr.write("module load mpi-mvapich2/1.4.1/intel\n\n")

  fPtr.write("#$PBS_O_WORKDIR\n\n")

  fPtr.write("LOGSTDOUT=$HDIR/climate_ivy_serial_%d.stdout\n" % rank)
  fPtr.write("LOGSTDERR=$HDIR/climate_ivy_serial_%d.stderr\n\n" % rank)

  fPtr.write("mpiexec -pernode -comm none -np $CORES $HDIR/wrapper.sh $SUBDIRNUM $HDIR 1> $LOGSTDOUT 2> $LOGSTDERR\n")

  fPtr.close()

if __name__ == '__main__':
   main()

################################ #

#python /nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/runSetup/inputCreate_stage_4.py /nobackupp6/aguzman4/climateLayers/inputLayers/regionCSV/list_reg4.csv /nobackupp6/aguzman4/climateLayers/out/reg4/ /nobackupp6/aguzman4/climateLayers/out/reg4/serialRuns/2014/ 10 0 250000 /nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/master_script_stage_4.R r4 2014

#parameters
#python script:/nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/runSetup/inputCreate_stage_4.py 
#inLatLonList: /nobackupp6/aguzman4/climateLayers/inputLayers/regionCSV/list_reg4.csv 
#outputFolder: /nobackupp6/aguzman4/climateLayers/out/reg4/ 
#outDirForShFiles: /nobackupp6/aguzman4/climateLayers/out/reg4/serialRuns/2014/ 
#filePerDir: 10 
#lowStatCount: 0 
#higStatCount: 250000 
#scriptFn: /nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/master_script_stage_4.R 
#r4 
#2014

#    inLatLonList: A text file that includes lat_lon of the tiles you want to process,one per line in format 10.0_-120.0 
#    outputFolder:  The path to the dir where you want to place the output files (e.g. /nobackupp4/climateLayers/output20Deg/reg4/
#    outDirForShFiles: The directory where you want to place the shell script serial files.
#    filePerDir: How many sh files to put in each subdirectory, this is useful if you want to run a subset.
#    lowStatCount: What's the minimum station count you want to use
#    higStatCount: The max station count you want to use
#    scriptFn: The filename of the master script you want to use to run stage4. 
#              For example 'Rscript /climateCode/environmental-layers/climate/research/world/interpolation/master_script_09012014_normal.R'


#python /nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/runSetup/inputCreate_stage_4.py /nobackupp6/aguzman4/climateLayers/inputLayers/regionCSV/list_reg4.csv /nobackupp6/aguzman4/climateLayers/out/reg4/ /nobackupp6/aguzman4/climateLayers/out/reg4/serialStage6/2010/ 10 0 250000 /nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/master_script_stage_4.R r4 2010
#python /nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/runSetup/master_script_stage_6_01182016.R /nobackupp6/aguzman4/climateLayers/inputLayers/regionCSV/list_reg4.csv /nobackupp6/aguzman4/climateLayers/out/reg4/ /nobackupp6/aguzman4/climateLayers/out/reg4/serialStage6/2010/ 10 0 250000 /nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/master_script_stage_6_01182016.R r4 2010

#Calculate times it takes by dates and region? 
#Add all bin count together for each year
#loop over year?

#inFn = "/nobackupp6/aguzman4/climateLayers/inputLayers/regionCSV/list_reg4.csv"
#outputFolder = "/nobackupp6/aguzman4/climateLayers/out/reg4/" 
#outDir = "/nobackupp6/aguzman4/climateLayers/out/reg4/serialStage6/2010/"
#fPer = 10 #subdivide the job for NEX
lowTresh = 0
hiTresh = int(250000)
scriptFn = "/nobackupp6/aguzman4/climateLayers/finalCode/environmental-layers/climate/research/world/interpolation/master_script_stage_4.R" 
prefix = r4
year = 2010 #give a range of date here
var = "TMAX" # variable being interpolated #param 1, arg 1
in_dir1 = "/nobackupp6/aguzman4/climateLayers/out/" #param 5, arg 2
region_name = "reg4" #param 6, arg 3
out_prefix = "run_global_analyses_pred_12282015" #param 7, arg 4
out_dir = "/nobackupp8/bparmen1/" #param 8, arg 5
create_out_dir_param = "TRUE" #param 9, arg 6
list_year_predicted = 2010 # param 10, arg 7
num_cores = 6 #number of cores used # param 13, arg 8
max_mem = 1e+07 #param 21, arg 9


inPtr = open(inFn,"r") 
if inPtr is None:
  print "Error: %s doesnt' exist" % inFn
lines = inPtr.readlines()
inPtr.close()
  
found = 0
notFound = 0
  
if os.path.exists(outDir) is False:
    print "Warning: %s doesn't exists, creating" % outDir
    os.mkdir(outDir)
    
wrap = "%s/wrapper.sh" % outDir #this is the wrapper to launch all the codes together
if os.path.exists(wrap) is False:
    print "Warning: %s doesn't exists, creating" % wrap
    wrapFn = os.path.abspath(os.path.dirname(sys.argv[0]))+"/wrapper.sh"
    shutil.copyfile(wrapFn,wrap)
    os.chmod(wrap,0750)

qs = "%s/serial_TEMPLATE.qsub" % outDir
if os.path.exists(qs) is False:
    print "Warning: %s doesn't exists, creating" % qs
    qsFn = os.path.abspath(os.path.dirname(sys.argv[0]))+"/serial_TEMPLATE.qsub"
    shutil.copyfile(qsFn,qs)
    os.chmod(wrap,0750)

