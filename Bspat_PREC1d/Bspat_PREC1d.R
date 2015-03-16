# --~- Bspat_PREC1d.R -~--
# Bayesian Spatial Interpolation of daily cumulated precipitation
#..............................................................................
# == Command line ==
#  $>R --vanilla yyyy.mm.dd yyyy.mm.dd 
#                   begin ---> end accumulation period 
#
# == Time specification ==
# Timezone is UTC: hour [0,23]. timestamp = end of accumulation period.
#  (i.e. 2014/09/01 12 -> precipitation sum 2014/09/01 11:01 2014/09/01 12:00)
#
# == Grid specifications ==
# 2 Grids:
# 1. FG -> Fine Grid, which is the high-resolution output grid masked on Norwegian mainland
# 2. CG -> Coarse Grid, which is the low-resolution grid larger than Norway
# High-Resolution Fine Grid 1Km grid (FG)
#> print(orog)
# class       : RasterLayer 
# dimensions  : 1550, 1195, 1852250  (nrow, ncol, ncell)
# resolution  : 1000, 1000  (x, y)
# extent      : -75000, 1120000, 6450000, 8e+06  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 
# data source : in memory
# names       : fenno_dem_u33 
# values      : 0, 2250  (min, max)
# 
# Coarse Grid (CG)
#
#> print(cbind(xy.CG[1:3,1],xy.CG[1:3,2]))
#       [,1]    [,2]
#[1,] -72500 7997500 coordinates of grid.CG[1,1] -> gridpoint is the center of the cell
#[2,] -67500 7997500 coordinates of grid.CG[1,2]
#[3,] -62500 7997500 coordinates of grid.CG[1,3]
#
#
#> print(cbind(xy.CG[1:3,1],xy.CG[1:3,2],rc.CG[1:3,1],rc.CG[1:3,2]))
#       [,1]    [,2] [,3] [,4]
#[1,] -72500 7997500    1    1
#[2,] -67500 7997500    1    2
#[3,] -62500 7997500    1    3
#
# row(y)  col(x)->      1       2    ... nx
# 1                   (1,1)   (1,2)  ... (1,nx)
# 2                   (2,1)   (2,2)  ... (2,nx)            
# ..
# ny                  (ny,1)  (ny,2)  ... (ny,nx)
#
# upper-left corner (xmin,ymax) -> upper-left corner of grid-point (1,1)
#> print(cbind(xy[1:3,1],xy[1:3,2]))
#        [,1]    [,2]
# [1,] -74500 7999500
# [2,] -73500 7999500
# [3,] -72500 7999500
#
# 
# == OUTPUT:  netCDF files + text files ==
# NB loof for @@@@ in the code and you'll get the output
#
# name: seNorge_v2test_PREC_grid_yyyymmdd.b_yyyymmdd.e.nc
#       [var=out.file.grd.ana] -> analysis field
#        
# name: seNorge_v2test_PREC_grid_normIDI_yyyymmdd.b_yyyymmdd.e.nc
#       [var=out.file.grd.idi] -> normalized IDI field
#
# name: seNorge_v2test_PREC_station_yyyymmdd.b_yyyymmdd.e.txt
#       [var=out.file.stn]
# ASCII file with data in columns separated by ";". Header:
# year;month;day;hour;nhour;stid;x;y;z;...
#  1     2    3    4    5     6  7 8 9
# ... eve;yo;yb;ya;yav;yidi;yidiv;dqcflag
#      10 11 12 13  14  15   16     17
# description:
# year/month/day hour => timestamp
# nhour => number of hours accumulated
# stid => station identifier
# x,y => geographical coordinates (UTM 33, Km)
# z => elevation [m]
# eve => Contiguous Rain Area
# yo => observation [mm]
# yb => Background value used in the last iteration [mm]
# ya => analysis [mm] (prediction)
# yav => Cross-Validated (CV) analysis [mm]
# yidi => Integral Data Influence (IDI) (accumulated and normalized)
# yidiv => CV IDI (accumulated and normalized)
# dqcflag => Data Quality Control (DQC) Flag
#
# name: seNorge_v2test_PREC_event_yyyymmdd.b_yyyymmdd.e.txt
#       [var=out.file.eve]
# ASCII file with data in columns separated by ";". Header:
# year;month;day;hour;nhour;eve.lab;nobs;...
#   1    2    3    4    5      6      7 
# ... area;volume;mean.idi.x;mean.idi.y;mean.idiv.y;
#      8       9     10         11        12   
# ... mean.rain;max.rain.x;max.rain.yo;max.rain.ya;max.rain.yav;
#        13        14          15            16        17       
# ... x.loc;y.loc;s.maj.ax;s.min.ax;dir.s.maj.ax;...
#      18     19      20       21         22
# ... cv.rel.all;cv.bias.all;cv.rmse.all;cv.made.all;...
#            23             24        25        26
# ... cv.rel.q50;cv.bias.q50;cv.rmse.q50;cv.made.q50;mean.idiv.y.q50;n.q50;
#            27             28        29        30              31         32
# ... cv.rel.q75;cv.bias.q75;cv.rmse.q75;cv.made.q75;mean.idiv.y.q75;n.q75;
#            33             34        35        36              37         38
# ... idi.norm.fac;
#           39 
# description:
# year/month/day hour => timestamp
# nhour => number of hours accumulated
# eve.lab => eve label
# nobs => number of observations in the eve
# area => eve areas (Km**2 == number of gridpoints if gridcell dimensions are 1x1 Km)
# volume => eve areas * eve rain field [mm * Km**2]
# mean.idi.x => average of the IDI field
# mean.idi.y => average of IDI at station locations
# mean.idiv.y => average of CV IDI at station locations
# mean.rain => average of rain field
# max.rain.x => maximum rain intensity in the rain field
# max.rain.yo => maximum rain intensity in the observations
# max.rain.ya => maximum rain intensity in the analysis 
# max.rain.yav => maximum rain intensity in the CV analysis 
# x.loc => x coordinate of the ellipsoid hull 
# “ellipsoid hull” or “spanning ellipsoid”, i.e. the ellipsoid of minimal volume
#  (‘area’ in 2D) such that all given points lie just inside or on the boundary
#  of the ellipsoid) center
# y.loc => y coordinate of the ellipsoid hull center
# s.maj.ax => semi-major axis of the ellipsoid hull
# s.min.ax => semi-minor axis of the ellipsoid hull
# dir.s.maj.ax => semi-major axis direction
#   dir=0 N-S orientation; dir=45 NE-SW; dir=90 E-W; dir=135 NW-SE; dir=180 N-S
# cv.rel.all  => mean(Cross-Validated relative error) set=all observation in that OAP
# cv.bias.all => mean(CV bias) set=all observation in that OAP
# cv.rmse.all => mean(CV root mean square error) set=all observation in that OAP
# cv.made.all => 1.4826 * Median Absolut Deviation set=all observation in that OAP
# cv.rel.q50  => mean(Cross-Validated relative error) set=observation in OAP >=q50.daily
# cv.bias.q50 => mean(CV bias) set=observation in OAP >=q50.daily
# cv.rmse.q50 => mean(CV root mean square error) set=observation in OAP >=q50.daily
# mean.idiv.y.q50 => mean(CV IDI for observation in OAP >=q50.daily)
# n.q50       => number of observation in OAP >=q50.daily
# cv.made.q50 => 1.4826 * Median Absolut Deviation set=observation in OAP >=q50.daily
# cv.rel.q75  => mean(Cross-Validated relative error) set=observation in OAP >=q75.daily
# cv.bias.q75 => mean(CV bias) set=observation in OAP >=q75.daily
# cv.rmse.q75 => mean(CV root mean square error) set=observation in OAP >=q75.daily
# cv.made.q75 => 1.4826 * Median Absolut Deviation set=observation in OAP >=q75.daily
# mean.idiv.y.q50 => mean(CV IDI for observation in OAP >=q75.daily)
# n.q50       => number of observation in OAP >=q75.daily
# idi.norm.fac => IDI normalization factor for that event 
# 
# == DQC flag ==
# pos  value  description
#
# History:
# 03.12.2014 - Cristian Lussana. Original code.
# 26.02.2015 - Cristian Lussana. original code from Bspat_TAMRR_v1_0.R
#  change log: 
#  - allow for the use of observations outside Norway
#  - revisied queries to KDVH
#  - geographical information on seNorge2_dem_UTM33.nc
#  - definition of PREC1d 
#  - definition of a new directory tree
#
# DQCflag
# NA  observation is NA
# -1  missing DQC info
# 0   good observation
# 100 bad  observation: KDVH flag > 2 | observation not good in external DQC | 
#                       observed value not plausible | station in blacklist/s  
# 200 bad  observation: dry-station surrounded only by wet-stations (close enough)
# 300 bad  observation: wet-stations surrounded only by dry-stations (close enough)
# 400 bad  observation: dry observation is (1) not included in a dry area
#                       (2) is in Norway 
# 500 bad  observation: wet observation is (1) not included in an event (2) in Norway
# -----------------------------------------------------------------------------
rm(list=ls())
# Libraries
library(raster)
rasterOptions(chunksize=2e+07,maxmemory=3e+08)
library(igraph)
library(rgdal)
library(ncdf)
require(tripack)
require(cluster)
# CRS strings
proj4.wgs84<-"+proj=longlat +datum=WGS84"
proj4.ETRS_LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
#..............................................................................
# Functions
# + st.log:
st.log<-function(x,lc=NULL) {
  if (is.null(lc)) lc<-1.5
  if (is.vector(x)) {
    res<-vector(mode="numeric",length=length(x))
    aux<-which(x>lc)
    if (length(aux)>0) {
      res[aux]<-log10(x[aux])
      res[-aux]<-log10(lc)+(x[-aux]-lc)/(lc*log(10))
    } else {
      res<-log10(lc)+(x-lc)/(lc*log(10))
    }
  } else {
    if (x>lc) {
      res<-log10(x)
    } else {
      res<-log10(lc)+(x-lc)/(lc*log(10))
    }
  }
  return(res)
}
# Gaussian filter for square cells
fgauss <- function(sigma, n) {
#  sigma<-7
#  n<-21
  m <- matrix(nc=n, nr=n)
  col <- rep(1:n, n)
  row <- rep(1:n, each=n)
  x <- col - ceiling(n/2)
  y <- row - ceiling(n/2)
  # according to http://en.wikipedia.org/wiki/Gaussian_filter
  m[cbind(row, col)] <- 1/(2*pi*sigma^2) * exp(-(x^2+y^2)/(2*sigma^2))
  # sum of weights should add up to 1     
  m / sum(m)
}
# manage fatal error
error_exit<-function(str=NULL) {
  print("Fatal Error:")
  if (!is.null(str)) print(str)
  quit(status=1)
}
#-------------------------------------------------------------------
# [] Setup parameters
# eps2 meaning: believe in observations
eps2<-0.1
# rain/norain threshold
rr.inf<-0.1
# Contiguous NO-Rain and Rain Area (first step) 
Lsubsample.max<-20
Lsubsample.DHmax<-150
# identify adjacent nodes - Establish (local-) triangulation (Delauney)
n.sector<-16
sector.angle<-360/n.sector
# DQC thresholds 
yo.dqc.plausible.min<-0 # mm/h
yo.dqc.plausible.max<-500 # mm/h
# evaluation thresholds
#q50.hourly<-0.4 #mm/h (weak/normal)
#q75.hourly<-1.0 #mm/h (normal/moderate-strong)
q50.daily<-3 #mm/daily (weak/normal)
q75.daily<-8 #mm/daily (normal/moderate-strong)
# definition of Dh sequence - horizontal de-correlation length scale
min.Dh.seq.allowed<-10 # Km
#Dh.seq.ref<-100 # Km
Dh.seq.ref<-50 # Km
#Dh.seq.reference<-c(5000,4000,3000,2500,2000,1750,1500,1400,1300,1200,
#                    1100,1000, 900, 800, 700, 650, 600, 550, 500, 450,
#                     400, 350, 300, 260, 230, 200, 180, 160, 140, 120,
#                     100,  90,  80,  70,  60,  50,  40,  30,  20,  10)
Dh.seq.reference<-c(5000,4000,3000,2000,1000, 900, 800, 700, 600, 500,
                     400, 300, 200, 150, 100,  80,  60,  40,  20,  10)
#Dh.seq.reference<-c(500, 100,  40,  20,  10)
Dz.seq.gt50<-c(5000,1000) # m
Dz.seq.le50<-c(1000,500) # m
# analysis cycle on FG
ndim.FG.iteration<-10000
#ndim.FG.iteration<-5000
# stations outside Norway
max.Km.stnINdomain<-300
# MAIN ========================================================================
# [] Read command line arguments
arguments <- commandArgs()
print("Arguments")
print(arguments)
date.b.string<-arguments[3]
date.e.string<-arguments[4]
file_blacklist_current<-arguments[5]
file_blacklist_never<-arguments[6]
file_errobs<-arguments[7]
config_file<-arguments[8]
config_par<-arguments[9]
if (length(arguments)!=9) 
  ext<-error_exit(paste("Error in command line arguments: \n",
  " R --vanilla yyyy.mm.dd yyyy.mm.dd blacklist_current blacklist_never errobs configFILE configPAR \n",
  "             begin -----> end of the accumulation period",sep=""))
# [] define/check paths
if (!file.exists(config_file)) ext<-error_exit("Fatal Error: configuration file not found")
source(config_file)
for (p in 1:length(config_list)) {
  if (config_list[[p]]$pname==config_par) break
}
if (p==length(config_list) & (config_list[[p]]$pname!=config_par) )  
  ext<-error_exit("Fatal Error: configuration parameter not in configuration file")
main.path<-config_list[[p]]$opt$main.path
main.path.output<-config_list[[p]]$opt$main.path.output
testmode<-config_list[[p]]$opt$testmode
if ( !(file.exists(main.path)) | !(file.exists(main.path.output)) ) 
  ext<-error_exit("Fatal Error: path not found")
# geographical information
main.path.geoinfo<-paste(main.path,"/geoinfo",sep="")
filenamedem<-paste(main.path.geoinfo,"/seNorge2_dem_UTM33.nc",sep="")
filenamedem.CG<-paste(main.path.geoinfo,"/fennodem_utm33.nc",sep="")
if (!file.exists(paste(main.path.geoinfo,"/seNorge2_dem_UTM33.nc",sep=""))) 
  ext<-error_exit(paste("File not found:",main.path.geoinfo,"/seNorge2_dem_UTM33.nc"))
if (!file.exists(paste(main.path.geoinfo,"/fennodem_utm33.nc",sep=""))) 
  ext<-error_exit(paste("File not found:",main.path.geoinfo,"/fennodem_utm33.nc"))
# common libs and etcetera
path2lib.com<-paste(main.path,"/lib",sep="")
path2etc.com<-paste(main.path,"/etc",sep="")
if (!file.exists(paste(path2lib.com,"/nogrid.ncout.R",sep=""))) 
  ext<-error_exit(paste("File not found:",path2lib.com,"/nogrid.ncout.R"))
if (!file.exists(paste(path2lib.com,"/ncout.spec.list.r",sep=""))) 
  ext<-error_exit(paste("File not found:",path2lib.com,"/ncout.spec.list.r"))
if (!file.exists(paste(path2lib.com,"/getStationData.R",sep=""))) 
  ext<-error_exit(paste("File not found:",path2lib.com,"/getStationData.R"))
source(paste(path2lib.com,"/nogrid.ncout.R",sep=""))
source(paste(path2lib.com,"/ncout.spec.list.r",sep=""))
source(paste(path2lib.com,"/getStationData.R",sep=""))
# test mode
print(testmode)
if (testmode) {
  print("TESTMODE TESTMODE TESTMODE")
  if (file.exists(paste(main.path,"/Bspat_PREC1d/testbed",sep=""))) {
    testbed<-paste(main.path,"/Bspat_PREC1d/testbed",sep="")
    station.info<-paste(testbed,"/station_data.csv",sep="")
    observed.data<-paste(testbed,"/observed_data.csv",sep="")
  } else {
    ext<-error_exit(paste("testbed not found"))
  }
}
# output directories
dir.create(file.path(main.path.output,"seNorge2"), showWarnings = FALSE)
dir.create(file.path(main.path.output,"seNorge2_addInfo"), showWarnings = FALSE)
path2output.main<-paste(main.path.output,"/seNorge2/PREC1d",sep="")
path2output.main.stn<-paste(path2output.main,"/station_dataset",sep="")
path2output.main.grd<-paste(path2output.main,"/gridded_dataset",sep="")
path2output.add<-paste(main.path.output,"/seNorge2_addInfo/PREC1d",sep="")
path2output.add.grd<-paste(path2output.add,"/gridded_dataset",sep="")
path2output.add.eve<-paste(path2output.add,"/event_dataset",sep="")
if (!(file.exists(path2output.main)))     dir.create(path2output.main,showWarnings=F) 
if (!(file.exists(path2output.main.stn))) dir.create(path2output.main.stn,showWarnings=F) 
if (!(file.exists(path2output.main.grd))) dir.create(path2output.main.grd,showWarnings=F) 
if (!(file.exists(path2output.add)))      dir.create(path2output.add,showWarnings=F) 
if (!(file.exists(path2output.add.grd)))  dir.create(path2output.add.grd,showWarnings=F) 
if (!(file.exists(path2output.add.eve)))  dir.create(path2output.add.eve,showWarnings=F) 
# netcdf fixed parameters
grid.type <- "utm33"
source.nc<-"daily precipitation from station data"
var.version.xa <- "1.0"
prod.date <- substr(Sys.time(),1,10)
pname.xa<-"RR"
for (p in 1:length(ncout.spec.list)) {
  if (ncout.spec.list[[p]]$pname==pname.xa) break
}
if (p==length(ncout.spec.list) & (ncout.spec.list[[p]]$pname!=pname.xa) ) {
  flag.write.xa<-F
} else {
  flag.write.xa<-T
  var.name.xa<-ncout.spec.list[[p]]$opts$var.name
  var.longname.xa<-ncout.spec.list[[p]]$opts$var.longname
  var.unit.xa<-ncout.spec.list[[p]]$opts$var.unit
  var.mv.xa<-as.numeric(ncout.spec.list[[p]]$opts$var.mv)
  times.unit.xa <-ncout.spec.list[[p]]$opts$t.unit
  times.ref.xa <-ncout.spec.list[[p]]$opts$t.ref
  reference.xa <- ncout.spec.list[[p]]$opts$reference
}
# [] set Time-related variables
start.string<-paste(date.b.string,sep="")
end.string<-paste(date.e.string,sep="")
start.string.day<-paste(substr(start.string,1,10),".01",sep="")
end.string.day  <-paste(substr(end.string,1,10),".24",sep="")
start <- strptime(start.string,"%Y.%m.%d","UTC")
end <- strptime(end.string,"%Y.%m.%d","UTC")
timeseq<-as.POSIXlt(seq(as.POSIXlt(start),as.POSIXlt(end),by="1 hour"),"UTC")
nhour<-length(timeseq)
start.day <- strptime(start.string.day,"%Y.%m.%d.%H","UTC")
end.day <- strptime(end.string.day,"%Y.%m.%d.%H","UTC")
dayseq<-as.POSIXlt(seq(as.POSIXlt(start.day),as.POSIXlt(end.day),by="1 day"),"UTC")
nday<-length(dayseq)
yyyy.b<-timeseq$year[1]+1900
mm.b<-timeseq$mon[1]+1
dd.b<-timeseq$mday[1]
yyyy.e<-timeseq$year[nhour]+1900
mm.e<-timeseq$mon[nhour]+1
dd.e<-timeseq$mday[nhour]
date.b<-paste(dd.b,".",mm.b,".",yyyy.b,sep="")
date.e<-paste(dd.e,".",mm.e,".",yyyy.e,sep="")
yyyymmdd.b<-paste(yyyy.b,formatC(mm.b,width=2,flag="0"),formatC(dd.b,width=2,flag="0"),sep="")
yyyymm.b<-paste(yyyy.b,formatC(mm.b,width=2,flag="0"),sep="")
datestring.b<-paste(yyyy.b,".",mm.b,".",dd.b," ",sep="")
yyyymmdd.e<-paste(yyyy.e,formatC(mm.e,width=2,flag="0"),formatC(dd.e,width=2,flag="0"),sep="")
yyyymm.e<-paste(yyyy.e,formatC(mm.e,width=2,flag="0"),sep="")
datestring.e<-paste(yyyy.e,".",mm.e,".",dd.e," ",sep="")
# Setup output files 
dir.create(paste(path2output.main.stn,"/",yyyymm.b,sep=""),showWarnings=F)
dir.create(paste(path2output.main.grd,"/",yyyymm.b,sep=""),showWarnings=F)
dir.create(paste(path2output.add.grd,"/",yyyymm.b,sep=""),showWarnings=F)
dir.create(paste(path2output.add.eve,"/",yyyymm.b,sep=""),showWarnings=F)
out.file.stn<- paste(path2output.main.stn,"/",yyyymm.b,
                     "/seNorge_v2_0_PREC1d_station_",
                     yyyymmdd.b,"_",yyyymmdd.e,".txt",sep="")
out.file.eve<- paste(path2output.add.eve,"/",yyyymm.b,
                     "/seNorge_v2_0_PREC1d_event_",
                     yyyymmdd.b,"_",yyyymmdd.e,".txt",sep="")
out.file.grd.ana<- paste(path2output.main.grd,"/",yyyymm.b,
                         "/seNorge_v2_0_PREC1d_grid_",
                         yyyymmdd.b,"_",yyyymmdd.e,".nc",sep="")
out.file.grd.idi<- paste(path2output.add.grd,"/",yyyymm.b,
                         "/seNorge_v2_0_PREC1d_grid_normIDI_",
                         yyyymmdd.b,"_",yyyymmdd.e,".nc",sep="")
#
print("Output files:")
print("analysis on the grid (netcdf)")
print(out.file.grd.ana)
print("event-normalized idi on the grid (netcdf)")
print(out.file.grd.idi)
print("station outputs (text)")
print(out.file.stn)
print("event outputs (text)")
print(out.file.eve)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# define header for the station data output file
cat(paste("year","month","day","nday","stid",
          "x","y","z","eve.lab","yo",
          "yb","ya","yav","yidi","yidiv","dqcflag","\n",sep=";"),
          file=out.file.stn,append=F)
cat(paste("year","month","day","nday","eve.lab","nobs",
          "area","volume",
          "mean.idi.x","mean.idi.y","mean.idiv.y",
          "mean.rain","max.rain.x",
          "max.rain.yo","max.rain.ya","max.rain.yav",
          "x.loc","y.loc","s.maj.ax","s.min.ax","dir.s.maj.ax",
          "cv.rel.all","cv.bias.all","cv.rmse.all","cv.made.all",
          "cv.rel.q50","cv.bias.q50","cv.rmse.q50","cv.made.q50",
          "mean.idiv.y.q50","n.q50",
          "cv.rel.q75","cv.bias.q75","cv.rmse.q75","cv.made.q75",
          "mean.idiv.y.q75","n.q75",
          "idi.norm.fac",
          "\n",sep=";"),
          file=out.file.eve,append=F)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [] Grid
# CRS Coordinate Reference System
r.orog.FG<-trim(raster(filenamedem))
nx.FG<-ncol(r.orog.FG)
ny.FG<-nrow(r.orog.FG)
dx.FG<-xres(r.orog.FG)
dy.FG<-yres(r.orog.FG)
area.1cell.FG<-dx.FG*dy.FG/(1000*1000) # Area Km**2
# 4 borders point (SW corner (xmn,ymn); NE corner (xmx,ymx))
xmn.FG<-xmin(r.orog.FG)
xmx.FG<-xmax(r.orog.FG)
ymn.FG<-ymin(r.orog.FG)
ymx.FG<-ymax(r.orog.FG)
#
r.orog.CG<-trim(raster(filenamedem.CG))
nx.CG<-ncol(r.orog.CG)
ny.CG<-nrow(r.orog.CG)
dx.CG<-xres(r.orog.CG)
dy.CG<-yres(r.orog.CG)
xmn.CG<-xmin(r.orog.CG)
xmx.CG<-xmax(r.orog.CG)
ymn.CG<-ymin(r.orog.CG)
ymx.CG<-ymax(r.orog.CG)
# extract all the cell values: zvalues[1] contains the orog[1,1] value
# Raster: cell numbers start at 1 in the upper left corner,
# and increase from left to right, and then from top to bottom
zvalues.FG<-extract(r.orog.FG,1:ncell(r.orog.FG))
xy<-xyFromCell(r.orog.FG,1:ncell(r.orog.FG))
x.FG<-sort(unique(xy[,1]))
y.FG<-sort(unique(xy[,2]),decreasing=T)
rc<-rowColFromCell(r.orog.FG,1:ncell(r.orog.FG))
aux<-as.vector(zvalues.FG)
mask.FG<-which(!is.na(aux))
zgrid<-aux[mask.FG]
xgrid<-xy[mask.FG,1]
ygrid<-xy[mask.FG,2]
rowgrid<-rc[mask.FG,1]
colgrid<-rc[mask.FG,2]
Lgrid.FG<-length(mask.FG)
# CG
zvalues.CG<-extract(r.orog.CG,1:ncell(r.orog.CG))
xy.CG<-xyFromCell(r.orog.CG,1:ncell(r.orog.CG))
rc.CG<-rowColFromCell(r.orog.CG,1:ncell(r.orog.CG))
aux.CG<-as.vector(zvalues.CG)
mask.CG<-which(!is.na(aux.CG))
zgrid.CG<-aux.CG[mask.CG]
xgrid.CG<-xy.CG[mask.CG,1]
ygrid.CG<-xy.CG[mask.CG,2]
rowgrid.CG<-rc.CG[mask.CG,1]
colgrid.CG<-rc.CG[mask.CG,2]
Lgrid.CG<-length(mask.CG)
fg<-fgauss(sigma=3.5,n=11)
# clean memory
rm(zvalues.CG,zvalues.FG)
rm(xy,rc,aux,rowgrid,colgrid)
rm(xy.CG,rc.CG,aux.CG)
# debug info
print("grid parameters")
print(paste("nx ny dx dy",as.integer(nx.FG),as.integer(ny.FG),round(dx.FG,2),round(dy.FG,2)))
print(paste("xmn xmx ymn ymx",round(xmn.FG,2),round(xmx.FG,2),round(ymn.FG,2),round(ymx.FG,2)))
print(paste("Lgrid.FG",as.integer(Lgrid.FG)))
print("grid.CG parameters")
print(paste("nx.CG ny.CG dx.CG dy.CG",as.integer(nx.CG),as.integer(ny.CG),round(dx.CG,2),round(dy.CG,2)))
print(paste("xmn.CG xmx.CG ymn.CG ymx.CG",round(xmn.CG,2),round(xmx.CG,2),round(ymn.CG,2),round(ymx.CG,2)))
print(paste("Lgrid.CG",as.integer(Lgrid.CG)))
#------------------------------------------------------------------------------
# [] Read Station Information 
# conditions:
# 1. stations in KDVH
# 2. stations in CG
if (!testmode) {
  stations.tmp<-getStationMetadata(from.year=yyyy.b,to.year=yyyy.b,
                                   max.Km=max.Km.stnINdomain)
} else {
  stations.tmp<-read.csv(file=station.info)
}
# check if in CG
n.tmp<-length(stations.tmp$stnr)
aux<-extract(r.orog.CG,cbind(stations.tmp$x,stations.tmp$y),na.rm=T)
stn.out.CG<-vector(length=n.tmp)
stn.out.CG[1:n.tmp]<-F
stn.out.CG[which(is.na(aux))]<-T
# definitive station list
L.y.tot<-length(which(!stn.out.CG))
stations<-data.frame(matrix(nrow=L.y.tot,ncol=5))
names(stations)<-c("stnr","z","x","y","NO")
stations$stnr<-stations.tmp$stnr[which(!stn.out.CG)]
stations$NO<-stations.tmp$NO[which(!stn.out.CG)]
stations$z<-stations.tmp$z[which(!stn.out.CG)]
stations$y<-stations.tmp$y[which(!stn.out.CG)]
stations$x<-stations.tmp$x[which(!stn.out.CG)]
# [] define Vectors and Matrices
VecX<-vector(mode="numeric",length=L.y.tot)
VecY<-vector(mode="numeric",length=L.y.tot)
VecZ<-vector(mode="numeric",length=L.y.tot)
VecS<-vector(mode="numeric",length=L.y.tot)
Disth<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
Distz<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
n.yo<-vector(mode="numeric",length=L.y.tot)
yo<-vector(mode="numeric",length=L.y.tot)
y.eve<-vector(mode="numeric",length=L.y.tot)
#yo.n<-vector(mode="numeric",length=L.y.tot)
xb<-vector(mode="numeric",length=Lgrid.FG)
ya<-vector(mode="numeric",length=L.y.tot)
yav<-vector(mode="numeric",length=L.y.tot)
ydqc<-vector(mode="numeric",length=L.y.tot)
ydqc.flag<-vector(mode="numeric",length=L.y.tot)
y.mask<-vector(mode="numeric",length=L.y.tot)
Lsubsample.vec<-vector(mode="numeric",length=L.y.tot)
D<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
S<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
# fill VecX/Y/Z/S
VecX<-as.numeric(as.vector(stations$x))
VecY<-as.numeric(as.vector(stations$y))
VecZ<-as.numeric(as.vector(stations$z))
VecS<-as.numeric(as.vector(stations$stnr))
# identify stations without dem (FG)
z.FG<-extract(r.orog.FG,cbind(VecX,VecY),na.rm=T)
stn.out.FG<-vector(length=L.y.tot)
stn.out.FG[1:L.y.tot]<-F
stn.out.FG[which(is.na(z.FG))]<-T
# [] compute Disth and Distz (symmetric) matrices: 
#  Disth(i,j)=horizontal distance between i-th station and j-th station [Km]
#  Distz(i,j)=elevation difference between i-th station and j-th station [m]
Disth<-(outer(VecY,VecY,FUN="-")**2.+outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
Distz<-abs(outer(VecZ,VecZ,FUN="-"))
print("list of station ids")
print(cbind(1:L.y.tot,VecS))
# Stations on the Norwegian mainland
stn.NO<-which(stations$NO & !stn.out.FG)
n.stn.NO<-length(stn.NO)
#------------------------------------------------------------------------------
# Elaborations
# Get and sum daily observations from KDVH
n.yo[]<-NA
yo[]<-NA
# set DQC flag to 0 for every station/observation
ydqc.flag[]<-NA
# cycle over days to get data from KDVH
if (!testmode) {
  if (nday==1) {
    data<-getStationData(var="RR", from.yyyy=yyyy.b, from.mm=mm.b, from.dd=dd.b,
                         to.yyyy=yyyy.e, to.mm=mm.e, to.dd=dd.e,
                         qa=NULL, statlist=stations, outside.Norway=T,
                         err.file=file_errobs, blist.perm=file_blacklist_never,
                         blist.curr=file_blacklist_current, verbose=T,
                         val.min.allowed=yo.dqc.plausible.min,
                         val.max.allowed=yo.dqc.plausible.max)
  } else {
    data<-getStationData(var="RR", from.yyyy=yyyy.b, from.mm=mm.b, from.dd=dd.b,
                         to.yyyy=yyyy.e, to.mm=mm.e, to.dd=dd.e,
                         qa=NULL, statlist=stations, outside.Norway=T,
                         err.file=file_errobs, blist.perm=file_blacklist_never, 
                         blist.curr=file_blacklist_current, verbose=T,
                         val.min.allowed=yo.dqc.plausible.min, 
                         val.max.allowed=yo.dqc.plausible.max,fun="sum")
  }
} else {
  data<-read.csv(file=observed.data)
}
# debug: start
print("Observed data (alligned with station info):")
print(data)
# debug: end
yo<-data$value
y.notNA<-which(!is.na(yo))
L.y.notNA<-length(y.notNA)
for (i in 1:L.y.tot) {
  if (is.na(data$value[i])) next
  ydqc.flag[i]<--1
  if (!is.na(data$KDVHflag[i])) if (data$KDVHflag[i]>2) ydqc.flag[i]<-100
  if (!is.na(data$KDVHflag[i])) if (data$KDVHflag[i]==0 | data$KDVHflag[i]==1 | data$KDVHflag[i]==2) ydqc.flag[i]<-0
  if (data$ntime[i]!=data$nvalue[i] |
      !data$plausible[i] |
      data$err.ext[i] |
      data$blist.perm[i] |
      data$blist.curr[i]) ydqc.flag[i]<-100
}
print(data)
#data.names<-c("stnr","year","month","day","hour","ntime",
# "value","nvalue","DQC","KDVHflag","plausible","err.ext","blist.perm","blist.curr")
print(paste("number of station having at least one available observation =",L.y.notNA,"(tot=",L.y.tot,")"))
print(paste("number of station having good observation           (so far)=",length(which(ydqc.flag<=0 & !is.na(yo)))))
print(paste("  # station having at least one erroneous observation (plausibility check) =",length(which(!data$plausible & !is.na(yo)))))
print(paste("  # station having at least one erroneous observation           (KDVH DQC) =",length(which(data$KDVHflag>2 & !is.na(yo)))))
print(paste("  # station having at least one erroneous observation       (external DQC) =",length(which(data$err.ext & !is.na(yo)))))
print(paste("  # station blacklisted for the time period considered =",length(which(data$blist.curr & !is.na(yo)))))
print(paste("  # station blacklisted (permanently) =",length(which(data$blist.perm & !is.na(yo)))))
print(paste("  # station in masked areas =",length(which(stn.out.CG & !is.na(yo)))))
# define Vectors and Matrices
xidi.CG.wet<-vector(mode="numeric",length=Lgrid.CG)
xidi.CG.dry<-vector(mode="numeric",length=Lgrid.CG)
# loop for DQC
yo.ok.pos<-which(ydqc.flag<=0 & !is.na(yo))
L.yo.ok<-length(yo.ok.pos)
while (L.yo.ok>0) {
# vector with the positions (pointers to VecS) of presumably good observations 
# daily observations [L.yo.ok]
  yo.ok.pos<-which(ydqc.flag<=0 & !is.na(yo))
  L.yo.ok<-length(yo.ok.pos)
# vectors with the positions (pointers to VecS) of valid (i.e. not NAs & DQC-ok)
# daily observations and "wet" or "dry" [L.yo.ok.wet,L.yo.ok.dry]
  yo.ok.wet<-ydqc.flag<=0 & !is.na(yo) & yo>=rr.inf
  yo.ok.dry<-ydqc.flag<=0 & !is.na(yo) & yo<rr.inf
  yo.ok.pos.wet<-which(yo.ok.wet)
  yo.ok.pos.dry<-which(yo.ok.dry)
  L.yo.ok.wet<-length(yo.ok.pos.wet)
  L.yo.ok.dry<-length(yo.ok.pos.dry)
  print(paste("observation not NAs, presumably good =",
              L.yo.ok,"(wet=",L.yo.ok.wet,"dry=",L.yo.ok.dry,")"))
#  print(cbind(VecS[yo.ok.pos.wet],yo[yo.ok.pos.wet]))
# case of no-rain over the whole domain
  if (L.yo.ok.wet==0) {
    print("no rain over the whole domain")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    cat(paste(yyyy.b,mm.b,dd.b,nday,
              round(VecS[stn.NO],0),
              round(VecX[stn.NO],0),
              round(VecY[stn.NO],0),
              round(VecZ[stn.NO],0),
              rep(NA,n.stn.NO), #eve.lab
              round(yo[stn.NO],1),
              rep(NA,n.stn.NO), #yb
              rep(0,n.stn.NO), #ya
              rep(NA,n.stn.NO), #yav
              rep(NA,n.stn.NO), #yidi.eve
              rep(NA,n.stn.NO), #yidiv.eve
              round(ydqc.flag[stn.NO],2),
              "\n",sep=";"),file=out.file.stn,append=T)
    # Figures
    r.aux.FG <-raster(ncol=nx.FG, nrow=ny.FG, xmn=xmn.FG, xmx=xmx.FG,
                  ymn=ymn.FG, ymx=ymx.FG, crs=proj4.utm33)
    xa.FG<-vector(mode="numeric",length=Lgrid.FG)
    xa.FG[]<-0
    r.aux.FG[mask.FG]<-round(xa.FG,1)
    nogrid.ncout(grid=t(as.matrix(r.aux.FG)),
                 x=x.FG,y=y.FG,grid.type=grid.type,
                 file.name=out.file.grd.ana,
                 var.name=var.name.xa,
                 var.longname=var.longname.xa,
                 var.unit=var.unit.xa,
                 var.mv=var.mv.xa,
                 var.version=var.version.xa,
                 times=c(paste(yyyymmdd.b,"0000",sep="")),times.unit=times.unit.xa,
                 times.ref=times.ref.xa,
                 prod.date=prod.date,
                 reference=reference.xa,
                 proj4.string="+proj=utm +zone=33 +ellps=WGS84",
                 source.string=source.nc)
    quit(status=0)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  } # end case of no-rain over the whole domain
  # [] Contiguous NO-Rain Area: First Guess
  Lsubsample.vec[]<-NA
  # vectors used in the identification of NO-Rain (nor) - first guess
  # ...aux-> temporary vectors
  # ...vec-> final vectors
  # nnor.1st.vec-> number of nor.1st
  # nor.1st.vec-> list of stations (i.e. pointers to VecS) inside nor.1st
  #           matrix(i,j) i=1,nnor.1st.vec j=1,Lnor.1st.vec[i]
  #   note: an nor.1st could contain both wet and dry stations
  # Lnor.1st.vec-> vector(i) (i=1,..,nnor.1st.vec) number of stations inside i-th nor.1st 
  nor.1st.vec<-matrix(data=NA,ncol=L.yo.ok,nrow=L.yo.ok)
  Lnor.1st.vec<-vector(mode="numeric",length=L.yo.ok)
  Lnor.1st.vec[]<-NA
  nnor.1st.vec<-0
  nor.1st.aux<-matrix(data=NA,ncol=L.yo.ok,nrow=L.yo.ok)
  Lnor.1st.aux<-vector(mode="numeric",length=L.yo.ok)
  Lnor.1st.aux[]<-NA
  nnor.1st.aux<-0
  # Contiguous NO-Rain Area: First Guess. Step I
  for (b in yo.ok.pos.dry) {  # START: Cycle over dry stations 
    # a. identify the closest stations to the b-th dry station
    #  outputs: -close2b-> position (i.e. pointer to VecS) of the closest stations.
    #            constraints: distance must be less than Lsubsample.DHmax 
    #                         max number of stations allowed is Lsubsample.max
    #           -Lsubsample.vec[b]-> actual number of stations used
    close2b.aux<-order(Disth[b,],decreasing=F)
    close2b.au1<-close2b.aux[which(close2b.aux%in%yo.ok.pos)][1:Lsubsample.max]
    Lsubsample.vec[b]<-Lsubsample.max
    if (Disth[b,close2b.au1[Lsubsample.max]]>Lsubsample.DHmax) {
      if (any(Disth[b,close2b.au1]<=Lsubsample.DHmax)) {
        Lsubsample.vec[b]<-max(which(Disth[b,close2b.au1]<=Lsubsample.DHmax))
      } else {
        Lsubsample.vec[b]<-0
      }
    }
    if (Lsubsample.vec[b]<3) Lsubsample.vec[b]<-Lsubsample.max
    close2b<-close2b.au1[1:Lsubsample.vec[b]]
    # b. setup vectors 
    VecX.b<-vector(mode="numeric",length=Lsubsample.vec[b])
    VecY.b<-vector(mode="numeric",length=Lsubsample.vec[b])
    VecZ.b<-vector(mode="numeric",length=Lsubsample.vec[b])
    yo.b<-vector(mode="numeric",length=Lsubsample.vec[b])
    VecX.b<-VecX[close2b]
    VecY.b<-VecY[close2b]
    VecZ.b<-VecZ[close2b]
    wet.b.aux<-close2b %in% yo.ok.pos.wet
    dry.b.aux<-close2b %in% yo.ok.pos.dry
    wet.b<-close2b[wet.b.aux]
    dry.b<-close2b[dry.b.aux]
    nwet<-length(wet.b)
    ndry<-length(dry.b)
    # c. Establish (local-) triangulation (Delauney)
    #    note: all stations are used (wet and dry)
    tri.rr<-tri.mesh(VecX.b,VecY.b)
    # d. identify all the stations (wheter wet or dry) which belongs
    #    to adjacent nodes (respect to the b-th dry station node)
    #  procedure: tri.find-> returns the three vertex indexes of triangle
    #   containing the specified point. To find all the neighbouring 
    #   triangles just span the area surrounding the b-th dry station
    #   (circle of radius=1Km)
    #  output: nodes-> station position respect to VecS
    nodes<-vector(mode="numeric")
    tri.loc<-tri.find(tri.rr,VecX[b],VecY[b])
    # note: due to the fact that (VecX[b],VecY[b]) belongs to the mesh
    # used to establish the triangulation, ...i1,i2,i3 must be >0
    nodes<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
    for (cont.clock in 1:n.sector) {
      x.aux<-VecX[b]+sin(sector.angle*(cont.clock-1)*pi/180.)*1000
      y.aux<-VecY[b]+cos(sector.angle*(cont.clock-1)*pi/180.)*1000
      tri.loc<-tri.find(tri.rr,x.aux,y.aux)
      nodes.aux<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
      nodes<-c(nodes,nodes.aux[which( (!(nodes.aux %in% nodes)) & (nodes.aux>0) )])
    }
    lnodes<-length(nodes)
    # DEBUG END
    # e. update the temporary structure used to identify nor 
    #    merge in a (new) nor all the (old) temporary nor (if any) 
    #    which contain at least one station present in nodes (avoiding
    #    repetition); erase (i.e. set to NAs) the (old) temporary nor merged
    aux<-vector()
    if (nnor.1st.aux>0) {
      for (nor.1st in 1:nnor.1st.aux) {
        if (Lnor.1st.aux[nor.1st]==0) next
        if ( any(nodes %in% nor.1st.aux[nor.1st,1:Lnor.1st.aux[nor.1st]]) ) {
          aux<-c(aux[which(!(aux %in% nor.1st.aux[nor.1st,1:Lnor.1st.aux[nor.1st]]))],nor.1st.aux[nor.1st,1:Lnor.1st.aux[nor.1st]])
          nor.1st.aux[nor.1st,]<-NA
          Lnor.1st.aux[nor.1st]<-0
        }
      }
    }
    nnor.1st.aux<-nnor.1st.aux+1
    Lnor.1st.aux[nnor.1st.aux]<-length(c(aux,nodes[which(!(nodes%in%aux))]))
    nor.1st.aux[nnor.1st.aux,1:Lnor.1st.aux[nnor.1st.aux]]<-c(aux,nodes[which(!(nodes%in%aux))])
    #DEBUG START
#    print("new dry-station group identified")
#    print(nor.1st.aux[nnor.1st.aux,1:Lnor.1st.aux[nnor.1st.aux]])
    #DEBUG END      
  } # END: Cycle over the dry stations
# Contiguous Rain Area: First Guess. Step II
# reorganise nor labels
  y.nor.1st<-vector(length=L.y.tot)
  y.nor.1st[1:L.y.tot]<-NA
  nnor.1st.vec<-0
  for (nor.1st in 1:nnor.1st.aux) {
    if (Lnor.1st.aux[nor.1st]==0) next
    nnor.1st.vec<-nnor.1st.vec+1
    Lnor.1st.vec[nnor.1st.vec]<-Lnor.1st.aux[nor.1st]
    nor.1st.vec[nnor.1st.vec,1:Lnor.1st.vec[nnor.1st.vec]]<-nor.1st.aux[nor.1st,1:Lnor.1st.aux[nor.1st]]
    y.nor.1st[nor.1st.vec[nnor.1st.vec,1:Lnor.1st.vec[nnor.1st.vec]]]<-nnor.1st.vec
  }
# DQC DQC DQC
# dry-stations surrounded only by wet-stations (close enough)
  flag.dry.dqc<-F
  for (j in 1:nnor.1st.vec) {
    nor.1st.j<-nor.1st.vec[j,1:Lnor.1st.vec[j]]
    indx.dry.stns<-which(nor.1st.j %in% yo.ok.pos.dry)
    n.dry<-length(indx.dry.stns)
    if (n.dry>1) next
    nor.1st.j.dry<-nor.1st.j[indx.dry.stns]
    aux.dist<-Disth[nor.1st.j.dry,nor.1st.j]
    min.dist.from.dry.stn<-min(aux.dist[aux.dist>0])
#    if (min.dist.from.dry.stn<min.dist.allowed) {
    if (min.dist.from.dry.stn<=20) {
      flag.dry.dqc<-T
      ydqc.flag[nor.1st.j.dry]<-200
    }
  }
  if (flag.dry.dqc) next
# DQC DQC end
  # [] Contiguous Rain Area: First Guess
  Lsubsample.vec[]<-NA
  # vectors used in the identification of eve - first guess
  # ...aux-> temporary vectors
  # ...vec-> final vectors
  # neve.1st.vec-> number of eve.1st
  # eve.1st.vec-> list of stations (i.e. pointers to VecS) inside eve.1st
  #           matrix(i,j) i=1,neve.1st.vec j=1,Leve.1st.vec[i]
  #   note: an eve.1st could contain both wet and dry stations
  # Leve.1st.vec-> vector(i) (i=1,..,neve.1st.vec) number of stations inside i-th eve.1st 
  eve.1st.vec<-matrix(data=NA,ncol=L.yo.ok,nrow=L.yo.ok)
  Leve.1st.vec<-vector(mode="numeric",length=L.yo.ok)
  Leve.1st.vec[]<-NA
  neve.1st.vec<-0
  eve.1st.aux<-matrix(data=NA,ncol=L.yo.ok,nrow=L.yo.ok)
  Leve.1st.aux<-vector(mode="numeric",length=L.yo.ok)
  Leve.1st.aux[]<-NA
  neve.1st.aux<-0
  # Contiguous Rain Area: First Guess. Step I
  for (b in yo.ok.pos.wet) {  # START: Cycle over wet stations 
    # a. identify the closest stations to the b-th wet station
    #  outputs: -close2b-> position (i.e. pointer to VecS) of the closest stations.
    #            constraints: distance must be less than Lsubsample.DHmax 
    #                         max number of stations allowed is Lsubsample.max
    #           -Lsubsample.vec[b]-> actual number of stations used
    close2b.aux<-order(Disth[b,],decreasing=F)
    close2b.au1<-close2b.aux[which(close2b.aux%in%yo.ok.pos)][1:Lsubsample.max]
    Lsubsample.vec[b]<-Lsubsample.max
    if (Disth[b,close2b.au1[Lsubsample.max]]>Lsubsample.DHmax) {
      if (any(Disth[b,close2b.au1]<=Lsubsample.DHmax)) {
        Lsubsample.vec[b]<-max(which(Disth[b,close2b.au1]<=Lsubsample.DHmax))
      } else {
        Lsubsample.vec[b]<-0
      }
    }
    if (Lsubsample.vec[b]<3) Lsubsample.vec[b]<-Lsubsample.max
    close2b<-close2b.au1[1:Lsubsample.vec[b]]
    # b. setup vectors 
    VecX.b<-vector(mode="numeric",length=Lsubsample.vec[b])
    VecY.b<-vector(mode="numeric",length=Lsubsample.vec[b])
    VecZ.b<-vector(mode="numeric",length=Lsubsample.vec[b])
    yo.b<-vector(mode="numeric",length=Lsubsample.vec[b])
    VecX.b<-VecX[close2b]
    VecY.b<-VecY[close2b]
    VecZ.b<-VecZ[close2b]
    wet.b.aux<-close2b %in% yo.ok.pos.wet
    dry.b.aux<-close2b %in% yo.ok.pos.dry
    wet.b<-close2b[wet.b.aux]
    dry.b<-close2b[dry.b.aux]
    nwet<-length(wet.b)
    ndry<-length(dry.b)
    # c. Establish (local-) triangulation (Delauney)
    #    note: all stations are used (wet and dry)
    tri.rr<-tri.mesh(VecX.b,VecY.b)
    # d. identify all the stations (wheter wet or dry) which belongs
    #    to adjacent nodes (respect to the b-th wet station node)
    #  procedure: tri.find-> returns the three vertex indexes of triangle
    #   containing the specified point. To find all the neighbouring 
    #   triangles just span the area surrounding the b-th wet station
    #   (circle of radius=1Km)
    #  output: nodes-> station position respect to VecS
    nodes<-vector(mode="numeric")
    tri.loc<-tri.find(tri.rr,VecX[b],VecY[b])
    # note: due to the fact that (VecX[b],VecY[b]) belongs to the mesh
    # used to establish the triangulation, ...i1,i2,i3 must be >0
    nodes<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
    for (cont.clock in 1:n.sector) {
      x.aux<-VecX[b]+sin(sector.angle*(cont.clock-1)*pi/180.)*1000
      y.aux<-VecY[b]+cos(sector.angle*(cont.clock-1)*pi/180.)*1000
      tri.loc<-tri.find(tri.rr,x.aux,y.aux)
      nodes.aux<-close2b[c(tri.loc$i1,tri.loc$i2,tri.loc$i3)]
      nodes<-c(nodes,nodes.aux[which( (!(nodes.aux %in% nodes)) & (nodes.aux>0) )])
    }
    lnodes<-length(nodes)
    # DEBUG END
    # e. update the temporary structure used to identify eve 
    #    merge in a (new) eve all the (old) temporary eve (if any) 
    #    which contain at least one station present in nodes (avoiding
    #    repetition); erase (i.e. set to NAs) the (old) temporary eve merged
    aux<-vector()
    if (neve.1st.aux>0) {
      for (eve.1st in 1:neve.1st.aux) {
        if (Leve.1st.aux[eve.1st]==0) next
        if ( any(nodes %in% eve.1st.aux[eve.1st,1:Leve.1st.aux[eve.1st]]) ) {
          aux<-c(aux[which(!(aux %in% eve.1st.aux[eve.1st,1:Leve.1st.aux[eve.1st]]))],eve.1st.aux[eve.1st,1:Leve.1st.aux[eve.1st]])
          eve.1st.aux[eve.1st,]<-NA
          Leve.1st.aux[eve.1st]<-0
        }
      }
    }
    neve.1st.aux<-neve.1st.aux+1
    Leve.1st.aux[neve.1st.aux]<-length(c(aux,nodes[which(!(nodes%in%aux))]))
    eve.1st.aux[neve.1st.aux,1:Leve.1st.aux[neve.1st.aux]]<-c(aux,nodes[which(!(nodes%in%aux))])
    #DEBUG START
#    print("new wet-station group identified")
#    print(eve.1st.aux[neve.1st.aux,1:Leve.1st.aux[neve.1st.aux]])
    #DEBUG END      
  } # END: Cycle over the wet stations
# Contiguous Rain Area: First Guess. Step II
# reorganise eve labels
  y.eve.1st<-vector(length=L.y.tot)
  y.eve.1st[1:L.y.tot]<-NA
  neve.1st.vec<-0
  for (eve.1st in 1:neve.1st.aux) {
    if (Leve.1st.aux[eve.1st]==0) next
    neve.1st.vec<-neve.1st.vec+1
    Leve.1st.vec[neve.1st.vec]<-Leve.1st.aux[eve.1st]
    eve.1st.vec[neve.1st.vec,1:Leve.1st.vec[neve.1st.vec]]<-eve.1st.aux[eve.1st,1:Leve.1st.aux[eve.1st]]
    y.eve.1st[eve.1st.vec[neve.1st.vec,1:Leve.1st.vec[neve.1st.vec]]]<-neve.1st.vec
  }
# DQC DQC DQC
# wet-stations surrounded only by dry-stations (close enough)
  flag.wet.dqc<-F
  for (j in 1:neve.1st.vec) {
    eve.1st.j<-eve.1st.vec[j,1:Leve.1st.vec[j]]
    indx.wet.stns<-which(eve.1st.j %in% yo.ok.pos.wet)
    n.wet<-length(indx.wet.stns)
    if (n.wet>1) next
    eve.1st.j.wet<-eve.1st.j[indx.wet.stns]
    aux.dist<-Disth[eve.1st.j.wet,eve.1st.j]
    min.dist.from.wet.stn<-min(aux.dist[aux.dist>0])
#    if (min.dist.from.wet.stn<min.dist.allowed) {
    if (min.dist.from.wet.stn<=20) {
      flag.wet.dqc<-T
      ydqc.flag[eve.1st.j.wet]<-300
    }
  }
  if (flag.wet.dqc) next
# GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID GRID
# [] Contiguous NO-Rain Area: First Guess on the grid.CG
  print("++ define nor.1st on the grid.CG")
  lab.nor.CG.1st<-vector(mode="numeric",length=Lgrid.CG)
  lab.nor.CG.1st[]<-NA
  # do the same for Contiguous Rain Area [to be used later]
  lab.eve.CG.1st<-vector(mode="numeric",length=Lgrid.CG)
  lab.eve.CG.1st[]<-NA
#  tri.rr<-tri.mesh(VecX[yo.ok.pos],VecY[yo.ok.pos])
  for (i in 1:Lgrid.CG) {
    Disth.i<-( (ygrid.CG[i]-VecY[yo.ok.pos])**2+(xgrid.CG[i]-VecX[yo.ok.pos])**2 )**0.5 / 1000.
    close2i<-yo.ok.pos[which.min(Disth.i)]
    if (!is.na(y.nor.1st[close2i])) {
      lab.nor.CG.1st[i]<-y.nor.1st[close2i]
    } else {
      lab.nor.CG.1st[i]<-0
    }
    # do the same for Contiguous Rain Area [to be used later]
    if (!is.na(y.eve.1st[close2i])) {
      lab.eve.CG.1st[i]<-y.eve.1st[close2i]
    } else {
      lab.eve.CG.1st[i]<-0
    }
  }
  rm(Disth.i)
  r.CG <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
                        ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
  r.CG[]<-NA
  r.CG[mask.CG]<-lab.nor.CG.1st
# [] Refine identification (Identify Wet and Dry Areas within the nor first guess)
  print("++ Refine nor identification (Identify Wet and Dry Areas within the nor first guess)")
  xa.Dh.nor.1st<-vector(mode="numeric",length=Lgrid.CG)
  xa.Dz.nor.1st<-vector(mode="numeric",length=Lgrid.CG)
  x.CG.nor.dry<-vector(mode="numeric",length=Lgrid.CG)
  n.nor<-0
# x.CG.nor.dry -> =0 gridpoint not classified as dry; =1 grid-point classified as wet
  x.CG.nor.dry[]<-0
  for (nor.1st in 1:nnor.1st.vec) { # START: cycle over nor.1st
#    print(paste(nor.1st,"nor.1st +++++++++++++++++++++++++++++++++"))
    # define vectors / variables / indexes 
    Lnor.1st<-Lnor.1st.vec[nor.1st]
    VecX.nor.1st<-vector(mode="numeric",length=Lnor.1st)
    VecY.nor.1st<-vector(mode="numeric",length=Lnor.1st)
    VecZ.nor.1st<-vector(mode="numeric",length=Lnor.1st)
    VecS.nor.1st<-vector(mode="numeric",length=Lnor.1st)
    yo.nor.1st<-vector(mode="numeric",length=Lnor.1st)
    VecX.nor.1st<-VecX[nor.1st.vec[nor.1st,1:Lnor.1st]]
    VecY.nor.1st<-VecY[nor.1st.vec[nor.1st,1:Lnor.1st]]
    VecZ.nor.1st<-VecZ[nor.1st.vec[nor.1st,1:Lnor.1st]]
    VecS.nor.1st<-VecS[nor.1st.vec[nor.1st,1:Lnor.1st]]
    yo.nor.1st<-yo[nor.1st.vec[nor.1st,1:Lnor.1st]]
    yoindx.nor.1st.wet<-which(yo.nor.1st>=rr.inf & !is.na(yo.nor.1st))
    yoindx.nor.1st.dry<-which(yo.nor.1st<rr.inf & !is.na(yo.nor.1st))
    yoindx.wet<-nor.1st.vec[nor.1st,yoindx.nor.1st.wet]
    yoindx.dry<-nor.1st.vec[nor.1st,yoindx.nor.1st.dry]
    n.nor.1st.wet<-length(yoindx.nor.1st.wet)
    n.nor.1st.dry<-length(yoindx.nor.1st.dry)
    # the next line should not happen...
    if (n.nor.1st.dry==0) next
    # case of an isolated dry station...
    if (n.nor.1st.wet==0) {
      x.CG.nor.dry[which(lab.nor.CG.1st==nor.1st)]<-1
      next
    }
    xindx.nor.CG<-which(lab.nor.CG.1st==nor.1st)
    Lgrid.nor.1st<-length(xindx.nor.CG)
    # define Dz for the current nor (min allowed value is 500 m)
    xa.Dz.nor.1st[xindx.nor.CG]<-max(500,max(zgrid.CG[xindx.nor.CG]))
    for (i in xindx.nor.CG) {    #START: cycle over xindx.nor.CG
      Disth.nor.1st<-( (ygrid.CG[i]-VecY.nor.1st)**2+(xgrid.CG[i]-VecX.nor.1st)**2 )**0.5 / 1000.
      Distz.nor.1st<-abs(zgrid.CG[i]-VecZ.nor.1st)
      close2i.ord<-order(Disth.nor.1st,decreasing=F)
#      xa.Dh.nor.1st[i] <- max(20,Disth.nor.1st[close2i.ord[1]])
      xa.Dh.nor.1st[i] <- Disth.nor.1st[close2i.ord[1]]
      notsofar.wet.i<-which(Disth.nor.1st[yoindx.nor.1st.wet]<(5*xa.Dh.nor.1st[i]))
      notsofar.dry.i<-which(Disth.nor.1st[yoindx.nor.1st.dry]<(5*xa.Dh.nor.1st[i]))
      n.notsofar.wet.i<-length(notsofar.wet.i)
      n.notsofar.dry.i<-length(notsofar.dry.i)
      if (n.notsofar.wet.i==0) {
        x.CG.nor.dry[i]<-1
      } else if (n.notsofar.dry.i>0) {
# note: if (length(notsofar.dry.i)==0) next is implicit
        # IDI wet stations
        G.wet<-exp(-0.5*(Disth.nor.1st[yoindx.nor.1st.wet[notsofar.wet.i]]/xa.Dh.nor.1st[i])**2.
                   -0.5*(Distz.nor.1st[yoindx.nor.1st.wet[notsofar.wet.i]]/xa.Dz.nor.1st[i])**2.)
        if (n.notsofar.wet.i==1) {
          K.wet<-G.wet/(1+eps2)
        } else {
          D.wet<-as.matrix(exp(-0.5*(Disth[yoindx.wet[notsofar.wet.i],
                                           yoindx.wet[notsofar.wet.i]]/xa.Dh.nor.1st[i])**2.
                               -0.5*(Distz[yoindx.wet[notsofar.wet.i],
                                           yoindx.wet[notsofar.wet.i]]/xa.Dz.nor.1st[i])**2. ))
          D.wet[row(D.wet)==col(D.wet)]<-D.wet[row(D.wet)==col(D.wet)]+eps2
          InvD.wet<-solve(D.wet)
          K.wet<-tcrossprod(G.wet,InvD.wet)
          rm(D.wet,InvD.wet)
        }
        rm(G.wet)
        xidi.CG.wet[i]<-sum(K.wet)
        rm(K.wet)
        # IDI dry stations
        G.dry<-exp(-0.5*(Disth.nor.1st[yoindx.nor.1st.dry[notsofar.dry.i]]/xa.Dh.nor.1st[i])**2.-
                    0.5*(Distz.nor.1st[yoindx.nor.1st.dry[notsofar.dry.i]]/xa.Dz.nor.1st[i])**2.)
        if (n.notsofar.dry.i==1) {
          K.dry<-G.dry/(1+eps2)
        } else {
          D.dry<-as.matrix(exp(-0.5*(Disth[yoindx.dry[notsofar.dry.i],
                                           yoindx.dry[notsofar.dry.i]]/xa.Dh.nor.1st[i])**2.+
                               -0.5*(Distz[yoindx.dry[notsofar.dry.i],
                                           yoindx.dry[notsofar.dry.i]]/xa.Dz.nor.1st[i])**2. ))
          D.dry[row(D.dry)==col(D.dry)]<-D.dry[row(D.dry)==col(D.dry)]+eps2
          InvD.dry<-solve(D.dry)
          K.dry<-tcrossprod(G.dry,InvD.dry)
          rm(D.dry,InvD.dry)
        }
        rm(G.dry)
        xidi.CG.dry[i]<-sum(K.dry)
        rm(K.dry)
        if (xidi.CG.dry[i]>xidi.CG.wet[i]) x.CG.nor.dry[i]<-1
      } # END: compare IDIs
    } # END: cycle over xindx.nor.CG
  } # END: cycle over nor.1st
# DQC DQC DQC
# to revise: if a dry observation is (1) not included in a dry area and (2) is in Norway then
# flag it as presumably bad observation and repeat everything
  r.CG.nor.dry <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
                        ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
  r.CG.nor.dry[]<-NA
  r.CG.nor.dry[mask.CG]<-round(x.CG.nor.dry,0)
  stn.nor.dry<-extract(r.CG.nor.dry,cbind(VecX,VecY),na.rm=T)
  yo.ok.dry.fail<-!(yo.ok.pos.dry %in% which(stn.nor.dry==1 | is.na(stn.nor.dry)))
  if (any(yo.ok.dry.fail)) {
    ydqc.flag[yo.ok.pos.dry[which(yo.ok.dry.fail)]]<-400
    next
  }
#  print(paste(yo,stn.nor.dry,stn.nor.dry,ydqc.flag,"\n"))
# [] Contiguous Rain Area: First Guess on the grid.CG
# [] Refine eve identification (Identify Wet and Dry Areas within the eve first guess)
  print("++ Refine eve identification (Identify Wet and Dry Areas within the eve first guess)")
  xa.Dh.eve.1st<-vector(mode="numeric",length=Lgrid.CG)
  xa.Dz.eve.1st<-vector(mode="numeric",length=Lgrid.CG)
  x.CG.wet<-vector(mode="numeric",length=Lgrid.CG)
  lab.eve.CG<-vector(mode="numeric",length=Lgrid.CG)
  n.eve<-0
# x.CG.@@@ -> =0 gridpoint not classified as @@@; =1 grid-point classified as @@@
# lab.eve.CG -> =0 gridpoint "dry"; =-1 gridpoint "wet", waiting to labelled;
#             n>0 gridpoint "wet" belonging to eve n (only for isolated wet station)
  x.CG.wet[]<-0
  lab.eve.CG[]<-0
  for (eve.1st in 1:neve.1st.vec) { # START: cycle over eve.1st
#    print(paste(eve.1st,"eve.1st +++++++++++++++++++++++++++++++++"))
    # define vectors / variables / indexes 
    Leve.1st<-Leve.1st.vec[eve.1st]
    VecX.eve.1st<-vector(mode="numeric",length=Leve.1st)
    VecY.eve.1st<-vector(mode="numeric",length=Leve.1st)
    VecZ.eve.1st<-vector(mode="numeric",length=Leve.1st)
    VecS.eve.1st<-vector(mode="numeric",length=Leve.1st)
    yo.eve.1st<-vector(mode="numeric",length=Leve.1st)
    VecX.eve.1st<-VecX[eve.1st.vec[eve.1st,1:Leve.1st]]
    VecY.eve.1st<-VecY[eve.1st.vec[eve.1st,1:Leve.1st]]
    VecZ.eve.1st<-VecZ[eve.1st.vec[eve.1st,1:Leve.1st]]
    VecS.eve.1st<-VecS[eve.1st.vec[eve.1st,1:Leve.1st]]
    yo.eve.1st<-yo[eve.1st.vec[eve.1st,1:Leve.1st]]
    yoindx.eve.1st.wet<-which(yo.eve.1st>=rr.inf & !is.na(yo.eve.1st))
    yoindx.eve.1st.dry<-which(yo.eve.1st<rr.inf & !is.na(yo.eve.1st))
    yoindx.wet<-eve.1st.vec[eve.1st,yoindx.eve.1st.wet]
    yoindx.dry<-eve.1st.vec[eve.1st,yoindx.eve.1st.dry]
    n.eve.1st.wet<-length(yoindx.eve.1st.wet)
    n.eve.1st.dry<-length(yoindx.eve.1st.dry)
    # the next line should not happen...
    if (n.eve.1st.wet==0) next
    # case of an isolated wet station...
    if (n.eve.1st.dry==0) {
      x.CG.eve.wet[which(lab.eve.CG.1st==eve.1st)]<-1
      next
    }
    xindx.eve.CG<-which(lab.eve.CG.1st==eve.1st)
    Lgrid.eve.1st<-length(xindx.eve.CG)
    # more than one single wet-station
    # define Dz for the current eve (min allowed value is 500 m)
    xa.Dz.eve.1st[xindx.eve.CG]<-max(500,max(zgrid.CG[xindx.eve.CG]))
    for (i in xindx.eve.CG) {    #START: cycle over xindx.eve.CG
      Disth.eve.1st<-( (ygrid.CG[i]-VecY.eve.1st)**2+(xgrid.CG[i]-VecX.eve.1st)**2 )**0.5 / 1000.
      Distz.eve.1st<-abs(zgrid.CG[i]-VecZ.eve.1st)
      close2i.ord<-order(Disth.eve.1st,decreasing=F)
#      xa.Dh.eve.1st[i] <- max(20,Disth.eve.1st[close2i.ord[1]])
      xa.Dh.eve.1st[i] <- Disth.eve.1st[close2i.ord[1]]
      notsofar.wet.i<-which(Disth.eve.1st[yoindx.eve.1st.wet]<(5*xa.Dh.eve.1st[i]))
      notsofar.dry.i<-which(Disth.eve.1st[yoindx.eve.1st.dry]<(5*xa.Dh.eve.1st[i]))
      n.notsofar.wet.i<-length(notsofar.wet.i)
      n.notsofar.dry.i<-length(notsofar.dry.i)
      if (n.notsofar.dry.i==0) {
        x.CG.wet[i]<-1
      } else if (n.notsofar.wet.i>0) {
# note: if (length(notsofar.wet.i)==0) next is implicit
        # IDI wet stations
        G.wet<-exp(-0.5*(Disth.eve.1st[yoindx.eve.1st.wet[notsofar.wet.i]]/xa.Dh.eve.1st[i])**2.
                   -0.5*(Distz.eve.1st[yoindx.eve.1st.wet[notsofar.wet.i]]/xa.Dz.eve.1st[i])**2.)
        if (n.notsofar.wet.i==1) {
          K.wet<-G.wet/(1+eps2)
        } else {
          D.wet<-as.matrix(exp(-0.5*(Disth[yoindx.wet[notsofar.wet.i],
                                           yoindx.wet[notsofar.wet.i]]/xa.Dh.eve.1st[i])**2.
                               -0.5*(Distz[yoindx.wet[notsofar.wet.i],
                                           yoindx.wet[notsofar.wet.i]]/xa.Dz.eve.1st[i])**2. ))
          D.wet[row(D.wet)==col(D.wet)]<-D.wet[row(D.wet)==col(D.wet)]+eps2
          InvD.wet<-solve(D.wet)
          K.wet<-tcrossprod(G.wet,InvD.wet)
          rm(D.wet,InvD.wet)
        }
        rm(G.wet)
        xidi.CG.wet[i]<-sum(K.wet)
        rm(K.wet)
        # IDI dry stations
        G.dry<-exp(-0.5*(Disth.eve.1st[yoindx.eve.1st.dry[notsofar.dry.i]]/xa.Dh.eve.1st[i])**2.-
                    0.5*(Distz.eve.1st[yoindx.eve.1st.dry[notsofar.dry.i]]/xa.Dz.eve.1st[i])**2.)
        if (n.notsofar.dry.i==1) {
          K.dry<-G.dry/(1+eps2)
        } else {
          D.dry<-as.matrix(exp(-0.5*(Disth[yoindx.dry[notsofar.dry.i],
                                           yoindx.dry[notsofar.dry.i]]/xa.Dh.eve.1st[i])**2.+
                               -0.5*(Distz[yoindx.dry[notsofar.dry.i],
                                           yoindx.dry[notsofar.dry.i]]/xa.Dz.eve.1st[i])**2. ))
          D.dry[row(D.dry)==col(D.dry)]<-D.dry[row(D.dry)==col(D.dry)]+eps2
          InvD.dry<-solve(D.dry)
          K.dry<-tcrossprod(G.dry,InvD.dry)
          rm(D.dry,InvD.dry)
        }
        rm(G.dry)
        xidi.CG.dry[i]<-sum(K.dry)
        rm(K.dry)
        if (xidi.CG.wet[i]>=xidi.CG.dry[i]) x.CG.wet[i]<-1
      }
    } # END: cycle over xindx.eve.CG
  } # END: cycle over eve.1st
# DQC DQC DQC
# to revise: if a wet observation is not included in an event and it is inside Norway then
# flag it as presumably bad observation and repeat everything
  r.CG.eve.wet <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
                        ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
  r.CG.eve.wet[]<-NA
  r.CG.eve.wet[mask.CG]<-round(x.CG.wet,0)
  stn.eve.wet<-extract(r.CG.eve.wet,cbind(VecX,VecY),na.rm=T)
  yo.ok.wet.fail<-!(yo.ok.pos.wet %in% which(stn.eve.wet==1 | is.na(stn.eve.wet)))
  if (any(yo.ok.wet.fail)) {
    ydqc.flag[yo.ok.pos.wet[which(yo.ok.wet.fail)]]<-500
    next
  }
# [] eve Labelling (Identify final eve on the grid)
# lab.eve.CG -> =0 gridpoint "dry"; =-1 gridpoint "wet", waiting to labelled;
#             n>0 gridpoint "wet" belonging to eve n (only for isolated wet station)
  r.lab.eve.CG<-clump(r.CG.eve.wet,directions=8)
  f.lab<-freq(r.lab.eve.CG)
  aux<-which(!is.na(f.lab[,1]))
  f.lab.val<-f.lab[aux,1]
  f.lab.n<-f.lab[aux,2]
  l.CG<-extract(r.lab.eve.CG,1:ncell(r.lab.eve.CG))
  lab.eve.CG<-l.CG[mask.CG]
  rm(f.lab,aux)
# [] assign stations at eve
#  remarks: 1. each eve must include at least 1 observation; 
#   2. it is possible to have wet-observations not included in any of the eve 
#      (for example, in areas where there are a lot of observations and most of them are "dry")
  y.eve<-extract(r.lab.eve.CG,cbind(VecX,VecY),na.rm=T)
  eve.labels<-as.vector(na.omit(unique(y.eve[yo.ok.pos.wet])))
  n.eve<-length(eve.labels)
  eve.nostn.pos<-which(!(f.lab.val %in% eve.labels))
  if (length(eve.nostn.pos)>0) {
    r <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
               ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
    x.aux<-vector(mode="numeric",length=Lgrid.CG)
    for (i in eve.nostn.pos) {
      aux<-which(l.CG==f.lab.val[i])
      x.aux[]<-NA
      x.aux[aux]<-l.CG[aux]
      r[]<-NA
      r[mask.CG]<-round(x.aux,0)
      dist.min<-NA
      for (e in eve.labels) {
        dist<-cellStats(distanceFromPoints(r,cbind(VecX[which(y.eve==e)],VecY[which(y.eve==e)])),mean)
        if (is.na(dist.min)) {
          dist.min<-dist
          new.lab<-e
        } else {
          if (dist<dist.min) {
            dist.min<-dist
            new.lab<-e
          }
        }
      }
      new.lab.pos<-which(f.lab.val==new.lab)
      l.CG[aux]<-new.lab
      f.lab.val[i]<-NA
      f.lab.n[new.lab.pos]<-f.lab.n[new.lab.pos]+f.lab.n[i]
      f.lab.n[i]<-NA
      r.lab.eve.CG<-l.CG
      lab.eve.CG<-l.CG[mask.CG]
    }
    y.eve<-extract(r.lab.eve.CG,cbind(VecX,VecY),na.rm=T)
    eve.labels<-as.vector(na.omit(unique(y.eve[yo.ok.pos.wet])))
    n.eve<-length(eve.labels)
  }
  rm(l.CG)
  if (exists("r")) rm(r,x.aux,aux,dist)
# finally, we're interested only in events occurring in Norway
  lab.eve.FG<-vector(mode="numeric",length=Lgrid.FG)
  lab.eve.FG[]<-0
  r.lab.eve.FG<-resample(r.lab.eve.CG,r.orog.FG,method="ngb")
  aux<-extract(r.lab.eve.FG,1:ncell(r.lab.eve.FG))
  lab.eve.FG<-as.integer(aux[mask.FG])
  rm(aux)
  aux<-which(eve.labels %in% lab.eve.FG)
  if (length(aux)<n.eve) {
    tmp<-eve.labels
    eve.labels<-NA
    eve.labels<-tmp[aux]
    n.eve<-length(eve.labels)
  }
  r.lab.eve.FG[]<-NA
  r.lab.eve.FG[mask.FG]<-lab.eve.FG
# ANALYSIS ANALYSIS ANALYSIS ANALYSIS ANALYSIS ANALYSIS ANALYSIS ANALYSIS
# [] Analysis
  print("++ events Analysis")
# + data structures definition
# FG -> Fine Grid
  xa.FG<-vector(mode="numeric",length=Lgrid.FG)
  xb.FG<-vector(mode="numeric",length=Lgrid.FG)
  xidi.FG<-vector(mode="numeric",length=Lgrid.FG)
  xidi.eve.FG<-vector(mode="numeric",length=Lgrid.FG)
# CG -> Coarse Grid
  xa.CG<-vector(mode="numeric",length=Lgrid.CG)
# Y vectors (station locations)
  ya<-vector(mode="numeric",length=L.y.tot)
  yb<-vector(mode="numeric",length=L.y.tot)
  ybv.mat<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
  yidi.eve<-vector(mode="numeric",length=L.y.tot)
  yidiv.eve<-vector(mode="numeric",length=L.y.tot)
  yo.n<-vector(mode="numeric",length=L.y.tot)
  yav<-vector(mode="numeric",length=L.y.tot)
  ya.tmp<-vector(mode="numeric",length=L.y.tot)
  yav.tmp<-vector(mode="numeric",length=L.y.tot)
#  eps2.eve.smallest<-vector(mode="numeric",length=n.eve)
# events descriptive vectors
  n.y.eve<-vector(mode="numeric",length=n.eve)
  n.ya.eve<-vector(mode="numeric",length=n.eve)
  area.eve<-vector(mode="numeric",length=n.eve)
  volume.eve<-vector(mode="numeric",length=n.eve)
  meanrain.eve<-vector(mode="numeric",length=n.eve)
  maxrain.x.eve<-vector(mode="numeric",length=n.eve)
  maxrain.yo.eve<-vector(mode="numeric",length=n.eve)
  maxrain.ya.eve<-vector(mode="numeric",length=n.eve)
  maxrain.yav.eve<-vector(mode="numeric",length=n.eve)
  meanidi.x.eve<-vector(mode="numeric",length=n.eve)
  meanidi.y.eve<-vector(mode="numeric",length=n.eve)
  meanidiv.y.eve<-vector(mode="numeric",length=n.eve)
  meanidiv.y.eve.q50<-vector(mode="numeric",length=n.eve)
  meanidiv.y.eve.q75<-vector(mode="numeric",length=n.eve)
  # ell -> ellipse
  ell.locx.eve<-vector(mode="numeric",length=n.eve)
  ell.locy.eve<-vector(mode="numeric",length=n.eve)
  ell.smajor.eve<-vector(mode="numeric",length=n.eve)
  ell.sminor.eve<-vector(mode="numeric",length=n.eve)
  ell.smadir.eve<-vector(mode="numeric",length=n.eve)
  # cv -> cross-validated
  cv.rel.eve.all<-vector(mode="numeric",length=n.eve)
  cv.bias.eve.all<-vector(mode="numeric",length=n.eve)
  cv.rmse.eve.all<-vector(mode="numeric",length=n.eve)
  cv.made.eve.all<-vector(mode="numeric",length=n.eve)
  cv.rel.eve.q50<-vector(mode="numeric",length=n.eve)
  cv.bias.eve.q50<-vector(mode="numeric",length=n.eve)
  cv.rmse.eve.q50<-vector(mode="numeric",length=n.eve)
  cv.made.eve.q50<-vector(mode="numeric",length=n.eve)
  cv.rel.eve.q75<-vector(mode="numeric",length=n.eve)
  cv.bias.eve.q75<-vector(mode="numeric",length=n.eve)
  cv.rmse.eve.q75<-vector(mode="numeric",length=n.eve)
  cv.made.eve.q75<-vector(mode="numeric",length=n.eve)
  idi.norm.fac<-vector(mode="numeric",length=n.eve)
  n.q50<-vector(mode="numeric",length=n.eve)
  n.q75<-vector(mode="numeric",length=n.eve)
# + initialization
  xa.FG[]<-NA
  xb.FG[]<-NA
  xidi.FG[]<-0
  xidi.eve.FG[]<-0
  xa.CG[]<-NA
  yo.n[]<-NA
  yo.n[yo.ok.pos]<-0
  yb[]<-NA
  yb[yo.ok.pos]<-0
  yidi.eve[]<-NA
  yidi.eve[yo.ok.pos]<-0
  yidiv.eve[]<-NA
  yidiv.eve[yo.ok.pos]<-0
  ya[]<-NA
  ya[yo.ok.pos]<-0
  yav[]<-NA
  yav[yo.ok.pos]<-0
  ell.locx.eve[]<-NA
  ell.locy.eve[]<-NA
  ell.smajor.eve[]<-NA
  ell.sminor.eve[]<-NA
  ell.smadir.eve[]<-NA
# ELLIPSOID HULLS  ELLIPSOID HULLS  ELLIPSOID HULLS  ELLIPSOID HULLS  ELLIPSOID HULLS 
# + Ellipsoid hulls
# ellipsoid hulls computation could fail (i.e. small events) and return NAs. 
  # debug: start
#  png(file="../../seNorge2_scratch/Bspat_PREC1d/elli.png",width=1200,height=1200)
#  mx<-as.integer(max(eve.labels,na.rm=T))
#  cols<-c("gray",rainbow((mx-1)))
#  plot(r.lab.eve.FG,breaks=seq(0.5,(mx+0.5),by=1),col=cols)
  # debug: end
  for (n in 1:n.eve) {
#    print(paste("n eve.labels[n]",n,eve.labels[n]))
    xindx.eve.CG<-which(lab.eve.CG==eve.labels[n])
    Lgrid.eve<-length(xindx.eve.CG)
#    print(Lgrid.eve)
    xy<-cbind(xgrid.CG[xindx.eve.CG],ygrid.CG[xindx.eve.CG])
    ell<-ellipsoidhull(xy)
    ell.locx.eve[n]<-ell$loc[1]
    ell.locy.eve[n]<-ell$loc[2]
#   debug plot: lines(predict(ell), col="blue")
#    lines(predict(ell), col="blue")
#    points(ell.locx.eve[n],ell.locy.eve[n],pch=19,col="black",cex=1.5)
    if ( is.na(ell$cov[1,1]) | is.na(ell$cov[2,2]) | is.na(ell$cov[1,2]) | 
         ell$cov[1,1]==0 | ell$cov[2,2]==0 | ell$cov[1,2]==0) {
      ell.smajor.eve[n] <- 2*min.Dh.seq.allowed
      ell.sminor.eve[n] <- 2*min.Dh.seq.allowed
      ell.smadir.eve[n] <- NA
      next
    }
#    eigen -> values sorted in decreasing order; 
#             vectors matrix whose columns contain the eigenvectors
    eigenval<-eigen(ell$cov)$values
    eigenvec<-eigen(ell$cov)$vectors
    e <- sqrt(eigenval)
    ell.smajor.eve[n] <- sqrt(ell$d2) * e[1] /1000  # semi-major axis [Km]
    ell.sminor.eve[n] <- sqrt(ell$d2) * e[2] /1000  # semi-minor axis [Km]
    # dir=0 N-S orientation; dir=45 NE-SW; dir=90 E-W; dir=135 NW-SE; dir=180 N-S
    ell.smadir.eve[n]<-atan2(eigenvec[1,1],eigenvec[2,1])
    if (ell.smadir.eve[n]<0) ell.smadir.eve[n]<-ell.smadir.eve[n]+pi
    ell.smadir.eve[n]<-ell.smadir.eve[n]/pi*180.
  } # end cycle: compute ellipsoid hulls
# debug: close plot session
#  dev.off()
  rm(xindx.eve.CG,xy,ell,eigenval,eigenvec,e)
# ANALYSIS OVER EVENTS  ANALYSIS OVER EVENTS  ANALYSIS OVER EVENTS  ANALYSIS OVER EVENTS  
# + Analysis cycle over events 
  print("++ Analysis cycle over events")
  for (n in 1:n.eve) {
    # initialization
    r.xb.CG <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
                     ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
    r.xidi.CG <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
                       ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
    r.xb.FG <-raster(ncol=nx.FG, nrow=ny.FG, xmn=xmn.FG, xmx=xmx.FG,
                     ymn=ymn.FG, ymx=ymx.FG, crs=proj4.utm33)
    r.xidi.FG <-raster(ncol=nx.FG, nrow=ny.FG, xmn=xmn.FG, xmx=xmx.FG,
                       ymn=ymn.FG, ymx=ymx.FG, crs=proj4.utm33)
    r.aux.CG <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
                      ymn=ymn.CG, ymx=ymx.CG, crs=proj4.utm33)
    r.aux.FG <-raster(ncol=nx.FG, nrow=ny.FG, xmn=xmn.FG, xmx=xmx.FG,
                      ymn=ymn.FG, ymx=ymx.FG, crs=proj4.utm33)
    xb.CG<-vector(mode="numeric",length=Lgrid.CG)
    xidi.eve.CG<-vector(mode="numeric",length=Lgrid.CG)
    #
    r.xb.CG[]<-NA
    r.xb.FG[]<-NA
    r.xidi.CG[]<-NA
    r.xidi.FG[]<-NA
    r.aux.CG[]<-NA
    r.aux.FG[]<-NA
    xb.CG[]<-NA
    xidi.eve.CG[]<-0
    # identify event extension on the grid
    xindx.eve.CG<-which(lab.eve.CG==eve.labels[n])
    xindx.eve.FG<-which(lab.eve.FG==eve.labels[n])
    Lgrid.eve.CG<-length(xindx.eve.CG)
    Lgrid.eve.FG<-length(xindx.eve.FG)
    ya.indx<-which(y.eve==eve.labels[n] & !(yo.ok.wet))
    yindx<-which(y.eve==eve.labels[n] & yo.ok.wet)
    eve.rx<-range(xgrid.CG[xindx.eve.CG])
    eve.ry<-range(ygrid.CG[xindx.eve.CG])
    #
    n.y.eve[n]<-length(yindx)
    n.ya.eve[n]<-length(ya.indx)
    #
    print(paste("+ eve lab >",eve.labels[n],
                "(",n,"/",n.eve,")",
                " #obs=",n.y.eve[n]," #stns=",(n.y.eve[n]+n.ya.eve[n]),
                sep=""))
    print(paste("event range x / y=",eve.rx,"/",eve.ry))
#   ---------------------------------------------------------------------------
#   + event observed by a single (isolated) station
    if (n.y.eve[n]==1) {
      # Improve the identification of event area (from CG to FG) and smooth the borders
      # note: xidi=-1; yidi/v=-1; yb=NA; yav=0
      r.xb.CG[]<-NA
      r.xb.CG[mask.CG[xindx.eve.CG]]<-yo[yindx]
      r.xb.CG<-focal(r.xb.CG,w=matrix(1/9, nr=3, nc=3),na.rm=T)
      r.xb.CG[mask.CG[xindx.eve.CG]]<-yo[yindx]
      r.xb.CG<-trim(r.xb.CG)
      r.aux.FG<-crop(r.orog.FG,r.xb.CG)
      r.aux.FG[]<-NA
      r.xb.FG<-resample(r.xb.CG,r.aux.FG,method="bilinear")
      r.xb.FG<-extend(r.xb.FG,r.orog.FG)
      xb.FG<-extract(r.xb.FG,mask.FG)
      xindx.eve.FG<-which(!is.na(xb.FG))
      Lgrid.eve.FG<-length(xindx.eve.FG)
      lab.eve.FG[xindx.eve.FG]<-eve.labels[n]
#
##      r.aux.FG[mask.FG[xindx.eve.FG]]<-rep(1,length=Lgrid.eve.FG)
#      r.aux.FG<-trim(r.aux.FG)
#      r.xb.FG<-resample(r.xb.CG,r.aux.FG,method="bilinear")
#      r.xb.FG<-extend(r.xb.FG,r.orog.FG)
#      xb.FG<-extract(r.xb.FG,mask.FG)
#      xb.FG[is.na(xb.FG)]<-0
#      r.xb.FG[mask.FG]<-xb.FG
#      r.xb.FG<-focal(r.xb.FG,w=fg,na.rm=T)
#      xb.FG<-extract(r.xb.FG,mask.FG)
#      aux<-which(xb.FG>=rr.inf & !is.na(xb.FG))
#      lab.eve.FG[aux]<-eve.labels[n]
#      xindx.eve.FG<-which(lab.eve.FG==eve.labels[n])
#      Lgrid.eve.FG<-length(xindx.eve.FG)
#      rm(r.aux.FG)
      xa.FG[xindx.eve.FG]<-xb.FG[xindx.eve.FG]
      ya[yindx]<-yo[yindx]
      xb.FG[xindx.eve.FG]<-NA
      xidi.FG[xindx.eve.FG]<-(-1)
      yb[yindx]<-NA
      yav[yindx]<-0
      yidi.eve[yindx]<-(-1)
      yidiv.eve[yindx]<-(-1)
      ya[ya.indx]<-extract(r.xb.FG,cbind(VecX[ya.indx],VecY[ya.indx]),na.rm=T)
#      aux<-which(ya[ya.indx]<rr.inf & !is.na(ya[ya.indx]) )
#      ya[ya.indx][aux]<-rr.inf
      yb[ya.indx]<-NA
      yav[ya.indx]<-0
      yidi.eve[ya.indx]<-(-1)
      yidiv.eve[ya.indx]<-(-1)
      # event properties 
      area.eve[n]<-Lgrid.eve.FG*area.1cell.FG # Area Km**2
      volume.eve[n]<-sum(xa.FG[xindx.eve.FG])
      meanidi.x.eve[n]<-NA
      meanidi.y.eve[n]<-NA
      meanidiv.y.eve[n]<-NA
      meanidiv.y.eve.q50[n]<-NA
      meanidiv.y.eve.q75[n]<-NA
      meanrain.eve[n]<-mean(xa.FG[xindx.eve.FG])
      maxrain.x.eve[n]<-max(xa.FG[xindx.eve.FG])
      maxrain.yo.eve[n]<-max(yo[yindx])
      maxrain.ya.eve[n]<-max(ya[yindx])
      maxrain.yav.eve[n]<-NA
      cv.rel.eve.all[n]<-NA
      cv.bias.eve.all[n]<-NA
      cv.rmse.eve.all[n]<-NA
      cv.made.eve.all[n]<-NA
      cv.rel.eve.q50[n]<-NA
      cv.bias.eve.q50[n]<-NA
      cv.rmse.eve.q50[n]<-NA
      cv.made.eve.q50[n]<-NA
      cv.rel.eve.q75[n]<-NA
      cv.bias.eve.q75[n]<-NA
      cv.rmse.eve.q75[n]<-NA
      cv.made.eve.q75[n]<-NA
      idi.norm.fac[n]<-NA
      n.q50[n]<-NA
      n.q75[n]<-NA
      next
    }
#   ---------------------------------------------------------------------------
#   + event including more than one observation
#   + Larger-scale Background = mode(distribution of observations)
    histog<-hist(yo[yindx],breaks=seq(0,round((max(yo[yindx])+1),0)),plot=F)
    mode<-which.max(histog$counts)-0.5
    if (mode<rr.inf) mode<-rr.inf
    xb.CG[]<-NA
    xb.CG[xindx.eve.CG]<-mode
    yb[yindx]<-mode
    ybv.mat[yindx,yindx]<-mode
    if (n.ya.eve[n]>0) yb[ya.indx]<-mode
#   + define Dh.seq: sequence of horizontal decorrelation length scales
#                    (coarser to finer spatial scales)
#   general case:
#    from ellipsoid hull major semiaxis (coarser scale)
#    to a value proportional to the minimum distance between two stations (finer scale)
#   details: particular cases (i.e. semiaxis is NA, semiaxis is smaller then min distance);
    Disth.aux<-Disth[yindx,yindx]
    # minimum distance between two eve stations
    min.Disth<-min(Disth.aux[row(Disth.aux)!=col(Disth.aux)])/4.
#    min.Dh.seq<-as.integer(max(min.Dh.seq.allowed,ceiling(min.Disth)))
    min.Dh.seq<-as.integer(min.Dh.seq.allowed)
    max.Dh.seq<-as.integer(max(5*min.Dh.seq.allowed,ceiling(ell.smajor.eve[n])))
    aux<-which(Dh.seq.reference>=min.Dh.seq & Dh.seq.reference<=max.Dh.seq)
    Dh.seq<-c(Dh.seq.reference[aux],min.Dh.seq.allowed)
    n.Dh.seq<-length(Dh.seq)
    print(paste("Dh.seq [Km] ->",Dh.seq))
#   + pre-elaboration to obtain err-cov matrices on CG
    G.CG<-matrix(ncol=n.y.eve[n],nrow=Lgrid.eve.CG,data=0.)
    aux.CG<-(outer(ygrid.CG[xindx.eve.CG],VecY[yindx],FUN="-")**2. +
             outer(xgrid.CG[xindx.eve.CG],VecX[yindx],FUN="-")**2.)**0.5/1000.
    auxz.CG<-abs(outer(zgrid.CG[xindx.eve.CG],VecZ[yindx],FUN="-"))
    InvD<-matrix(data=0,ncol=n.y.eve[n],nrow=n.y.eve[n])
#   + cycle on horizontal decorrelation length scales
    flag.firsttime<-T
    n.iter<-0
    for (Dh.test in Dh.seq) {
      n.iter<-n.iter+1
      D.test.part<-exp(-0.5*(Disth[yindx,yindx]/Dh.test)**2.)
      J<-1000000
      if (Dh.test>=300) Dz.seq<-c(10000)
      if (Dh.test>50 & Dh.test<300) Dz.seq<-Dz.seq.gt50
      if (Dh.test<=50) Dz.seq<-Dz.seq.le50
      if (n.iter==n.Dh.seq) {
          yo.n.tmp<-vector(mode="numeric",length=L.y.tot)
          yo.n.tmp[yindx]<-yo[yindx]
      }
#     + cycle on vertical decorrelation length scales
      for (Dz.test in Dz.seq) {
        # compute superobservations according to decorrelation lenght-scales
        if (n.iter!=n.Dh.seq) {
          yo.n.tmp<-vector(mode="numeric",length=L.y.tot)
          for (b in yindx) {
            aux.yo<-which(Disth[b,yindx]<=(Dh.test) & Distz[b,yindx]<=(Dz.test))
            yo.n.tmp[b]<-mean(yo[yindx[aux.yo]])
          }
        }
        # analysis at eve station locations
        D.test<-D.test.part*exp(-0.5*(Distz[yindx,yindx]/Dz.test)**2.)
        S.test<-D.test
        D.test[row(D.test)==col(D.test)]<-D.test[row(D.test)==col(D.test)]+eps2
        InvD.test<-solve(D.test)
        W.test<-tcrossprod(S.test,InvD.test)
        aux<-t(yo.n.tmp[yindx]-yb[yindx])
        ya.tmp[yindx]<-tcrossprod(W.test,aux)
        ya.tmp[yindx]<-yb[yindx] + ya.tmp[yindx]
        rm(aux)
        yav.tmp[yindx]<-yo.n.tmp[yindx] +
                           1./(1.-diag(W.test)) * (ya.tmp[yindx]-yo.n.tmp[yindx])
        aux<-which(ya.tmp[yindx]<rr.inf & !is.na(ya.tmp[yindx]))
        ya.tmp[yindx][aux]<-rr.inf
        yav.tmp[yindx][aux]<-rr.inf
        # cost function
        J.tmp<-(1/n.y.eve[n]*sum((st.log(yav.tmp[yindx])-st.log(yo.n.tmp[yindx]))**2))**0.5
        # J.tmp<-(1/n.y.eve[n]*sum((yav.tmp[yindx]-yo.n.tmp[yindx])**2))**0.5
        if (Dz.test==Dz.seq[1] | J.tmp<J) {
          # save the best values (minimum value of the cost function)
          J<-J.tmp
          InvD<-InvD.test
          D<-D.test
          S<-S.test
          W<-W.test
          Dz.choice<-Dz.test
          ya[yindx]<-ya.tmp[yindx]
          yav[yindx]<-yav.tmp[yindx]
          yo.n[yindx]<-yo.n.tmp[yindx]
          yidi.eve[yindx]<-yidi.eve[yindx]+rowSums(W)
          yidiv.eve[yindx]<-yidiv.eve[yindx]+
                            rep(1,n.y.eve[n])+ 1./(1.-diag(W)) * (rowSums(W)-rep(1,n.y.eve[n]))
        }
      } # end cycle on vertical decorrelation lenght scales
      print(paste("Dh=",Dh.test,"Km ==> Dz=",Dz.choice,"m"))
      rm(D.test.part,D.test,S.test,InvD.test,W.test)
      t.d<-t(yo.n[yindx]-yb[yindx])
      if ( !(any(abs(t.d)>0.01)) & n.iter!=n.Dh.seq) next
#     + CrossValidated-analysis 
      i<-0
      for (b in yindx) {
        i<-i+1
        # InvD.1<-solve(D[-i,-i])
        InvD.1<-InvD[-i,-i]-1/InvD[i,i]*outer(InvD[-i,i],InvD[i,-i])
        # W.1<-S[,-i]%*%InvD.1
        # ya.tmp[yindx]<-ybv.mat[yindx,b] + W.1 %*% (yo.n[yindx[-i]]-ybv.mat[yindx[-i],b])
        W.1<-tcrossprod(S[,-i],t(InvD.1))
        aux<-t(yo.n[yindx[-i]]-ybv.mat[yindx[-i],b])
        ya.tmp[yindx]<-tcrossprod(W.1, aux )
        ya.tmp[yindx]<-ybv.mat[yindx,b] + ya.tmp[yindx]
        ybv.mat[yindx,b]<-ya.tmp[yindx]
      }
      rm(W.1,InvD.1,D,W,S,aux)
#     + Analysis on other station points within the event areas
      if (n.ya.eve[n]>0) {
        G.ya<-matrix(ncol=n.y.eve[n],nrow=n.ya.eve[n],data=0.)
        G.ya<-exp(-0.5*(Disth[ya.indx,yindx]/Dh.test)**2.
                  -0.5*(Distz[ya.indx,yindx]/Dz.choice)**2.)
        K<-tcrossprod(G.ya,InvD)
        yidi.eve[ya.indx]<-yidi.eve[ya.indx]+rowSums(K)
        yidiv.eve[ya.indx]<-yidi.eve[ya.indx]
        ya[ya.indx]<-tcrossprod(K,t.d)
        ya[ya.indx]<-ya[ya.indx] + yb[ya.indx]
        aux<-which(ya[ya.indx]<rr.inf & !is.na(ya[ya.indx]))
        ya[ya.indx][aux]<-rr.inf
        yav[ya.indx]<-ya[ya.indx]
        rm(K,G.ya)
      }
#     + Dh>ref.value ->analysis on CG; Dh<ref.value ->analysis on FG
#     + Analysis on CG
      if (Dh.test>Dh.seq.ref) { 
        G.CG<-exp(-0.5*(aux.CG/Dh.test)**2.-0.5*(auxz.CG/Dz.choice)**2.)
        K.CG<-tcrossprod(G.CG,InvD)
        rm(G.CG)
        # update analysis
        xa.CG[xindx.eve.CG]<-xb.CG[xindx.eve.CG]+tcrossprod(K.CG,t.d)
        aux<-which(xa.CG[xindx.eve.CG]<rr.inf & !is.na(xa.CG[xindx.eve.CG]))
        xa.CG[xindx.eve.CG][aux]<-rr.inf
        # update IDI (which is cumulative)
        xidi.eve.CG[xindx.eve.CG]<-xidi.eve.CG[xindx.eve.CG]+rowSums(K.CG)
        rm(K.CG)
        # next background is the current analysis, only a little bit smoothed
        # smoothing is for realistic border effects
        xb.CG[]<-NA
        xb.CG[xindx.eve.CG]<-xa.CG[xindx.eve.CG]
      } else { 
#     + Analysis on FG
        # get background from CG, if it is needed, otherwise the background is the "mode"
        # note: if Dh.test==Dh.seq[1] no CG background is present
        if (flag.firsttime) {
          flag.firsttime<-F
          #
          r.xidi.CG[]<-NA
          r.xidi.CG[mask.CG[xindx.eve.CG]]<-xidi.eve.CG[xindx.eve.CG]
          r.xidi.CG<-focal(r.xidi.CG,w=matrix(1/9, nr=3, nc=3),na.rm=T)
          r.xidi.CG[mask.CG[xindx.eve.CG]]<-xidi.eve.CG[xindx.eve.CG]
          r.xidi.CG<-trim(r.xidi.CG)
          #          
          r.xb.CG[]<-NA
          r.xb.CG[mask.CG[xindx.eve.CG]]<-xb.CG[xindx.eve.CG]
          r.xb.CG<-focal(r.xb.CG,w=matrix(1/9, nr=3, nc=3),na.rm=T)
          r.xb.CG[mask.CG[xindx.eve.CG]]<-xb.CG[xindx.eve.CG]
          r.xb.CG<-trim(r.xb.CG)
          #
          r.aux.FG<-crop(r.orog.FG,r.xb.CG)
          r.aux.FG[]<-NA
          r.xb.FG<-resample(r.xb.CG,r.aux.FG,method="bilinear")
          r.xb.FG<-extend(r.xb.FG,r.orog.FG)
          xb.FG<-extract(r.xb.FG,mask.FG)
#          xindx.eve.FG<-which(!is.na(xb.FG) & xb.FG>=rr.inf)
          xindx.eve.FG<-which(!is.na(xb.FG))
          Lgrid.eve.FG<-length(xindx.eve.FG)
          lab.eve.FG[xindx.eve.FG]<-eve.labels[n]
          #
          r.aux.FG[]<-NA
          r.xidi.FG<-resample(r.xidi.CG,r.aux.FG,method="bilinear")
          r.xidi.FG<-extend(r.xidi.FG,r.orog.FG)
          xidi.eve.FG<-extract(r.xidi.FG,mask.FG)
          xidi.FG[xindx.eve.FG]<-xidi.eve.FG[xindx.eve.FG]
          #
#          rm(r.aux.FG)
#
#          xindx.eve.CG<-which(!is.na(xb.CG))
#          Lgrid.eve.CG<-length(xindx.eve.CG)
#
#          r.aux.FG[]<-NA
#          r.aux.FG[mask.FG[xindx.eve.FG]]<-rep(1,length=Lgrid.eve.FG)
#          r.aux.FG<-trim(r.aux.FG)
#
#          r.xidi.FG<-resample(r.xidi.CG,r.aux.FG,method="bilinear")
#          r.xidi.FG<-extend(r.xidi.FG,r.orog.FG)
#          xidi.eve.FG<-extract(r.xidi.FG,mask.FG)
#          xidi.FG[xindx.eve.FG]<-xidi.eve.FG[xindx.eve.FG]
#          r.xb.FG<-resample(r.xb.CG,r.aux.FG,method="bilinear")
#          r.xb.FG<-extend(r.xb.FG,r.orog.FG)
#          xb.FG<-extract(r.xb.FG,mask.FG)
#          # refine event area
#          xb.FG[is.na(xb.FG)]<-0
#          r.xb.FG[mask.FG]<-xb.FG
#          r.xb.FG<-focal(r.xb.FG,w=fg,na.rm=T)
#          xb.FG<-extract(r.xb.FG,mask.FG)
#          aux<-which(is.na(xb.FG[xindx.eve.FG]))
#          if (length(aux>0)) xb.FG[xindx.eve.FG][aux]<-rr.inf
#          aux<-which(xb.FG[xindx.eve.FG]<rr.inf)
#          if (length(aux>0)) xb.FG[xindx.eve.FG][aux]<-rr.inf
#          aux<-which(xb.FG>=rr.inf & !is.na(xb.FG))
#          lab.eve.FG[aux]<-eve.labels[n]
#          xindx.eve.FG<-which(lab.eve.FG==eve.labels[n])
#          Lgrid.eve.FG<-length(xindx.eve.FG)
        }
        # analysis on FG is done iteratively (save memory)
        i<-0
        while ((i*ndim.FG.iteration)<Lgrid.eve.FG) {
          start<-i*ndim.FG.iteration+1
          end<-min( (i+1)*ndim.FG.iteration, Lgrid.eve.FG)
          ndimaux.FG<-end-start+1
          # print(paste("start end ndimaux.FG",start,end,ndimaux.FG))
          pos<-xindx.eve.FG[start:end]
          aux.FG<-matrix(ncol=n.y.eve[n],nrow=ndimaux.FG,data=0.)
          auxz.FG<-matrix(ncol=n.y.eve[n],nrow=ndimaux.FG,data=0.)
          aux.FG<-(outer(ygrid[pos],VecY[yindx],FUN="-")**2. +
                outer(xgrid[pos],VecX[yindx],FUN="-")**2.)**0.5/1000.
          auxz.FG<-abs(outer(zgrid[pos],VecZ[yindx],FUN="-"))
          G.FG<-matrix(ncol=n.y.eve[n],nrow=ndimaux.FG,data=0.)
          G.FG<-exp(-0.5*(aux.FG/Dh.test)**2.-0.5*(auxz.FG/Dz.choice)**2.)
          K.FG<-tcrossprod(G.FG,InvD)
          rm(G.FG)
          xa.FG[pos]<-tcrossprod(K.FG,t.d)
          xa.FG[pos]<-xa.FG[pos]+xb.FG[pos]
          xidi.FG[pos]<-xidi.FG[pos]+rowSums(K.FG)
          rm(aux.FG,auxz.FG,K.FG)
          i<-i+1
        }
        aux<-which(xa.FG[xindx.eve.FG]<rr.inf & !is.na(xa.FG[xindx.eve.FG]))
        if (length(aux)>0) xa.FG[xindx.eve.FG][aux]<-rr.inf
        # current analysis is next iteration background
        xb.FG[]<-NA
        xb.FG[xindx.eve.FG]<-xa.FG[xindx.eve.FG]
      } # end: analysis on FG
      # current analysis is next iteration background
      yb[yindx]<-ya[yindx]
    } # end cycle on horizontal decorrelation lenght scales
    if (exists("K.CG")) rm(K.CG,G.CG)
    #
    idi.norm.fac[n]<-max(xidi.FG[xindx.eve.FG],yidi.eve[yindx],yidiv.eve[yindx])
    xidi.FG[xindx.eve.FG]<-xidi.FG[xindx.eve.FG]/idi.norm.fac[n]
    yidi.eve[yindx]<-yidi.eve[yindx]/idi.norm.fac[n]
    yidiv.eve[yindx]<-yidiv.eve[yindx]/idi.norm.fac[n]
    area.eve[n]<-Lgrid.eve.FG*area.1cell.FG # Area Km**2
    volume.eve[n]<-sum(xa.FG[xindx.eve.FG])
    meanidi.x.eve[n]<-mean(xidi.FG[xindx.eve.FG])
    meanidi.y.eve[n]<-mean(yidi.eve[yindx])
    meanidiv.y.eve[n]<-mean(yidiv.eve[yindx])
    meanrain.eve[n]<-mean(xa.FG[xindx.eve.FG])
    maxrain.x.eve[n]<-max(xa.FG[xindx.eve.FG])
    yo.aux<-yo[yindx]
    ya.aux<-ya[yindx]
    yav.aux<-yav[yindx]
    yidiv.aux<-yidiv.eve[yindx]
    maxrain.yo.eve[n]<-max(yo.aux)
    maxrain.ya.eve[n]<-max(ya.aux)
    maxrain.yav.eve[n]<-max(yav.aux)
    cv.rel.eve.all[n]<-mean(yav.aux/yo.aux)
    aux<-yav.aux-yo.aux
    cv.bias.eve.all[n]<-mean(aux)
    cv.rmse.eve.all[n]<-sqrt(mean(aux**2.))
    cv.made.eve.all[n]<-1.4826*median(abs(aux))
    q50.pos<-which(yo.aux>=q50.daily)
    q75.pos<-which(yo.aux>=q75.daily)
    n.q50[n]<-length(q50.pos)
    n.q75[n]<-length(q75.pos)
    if (n.q50[n]>0) {
      cv.rel.eve.q50[n]<-mean(yav.aux[q50.pos]/yo.aux[q50.pos])
      cv.bias.eve.q50[n]<-mean(aux[q50.pos])
      cv.rmse.eve.q50[n]<-sqrt(mean(aux[q50.pos]**2.))
      cv.made.eve.q50[n]<-1.4826*median(abs(aux[q50.pos]))
      meanidiv.y.eve.q50[n]<-mean(yidiv.aux[q50.pos])
    } else {
      cv.rel.eve.q50[n]<-NA
      cv.bias.eve.q50[n]<-NA
      cv.rmse.eve.q50[n]<-NA
      cv.made.eve.q50[n]<-NA
      meanidiv.y.eve.q50[n]<-NA
    }
    if (n.q75[n]>0) {
      cv.rel.eve.q75[n]<-mean(yav.aux[q75.pos]/yo.aux[q75.pos])
      cv.bias.eve.q75[n]<-mean(aux[q75.pos])
      cv.rmse.eve.q75[n]<-sqrt(mean(aux[q75.pos]**2.))
      cv.made.eve.q75[n]<-1.4826*median(abs(aux[q75.pos]))
      meanidiv.y.eve.q75[n]<-mean(yidiv.aux[q75.pos])
    } else {
      cv.rel.eve.q75[n]<-NA
      cv.bias.eve.q75[n]<-NA
      cv.rmse.eve.q75[n]<-NA
      cv.made.eve.q75[n]<-NA
      meanidiv.y.eve.q75[n]<-NA
    }
  } # END CYCLE over events
  if (exists("r.xb.CG")) rm(r.xidi.CG,r.xb.CG,xidi.eve.CG,xa.CG,xindx.eve.CG,r.aux.CG,xb.CG)
  if (exists("r.xb.FG")) rm(r.xidi.FG,r.xb.FG,r.aux.FG)
  break
} # end of DQC loop
#
xa.FG[is.na(xa.FG)]<-0.
#
r.aux.FG <-raster(ncol=nx.FG, nrow=ny.FG, xmn=xmn.FG, xmx=xmx.FG,
                  ymn=ymn.FG, ymx=ymx.FG, crs=proj4.utm33)
r.aux.FG[]<-NA
r.aux.FG[mask.FG]<-xa.FG
#r.aux.FG<-trim(r.aux.FG)
ya.tmp<-extract(r.aux.FG,cbind(VecX,VecY))
y.eve[is.na(ya)]<-NA
ya[is.na(ya)]<-round(ya.tmp[is.na(ya)],1)
yb[is.na(y.eve)]<-NA
yav[is.na(y.eve)]<-NA
yidi.eve[is.na(y.eve)]<-NA
yidiv.eve[is.na(y.eve)]<-NA
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
print("++ Output")
# Station Points - Write output on file 
cat(paste(yyyy.b,mm.b,dd.b,nday,
          formatC(VecS[stn.NO],format="f",digits=0),
          formatC(VecX[stn.NO],format="f",digits=0),
          formatC(VecY[stn.NO],format="f",digits=0),
          formatC(VecZ[stn.NO],format="f",digits=0),
          formatC(y.eve[stn.NO],format="f",digits=0),
          formatC(yo[stn.NO],format="f",digits=1),
          formatC(yb[stn.NO],format="f",digits=2),
          formatC(ya[stn.NO],format="f",digits=2),
          formatC(yav[stn.NO],format="f",digits=2),
          formatC(100*yidi.eve[stn.NO],format="f",digits=2),
          formatC(100*yidiv.eve[stn.NO],format="f",digits=2),
          formatC(ydqc.flag[stn.NO],format="f",digits=0),
          "\n",sep=";"),file=out.file.stn,append=T)
cat(paste(yyyy.b,mm.b,dd.b,nday,
          eve.labels,
          n.y.eve,
          formatC(area.eve,format="f",digits=0),
          formatC(volume.eve,format="f",digits=2),
          formatC(100*meanidi.x.eve,format="f",digits=2),
          formatC(100*meanidi.y.eve,format="f",digits=2),
          formatC(100*meanidiv.y.eve,format="f",digits=2),
          formatC(meanrain.eve,format="f",digits=2),
          formatC(maxrain.x.eve,format="f",digits=2),
          formatC(maxrain.yo.eve,format="f",digits=2),
          formatC(maxrain.ya.eve,format="f",digits=2),
          formatC(maxrain.yav.eve,format="f",digits=2),
          formatC(ell.locx.eve,format="f",digits=1),
          formatC(ell.locy.eve,format="f",digits=1),
          formatC(ell.smajor.eve,format="f",digits=1),
          formatC(ell.sminor.eve,format="f",digits=1),
          formatC(ell.smadir.eve,format="f",digits=1),
          formatC(cv.rel.eve.all,format="f",digits=4),
          formatC(cv.bias.eve.all,format="f",digits=4),
          formatC(cv.rmse.eve.all,format="f",digits=4),
          formatC(cv.made.eve.all,format="f",digits=4),
          formatC(cv.rel.eve.q50,format="f",digits=4),
          formatC(cv.bias.eve.q50,format="f",digits=4),
          formatC(cv.rmse.eve.q50,format="f",digits=4),
          formatC(cv.made.eve.q50,format="f",digits=4),
          formatC(100*meanidiv.y.eve.q50,format="f",digits=2),
          formatC(n.q50,format="f",digits=0),
          formatC(cv.rel.eve.q75,format="f",digits=4),
          formatC(cv.bias.eve.q75,format="f",digits=4),
          formatC(cv.rmse.eve.q75,format="f",digits=4),
          formatC(cv.made.eve.q75,format="f",digits=4),
          formatC(100*meanidiv.y.eve.q75,format="f",digits=2),
          formatC(n.q75,format="f",digits=0),
          formatC(idi.norm.fac,format="f",digits=5),
          "\n",sep=";"),file=out.file.eve,append=T)
# Figure: eve
r.aux.FG[]<-NA
r.aux.FG[mask.FG]<-round(100*xidi.FG,1)
rnc <- writeRaster(r.aux.FG,
                   filename=out.file.grd.idi,
                   format="CDF", overwrite=TRUE)
# Figure: Analysis on high-resolution grid
r.aux.FG[mask.FG]<-round(xa.FG,1)
if (flag.write.xa) {
  nogrid.ncout(grid=t(as.matrix(r.aux.FG)),
               x=x.FG,y=y.FG,grid.type=grid.type,
               file.name=out.file.grd.ana,
               var.name=var.name.xa,
               var.longname=var.longname.xa,
               var.unit=var.unit.xa,
               var.mv=var.mv.xa,
               var.version=var.version.xa,
               times=c(paste(yyyymmdd.b,"0000",sep="")),times.unit=times.unit.xa,
               times.ref=times.ref.xa,
               prod.date=prod.date,
               reference=reference.xa,
               proj4.string="+proj=utm +zone=33 +ellps=WGS84",
               source.string=source.nc)
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Success: exit!
quit(status=0)
