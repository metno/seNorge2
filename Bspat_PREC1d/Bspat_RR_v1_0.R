# --~- Bspat_RR_v1_0.R -~--
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
# High-Resolution (final) 1Km grid
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
# intermediate 5Km (coarse) grid (grid.CG)
# aggregate, upscale orog field with a factor of "fact"(5) (faster computations)
# orog.CG[1,1]<-sum(orog[i,j],i=1:fact,j=1:fact)/(fact**2)
# orog.CG[1,2]<-sum(orog[i,j],i=1:fact,j=(fact+1):(2*fact))/(fact**2)
# orog.CG[2,1]<-sum(orog[i,j],i=(fact+1):(2*fact),j=1:fact)/(fact**2)
# > print(orog.CG)
# class       : RasterLayer 
# dimensions  : 310, 239, 74090  (nrow, ncol, ncell)
# resolution  : 5000, 5000  (x, y)
# extent      : -75000, 1120000, 6450000, 8e+06  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 
# data source : in memory
# names       : fenno_dem_u33 
# values      : 0, 1984.88  (min, max)
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
# 1    2**0   missing
# 2    2**1   good 
# 3    2**2   err: 1h plausibility check failed 
# 4    2**3   not enough available observation 
# 5    2**4   err: external DQC 
# 6    2**5   masked station location (not enough geographical info) 
# 7    2**6   blacklist 1 
# 8    2**7   blacklist 2 
#
# Abbreviations:
# eve-> Event (described as a contiguous Observed Area of Precipitation, OAP)
# FG -> Higher Resolution grid (fine grid)
# CG -> Lower Resolution grid (coarse grid)
# q50 -> quantile(50) estimate for the daily precip in Norway (maded using 2013.09->2014.09 data)
# q75 -> quantile(50) estimate for the daily precip in Norway (maded using 2013.09->2014.09 data)
#
# History:
# 03.12.2014 - Cristian Lussana. Original code.
# -----------------------------------------------------------------------------
rm(list=ls())
# Libraries
library(raster)
library(rgdal)
library(ncdf)
require(tripack)
require(cluster)
# paths
#[DEVELOPMENT]
main.path.prog<-"/home/cristianl/projects/Bspat"
main.path.output<-"/home/cristianl"
main.path.geoinfo<-"/home/cristianl/geoinfo"
# common libs and etcetera
path2lib.com<-paste(main.path.prog,"/lib",sep="")
path2etc.com<-paste(main.path.prog,"/etc",sep="")
# libs and etcetera
path2lib.loc<-paste(main.path.prog, "/Bspat_RR/lib",sep="")
path2etc.loc<-paste(main.path.prog, "/Bspat_RR/etc",sep="")
#
path2output.main<-paste(main.path.output,"/seNorge2/PREC_daily",sep="")
path2output.main.stn<-paste(path2output.main,"/station_dataset",sep="")
path2output.main.grd<-paste(path2output.main,"/gridded_dataset",sep="")
path2output.add<-paste(main.path.output,"/seNorge2_addInfo/PREC_daily",sep="")
path2output.add.grd<-paste(path2output.add,"/gridded_dataset",sep="")
path2output.add.eve<-paste(path2output.add,"/event_dataset",sep="")
# External Functions
source(paste(path2lib.com,"/nogrid.ncout.R",sep=""))
source(paste(path2lib.com,"/ncout.spec.list.r",sep=""))
# Read Geographical Information
filenamedem<-paste(main.path.geoinfo,"/fenno_dem_u33.asc",sep="")
filenamemask<-paste(main.path.geoinfo,"/maskgrid.bil",sep="")
# Aux files
file_CG2FG<-paste(path2etc.loc,"/CG2FG.bin",sep="")
# Blacklist files
#file_errobs<-paste(path2etc.loc,"/SpInt_PREC_daily_ERRONEOUSOBSERVATIONS.txt",sep="")
#file_blacklist_current<-paste(path2etc.loc,"/SpInt_PREC_daily_BLACKLIST_current.txt",sep="")
#file_blacklist_never<-paste(path2etc.loc,"/SpInt_PREC_daily_BLACKLIST_never.txt",sep="")
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
#-------------------------------------------------------------------
# [] Setup parameters
# eps2 meaning: believe in observations
eps2<-0.1
# rain/norain threshold
rr.inf<-0.1
# upscale 1Km orog field with a factor of "fact" (faster computations)  
fact<-5
#
Lsubsample.max<-20
Lsubsample.DHmax<-150
# identify adjacent nodes
n.sector<-16
sector.angle<-360/n.sector
# 
superobs.radius<-30 #Km
# DQC thresholds 
yo.dqc.plausible.min<-0 # mm/h
yo.dqc.plausible.max<-500 # mm/h
# evaluation thresholds
#q50.hourly<-0.4 #mm/h (weak/normal)
#q75.hourly<-1.0 #mm/h (normal/moderate-strong)
q50.daily<-3 #mm/daily (weak/normal)
q75.daily<-8 #mm/daily (normal/moderate-strong)
# 
n.Dh.seq<-10
min.Dh.seq.allowed<-10 # Km
Dh.seq.ref<-100 # Km
#Dh.seq.reference<-c(5000,4000,3000,2500,2000,1750,1500,1400,1300,1200,
#                    1100,1000, 900, 800, 700, 650, 600, 550, 500, 450,
#                     400, 350, 300, 260, 230, 200, 180, 160, 140, 120,
#                     100,  90,  80,  70,  60,  50,  40,  30,  20,  10)
Dh.seq.reference<-c(5000,4000,3000,2000,1000, 900, 800, 700, 600, 500,
                     400, 300, 200, 150, 100,  80,  60,  40,  20,  10)
Dz.seq.gt50<-c(5000,2000,1000,500) # m
Dz.seq.le50<-c(1000,500) # m
#
ndim.FG.iteration<-10000
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
# MAIN ========================================================================
# [] Read command line arguments
arguments <- commandArgs()
arguments
date.b.string<-arguments[3]
date.e.string<-arguments[4]
file_blacklist_current<-arguments[5]
file_blacklist_never<-arguments[6]
file_errobs<-arguments[7]
if (length(arguments)!=7) {
  print("Error in command line arguments:")
  print("R --vanilla yyyy.mm.dd yyyy.mm.dd blacklist_current blacklist_never errobs")
  print("             begin -----> end of the accumulation period")
  quit(status=1)
}
# mode.flag=0 -> operational mode
#
# date.x.string<-"YYYY.MM.DD"
#                 1234567890
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
monthseq<-as.POSIXlt(seq(as.POSIXlt(start.day),as.POSIXlt(end.day),by="1 month"),"UTC")
nmonth<-length(monthseq)
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
                     "/seNorge_v2test_PREC_station_",
                     yyyymmdd.b,"_",yyyymmdd.e,".txt",sep="")
out.file.eve<- paste(path2output.add.eve,"/",yyyymm.b,
                     "/seNorge_v2test_PREC_event_",
                     yyyymmdd.b,"_",yyyymmdd.e,".txt",sep="")
out.file.grd.ana<- paste(path2output.main.grd,"/",yyyymm.b,
                         "/seNorge_v2test_PREC_grid_",
                         yyyymmdd.b,"_",yyyymmdd.e,".nc",sep="")
out.file.grd.idi<- paste(path2output.add.grd,"/",yyyymm.b,
                         "/seNorge_v2test_PREC_grid_normIDI_",
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
#------------------------------------------------------------------------------
# [] Grid - same as the DEM grid 
# CRS Coordinate Reference System
projstr<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
proj_format_str<-"PROJ.4"
stackGeoGrid<-raster(filenamedem)
mask.FG<-raster(filenamemask)
projection(stackGeoGrid)<-projstr
projection(mask.FG)<-projstr
nx.FG<-ncol(stackGeoGrid)
ny.FG<-nrow(stackGeoGrid)
dx.FG<-xres(stackGeoGrid)
dy.FG<-yres(stackGeoGrid)
area.1cell.FG<-dx.FG*dy.FG/(1000*1000) # Area Km**2
# 4 borders point (SW corner (xmn,ymn); NE corner (xmx,ymx))
xmn.FG<-xmin(stackGeoGrid)
xmx.FG<-xmax(stackGeoGrid)
ymn.FG<-ymin(stackGeoGrid)
ymx.FG<-ymax(stackGeoGrid)
# Extract orography on unmasked gridpoints only
# orog has the same dimensions of stackGeoGrid but it is masked (NAs) outside Norway
orog.FG<-mask(stackGeoGrid,mask.FG,maskvalue=0)
# LowRes CG orography
orog.CG<-aggregate(orog.FG,fact=fact,fun=mean,na.rm=T)
nx.CG<-ncol(orog.CG)
ny.CG<-nrow(orog.CG)
dx.CG<-xres(orog.CG)
dy.CG<-yres(orog.CG)
xmn.CG<-xmin(orog.CG)
xmx.CG<-xmax(orog.CG)
ymn.CG<-ymin(orog.CG)
ymx.CG<-ymax(orog.CG)
# extract all the cell values: zvalues[1] contains the orog[1,1] value
# Raster: cell numbers start at 1 in the upper left corner,
# and increase from left to right, and then from top to bottom
zvalues.FG<-extract(orog.FG,1:ncell(orog.FG))
xy<-xyFromCell(orog.FG,1:ncell(orog.FG))
x.FG<-sort(unique(xy[,1]))
y.FG<-sort(unique(xy[,2]),decreasing=T)
rc<-rowColFromCell(orog.FG,1:ncell(orog.FG))
aux<-as.vector(zvalues.FG)
gridmask.notNAs<-which(!is.na(aux))
zgrid<-aux[gridmask.notNAs]
xgrid<-xy[gridmask.notNAs,1]
ygrid<-xy[gridmask.notNAs,2]
rowgrid<-rc[gridmask.notNAs,1]
colgrid<-rc[gridmask.notNAs,2]
Lgrid.FG<-length(xgrid)
#
zvalues.CG<-extract(orog.CG,1:ncell(orog.CG))
xy.CG<-xyFromCell(orog.CG,1:ncell(orog.CG))
rc.CG<-rowColFromCell(orog.CG,1:ncell(orog.CG))
aux.CG<-as.vector(zvalues.CG)
gridmask.notNAs.CG<-which(!is.na(aux.CG))
zgrid.CG<-aux.CG[gridmask.notNAs.CG]
xgrid.CG<-xy.CG[gridmask.notNAs.CG,1]
ygrid.CG<-xy.CG[gridmask.notNAs.CG,2]
rowgrid.CG<-rc.CG[gridmask.notNAs.CG,1]
colgrid.CG<-rc.CG[gridmask.notNAs.CG,2]
Lgrid.CG<-length(xgrid.CG)
# variables facilitating the matrix to/from vector conversion
mat.CG.eve<-matrix(ncol=nx.CG,nrow=ny.CG,data=NA)
mat2vec.CG<-matrix(ncol=nx.CG,nrow=ny.CG,data=NA)
CG2FG<-matrix(ncol=(fact**2),nrow=Lgrid.CG,data=NA)
#to.write <- file("CG2FG.bin", "wb")
for (i in 1:Lgrid.CG) {
  mat2vec.CG[rowgrid.CG[i],colgrid.CG[i]]<-i
# following lines are used to create matrix CG2FG
#  aux<-which(rowgrid>=(fact*rowgrid.CG[i]-(fact-1)) & rowgrid<=(fact*rowgrid.CG[i]) &
#             colgrid>=(fact*colgrid.CG[i]-(fact-1)) & colgrid<=(fact*colgrid.CG[i]))
#  CG2FG[i,1:length(aux)]<-aux
#  writeBin(as.numeric(CG2FG[i,]),to.write,size=4)
}
#close(to.write)
#q()
if (!file.exists(file_CG2FG)) {
  print("FATAL ERROR: file not found!")
  print(file_CG2FG)
  quit(status=1)  
}
to.read<-file(file_CG2FG, "rb")
CG2FG<-t(matrix(readBin(to.read,numeric(),size=4,Lgrid.CG*fact**2),fact**2,Lgrid.CG))
close(to.read)
# clean memory
rm(xy,rc,aux,stackGeoGrid,mask.FG,rowgrid,colgrid)
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
myurl <- paste("http://klapp.oslo.dnmi.no/metnopub/production/metno?",
               "re=16&nod=NA&ct=text/plain&ddel=dot&del=semicolon",
  "&fy=",timeseq$year[1]+1900,formatC(timeseq$mon[1]+1,width=2,flag="0"),
                              formatC(timeseq$mday[1],width=2,flag="0"),
  "&ty=",timeseq$year[nhour]+1900,formatC(timeseq$mon[nhour]+1,width=2,flag="0"),
                          formatC(timeseq$mday[nhour],width=2,flag="0"),sep="")
#print(myurl)
o.cont<-1
while (o.cont<=10) {
  stataux<-NULL
  try(stataux <-read.table(myurl, header = TRUE,  sep = ";",
                         stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
                         encoding = "UTF-8", quote = "",na.string=-999))
  if (length(stataux)<10) {
    print("exit with error in command:")
        print(myurl)
        o.cont<-o.cont+1
        Sys.sleep(5)
#      q(status=1)
  } else {
    break
  }
}
if (o.cont>10) {
  print("Fatal Error in command:")
  print(myurl)
  print(stataux)
  q(status=1)
}
#if (is.null(stations)) {
#}
# stataux column names
# DEPARTMENT;DEPT_NO;MUNICIPALITY;MUNI_NO;ST_NAME;STNR;UTM_E;UTM_N;AMSL;LAT_DEC;LON_DEC;WMO_NO
# Select stations having the geographical information needed
# Find stations in Norway (Mainland).
#indx<-which( (as.numeric(stataux$STNR)<99710)  & 
#             (as.numeric(stataux$STNR)!=76900) &
#             (as.numeric(stataux$STNR)!=76933) &
#             (as.numeric(stataux$STNR)!=76930) &
#             (as.numeric(stataux$STNR)!=76928) &
#             (as.numeric(stataux$UTM_E)>=-65000) &
#             (as.numeric(stataux$UTM_N)<=8000000) &
#             !is.na(as.numeric(stataux$UTM_E)) & 
#             !is.na(as.numeric(stataux$UTM_N)) & 
#             !is.na(as.numeric(stataux$AMSL)) )
indx<-which( !is.na(as.numeric(stataux$UTM_E)) & 
             !is.na(as.numeric(stataux$UTM_N)) & 
             !is.na(as.numeric(stataux$AMSL)) )
stations<-data.frame(matrix(nrow=length(indx),ncol=4))
stations<-data.frame(matrix(nrow=length(indx),ncol=4))
names(stations)<-c("Stnr","z","x","y")
stations$Stnr<-as.numeric(stataux$STNR[indx])
stations$z<-stataux$AMSL[indx]
stations$x<-stataux$UTM_E[indx]
stations$y<-stataux$UTM_N[indx]
L.y.tot<-length(stations$Stnr)
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
#yo.superobs<-vector(mode="numeric",length=L.y.tot)
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
VecS<-as.numeric(as.vector(stations$Stnr))
# [] Read blacklists/errobs
if (file.exists(file_errobs)) {
  errobs<-read.table(file=file_errobs,header=T,sep=";")
  names(errobs)<-c("Stnr","Year","Month","Day","Hour","RR_1")
  err.timeseq<-as.POSIXlt(strptime(paste(errobs$Year,
                          formatC(errobs$Month,width=2,flag="0"),
                          formatC(errobs$Day,width=2,flag="0"),
                          formatC(errobs$Hour,width=2,flag="0"),sep=""),"%Y%m%d%H"),"UTC")
  err.Stnr<-unique(errobs$Stnr)
} else {
  print("Warning: file not found")
  print(file_errobs)
}
if (file.exists(file_blacklist_current)) {
  blacklist_current<-read.table(file=file_blacklist_current,header=T,sep=";")
  names(blacklist_current)<-c("Stnr")
} else {
  print("Warning: file not found")
  print(file_blacklist_current)
}
if (file.exists(file_blacklist_never)) {
  blacklist_never<-read.table(file=file_blacklist_never,header=T,sep=";")
  names(blacklist_never)<-c("Stnr")
} else {
  print("Warning: file not found")
  print(file_blacklist_never)
}
# [] compute Disth and Distz (symmetric) matrices: 
#  Disth(i,j)=horizontal distance between i-th station and j-th station [Km]
#  Distz(i,j)=elevation difference between i-th station and j-th station [m]
Disth<-(outer(VecY,VecY,FUN="-")**2.+outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
Distz<-abs(outer(VecZ,VecZ,FUN="-"))
print("list of station ids")
print(cbind(1:L.y.tot,VecS))
#------------------------------------------------------------------------------
# Elaborations
# xx/xx.CG raster structures usefull for map production
xx <-raster(ncol=nx.FG, nrow=ny.FG, xmn=xmn.FG, xmx=xmx.FG, ymn=ymn.FG, ymx=ymx.FG,
            crs=projstr)
xx[]<-NA
xx.CG <-raster(ncol=nx.CG, nrow=ny.CG, xmn=xmn.CG, xmx=xmx.CG,
               ymn=ymn.CG, ymx=ymx.CG, crs=projstr)
xx.CG[]<-NA
# Get and sum daily observations from KDVH
n.yo[]<-NA
yo[]<-NA
# set DQC flag to 0 for every station/observation
ydqc.flag[]<-0
# cycle over days to get data from KDVH
for (d in 1:nday) {
  yyyy<-dayseq$year[d]+1900
  mm<-dayseq$mon[d]+1
  dd<-dayseq$mday[d]
  date<-paste(dd,".",mm,".",yyyy,sep="")
  ulric<-paste("http://klapp.oslo.dnmi.no/metnopub/production/",
               "metno?re=14&p=RR&fd=",date,"&td=",date,
               "&nob=0.0&ddel=dot&del=semicolon&nmt=0",
               "&ct=text/plain&split=1&nod=-999&qa=2",sep="")
  print(ulric)
  o.cont<-1
  while (o.cont<=10) {
    o<-NULL
    try(o <- read.table(ulric, header = TRUE,  sep = ";", #nrows = nrows,
            stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
 	            encoding = "UTF-8", quote = "",na.strings="-999"))
    indx<-which(!is.na(o$RR))
    L.y.fromKDVH<-length(indx)
    if (L.y.fromKDVH<10) {
      print("exit with error in command:")
      print(ulric)
      o.cont<-o.cont+1
      Sys.sleep(5)
#      q(status=1)
    } else {
      break
    }
  }
  if (o.cont>10) {
    print("Fatal Error in command:")
    print(ulric)
    q(status=1)
  }
  if (exists("OBS")) rm(OBS)
  OBS<-data.frame(matrix(nrow=L.y.fromKDVH,ncol=7))
  names(OBS)<-c("Stnr","Year","Month","Day","RR","datetime")
  OBS$Stnr<-as.numeric(o$Stnr[indx])
  OBS$Year<-o$Year[indx]
  OBS$Month<-o$Month[indx]
  OBS$Day<-o$Day[indx]
  OBS$RR<-as.numeric(o$RR[indx])
  rm(indx)
  OBS$datetime<-as.POSIXlt(strptime(paste(OBS$Year,
                                    formatC(OBS$Month,width=2,flag="0"),
                                    formatC(OBS$Day,width=2,flag="0"),
                                    sep=""),"%Y%m%d"),"UTC")
  # obs having complete geo info 
  indx<-which(OBS$Stnr %in% stations$Stnr)
  L.y.fromKDVH.ofInterest<-length(indx)
  if (L.y.fromKDVH.ofInterest<=0) {
    print("Error: NO valid observations found")
    next
  }
  #
  obs.match.stn<-match(OBS$Stnr[indx],stations$Stnr)
  RR<-as.numeric(as.vector(OBS$RR[indx]))
  obs.match.stn.unique<-unique(obs.match.stn)
  n.obs.match.stn.unique<-length(obs.match.stn.unique)
  for (i in 1:n.obs.match.stn.unique) {
    indx<-which(obs.match.stn==obs.match.stn.unique[i])
    n.indx<-length(indx)
    if (n.indx==0) next 
    if (!any(!is.na(RR[indx]))) next
    # DQC check: plausibility check
    if (any(RR[indx]<yo.dqc.plausible.min) | any(RR[indx]>yo.dqc.plausible.max)) {
      flag.aux<-sapply(ydqc.flag[obs.match.stn.unique[i]],function(x){ as.integer(intToBits(x))})
      if (flag.aux[3,1]==0) ydqc.flag[obs.match.stn.unique[i]]<-ydqc.flag[obs.match.stn.unique[i]]+2**2
    }
    # DQC External
    Stnr.cur<-stations$Stnr[obs.match.stn.unique[i]]
    if ( Stnr.cur %in% err.Stnr ) {
      aux<-which(errobs$Stnr==Stnr.cur)
      err.timeseq<-as.POSIXlt(strptime(paste(errobs$Year[aux],
                              formatC(errobs$Month[aux],width=2,flag="0"),
                              formatC(errobs$Day[aux],width=2,flag="0"),
                              formatC(errobs$Hour[aux],width=2,flag="0"),sep=""),"%Y%m%d%H"),"UTC")
      indx.1<-which(OBS$Stnr==Stnr.cur)
      if (any(OBS$datetime[indx.1] %in% err.timeseq)) {
        flag.aux<-sapply(ydqc.flag[obs.match.stn.unique[i]],function(x){ as.integer(intToBits(x))})
        if (flag.aux[5,1]==0) ydqc.flag[obs.match.stn.unique[i]]<-ydqc.flag[obs.match.stn.unique[i]]+2**4
      }
    }
    #
    if (is.na(n.yo[obs.match.stn.unique[i]])) {
      n.yo[obs.match.stn.unique[i]]<-length(indx)
      yo[obs.match.stn.unique[i]]<-sum(RR[indx])
    } else {
      n.yo[obs.match.stn.unique[i]]<-n.yo[obs.match.stn.unique[i]]+length(indx)
      yo[obs.match.stn.unique[i]]<-yo[obs.match.stn.unique[i]]+sum(RR[indx])
    }
  }
} # cycle over days to get data from KDVH 
#
y.notNA<-which(!is.na(yo))
L.y.notNA<-length(y.notNA)
# DQC
aux<-which(n.yo[y.notNA]!=nday)
if (length(aux>0)) ydqc.flag[y.notNA][aux]<-ydqc.flag[y.notNA][aux]+2**3
# DQC blacklists
if (exists("blacklist_current")) {
  aux<-which(VecS %in% blacklist_current$Stnr)
  if (length(aux)>0) {
#    print("blacklist (current):")
#    print(cbind(VecS[aux],yo[aux]))
#    yo[aux]<-NA
    ydqc.flag[aux]<-ydqc.flag[aux]+2**6
  }
}
if (exists("blacklist_never")) {
  aux<-which(VecS %in% blacklist_never$Stnr)
  if (length(aux)>0) {
#    print("blacklist (never):")
#    print(cbind(VecS[aux],yo[aux]))
#    yo[aux]<-NA
    ydqc.flag[aux]<-ydqc.flag[aux]+2**7
  }
}
# mask observations that are outside the grid
# ydqc.flag=100 masked
for (b in 1:L.y.notNA) {
  r.b<-1+as.integer(abs(VecY[y.notNA[b]]-ymx.CG)/dy.CG)
  c.b<-1+as.integer((VecX[y.notNA[b]]-xmn.CG)/dx.CG)
  flag<-F
  for (ii in -1:1) {
    for (jj in -1:1) {
      r.ii<-r.b+ii
      c.jj<-c.b+jj
      if (r.ii<1 | r.ii>ny.CG | c.jj<1 | c.jj>nx.CG) {
        next
      }
      if (!is.na(mat2vec.CG[r.ii,c.jj])) flag<-T 
    }
  }
  if (!flag) ydqc.flag[y.notNA[b]]<-ydqc.flag[y.notNA[b]]+2**5
}
# DQC flag set to missing
ydqc.flag[y.notNA][ydqc.flag[y.notNA]==0]<-2**0
#
flag.aux<-sapply(ydqc.flag,function(x){ as.integer(intToBits(x))})
print(paste("number of station having at least one available observation =",L.y.notNA,"(tot=",L.y.tot,")"))
print(paste("  # station having at least one erroneous observation (plausibility check) =",length(which(flag.aux[3,]==1))))
print(paste("  # station having at least one erroneous observation       (external DQC) =",length(which(flag.aux[5,]==1))))
print(paste("  # station blacklisted for the time period considered =",length(which(flag.aux[7,]==1))))
print(paste("  # station blacklisted (permanently) =",length(which(flag.aux[8,]==1))))
print(paste("  # station in masked areas =",length(which(flag.aux[6,]==1))))
# define Vectors and Matrices
xidi.CG.wet<-vector(mode="numeric",length=Lgrid.CG)
xidi.CG.dry<-vector(mode="numeric",length=Lgrid.CG)
# loop for DQC
while (TRUE) {
  flag.aux<-sapply(ydqc.flag,function(x){ as.integer(intToBits(x))})
# vector with the positions (pointers to VecS) of: not NAs, not-in-masked-areas, valid observations 
# daily observations [L.yo.ok]
  yo.ok.pos<-which(flag.aux[1,]==1 | flag.aux[2,]==1)
  L.yo.ok<-length(yo.ok.pos)
# vectors with the positions (pointers to VecS) of valid (i.e. not NAs & DQC-ok)
# daily observations and "wet" or "dry" [L.yo.ok.wet,L.yo.ok.dry]
  yo.ok.pos.wet<-which((flag.aux[1,]==1 | flag.aux[2,]==1) & yo>=rr.inf)
  yo.ok.pos.dry<-which((flag.aux[1,]==1 | flag.aux[2,]==1) & yo<rr.inf)
  L.yo.ok.wet<-length(yo.ok.pos.wet)
  L.yo.ok.dry<-length(yo.ok.pos.dry)
  print(paste("observation not NAs, not-in-masked-areas, presumably good =",
              L.yo.ok,"(wet=",L.yo.ok.wet,"dry=",L.yo.ok.dry,")"))
#  print(cbind(VecS[yo.ok.pos.wet],yo[yo.ok.pos.wet]))
# case of no-rain over the whole domain
  if (L.yo.ok.wet==0) {
    print("no rain over the whole domain")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    cat(paste(yyyy.b,mm.b,dd.b,nday,
              round(VecS[yo.ok.pos],0),
              round(VecX[yo.ok.pos],0),
              round(VecY[yo.ok.pos],0),
              round(VecZ[yo.ok.pos],0),
              rep(NA,L.yo.ok), #eve.lab
              round(yo[yo.ok.pos],1),
              rep(NA,L.yo.ok), #yb
              rep(0,L.yo.ok), #ya
              rep(NA,L.yo.ok), #yav
              rep(NA,L.yo.ok), #yidi.eve
              rep(NA,L.yo.ok), #yidiv.eve
              round(ydqc.flag[yo.ok.pos],2),
              "\n",sep=";"),file=out.file.stn,append=T)
    # Figures
    xa.FG<-vector(mode="numeric",length=Lgrid.FG)
    xa.FG[]<-0
    xx[gridmask.notNAs]<-round(xa.FG,1)
    nogrid.ncout(grid=t(as.matrix(xx)),
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
    quit()
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  } # end case of no-rain over the whole domain
  # [] Contiguous Rain Area: First Guess
  Lsubsample.vec[]<-NA
  # vectors used in the identification of eve - first guess
  # ...aux-> temporary vectors
  # ...vec-> final vectors
  # neve.firstguess.vec-> number of eve.firstguess
  # eve.firstguess.vec-> list of stations (i.e. pointers to VecS) inside eve.firstguess
  #           matrix(i,j) i=1,neve.firstguess.vec j=1,Leve.firstguess.vec[i]
  #   note: an eve.firstguess could contain both wet and dry stations
  # Leve.firstguess.vec-> vector(i) (i=1,..,neve.firstguess.vec) number of stations inside i-th eve.firstguess 
  eve.firstguess.vec<-matrix(data=NA,ncol=L.yo.ok,nrow=L.yo.ok)
  Leve.firstguess.vec<-vector(mode="numeric",length=L.yo.ok)
  Leve.firstguess.vec[]<-NA
  neve.firstguess.vec<-0
  eve.firstguess.aux<-matrix(data=NA,ncol=L.yo.ok,nrow=L.yo.ok)
  Leve.firstguess.aux<-vector(mode="numeric",length=L.yo.ok)
  Leve.firstguess.aux[]<-NA
  neve.firstguess.aux<-0
  # [] Contiguous Rain Area: First Guess. Step I
  b.inc<-0
  for (b in yo.ok.pos.wet) {  # START: Cycle over wet stations 
    b.inc<-b.inc+1
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
    if (neve.firstguess.aux>0) {
      for (eve.firstguess in 1:neve.firstguess.aux) {
        if (Leve.firstguess.aux[eve.firstguess]==0) next
        if ( any(nodes %in% eve.firstguess.aux[eve.firstguess,1:Leve.firstguess.aux[eve.firstguess]]) ) {
          aux<-c(aux[which(!(aux %in% eve.firstguess.aux[eve.firstguess,1:Leve.firstguess.aux[eve.firstguess]]))],eve.firstguess.aux[eve.firstguess,1:Leve.firstguess.aux[eve.firstguess]])
          eve.firstguess.aux[eve.firstguess,]<-NA
          Leve.firstguess.aux[eve.firstguess]<-0
        }
      }
    }
    neve.firstguess.aux<-neve.firstguess.aux+1
    Leve.firstguess.aux[neve.firstguess.aux]<-length(c(aux,nodes[which(!(nodes%in%aux))]))
    eve.firstguess.aux[neve.firstguess.aux,1:Leve.firstguess.aux[neve.firstguess.aux]]<-c(aux,nodes[which(!(nodes%in%aux))])
    #DEBUG START
#    print("new wet-station group identified")
#    print(eve.firstguess.aux[neve.firstguess.aux,1:Leve.firstguess.aux[neve.firstguess.aux]])
    #DEBUG END      
  } # END: Cycle over the wet stations
# Contiguous Rain Area: First Guess. Step II
# reorganise eve labels
  neve.firstguess.vec<-0
  for (eve.firstguess in 1:neve.firstguess.aux) {
    if (Leve.firstguess.aux[eve.firstguess]==0) next
    neve.firstguess.vec<-neve.firstguess.vec+1
    Leve.firstguess.vec[neve.firstguess.vec]<-Leve.firstguess.aux[eve.firstguess]
    eve.firstguess.vec[neve.firstguess.vec,1:Leve.firstguess.vec[neve.firstguess.vec]]<-eve.firstguess.aux[eve.firstguess,1:Leve.firstguess.aux[eve.firstguess]]
  }
# [] Contiguous Rain Area: First Guess on the grid.CG
  print("++ define eve.firstguess on the grid.CG")
  lab.eve.CG.firstguess<-vector(mode="numeric",length=Lgrid.CG)
  close2i<-vector(mode="numeric",length=Lgrid.CG)
  close2i.wet<-vector(mode="numeric",length=Lgrid.CG)
  Disth2i.min.wet<-vector(mode="numeric",length=Lgrid.CG)
  lab.eve.CG.firstguess[]<-NA
#  tri.rr<-tri.mesh(VecX[yo.ok.pos],VecY[yo.ok.pos])
  for (i in 1:Lgrid.CG) {
    Disth.i<-( (ygrid.CG[i]-VecY[yo.ok.pos])**2+(xgrid.CG[i]-VecX[yo.ok.pos])**2 )**0.5 / 1000.
    close2i[i]<-yo.ok.pos[which.min(Disth.i)]
    Disth.i.wet<-( (ygrid.CG[i]-VecY[yo.ok.pos.wet])**2+(xgrid.CG[i]-VecX[yo.ok.pos.wet])**2 )**0.5 / 1000.
    close2i.wet[i]<-yo.ok.pos.wet[which.min(Disth.i.wet)]
    Disth2i.min.wet[i]<-min(Disth.i.wet)
    for (eve.firstguess in 1:neve.firstguess.vec) {
      if (close2i[i] %in% eve.firstguess.vec[eve.firstguess,1:Leve.firstguess.vec[eve.firstguess]]) {
        lab.eve.CG.firstguess[i]<-eve.firstguess
        break
      }
    }
    if (is.na(lab.eve.CG.firstguess[i])) lab.eve.CG.firstguess[i]<-0
  }
  rm(Disth.i,Disth.i.wet)
# [] Refine eve identification (Identify Wet and Dry Areas within the eve first guess)
  print("++ Refine eve identification (Identify Wet and Dry Areas within the eve first guess)")
  xa.Dh.eve.firstguess<-vector(mode="numeric",length=Lgrid.CG)
  xa.Dz.eve.firstguess<-vector(mode="numeric",length=Lgrid.CG)
  x.CG.wet<-vector(mode="numeric",length=Lgrid.CG)
  x.CG.dry<-vector(mode="numeric",length=Lgrid.CG)
  lab.eve.CG<-vector(mode="numeric",length=Lgrid.CG)
  n.eve<-0
# x.CG.@@@ -> =0 gridpoint not classified as @@@; =1 grid-point classified as @@@
# lab.eve.CG -> =0 gridpoint "dry"; =-1 gridpoint "wet", waiting to labelled;
#             n>0 gridpoint "wet" belonging to eve n (only for isolated wet station)
  x.CG.wet[]<-0
  x.CG.dry[]<-0
  lab.eve.CG[]<-0
  for (eve.firstguess in 1:neve.firstguess.vec) { # START: cycle over eve.firstguess
#    print(paste(eve.firstguess,"eve.firstguess +++++++++++++++++++++++++++++++++"))
    # define vectors / variables / indexes 
    Leve.firstguess<-Leve.firstguess.vec[eve.firstguess]
    VecX.eve.firstguess<-vector(mode="numeric",length=Leve.firstguess)
    VecY.eve.firstguess<-vector(mode="numeric",length=Leve.firstguess)
    VecZ.eve.firstguess<-vector(mode="numeric",length=Leve.firstguess)
    VecS.eve.firstguess<-vector(mode="numeric",length=Leve.firstguess)
    yo.eve.firstguess<-vector(mode="numeric",length=Leve.firstguess)
    VecX.eve.firstguess<-VecX[eve.firstguess.vec[eve.firstguess,1:Leve.firstguess]]
    VecY.eve.firstguess<-VecY[eve.firstguess.vec[eve.firstguess,1:Leve.firstguess]]
    VecZ.eve.firstguess<-VecZ[eve.firstguess.vec[eve.firstguess,1:Leve.firstguess]]
    VecS.eve.firstguess<-VecS[eve.firstguess.vec[eve.firstguess,1:Leve.firstguess]]
    yo.eve.firstguess<-yo[eve.firstguess.vec[eve.firstguess,1:Leve.firstguess]]
    yoindx.eve.firstguess.wet<-which(yo.eve.firstguess>=rr.inf)
    yoindx.eve.firstguess.dry<-which(yo.eve.firstguess<rr.inf)
    yoindx.wet<-eve.firstguess.vec[eve.firstguess,yoindx.eve.firstguess.wet]
    yoindx.dry<-eve.firstguess.vec[eve.firstguess,yoindx.eve.firstguess.dry]
    n.eve.firstguess.wet<-length(yoindx.eve.firstguess.wet)
    n.eve.firstguess.dry<-length(yoindx.eve.firstguess.dry)
    xindx.eve.CG<-which(lab.eve.CG.firstguess==eve.firstguess)
    Lgrid.eve.firstguess<-length(xindx.eve.CG)
    # isolated wet-station
    if (n.eve.firstguess.wet==1) n.eve<-n.eve+1
    # more than one single wet-station
    # define Dz for the current eve (min allowed value is 500 m)
    xa.Dz.eve.firstguess[xindx.eve.CG]<-max(500,max(zgrid.CG[xindx.eve.CG]))
    for (i in xindx.eve.CG) {    #START: cycle over xindx.eve.CG
      Disth.eve.firstguess<-( (ygrid.CG[i]-VecY.eve.firstguess)**2+(xgrid.CG[i]-VecX.eve.firstguess)**2 )**0.5 / 1000.
      Distz.eve.firstguess<-abs(zgrid.CG[i]-VecZ.eve.firstguess)
      close2i.ord<-order(Disth.eve.firstguess,decreasing=F)
#      xa.Dh.eve.firstguess[i] <- max(20,Disth.eve.firstguess[close2i.ord[1]])
      xa.Dh.eve.firstguess[i] <- Disth.eve.firstguess[close2i.ord[1]]
      notsofar.wet.i<-which(Disth.eve.firstguess[yoindx.eve.firstguess.wet]<(5*xa.Dh.eve.firstguess[i]))
      notsofar.dry.i<-which(Disth.eve.firstguess[yoindx.eve.firstguess.dry]<(5*xa.Dh.eve.firstguess[i]))
      if (length(notsofar.wet.i)==0) {
        x.CG.dry[i]<-1
      } else if (length(notsofar.dry.i)==0) {
        x.CG.wet[i]<-1
        if (n.eve.firstguess.wet==1) {
          lab.eve.CG[i]<-n.eve
        } else {
          lab.eve.CG[i]<--1
        }
      } else {
        G.wet<-exp(-0.5*(Disth.eve.firstguess[yoindx.eve.firstguess.wet[notsofar.wet.i]]/xa.Dh.eve.firstguess[i])**2.
                   -0.5*(Distz.eve.firstguess[yoindx.eve.firstguess.wet[notsofar.wet.i]]/xa.Dz.eve.firstguess[i])**2.)
        D.wet<-as.matrix(exp(-0.5*(Disth[yoindx.wet[notsofar.wet.i],
                                         yoindx.wet[notsofar.wet.i]]/xa.Dh.eve.firstguess[i])**2.
                             -0.5*(Distz[yoindx.wet[notsofar.wet.i],
                                         yoindx.wet[notsofar.wet.i]]/xa.Dz.eve.firstguess[i])**2. ))
        D.wet[row(D.wet)==col(D.wet)]<-D.wet[row(D.wet)==col(D.wet)]+eps2
        if (n.eve.firstguess.wet>1) {
          InvD.wet<-solve(D.wet)
          K.wet<-G.wet%*%InvD.wet
          rm(InvD.wet)
        } else {
          K.wet<-G.wet/(1+eps2)
        }
        xidi.CG.wet[i]<-sum(K.wet)
        # dry
        G.dry<-exp(-0.5*(Disth.eve.firstguess[yoindx.eve.firstguess.dry[notsofar.dry.i]]/xa.Dh.eve.firstguess[i])**2.-
                    0.5*(Distz.eve.firstguess[yoindx.eve.firstguess.dry[notsofar.dry.i]]/xa.Dz.eve.firstguess[i])**2.)
        D.dry<-as.matrix(exp(-0.5*(Disth[yoindx.dry[notsofar.dry.i],
                                         yoindx.dry[notsofar.dry.i]]/xa.Dh.eve.firstguess[i])**2.+
                             -0.5*(Distz[yoindx.dry[notsofar.dry.i],
                                         yoindx.dry[notsofar.dry.i]]/xa.Dz.eve.firstguess[i])**2. ))
        D.dry[row(D.dry)==col(D.dry)]<-D.dry[row(D.dry)==col(D.dry)]+eps2
        InvD.dry<-solve(D.dry)
        K.dry<-G.dry%*%InvD.dry
        rm(InvD.dry)
        xidi.CG.dry[i]<-sum(K.dry)
        if (xidi.CG.wet[i]>xidi.CG.dry[i]) {
          x.CG.wet[i]<-1
          if (n.eve.firstguess.wet==1) {
            lab.eve.CG[i]<-n.eve
          } else {
            lab.eve.CG[i]<--1
          }
        } else {
          x.CG.dry[i]<-1
        }
      }
    } # END: cycle over xindx.eve.CG
    if (n.eve.firstguess.wet==1) {
      aux<-which(lab.eve.CG==n.eve)
      if (length(aux)==0) {
        r.b<-1+as.integer(abs(VecY[yoindx.wet]-ymx.CG)/dy.CG)
        c.b<-1+as.integer((VecX[yoindx.wet]-xmn.CG)/dx.CG)
        i<-mat2vec.CG[r.b,c.b]
        for (ii in -1:1) {
          for (jj in -1:1) {
            if ((r.b+ii)>ny.CG | (r.b+ii)<1) next
            if ((c.b+jj)>nx.CG | (c.b+jj)<1) next
            i<-mat2vec.CG[(r.b+ii),(c.b+jj)]
            x.CG.wet[i]<-1 
            x.CG.dry[i]<-0
            lab.eve.CG[i]<-n.eve
          }
        }
      }
    }
  } # END: cycle over eve.firstguess
# assure that a wet observation fell in a eve
  for (b in yo.ok.pos.wet) {
    r.b<-1+as.integer(abs(VecY[b]-ymx.CG)/dy.CG)
    c.b<-1+as.integer((VecX[b]-xmn.CG)/dx.CG)
    i<-mat2vec.CG[r.b,c.b]
    flag<-F
    if (is.na(i)) flag<-T
    if (!flag) {
      if (x.CG.wet[i]!=1) flag<-T
    }
    if (flag) {
      for (ii in -1:1) {
        for (jj in -1:1) {
          if ((r.b+ii)>ny.CG | (r.b+ii)<1) next
          if ((c.b+jj)>nx.CG | (c.b+jj)<1) next
          i<-mat2vec.CG[(r.b+ii),(c.b+jj)]
          if (is.na(i)) next
          x.CG.wet[i]<-1 
          x.CG.dry[i]<-0
          lab.eve.CG[i]<--1
        }
      }
    }
  }
# [] eve Labelling (Identify final eve on the grid)
  print("++ eve Labelling (Identify final eve on the grid)")
#   assign neighbouring wet grid-points to the same eve
  # ...pre = preliminary
  n.eve.pre<-n.eve
  for (i in which(x.CG.wet==1)) {
    # cycle over un-labelled wet gridpoints
    if (lab.eve.CG[i]!=-1) next
    aux<-vector()
    n.aux<-0
    for (j in which(x.CG.wet==1)) {
      if (abs(rowgrid.CG[i]-rowgrid.CG[j])<=1 & abs(colgrid.CG[i]-colgrid.CG[j])<=1 ) {
        n.aux<-n.aux+1
        aux[n.aux]<-j
      }
    }
    n.eve.pre<-n.eve.pre+1
    for (j in 1:n.aux) {
      if (lab.eve.CG[aux[j]]>0) lab.eve.CG[which(lab.eve.CG==lab.eve.CG[aux[j]])]<-n.eve.pre
      if (lab.eve.CG[aux[j]]==-1) lab.eve.CG[aux[j]]<-n.eve.pre
    }
  }
# reorganise eve labels
  if (n.eve.pre>n.eve) {
    for (n in (n.eve+1):n.eve.pre) {
      aux<-which(lab.eve.CG==n)
      if (length(aux)>0) {
        n.eve<-n.eve+1
        lab.eve.CG[aux]<-n.eve
      }
    }
  }
# [] assign stations at eve
#  remarks: 1. each eve must include at least 1 observation; 
#   2. it is possible to have wet-observations not included in any of the eve 
#      (for example, in areas where there are a lot of observations and most of them are "dry")
# eve on the grid: from vector to matrix
  for (i in which(x.CG.wet==1)) {
    mat.CG.eve[rowgrid.CG[i],colgrid.CG[i]]<-lab.eve.CG[i]
#      print(paste(rowgrid.CG[i],colgrid.CG[i],mat.CG.eve[rowgrid.CG[i],colgrid.CG[i]]))
  }
# eve at wet-station locations
  y.eve[]<-NA
  for (b in yo.ok.pos.wet) {
    r.b<-1+as.integer(abs(VecY[b]-ymx.CG)/dy.CG)
    c.b<-1+as.integer((VecX[b]-xmn.CG)/dx.CG)
    y.eve[b]<-mat.CG.eve[r.b,c.b]
    # case: wet station on the edge (but outside) of eve
    if (is.na(y.eve[b])) {
      for (ii in -1:1) {
        for (jj in -1:1) {
          if ((r.b+ii)>ny.CG | (r.b+ii)<1) next
          if ((c.b+jj)>nx.CG | (c.b+jj)<1) next
          if (!is.na(mat.CG.eve[r.b+ii,c.b+jj])) y.eve[b]<-mat.CG.eve[r.b+ii,c.b+jj]
        }
      }
    }
#  print(paste(xmn.CG,ymx.CG,VecS[b],yo[b],VecX[b],VecY[b],r.b,c.b,y.eve[b]))
  }
# search for possible eves on the grid not including any station and assign
# them to the "nearest" eve which includes stations.
  aux.y.eve<-unique(y.eve[which(!is.na(y.eve))])
  aux.x.eve<-unique(lab.eve.CG[which(!is.na(lab.eve.CG) & lab.eve.CG!=0)])
#  print(aux.y.eve)
#  print(aux.x.eve)
#  print(!(aux.x.eve %in% aux.y.eve))
#  print(any(!(aux.x.eve %in% aux.y.eve)))
  if (any( !(aux.x.eve %in% aux.y.eve) )) {
    aux.indx<-which( !(aux.x.eve %in% aux.y.eve) )
#    print(aux.indx)
    for (a in aux.indx) {
      indx<-which(lab.eve.CG==aux.x.eve[a])
      new.eve<-y.eve[close2i.wet[indx[1]]]
      lab.eve.CG[which(lab.eve.CG==aux.x.eve[a])]<-new.eve
    }
  }
  eve.labels<-sort(unique(lab.eve.CG[which(!is.na(lab.eve.CG) & lab.eve.CG!=0)]))
  n.eve<-length(eve.labels)
#  print(paste("n.eve",n.eve))
#  print("eve.labels")
#  print(eve.labels)
# eve on the grid: from vector to matrix
  for (i in which(x.CG.wet==1)) {
    mat.CG.eve[rowgrid.CG[i],colgrid.CG[i]]<-lab.eve.CG[i]
  }
# [] Analysis
  print("++ eves Analysis")
# + data structures definition
# FG -> Fine Grid
  xa.FG<-vector(mode="numeric",length=Lgrid.FG)
  xb.FG<-vector(mode="numeric",length=Lgrid.FG)
  xidi.eve.FG<-vector(mode="numeric",length=Lgrid.FG)
  lab.eve.FG<-vector(mode="numeric",length=Lgrid.FG)
# CG -> Coarse Grid
  xa.CG<-vector(mode="numeric",length=Lgrid.CG)
  xb.CG<-vector(mode="numeric",length=Lgrid.CG)
  xidi.eve.CG<-vector(mode="numeric",length=Lgrid.CG)
# Y vectors (station locations)
  ya<-vector(mode="numeric",length=L.y.tot)
  yb<-vector(mode="numeric",length=L.y.tot)
  ybv.mat<-matrix(ncol=L.y.tot,nrow=L.y.tot,data=0.)
  yidi.eve<-vector(mode="numeric",length=L.y.tot)
  yidiv.eve<-vector(mode="numeric",length=L.y.tot)
  yo.superobs<-vector(mode="numeric",length=L.y.tot)
  yav<-vector(mode="numeric",length=L.y.tot)
  ya.tmp<-vector(mode="numeric",length=L.y.tot)
  yav.tmp<-vector(mode="numeric",length=L.y.tot)
# vectors used in the last iteration of the analysis procedure
  Dh.eve.smallest<-vector(mode="numeric",length=n.eve)
  Dz.eve.smallest<-vector(mode="numeric",length=n.eve)
#  eps2.eve.smallest<-vector(mode="numeric",length=n.eve)
# eves descriptive vectors
  n.y.eve<-vector(mode="numeric",length=n.eve)
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
  xidi.eve.FG[]<-0
  lab.eve.FG[]<-0
  for (i in which(x.CG.wet==1)) {
    aux<-CG2FG[i,which(!is.na(CG2FG[i,]))]
    lab.eve.FG[aux]<-lab.eve.CG[i]
  }
  xidi.eve.CG[]<-0
  yo.superobs[]<-NA
  yo.superobs[yo.ok.pos]<-0
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
  Dh.eve.smallest[]<-NA
  Dz.eve.smallest[]<-NA
  ell.locx.eve[]<-NA
  ell.locy.eve[]<-NA
  ell.smajor.eve[]<-NA
  ell.sminor.eve[]<-NA
  ell.smadir.eve[]<-NA
#  eps2.eve.smallest[]<-NA
# + Ellipsoid hulls
# ellipsoid hulls computation could fail (i.e. small eves) and return NAs. 
  for (n in 1:n.eve) {
    xindx.eve.CG<-which(lab.eve.CG==eve.labels[n])
    Lgrid.eve<-length(xindx.eve.CG)
    xy<-cbind(xgrid.CG[xindx.eve.CG],ygrid.CG[xindx.eve.CG])
    ell<-ellipsoidhull(xy)
#    to plot: lines(predict(ell), col="blue")
    ell.locx.eve[n]<-ell$loc[1]
    ell.locy.eve[n]<-ell$loc[2]
    if (is.na(ell$cov[1,1]) | is.na(ell$cov[2,2]) | is.na(ell$cov[1,2])) next
    if (ell$cov[1,1]==0 | ell$cov[2,2]==0 | ell$cov[1,2]==0) next
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
  rm(xindx.eve.CG,xy,ell,eigenval,eigenvec,e)
# + Analysis cycle over eves 
  print("++ Analysis cycle over eves")
  for (n in 1:n.eve) {
    xindx.eve.CG<-which(lab.eve.CG==eve.labels[n])
    xindx.eve.FG<-which(lab.eve.FG==eve.labels[n])
    Lgrid.eve.CG<-length(xindx.eve.CG)
    Lgrid.eve.FG<-length(xindx.eve.FG)
    yindx<-which(y.eve==eve.labels[n])
    n.y.eve[n]<-length(yindx)
    print(paste("+ eve lab >",eve.labels[n],
                "(",n,"/",n.eve,")",
                " #obs=",n.y.eve[n],sep=""))
#   + eve from isolated station
    if (n.y.eve[n]==1) {
      xa.FG[xindx.eve.FG]<-yo[yindx]
      xa.CG[xindx.eve.CG]<-yo[yindx]
      ya[yindx]<-yo[yindx]
      xb.FG[xindx.eve.FG]<-NA
      xidi.eve.FG[xindx.eve.FG]<--1
      xb.CG[xindx.eve.CG]<-NA
      xidi.eve.CG[xindx.eve.CG]<--1
      yb[yindx]<-NA
      yav[yindx]<-NA
      yidi.eve[yindx]<-NA
      yidiv.eve[yindx]<-NA
      next
    }
#   + eve including more than one station
#     + Larger-scale Background = mode(distribution of observations)
    histog<-hist(yo[yindx],breaks=seq(0,round((max(yo[yindx])+1),0)),plot=F)
    mode<-which.max(histog$counts)-0.5
    xb.CG[xindx.eve.CG]<-mode
    xb.FG[xindx.eve.FG]<-mode
    yb[yindx]<-mode
    ybv.mat[yindx,yindx]<-mode
#     + sequence of horizontal decorrelation length scales Dh.seq: largest to smallest
#   general case (a.2):
#    from ellipsoid hull major semiaxis (large-scale reference value)
#    to the minimum distance between two stations (small-scale reference value)
#   details: particular cases (i.e. semiaxis is NA, semiaxis is smaller then min distance);
#            two different steps are chosen in defining the sequence:
#            above a Dh refence value (100 Km) the step is 20 Km
#            below is 10 Km (reduction of computational time) 
    Disth.aux<-Disth[yindx,yindx]
    # minimum distance between two eve stations
    min.Disth<-min(Disth.aux[row(Disth.aux)!=col(Disth.aux)])
    min.Dh.seq<-round(max(min.Dh.seq.allowed,min.Disth),0)
    if (!is.na(ell.smajor.eve[n])) {                              # (a)
      max.Dh.seq<-(round(ell.smajor.eve[n]/10,0)+1)*10
      if (min.Dh.seq>max.Dh.seq) {                                # (a.0)
        max.Dh.seq<-min.Dh.seq
        min.Dh.seq<-(round(ell.smajor.eve[n]/10,0)+1)*10
      }
      max.ind<-max(which(Dh.seq.reference>max.Dh.seq))
      if (min.Dh.seq<=min.Dh.seq.allowed) {
        min.ind<-length(Dh.seq.reference)
      } else {
        min.ind<-min(which(Dh.seq.reference<min.Dh.seq))
      }
      if (max.ind==min.ind) max.ind<-min.ind+1
      Dh.seq<-Dh.seq.reference[max.ind:min.ind]
      n.Dh.seq<-length(Dh.seq)
    } else {                                                      # (b)
      Dh.seq<-seq(2*min.Dh.seq,min.Dh.seq,length=n.Dh.seq)
    }
    print("Dh.seq")
    print(Dh.seq)
#     + define matrices on Coarse Grid
    G.CG<-matrix(ncol=n.y.eve[n],nrow=Lgrid.eve.CG,data=0.)
    aux.CG<-(outer(ygrid.CG[xindx.eve.CG],VecY[yindx],FUN="-")**2. +
             outer(xgrid.CG[xindx.eve.CG],VecX[yindx],FUN="-")**2.)**0.5/1000.
    auxz.CG<-abs(outer(zgrid.CG[xindx.eve.CG],VecZ[yindx],FUN="-"))
    InvD<-matrix(data=0,ncol=n.y.eve[n],nrow=n.y.eve[n])
#     + cycle on horizontal decorrelation length scales
    flag.firsttime<-T
    for (Dh.test in Dh.seq) {
      print(Dh.test)
      D.test.part<-exp(-0.5*(Disth[yindx,yindx]/Dh.test)**2.)
      J<-1000000
      if (Dh.test>50) Dz.seq<-Dz.seq.gt50
      if (Dh.test<=50) Dz.seq<-Dz.seq.le50
#     + cycle on vertical decorrelation length scales
      for (Dz.test in Dz.seq) {
        # compute superobservations according to decorrelation lenght-scales
        yo.superobs.tmp<-vector(mode="numeric",length=L.y.tot)
        for (b in yindx) {
          aux.yo<-which(Disth[b,yindx]<=(Dh.test) & Distz[b,yindx]<=(Dz.test))
          yo.superobs.tmp[b]<-mean(yo[yindx[aux.yo]])
        }
        # analysis at eve station locations
        D.test<-D.test.part*exp(-0.5*(Distz[yindx,yindx]/Dz.test)**2.)
        S.test<-D.test
        D.test[row(D.test)==col(D.test)]<-D.test[row(D.test)==col(D.test)]+eps2
        InvD.test<-solve(D.test)
        W.test<-S.test%*%InvD.test
        ya.tmp[yindx]<-yb[yindx] + W.test %*% (yo.superobs.tmp[yindx]-yb[yindx])
        yav.tmp[yindx]<-yo.superobs.tmp[yindx] +
                           1./(1.-diag(W.test)) * (ya.tmp[yindx]-yo.superobs.tmp[yindx])
        ya.tmp[yindx][ya.tmp[yindx]<rr.inf]<-rr.inf
        yav.tmp[yindx][yav.tmp[yindx]<rr.inf]<-rr.inf
        # cost function
        J.tmp<-(1/n.y.eve[n]*sum((st.log(yav.tmp[yindx])-st.log(yo.superobs.tmp[yindx]))**2))**0.5
#        J.tmp<-(1/n.y.eve[n]*sum((yav.tmp[yindx]-yo.superobs.tmp[yindx])**2))**0.5
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
          yo.superobs[yindx]<-yo.superobs.tmp[yindx]
        }
      } # end cycle on vertical decorrelation lenght scales
#     + Dh>100Km then analysis is on CG; Dh<100Km analysis on FG
      if (Dh.test>Dh.seq.ref) { # analysis on CG
        G.CG<-exp(-0.5*(aux.CG/Dh.test)**2.-0.5*(auxz.CG/Dz.choice)**2.)
        K.CG<-G.CG%*%InvD
        xa.CG[xindx.eve.CG]<-xb.CG[xindx.eve.CG]+K.CG%*%(yo.superobs[yindx]-yb[yindx])
        xa.CG[xindx.eve.CG][xa.CG[xindx.eve.CG]<rr.inf]<-rr.inf
        xidi.eve.CG[xindx.eve.CG]<-xidi.eve.CG[xindx.eve.CG]+rowSums(K.CG)
        xb.CG[xindx.eve.CG]<-xa.CG[xindx.eve.CG]
      } else { # analysis on FG
        # get background from CG, if it is needed
        if (Dh.test!=Dh.seq[1] & flag.firsttime) {
          for (i in xindx.eve.CG) {
            aux.CG2FG<-CG2FG[i,which(!is.na(CG2FG[i,]))]
            xb.FG[aux.CG2FG]<-xb.CG[i]
            xidi.eve.FG[aux.CG2FG]<-xidi.eve.CG[i]
            rm(aux.CG2FG)
          }
        }
        flag.firsttime<-F
        # analysis on FG is done iteratively (limit memory use)
        i<-0
        while ((i*ndim.FG.iteration)<Lgrid.eve.FG) {
          start<-i*ndim.FG.iteration+1
          end<-(i+1)*ndim.FG.iteration
          if (end>Lgrid.eve.FG) {
            end<-Lgrid.eve.FG
          }
          ndimaux.FG<-end-start+1
          aux.FG<-matrix(ncol=n.y.eve[n],nrow=ndimaux.FG,data=0.)
          auxz.FG<-matrix(ncol=n.y.eve[n],nrow=ndimaux.FG,data=0.)
          aux.FG<-(outer(ygrid[xindx.eve.FG[start:end]],VecY[yindx],FUN="-")**2. +
                outer(xgrid[xindx.eve.FG[start:end]],VecX[yindx],FUN="-")**2.)**0.5/1000.
          auxz.FG<-abs(outer(zgrid[xindx.eve.FG[start:end]],VecZ[yindx],FUN="-"))
          G.FG<-matrix(ncol=n.y.eve[n],nrow=ndimaux.FG,data=0.)
          G.FG<-exp(-0.5*(aux.FG/Dh.test)**2.-0.5*(auxz.FG/Dz.choice)**2.)
          K.FG<-G.FG%*%InvD
          xa.FG[xindx.eve.FG[start:end]]<-xb.FG[xindx.eve.FG[start:end]]+K.FG%*%(yo.superobs[yindx]-yb[yindx])
          xa.FG[xindx.eve.FG[start:end]][xa.FG[xindx.eve.FG[start:end]]<rr.inf]<-rr.inf
          xidi.eve.FG[xindx.eve.FG[start:end]]<-xidi.eve.FG[xindx.eve.FG[start:end]]+rowSums(K.FG)
          rm(aux.FG,auxz.FG,G.FG,K.FG)
          i<-i+1
        }
        # current analysis is next iteration background
        xb.FG[xindx.eve.FG]<-xa.FG[xindx.eve.FG]
      } # end: analysis on FG
#     + update yidi, yidiv and ybv.mat
      yidi.eve[yindx]<-yidi.eve[yindx]+rowSums(W)
      yidiv.eve[yindx]<-yidiv.eve[yindx]+
                         rep(1,n.y.eve[n])+ 1./(1.-diag(W)) * (rowSums(W)-rep(1,n.y.eve[n]))
      i<-0
      for (b in yindx) {
        i<-i+1
#        InvD.1<-solve(D[-i,-i])
        InvD.1<-InvD[-i,-i]-1/InvD[i,i]*outer(InvD[-i,i],InvD[i,-i])
        W.1<-S[,-i]%*%InvD.1
        ya.tmp[yindx]<-ybv.mat[yindx,b] + W.1 %*% (yo.superobs[yindx[-i]]-ybv.mat[yindx[-i],b])
        ybv.mat[yindx,b]<-ya.tmp[yindx]
      }
      # current analysis is next iteration background
      yb[yindx]<-ya[yindx]
    } # end cycle on horizontal decorrelation lenght scales
    rm(D.test.part,D.test,S.test,W.test,InvD.test,InvD,D,S,W)
    if (exists("K.CG")) rm(K.CG,G.CG)
    rm(InvD.1,W.1)
    Dh.eve.smallest[n]<-max(Dh.seq[n.Dh.seq]/2,min.Disth)
    Dz.eve.smallest[n]<-Dz.choice
#    eps2.eve.smallest[n]<-0.1
  } # END CYCLE over eves
  break
} # end of DQC loop
# Analysis on high-res grid
print("++ Analysis on the high-resolution grid")
# Fine Grid Analysis: cycle over eves
for (n in 1:n.eve) {
  xindx.eve.FG<-which(lab.eve.FG==eve.labels[n])
  Lgrid.eve.FG<-length(xindx.eve.FG)
  area.eve[n]<-Lgrid.eve.FG*area.1cell.FG # Area Km**2
  yindx<-which(y.eve==eve.labels[n])
  if (n.y.eve[n]==1) {
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
  if (Dh.eve.smallest[n]>Dh.seq.ref) {
    xindx.eve.CG<-which(lab.eve.CG==eve.labels[n])
    for (i in xindx.eve.CG) {
      aux.CG2FG<-CG2FG[i,which(!is.na(CG2FG[i,]))]
      xb.FG[aux.CG2FG]<-xb.CG[i]
      xidi.eve.FG[aux.CG2FG]<-xidi.eve.CG[i]
      rm(aux.CG2FG)
    }
  }
  D<-exp(-0.5*(Disth[yindx,yindx]/Dh.eve.smallest[n])**2.
         -0.5*(Distz[yindx,yindx]/Dz.eve.smallest[n])**2.)
  S<-D
  D[row(D)==col(D)]<-D[row(D)==col(D)]+eps2
  InvD<-solve(D)
  W<-S%*%InvD
  ya[yindx]<-yb[yindx]+W%*%(yo[yindx]-yb[yindx])
  ya[yindx][ya[yindx]<rr.inf]<-rr.inf
  i<-0
  for (b in yindx) {
    i<-i+1
#    InvD.1<-solve(D[-i,-i])
    InvD.1<-InvD[-i,-i]-1/InvD[i,i]*outer(InvD[-i,i],InvD[i,-i])
    W.1<-S[,-i]%*%InvD.1
    ya.tmp[yindx]<-ybv.mat[yindx,b] + W.1 %*% (yo[yindx[-i]]-ybv.mat[yindx[-i],b])
    yav[yindx[i]]<-ya.tmp[yindx[i]]
    if (yav[yindx[i]]<rr.inf) yav[yindx[i]]<-rr.inf
  }
  yidi.eve[yindx]<-yidi.eve[yindx]+rowSums(W)
  yidiv.eve[yindx]<-yidiv.eve[yindx]+
                     rep(1,n.y.eve[n])+ 1./(1.-diag(W)) * (rowSums(W)-rep(1,n.y.eve[n]))
  rm(W,S,D)
  i<-0
  while ((i*ndim.FG.iteration)<Lgrid.eve.FG) {
    start<-i*ndim.FG.iteration+1
    end<-(i+1)*ndim.FG.iteration
    if (end>Lgrid.eve.FG) {
      end<-Lgrid.eve.FG
    }
    ndimaux<-end-start+1
    aux<-matrix(ncol=n.y.eve[n],nrow=ndimaux,data=0.)
    auxz<-matrix(ncol=n.y.eve[n],nrow=ndimaux,data=0.)
    aux<-(outer(ygrid[xindx.eve.FG[start:end]],VecY[yindx],FUN="-")**2. +
          outer(xgrid[xindx.eve.FG[start:end]],VecX[yindx],FUN="-")**2.)**0.5/1000.
    auxz<-abs(outer(zgrid[xindx.eve.FG[start:end]],VecZ[yindx],FUN="-"))
    G<-matrix(ncol=n.y.eve[n],nrow=ndimaux,data=0.)
    G<-exp(-0.5*(aux/Dh.eve.smallest[n])**2.-0.5*(auxz/Dz.eve.smallest[n])**2.)
    K<-G%*%InvD
    xa.FG[xindx.eve.FG[start:end]]<-xb.FG[xindx.eve.FG[start:end]]+K%*%(yo[yindx]-yb[yindx])
    xa.FG[xindx.eve.FG[start:end]][xa.FG[xindx.eve.FG[start:end]]<rr.inf]<-rr.inf
    xidi.eve.FG[xindx.eve.FG[start:end]]<-xidi.eve.FG[xindx.eve.FG[start:end]]+rowSums(K)
    rm(aux,auxz,G,K)
    i<-i+1
  }
  idi.norm.fac[n]<-max(xidi.eve.FG[xindx.eve.FG],yidi.eve[yindx],yidiv.eve[yindx])
  xidi.eve.FG[xindx.eve.FG]<-xidi.eve.FG[xindx.eve.FG]/idi.norm.fac[n]
  yidi.eve[yindx]<-yidi.eve[yindx]/idi.norm.fac[n]
  yidiv.eve[yindx]<-yidiv.eve[yindx]/idi.norm.fac[n]
  volume.eve[n]<-sum(xa.FG[xindx.eve.FG])
  meanidi.x.eve[n]<-mean(xidi.eve.FG[xindx.eve.FG])
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
} # Fine Grid Analysis: cycle over eves
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
print("++ Output")
# Station Points - Write output on file 
cat(paste(yyyy.b,mm.b,dd.b,nday,
          round(VecS[yo.ok.pos],0),
          round(VecX[yo.ok.pos],0),
          round(VecY[yo.ok.pos],0),
          round(VecZ[yo.ok.pos],0),
          round(y.eve[yo.ok.pos],0),
          round(yo[yo.ok.pos],1),
          round(yb[yo.ok.pos],2),
          round(ya[yo.ok.pos],2),
          round(yav[yo.ok.pos],2),
          round(100*yidi.eve[yo.ok.pos],2),
          round(100*yidiv.eve[yo.ok.pos],2),
          round(ydqc.flag[yo.ok.pos],2),
          "\n",sep=";"),file=out.file.stn,append=T)
cat(paste(yyyy.b,mm.b,dd.b,nday,
          eve.labels,
          n.y.eve,
          round(area.eve,0),
          round(volume.eve,2),
          round(100*meanidi.x.eve,2),
          round(100*meanidi.y.eve,2),
          round(100*meanidiv.y.eve,2),
          round(meanrain.eve,2),
          round(maxrain.x.eve,2),
          round(maxrain.yo.eve,2),
          round(maxrain.ya.eve,2),
          round(maxrain.yav.eve,2),
          round(ell.locx.eve,1),
          round(ell.locy.eve,1),
          round(ell.smajor.eve,1),
          round(ell.sminor.eve,1),
          round(ell.smadir.eve,1),
          round(cv.rel.eve.all,4),
          round(cv.bias.eve.all,4),
          round(cv.rmse.eve.all,4),
          round(cv.made.eve.all,4),
          round(cv.rel.eve.q50,4),
          round(cv.bias.eve.q50,4),
          round(cv.rmse.eve.q50,4),
          round(cv.made.eve.q50,4),
          round(100*meanidiv.y.eve.q50,2),
          round(n.q50,0),
          round(cv.rel.eve.q75,4),
          round(cv.bias.eve.q75,4),
          round(cv.rmse.eve.q75,4),
          round(cv.made.eve.q75,4),
          round(100*meanidiv.y.eve.q75,2),
          round(n.q75,0),
          round(idi.norm.fac,5),
          "\n",sep=";"),file=out.file.eve,append=T)
# Figure: eve
xx[gridmask.notNAs]<-round(100*xidi.eve.FG,1)
rnc <- writeRaster(xx,
                   filename=out.file.grd.idi,
                   format="CDF", overwrite=TRUE)
# Figure: Analysis on high-resolution grid
xx[gridmask.notNAs]<-round(xa.FG,1)
if (flag.write.xa) {
  nogrid.ncout(grid=t(as.matrix(xx)),
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
#
quit(status=0)
