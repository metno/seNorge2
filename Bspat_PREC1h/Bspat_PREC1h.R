# --~- SpInt_daily2monthly_v01.R  -~--
# Compute monthly average from daily averages of temperature
# 30.05.2014 - Cristian Lussana
# -----------------------------------------------------------------------------
rm(list=ls())
# Libraries
library(raster)
library(ncdf)
#
proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
#------------------------------------------------------------------------------
# Read command line arguments
arguments <- commandArgs()
arguments
yyyy.mm.dd<-arguments[3]
config_file<-arguments[4]
config_par<-arguments[5]
if (length(arguments)!=5) 
  ext<-error_exit(paste("Error in command line arguments: \n",
  " R --vanilla yyyy.mm.dd configFILE configPAR \n",sep=""))
# define/check paths
if (!file.exists(config_file)) 
  ext<-error_exit("Fatal Error: configuration file not found")
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
  if (file.exists(paste(main.path,"/Bspat_PREC1h/testbed",sep=""))) {
    testbed<-paste(main.path,"/Bspat_PREC1h/testbed",sep="")
    station.info<-paste(testbed,"/station_data.csv",sep="")
    observed.data<-paste(testbed,"/observed_data.csv",sep="")
  } else {
    ext<-error_exit(paste("testbed not found"))
  }
}
# netcdf fixed parameters
grid.type <- "utm33"
prod.date <- substr(Sys.time(),1,10)
xa.source.nc<-"hourly precipitation from station data"
xa.var.version <- "1.0"
xa.pname<-"PREC1h"
for (p in 1:length(ncout.spec.list)) {
  if (ncout.spec.list[[p]]$pname==xa.pname) break
}
if (p==length(ncout.spec.list) & (ncout.spec.list[[p]]$pname!=xa.pname) ) {
  xa.flag.write<-F
} else {
  xa.flag.write<-T
  xa.var.name<-ncout.spec.list[[p]]$opts$var.name
  xa.var.longname<-ncout.spec.list[[p]]$opts$var.longname
  xa.var.unit<-ncout.spec.list[[p]]$opts$var.unit
  xa.var.mv<-as.numeric(ncout.spec.list[[p]]$opts$var.mv)
  xa.times.unit <-ncout.spec.list[[p]]$opts$t.unit
  xa.times.ref <-ncout.spec.list[[p]]$opts$t.ref
  xa.reference <- ncout.spec.list[[p]]$opts$reference
}
# set Time-related variables
yyyy<-substr(yyyy.mm.dd,1,4)
mm<-substr(yyyy.mm.dd,6,7)
dd<-substr(yyyy.mm.dd,9,10)
yyyymm<-paste(yyyy,mm,sep="")
yyyymmdd<-paste(yyyy,mm,dd,sep="")
# set Time variables
end.string<-paste(yyyy.mm.dd,".06",sep="")
end<-strptime(end.string,"%Y.%m.%d.%H","UTC")
timeseq<-as.POSIXlt(seq(as.POSIXlt(end),length=24,by="-1 hour"),"UTC")
timeseq<-rev(timeseq)
yyyymm.v<-paste(formatC(timeseq$year+1900,width=4,flag="0"),
                formatC(timeseq$mon+1,width=2,flag="0"),sep="")
yyyymmddhh.v<-paste(formatC(timeseq$year+1900,width=4,flag="0"),
                    formatC(timeseq$mon+1,width=2,flag="0"),
                    formatC(timeseq$mday,width=2,flag="0"),
                    formatC(timeseq$hour,width=2,flag="0"),sep="")
print(timeseq)
print(yyyymmddhh.v)
nt<-length(yyyymmddhh.v)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# [] Grid
# CRS Coordinate Reference System
r.orog.FG<-raster(filenamedem)
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
# extract all the cell values: zvalues[1] contains the orog[1,1] value
# Raster: cell numbers start at 1 in the upper left corner,
# and increase from left to right, and then from top to bottom
zvalues.FG<-getValues(r.orog.FG)
storage.mode(zvalues.FG)<-"numeric"
xy.FG<-xyFromCell(r.orog.FG,1:ncell(r.orog.FG))
x.FG<-sort(unique(xy.FG[,1]))
y.FG<-sort(unique(xy.FG[,2]),decreasing=T)
mask.FG<-which(!is.na(zvalues.FG))
zgrid<-zvalues.FG[mask.FG]
xgrid<-xy.FG[mask.FG,1]
ygrid<-xy.FG[mask.FG,2]
Lgrid.FG<-length(mask.FG)
# clean memory
rm(zvalues.FG)
# debug info
print("grid parameters")
print(paste("nx ny dx dy",as.integer(nx.FG),as.integer(ny.FG),round(dx.FG,2),round(dy.FG,2)))
print(paste("xmn xmx ymn ymx",round(xmn.FG,2),round(xmx.FG,2),round(ymn.FG,2),round(ymx.FG,2)))
print(paste("Lgrid.FG",as.integer(Lgrid.FG)))
#..............................................................................
# input directories
# daily precipitation data
d.path.main<-paste(main.path.output,"/seNorge2/PREC1d",sep="")
d.path.main.stn<-paste(d.path.main,"/station_dataset",sep="")
d.path.main.grd<-paste(d.path.main,"/gridded_dataset",sep="")
d.path.add<-paste(main.path.output,"/seNorge2_addInfo/PREC1d",sep="")
d.path.add.grd<-paste(d.path.add,"/gridded_dataset",sep="")
d.path.add.eve<-paste(d.path.add,"/event_dataset",sep="")
d.file.stn<- paste(d.path.main.stn,"/",yyyymm,
                       "/seNorge_v2_0_PREC1d_station_",
                       yyyymmdd,"_",yyyymmdd,".txt",sep="")
d.file.eve<- paste(d.path.add.eve,"/",yyyymm,
                       "/seNorge_v2_0_PREC1d_event_",
                       yyyymmdd,"_",yyyymmdd,".txt",sep="")
d.file.grd<- paste(d.path.main.grd,"/",yyyymm,
                           "/seNorge_v2_0_PREC1d_grid_",
                           yyyymmdd,"_",yyyymmdd,".nc",sep="")
d.file.idi<- paste(d.path.add.grd,"/",yyyymm,
                           "/seNorge_v2_0_PREC1d_grid_normIDI_",
                           yyyymmdd,"_",yyyymmdd,".nc",sep="")
d.file.check<-file.exists(d.file.grd)
# hourly precipitation data
hrt.path.main<-paste(main.path.output,"/seNorge2/PREC1hRT",sep="")
hrt.path.main.stn<-paste(hrt.path.main,"/station_dataset",sep="")
hrt.path.main.grd<-paste(hrt.path.main,"/gridded_dataset",sep="")
hrt.path.add<-paste(main.path.output,"/seNorge2_addInfo/PREC1hRT",sep="")
hrt.path.add.grd<-paste(hrt.path.add,"/gridded_dataset",sep="")
hrt.path.add.eve<-paste(hrt.path.add,"/event_dataset",sep="")
hrt.files.grd<-vector()
hrt.files.grd<-paste(hrt.path.main.grd,"/",yyyymm.v,
                           "/seNorge_v2_0_PREC1hRT_grid_",
                           yyyymmddhh.v,"_",yyyymmddhh.v,".nc",sep="")
hrt.files.check<-file.exists(hrt.files.grd)
#..............................................................................
# output directories
# hourly precipitation data
dir.create(file.path(main.path.output,"seNorge2"), showWarnings = FALSE)
dir.create(file.path(main.path.output,"seNorge2_addInfo"), showWarnings = FALSE)
h.path.main<-paste(main.path.output,"/seNorge2/PREC1h",sep="")
h.path.main.stn<-paste(h.path.main,"/station_dataset",sep="")
h.path.main.grd<-paste(h.path.main,"/gridded_dataset",sep="")
h.path.add<-paste(main.path.output,"/seNorge2_addInfo/PREC1h",sep="")
h.path.add.grd<-paste(h.path.add,"/gridded_dataset",sep="")
h.path.add.eve<-paste(h.path.add,"/event_dataset",sep="")
if (!(file.exists(h.path.main)))     dir.create(h.path.main,showWarnings=F) 
if (!(file.exists(h.path.main.stn))) dir.create(h.path.main.stn,showWarnings=F) 
if (!(file.exists(h.path.main.grd))) dir.create(h.path.main.grd,showWarnings=F) 
if (!(file.exists(h.path.add)))      dir.create(h.path.add,showWarnings=F) 
if (!(file.exists(h.path.add.grd)))  dir.create(h.path.add.grd,showWarnings=F) 
if (!(file.exists(h.path.add.eve)))  dir.create(h.path.add.eve,showWarnings=F) 
# Setup output files
for (i in 1:length(yyyymm.v)) {
  dir.create(paste(h.path.main.stn,"/",yyyymm.v[i],sep=""),showWarnings=F)
  dir.create(paste(h.path.main.grd,"/",yyyymm.v[i],sep=""),showWarnings=F)
  dir.create(paste(h.path.add.grd,"/",yyyymm.v[i],sep=""),showWarnings=F)
  dir.create(paste(h.path.add.eve,"/",yyyymm.v[i],sep=""),showWarnings=F)
}
h.files.stn<- paste(h.path.main.stn,"/",yyyymm.v,
                    "/seNorge_v2_0_PREC1h_station_",
                    yyyymmddhh.v,"_",yyyymmddhh.v,".txt",sep="")
h.files.eve<- paste(h.path.add.eve,"/",yyyymm.v,
                    "/seNorge_v2_0_PREC1h_event_",
                    yyyymmddhh.v,"_",yyyymmddhh.v,".txt",sep="")
h.files.grd<- paste(h.path.main.grd,"/",yyyymm.v,
                    "/seNorge_v2_0_PREC1h_grid_",
                    yyyymmddhh.v,"_",yyyymmddhh.v,".nc",sep="")
h.files.idi<- paste(h.path.add.grd,"/",yyyymm.v,
                    "/seNorge_v2_0_PREC1h_grid_normIDI_",
                    yyyymmddhh.v,"_",yyyymmddhh.v,".nc",sep="")
#..............................................................................
print("Input files:")
print("PREC1d:")
print(cbind(d.file.grd,d.file.check))
print("PREC1hRT:")
print(cbind(hrt.files.grd,hrt.files.check))
print("PREC1h Output files:")
print("PREC1h analysis on the grid (netcdf)")
print(h.files.grd)
print("PREC1h event-normalized idi on the grid (netcdf)")
print(h.files.idi)
print("PREC1h station outputs (text)")
print(h.files.stn)
print("PREC1h event outputs (text)")
print(h.files.eve)
#..............................................................................
# daily precipitation: open/read/close netcdf file
nc<-open.ncdf(d.file.grd)
data<-get.var.ncdf(nc)
#aux<-att.get.ncdf(nc,"UTM_Zone_33","proj4")
#projstr<-aux$value
#dx<-nc$dim$X$vals[2]-nc$dim$X$vals[1]
#ex.xmin<-min(nc$dim$X$vals)-dx/2
#ex.xmax<-max(nc$dim$X$vals)+dx/2
#dy<-nc$dim$Y$vals[2]-nc$dim$Y$vals[1]
#ex.ymin<-min(nc$dim$Y$vals)-dy/2
#ex.ymax<-max(nc$dim$Y$vals)+dy/2
#nx<-nc$dim$X$len
#ny<-nc$dim$Y$len
close.ncdf(nc)
d.ra<-raster(ncol=nx.FG, nrow=ny.FG,
             xmn=xmn.FG, xmx=xmx.FG,
             ymn=ymn.FG, ymx=ymx.FG,
             crs=proj4.utm33)
d.ra[]<-NA
d.ra[]<-t(data)
#d.xa<-extract(d.ra,1:ncell(d.ra))
#d.indx.xa.noNA<-which(!is.na(d.xa))
#d.indx.xa.no0<-which(d.xa>=0.1)
#..............................................................................
# cycle: filter the sub-daily precipitation field; 
#        obtain the daily accumulated prec field from sub-daily prec fields
# create raster hrt.ra.i
hrt.ra.i<-d.ra
projection(hrt.ra.i)<-proj4.utm33
i<-0
for (hrt.file in hrt.files.grd) {
  i<-i+1
  print(hrt.file)
  hrt.ra.i[]<-NA
# sub-daily precipitation: open/read/close netcdf file
  nc<-open.ncdf(hrt.file)
  data<-get.var.ncdf(nc)
  close.ncdf(nc)
  hrt.ra.i[]<-t(data)
  hrt.ra.i.filt<-focal(hrt.ra.i,w=matrix(1,nc=11,nr=11),fun=mean,na.rm=T)
  projection(hrt.ra.i.filt)<-proj4.utm33
  if (i==1) {
    hrt.ra.sum<-hrt.ra.i.filt
    hrt.ra.filt<-hrt.ra.i.filt
    projection(hrt.ra.sum)<-proj4.utm33
    projection(hrt.ra.filt)<-proj4.utm33
#    hrt.ra.orig<-hrt.ra.i
  } else {
    hrt.ra.sum<-hrt.ra.sum+hrt.ra.i.filt
    hrt.ra.filt<-stack(hrt.ra.filt,hrt.ra.i.filt)
#    hrt.ra.orig<-stack(hrt.ra.orig,hrt.ra.i)
  }
}
# Disaggregate daily precipitation according to sub-daily fraction
h.ra<-d.ra*hrt.ra.filt/hrt.ra.sum
# Output Session
for (i in 1:nt) {
  b2r<-raster(h.ra,layer=i)
  projection(b2r)<-proj4.utm33
  h.xa<-extract(b2r,1:ncell(b2r))
  h.xa[mask.FG][is.na(h.xa[mask.FG])]<-0
  b2r[mask.FG]<-h.xa[mask.FG]
  nogrid.ncout(file.name=h.files.grd[i],
               grid=t(as.matrix(b2r)),
               x=x.FG,y=y.FG,grid.type=grid.type,
               times=c(paste(yyyymmddhh.v[i],"00",sep="")),
               prod.date=prod.date,
               proj4.string=proj4.utm33,
               var.name=xa.var.name,
               var.longname=xa.var.longname,
               var.unit=xa.var.unit,
               var.mv=xa.var.mv,
               var.version=xa.var.version,
               times.unit=xa.times.unit,
               times.ref=xa.times.ref,
               reference=xa.reference,
               source.string=xa.source.nc)
}
# Exit with success
quit(status=0)
