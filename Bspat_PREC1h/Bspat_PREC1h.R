# --~- SpInt_daily2monthly_v01.R  -~--
# Compute monthly average from daily averages of temperature
# 30.05.2014 - Cristian Lussana
# -----------------------------------------------------------------------------
rm(list=ls())
# Libraries
library(raster)
#------------------------------------------------------------------------------
dir.in.RR.main<-"/lustre/mnt/cristianl/seNorge2/PREC_daily/gridded_dataset_Bs1.1"
dir.in.RR1.main<-"/lustre/mnt/cristianl/seNorge2/PREC_hourly/gridded_dataset"
dir.in.RR1.idi.main<-"/lustre/mnt/cristianl/seNorge2_addInfo/PREC_hourly/gridded_dataset"
#/lustre/mnt/cristianl/seNorge2_addInfo/PREC_hourly/gridded_dataset/201407/seNorge_v2test_PREC_grid_normIDI_2014070715_2014070715.nc
# Read command line arguments
arguments <- commandArgs()
arguments
yyyy.mm.dd<-arguments[3]
config_file<-arguments[4]
config_par<-arguments[5]
if (length(arguments)!=5) 
  ext<-error_exit(paste("Error in command line arguments: \n",
  " R --vanilla yyyy.mm.dd configFILE configPAR \n",sep=""))
# [] define/check paths
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
#
dir.in.RR<-paste(dir.in.RR.main,"/",yyyymm,sep="")
file.RR<-paste(dir.in.RR,"/seNorge_v2test_PREC_grid_",yyyymmdd,"_",yyyymmdd,".nc",sep="")
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
# input directories
# daily precipitation data
path2input.1d.main<-paste(main.path.output,"/seNorge2/PREC1d",sep="")
path2input.1d.main.stn<-paste(path2output.main,"/station_dataset",sep="")
path2input.1d.main.grd<-paste(path2output.main,"/gridded_dataset",sep="")
path2input.1d.add<-paste(main.path.output,"/seNorge2_addInfo/PREC1d",sep="")
path2input.1d.add.grd<-paste(path2output.add,"/gridded_dataset",sep="")
path2input.1d.add.eve<-paste(path2output.add,"/event_dataset",sep="")
in.1d.file.stn<- paste(path2input.1d.main.stn,"/",yyyymm,
                       "/seNorge_v2_0_PREC1d_station_",
                       yyyymmdd,"_",yyyymmdd,".txt",sep="")
in.1d.file.eve<- paste(path2input.1d.add.eve,"/",yyyymm.b,
                       "/seNorge_v2_0_PREC1d_event_",
                       yyyymmdd,"_",yyyymmdd,".txt",sep="")
in.1d.file.grd.ana<- paste(path2input.1d.main.grd,"/",yyyymm.b,
                           "/seNorge_v2_0_PREC1d_grid_",
                           yyyymmdd,"_",yyyymmdd,".nc",sep="")
in.1d.file.grd.idi<- paste(path2input.1d.add.grd,"/",yyyymm.b,
                           "/seNorge_v2_0_PREC1d_grid_normIDI_",
                           yyyymmdd,"_",yyyymmdd,".nc",sep="")
# hourly precipitation data
path2input.1h.main<-paste(main.path.output,"/seNorge2/PREC1h",sep="")
path2input.1h.main.stn<-paste(path2output.main,"/station_dataset",sep="")
path2input.1h.main.grd<-paste(path2output.main,"/gridded_dataset",sep="")
path2input.1h.add<-paste(main.path.output,"/seNorge2_addInfo/PREC1h",sep="")
path2input.1h.add.grd<-paste(path2output.add,"/gridded_dataset",sep="")
path2input.1h.add.eve<-paste(path2output.add,"/event_dataset",sep="")
#in.1h.file.stn<- paste(path2input.1h.main.stn,"/",yyyymm,
#                       "/seNorge_v2_0_PREC1h_station_",
#                       yyyymmdd,"_",yyyymmdd,".txt",sep="")
#in.1h.file.eve<- paste(path2input.1h.add.eve,"/",yyyymm.b,
#                       "/seNorge_v2_0_PREC1h_event_",
#                       yyyymmdd,"_",yyyymmdd,".txt",sep="")
#in.1h.file.grd.ana<- paste(path2input.1h.main.grd,"/",yyyymm.b,
#                           "/seNorge_v2_0_PREC1h_grid_",
#                           yyyymmdd,"_",yyyymmdd,".nc",sep="")
#in.1h.file.grd.idi<- paste(path2input.1h.add.grd,"/",yyyymm.b,
#                           "/seNorge_v2_0_PREC1h_grid_normIDI_",
#                           yyyymmdd,"_",yyyymmdd,".nc",sep="")
in.1h.files.grd.ana<-vector()
in.1h.files.grd.ana<-paste(path2input.1h.main.grd,"/",yyyymm.v,
                           "/seNorge_v2_0_PREC1h_grid_",
                           yyyymmddhh.v,"_",yyyymmddhh.v,".nc",sep="")
files.check<-file.exists(in.1h.files.grd.ana)
i<-0
for (file.RR1 in files.RR1) {
  i<-i+1
  print(file.RR1)
  xa.RR1.i<-raster(file.RR1)
  projection(xa.RR1.i)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
  val.RR1<-extract(xa.RR1.i,1:ncell(xa.RR1.i))
#  indx1<-which(val.RR1<0.1)
  indx2<-which(is.na(val.RR1))
  print("before")
#  xa.RR1.filt<-focal(xa.RR1,w=gf)
  xa.RR1.i.filt<-focal(xa.RR1.i,w=matrix(1,nc=11,nr=11),fun=mean,na.rm=T)
  xa.RR1.i.filt[indx2]<-NA
#  xa.RR1.filt[indx1]<-0
  projection(xa.RR1.i.filt)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
  print("after")
  rnc<-writeRaster(xa.RR1.i.filt,filename=paste("./seNorge_v2test_filt",yyyymmddhh.v[i],".nc",sep=""),format="CDF",overwrite=TRUE)
  if (i==1) {
    xa.RR1.sum<-xa.RR1.i.filt
    projection(xa.RR1.sum)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
    xa.RR1<-xa.RR1.i.filt
    xa.RR1.orig<-xa.RR1.i
  } else {
    xa.RR1.sum<-xa.RR1.sum+xa.RR1.i.filt
    xa.RR1<-stack(xa.RR1,xa.RR1.i.filt)
    xa.RR1.orig<-stack(xa.RR1.orig,xa.RR1.i)
  }
}




# hourly precipitation data
# output directories
dir.create(file.path(main.path.output,"seNorge2"), showWarnings = FALSE)
dir.create(file.path(main.path.output,"seNorge2_addInfo"), showWarnings = FALSE)
path2output.main<-paste(main.path.output,"/seNorge2/PREC1h",sep="")
path2output.main.stn<-paste(path2output.main,"/station_dataset",sep="")
path2output.main.grd<-paste(path2output.main,"/gridded_dataset",sep="")
path2output.add<-paste(main.path.output,"/seNorge2_addInfo/PREC1h",sep="")
path2output.add.grd<-paste(path2output.add,"/gridded_dataset",sep="")
path2output.add.eve<-paste(path2output.add,"/event_dataset",sep="")
if (!(file.exists(path2output.main)))     dir.create(path2output.main,showWarnings=F) 
if (!(file.exists(path2output.main.stn))) dir.create(path2output.main.stn,showWarnings=F) 
if (!(file.exists(path2output.main.grd))) dir.create(path2output.main.grd,showWarnings=F) 
if (!(file.exists(path2output.add)))      dir.create(path2output.add,showWarnings=F) 
if (!(file.exists(path2output.add.grd)))  dir.create(path2output.add.grd,showWarnings=F) 
if (!(file.exists(path2output.add.eve)))  dir.create(path2output.add.eve,showWarnings=F) 
# Setup output files 
dir.create(paste(path2output.main.stn,"/",yyyymm.b,sep=""),showWarnings=F)
dir.create(paste(path2output.main.grd,"/",yyyymm.b,sep=""),showWarnings=F)
dir.create(paste(path2output.add.grd,"/",yyyymm.b,sep=""),showWarnings=F)
dir.create(paste(path2output.add.eve,"/",yyyymm.b,sep=""),showWarnings=F)
out.file.stn<- paste(path2output.main.stn,"/",yyyymm.b,
                     "/seNorge_v2_0_PREC1h_station_",
                     yyyymmdd.b,"_",yyyymmdd.e,".txt",sep="")
out.file.eve<- paste(path2output.add.eve,"/",yyyymm.b,
                     "/seNorge_v2_0_PREC1h_event_",
                     yyyymmdd.b,"_",yyyymmdd.e,".txt",sep="")
out.file.grd.ana<- paste(path2output.main.grd,"/",yyyymm.b,
                         "/seNorge_v2_0_PREC1h_grid_",
                         yyyymmdd.b,"_",yyyymmdd.e,".nc",sep="")
out.file.grd.idi<- paste(path2output.add.grd,"/",yyyymm.b,
                         "/seNorge_v2_0_PREC1h_grid_normIDI_",
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






#
print(files.RR1)
print(files.RR1.check)

xa.RR<-raster(file.RR)
projection(xa.RR)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
val.RR<-extract(xa.RR,1:ncell(xa.RR))
indx.RR.noNA<-which(!is.na(val.RR))
indx.RR.no0<-which(val.RR>=0.1)
#
i<-0
for (file.RR1 in files.RR1) {
  i<-i+1
  print(file.RR1)
  xa.RR1.i<-raster(file.RR1)
  projection(xa.RR1.i)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
  val.RR1<-extract(xa.RR1.i,1:ncell(xa.RR1.i))
#  indx1<-which(val.RR1<0.1)
  indx2<-which(is.na(val.RR1))
  print("before")
#  xa.RR1.filt<-focal(xa.RR1,w=gf)
  xa.RR1.i.filt<-focal(xa.RR1.i,w=matrix(1,nc=11,nr=11),fun=mean,na.rm=T)
  xa.RR1.i.filt[indx2]<-NA
#  xa.RR1.filt[indx1]<-0
  projection(xa.RR1.i.filt)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
  print("after")
  rnc<-writeRaster(xa.RR1.i.filt,filename=paste("./seNorge_v2test_filt",yyyymmddhh.v[i],".nc",sep=""),format="CDF",overwrite=TRUE)
  if (i==1) {
    xa.RR1.sum<-xa.RR1.i.filt
    projection(xa.RR1.sum)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
    xa.RR1<-xa.RR1.i.filt
    xa.RR1.orig<-xa.RR1.i
  } else {
    xa.RR1.sum<-xa.RR1.sum+xa.RR1.i.filt
    xa.RR1<-stack(xa.RR1,xa.RR1.i.filt)
    xa.RR1.orig<-stack(xa.RR1.orig,xa.RR1.i)
  }
}
#
#xa.RR1<-lapply(seq(files.RR1),function(i) {raster(files.RR1[i]) })
#xa.RR1.stack<-stack(xa.RR1)
#projection(xa.RR1.stack)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
#rg.xa.RR1<-focal(xa.RR1.stack,w=gf)
#print(rg.xa.RR1)
#xa.RR1.sum<-sum(xa.RR1.stack)

#val.RR1.sum<-extract(xa.RR1.sum,1:ncell(xa.RR1.sum))
#indx.weird<-which(val.RR>=0.1 & val.RR1.sum<0.1)
#val.RR1.sum[indx.weird]<-(-1)
#
#xa.RR1.idi<-lapply(seq(files.RR1.idi),function(i) {raster(files.RR1.idi[i]) })
#xa.RR1.idi.stack<-stack(xa.RR1.idi)
#projection(xa.RR1.idi.stack)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
#xa.RR1.idi.sum<-sum(xa.RR1.idi.stack)
#val.RR1.idi.sum<-extract(xa.RR1.idi.sum,1:ncell(xa.RR1.idi.sum))
#
rnc<-writeRaster(xa.RR,filename=paste("./seNorge_v2test_grid_RRsum.nc",sep=""),format="CDF",overwrite=TRUE)
rnc<-writeRaster(xa.RR1.sum,filename=paste("./seNorge_v2test_grid_RR1sum.nc",sep=""),format="CDF",overwrite=TRUE)
xa.RR1.dis<-xa.RR*xa.RR1/xa.RR1.sum
rnc<-writeRaster(xa.RR1.dis,filename=paste("./seNorge_v2test_grid_RR1d.nc",sep=""),format="CDF",overwrite=TRUE)
rnc<-writeRaster(xa.RR1.orig,filename=paste("./seNorge_v2test_grid_RR1o.nc",sep=""),format="CDF",overwrite=TRUE)
q()
#rnc<-writeRaster(xa.RR1.idi.sum,filename=paste("./seNorge_v2test_grid_RR1idisum.nc",sep=""),format="CDF",overwrite=TRUE)
rm(xa.RR1,xa.RR1.stack)
#val.RR1.dis<-vector(length=ncell(xa.RR))
#val.RR1.dis[]<-NA
i<-0
for (file.RR1 in files.RR1) {
  i<-i+1
  print(file.RR1)
  xa.RR1<-raster(file.RR1)
  val.RR1<-extract(xa.RR1,1:ncell(xa.RR1))
  val.RR1.dis[indx.RR.noNA]<-0
  val.RR1.dis[indx.RR.no0]<-val.RR[indx.RR.no0] * val.RR1[indx.RR.no0] / val.RR1.sum[indx.RR.no0]
  val.RR1.dis[indx.weird]<-val.RR[indx.weird] / 24
  xa.RR1[indx.RR.noNA]<-round(val.RR1.dis[indx.RR.noNA],2)
  rnc<-writeRaster(xa.RR1,filename=paste("./seNorge_v2test_grid_",yyyymmddhh.v[i],".nc",sep=""),format="CDF",overwrite=TRUE)
}

q()

#dir.in.stat<-"/vol/klimagrid/test/temp_cristianl/SpInt_TEMP_daily/stationpoints"
i<-0
fnames.tot<-vector()
dir<-paste(dir.in.grid,"/",yyyymm,"/",sep="")
fname.aux<-paste("seNorge_v2test_TEMP_grid_",yyyymm,"*.nc",sep="")
fnames.xa<-list.files(dir,pattern=glob2rx(fname.aux),recursive=T)
for (f in 1:length(fnames.xa)) {
  i<-i+1
  fnames.tot[i]<-paste(dir,fnames.xa[f],sep="")
}
print(fnames.tot) 
x<-raster(fnames.tot[1])
# read input data into raster format
xa.day<-lapply(seq(fnames.tot),function(i) {raster(fnames.tot[i]) })
print(xa.day) 
xa.day.stack<-stack(xa.day)
projection(xa.day.stack)<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
# 
xa.month<-mean(xa.day.stack)
rnc<-writeRaster(xa.month,filename=paste(dir.out,"/seNorge_v2test_TEMPdaily_grid_",yyyymm,".nc",sep=""),format="CDF",overwrite=TRUE)
#   
quit(status=0)
