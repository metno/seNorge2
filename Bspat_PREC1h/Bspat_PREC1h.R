# --~- SpInt_daily2monthly_v01.R  -~--
# Compute monthly average from daily averages of temperature
# 30.05.2014 - Cristian Lussana
# -----------------------------------------------------------------------------
rm(list=ls())
# Libraries
library(raster)
library(ncdf)
library(cluster)
#-----------------------------------------------------------------------------
# + manage fatal error
error_exit<-function(str=NULL) {
  print("Fatal Error:")
  if (!is.null(str)) print(str)
  quit(status=1)
}
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
d.file.grd.check<-file.exists(d.file.grd)
d.file.stn.check<-file.exists(d.file.stn)
# hourly precipitation data
hrt.path.main<-paste(main.path.output,"/seNorge2/PREC1hRT",sep="")
hrt.path.main.stn<-paste(hrt.path.main,"/station_dataset",sep="")
hrt.path.main.grd<-paste(hrt.path.main,"/gridded_dataset",sep="")
hrt.path.add<-paste(main.path.output,"/seNorge2_addInfo/PREC1hRT",sep="")
hrt.path.add.grd<-paste(hrt.path.add,"/gridded_dataset",sep="")
hrt.path.add.eve<-paste(hrt.path.add,"/event_dataset",sep="")
hrt.files.grd<-vector()
hrt.files.stn<-vector()
hrt.files.grd<-paste(hrt.path.main.grd,"/",yyyymm.v,
                     "/seNorge_v2_0_PREC1hRT_grid_",
                     yyyymmddhh.v,"_",yyyymmddhh.v,".nc",sep="")
hrt.files.stn<-paste(hrt.path.main.stn,"/",yyyymm.v,
                     "/seNorge_v2_0_PREC1hRT_station_",
                     yyyymmddhh.v,"_",yyyymmddhh.v,".txt",sep="")
hrt.files.grd.check<-file.exists(hrt.files.grd)
hrt.files.stn.check<-file.exists(hrt.files.stn)
#..............................................................................
# check files
  if (!d.file.grd.check | !d.file.stn.check) {
    print("Fatal Error: Error reading daily precipitation files")
    print(cbind(d.file.grd,d.file.grd.check))
    print(cbind(d.file.stn,d.file.stn.check))
    error_exit()
  }
  if (any(!hrt.files.grd.check) | any(!hrt.files.stn.check)) {
    print("Fatal Error: Error reading sud-daily precipitation files")
    print(cbind(hrt.files.grd,hrt.files.grd.check))
    print(cbind(hrt.files.stn,hrt.files.stn.check))
    error_exit()
  }
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
print("PREC1d stations:")
print(cbind(d.file.stn,d.file.stn.check))
print("PREC1d grid:")
print(cbind(d.file.grd,d.file.grd.check))
print("PREC1hRT stations:")
print(cbind(hrt.files.stn,hrt.files.stn.check))
print("PREC1hRT grid:")
print(cbind(hrt.files.grd,hrt.files.grd.check))
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
# daily prec at station locations
#year,month,day,nday,stid,x,y,z,eve.lab,yo,yb,ya,yav,yidi,yidiv,dqcflag
d.y<-read.csv(file=d.file.stn,header=T,sep=";")
d.nstn<-length(d.y$stid)
print("daily prec at station locations")
print(d.y)
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
#
ya.aux<-extract(d.ra,cbind(d.y$x,d.y$y),na.rm=T)
indx<-which( is.na(d.y$ya) & (!is.na(ya.aux)))
if (length(indx)>0) {
  d.y$ya[indx]<-ya.aux[indx]
  d.y$dqcflag[indx]<-1000
}

#d.xa<-extract(d.ra,1:ncell(d.ra))
#d.indx.xa.noNA<-which(!is.na(d.xa))
#d.indx.xa.no0<-which(d.xa>=0.1)
#..............................................................................
# cycle: filter the sub-daily precipitation field; 
#        obtain the daily accumulated prec field from sub-daily prec fields
#
# create raster hrt.ra.i
hrt.ra.i<-d.ra
projection(hrt.ra.i)<-proj4.utm33
hrt.ra.i[]<-NA
#
hrt.yasum<-matrix(nrow=d.nstn)
hrt.yavsum<-matrix(nrow=d.nstn)
hrt.yo<-matrix(ncol=nt,nrow=d.nstn)
hrt.ya<-matrix(ncol=nt,nrow=d.nstn)
hrt.yav<-matrix(ncol=nt,nrow=d.nstn)
hrt.yidi<-matrix(ncol=nt,nrow=d.nstn)
hrt.yidiv<-matrix(ncol=nt,nrow=d.nstn)
hrt.eve.lab<-matrix(ncol=nt,nrow=d.nstn)
hrt.dqcflag<-matrix(ncol=nt,nrow=d.nstn)
hrt.yasum[]<-0
hrt.yavsum[]<-0
hrt.yo[]<-NA
hrt.ya[]<-NA
hrt.yav[]<-NA
hrt.yidi[]<-NA
hrt.yidiv[]<-NA
hrt.eve.lab[]<-NA
hrt.dqcflag[]<-NA
for (i in 1:nt) {
  # sub-daily prec grid 
  hrt.file<-hrt.files.grd[i]
  # sub-daily precipitation: open/read/close netcdf file
  nc<-open.ncdf(hrt.file)
  data<-get.var.ncdf(nc)
  close.ncdf(nc)
  hrt.ra.i[]<-NA
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
  # sub-daily prec station locations
  #year,month,day,hour,nhour,stid,x,y,z,eve.lab,yo,yb,ya,yav,yidi,yidiv,dqcflag
  hrt.y<-read.csv(file=hrt.files.stn[i],header=T,sep=";")
  print("sub-daily prec at station locations")
  print(hrt.y)
  print(hrt.file)
  match<-match(d.y$stid,hrt.y$stid)
  hrt.yo[,i]<-hrt.y$yo[match]
  hrt.ya[,i]<-hrt.y$ya[match]
  hrt.yav[,i]<-hrt.y$yav[match]
  hrt.yidi[,i]<-hrt.y$yidi[match]
  hrt.yidiv[,i]<-hrt.y$yidiv[match]
  hrt.dqcflag[,i]<-hrt.y$dqcflag[match]
  print(cbind(d.y$stid,d.y$yo,d.y$ya,d.y$dqcflag,hrt.yo[,i],hrt.ya[,i],hrt.dqcflag[,i]))
  ya.aux<-extract(hrt.ra.i,cbind(hrt.y$x,hrt.y$y),na.rm=T)
  indx<-which( is.na(hrt.ya[,i]) & (!is.na(ya.aux)) )
  if (length(indx)>0) {
    hrt.ya[indx,i]<-ya.aux[indx]
    hrt.yav[indx,i]<-ya.aux[indx]
    hrt.dqcflag[indx,i]<-1000
  }
  hrt.yasum<-hrt.yasum+hrt.ya
  hrt.yavsum<-hrt.yavsum+hrt.yav
  if (i==2) break
}
#..............................................................................
# Disaggregate daily precipitation according to sub-daily fraction
# grid
h.ra<- d.ra * hrt.ra.filt/hrt.ra.sum
# stations
h.ya<-matrix(ncol=nt,nrow=d.nstn)
h.yav<-matrix(ncol=nt,nrow=d.nstn)
h.ya[]<-NA
h.yav[]<-NA
for (i in 1:nt) {
  h.ya[,i]<- d.y$ya * hrt.ya[,i]/hrt.yasum
  indx<-which(d.y$ya==0 | hrt.yasum==0 | hrt.ya[,i]==0)
  h.ya[indx,i]<-0.
  h.yav[,i]<- d.y$yav * hrt.yav[,i]/hrt.yavsum
  indx<-which(d.y$yav==0 | hrt.yavsum==0 | hrt.yav[,i]==0)
  h.yav[indx,i]<-0.
}
#..............................................................................
# Output Session
for (i in 1:nt) {
  # sub-daily prec
  #year,month,day,hour,nhour,stid,x,y,z,eve.lab,yo,yb,ya,yav,yidi,yidiv,dqcflag
  hrt.y<-read.csv(file=hrt.files.stn[i],header=T,sep=";")
  print("sub-daily prec at station locations")
  print(hrt.y)
  #
  b2r<-raster(h.ra,layer=i)
  projection(b2r)<-proj4.utm33
  h.xa<-extract(b2r,1:ncell(b2r))
  h.xa[mask.FG][is.na(h.xa[mask.FG])]<-0
  b2r[mask.FG]<-h.xa[mask.FG]
  # Detect clumps (patches) of connected cells. Each clump gets a unique ID.
  h.reve<-clump(b2r,directions=8,gaps=F)
  h.xeve<-getValues(h.reve)[mask.FG]
  mask.clu.FG<-which(!is.na(h.xeve))
  NAmask.clu.FG<-which(is.na(h.xeve))
  f.lab<-freq(h.reve)
  aux<-which(!is.na(f.lab[,1]))
  eve.labels<-f.lab[aux,1]
  eve.labels.n<-f.lab[aux,2]
  n.eve<-length(eve.labels)
  print(n.eve)
  print(cbind(eve.labels,eve.labels.n))
  rnc<-writeRaster(h.reve,filename=paste("clump.nc",sep=""),format="CDF",overwrite=TRUE)
  # ell -> ellipse
  ell.locx.eve<-vector(mode="numeric",length=n.eve)
  ell.locy.eve<-vector(mode="numeric",length=n.eve)
  ell.smajor.eve<-vector(mode="numeric",length=n.eve)
  ell.sminor.eve<-vector(mode="numeric",length=n.eve)
  ell.smadir.eve<-vector(mode="numeric",length=n.eve)
#------------------------------------------------------------------------------
# ELLIPSOID HULLS  ELLIPSOID HULLS  ELLIPSOID HULLS  ELLIPSOID HULLS  ELLIPSOID HULLS 
# + Ellipsoid hulls
# ellipsoid hulls computation could fail (i.e. small events) and return NAs. 
  # debug: start
#  png(file="../../seNorge2_scratch/Bspat_PREC1hRT/elli.png",width=1200,height=1200)
#  mx<-as.integer(max(eve.labels,na.rm=T))
#  cols<-c("gray",rainbow((mx-1)))
#  plot(r.eve.FG,breaks=seq(0.5,(mx+0.5),by=1),col=cols)
  # debug: end
  for (n in 1:n.eve) {
#    print(paste("n eve.labels[n]",n,eve.labels[n]))
    xindx.eve<-which(h.xeve==eve.labels[n])
    Lgrid.eve<-length(xindx.eve)
#    print(Lgrid.eve)
    min.Dh.seq.allowed<-10
    if (Lgrid.eve<=(min.Dh.seq.allowed)**2.) {
      ell.smajor.eve[n] <- 2*min.Dh.seq.allowed
      ell.sminor.eve[n] <- 2*min.Dh.seq.allowed
      ell.smadir.eve[n] <- NA
      next
    }
    xy<-cbind(xgrid[xindx.eve],ygrid[xindx.eve])
    ell<-ellipsoidhull(xy)
    ell.locx.eve[n]<-ell$loc[1]
    ell.locy.eve[n]<-ell$loc[2]
    if ( is.na(ell$cov[1,1]) | is.na(ell$cov[2,2]) | is.na(ell$cov[1,2]) | 
         ell$cov[1,1]==0 | ell$cov[2,2]==0 | ell$cov[1,2]==0) {
      ell.smajor.eve[n] <- 2*min.Dh.seq.allowed
      ell.sminor.eve[n] <- 2*min.Dh.seq.allowed
      ell.smadir.eve[n] <- NA
      next
    }
#   debug plot: lines(predict(ell), col="blue")
#    lines(predict(ell), col="blue")
#    points(ell.locx.eve[n],ell.locy.eve[n],pch=19,col="black",cex=1.5)
#    eigen -> values sorted in decreasing order; 
#             vectors matrix whose columns contain the eigenvectors
    eigenval<-eigen(ell$cov)$values
    eigenvec<-eigen(ell$cov)$vectors
    e <- sqrt(eigenval)
    ell.smajor.eve[n] <- sqrt(ell$d2) * e[1] /1000  # semi-major axis [Km]
    ell.sminor.eve[n] <- sqrt(ell$d2) * e[2] /1000  # semi-minor axis [Km]
    # orientation of the ellipse: angle between the major axis and the y-axis
    # dir=0 N-S orientation; dir=45 NE-SW; dir=90 E-W; dir=135 NW-SE; dir=180 N-S
    ell.smadir.eve[n]<-atan2(eigenvec[1,1],eigenvec[2,1])
    if (ell.smadir.eve[n]<0) ell.smadir.eve[n]<-ell.smadir.eve[n]+pi
    ell.smadir.eve[n]<-ell.smadir.eve[n]/pi*180.
#    print(paste(eve.labels[n]))
#    print(paste(ell.locx.eve[n],ell.locy.eve[n]))
#    print(paste(ell.smajor.eve[n],ell.sminor.eve[n]))
#    print(paste(ell.smadir.eve[n]))
  } # end cycle: compute ellipsoid hulls
# debug: close plot session
#  dev.off()
  if (exists("xindx.eve")) rm(xindx.eve)
  if (exists("eigenval")) rm(xy,ell,eigenval,eigenvec,e)

# 
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
  if (i==2) break
}
# Exit with success
quit(status=0)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# HEADER HEADER HEADER HEADER HEADER HEADER HEADER HEADER
# define header for the station data output file
cat(paste("year","month","day","hour","nhour","stid",
          "x","y","z","eve.lab","yo",
          "yb","ya","yav","yidi","yidiv","dqcflag","\n",sep=";"),
          file=out.file.stn,append=F)
cat(paste("year","month","day","hour","nhour","eve.lab","nobs",
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


print("++ Output")
# Station Points - Write output on file 
cat(paste(yyyy.b,mm.b,dd.b,hh.b,nhour,
          formatC(VecS[stn.output],format="f",digits=0),
          formatC(VecX[stn.output],format="f",digits=0),
          formatC(VecY[stn.output],format="f",digits=0),
          formatC(VecZ[stn.output],format="f",digits=0),
          formatC(y.eve[stn.output],format="f",digits=0),
          formatC(yo[stn.output],format="f",digits=1),
          formatC(yb[stn.output],format="f",digits=2),
          formatC(ya[stn.output],format="f",digits=2),
          formatC(yav[stn.output],format="f",digits=2),
          formatC(yidi.eve[stn.output],format="f",digits=4),
          formatC(yidiv.eve[stn.output],format="f",digits=4),
          formatC(ydqc.flag[stn.output],format="f",digits=0),
          "\n",sep=";"),file=out.file.stn,append=T)
cat(paste(yyyy.b,mm.b,dd.b,hh.b,nhour,
          eve.labels,
          n.y.eve,
          formatC(area.eve,format="f",digits=0),
          formatC(volume.eve,format="f",digits=2),
          formatC(meanidi.x.eve,format="f",digits=4),
          formatC(meanidi.y.eve,format="f",digits=4),
          formatC(meanidiv.y.eve,format="f",digits=4),
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
          formatC(meanidiv.y.eve.q50,format="f",digits=4),
          formatC(n.q50,format="f",digits=0),
          formatC(cv.rel.eve.q75,format="f",digits=4),
          formatC(cv.bias.eve.q75,format="f",digits=4),
          formatC(cv.rmse.eve.q75,format="f",digits=4),
          formatC(cv.made.eve.q75,format="f",digits=4),
          formatC(meanidiv.y.eve.q75,format="f",digits=4),
          formatC(n.q75,format="f",digits=0),
          formatC(idi.norm.fac,format="f",digits=5),
          "\n",sep=";"),file=out.file.eve,append=T)

    #
    idi.norm.fac[n]<-max(xidi.FG[xindx.eve.FG],yidi.eve[yindx],yidiv.eve[yindx])
    xidi.FG[xindx.eve.FG]<-xidi.FG[xindx.eve.FG]/idi.norm.fac[n]
    yidi.eve[yindx]<-yidi.eve[yindx]/idi.norm.fac[n]
    yidiv.eve[yindx]<-yidiv.eve[yindx]/idi.norm.fac[n]
    if (n.ya.eve[n]>0) {
      yidi.eve[ya.indx]<-yidi.eve[ya.indx]/idi.norm.fac[n]
      yidiv.eve[ya.indx]<-yidiv.eve[ya.indx]/idi.norm.fac[n]
    }
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
    q50.pos<-which(yo.aux>=q50)
    q75.pos<-which(yo.aux>=q75)
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

