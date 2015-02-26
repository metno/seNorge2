# --~- Bspat_TEMP1h_v1_1.R  -~--
# Bayesian spatial interpolation of TA.
# TA = istantaneous air temperature (sampled with hourly frequency)
# Spatial Consistency Test (SCT) is included.
#
# Outputs: look for "@@@@@@@@@" in the code and you'll get the ouput formats
#
# History:
# 26.02.2015 - Cristian Lussana. original code from Bspat_TAMRR_v1_0.R
#  change log: 
#  - allow for the use of observations outside Norway
#  - revisied queries to KDVH
#  - geographical information on seNorge2_dem_UTM33.nc
#  - definition of TEMP1h
#  - definition of a new directory tree
# ==~==========================================================================
rm(list=ls())
# Libraries
library(raster)
library(rgdal)
library(ncdf)
# paths
#[DEVELOPMENT]
main.path<-"/disk1/projects/seNorge2"
main.path.output<-"/disk1"
#
main.path.prog<-paste(main.path,"/Bspat_TEMP1h",sep="")
main.path.geoinfo<-paste(main.path,"/geoinfo",sep="")
# common libs and etcetera
path2lib.com<-paste(main.path,"/lib",sep="")
path2etc.com<-paste(main.path,"/etc",sep="")
#
path2output.main<-paste(main.path.output,"/seNorge2/TEMP1h",sep="")
path2output.main.stn<-paste(path2output.main,"/station_dataset",sep="")
path2output.main.grd<-paste(path2output.main,"/gridded_dataset",sep="")
path2output.add<-paste(main.path.output,"/seNorge2_addInfo/TEMP1h",sep="")
path2output.add.grd<-paste(path2output.add,"/gridded_dataset",sep="")
# External Functions
source(paste(path2lib.com,"/SpInt_PseudoBackground.R",sep=""))
source(paste(path2lib.com,"/nogrid.ncout.R",sep=""))
source(paste(path2lib.com,"/ncout.spec.list.r",sep=""))
# Read Geographical Information
filenamedem<-paste(main.path.geoinfo,"/seNorge2_dem_UTM33.nc",sep="")
# CRS strings
proj4.wgs84<-"+proj=longlat +datum=WGS84"
proj4.ETRS_LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
#
max.Km.stnINdomain<-50
#-------------------------------------------------------------------
# [] Setup OI parameters
sig2o<-3.0
eps2<-0.5
Dh<-60
Dz<-600
T2<-20
#
print("ANALYSIS and DQC parameters")
print(paste("EPS2 Dh[Km] Dz[m] > ", round(eps2,3),
                                    round(Dh,3),
                                    round(Dz,3),
                                    sep=" "))
print(paste("VARobs[C^2] T^2 = ", round(sig2o,3), round(T2,3),sep=" "))
# Background - parameters
Dh.b<-70.
eps2.b<-0.5
Lsubsample<-50
Lsubsample.max<-50
Lsubsample.DHmax<-200
Lsubsample.vec<-vector()
# netcdf fixed parameters
grid.type <- "utm33"
source.nc<-"air temperature from station data"
var.version.xa <- "1.0"
prod.date <- substr(Sys.time(),1,10)
pname.xa<-"TEMP1h"
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
# Read command line arguments
arguments <- commandArgs()
arguments
date.string<-arguments[3]
#1234567890123
#yyyy.mm.dd.hh
yyyy<-substr(date.string,1,4)
mm<-substr(date.string,6,7)
dd<-substr(date.string,9,10)
hh<-substr(date.string,12,13)
# useful only if nmt=1 in database query
#if (as.numeric(hh)==0) {
#  print(paste("Warning! Timestamp will be modified: from ",yyyy,".",mm,".",dd," ",hh,sep=""))
#  aux.date <- strptime(paste(yyyy,mm,dd,hh,sep="."),"%Y.%m.%d.%H","UTC")
#  date.minus.1h<-as.POSIXlt(seq(as.POSIXlt(aux.date),length.out=2,by="-1 hour"),"UTC")
#  aux.date<-date.minus.1h[2]
#  yyyy<-aux.date$year+1900
#  mm<-formatC(aux.date$mon+1,width=2,flag="0")
#  dd<-formatC(aux.date$mday,width=2,flag="0")
#  hh<-24
#  print(paste("to ",yyyy,".",mm,".",dd," ",hh," (UTC+1)",sep=""))
#}
date.dot<-paste(dd,".",mm,".",yyyy,sep="")
yyyymm<-paste(yyyy,mm,sep="")
yyyymmdd<-paste(yyyymm,dd,sep="")
yyyymmddhh<-paste(yyyymmdd,hh,sep="")
h<-as.numeric(hh)
#
dir.create(paste(path2output.main.stn,"/",yyyymm,sep=""),showWarnings=F)
dir.create(paste(path2output.main.grd,"/",yyyymm,sep=""),showWarnings=F)
dir.create(paste(path2output.add.grd,"/",yyyymm,sep=""),showWarnings=F)
out.file.stn<- paste(path2output.main.stn,"/",yyyymm,
                "/seNorge_v2_0_TEMP1h_station_",yyyymmddhh,".txt",sep="")
out.file.grd.ana<- paste(path2output.main.grd,"/",yyyymm,
                   "/seNorge_v2_0_TEMP1h_grid_",yyyymmddhh,".nc",sep="")
out.file.grd.bck<- paste(path2output.add.grd,"/",yyyymm,
                   "/seNorge_v2_0_TEMP1h_grid_background_",yyyymmddhh,".nc",sep="")
out.file.grd.idi<- paste(path2output.add.grd,"/",yyyymm,
                   "/seNorge_v2_0_TEMP1h_grid_idi_",yyyymmddhh,".nc",sep="")
#
print("Output files:")
print("analysis on the grid (netcdf)")
print(out.file.grd.ana)
print("background on the grid (netcdf)")
print(out.file.grd.bck)
print("idi on the grid (netcdf)")
print(out.file.grd.idi)
print("station outputs (text)")
print(out.file.stn)
#------------------------------------------------------------------------------
# Grid - it is defined by the DEM file
# CRS Coordinate Reference System
stackGeoGrid<-raster(filenamedem)
nx<-ncol(stackGeoGrid)
ny<-nrow(stackGeoGrid)
dx<-xres(stackGeoGrid)
dy<-yres(stackGeoGrid)
# 4 borders point
xmn<-xmin(stackGeoGrid)
xmx<-xmax(stackGeoGrid)
ymn<-ymin(stackGeoGrid)
ymx<-ymax(stackGeoGrid)
# South-West Point Coordinates
#Xnodesw<-xmn+dx/2
#Ynodesw<-ymn+dy/2
Xnodesw<-xmn
Ynodesw<-ymn
Xnode<-Xnodesw+(0:nx-1)*dx
Ynode<-Ynodesw+(0:ny-1)*dy
# Extract orography and urban-weight on unmasked gridpoints only
orog<-stackGeoGrid
# extract all the cell values: cells[1] contains the orog[1,1] value
# Raster: cell numbers start at 1 in the upper left corner,
# and increase from left to right, and then from top to bottom
cells<-extract(orog,1:ncell(orog))
xy<-xyFromCell(orog,1:ncell(orog))
x.G<-sort(unique(xy[,1]))
y.G<-sort(unique(xy[,2]),decreasing=T)
rc<-rowColFromCell(orog,1:ncell(orog))
aux<-as.vector(cells)
ccgrid<-which(!is.na(aux))
zgrid<-aux[ccgrid]
xgrid<-xy[ccgrid,1]
ygrid<-xy[ccgrid,2]
rowgrid<-rc[ccgrid,1]
colgrid<-rc[ccgrid,2]
Lgrid<-length(xgrid)
print(Lgrid)
rm(xy,rc,rowgrid,colgrid,cells,aux,stackGeoGrid)
#------------------------------------------------------------------------------
# Read Station Information 
myurl <- paste("http://klapp.oslo.dnmi.no/metnopub/production/metno?",
               "re=16&nod=NA&ct=text/plain&ddel=dot&del=semicolon",
               "&fy=",yyyy,"&ty=",yyyy,sep="")
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
# first step: lat, lon and elevation must be present
lat_dec<-as.numeric(stataux$LAT_DEC)
lon_dec<-as.numeric(stataux$LON_DEC)
z<-as.numeric(stataux$AMSL)
#indx<-which( !is.na(lat_dec) & !is.na(lon_dec) & !is.na(z) & lat_dec>40 & lat_dec<85 & lon_dec>0 & lon_dec<60 )
indx<-which( !is.na(lat_dec) & !is.na(lon_dec) & !is.na(z) )
# second step: the location must be in Norway or on the border (lee than max.Km.stnINdomain)
#  intermediate step: transformation in Km-coordinates ETRS_LAEA, which has a transformation 
#    less problematic than UTM33
coord<-SpatialPoints(cbind(lon_dec[indx],lat_dec[indx]), proj4string=CRS(proj4.wgs84))
coord.new<-spTransform(coord, CRS(proj4.ETRS_LAEA))
xy.RR<-coordinates(coord.new)
x<-round(xy.RR[,1],0)
y<-round(xy.RR[,2],0)
stations.tmp<-data.frame(matrix(nrow=length(indx),ncol=7))
names(stations.tmp)<-c("stnr","dept_no","lat_dec","lon_dec","z","x","y")
stations.tmp$stnr<-as.numeric(stataux$STNR[indx])
stations.tmp$dept_no<-as.numeric(stataux$DEPT_NO[indx])
stations.tmp$z<-z[indx]
stations.tmp$y<-y
stations.tmp$x<-x
stations.tmp$lat_dec<-lat_dec[indx]
stations.tmp$lon_dec<-lon_dec[indx]
n.stn<-length(stations.tmp$stnr)
indx.no<-which(stations.tmp$dept_no>=1 & stations.tmp$dept_no<=20)
Disth<-matrix(ncol=n.stn,nrow=n.stn,data=0.)
Disth<-(outer(stations.tmp$y,stations.tmp$y,FUN="-")**2.+outer(stations.tmp$x,stations.tmp$x,FUN="-")**2.)**0.5/1000.
aux.in<-vector(length=n.stn)
aux.in[1:n.stn]<-F
for (s in indx.no) {
  aux.in[which(Disth[s,]<max.Km.stnINdomain)]<-T 
}
indx.in<-which(aux.in)
# final step: create the structure stations containing the selected stations in utm33 coordinates
# note: for each of the selected we will obtain an analysis value
LOBS<-length(indx.in)
stations<-data.frame(matrix(nrow=LOBS,ncol=4))
names(stations)<-c("stnr","z","x","y")
lat_dec<-stations.tmp$lat_dec[indx.in]
lon_dec<-stations.tmp$lon_dec[indx.in]
coord<-SpatialPoints(cbind(lon_dec,lat_dec), proj4string=CRS(proj4.wgs84))
coord.new<-spTransform(coord, CRS(proj4.utm33))
xy.RR<-coordinates(coord.new)
stations$x<-round(xy.RR[,1],0)
stations$y<-round(xy.RR[,2],0)
stations$z<-stations.tmp$z[indx.in]
stations$stnr<-stations.tmp$stnr[indx.in]
rm(lat_dec,lon_dec,z,indx,coord,coord.new,xy.RR,x,y)
rm(stations.tmp,n.stn,indx.no,Disth,aux.in,indx.in)
# OLD code
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
print(LOBS)
# define Vectors and Matrices
VecX<-vector(mode="numeric",length=LOBS)
VecY<-vector(mode="numeric",length=LOBS)
VecZ<-vector(mode="numeric",length=LOBS)
VecS<-vector(mode="numeric",length=LOBS)
Disth<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
Distz<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
yo<-vector(mode="numeric",length=LOBS)
yb<-vector(mode="numeric",length=LOBS)
ya<-vector(mode="numeric",length=LOBS)
yav<-vector(mode="numeric",length=LOBS)
yidi<-vector(mode="numeric",length=LOBS)
yidiv<-vector(mode="numeric",length=LOBS)
ydqc<-vector(mode="numeric",length=LOBS)
ydqc.flag<-vector(mode="numeric",length=LOBS)
D<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
S<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
VecX<-as.numeric(as.vector(stations$x))
VecY<-as.numeric(as.vector(stations$y))
VecZ<-as.numeric(as.vector(stations$z))
VecS<-as.numeric(as.vector(stations$stnr))
# compute S and D=S+R matrices (R is assumed to be diagonal = sigma_obs**2*I)
# Disth and Distz are the (symmetric) matrices where 
#      Disth(i,j)=horizontal distance between i-th station and j-th station [Km]
#      Distz(i,j)=elevation difference between i-th station and j-th station [m]
Disth<-(outer(VecY,VecY,FUN="-")**2.+outer(VecX,VecX,FUN="-")**2.)**0.5/1000.
Distz<-abs(outer(VecZ,VecZ,FUN="-"))
D<-exp(-0.5*(Disth/Dh)**2.-0.5*(Distz/Dz)**2.)
S<-D
D[row(D)==col(D)]<-D[row(D)==col(D)]+eps2
#
D.b<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
S.b<-matrix(ncol=LOBS,nrow=LOBS,data=0.)
# NOTE: for the background Dh.b replaces Dh and Dz.b replaces Dz
# CL: is it eps2.b needed??
#D.b<-exp(-0.5*(Disth/Dh.b)**2.-0.5*(Distz/Dz.b)**2.)
D.b<-exp(-0.5*(Disth/Dh.b)**2.)
S.b<-D.b
D.b[row(D.b)==col(D.b)]<-D.b[row(D.b)==col(D.b)]+eps2.b
#
# DEBUG
#print(stations)
#write.table(file="stations.txt",stations)
#stations<-read.table(file="stations.txt")
#print(stations[1,])
#stations<-read.table("stations.txt",header=TRUE)
#print(stations[1,])
#  stnr z      x       y fennomean fenno_min4 fenno_long fenno_lat
# number of days
#print(ndays)
#------------------------------------------------------------------------------
# Elaborations
# define header for the station data output file
cat(paste("year","month","day","hour","stid","x","y","z","yo",
          "yb","ya","yav","yidi","yidiv","dqcflag","\n",sep=";"),
          file=out.file.stn,append=F)
# xx is the raster structure used for map production
xx <-raster(ncol=nx, nrow=ny, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx,
            crs=proj4.utm33)
xx[]<-NA
#
print(paste("+ ",yyyymmddhh," -~-----------------",sep=""))
ulric<-paste("http://klapp.oslo.dnmi.no/metnopub/production/",
             "metno?re=17&p=TA&fd=",date.dot,"&td=",date.dot,"&h=",h,
             "&nob=0.0&ddel=dot&del=semicolon&nmt=0",
             "&ct=text/plain&split=1&nod=-999",sep="")
#print(ulric)
o.cont<-1
while (o.cont<=10) {
  o<-NULL
  try(o <- read.table(ulric, header = TRUE,  sep = ";", #nrows = nrows,
          stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
 	          encoding = "UTF-8", quote = "",na.string=-999))
# names(o) -->  stnr;Year;Month;Day;TA
#    print(o)
#    print(length(o))
  value<-as.numeric(o$TA)
  stnr<-as.numeric(o$Stnr)
  hour.seq<-as.integer(o$Time.UTC.)
  value[value==-999]<-NA
  indx<-which( !is.na(value) & (stnr %in% stations$stnr) & (hour.seq==h) )
  LOBSt<-length(indx)
  if (LOBSt<10) {
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
# note: LOBSt always greater than 10
print(paste("  Total number of observations [not NA] =",LOBSt))
# OBS: "d"-day observations without NAs
OBS<-data.frame(matrix(nrow=LOBSt,ncol=7))
names(OBS)<-c("stnr","year","month","day","hour","value","DQC")
OBS$stnr<-stnr[indx]
OBS$year<-o$Year[indx]
OBS$month<-o$Month[indx]
OBS$day<-o$Day[indx]
OBS$hour<-hour.seq[indx]
OBS$value<-value[indx]
rm(stnr,o,indx,value,hour.seq)
# DQC flag: -1 missing, 0 good obs, 1 erroneous obs
OBS$DQC<-rep(-1,LOBSt)
#  write.table(file="o.txt",o)
stn.match.obs<-match(stations$stnr,OBS$stnr)
obs.match.stn<-match(OBS$stnr,stations$stnr)
yo[]<-NA
yo[obs.match.stn]<-as.numeric(as.vector(OBS$value))
yo.h.pos<-which(!is.na(yo))
rm(OBS)
# BACKGROUND AT STATION LOCATIONS
# For each station, compute a non-linear vertical profile using 
# the Lsubsample (closest) surrounding stations.
# The non-linear profile allows for a lat-lon dependance and
# for an inversion in the vertical profile.
#
# for each iteration, the yb is computed using only Lsubsample stations
# but over all stations (also the distant ones). A weight is then assigned
# for all the stations, giving more weight to the stations closest to the 
# actually used stations. This procedure should ensure a continuous background
yb.set<-matrix(data=NA,ncol=LOBS,nrow=LOBS)
ybweights.set<-matrix(data=0,ncol=LOBS,nrow=LOBS)
yb.param<-matrix(data=NA,ncol=12,nrow=LOBS)
yb.param1.mat<-matrix(data=NA,ncol=12,nrow=LOBS)
yb.param2.mat<-matrix(data=NA,ncol=12,nrow=LOBS)
VecS.set<-matrix(data=NA,ncol=Lsubsample.max,nrow=LOBS)
VecS.set.pos<-matrix(data=NA,ncol=Lsubsample.max,nrow=LOBS)
yb.h.pos<-vector()
# cycle over the available stations
b.inc<-0
for (b in yo.h.pos) {
  if (b.inc>0) {
    finoa<-as.integer(Lsubsample.max*1/4)
    if (VecY[b]<7750000 & VecZ[b]<800) finoa<-as.integer(Lsubsample.max*2/4)
    if (VecY[b]<7100000 & VecZ[b]<800) finoa<-as.integer(Lsubsample.max*3/4)
    if (VecS[b]%in%VecS.set[yb.h.pos[1:b.inc],2:finoa]) {
#              print(VecS[b])
#              print(VecS.set[yo.h.pos[1:(which(b==yo.h.pos)-1)],2:Lsubsample])
      next
    }
####        if (VecY[b]<7750000) {
####          if (VecY[b]<7100000) {
####            if (b.inc>0) {
####              aux.sort<-sort(Disth[b,yo.h.pos])
####              if (VecS[b]%in%VecS.set[yb.h.pos[1:b.inc],2:as.integer(Lsubsample.max/2)] &
####                  min(Disth[b,yb.h.pos[1:b.inc]])<(2*Dh.b) & aux.sort[2]<Dh.b) {
#####              print(VecS[b])
#####              print(VecS.set[yo.h.pos[1:(which(b==yo.h.pos)-1)],2:Lsubsample])
####                next
####              }
####            }
####        # VecS>=70000 & VecS<91000
####          } else {
####            if (b.inc>0) {
####              aux.sort<-sort(Disth[b,yo.h.pos])
####              if (VecS[b]%in%VecS.set[yb.h.pos[1:b.inc],2:as.integer(Lsubsample.max/4)] &
####                  min(Disth[b,yb.h.pos[1:b.inc]])<(2*Dh.b) & aux.sort[2]<Dh.b) {
#####              print(VecS[b])
#####              print(VecS.set[yo.h.pos[1:(which(b==yo.h.pos)-1)],2:Lsubsample])
####                next
####              }
 ####           }
####          }
####        # VecS>=91000
####        } else  {
####          if (b.inc>0) {
####            if (VecS[b]%in%VecS.set[yb.h.pos[1:b.inc],2:4]) {
#####              print(VecS[b])
#####              print(VecS.set[yo.h.pos[1:(which(b==yo.h.pos)-1)],2:Lsubsample])
####              next
####            }
 ####         }
####        }
  }
# select the closest (respect to b-th station) Lsubsample stations
  close2b.aux<-order(Disth[b,],decreasing=F)
  close2b.au1<-close2b.aux[which(close2b.aux%in%yo.h.pos)][1:Lsubsample]
  if (Disth[b,close2b.au1[2]]>(2*Dh.b)) next
  Lsubsample.vec[b]<-Lsubsample.max
  if (Disth[b,close2b.au1[Lsubsample.max]]>Lsubsample.DHmax) {
#        print(Disth[b,close2b.au1])
    Lsubsample.vec[b]<-max(which(Disth[b,close2b.au1]<=Lsubsample.DHmax))
#        print(Lsubsample.vec[b])
  }
  close2b<-close2b.au1[1:Lsubsample.vec[b]]
  b.inc<-b.inc+1
  yb.h.pos[b.inc]<-b
  yb.set0<-matrix(data=NA,nrow=LOBS)
  yb.set1<-matrix(data=NA,nrow=LOBS)
  yb.set2<-matrix(data=NA,nrow=LOBS)
  VecX.b<-vector(mode="numeric",length=Lsubsample.vec[b])
  VecY.b<-vector(mode="numeric",length=Lsubsample.vec[b])
  VecZ.b<-vector(mode="numeric",length=Lsubsample.vec[b])
  yo.b<-vector(mode="numeric",length=Lsubsample.vec[b])
  VecX.b<-VecX[close2b]
  VecY.b<-VecY[close2b]
  VecZ.b<-VecZ[close2b]
  VecS.set[b,]<-NA
  VecS.set.pos[b,]<-NA
  VecS.set[b,1:Lsubsample.vec[b]]<-VecS[close2b]
  VecS.set.pos[b,1:Lsubsample.vec[b]]<-close2b
  yo.b<-yo[close2b]
# NOTE: this is D.b! not D...
  ide.b<-matrix(data=0,ncol=Lsubsample.vec[b],nrow=Lsubsample.vec[b])
  ide.b[row(ide.b)==col(ide.b)]=1
  InvD.b<-solve(D.b[close2b,close2b],ide.b)
# Background 0  
  yb.param0<-XYZnoinv_step0(VecX.b,VecY.b,VecZ.b,yo.b)
  # 1.mean(z)[m];2.NA;3.mean(yo)[C];4.gamma[C/m];5.NA
  # 6. alpha[C/m];7.NA;8.Beta[C/m];9NA
  yb.set0[]<-yb.param0[3]+ yb.param0[6]*(VecX-yb.param0[10])+
                           yb.param0[8]*(VecY-yb.param0[11])+
                           yb.param0[4]* VecZ
  J0<-sum((yb.set0[close2b]-yo.b)**2.)
# Background 1  
  back1.flag<-T
  par.aux<-c(mean(VecZ.b),50,yb.param0[3],
             yb.param0[4],yb.param0[4],
             yb.param0[6],yb.param0[6],
             yb.param0[8],yb.param0[8])
  yb.param1<-XYZinv_step1(par.aux,VecX.b,VecY.b,VecZ.b,yo.b)
  yb.param1.mat[b,1:11]<-yb.param1
  if (!is.na(yb.param1[11])) {
    # 1.zinv[m];2.dz[m];3.Tinv[C];4.gammaAbove[C/m];5.gammaBelow[C]
    # 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
    zinv<-yb.param1[1]
    zabov<-zinv+yb.param1[2]
    zbelo<-zinv-yb.param1[2]
    bfabov<-yb.param1[3]+ yb.param1[6]*(VecX-yb.param1[10])+
                          yb.param1[8]*(VecY-yb.param1[11])+
                          yb.param1[4]*(VecZ-zinv)
    bfbelo<-yb.param1[3]+ yb.param1[7]*(VecX-yb.param1[10])+
                          yb.param1[9]*(VecY-yb.param1[11])+ 
                          yb.param1[5]*(VecZ-zinv)
    yb.set1[VecZ>zabov]<-bfabov[VecZ>zabov]
    yb.set1[VecZ<=zbelo]<-bfbelo[VecZ<=zbelo]
    aux<-(VecZ>zbelo)&(VecZ<=zabov)
    yb.set1[aux]<-(bfabov[aux]*(VecZ[aux]-zbelo) + bfbelo[aux]*(zabov-VecZ[aux]) ) / (zabov-zbelo)
    J1<-sum((yb.set1[close2b]-yo.b)**2.)
  } else {
    back1.flag<-F
    J1<-J0+1
  }
# Background 2  
# NOTE: yb.set* is defined for all the LOBSt stations
# 1.h0[m];2.h1-h0[m];3.T0[C];4.gamma[C/m];5.a[C]
# 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
  back2.flag<-T
  if (back1.flag) {
    par.aux<-c(zinv,
               (max(VecZ.b)-min(VecZ.b))/10,
               yb.param1[3]+yb.param1[4]*(-zinv),
               yb.param1[4],
               0,
               yb.param1[6],yb.param1[7],
               yb.param1[8],yb.param1[9])
  } else {
    par.aux<-c(mean(VecZ.b),
               (max(VecZ.b)-min(VecZ.b))/10,
               yb.param0[3],
               yb.param0[4],
               0,
               yb.param0[6],yb.param0[6],
               yb.param0[8],yb.param0[8])
  }
  yb.param2<-XYZinv_step2(par.aux,VecX.b,VecY.b,VecZ.b,yo.b)
  yb.param2.mat[b,1:11]<-yb.param2
  if (!is.na(yb.param2[11])) {
    yb.param2[2]<-abs(yb.param2[2])
    ab<-which( VecZ>=(yb.param2[1]+yb.param2[2]) )
    bl<-which( VecZ<=yb.param2[1] )
    bw<-which( (VecZ>yb.param2[1]) & (VecZ<(yb.param2[1]+yb.param2[2])) )
    yb.set2[ab]<-yb.param2[3] + yb.param2[4]*VecZ[ab] + 
                   yb.param2[6]*(VecX[ab]-yb.param2[10]) + 
                   yb.param2[8]*(VecY[ab]-yb.param2[11])
    if (length(bl)>0)
      yb.set2[bl]<-yb.param2[3] + yb.param2[4]*VecZ[bl] - 
                     yb.param2[5] +
                     yb.param2[7]*(VecX[bl]-yb.param2[10]) + 
                     yb.param2[9]*(VecY[bl]-yb.param2[11])
    if (length(bw)>0)
      yb.set2[bw]<-yb.param2[3] + yb.param2[4]*VecZ[bw] - 
                     yb.param2[5]/2.*(1+cos(pi*(VecZ[bw]-yb.param2[1])/yb.param2[2])) + 
                    ( (yb.param2[6]*(VecX[bw]-yb.param2[10])+yb.param2[8]*(VecY[bw]-yb.param2[11]))*(VecZ[bw]-yb.param2[1]) + 
                      (yb.param2[7]*(VecX[bw]-yb.param2[10])+yb.param2[9]*(VecY[bw]-yb.param2[11]))*(yb.param2[1]+yb.param2[2]-VecZ[bw]) ) / yb.param2[2]
    J2<-sum((yb.set2[close2b]-yo.b)**2.)
  } else {
    back2.flag<-F
    J2<-J0+1
  }
# Best background
  aux<-order(c(J0,J1,J2))
  best<-aux[1]-1
  if (best==0) {
    yb.set[,b]<-yb.set0[]
    yb.param[b,]<-c(yb.param0,0)
  }
  if (best==1) {
    yb.set[,b]<-yb.set1[]
    yb.param[b,]<-c(yb.param1,1)
  }
  if (best==2) {
    yb.set[,b]<-yb.set2[]
    yb.param[b,]<-c(yb.param2,2)
  }
# DEBUG start
  zero.str<-"0."
  uno.str<-"1."
  due.str<-"2."
  if (best==0) zero.str<-"+0."
  if (best==1) uno.str<-"+1."
  if (best==2) due.str<-"+2."
  print(paste("@@",b.inc,". id pos Zmn/x DisthMAX / J0 J1 J2 / #stn:",VecS[b],b,
                           round(min(VecZ.b),0),round(max(VecZ.b),0),
                           round(Disth[b,close2b[Lsubsample.vec[b]]],0),"/",
                           round(J0,0),round(J1,0),round(J2,0),"/",Lsubsample.vec[b]))
  print(paste(zero.str,"alpha beta gamma:",round(yb.param0[6],6),round(yb.param0[8],6),round(yb.param0[4],4)))
  print(paste(uno.str,"zinv dz tinv aA aB bA bB ga gb:",
        round(yb.param1[1],1), round(yb.param1[2],0), round(yb.param1[3],1),
        round(yb.param1[6],6), round(yb.param1[7],6), round(yb.param1[8],6),
        round(yb.param1[9],6), round(yb.param1[4],4), round(yb.param1[5],4)))
  print(paste(due.str,"h0 h1-h0 t0 a aA aB bA bB gamma:",
        round(yb.param2[1],1), round(yb.param2[2],0), round(yb.param2[3],1),
        round(yb.param2[5],6), round(yb.param2[6],6), round(yb.param2[7],6),
        round(yb.param2[8],6), round(yb.param2[9],6), round(yb.param2[4],4)))
# DEBUG stop
# G1 is LOBStOK x Lsubsample matrix to interpolate the Lsubsample values
# on the whole LOBStOK station dataset
  G1.b<-S.b[,close2b]
  K.b<-G1.b%*%InvD.b
  rm(G1.b)
  W.b<-S.b[close2b,close2b]%*%InvD.b
  # this is idi (sort of)
  ybweights.set[,b]<-rowSums(K.b)
  if (any(is.na(ybweights.set[,b]))) {
    print("Error: NA value presents in ybweights.set vector. Unexpected.")
    print(paste(round(yb.set[,b],2),round(yo,2),round(ybweights.set[,b],3),"\n"))
    q(status=1)
  }
  rm(K.b,W.b)
} # end cycle LOBSt
LBAKh<-b.inc
print(paste("# stations used in background elaborations=",LBAKh))
# At this point I've this 5 outputs
# 1. yb.set<-matrix(data=NA,ncol=btimes,nrow=LOBSt)
# 2. ybweights.set<-matrix(data=NA,ncol=btimes,nrow=LOBSt)
# 3. yb.param<-matrix(data=NA,ncol=12,nrow=LOBSt)
# 4. VecS.set<-matrix(data=NA,ncol=Lsubsample,nrow=LOBSt)
# 5. VecS.set.pos[b,]<-close2b
#
#====================================================================
### normalization of the yb weights such that their sum equals to one
##    ybweights.norm<-ybweights.set[,yo.h.pos] / 
##                    rowSums(ybweights.set[,yo.h.pos])
##    print(ncol(ybweights.norm))
##    print(nrow(ybweights.norm))
### background on station points
##    for (m in 1:LOBS) {
##      yb[m]<-ybweights.norm[m,] %*% yb.set[m,yo.h.pos]
##    }
##    print(paste("id yo yb",VecS,
##                round(yo,1),
##                round(yb,2),"\n"))
##    q()
#====================================================================
# Station Analysis cycle (with SCT!)
yoBad.id<-NA
yoBad.pos<-NA
isct<-0
ydqc.flag[]<-rep(0,LOBS)
ydqc.flag[-yo.h.pos]<-(-1)
yo.OKh.pos<-which(!is.na(yo) & ydqc.flag!=1)
LOBStOK<-length(yo.OKh.pos)
while (TRUE) {
  print(paste(">> Total number of observations [not NA & not ERR(so far)] =",
              LOBStOK))
# Station points - Background
  if (!is.na(yoBad.id)) {
    isct<-isct+1
    for (b in yb.h.pos) {
      if (any(VecS.set[b,1:Lsubsample.vec[b]]==yoBad.id)) {
        yb.set0<-matrix(data=NA,nrow=LOBS)
        yb.set1<-matrix(data=NA,nrow=LOBS)
        yb.set2<-matrix(data=NA,nrow=LOBS)
        Lsubsample.vec[b]<-Lsubsample.vec[b]-1
        VecX.b<-vector(mode="numeric",length=Lsubsample.vec[b])
        VecY.b<-vector(mode="numeric",length=Lsubsample.vec[b])
        VecZ.b<-vector(mode="numeric",length=Lsubsample.vec[b])
        yo.b<-vector(mode="numeric",length=Lsubsample.vec[b])
# select the closest (respect to indxb stations) Lsubsample stations
        close2b.aux<-order(Disth[b,],decreasing=F)
        close2b<-close2b.aux[which(close2b.aux%in%yo.OKh.pos)][1:Lsubsample.vec[b]]
        VecX.b<-VecX[close2b]
        VecY.b<-VecY[close2b]
        VecZ.b<-VecZ[close2b]
        VecS.set[b,]<-NA
        VecS.set.pos[b,]<-NA
        VecS.set[b,1:Lsubsample.vec[b]]<-VecS[close2b]
        VecS.set.pos[b,1:Lsubsample.vec[b]]<-close2b
        yo.b<-yo[close2b]
# NOTE: this is D.b! not D...
        ide.b<-matrix(data=0,ncol=Lsubsample.vec[b],nrow=Lsubsample.vec[b])
        ide.b[row(ide.b)==col(ide.b)]<-1
        InvD.b<-solve(D.b[close2b,close2b],ide.b)
# Background 0  
        yb.param0<-XYZnoinv_step0(VecX.b,VecY.b,VecZ.b,yo.b)
        # 1.mean(z)[m];2.NA;3.mean(yo)[C];4.gamma[C/m];5.NA
        # 6. alpha[C/m];7.NA;8.Beta[C/m];9NA
        yb.set0[]<-yb.param0[3]+ yb.param0[6]*(VecX-yb.param0[10])+
                                 yb.param0[8]*(VecY-yb.param0[11])+
                                 yb.param0[4]*VecZ
        J0<-sum((yb.set0[close2b]-yo.b)**2.)
# Background 1  
        back1.flag<-T
        if (!is.na(yb.param1.mat[b,11])) {
          par.aux<-yb.param1.mat[b,1:9]
        } else {
          par.aux<-c(mean(VecZ.b),50,yb.param0[3],
                     yb.param0[4],yb.param0[4],
                     yb.param0[6],yb.param0[6],
                     yb.param0[8],yb.param0[8])
        }
        yb.param1<-XYZinv_step1(par.aux,VecX.b,VecY.b,VecZ.b,yo.b)
        yb.param1.mat[b,1:11]<-yb.param1
        if (!is.na(yb.param1[11])) {
          # 1.zinv[m];2.dz[m];3.Tinv[C];4.gammaAbove[C/m];5.gammaBelow[C]
          # 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
          zinv<-yb.param1[1]
          zabov<-zinv+yb.param1[2]
          zbelo<-zinv-yb.param1[2]
          bfabov<-yb.param1[3]+ yb.param1[6]*(VecX-yb.param1[10])+
                                yb.param1[8]*(VecY-yb.param1[11])+
                                yb.param1[4]*(VecZ-zinv)
          bfbelo<-yb.param1[3]+ yb.param1[7]*(VecX-yb.param1[10])+
                                yb.param1[9]*(VecY-yb.param1[11])+ 
                                yb.param1[5]*(VecZ-zinv)
          yb.set1[VecZ>zabov]<-bfabov[VecZ>zabov]
          yb.set1[VecZ<=zbelo]<-bfbelo[VecZ<=zbelo]
          aux<-(VecZ>zbelo)&(VecZ<=zabov)
          yb.set1[aux]<-(bfabov[aux]*(VecZ[aux]-zbelo) + bfbelo[aux]*(zabov-VecZ[aux]) ) / (zabov-zbelo)
          J1<-sum((yb.set1[close2b]-yo.b)**2.)
        } else {
          back1.flag<-F
          J1<-J0+1
        }
# Background 2  
# NOTE: yb.set* is defined for all the LOBSt stations
# 1.h0[m];2.h1-h0[m];3.T0[C];4.gamma[C/m];5.a[C]
# 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
        back2.flag<-T
        if (!is.na(yb.param2.mat[b,11])) {
          par.aux<-yb.param2.mat[b,1:9]
        } else {
          if (back1.flag) {
            par.aux<-c(zinv,
                   (max(VecZ.b)-min(VecZ.b))/10,
                   yb.param1[3]+yb.param1[4]*(-zinv),
                   yb.param1[4],
                   0,
                   yb.param1[6],yb.param1[7],
                   yb.param1[8],yb.param1[9])
          } else {
            par.aux<-c(mean(VecZ.b),
                       (max(VecZ.b)-min(VecZ.b))/10,
                       yb.param0[3],
                       yb.param0[4],
                       0,
                       yb.param0[6],yb.param0[6],
                       yb.param0[8],yb.param0[8])
          }
        }
        yb.param2<-XYZinv_step2(par.aux,VecX.b,VecY.b,VecZ.b,yo.b)
        yb.param2.mat[b,1:11]<-yb.param2
        if (!is.na(yb.param2[11])) {
          yb.param2[2]<-abs(yb.param2[2])
          ab<-which( VecZ>=(yb.param2[1]+yb.param2[2]) )
          bl<-which( VecZ<=yb.param2[1] )
          bw<-which( (VecZ>yb.param2[1]) & (VecZ<(yb.param2[1]+yb.param2[2])) )
          yb.set2[ab]<-yb.param2[3] + yb.param2[4]*VecZ[ab] + 
                         yb.param2[6]*(VecX[ab]-yb.param2[10]) + 
                         yb.param2[8]*(VecY[ab]-yb.param2[11])
          if (length(bl)>0)
            yb.set2[bl]<-yb.param2[3] + yb.param2[4]*VecZ[bl] - 
                           yb.param2[5] +
                           yb.param2[7]*(VecX[bl]-yb.param2[10]) + 
                           yb.param2[9]*(VecY[bl]-yb.param2[11])
          if (length(bw)>0)
            yb.set2[bw]<-yb.param2[3] + yb.param2[4]*VecZ[bw] - 
                           yb.param2[5]/2.*(1+cos(pi*(VecZ[bw]-yb.param2[1])/yb.param2[2])) + 
                          ( (yb.param2[6]*(VecX[bw]-yb.param2[10])+yb.param2[8]*(VecY[bw]-yb.param2[11]))*(VecZ[bw]-yb.param2[1]) + 
                            (yb.param2[7]*(VecX[bw]-yb.param2[10])+yb.param2[9]*(VecY[bw]-yb.param2[11]))*(yb.param2[1]+yb.param2[2]-VecZ[bw]) ) / yb.param2[2]
          J2<-sum((yb.set2[close2b]-yo.b)**2.)
        } else {
          back2.flag<-F
          J2<-J0+1
        }
# Best background
        aux<-order(c(J0,J1,J2))
        best<-aux[1]-1
        if (best==0) {
          yb.set[,b]<-yb.set0[]
          yb.param[b,]<-c(yb.param0,0)
        }
        if (best==1) {
          yb.set[,b]<-yb.set1[]
          yb.param[b,]<-c(yb.param1,1)
        }
        if (best==2) {
          yb.set[,b]<-yb.set2[]
          yb.param[b,]<-c(yb.param2,2)
        }
# DEBUG start
        zero.str<-"0."
        uno.str<-"1."
        due.str<-"2."
        if (best==0) zero.str<-"*0."
        if (best==1) uno.str<-"*1."
        if (best==2) due.str<-"*2."
        print(paste("@@",isct,". id pos Zmn/x DisthMAX / J0 J1 J2 / #stn:",VecS[b],b,
                                 round(min(VecZ.b),0),round(max(VecZ.b),0),
                                 round(Disth[b,close2b[Lsubsample.vec[b]]],0),"/",
                                 round(J0,0),round(J1,0),round(J2,0),"/",Lsubsample.vec[b]))
        print(paste(zero.str,"alpha beta gamma:",round(yb.param0[6],6),round(yb.param0[8],6),round(yb.param0[4],4)))
        print(paste(uno.str,"zinv dz tinv aA aB bA bB ga gb:",
              round(yb.param1[1],1), round(yb.param1[2],0), round(yb.param1[3],1),
              round(yb.param1[6],6), round(yb.param1[7],6), round(yb.param1[8],6),
              round(yb.param1[9],6), round(yb.param1[4],4), round(yb.param1[5],4)))
        print(paste(due.str,"h0 h1-h0 t0 a aA aB bA bB gamma:",
              round(yb.param2[1],1), round(yb.param2[2],0), round(yb.param2[3],1),
              round(yb.param2[5],6), round(yb.param2[6],6), round(yb.param2[7],6),
              round(yb.param2[8],6), round(yb.param2[9],6), round(yb.param2[4],4)))
# DEBUG stop
# G1 is LOBStOK x Lsubsample matrix to interpolate the Lsubsample values
# on the whole LOBStOK station dataset
        G1.b<-S.b[,close2b]
        K.b<-G1.b%*%InvD.b
        rm(G1.b)
        W.b<-S.b[close2b,close2b]%*%InvD.b
        # this is idi (sort of)
        ybweights.set[,b]<-rowSums(K.b)
        if (any(is.na(ybweights.set[,b]))) {
          print("Error: NA value presents in ybweights.set vector. Unexpected.")
          print(paste(round(yb.set[,b],2),round(yo,2),round(ybweights.set[,b],3),"\n"))
          q(status=1)
        }
        rm(K.b,W.b)
      }
    } # end cycle btimes (no optimization: over the station number)
  }
# normalization of the yb weights such that their sum equals to one
  ybweights.norm<-ybweights.set[,yb.h.pos] / 
                  rowSums(ybweights.set[,yb.h.pos])
# background on station points
  for (m in 1:LOBS) {
    yb[m]<-ybweights.norm[m,] %*% yb.set[m,yb.h.pos]
  }
# deallocate memory
#      rm(D.b,S.b,yb.set,ybweights.set,ybweights.norm)
#      rm(VecY.b,VecX.b,VecZ.b,yo.b,close2b,y1.b)
# Station (CV)Analysis/(CV)IDI
  ide<-matrix(data=0,ncol=LOBStOK,nrow=LOBStOK)
  ide[row(ide)==col(ide)]=1
  InvD<-solve(D[yo.OKh.pos,yo.OKh.pos],ide)
  W<-S[yo.OKh.pos,yo.OKh.pos] %*% InvD
  G1<-S[,yo.OKh.pos]
  K<-G1 %*% InvD
  rm(G1)
  ya<-yb + K %*% (yo[yo.OKh.pos]-yb[yo.OKh.pos])
  yav<-ya
  yav[yo.OKh.pos]<-yo[yo.OKh.pos] + 1./(1.-diag(W)) * (ya[yo.OKh.pos]-yo[yo.OKh.pos])
  yidi<-rowSums(K)
  yidiv<-yidi
  yidiv[yo.OKh.pos]<-rep(1,LOBStOK) + 1./(1.-diag(W)) * (yidi[yo.OKh.pos]-rep(1,LOBStOK))
# DQC CHECK - Spatial Continuity Check
  ydqc[]<--9999.
  ydqc[yo.OKh.pos]<-(yo[yo.OKh.pos]-yav[yo.OKh.pos]) *
                               (yo[yo.OKh.pos]-ya[yo.OKh.pos]) / sig2o
# DQC test: T2 is the SCT threshold, if any(ydqc>T2) is true then reject the 
# station with highest ydqc value. This mean that the station is not used to
# compute the analysis
  if (max(ydqc)<=T2) {
    break
  }
  aux<-which.max(ydqc)
  yoBad.id<-VecS[aux]
  yoBad.pos<-aux
  ydqc.flag[yoBad.pos]<-1
  yo.OKh.pos<-which(!is.na(yo) & ydqc.flag!=1)
  LOBStOK<-length(yo.OKh.pos)
# DQC flag "1" means erroneous observation
##      OBS$DQC[indx[aux]]<-1
  print(paste("SCT found GE-> id yo yb ya yav ydqc",VecS[aux],
              round(yo[aux],1),
              round(yb[aux],2),
              round(ya[aux],2),
              round(yav[aux],2),
              round(ydqc[aux],2),"\n"))
# write output file
} # End of Station Analysis cycle (with SCT!)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Station Points - Write output on file 
cat(paste(yyyy,mm,dd,h,
          round(VecS,0),round(VecX,0),round(VecY,0),round(VecZ,0),
          round(yo,1),round(yb,2),round(ya,2),round(yav,2),
          round(yidi,3),round(yidiv,3),round(ydqc.flag,2),
          "\n",sep=";"),
          file=out.file.stn,append=T)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# if something strange take place in the Analysis/SCT cycle the pass to the next iteration
if (LOBStOK<=0) next
# DQC flag "0" means good observations
print(paste(">>>>Total number of observations [not NA & good] =",LOBStOK))
#------------------------------------------------------------------------------
# Grid - Background
xb<-vector(mode="numeric",length=Lgrid)
xb.set<-matrix(data=0,ncol=LBAKh,nrow=Lgrid)
xbweights.set<-matrix(data=NA,ncol=LBAKh,nrow=Lgrid)
##    for (b in 1:btimes) {
b.aux<-0
print("++ Grid - Background elaborations\n")
print("#/tot stnid #grid.points/#grid.points.tot")
for (b in yb.h.pos) {
  b.aux<-b.aux+1
#      xindx<-which( (abs(xgrid-VecX[b])<(2*(Dh.b*1000))) & (abs(ygrid-VecY[b])<(2*(Dh.b*1000))) )
  xindx<-which( (xgrid-min(VecX[VecS.set.pos[b,1:Lsubsample.vec[b]]]))>(-3*Dh.b*1000) & (xgrid-max(VecX[VecS.set.pos[b,1:Lsubsample.vec[b]]]))<(3*Dh.b*1000) &
                (ygrid-min(VecY[VecS.set.pos[b,1:Lsubsample.vec[b]]]))>(-3*Dh.b*1000) & (ygrid-max(VecY[VecS.set.pos[b,1:Lsubsample.vec[b]]]))<(3*Dh.b*1000) )
  Lgrid.b<-length(xindx)
  if (Lgrid.b==0) {
    xb.set[,b.aux]<-0
    xbweights.set[,b.aux]<-0.
    next
  }
  ide.b<-matrix(data=0,ncol=Lsubsample.vec[b],nrow=Lsubsample.vec[b])
  ide.b[row(ide.b)==col(ide.b)]<-1
  InvD.b<-solve(D.b[VecS.set.pos[b,1:Lsubsample.vec[b]],VecS.set.pos[b,1:Lsubsample.vec[b]]],ide.b)
# G matrix
  Disth.b<-matrix(ncol=Lsubsample.vec[b],nrow=Lgrid.b,data=0.)
#      Distz.b<-matrix(ncol=Lsubsample.vec[b],nrow=Lgrid.b,data=0.)
  G.b<-matrix(ncol=Lsubsample.vec[b],nrow=Lgrid.b,data=0.)
  Disth.b<-(outer(ygrid[xindx],VecY[VecS.set.pos[b,1:Lsubsample.vec[b]]],FUN="-")**2.+
            outer(xgrid[xindx],VecX[VecS.set.pos[b,1:Lsubsample.vec[b]]],FUN="-")**2.)**0.5/1000.
#      Distz.b<-abs(outer(zgrid[xindx],VecZ[VecS.set.pos[b,1:Lsubsample.vec[b]]],FUN="-"))
#      G.b<-exp(-0.5*(Disth.b/Dh.b)**2.-0.5*(Distz.b/Dz.b)**2.)
  G.b<-exp(-0.5*(Disth.b/Dh.b)**2.)
#      rm(Disth.b,Distz.b)
  rm(Disth.b)
#  compute analysis/idi over grid 
  K.b<-G.b%*%InvD.b
  xbweights.set[,b.aux]<-0
  xbweights.set[xindx,b.aux]<-rowSums(K.b)
  xindx1<-which(xbweights.set[,b.aux]>0.05)
  Lgrid.b1<-length(xindx1)
  print(paste(b.aux,"/",LBAKh," ",VecS[b]," ",Lgrid.b,"-->",Lgrid.b1,"/",Lgrid,"\n",sep=""))
  rm(G.b,K.b)
  if (yb.param[b,12]==0) {
    # 1.mean(z)[m];2.NA;3.mean(yo)[C];4.gamma[C/m];5.NA
    # 6. alpha[C/m];7.NA;8.Beta[C/m];9NA
    xb.set[xindx1,b.aux]<-yb.param[b,3]+ yb.param[b,6]*(xgrid[xindx1]-yb.param[b,10])+
                                        yb.param[b,8]*(ygrid[xindx1]-yb.param[b,11])+
                                        yb.param[b,4]* zgrid[xindx1]
  }
  if (yb.param[b,12]==1) {
    # 1.zinv[m];2.dz[m];3.Tinv[C];4.gammaAbove[C/m];5.gammaBelow[C]
    # 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
    zinv<-yb.param[b,1]
    zabov<-zinv+yb.param[b,2]
    zbelo<-zinv-yb.param[b,2]
    bfabov<-yb.param[b,3]+ yb.param[b,6]*(xgrid[xindx1]-yb.param[b,10])+
                           yb.param[b,8]*(ygrid[xindx1]-yb.param[b,11])+
                           yb.param[b,4]*(zgrid[xindx1]-zinv)
    bfbelo<-yb.param[b,3]+ yb.param[b,7]*(xgrid[xindx1]-yb.param[b,10])+
                           yb.param[b,9]*(ygrid[xindx1]-yb.param[b,11])+ 
                           yb.param[b,5]*(zgrid[xindx1]-zinv)
    aux.ab<-which(zgrid[xindx1]>zabov)
    aux.bl<-which(zgrid[xindx1]<=zbelo)
    aux.bw<-which((zgrid[xindx1]>zbelo)&(zgrid[xindx1]<=zabov))
    xb.set[xindx1[aux.ab],b.aux]<-bfabov[aux.ab]
    xb.set[xindx1[aux.bl],b.aux]<-bfbelo[aux.bl]
    xb.set[xindx1[aux.bw],b.aux]<-(bfabov[aux.bw]*(zgrid[xindx1[aux.bw]]-zbelo) + bfbelo[aux.bw]*(zabov-zgrid[xindx1[aux.bw]]) ) / (zabov-zbelo)
    rm(aux.ab,aux.bl,aux.bw)
  }
  if (yb.param[b,12]==2) {
# 1.h0[m];2.h1-h0[m];3.T0[C];4.gamma[C/m];5.a[C]
# 6. alphaAbove[C/m];7.AlphaBelow[C/m];8.BetaAbove[C/m];9BetaBelow[C/m]
    h0<-yb.param[b,1]
    h1<-yb.param[b,1]+yb.param[b,2]
    aux.ab<-which(zgrid[xindx1]>=h1)
    aux.bl<-which(zgrid[xindx1]<=h0)
    aux.bw<-which((zgrid[xindx1]>h0) & (zgrid[xindx1]<h1))
    xb.set[xindx1[aux.ab],b.aux]<-yb.param[b,3] + 
                             yb.param[b,4]* zgrid[xindx1[aux.ab]] + 
                             yb.param[b,6]*(xgrid[xindx1[aux.ab]]-yb.param[b,10]) + 
                             yb.param[b,8]*(ygrid[xindx1[aux.ab]]-yb.param[b,11])
    xb.set[xindx1[aux.bl],b.aux]<-yb.param[b,3] + 
                             yb.param[b,4]*zgrid[xindx1[aux.bl]] - 
                             yb.param[b,5] +
                             yb.param[b,7]*(xgrid[xindx1[aux.bl]]-yb.param[b,10]) + 
                             yb.param[b,9]*(ygrid[xindx1[aux.bl]]-yb.param[b,11])
    xb.set[xindx1[aux.bw],b.aux]<-yb.param[b,3] +
                             yb.param[b,4]*zgrid[xindx1[aux.bw]] - 
                             yb.param[b,5]/2.*(1+cos(pi*(zgrid[xindx1[aux.bw]]-h0)/yb.param[b,2])) + 
     ( (yb.param[b,6]*(xgrid[xindx1[aux.bw]]-yb.param[b,10])+yb.param[b,8]*(ygrid[xindx1[aux.bw]]-yb.param[b,11]))*(zgrid[xindx1[aux.bw]]-h0) + 
       (yb.param[b,7]*(xgrid[xindx1[aux.bw]]-yb.param[b,10])+yb.param[b,9]*(ygrid[xindx1[aux.bw]]-yb.param[b,11]))*(h1-zgrid[xindx1[aux.bw]]) ) / yb.param[b,2]
    rm(aux.ab,aux.bl,aux.bw)
  }
# G matrix
#      Disth.b<-matrix(ncol=Lsubsample.vec[b],nrow=Lgrid.b,data=0.)
#      Distz.b<-matrix(ncol=Lsubsample.vec[b],nrow=Lgrid.b,data=0.)
#      G.b<-matrix(ncol=Lsubsample.vec[b],nrow=Lgrid.b,data=0.)
#      Disth.b<-(outer(ygrid[xindx],VecY[VecS.set.pos[b,1:Lsubsample.vec[b]]],FUN="-")**2.+
#                outer(xgrid[xindx],VecX[VecS.set.pos[b,1:Lsubsample.vec[b]]],FUN="-")**2.)**0.5/1000.
#      Distz.b<-abs(outer(zgrid[xindx],VecZ[VecS.set.pos[b,1:Lsubsample.vec[b]]],FUN="-"))
##      G.b<-exp(-0.5*(Disth.b/Dh.b)**2.-0.5*(Distz.b/Dz.b)**2.)
#      G.b<-exp(-0.5*(Disth.b/Dh.b)**2.)
#      rm(Disth.b,Distz.b)
##  compute analysis/idi over grid 
#      K.b<-G.b%*%InvD.b
#      xbweights.set[,b.aux]<-0
#      xbweights.set[xindx,b.aux]<-rowSums(K.b)
#      print(paste(b.aux,"/",LBAKh," ",VecS[b]," ",Lgrid.b,"/",Lgrid,"\n",sep=""))
#      rm(G.b,K.b)
}
xbweights.norm<-xbweights.set/rowSums(xbweights.set)
for (i in 1:Lgrid) {
  xb[i]<-xbweights.norm[i,] %*% xb.set[i,]
}
rm(xb.set,xbweights.norm,xbweights.set,xindx,xindx1)
#------------------------------------------------------------------------------
# Gridded Analysis/IDI    
xa<-vector(mode="numeric",length=Lgrid)
xidi<-vector(mode="numeric",length=Lgrid)
ndim<-10000
i<-0
print("++ Grid - Analysis and IDI elaborations\n")
print("# gridpoint.start gridpoint.end gridpoint.total\n")
#to.write <- file("aux.bin", "wb")
#to.read<-file("aux.bin", "rb")
while ((i*ndim)<Lgrid) {
  start<-i*ndim+1
  end<-(i+1)*ndim
  if (end>Lgrid) {
    end<-Lgrid
  }
  ndimaux<-end-start+1
  print(paste(round(i),round(start,0),round(end,0),round(Lgrid,0)))
  aux<-matrix(ncol=LOBStOK,nrow=ndimaux,data=0.)
  auxz<-matrix(ncol=LOBStOK,nrow=ndimaux,data=0.)
  G<-matrix(ncol=LOBStOK,nrow=ndimaux,data=0.)
  aux<-(outer(ygrid[start:end],VecY[yo.OKh.pos],FUN="-")**2. +
        outer(xgrid[start:end],VecX[yo.OKh.pos],FUN="-")**2.)**0.5/1000.
#  for (j in 1:Linfo) writeBin(as.numeric(aux[,j]),to.write,size=4)
#  aux<-matrix(readBin(to.read,numeric(),size=4,ndimaux*Linfo),ndimaux,Linfo)
# vertical distance
  auxz<-abs(outer(zgrid[start:end],VecZ[yo.OKh.pos],FUN="-"))
  G<-exp(-0.5*(aux/Dh)**2.-0.5*(auxz/Dz)**2.)
  rm(aux,auxz)
#  Gt<-exp(-0.5*(aux/Dh)**2.)
  K<-G%*%InvD
  rm(G)
  xa[start:end]<-xb[start:end]+K%*%(yo[yo.OKh.pos]-yb[yo.OKh.pos])
  xidi[start:end]<-rowSums(K)
#  Disth[start:end,1:Linfo] <-aux[1:ndimaux,1:Linfo]
  rm(K)
#  print("finito conti")
  i<-i+1
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# output - ANALYSIS
xx[]<-NA
# 
xx[ccgrid]<-round(xa,1)
nogrid.ncout(grid=t(as.matrix(xx)),
             x=x.G,y=y.G,grid.type=grid.type,
             file.name=out.file.grd.ana,
             var.name=var.name.xa,
             var.longname=var.longname.xa,
             var.unit=var.unit.xa,
             var.mv=var.mv.xa,
             var.version=var.version.xa,
             times=c(paste(yyyymmddhh,"00",sep="")),times.unit=times.unit.xa,
             times.ref=times.ref.xa,
             prod.date=prod.date,
             reference=reference.xa,
             proj4.string="+proj=utm +zone=33 +ellps=WGS84",
             source.string=source.nc)
# output - background
xx[]<-NA
xx[ccgrid]<-round(xb,1)
nogrid.ncout(grid=t(as.matrix(xx)),
             x=x.G,y=y.G,grid.type=grid.type,
             file.name=out.file.grd.bck,
             var.name=var.name.xa,
             var.longname=var.longname.xa,
             var.unit=var.unit.xa,
             var.mv=var.mv.xa,
             var.version=var.version.xa,
             times=c(paste(yyyymmddhh,"00",sep="")),times.unit=times.unit.xa,
             times.ref=times.ref.xa,
             prod.date=prod.date,
             reference=reference.xa,
             proj4.string="+proj=utm +zone=33 +ellps=WGS84",
             source.string=source.nc)
# output - IDI
xx[]<-NA
xx[ccgrid]<-round(xidi,3)
rnc <- writeRaster(xx, filename=out.file.grd.idi,format="CDF",overwrite=TRUE)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Exit - Success
q(status=0)

