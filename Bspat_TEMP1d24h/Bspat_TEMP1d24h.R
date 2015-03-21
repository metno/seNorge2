# --~- Bspat_TEMP1d24h.R  -~--
# Bayesian spatial interpolation of daily mean air temperature (06-06 UTC)
# Spatial Consistency Test (SCT) is included.
# Observation: TAMRR from KDVH
# Background: Spatial Interpolation of TA
#
# Outputs: look for "@@@@@@@@@" in the code and you'll get the ouput formats
#
# History:
# 02.12.2014 - Cristian Lussana. Original code.
# 25.02.2015 - Cristian Lussana. version 1.0.
#  change log: 
#  - allow for the use of observations outside Norway
#  - revisied queries to KDVH
#  - geographical information on seNorge2_dem_UTM33.nc
#  - definition of TEMP1d24h
#  - definition of a new directory tree
# ==~==========================================================================
rm(list=ls())
# Libraries
library(raster)
library(rgdal)
library(colorspace)
library(ncdf)
#-------------------------------------------------------------------
# FUNCTIONS 
# manage fatal error
error_exit<-function(str=NULL) {
  print("Fatal Error:")
  if (!is.null(str)) print(str)
  quit(status=1)
}
#-------------------------------------------------------------------
# CRS strings
proj4.wgs84<-"+proj=longlat +datum=WGS84"
proj4.ETRS_LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
#
max.Km.stnINdomain<-50
#-----------------------------------------------------------------------------
# [] Setup OI parameters
sig2o<-3.0
eps2<-0.5
Dh<-60
Dz<-600
T2<-20
# Background - parameters
Dh.b<-70.
eps2.b<-0.5
Lsubsample<-50
Lsubsample.max<-50
Lsubsample.DHmax<-200
Lsubsample.vec<-vector()
#
print("ANALYSIS and DQC parameters")
print(paste("EPS2 Dh[Km] Dz[m] > ", round(eps2,3),
                                    round(Dh,3),
                                    round(Dz,3),
                                    sep=" "))
print(paste("VARobs[C^2] T^2 = ", round(sig2o,3), round(T2,3),sep=" "))
# MAIN ========================================================================
# Read command line arguments
arguments <- commandArgs()
arguments
date.string<-arguments[3]
config_file<-arguments[4]
config_par<-arguments[5]
if (length(arguments)!=5) 
  ext<-error_exit(paste("Error in command line arguments: \n",
  " R --vanilla yyyy.mm.dd configFILE configPAR \n",
  sep=""))
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
if (!file.exists(paste(main.path.geoinfo,"/seNorge2_dem_UTM33.nc",sep=""))) 
  ext<-error_exit(paste("File not found:",main.path.geoinfo,"/seNorge2_dem_UTM33.nc"))
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
source(paste(path2lib.com,"/SpInt_PseudoBackground.R",sep=""))
# test mode
print(testmode)
if (testmode) {
  print("TESTMODE TESTMODE TESTMODE")
  if (file.exists(paste(main.path,"/Bspat_TEMP1d24h/testbed",sep=""))) {
    testbed<-paste(main.path,"/Bspat_TEMP1d24h/testbed",sep="")
    station.info<-paste(testbed,"/station_data.csv",sep="")
    observed.data<-paste(testbed,"/observed_data.csv",sep="")
  } else {
    ext<-error_exit(paste("testbed not found"))
  }
}
# netcdf fixed parameters
grid.type <- "utm33"
prod.date <- substr(Sys.time(),1,10)
xa.source.nc<-"daily mean air temperature from station data (06-06 UTC)"
xa.var.version <- "1.0"
xa.pname<-"TEMP1d"
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
xidi.source.nc<-"IDI for daily mean air temperature from station data (06-06 UTC)"
xidi.var.version <- "1.0"
xidi.pname<-"IDI"
for (p in 1:length(ncout.spec.list)) {
  if (ncout.spec.list[[p]]$pname==xidi.pname) break
}
if (p==length(ncout.spec.list) & (ncout.spec.list[[p]]$pname!=xidi.pname) ) {
  xidi.flag.write<-F
} else {
  xidi.flag.write<-T
  xidi.var.name<-ncout.spec.list[[p]]$opts$var.name
  xidi.var.longname<-ncout.spec.list[[p]]$opts$var.longname
  xidi.var.unit<-ncout.spec.list[[p]]$opts$var.unit
  xidi.var.mv<-as.numeric(ncout.spec.list[[p]]$opts$var.mv)
  xidi.times.unit <-ncout.spec.list[[p]]$opts$t.unit
  xidi.times.ref <-ncout.spec.list[[p]]$opts$t.ref
  xidi.reference <- ncout.spec.list[[p]]$opts$reference
}
# set Time-related variables
yyyy<-substr(date.string,1,4)
mm<-substr(date.string,6,7)
dd<-substr(date.string,9,10)
date.dot<-paste(dd,".",mm,".",yyyy,sep="")
yyyymm<-paste(yyyy,mm,sep="")
yyyymmdd<-paste(yyyymm,dd,sep="")
end.string<-paste(substr(date.string,1,10),".06",sep="")
end<-strptime(end.string,"%Y.%m.%d.%H","UTC")
timeseq<-as.POSIXlt(seq(as.POSIXlt(end),length=24,by="-1 hour"),"UTC")
yyyy.v<-timeseq$year+1900
mm.v<-timeseq$mon+1
dd.v<-timeseq$mday
hh.v<-timeseq$hour
yyyymm.v<-paste(yyyy.v,formatC(mm.v,width=2,flag="0"),sep="")
yyyymmddhh.v<-paste(yyyy.v,
                  formatC(mm.v,width=2,flag="0"),
                  formatC(dd.v,width=2,flag="0"),
                  formatC(hh.v,width=2,flag="0"),sep="")
n.timeseq<-length(yyyymmddhh.v)
# output directories
dir.create(file.path(main.path.output,"seNorge2"), showWarnings = FALSE)
dir.create(file.path(main.path.output,"seNorge2_addInfo"), showWarnings = FALSE)
path2output.main<-paste(main.path.output,"/seNorge2/TEMP1d24h",sep="")
path2output.main.stn<-paste(path2output.main,"/station_dataset",sep="")
path2output.main.grd<-paste(path2output.main,"/gridded_dataset",sep="")
path2output.add<-paste(main.path.output,"/seNorge2_addInfo/TEMP1d24h",sep="")
path2output.add.grd<-paste(path2output.add,"/gridded_dataset",sep="")
if (!(file.exists(path2output.main)))     dir.create(path2output.main,showWarnings=F) 
if (!(file.exists(path2output.main.stn))) dir.create(path2output.main.stn,showWarnings=F) 
if (!(file.exists(path2output.main.grd))) dir.create(path2output.main.grd,showWarnings=F) 
if (!(file.exists(path2output.add)))      dir.create(path2output.add,showWarnings=F) 
if (!(file.exists(path2output.add.grd)))  dir.create(path2output.add.grd,showWarnings=F) 
# Setup output files 
dir.create(paste(path2output.main.stn,"/",yyyymm,sep=""),showWarnings=F)
dir.create(paste(path2output.main.grd,"/",yyyymm,sep=""),showWarnings=F)
dir.create(paste(path2output.add.grd,"/",yyyymm,sep=""),showWarnings=F)
#
out.file.stn<- paste(path2output.main.stn,"/",yyyymm,
                "/seNorge_v2_0_TEMP1d24h_station_",yyyymmdd,".txt",sep="")
out.file.grd.ana<- paste(path2output.main.grd,"/",yyyymm,
                   "/seNorge_v2_0_TEMP1d24h_grid_",yyyymmdd,".nc",sep="")
out.file.grd.bck<- paste(path2output.add.grd,"/",yyyymm,
                   "/seNorge_v2_0_TEMP1d24h_grid_background_",yyyymmdd,".nc",sep="")
out.file.grd.idi<- paste(path2output.add.grd,"/",yyyymm,
                   "/seNorge_v2_0_TEMP1d24h_grid_idi_",yyyymmdd,".nc",sep="")
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
# input directories
path2input.main<-paste(main.path.output,"/seNorge2/TEMP1h",sep="")
path2input.main.stn<-paste(path2input.main,"/station_dataset",sep="")
path2input.main.grd<-paste(path2input.main,"/gridded_dataset",sep="")
if (!(file.exists(path2input.main))| 
    !(file.exists(path2input.main))| 
    !(file.exists(path2input.main))) {
  ext<-error_exit(paste("input path not found"))
}
# Setup input files
in.files.stn<-vector(length=n.timeseq)
in.files.grd.ana<-vector(length=n.timeseq)
for (i in 1:n.timeseq) {
  in.files.stn[i]<- paste(path2input.main.stn,"/",yyyymm.v[i],
                       "/seNorge_v2_0_TEMP1h_station_",yyyymmddhh.v[i],".txt",sep="")
  in.files.grd.ana[i]<- paste(path2input.main.grd,"/",yyyymm.v[i],
                           "/seNorge_v2_0_TEMP1h_grid_",yyyymmddhh.v[i],".nc",sep="")
  if ( !(file.exists(in.files.stn[i]))|
       !(file.exists(in.files.grd.ana[i])) ) {
    ext<-error_exit(paste("input file/s not found:",
                    in.files.stn[i],in.files.grd.ana[i]))
  }
}
#------------------------------------------------------------------------------
# Grid - it is defined by the DEM file
# CRS Coordinate Reference System
orog<-raster(filenamedem)
nx<-ncol(orog)
ny<-nrow(orog)
dx<-xres(orog)
dy<-yres(orog)
# 4 borders point
xmn<-xmin(orog)
xmx<-xmax(orog)
ymn<-ymin(orog)
ymx<-ymax(orog)
# South-West Point Coordinates
#Xnodesw<-xmn+dx/2
#Ynodesw<-ymn+dy/2
Xnodesw<-xmn
Ynodesw<-ymn
Xnode<-Xnodesw+(0:nx-1)*dx
Ynode<-Ynodesw+(0:ny-1)*dy
# extract all the cell values: cells[1] contains the orog[1,1] value
# Raster: cell numbers start at 1 in the upper left corner,
# and increase from left to right, and then from top to bottom
cells<-extract(orog,1:ncell(orog))
xy<-xyFromCell(orog,1:ncell(orog))
x.G<-sort(unique(xy[,1]))
y.G<-sort(unique(xy[,2]),decreasing=T)
rc<-rowColFromCell(orog,1:ncell(orog))
aux<-as.vector(cells)
mask<-which(!is.na(aux))
zgrid<-aux[mask]
xgrid<-xy[mask,1]
ygrid<-xy[mask,2]
rowgrid<-rc[mask,1]
colgrid<-rc[mask,2]
Lgrid<-length(xgrid)
print(Lgrid)
rm(xy,rc,rowgrid,colgrid,cells,aux)
#------------------------------------------------------------------------------
# [] Read Station Information 
# conditions:
# 1. stations in KDVH having: lat, lon and elevation. Note that UTM33 is 
#    obtained from lat,lon. Furthermore, the location must be in Norway or on
#    the border (less than max.Km.stnINdomain)
# 2. stations in CG
if (!testmode) {
  stations<-getStationMetadata(from.year=yyyy,to.year=yyyy,
                                   max.Km=max.Km.stnINdomain)
} else {
  stations<-read.csv(file=station.info)
}
LOBS<-length(stations$stnr)
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
z.G<-extract(orog,cbind(VecX,VecY))
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
cat(paste("year","month","day","stid","x","y","z","yo",
          "yb","ya","yav","yidi","yidiv","dqcflag","\n",sep=";"),
    file=out.file.stn,append=F)
# r is the raster structure used for map production
r <-raster(ncol=nx, nrow=ny, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx,
            crs=proj4.utm33)
r[]<-NA
# get data from KDVH
if (!testmode) {
  data<-getStationData(var="TAMRR", from.yyyy=yyyy, from.mm=mm, from.dd=dd,
                       to.yyyy=yyyy, to.mm=mm, to.dd=dd,
                       qa=NULL, statlist=stations, outside.Norway=T,
                       verbose=T)
} else {
  data<-read.csv(file=observed.data)
}
yo<-as.numeric(data$value)
yo.h.pos<-which(!is.na(yo))
#yo.pos<-which(!is.na(yo))
#------------------------------------------------------------------------------
# BACKGROUND AT STATION LOCATIONS/GRIDPOINTS
yb.tab<-matrix(nrow=LOBS,ncol=n.timeseq)
yb<-vector(length=LOBS)
for (i in 1:n.timeseq) {
# open/read/close netcdf file
  nc <- open.ncdf(in.files.grd.ana[i])
  x.data <- get.var.ncdf( nc )
  aux<-att.get.ncdf( nc, "UTM_Zone_33","proj4" )
  projstr<-aux$value
  dx<-nc$dim$X$vals[2]-nc$dim$X$vals[1]
  ex.xmin<-min(nc$dim$X$vals)-dx/2
  ex.xmax<-max(nc$dim$X$vals)+dx/2
  dy<-nc$dim$Y$vals[2]-nc$dim$Y$vals[1]
  ex.ymin<-min(nc$dim$Y$vals)-dy/2
  ex.ymax<-max(nc$dim$Y$vals)+dy/2
  nx<-nc$dim$X$len
  ny<-nc$dim$Y$len
  close.ncdf(nc)
# Define raster variable "xx"
  r <-raster(ncol=nx, nrow=ny,
              xmn=ex.xmin, xmx=ex.xmax, ymn=ex.ymin, ymx=ex.ymax,
              crs=projstr)
  r[]<-NA
# put data on raster variable (t=transpose)
  r[]<-t(x.data)
  if (i==1) st<-stack(r)
  if (i>1)  st<-stack(st,r)
  yb.tmp<-extract(r,cbind(VecX,VecY))
  table<-read.table(file=in.files.stn[i],header=T,sep=";",stringsAsFactors=F)
  table$stid<-as.integer(table$stid)
  table$ya<-as.numeric(table$ya)
  yb.tab[,i]<-table$ya[match(VecS,table$stid)]
  # if yb is NA take yb.tab from the grid
  aux<-which(is.na(yb.tab[,i]))
  gamma<--0.00649 #C/m
  yb.tab[aux,i]<-yb.tmp[aux]+gamma*(VecZ[aux]-z.G[aux])
#  for (s in yo.h.pos) {
#    pos<-which(table$stid==VecS[s])
#    if (length(pos)==1) yb.tab[s,i]<-table$ya[i]
#  }
  rm(r)
}
#> x
# [1]  1  2  3  4  5  6  7  8  9 10
#> y
#[1]  1 11 15  5  6 21 22  2
#> a
#[1] 0.1 1.1 1.5 0.5 0.6 2.1 2.2 0.2
#> a[match(x,y)]
# [1] 0.1 0.2  NA  NA 0.5 0.6  NA  NA  NA  NA
yb<-rowMeans(yb.tab)
r.b<-mean(st)
xb<-extract(r.b,mask)
#====================================================================
# Station Analysis cycle (with SCT!)
yoBad.id<-NA
yoBad.pos<-NA
isct<-0
ydqc.flag[]<-rep(0,LOBS)
ydqc.flag[-yo.h.pos]<-(-1)
yo.OKh.pos<-which(!is.na(yo) & ydqc.flag!=1)
LOBStOK<-length(yo.OKh.pos)
while (LOBStOK>1) {
  print(paste(">> Total number of observations [not NA & not ERR(so far)] =",
              LOBStOK))
# Station (CV)Analysis/(CV)IDI
  ide<-matrix(data=0,ncol=LOBStOK,nrow=LOBStOK)
  ide[row(ide)==col(ide)]=1
  InvD<-solve(D[yo.OKh.pos,yo.OKh.pos],ide)
  W<-tcrossprod(S[yo.OKh.pos,yo.OKh.pos],InvD)
  G1<-S[,yo.OKh.pos]
  K<-tcrossprod(G1,InvD)
  rm(G1)
#  ya<-yb + K %*% (yo[yo.OKh.pos]-yb[yo.OKh.pos])
  ya<-yb + tcrossprod(K,t(yo[yo.OKh.pos]-yb[yo.OKh.pos]))
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
cat(paste(yyyy,mm,dd,
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
  K<-tcrossprod(G,InvD)
  rm(G)
  xa[start:end]<-xb[start:end]+tcrossprod(K,t(yo[yo.OKh.pos]-yb[yo.OKh.pos]))
  xidi[start:end]<-rowSums(K)
  rm(K)
  i<-i+1
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# output - ANALYSIS
  r <-raster(ncol=nx, nrow=ny,
              xmn=ex.xmin, xmx=ex.xmax, ymn=ex.ymin, ymx=ex.ymax,
              crs=projstr)
r[]<-NA
r[mask]<-round(xa,1)
if (xa.flag.write) {
  nogrid.ncout(file.name=out.file.grd.ana,
               grid=t(as.matrix(r)),
               x=x.G,y=y.G,grid.type=grid.type,
               times=c(paste(yyyymmdd,"0000",sep="")),
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
# output - background
r[]<-NA
r[mask]<-round(xb,1)
if (xa.flag.write) {
  nogrid.ncout(file.name=out.file.grd.bck,
               grid=t(as.matrix(r)),
               x=x.G,y=y.G,grid.type=grid.type,
               times=c(paste(yyyymmdd,"0000",sep="")),
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
# output - IDI
r[]<-NA
r[mask]<-round(100*xidi,1)
if (xidi.flag.write) {
  nogrid.ncout(file.name=out.file.grd.idi,
               grid=t(as.matrix(r)),
               x=x.G,y=y.G,grid.type=grid.type,
               times=c(paste(yyyymmdd,"0000",sep="")),
               prod.date=prod.date,
               proj4.string=proj4.utm33,
               var.name=xidi.var.name,
               var.longname=xidi.var.longname,
               var.unit=xidi.var.unit,
               var.mv=xidi.var.mv,
               var.version=xidi.var.version,
               times.unit=xidi.times.unit,
               times.ref=xidi.times.ref,
               reference=xidi.reference,
               source.string=xidi.source.nc)
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Exit - Success
q(status=0)
