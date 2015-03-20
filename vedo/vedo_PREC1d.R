rm(list=ls())
# Libraries
library(raster)
library(rgdal)
library(ncdf)
#
# External Functions
source(paste("/disk1/projects/seNorge2/lib/SpInt_plots.R",sep=""))
#
# Graphic parameter
xlim.sw<--75000
xlim.ne<-1120000
ylim.sw<-6450000
ylim.ne<-8000000
# South Norway only
#xlim.sw<--75000
#xlim.ne<-500000
#ylim.sw<-6450000
#ylim.ne<-7150000
#-----------------------------------------------------------------------------
# [] Colors
tcol<-c("mediumorchid","plum","paleturquoise3","paleturquoise","palegreen",
        "green3","forestgreen","yellow2","darkorange","red")
tcol_ext<-c("gray","orangered4")
# day
bcol.daily<-c(0.05,0.5,3,7,10,15,20,30,40,50,60)
# month
bcol.monthly<-c(0.05,0.5,3,10,50,100,200,300,500,750,1500)
# annual
bcol.annual<-c(0.05,0.5,250,500,750,1000,1250,1500,2000,3000,4000)
bcol.idi<-c(0.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.)
#==============================================================================
arguments <- commandArgs()
arguments
file.nc<-arguments[3]
file.txt<-arguments[4]
file.out<-arguments[5]
if (length(arguments)!=5) {
  print("Error in command line arguments:")
  print("R --vanilla yyyy.mm.dd yyyy.mm.dd blacklist_current blacklist_never errobs")
  print("             begin -----> end of the accumulation period")
  quit(status=1)
}
#
# Read Geographical Information
filenamedem<-paste("/disk1/projects/seNorge2/geoinfo/seNorge2_dem_UTM33.nc",sep="")
fileborders<-paste("/disk1/projects/seNorge2/geoinfo/TM_WORLD_BORDERS_UTM33/TM_WORLD_BORDERS_UTM33-0.2.shp",sep="")
#
orog<-raster(filenamedem)
borders<-readOGR(fileborders,"TM_WORLD_BORDERS_UTM33-0.2")
#
# open/read/close netcdf file
nc <- open.ncdf(file.nc)
data <- get.var.ncdf( nc )
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
r[]<-t(data)
#
y.data<-read.table(file=file.txt,sep=";",header=T,stringsAsFactors=F)
y.data$dqcflag<-as.integer(y.data$dqcflag)
print(y.data$dqcflag)
#year;month;day;nday;stid;x;y;z;eve.lab;yo;yb;ya;yav;yidi;yidiv;dqcflag;
#
ee<-rainspatplot(x=y.data$x,#VecX[yo.ok.pos],
                 y=y.data$y,#VecY[yo.ok.pos],
                 yvar=y.data$yo,#yo[yo.ok.pos],
                 ydqc=y.data$dqcflag,
#                 xvar=x.data,
                 xvar=r,
                 xvar.orog=orog,
                 brk=bcol.daily,
                 col=tcol,
                 colext=tcol_ext,
                 legcex=2,
                 mtxt=paste("TEXT",sep=""),
                 namefileout=file.out,
                 cx=rep(1.8,(length(tcol)+2)),
                 bnd=borders,
                 xl=c(xlim.sw,xlim.ne),yl=c(ylim.sw,ylim.ne))
warnings()
#
quit()
