# get station metadata from KDVH. Only stations having defined UTM_E, UTM_N and AMSL.
getStationMetadata<-function(from.year,to.year,max.Km)
# from/to format => yyyy 
{
  require(raster)
  require(rgdal)
  # CRS strings
  proj4.wgs84<-"+proj=longlat +datum=WGS84"
  proj4.ETRS_LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
  proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
#------------------------------------------------------------------------------
# Read Station Information 
  myurl <- paste("http://klapp.oslo.dnmi.no/metnopub/production/metno?",
                 "re=16&nod=NA&ct=text/plain&ddel=dot&del=semicolon",
                 "&fy=",from.year,"&ty=",to.year,sep="")
#print(myurl)
  o.cont<-1
  while (o.cont<=10) {
    stataux<-NULL
    try(stataux <-read.table(myurl, header = TRUE,  sep = ";",
                           stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
                           encoding = "UTF-8", quote = "",na.string=-999))
# stataux column names
# DEPARTMENT;DEPT_NO;MUNICIPALITY;MUNI_NO;ST_NAME;STNR;UTM_E;UTM_N;AMSL;LAT_DEC;LON_DEC;WMO_NO
    if (length(stataux)<10) {
      print("exit with error in command:")
          print(myurl)
          o.cont<-o.cont+1
          Sys.sleep(5)
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
# Select stations having the geographical information needed
# first step: lat, lon and elevation must be present
  lat_dec<-suppressWarnings(as.numeric(stataux$LAT_DEC))
  lon_dec<-suppressWarnings(as.numeric(stataux$LON_DEC))
  z<-suppressWarnings(as.numeric(stataux$AMSL))
  indx<-which( !is.na(lat_dec) & !is.na(lon_dec) & !is.na(z) )
# second step: the location must be in Norway or on the border (lee than max.Km)
#  intermediate step: transformation in Km-coordinates ETRS_LAEA, which has a transformation 
#    less problematic than UTM33
  coord<-SpatialPoints(cbind(lon_dec[indx],lat_dec[indx]), proj4string=CRS(proj4.wgs84))
  coord.new<-spTransform(coord, CRS(proj4.ETRS_LAEA))
  xy.RR<-coordinates(coord.new)
  x<-round(xy.RR[,1],0)
  y<-round(xy.RR[,2],0)
  stations.tmp<-data.frame(matrix(nrow=length(indx),ncol=7))
  names(stations.tmp)<-c("stnr","dept_no","lat_dec","lon_dec","z","x","y")
  stations.tmp$stnr<-suppressWarnings(as.numeric(stataux$STNR[indx]))
  stations.tmp$dept_no<-suppressWarnings(as.numeric(stataux$DEPT_NO[indx]))
  stations.tmp$z<-z[indx]
  stations.tmp$y<-y
  stations.tmp$x<-x
  stations.tmp$lat_dec<-lat_dec[indx]
  stations.tmp$lon_dec<-lon_dec[indx]
  n.stn<-length(stations.tmp$stnr)
  # Norwegian stations
  aux.no<- stations.tmp$dept_no>=1 & stations.tmp$dept_no<=28
  # Norwegian stations -> Norwegian mainland!
  indx.no<-which(stations.tmp$dept_no>=1 & stations.tmp$dept_no<=20)
  Disth<-matrix(ncol=n.stn,nrow=n.stn,data=0.)
  Disth<-(outer(stations.tmp$y,stations.tmp$y,FUN="-")**2.+outer(stations.tmp$x,stations.tmp$x,FUN="-")**2.)**0.5/1000.
  aux.in<-vector(length=n.stn)
  aux.in[1:n.stn]<-F
  for (s in indx.no) {
    aux.in[which(Disth[s,]<max.Km)]<-T 
  }
  indx.in<-which(aux.in)
# final step: create the structure stations containing the selected stations in utm33 coordinates
# note: for each of the selected we will obtain an analysis value
  n<-length(indx.in)
  stations<-data.frame(matrix(nrow=n,ncol=5))
  names(stations)<-c("stnr","z","x","y","NO")
  lat_dec<-stations.tmp$lat_dec[indx.in]
  lon_dec<-stations.tmp$lon_dec[indx.in]
  coord<-SpatialPoints(cbind(lon_dec,lat_dec), proj4string=CRS(proj4.wgs84))
  coord.new<-spTransform(coord, CRS(proj4.utm33))
  xy.RR<-coordinates(coord.new)
  stations$x<-round(xy.RR[,1],0)
  stations$y<-round(xy.RR[,2],0)
  stations$z<-stations.tmp$z[indx.in]
  stations$stnr<-stations.tmp$stnr[indx.in]
  stations$NO<-aux.no[indx.in]
  rm(lat_dec,lon_dec,z,indx,coord,coord.new,xy.RR,x,y)
  rm(stations.tmp,n.stn,indx.no,Disth,aux.in,indx.in)
  return(stations)
}

# get station data from KDVH.
getStationData<-function(var=NULL, from.date, to.date, h=NULL, qa=NULL, statlist=NULL, outside.Norway=F)
# from.date/to.date -> dd.mm.yyyy
{
  str.var<-""
  str.qa<-""
  str.out<-""
  daily<-F
  if (var=="TAMRR" | var=="TAM" | var=="RR") daily<-T 
  if (daily) str.var<-"re=14"
  if (!daily) {
    str.var<-"re=17"
    if (!is.null(h)) {
      for (h.t in h) str.var<-paste(str.var,"&h=",h.t,sep="")
    }
  }
  str.var<-paste(str.var,"&p=",var,sep="")
  #date<-paste(dd,".",mm,".",yyyy,sep="")
  ulric.data<-paste("http://klapp.oslo.dnmi.no/metnopub/production/",
                    "metno?",str.var,"&fd=",from.date,"&td=",to.date,
                    "&nob=0.0&ddel=dot&del=semicolon&nmt=0",
                    "&ct=text/plain&split=1&nod=-999",sep="")
  ulric.flag<-paste("http://klapp.oslo.dnmi.no/metnopub/production/",
                    "metno?",str.var,"&fd=",from.date,"&td=",to.date,
                    "&nob=0.0&ddel=dot&del=semicolon&nmt=0",
                    "&ct=text/plain&split=1&nod=-999&flag=10",sep="")
  if (outside.Norway) {
    aux<-which(!statlist$NO)
    for (i in aux) str.out<-paste(str.out,"&s=",statlist$stnr[i],sep="")
    ulric.data.out<-paste("http://klapp.oslo.dnmi.no/metnopub/production/",
                          "metno?",str.var,"&fd=",from.date,"&td=",to.date,
                          "&nob=0.0&ddel=dot&del=semicolon&nmt=0",
                          "&ct=text/plain&split=1&nod=-999",
                          str.out,sep="")
    ulric.flag.out<-paste("http://klapp.oslo.dnmi.no/metnopub/production/",
                          "metno?",str.var,"&fd=",from.date,"&td=",to.date,
                          "&nob=0.0&ddel=dot&del=semicolon&nmt=0",
                          "&ct=text/plain&split=1&nod=-999&flag=10",
                          str.out,sep="")
  }
  print(ulric.data)
  print(ulric.flag)
  print(ulric.data.out)
  print(ulric.flag.out)
  o.cont<-1
  while (o.cont<=10) {
    o<-NULL
    q<-NULL
    o.out<-NULL
    q.out<-NULL
    try(o <- read.table(ulric.data, header = TRUE,  sep = ";",
                        stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
 	                    encoding = "UTF-8", quote = "",na.string=-999))
    try(q <- read.table(ulric.flag, header = TRUE,  sep = ";",
                        stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
 	                    encoding = "UTF-8", quote = "",na.string=-999))
    if (outside.Norway) {
      try(o.out <- read.table(ulric.data.out, header = TRUE,  sep = ";",
                              stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
 	                          encoding = "UTF-8", quote = "",na.string=-999))
      try(q.out <- read.table(ulric.flag.out, header = TRUE,  sep = ";",
                              stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
 	                          encoding = "UTF-8", quote = "",na.string=-999))
    }
    if (var=="TAMRR") value<-suppressWarnings(as.numeric(o$TAMRR))
    if (var=="TAM") value<-suppressWarnings(as.numeric(o$TAM))
    if (var=="RR") value<-suppressWarnings(as.numeric(o$RR))
    if (var=="TA") value<-suppressWarnings(as.numeric(o$TA))
    if (var=="RR_1") value<-suppressWarnings(as.numeric(o$RR_1))
    stnr.o<-suppressWarnings(as.numeric(o$Stnr))
    value[value==-999]<-NA
    indx<-which( !is.na(value) & (stnr.o %in% stations$stnr) )
    n.o<-length(indx)
    if (n.o<10) {
      print("getStationData: exit with error in command:")
      print(ulric.data)
      o.cont<-o.cont+1
      Sys.sleep(5)
    } else {
      break
    }
  }
  if (o.cont>10) {
    print("getStationData: Fatal Error in command:")
    print(ulric.data)
    q(status=1)
  }
  if (outside.Norway) {
    if (length(o.out$Stnr)>0) {
      o.tot<-merge(o,o.out,all=T)
    } else {
      o.tot<-o
    }
    if (length(q.out$Stnr)>0) {
      q.tot<-merge(q,q.out,all=T)
    } else {
      q.tot<-q
    }
  } else {
    o.tot<-o
    q.tot<-q
  }
  print("oooooooooooooooo")
  print(o)
  print("oooooooooooooooo.out")
  print(o.out)
  print("oooooooooooooooo.tot")
  print(o.tot)
  print("qqqqqqqqqqqqqqqq.tot")
  print(q.tot)
  if (var=="TAMRR") value<-suppressWarnings(as.numeric(o.tot$TAMRR))
  if (var=="TAM") value<-suppressWarnings(as.numeric(o.tot$TAM))
  if (var=="RR") value<-suppressWarnings(as.numeric(o.tot$RR))
  if (var=="TA") value<-suppressWarnings(as.numeric(o.tot$TA))
  if (var=="RR_1") value<-suppressWarnings(as.numeric(o.tot$RR_1))
  if (var=="TAMRR") flag<-suppressWarnings(as.numeric(q.tot$TAMRR))
  if (var=="TAM") flag<-suppressWarnings(as.numeric(q.tot$TAM))
  if (var=="RR") flag<-suppressWarnings(as.numeric(q.tot$RR))
  if (var=="TA") flag<-suppressWarnings(as.numeric(q.tot$TA))
  if (var=="RR_1") flag<-suppressWarnings(as.numeric(q.tot$RR_1))
  stnr.o<-suppressWarnings(as.numeric(o.tot$Stnr))
  stnr.q<-suppressWarnings(as.numeric(q.tot$Stnr))
  value[value==-999]<-NA
  indx.o<-which( !is.na(value) & (stnr.o %in% stations$stnr) )
  n.o<-length(indx.o)
  print(paste("  Total number of observations [not NA] =",n.o))
  # OBS: "d"-day observations without NAs
  if (daily) {
    datetime.o<-as.POSIXlt(strptime(paste(o.tot$Year,
                           formatC(o.tot$Month,width=2,flag="0"),
                           formatC(o.tot$Day,width=2,flag="0"),
                           sep=""),"%Y%m%d"),"UTC")
    datetime.q<-as.POSIXlt(strptime(paste(q.tot$Year,
                           formatC(q.tot$Month,width=2,flag="0"),
                           formatC(q.tot$Day,width=2,flag="0"),
                           sep=""),"%Y%m%d"),"UTC")
  } else {
    datetime.o<-as.POSIXlt(strptime(paste(o.tot$Year,
                           formatC(o.tot$Month,width=2,flag="0"),
                           formatC(o.tot$Day,width=2,flag="0"),
                           formatC(o.tot$Hour,width=2,flag="0"),
                           sep=""),"%Y%m%d%H"),"UTC")
    datetime.q<-as.POSIXlt(strptime(paste(q.tot$Year,
                           formatC(q.tot$Month,width=2,flag="0"),
                           formatC(q.tot$Day,width=2,flag="0"),
                           formatC(q.tot$Hour,width=2,flag="0"),
                           sep=""),"%Y%m%d%H"),"UTC")
  }
  print(datetime.o)
  print(datetime.q)
  datetime.seq<-as.POSIXlt(unique(datetime.o))
  n.t<-length(datetime.seq)
  n.stat<-length(statlist$stnr)
  n.data<-n.stat*n.t
  data<-data.frame(matrix(nrow=n.data,ncol=8))
  names(data)<-c("stnr","year","month","day","hour","value","DQC","KDVHflag")
  data$hour[1:n.data]<-NA
  data$value[1:n.data]<-NA
  data$DQC[1:n.data]<-NA
  data$KDVHflag[1:n.data]<-NA
  i<-0
  options(warn=1)
  for (t in 1:n.t) {
    yyyy.t<-datetime.seq$year[t]+1900
    mm.t<-datetime.seq$mon[t]+1
    dd.t<-datetime.seq$mday
    for (s in 1:n.stat) {
      i<-i+1
      data$stnr[i]<-statlist$stnr[s]
      data$year[i]<-yyyy.t
      data$month[i]<-mm.t
      data$day[i]<-dd.t
      indx.o<-which(stnr.o==data$stnr[i] & datetime.o==datetime.seq[t])
      if (length(indx.o)>0) data$value[i]<-value[indx.o]  
      indx.q<-which(stnr.q==data$stnr[i] & datetime.q==datetime.seq[t])
      if (length(indx.q)>0) {
        if (length(indx.q)>1) {
          print(paste("cazzo:",data$stnr[i]))
        } else {
          data$KDVHflag[i]<-flag[indx.q]
        }
      }
    }
  }
  print(data)
  q()
  return(data)
}
