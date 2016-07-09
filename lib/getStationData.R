# get station metadata from KDVH
getStationMetadata<-function(from.year,to.year,max.Km)
# interaction with KDVH using ULRIC
# from/to.year format => yyyy 
{
  require(raster)
  require(rgdal)
  #
  file.statauxNO<-"/lustre/storeB/users/cristianl/data/KDVHstationMetadata/KDVHstationMetadata_NorwegianStations.txt" 
  file.stataux<-"/lustre/storeB/users/cristianl/data/KDVHstationMetadata/KDVHstationMetadata_NotNorwegianStations.txt" 
  # CRS strings
  proj4.wgs84<-"+proj=longlat +datum=WGS84"
  proj4.ETRS_LAEA<-"+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
  proj4.utm33<-"+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"
  # domain definition
  lat.mn<-54.
  lat.mx<-72.
  lon.mn<-2.
  lon.mx<-33.
#------------------------------------------------------------------------------
# Read Station Information 
# Fylke code
# 1 ØSTFOLD, 2 AKERSHUS, 3 OSLO, 4 HEDMARK, 5 OPPLAND, 6 BUSKERUD, 7 VESTFOLD, 
# 8 TELEMARK, 9 AUST-AGDER, 10 VEST-AGDER, 11 ROGALAND, 12 HORDALAND, 
# 14 SOGN OG FJORDANE, 15 MØRE OG ROMSDAL, 16 SØR-TRØNDELAG, 
# 17 NORD-TRØNDELAG, 18 NORDLAND, 19 TROMS, 20 FINNMARK,
#  myurlno <- paste("http://klapp.oslo.dnmi.no/metnopub/production/metno?",
#                   "re=16&nod=NA&ct=text/plain&ddel=dot&del=semicolon",
#                   "&fy=",from.year,"&ty=",to.year,
#                   "&cnr=1&cnr=2&cnr=3&cnr=4&cnr=5&cnr=6&cnr=7&cnr=8&cnr=9",
#                   "&cnr=10&cnr=11&cnr=12&cnr=14&cnr=15&cnr=16&cnr=17",
#                   "&cnr=18&cnr=19&cnr=20&cnr=13",sep="")
#  myurl <- paste("http://klapp.oslo.dnmi.no/metnopub/production/metno?",
#                 "re=16&nod=NA&ct=text/plain&ddel=dot&del=semicolon",
#                 "&fy=",from.year,"&ty=",to.year,sep="")
# Read files and not make a query to the database because of the problem in getting
#  the data from KDVH.
  o.cont<-1
  while (o.cont<=10) {
    statauxNO<-NULL
    stataux<-NULL
    if (file.exists(file.statauxNO)) {
      print(paste("Sys.getlocale(\"LC_CTYPE\")=",Sys.getlocale("LC_CTYPE")))
      if (Sys.getlocale("LC_CTYPE")=="C") {
        Sys.setlocale(category = "LC_CTYPE", locale = "en_US.UTF-8")
        print(paste("Sys.getlocale(\"LC_CTYPE\")=",Sys.getlocale("LC_CTYPE")))
      }
      enc.str.list<-c("Latin1","ISO-8859-1","ISO-8859-15","latin-9")
      for (enc.str in enc.str.list) {
        print(enc.str)
        deb<-try( statauxNO <-read.table(file=file.statauxNO,
                                         header=TRUE,sep=";",
                                         stringsAsFactors=FALSE,
                                         na.strings="NA",strip.white=T,
                                         fileEncoding=enc.str,
                                         quote="",skipNul=T))
        if (attr(deb,"class")!="try-error") break
      }
      if (attr(deb,"class")=="try-error") quit(status=1)
    }
    if (file.exists(file.stataux)) {
#      stataux   <-read.table(file=file.stataux,
#                             header=TRUE,sep=";",
#                             stringsAsFactors=FALSE,
#                             fileEncoding="Latin1",
#                             quote="",na.strings="NA",
#                             strip.white=T)
      print(paste("Sys.getlocale(\"LC_CTYPE\")=",Sys.getlocale("LC_CTYPE")))
      if (Sys.getlocale("LC_CTYPE")=="C") {
        Sys.setlocale(category = "LC_CTYPE", locale = "en_US.UTF-8")
        print(paste("Sys.getlocale(\"LC_CTYPE\")=",Sys.getlocale("LC_CTYPE")))
      }
      enc.str.list<-c("Latin1","ISO-8859-1","ISO-8859-15","latin-9")
      for (enc.str in enc.str.list) {
        print(enc.str)
        deb<-try( stataux <-read.table(file=file.stataux,
                                         header=TRUE,sep=";",
                                         stringsAsFactors=FALSE,
                                         na.strings="NA",strip.white=T,
                                         fileEncoding=enc.str,
                                         quote="",skipNul=T))
        if (attr(deb,"class")!="try-error") break
      }
      if (attr(deb,"class")=="try-error") quit(status=1)
    }
# column names
# DEPARTMENT;DEPT_NO;MUNICIPALITY;MUNI_NO;ST_NAME;STNR;UTM_E;UTM_N;AMSL;LAT_DEC;LON_DEC;WMO_NO
    if (is.null(statauxNO) | is.null(stataux)) {
      print("Fatal error: KDVH Station Metadata files not found")
      print(paste("    Norwegian stations file = ",file.statauxNO))
      print(paste("Not Norwegian stations file = ",file.stataux))
      o.cont<-o.cont+1
      Sys.sleep(round(runif(1)*10))
      print(warnings())
    } else if (nrow(stataux)<10 | nrow(statauxNO)<10) {
      print("Fatal error: KDVH Station Metadata files found but empty")
      print(paste("    Norwegian stations file = ",file.statauxNO,"number of rows=",nrow(statauxNO)))
      print(paste("Not Norwegian stations file = ",file.stataux,"number of rows=",nrow(stataux)))
      o.cont<-o.cont+1
      Sys.sleep(round(runif(1)*10))
      print(warnings())
    } else {
      break
    }
  }
  if (o.cont>10) {
    print("Fatal Error: Can't read the KDVH metadata files")
    print(paste("    Norwegian stations file = ",file.statauxNO))
#    print("@@statauxNO")
#    print(statauxNO)
    print(paste("Not Norwegian stations file = ",file.stataux))
#    print("@@stataux")
#    print(stataux)
    print(warnings())
    q(status=1)
  }
# from this point on
  nn<-nrow(stataux)
  nn.no<-nrow(statauxNO)
insertRow <- function(existingDF, newrow, r) {
  if (nrow(existingDF)>=r) {
    existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
    existingDF[r,] <- newrow
  } else {
    existingDF[r,] <- newrow
  }
  existingDF
}
  for (i in 1:nn.no) {
    if (!(statauxNO$STNR[i] %in% stataux$STNR)) {
      stataux<-insertRow(stataux, statauxNO[i,], nn+1)
      nn<-nn+1
      if (is.na(stataux$DEPT_NO[nn])) stataux$DEPT_NO[nn]<-1
    }
  }
# Select stations having the geographical information needed
# first step: lat, lon and elevation must be present
  lat_dec<-suppressWarnings(as.numeric(stataux$LAT_DEC))
  lon_dec<-suppressWarnings(as.numeric(stataux$LON_DEC))
  z<-suppressWarnings(as.numeric(stataux$AMSL))
  indx1<-!is.na(lat_dec) & !is.na(lon_dec) & !is.na(z)
  indx2<-lon_dec>lon.mn & lon_dec<lon.mx & lat_dec>lat.mn & lat_dec<lat.mx
  indx2[is.na(indx2)]<-F
  indx<-which(indx1 & indx2)
  print(paste("total number of stations in KDVH=",length(lat_dec)))
  print(paste("total number of stations in the predefined domain=",length(indx)))
# debug
#  print(paste("lon_dec[indx],lat_dec[indx]\n"))
#  print(paste(lon_dec[indx],lat_dec[indx],"\n"))
# second step: the location must be in Norway or on the border (lee than max.Km)
#  intermediate step: transformation in Km-coordinates ETRS_LAEA, which has a transformation 
#    less problematic than UTM33
  coord<-SpatialPoints(cbind(lon_dec[indx],lat_dec[indx]), proj4string=CRS(proj4.wgs84))
  coord.new<-spTransform(coord, CRS(proj4.ETRS_LAEA))
  xy.RR<-coordinates(coord.new)
  x<-as.numeric(round(xy.RR[,1],0))
  y<-as.numeric(round(xy.RR[,2],0))
# debug
#  print(paste("x,y,dept\n"))
#  print(paste(x,y,stataux$DEPT_NO[indx],"\n"))
#
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
  print(paste("total number of Norwegian stations in the predefined domain=",length(indx.no)))
  Disth<-matrix(ncol=n.stn,nrow=n.stn,data=0.)
  Disth<-(outer(stations.tmp$y,stations.tmp$y,FUN="-")**2.+outer(stations.tmp$x,stations.tmp$x,FUN="-")**2.)**0.5/1000.
  aux.in<-vector(length=n.stn)
  aux.in[1:n.stn]<-F
  for (s in indx.no) {
    aux.in[which(Disth[s,]<=max.Km)]<-T 
  }
  indx.in<-which(aux.in)
# final step: create the structure stations containing the selected stations in utm33 coordinates
# note: for each of the selected we will obtain an analysis value
  n<-length(indx.in)
  print(paste("total number of stations within the domain of interest (predefined domain AND not too distant from a Norwegian station)=",n))
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

# get station data from KDVH by using Ulric web-based tool
getStationData<-function(var=NULL, from.dd, from.mm, from.yyyy, from.hh=NULL,
                         to.dd, to.mm, to.yyyy, to.hh,
                         h=NULL, qa=NULL, statlist=NULL, outside.Norway=F,
                         err.file=NULL, blist=NULL,
                         fun=NULL, verbose=F, 
                         val.min.allowed=NULL, val.max.allowed=NULL)
# Ulric doc, see (http://metklim.met.no/klima/userservices/urlinterface/brukerdok)
# >> input:
#  var: "TAMRR","TAM","TA","RR","RR_1"
#  statlist: station list, if NULL then call getStationMetadata with the current year and max.Km=0
#  from.date/to.date: from.yyyy from.mm from.dd / to.yyyy to.mm to.dd 
#  h: hours (vector) 
#  from.hh and to.hh: not used
#  qa: not used
#  outside.Norway: TRUE use observations outside Norway, FALSE do not use them
#  err.file: path to file with suspect/erroneous observations (supplementary DQC)
#  blist: path to file with a list of stations which must not enter the spatial interpolation procedure
#  fun: function applied to observed data. implemented: "sum"
#  val.min/max.allowed: DQC on plausible range of values
# >> output (follows the statlist order for stations):
#  stnr: station number
#  year month day: date
#  h: hour ("NA" in case of daily variables)
#  note: if fun="sum" then the date is the date of the first step (older)
#  ntime: number of timestamp requested
#  value: observed value
#  nvalue: number of observed values (used in fun)
#  DQC: empty, set to NA (used after the call at this function to create a definitive DQC flag)
#  KDVHflag: DQC flag from KDVH. in fun="sum" then KDVHflag is set to 100 if at least one of the flags is >2
#  plausible: plausibility-check result (T=plausible, F=not plausible)
#  err.ext: DQC flag from external files (T=err, F=good)
#  blist: T if station is blacklisted
#    note: if fun="sum" then 
#     plausible=F if at least one of the observation has plausible=F
#     err.ext,blist=T if at least one of the observation has err.ext,blist=T
#================================================================================
{
  max.try<-10
  sleep.int<-5
  names.daily<-c("stnr","year","month","day","value")
  names.hourly<-c("stnr","year","month","day","hour","value")
  data.names<-c("stnr","year","month","day","hour","ntime",
                "value","nvalue","DQC",
                "KDVHflag",
                "plausible","err.ext",
                "blist")
  data.names.col<-length(data.names)
# check input information
  from.date<-paste(from.dd,".",from.mm,".",from.yyyy,sep="")
  to.date<-paste(to.dd,".",to.mm,".",to.yyyy,sep="")
#
  valid.var<-c("TAMRR","TAM","TA","RR","RR_1")
  if (!var %in% valid.var) {
    data<-data.frame(matrix(nrow=1,ncol=data.names.col))
    names(data)<-data.names
    data$stnr[1]<--100
    return(data)
  }
  if (is.null(statlist)) {
    yyyy<-format(Sys.time(), "%Y")
    statlist<-getStationMetadata(from.year=yyyy,to.year=yyyy,max.Km=0)
  }
  err.file.ok<-F
  if (!is.null(err.file)) {
    if (file.exists(err.file)) {
      errobs<-read.table(file=err.file,header=T,sep=";")
      names(errobs)<-c("stnr","year","month","day","hour","value")
      err.timeseq<-as.POSIXlt(strptime(paste(errobs$year,
                              formatC(errobs$month,width=2,flag="0"),
                              formatC(errobs$day,width=2,flag="0"),
                              formatC(errobs$hour,width=2,flag="0"),sep=""),
                              "%Y%m%d%H","UTC"),"UTC")
      err.stnr<-unique(errobs$stnr)
      err.file.ok<-T
    } else {
      print("Warning: file not found")
      print(file_errobs)
    }
  }
  blist.file.ok<-F
  if (!is.null(blist)) {
    if (file.exists(blist)) {
      blacklist<-read.table(file=blist,header=T,sep=";")
      names(blacklist)<-c("stnr")
      blist.file.ok<-T
    } else {
      print("Warning: file not found")
      print(blist)
    }
  }
#
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
  o.cont<-1
  while (o.cont<=max.try) {
    o<-NULL
    q<-NULL
    o.out<-NULL
    q.out<-NULL
    if (verbose) print("query under processing: data inside Norway")
    if (verbose) print(ulric.data)
    try(o <- read.table(ulric.data, header = TRUE,  sep = ";",
                        stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
 	                    encoding = "UTF-8", quote = "",na.string=-999))
#    print(o)
    if (length(o$Stnr)==0) {
      if (verbose) {
        print("getStationData: exit with error in command:")
        print(ulric.data)
      }
      o.cont<-o.cont+1
      Sys.sleep(sleep.int)
      next
    }
    if (verbose) print("query under processing: quality flags inside Norway")
    if (verbose) print(ulric.flag)
    try(q <- read.table(ulric.flag, header = TRUE,  sep = ";",
                        stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
 	                    encoding = "UTF-8", quote = "",na.string=-999))
#    print(q)
    if (outside.Norway) {
      if (verbose) print("query under processing: data outside Norway")
      if (verbose) print(ulric.data.out)
      try(o.out <- read.table(ulric.data.out, header = TRUE,  sep = ";",
                              stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
 	                          encoding = "UTF-8", quote = "",na.string=-999))
#      print(o.out)
      if (verbose) print("query under processing: quality flags outside Norway")
      if (verbose) print(ulric.flag.out)
      try(q.out <- read.table(ulric.flag.out, header = TRUE,  sep = ";",
                              stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
 	                          encoding = "UTF-8", quote = "",na.string=-999))
#      print(q.out)
    }
    if (verbose) print("all queries done")
    if ( daily) {
      names(o)<-names.daily
      if (length(q$Stnr)>0) names(q)<-names.daily
      if (length(o.out$Stnr)>0) names(o.out)<-names.daily
      if (length(q.out$Stnr)>0) names(q.out)<-names.daily
    } else {
      names(o)<-names.hourly
      if (length(q$Stnr)>0) names(q)<-names.hourly
      if (length(o.out$Stnr)>0) names(o.out)<-names.hourly
      if (length(q.out$Stnr)>0) names(q.out)<-names.hourly
    }
    value<-suppressWarnings(as.numeric(o$value))
    stnr.o<-suppressWarnings(as.numeric(o$stnr))
    value[value==-999]<-NA
    indx<-which( !is.na(value) & (stnr.o %in% stations$stnr) )
    n.o<-length(indx)
    if (n.o<max.try) {
      if (verbose) {
        print("getStationData: exit with error in command:")
        print(ulric.data)
      }
      o.cont<-o.cont+1
      Sys.sleep(sleep.int)
    } else {
      break
    }
  }
  if (o.cont>max.try) {
    if (verbose) {
      print("getStationData: Fatal Error in command:")
      print(ulric.data)
    }
    data<-data.frame(matrix(nrow=1,ncol=data.names.col))
    names(data)<-data.names
    data$stnr[1]<--200
    return(data)
  }
  if (outside.Norway) {
    if (length(o.out$stnr)>0) {
      o.tot<-merge(o,o.out,all=T)
    } else {
      o.tot<-o
    }
    if (length(q.out$stnr)>0) {
      q.tot<-merge(q,q.out,all=T)
    } else {
      q.tot<-q
    }
  } else {
    o.tot<-o
    q.tot<-q
  }
#  if ( daily) {
#    names(o.tot)<-c("stnr","year","month","day","hour","value")
#    names(q.tot)<-c("stnr","year","month","day","hour","flag")
#  }
  value<-suppressWarnings(as.numeric(o.tot$value))
  flag<-suppressWarnings(as.numeric(q.tot$value))
  stnr.o<-suppressWarnings(as.numeric(o.tot$stnr))
  stnr.o.unique<-unique(stnr.o)
  indx.stn<-which(statlist$stnr %in% stnr.o)
  stnr.q<-suppressWarnings(as.numeric(q.tot$stnr))
  value[value==-999]<-NA
  indx.o<-which( !is.na(value) & (stnr.o %in% stations$stnr) )
  n.o<-length(indx.o)
  if (verbose) print(paste("  Total number of observations [not NA] =",n.o))
  # OBS: "d"-day observations without NAs
  if (daily) {
    datetime.o<-as.POSIXlt(strptime(paste(o.tot$year,
                           formatC(o.tot$month,width=2,flag="0"),
                           formatC(o.tot$day,width=2,flag="0"),
                           sep=""),"%Y%m%d","UTC"),"UTC")
    datetime.q<-as.POSIXlt(strptime(paste(q.tot$year,
                           formatC(q.tot$month,width=2,flag="0"),
                           formatC(q.tot$day,width=2,flag="0"),
                           sep=""),"%Y%m%d","UTC"),"UTC")
  } else {
    datetime.o<-as.POSIXlt(strptime(paste(o.tot$year,
                           formatC(o.tot$month,width=2,flag="0"),
                           formatC(o.tot$day,width=2,flag="0"),
                           formatC(o.tot$hour,width=2,flag="0"),
                           sep=""),"%Y%m%d%H","UTC"),"UTC")
    datetime.q<-as.POSIXlt(strptime(paste(q.tot$year,
                           formatC(q.tot$month,width=2,flag="0"),
                           formatC(q.tot$day,width=2,flag="0"),
                           formatC(q.tot$hour,width=2,flag="0"),
                           sep=""),"%Y%m%d%H","UTC"),"UTC")
  }
#  print("datetime.o")
#  print(datetime.o)
#  print("datetime.q")
#  print(datetime.q)
  datetime.seq<-as.POSIXlt(unique(as.POSIXct(datetime.o,"UTC")),"UTC")
#  print("datetime.seq")
#  print(datetime.seq)
  n.t<-length(datetime.seq)
  n.stat<-length(statlist$stnr)
  n.data<-n.stat*n.t
  data<-data.frame(matrix(nrow=n.data,ncol=data.names.col))
  names(data)<-data.names
  data$stnr[1:n.data]<-NA
  data$year[1:n.data]<-NA
  data$month[1:n.data]<-NA
  data$day[1:n.data]<-NA
  data$hour[1:n.data]<-NA
  if (is.null(fun)) data$ntime[1:n.data]<-1
  if (!is.null(fun)) data$ntime[1:n.data]<-NA
  data$nvalue[1:n.data]<-NA
  data$value[1:n.data]<-NA
  data$DQC[1:n.data]<-NA
  data$KDVHflag[1:n.data]<-NA
  data$plausible[1:n.data]<-T
  data$err.ext[1:n.data]<-F
  data$blist[1:n.data]<-F
  options(warn=1)
  for (t in 1:n.t) {
    yyyy.t<-datetime.seq$year[t]+1900
    mm.t<-datetime.seq$mon[t]+1
    dd.t<-datetime.seq$mday[t]
    if (!daily) hh.t<-datetime.seq$hour[t]
    aux.d.o<-datetime.o==datetime.seq[t]
    aux.d.q<-datetime.q==datetime.seq[t]
    if (err.file.ok) aux.d.e<-err.timeseq==datetime.seq[t]
    t.b<-(t-1)*n.stat+1
    t.e<-t*n.stat
    data$stnr[t.b:t.e]<-statlist$stnr
    data$year[t.b:t.e]<-yyyy.t
    data$month[t.b:t.e]<-mm.t
    data$day[t.b:t.e]<-dd.t
    if (!daily) data$hour[t.b:t.e]<-hh.t
    for (s in indx.stn) {
      s.t<-s+t.b-1
      indx.o<-which(stnr.o==data$stnr[s.t] & aux.d.o)
      data$nvalue[s.t]<-length(indx.o)
      if (data$nvalue[s.t]>0) { 
        if (data$nvalue[s.t]==1) {
          data$value[s.t]<-value[indx.o]
        } else {
          if (verbose) print(paste("anomaly in retrieving ",data$stnr[s.t]))
        }
      }
      indx.q<-which(stnr.q==data$stnr[s.t] & aux.d.q)
      if (length(indx.q)>0) {
        if (length(indx.q)==1) {
          data$KDVHflag[s.t]<-flag[indx.q]
        } else {
          if (verbose) print(paste("anomaly in retrieving ",data$KDVHflag[s.t]))
        }
      }
      if (err.file.ok) if (any(errobs$stnr==data$stnr[s.t] & aux.d.e)) data$err.ext[s.t]<-T
      if (blist.file.ok) if (any(blacklist$stnr==data$stnr[s.t])) data$blist[s.t]<-T
    }
  }
  if (!is.null(val.min.allowed)) data$plausible[data$value<val.min.allowed]<-F
  if (!is.null(val.max.allowed)) data$plausible[data$value>val.max.allowed]<-F
  if (!is.null(fun)) {
    if (fun=="sum") {
      data.fun<-data.frame(matrix(nrow=n.stat,ncol=data.names.col))
      names(data.fun)<-data.names
      data.fun$stnr[1:n.stat]<-NA
      data.fun$year[1:n.stat]<-from.yyyy
      data.fun$month[1:n.stat]<-from.mm
      data.fun$day[1:n.stat]<-from.dd
      if (!daily & !is.null(from.hh)) data.fun$hour[1:n.stat]<-from.hh
      if (!daily & is.null(from.hh)) data.fun$hour[1:n.stat]<-min(h)
      if (daily)  data.fun$hour[1:n.stat]<-NA
      data.fun$ntime[1:n.stat]<-n.t
      data.fun$nvalue[1:n.stat]<-NA
      data.fun$value[1:n.stat]<-NA
      data.fun$DQC[1:n.stat]<-NA
      data.fun$KDVHflag[1:n.stat]<-0
      data.fun$plausible[1:n.stat]<-T
      data.fun$err.ext[1:n.stat]<-F
      data.fun$blist[1:n.stat]<-F
      for (s in 1:n.stat) {
        data.fun$stnr[s]<-statlist$stnr[s]
        indx.o<-which(data$stnr==data.fun$stnr[s])
        data.fun$nvalue[s]<-length(which(!is.na(data$value[indx.o])))
        if (data.fun$nvalue[s]>0)  data.fun$value[s]<-sum(data$value[indx.o],na.rm=T)
        if (any(data$KDVHflag[indx.o]>2,na.rm=T)) data.fun$KDVHflag[s]<-100
        if (any(data$err.ext[indx.o],na.rm=T)) data.fun$err.ext[s]<-T
        if (any(data$blist[indx.o],na.rm=T)) data.fun$blist[s]<-T
        if (any(!data$plausible[indx.o],na.rm=T)) data.fun$plausible[s]<-F
      }
    }
    return(data.fun)
  } else {
    return(data)
  }
}
