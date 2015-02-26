# get station metadata from KDVH. Only stations having defined UTM_E, UTM_N and AMSL.
getStationMetadata<-function(from.year,to.year)
# from/to format => yyyy 
{
#------------------------------------------------------------------------------
# Read Station Information 
#myurl <- paste("http://klapp.oslo.dnmi.no/metnopub/production/metno?",
#               "re=16&nod=NA&ct=text/plain&ddel=dot&del=semicolon",
#               "&fy=",yyyymmdd,"&ty=",yyyymmdd,sep="")
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
# stataux column names
# DEPARTMENT;DEPT_NO;MUNICIPALITY;MUNI_NO;ST_NAME;STNR;UTM_E;UTM_N;AMSL;LAT_DEC;LON_DEC;WMO_NO
# Select stations having the geographical information needed
  indx<-which( !is.na(as.numeric(stataux$UTM_E)) & 
               !is.na(as.numeric(stataux$UTM_N)) & 
               !is.na(as.numeric(stataux$AMSL)) )
  stations<-data.frame(matrix(nrow=length(indx),ncol=4))
  names(stations)<-c("Stnr","z","x","y")
  stations$Stnr<-as.numeric(stataux$STNR[indx])
  stations$z<-stataux$AMSL[indx]
  stations$x<-stataux$UTM_E[indx]
  stations$y<-stataux$UTM_N[indx]
  return(stations)
}

getStnr.outside.Norway<-function(from.year, to.year, qa=NULL, Km=50) 
{
  myurl <- paste("http://klapp.oslo.dnmi.no/metnopub/production/metno?",
                 "re=16&nod=NA&ct=text/plain&ddel=dot&del=semicolon",
                 "&fy=",from.year,"&ty=",to.year,sep="")
#  print(myurl)
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
# stataux column names
# DEPARTMENT;DEPT_NO;MUNICIPALITY;MUNI_NO;ST_NAME;STNR;UTM_E;UTM_N;AMSL;LAT_DEC;LON_DEC;WMO_NO
# Select stations having the geographical information needed
  L<-length(stataux$STNR)
  stations<-data.frame(matrix(nrow=L,ncol=4))
  names(stations)<-c("Stnr","z","x","y")
  stations$Stnr<-as.numeric(stataux$STNR)
  stations$z<-stataux$AMSL
  stations$x<-stataux$UTM_E
  stations$y<-stataux$UTM_N



return(stations)

}

# get station data from KDVH.
getStationData<-function(vars=NULL, from.date, to.date, qa=NULL, outside.Norway=F, outside.Norway.sel=NULL)
# from.date/to.date -> dd.mm.yyyy
{
  str.var<-""
  str.qa<-""
  if (length(grep(vars,"TAMRR"))>0) str.var<-"re=14&p=TAMRR"
  if (!is.null(qa)) str.qa<-"&qa=2"
  #date<-paste(dd,".",mm,".",yyyy,sep="")
  ulric<-paste("http://klapp.oslo.dnmi.no/metnopub/production/",
               "metno?",str.var,"&fd=",from.date,"&td=",to.date,
               "&nob=0.0&ddel=dot&del=semicolon&nmt=0",
               "&ct=text/plain&split=1&nod=-999",str.qa,sep="")
  #print(ulric)
  o.cont<-1
  while (o.cont<=10) {
    o<-NULL
    try(o <- read.table(ulric, header = TRUE,  sep = ";", #nrows = nrows,
            stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1",
 	            encoding = "UTF-8", quote = "",na.string=-999))
    indx<-which(o$TA!=-999)
    LOBS<-length(indx==T)
    if (LOBS<10) {
      print("exit with error in command:")
      print(ulric)
      o.cont<-o.cont+1
      Sys.sleep(5)
    } else {
      break
    }
  }
  if (o.cont>10) {
    print("Fatal Error in command:")
    print(ulric)
    q(status=1)
  }
  # o: column names
  if (length(grep(vars,"TAMRR"))>0) indx<-which(o$TAMRR!=-999)
  LOBS<-length(indx)
  print(paste("  Total number of observations [not NA] =",LOBS))
  if (LOBS==0) return()
  # OBS: "d"-day observations without NAs
  data<-data.frame(matrix(nrow=LOBS,ncol=7))
  names(data)<-c("stnr","year","month","day","hour","value","DQC")
  data$stnr<-as.numeric(o$Stnr[indx])
  data$year<-o$Year[indx]
  data$month<-o$Month[indx]
  data$day<-o$Day[indx]
  if (length(grep(vars,"TAMRR"))>0) data$hour[]<-NA
  data$value<-as.numeric(o$TAM[indx])
  # DQC flag: -1 missing, 0 good obs, 1 erroneous obs
  data$DQC<-rep(-1,LOBS)
  if (outside.Norway) {
    ulric<-paste("http://klapp.oslo.dnmi.no/metnopub/production/",
                 "metno?",str.var,"&fd=",from.date,"&td=",to.date,
                 "&nob=0.0&ddel=dot&del=semicolon&nmt=0",
                 "&ct=text/plain&split=1&nod=-999",str.qa,
                 outside.Norway.sel,sep="")
  }
  return(data)
}
