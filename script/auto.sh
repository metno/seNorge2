# Variables
  R=/usr/bin/R
  Bspat=/home/cristianl/projects/Bspat/Bspat_RR1/Bspat_RR1_v1_0.R
  LOGDIR=/home/cristianl/seNorge2_driverscripts/PREC_hourly/log
#---------------------------------
# [] Read command line arguments 
#---------------------------------
  flagDATESTART=0
  flagDATEEND=0
  while getopts "s:e:" Option
  do
    case $Option in
    s ) flagDATESTART=1
    DATESTART=$OPTARG
    ;;
    e ) flagDATEEND=1
    DATEEND=$OPTARG
    ;;
    * ) echo " not recognized Option ";;
    esac
  done
# checks
  if [ "$flagDATESTART" -eq 0 ] || [ "$flagDATEEND" -eq 0 ]
  then
    echo "error in input date"
    exit 1
  fi
# log
  DATAC=`date +%Y%m%d%H%M`
  echo "seNorge2_PREC_hourly.sh "`date +%Y-%m-%d" "%H:%M`" > elaborations from "$DATESTART" UTC to "$DATEEND" UTC"
#--------------------------------------------------
# [] Clean directories
#--------------------------------------------------
  echo "seNorge2_PREC_hourly.sh "`date +%Y-%m-%d" "%H:%M`" > clean directories"
#-----------------------------------------------------
# [] Elaboration
#-----------------------------------------------------
  echo "seNorge2_PREC_hourly.sh "`date +%Y-%m-%d" "%H:%M`" > ==============================================================="
  echo "seNorge2_PREC_hourly.sh "`date +%Y-%m-%d" "%H:%M`" > begin main"
  YYYYbeg=${DATESTART:0:4}
  MMbeg=${DATESTART:5:2}
  DDbeg=${DATESTART:8:2}
  HHbeg=${DATESTART:11:2}
  YYYYend=${DATEEND:0:4}
  MMend=${DATEEND:5:2}
  DDend=${DATEEND:8:2}
  HHend=${DATEEND:11:2}
  SECbeg=`date +%s -d "$YYYYbeg-$MMbeg-$DDbeg $HHbeg:00:00"`
  SECend=`date +%s -d "$YYYYend-$MMend-$DDend $HHend:00:00"`
  SECcur=$SECbeg
  while (( "$SECcur" <= "$SECend" )) 
  do
    DATEcur=`date --date="1970-01-01 $SECcur sec UTC" +%Y.%m.%d.%H`
    echo "============================================================================="
    echo "$R --vanilla $DATEcur $DATEcur /disk1/projects/seNorge2/etc/SpInt_PREC_daily_BLACKLIST_never.txt /disk1/projects/seNorge2/etc/SpInt_PREC_daily_BLACKLIST_current.txt /disk1/projects/seNorge2/etc/SpInt_PREC_daily_ERRONEOUSOBSERVATIONS.txt /disk1/projects/seNorge2/etc/config_list.r pc4436 < /disk1/projects/seNorge2/Bspat_PREC1hRT/Bspat_PREC1hRT.R"
    $R --vanilla $DATEcur $DATEcur /disk1/projects/seNorge2/etc/SpInt_PREC_daily_BLACKLIST_never.txt /disk1/projects/seNorge2/etc/SpInt_PREC_daily_BLACKLIST_current.txt /disk1/projects/seNorge2/etc/SpInt_PREC_daily_ERRONEOUSOBSERVATIONS.txt /disk1/projects/seNorge2/etc/config_list.r pc4436 < /disk1/projects/seNorge2/Bspat_PREC1hRT/Bspat_PREC1hRT.R
    SECcur=$(( SECcur+3600 ))
  done
  echo "seNorge2_PREC_hourly.sh "`date +%Y-%m-%d" "%H:%M`" > end main"
#--------------------------
# [] Exit
#--------------------------
  echo "seNorge2_PREC_hourly.sh "`date +%Y-%m-%d" "%H:%M`" > success!"
  exit 0

