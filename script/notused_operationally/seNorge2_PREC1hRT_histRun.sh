#!/bin/bash
#===============================================================================
# <seNorge2_PREC1hRT.sh>
#
# DESCRIPTION:
# ===========
# Spatial Interpolation of hourly accumulated precipitation (Real-Time mode).
#
# COMMAND LINE:
# =============
#  >seNorge2_PREC1hRT.sh -s yyyy.mm.dd.hh  (date start)
#                        -e yyyy.mm.dd.hh  (date end)
#                        -c config.file (configuration file)
#                        -p config.par  (parameter name in configuration file)
#                        -l log.directory (log directory for Bspat)
#===============================================================================
function trim()
{
    local var=$1;
    var="${var#"${var%%[![:space:]]*}"}";   # remove leading whitespace characters
    var="${var%"${var##*[![:space:]]}"}";   # remove trailing whitespace characters
    var="${var%%,}";   # remove trailing comma characters
    echo -n "$var";
}
# whereis R?
  R=/usr/bin/R
#----------------------------
# Read command line arguments
#----------------------------
  flag_s=0
  flag_e=0
  flag_c=0
  flag_p=0
  flag_l=0
  while getopts "s:e:c:p:l:b:o:" Option
  do
    case $Option in
    s ) flag_s=1
    DATESTART=$OPTARG
    ;;
    e ) flag_e=1
    DATEEND=$OPTARG
    ;;
    c ) flag_c=1
    CONFIG_FILE=$OPTARG
    ;;
    p ) flag_p=1
    CONFIG_PAR=$OPTARG
    ;;
    l ) flag_l=1
    LOGDIR=$OPTARG
    ;;
    b ) BLACKL=$OPTARG 
    ;;
    o ) ERROBS=$OPTARG 
    ;;
    * ) echo " not recognized Option ";;
    esac
  done
#------------------------------------------------------------------------------
  # time-related operations
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
#------------------------------------------------------------------------------
# checks
  if [ "$flag_s" -eq 0 ] || [ "$flag_e" -eq 0 ] || [ "$SECbeg" -gt "$SECend" ]
  then
    echo "error in dates"
    exit 1
  fi
  if [ "$flag_c" -eq 0 ] || [ "$flag_p" -eq 0 ]
  then
    echo "error in configuration file specification"
    exit 1
  fi
  if [ "$flag_l" -eq 0 ]
  then
    echo "error in log file specification"
    exit 1
  fi
#------------------------------------------------------------------------------
# read config file and setup path variable
  flag_var=0
  while read line; do
    name=$(trim "${line%%=*}");
    value=$(trim "${line#*=}");
    aux="\"$CONFIG_PAR\""
    aux1="main.path"
    if [ "$value" == $aux ]; then
      flag_var=1
    fi
    if [ "$flag_var" -eq 1 ] && [ "$name" == $aux1 ]; then
        MAINDIR=$value
        break
      fi
  done < $CONFIG_FILE
  MAINDIR=${MAINDIR:1:-1}
#------------------------------------------------------------------------------
# log
  echo "seNorge2_PREC1hRT.sh "`date +%Y-%m-%d" "%H:%M`" > elaborations from "$DATESTART" UTC to "$DATEEND" UTC"
  echo "configuration file: "$CONFIG_FILE" configuration parameter:"$CONFIG_PAR
  echo "main directory:"$MAINDIR
  vis="vis-m1"
  if [ "$CONFIG_PAR" == "$vis" ]; then
    export R_LIBS=/home/cristianl/programs/spatial_interpolation/lib/rpacks
    echo $R_LIBS
  fi
#------------------------------------------------------------------------------
# Variables
  Bspat=$MAINDIR/Bspat_PREC1hRT/Bspat_PREC1hRT.R
#  BLACKL=$MAINDIR/etc/blacklists/seNorge2_PREC1h_blacklist.txt
#  ERROBS=$MAINDIR/etc/suspect_observations/seNorge2_PREC1h_suspect_observations.txt
#--------------------------------------------------
#  Clean temporary directories, if needed 
#--------------------------------------------------
#  echo "seNorge2_PREC1hRT.sh "`date +%Y-%m-%d" "%H:%M`" > clean directories"
#  find $OUTDIR/ -mtime +1 -exec rm -vf {} \;
#-----------------------------------------------------
# Statistical Interpolation 
#-----------------------------------------------------
  SECcur=$SECbeg
  while (( "$SECcur" <= "$SECend" )) 
  do
    DATEcur=`date --date="1970-01-01 $SECcur sec UTC" +%Y.%m.%d.%H`
    echo "================================================================================"
    echo "`date +%Y-%m-%d" "%H:%M` > $R --vanilla $DATEcur $DATEcur $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC1hRT_$DATEcur_$DATEcur.log 2>&1"
    $R --vanilla $DATEcur $DATEcur $BLACKL $ERROBS $CONFIG_FILE $CONFIG_PAR < $Bspat > $LOGDIR/Bspat_PREC1hRT_$DATEcur_$DATEcur.log 2>&1
    SECcur=$(( SECcur+3600 ))
  done
#--------------------------
# Exit 
#--------------------------
  echo "seNorge2_PREC1hRT.sh "`date +%Y-%m-%d" "%H:%M`" > success!"
  exit 0
