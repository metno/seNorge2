seNorge v 2.0
2015.02.26 - Cristian Lussana (cristian.lussana@met.no)

The objective is to create Temperature and Precipitation high-resolution (1Km) gridded datasets for the Norwegian mainland using observation from the KDVH (Climate Database).
The spatial interpolation method is based on Bayesian statical interpolation and the resulting gridded datasets constitutes the seNorge v 2.0 beta dataset.
Any feedback is strongly appreciated. Please, report any idea/comment/bug to cristian.lussana@met.no
!! Note that it is possible for observations affected by gross measurement errors to enter the spatial interpolation procedure.!!

Input variables from KDVH (through Ulric): TAMRR, RR, RR1, TA

seNorge v 2.0 Outputs:\\
    TEMP1d: Mean temperature (yesterday at 06 - today at 06 UTC); unit ºC
        input: TAMRR
        time range: 1957.01.01 to present
    TEMP1h: Air temperature at time of observation (hourly sampling rate); unit ºC
        input: TA
        time range: 2010.01.01 to present
    PREC1d: Daily total of precipitation (precipitation day definition: yesterday at 06 UTC - today at 06 UTC); unit mm
        input: RR
        time range: 1957.01.01 to present
    PRECcor1d: Daily total of precipitation with exposure correction (precipitation day definition: yesterday at 06 UTC - today at 06 UTC); unit mm
        input: RR, TEMP1d
        time range: 1957.01.01 to present
    PREC1hRT: Amount of precipitation in 1 hour (for real-time applications); unit mm
        input: RR1
        time range: 2010.01.01 to present
    PREC1h: Amount of precipitation in 1 hour (consistent with the daily total of precipitation); unit mm
        input: RR1, PREC1d, PREC1hRT
        time range: 2010.01.01 to present
    PREC3hRT: Amount of precipitation in 3 hour (for real-time applications); 00-03, 03-06, 06-09, …, 21-00; unit mm
        input: RR1
        time range: 2010.01.01 to present
    PREC3h: Amount of precipitation in 3 hour (consistent with the daily total of precipitation); unit mm
        input: RR1, PREC1d, PREC3hRT
        time range: 2010.01.01 to present
General info:
    Geographical information used: see /lustre/mnt/cristianl/seNorge2/geoinfo
    Horizontal resolution: 1 Km
    PROJ4 type description of the Coordinate Reference System: ”+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0”
    Timezone: UTC with 0-23 hour
    output filename format: netCDF
