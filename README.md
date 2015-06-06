Spatial Interpolation of in-situ observations
==============================================

The objective is to create Temperature and Precipitation high-resolution (1Km) gridded datasets for the Norwegian mainland using observation from the KDVH (Climate Database).
The spatial interpolation method is based on Bayesian statical interpolation and the resulting gridded datasets constitutes the seNorge v 2.0 beta dataset.

Notes:
1. Input data are retrieved from MET Norway Climate Database.
2. it is possible for observations affected by gross measurement errors to enter the spatial interpolation procedure.

seNorge v 2.0 Outputs:
----------------------
* TEMP1d: Mean temperature (yesterday at 06 - today at 06 UTC); unit ºC
* TEMP1h: Air temperature at time of observation (hourly sampling rate); unit ºC
* TEMP1d24h: Mean temperature (yesterday at 06 - today at 06 UTC) using TEMP1h as background; unit ºC
* PREC1d: Daily total of precipitation (precipitation day definition: yesterday at 06 UTC - today at 06 UTC); unit mm
* PRECcor1d: Daily total of precipitation with exposure correction (precipitation day definition: yesterday at 06 UTC - today at 06 UTC); unit mm
* REC1hRT: Amount of precipitation in 1 hour (for real-time applications); unit mm
* PREC1h: Amount of precipitation in 1 hour (consistent with the daily total of precipitation); unit mm
* PREC3hRT: Amount of precipitation in 3 hour (for real-time applications); 00-03, 03-06, 06-09, …, 21-00; unit mm
* PREC3h: Amount of precipitation in 3 hour (consistent with the daily total of precipitation); unit mm

Installation:
-------------
1. following libraries must be installed on your system:
  * gdal-bin
  * proj-bin
  * libgdal1-dev
  * libproj-dev
  * netcdf-bin
  * libnetcdf-dev
  * netcdf-bin

2. get the following packages from r-cran repository:
  * sp_1.0-9.tar.gz
  * raster_2.1-25.tar.gz
  * rgdal_0.8-9.tar.gz
  * tripack_1.3-6.tar.gz
  * ncdf_1.6.6.tar.gz
  * igraph_0.7.1.tar.gz

3. create directories (DIR1 should be any directory path on your system):

.. code-block:: bash

mkdir DIR1/log
mkdir DIR1/projects
mkdir DIR1/scripts
mkdir DIR1/seNorge2
mkdir DIR1/seNorge2_addInfo
mkdir DIR1/seNorge2_scratch

4. clone git hub repository (let's say you've done it under directory DIR2)

5. Edit/Create config file (in a standard installation should be DIR2/seNorge2/etc/config_list.r):

.. code-block:: bash

list( pname="your name here (let's call it NAMECONF)",
  opt = alist(
  main.path=DIR2,
  main.path.output=DIR1,
  testmode=FALSE)
),

testmode=FALSE get the data from MET Norway's Climate Database


Running the programs (examples):
--------------------------------
1. PREC1hRT from 2012.01.01.04 (yyyy.mm.dd.hh, UTC) to 2012.01.01.06:

.. code-block:: bash

DIR2/seNorge2/script/seNorge2_PREC1hRT.sh -s 2012.01.01.04 -e 2012.01.01.06 -c DIR2/seNorge2/etc/config_list.r -p NAMECONF -l DIR1/log

2. PREC3hRT from today at O1 UTC to O3 UTC:

.. code-block:: bash

DIR2/seNorge2/script/seNorge2_PREC3hRT_just3h.sh -s `date +\%Y.\%m.\%d`.01 -e `date +\%Y.\%m.\%d`.03 -c DIR2/seNorge2/etc/config_list.r -p NAMECONF -l DIR1/log

3. PREC1d for today:

.. code-block:: bash
DIR2/seNorge2/script/seNorge2_PREC1d.sh -s `date +\%Y.\%m.\%d` -e `date +\%Y.\%m.\%d` -c DIR2/seNorge2/etc/config_list.r -p NAMECONF -l DIR1/log
