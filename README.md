seNorge v 2.0 - 2015.02.26 - Cristian Lussana (cristian.lussana@met.no)

The objective is to create Temperature and Precipitation high-resolution (1Km) gridded datasets for the Norwegian mainland using observation from the KDVH (Climate Database).
The spatial interpolation method is based on Bayesian statical interpolation and the resulting gridded datasets constitutes the seNorge v 2.0 beta dataset.
Any feedback is strongly appreciated.
Note that it is possible for observations affected by gross measurement errors to enter the spatial interpolation procedure.

seNorge v 2.0 Outputs:
- TEMP1d: Mean temperature (yesterday at 06 - today at 06 UTC); unit ºC
- TEMP1h: Air temperature at time of observation (hourly sampling rate); unit ºC
- PREC1d: Daily total of precipitation (precipitation day definition: yesterday at 06 UTC - today at 06 UTC); unit mm
- PRECcor1d: Daily total of precipitation with exposure correction (precipitation day definition: yesterday at 06 UTC - today at 06 UTC); unit mm
- REC1hRT: Amount of precipitation in 1 hour (for real-time applications); unit mm
- PREC1h: Amount of precipitation in 1 hour (consistent with the daily total of precipitation); unit mm
- PREC3hRT: Amount of precipitation in 3 hour (for real-time applications); 00-03, 03-06, 06-09, …, 21-00; unit mm
- PREC3h: Amount of precipitation in 3 hour (consistent with the daily total of precipitation); unit mm
