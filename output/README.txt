These GR-Placefiles are produced by upscaling a 1-hour HRRR forecast to a 9-km grid-spacing 
(computational resources and data storage considerations) and running each grid point through
SHARPpy's parcel functions. While these fields will not exactly match those on the SPC 
mesoanalysis page, they should be close, especially for parameters that rely on effective 
inflow bases and tops (ESRH, EBWD, eSTP, etc.)

These files should be available by :01 after the hour, and are set to auto-update every minute
after they're loaded into GR. 

In the event of a NOMADS outage, this will fall back to searching the back-up FTPPRD site and 
finally a Google Cloud repository for realtime HRRR data. Data availablity is not gauranteed,
however, so be sure to check the timestamps in the GR Placefile Window!
