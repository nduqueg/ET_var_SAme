# Drivers of Evapotranspiration in tropical south America

This repository contains the set of scripts to reproduce the results of the paper:

Duque-Gardeazabal, Nicolas., Friedman, Andrew R., Brönnimann, Stefan. An Atlantic influence on evaporation in the Orinoco and Amazon basins (2024).

## prerequisites

To run the scripts: 
- Climate Data Operators (CDO) tool is necessary for preprocessing the data, specifically its aggregation from monthly to seasonal timescale
- they must be run in the sequence they were listed, starting by the preprocess directory and then some scripts in the main data subdirectories
- most necessary libraries are Raster, ncdf4, hydroTSM, metR, reshape, ggplot2, RColorBrewer, magrittr, gridExtra, ggrepel

## main data organization

As the research used several data sources, to reproduce the results is imperative to allocate the downloaded data to the appropriate subdirectory. The latter is:
- The [ERSST](https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/netcdf/) is downloaded automatically with an script and stored in the preprocessing directory
- CO_2 concentrations from Mauna Loa are already included in the repository
- Precipitation from [MSWEP](http://www.gloh2o.org/mswep/) must be downloaded and crop to the study region, and stored in the subdirectory "02_Ppt" inside the data directory. The same for [CHIRPS](https://data.chc.ucsb.edu/products/CHIRPS-2.0)
- ERA5 and ERA5-Land variables must be placed in subdiretories from 03 to 09, and the datasets are available from Copernicus Climate Data Store web portal https://cds.climate.copernicus.eu
- ESA CCI SM dataset must be placed in subdirectory "03_SM" and is available at: https://catalogue.ceda.ac.uk/uuid/ff890589c21f4033803aa550f52c980c
- GLEAM must be placed in subdirectory "04_Evap" and is available at: https://www.gleam.eu/
- EUMETSAT-CLARA must be placed in subdirectory "05_Rad" and is available at: https://wui.cmsaf.eu/safira/action/viewProduktDetails?fid=40&eid=22277_22492
