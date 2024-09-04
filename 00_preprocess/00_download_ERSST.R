rm(list = ls())
cat("\014")

library(raster)
library(ncdf4)

ersst_url <- "https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/netcdf/"
dist.dir <- "./individual/"

yr_start <-1854
yr_end <-2020

for (yr in seq(yr_start,yr_end,1)){
  print(paste("downloading year",yr))
  for(mo in c(paste0("0",seq(1,9,1)),"10","11","12")){
    download.file(paste0(ersst_url,"ersst.v5.",yr,mo,".nc"), paste0(dist.dir,"ersst.v5.",yr,mo,".nc"))
  }
}

print("stacking")
dir.create(dist.dir)
a<-list.files(path=dist.dir, pattern = ".nc")

SST <- stack(paste0(dist.dir,a),varname="sst")
SSTa <- stack(paste0(dist.dir,a),varname="ssta")

writeRaster(SST,filename = "1_ersst_v5_1854_2020.nc",format="CDF")
writeRaster(SSTa,filename = "Anom_ersst_v5_1854_2020.nc",format="CDF")

file.remove(dist.dir)