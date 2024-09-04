rm(list=ls())
cat("\014")

library(ggplot2)
library(RColorBrewer)
library(raster)

Dates  <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month")
date.ini <- format(Dates[1],format="%Y-%m-%d"); date.end <- format(Dates[length(Dates)],format="%Y-%m-%d")

# identifying trends ERA5-Land ----
ifile <- "../../../01_DataSets/04_Evap/02_ERA5L/ERA5L_TEvap_1950_2020.nc"
ofile <- "01_ERA5L_TEvap_1980_2020.nc"
cdo.cmd <- paste0("cdo seldate,",date.ini,",",date.end," ",ifile," ",ofile)
system(cdo.cmd)

ifile <- ofile
ofile.b <- "02_trend_b_ERA5L_TEvap_1980-2020.nc"; ofile.m <- "02_trend_m_ERA5L_TEvap_1980-2020.nc"
cdo.cmd <- paste0("cdo trend ",ifile," ", ofile.b," ", ofile.m)
system(cdo.cmd)

b <- brick(ofile.b); m <- brick(ofile.m)
b.t <- rasterToPoints(b); m.t <- rasterToPoints(m)

E <- cbind(b.t, m.t[,3]); colnames(E) <- c("lon","lat","b","m"); E <- as.data.frame(E)
E[,"m"] <- E[,"m"] * 12 * 10 # transform to mm/decade


# ploting
library(gridExtra)
Basins <- shapefile("../../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")
SA <- shapefile("../../../01_DataSets/South_America/South_America.shp")

at.m <- round(seq( min(E[,"b"]), max(E[,"b"]), length.out = 12))
p1 <- ggplot(E[,c("lon","lat","b")], aes(lon,lat)) + 
  geom_raster(aes(fill=b)) + scale_fill_stepsn(colours=brewer.pal(11,"RdYlGn"), breaks=at.m, guide=guide_colorsteps(barwidth=unit(10,"cm")),
                                                 limits=c(min(at.m),max(at.m)), # values=scales::rescale(at.v,from=range(at.m))
                                                 name="Intercept [mm]")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  coord_sf(xlim=c(-82,-35),ylim=c(-20,15))+labs(caption = "ERA5L")+
  theme_bw()+theme(legend.position = "bottom")

range <- max(abs(min(E[,"m"])), max(E[,"m"]))
at.m <- round(seq( -range, range, length.out = 12))
p2 <- ggplot(E[,c("lon","lat","m")], aes(lon,lat)) + 
  geom_raster(aes(fill=m)) + scale_fill_stepsn(colours=brewer.pal(11,"PiYG"), breaks=at.m, guide=guide_colorsteps(barwidth=unit(10,"cm")),
                                               limits=c(min(at.m),max(at.m)), # values=scales::rescale(at.v,from=range(at.m))
                                               name="slope/change [mm/decade]")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  coord_sf(xlim=c(-82,-35),ylim=c(-20,15))+labs(caption = "ERA5L")+
  theme_bw()+theme(legend.position = "bottom")

grid.arrange(p1,p2,ncol=2)

## detrend ERA5-Land & seasonal acummulation  ---- 
ifile <- "01_ERA5L_TEvap_1980_2020.nc"
ofile <- "03_detrended_ERA5L_TEvap_1980_2020.nc"
cdo.cmd <- paste0("cdo detrend"," ",ifile," ",ofile)
system(cdo.cmd)
file.remove(ifile)

ifile <- ofile; ofile <- "03_detrended_ET_ERA5L_seasonal_1980_2020.nc"
cdo.cmd <- paste0("cdo seassum ",ifile," ",ofile)
system(cdo.cmd)

# identifying trends GLEAM ----

Dates  <- seq(as.Date("1979-12-01"),as.Date("2020-12-01"), by="month")
date.ini <- format(Dates[1],format="%Y-%m-%d"); date.end <- format(Dates[length(Dates)],format="%Y-%m-%d")

ifile <- "../../../01_DataSets/04_Evap/01_GLEAM/GLEAM_E_1980-2020.nc"
ofile <- "01_GLEAM_TEvap_1980_2020.nc"
cdo.cmd <- paste0("cdo seldate,",date.ini,",",date.end," ",ifile," ",ofile)
system(cdo.cmd)

ifile <- ofile
ofile.b <- "02_trend_b_GLEAM_TEvap_1980-2020.nc"; ofile.m <- "02_trend_m_GLEAM_TEvap_1980-2020.nc"
cdo.cmd <- paste0("cdo trend ",ifile," ", ofile.b," ", ofile.m)
system(cdo.cmd)

b <- brick(ofile.b); m <- brick(ofile.m)
b.t <- rasterToPoints(b); m.t <- rasterToPoints(m)

E <- cbind(b.t, m.t[,3]); colnames(E) <- c("lon","lat","b","m"); E <- as.data.frame(E)
E[,"m"] <- E[,"m"] * 12 * 10 # transform to mm/decade


# ploting
library(gridExtra)
Basins <- shapefile("../../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")
SA <- shapefile("../../../01_DataSets/South_America/South_America.shp")

at.m <- round(seq( min(E[,"b"]), max(E[,"b"]), length.out = 12))
p1 <- ggplot(E[,c("lon","lat","b")], aes(lon,lat)) + 
  geom_raster(aes(fill=b)) + scale_fill_stepsn(colours=brewer.pal(11,"RdYlGn"), breaks=at.m, guide=guide_colorsteps(barwidth=unit(10,"cm")),
                                               limits=c(min(at.m),max(at.m)), # values=scales::rescale(at.v,from=range(at.m))
                                               name="Intercept [mm]")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  coord_sf(xlim=c(-82,-35),ylim=c(-20,15))+labs(caption = "GLEAM")+
  theme_bw()+theme(legend.position = "bottom")

range <- max(abs(min(E[,"m"])), max(E[,"m"]))
at.m <- round(seq( -range, range, length.out = 12))
p2 <- ggplot(E[,c("lon","lat","m")], aes(lon,lat)) + 
  geom_raster(aes(fill=m)) + scale_fill_stepsn(colours=brewer.pal(11,"PiYG"), breaks=at.m, guide=guide_colorsteps(barwidth=unit(10,"cm")),
                                               limits=c(min(at.m),max(at.m)), # values=scales::rescale(at.v,from=range(at.m))
                                               name="slope/change [mm/decade]")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  coord_sf(xlim=c(-82,-35),ylim=c(-20,15))+labs(caption = "GLEAM")+
  theme_bw()+theme(legend.position = "bottom")

grid.arrange(p1,p2,ncol=2)

## detrend GLEAM  & seasonal acummulation ---- 
ifile <- "01_GLEAM_TEvap_1980_2020.nc"
ofile <- "03_detrended_GLEAM_TEvap_1980_2020.nc"
cdo.cmd <- paste0("cdo detrend"," ",ifile," ",ofile)
system(cdo.cmd)
file.remove(ifile)

ifile <- ofile; ofile <- "03_detrended_ET_GLEAM_seasonal_1980_2020.nc"
cdo.cmd <- paste0("cdo seassum ",ifile," ",ofile)
system(cdo.cmd)