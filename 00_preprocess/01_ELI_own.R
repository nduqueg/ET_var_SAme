rm(list=ls())
cat("\014")

library(raster)
library(ncdf4)
library(hydroTSM)
library(ggplot2)
library(ggrepel)

# Dates ----
Dm.h <- seq(as.Date("1850-01-01"),as.Date("2022-11-30"),by="month")
Dm.e <- seq(as.Date("1854-01-01"),as.Date("2020-12-31"),by="month")
Dm.ana <- seq(as.Date("1885-01-01"),as.Date("2020-12-30"),by="month")

# load data ####

H.sst <- brick("./00_preprocess/HadSST.4.0.1.0_median.nc") # not with the anomalies, it's with the actual value
E.sst <- brick("./00_preprocess/1_ersst_v5_1854_2020.nc")

H.sst.a <- H.sst[[ which(Dm.h >= Dm.ana[1] & Dm.h <= Dm.ana[length(Dm.ana)] )]]
E.sst.a <- E.sst[[ which(Dm.e >= Dm.ana[1] & Dm.e <= Dm.ana[length(Dm.ana)] )]]

lat_Weight.Mean <- function(R, Dm.ana=Dm.ana){
  aux <- rasterToPoints(R)
  R.weight <- apply(
    aux[,3:length(Dm.ana)]*sqrt(cos(3.14159*aux[,2]/180)) 
    ,2, mean, na.rm=T )
  return( zoo(R.weight, order.by = Dm.ana))
} # Special function for calculated the mean with latitude weighting, also for the two datasets

## Convection Treshold ----
# H.sst.Tres <- lat_Weight.Mean( crop(H.sst.a,extent(c(-180,180,-5,5))), Dm.ana=Dm.ana)
# E.sst.Tres <- lat_Weight.Mean(crop(E.sst.a,extent(c(0,360,-5,5))), Dm.ana=Dm.ana)
H.sst.Tres <- cellStats( crop(H.sst.a,extent(c(-180,180,-5,5))), stat="mean", na.rm=T)
E.sst.Tres <- cellStats(crop(E.sst.a,extent(c(0,360,-5,5))), stat="mean",na.rm=T)
SST.Tres <- cbind(H.sst.Tres,E.sst.Tres)

### masked values ---
#Coord.H <- rasterToPoints(crop(H.sst.a[[1]],extent(c(-180,180,-5,5))))[,-3]; Coord.H[,1] <- Coord.H[,1] +360
Coord.E <- rasterToPoints(crop(E.sst.a[[1]],extent(c(0,360,-5,5))))[,-3]

#H.sst.c <- rasterToPoints(crop(H.sst.a,extent(c(-180,180,-5,5))))[,-c(1,2)]
E.sst.c <- rasterToPoints(crop(E.sst.a,extent(c(0,360,-5,5))))[,-c(1,2)]

#H.above <- matrix(FALSE, ncol=length(Dm.ana), nrow= nrow(H.sst.c))
E.above <- matrix(FALSE, ncol=length(Dm.ana), nrow= nrow(E.sst.c))
for (i in 1:length(Dm.ana)){
  #H.above[,i] <- H.sst.c[,i] >= as.numeric(H.sst.Tres[i])
  E.above[,i] <- E.sst.c[,i] >= as.numeric(E.sst.Tres[i])
}
#H.above[is.na(H.above)] <- FALSE

#### ELI index ----
#ELI.H <- zoo(NA, order.by = Dm.ana)
ELI.E <- zoo(NA, order.by = Dm.ana)
for (i in 1:length(Dm.ana)){
  aux <- Coord.E[ E.above[,i] ,1]
  aux <- aux[aux >=115 & aux <=360-70]
  ELI.E[i] <- mean(aux)
  #ELI.H[i] <- ELI.H[i]
}

ELI.E <- as.data.frame(ELI.E)
write.csv(ELI.E,"./01_data/01_SSTindices/00A_ELI_SST.csv")
