rm(list=ls())
cat("\014")

library(reshape)
library(raster)
library(ncdf4)
library(RColorBrewer)
library(hydroTSM)
library(ggplot2)
library(metR)

# DAtes ----
dates.m <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"),by="month")
Dm.e5 <- seq(as.Date("1950-01-01"),as.Date("2020-12-30"), by="month"); Dy.e5 <- format(Dm.e5,format="%Y")
month.Multi <- format(dates.m,format="%m")

years <- factor(format(dates.m,format="%Y")); Years <- unique(years); Years <- Years[-1]
Season.y <- paste0(years[-1],"-",time2season(dates.m[-length(dates.m)])); Season.y <- c(Season.y,Season.y[length(Season.y)])
Year.s <- unique(Season.y)

years.m.e5 <- format(Dm.e5, format="%Y"); Years.e5 <- unique(format(Dm.e5, format="%Y")); Years.e5 <- Years.e5[-1]
Season.y.e5 <- paste0(years.m.e5[-1], "-",time2season(Dm.e5)[-length(Dm.e5)])
Year.s.e5 <- unique(Season.y.e5)

seasons <- c("DJF","MAM","JJA","SON")

Dm.sst <- seq(as.Date("1885-01-01"),as.Date("2020-12-30"),by="month")
## loading files ####
dir.mf <- "../../../01_DataSets/10_VIMF/"

# cdo seasmean ERA5_VIMF_eastNorth_1950-2020.nc ERA5_VIMF_eastNorth_seasonal.nc
Mt.e.s <- brick(paste0(dir.mf,"ERA5_VIMF_eastNorth_seasonal.nc"),varname="p71.162"); Mt.e.s <- Mt.e.s[[match(Year.s,Year.s.e5)]]
Mt.n.s <- brick(paste0(dir.mf,"ERA5_VIMF_eastNorth_seasonal.nc"),varname="p72.162"); Mt.n.s <- Mt.n.s[[match(Year.s,Year.s.e5)]]

basins <- shapefile("../../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")

AMM <- read.zoo(read.csv("../../../01_DataSets/01_SST/04_Cal_Indices/00A_AMM_deltaSST.csv")[,-2], index.column = 1); AMM <- AMM[dates.m]

## seasonal preprocessing ----
#------------------------------------------------------------------ VIMF seasonal separation
{
Mt.e.s.t <- rasterToPoints(Mt.e.s)[,-c(1,2)]
Mt.n.s.t <- rasterToPoints(Mt.n.s)[,-c(1,2)]

Coord <- as.data.frame(rasterToPoints(Mt.e.s[[1]])[,c(1,2)])

Mt.E.seasons <- list(); Mt.N.seasons <- list()
for ( i in seasons){
  Mt.E.seasons[[i]] <- as.data.frame(matrix(NA, ncol=nrow(Mt.e.s.t), nrow=length(Years)))
  Mt.N.seasons[[i]] <- as.data.frame(matrix(NA, ncol=nrow(Mt.n.s.t), nrow=length(Years)))
  }

for(i in seasons){
  Mt.E.seasons[[i]][1:length(Years),] <- t(Mt.e.s.t[, which(substr(Year.s,6,8) == i)])
  Mt.N.seasons[[i]][1:length(Years),] <- t(Mt.n.s.t[, which(substr(Year.s,6,8) == i)])
}
}#------------------ VIMF separation

#------------------------------------------------------------------ seasonal mean
{
  Mean.MtE.seasons <- list(); Mean.MtN.seasons <- list()
  for ( i in seasons){
    Mean.MtE.seasons[[i ]] <- apply(Mt.E.seasons[[i ]],2, mean, na.rm=T)
    Mean.MtN.seasons[[i ]] <- apply(Mt.N.seasons[[i ]],2, mean, na.rm=T)
  }
  Mean.MtE.seasons <- simplify2array(Mean.MtE.seasons)
  Mean.MtN.seasons <- simplify2array(Mean.MtN.seasons)
  
  Mean.Angle.seasons <- Mean.MtE.seasons
  Mean.Mag.seasons <- Mean.MtE.seasons
  
  # for ( i in seasons){
  #   Mean.Angle.seasons[,i] <- atan2(dlat(Mean.MtN.seasons[,i]), dlon(Mean.MtE.seasons[,i],Coord[,"y"]))*180/pi
  #   Mean.Mag.seasons[,i] <- Mag(Mean.MtN.seasons[,i], Mean.MtE.seasons[,i])
  # }
}#----------- seasonal mean

#------------------------------------------------------------------ SST accummulation to seasons and separation
{
  # AMM.s <- aggregate(AMM[-length(AMM)], by= list(Season.y), FUN= mean)
  AMM.seasons <- read.csv("../AMM_std_1980-2020.csv")[,-1]
}#--------- SST accummulation to seasons and separation

#### identification of periods SST higher and lower than 1*SD --------
#--------------------------------------------------------------------------- initialization of datasets
{
  Mt.E.anom <- list(); Mt.N.anom <- list(); Mt.Mag.anom <- list(); Mt.Angle.anom <- list()
  matrix(NA,nrow=nrow(Coord),ncol=4) -> Mt.Mag.anom[["Pos"]] -> Mt.Mag.anom[["Neg"]] -> Mt.Angle.anom [["Pos"]] -> Mt.Angle.anom [["Neg"]] ->
    Mt.E.anom[["Pos"]] -> Mt.E.anom[["Neg"]] -> Mt.N.anom[["Pos"]] -> Mt.N.anom[["Neg"]]
  
  Mt.E.anom <- lapply(Mt.E.anom, 'colnames<-',seasons); Mt.N.anom <- lapply(Mt.N.anom, 'colnames<-',seasons)
  Mt.Mag.anom <- lapply(Mt.Mag.anom, 'colnames<-',seasons); Mt.Angle.anom <- lapply(Mt.Angle.anom, 'colnames<-',seasons)
  
  Mt.E.test <- list(); Mt.N.test <- list(); Mt.Mag.test <- list()
  matrix(NA,nrow=nrow(Coord),ncol=4) -> Mt.Mag.test[["Pos"]] -> Mt.Mag.test[["Neg"]] ->
    Mt.E.test[["Pos"]] -> Mt.E.test[["Neg"]] -> Mt.N.test[["Pos"]] -> Mt.N.test[["Neg"]]
  
  Mt.E.test <- lapply(Mt.E.test, 'colnames<-',seasons); Mt.N.test <- lapply(Mt.N.test, 'colnames<-',seasons)
  Mt.Mag.test <- lapply(Mt.Mag.test, 'colnames<-',seasons)
  }

AMM.bool <- list()
#---------------------------------------------------------------------- identification of periods with SST higher than 1*SD
{
  AMM.bool[["Pos"]] <- data.frame(matrix(FALSE, ncol=4, nrow= length(Years))); colnames(AMM.bool[["Pos"]]) <- seasons; rownames(AMM.bool[["Pos"]]) <- Years

AMM.bool[["Pos"]][which(AMM.seasons >=1, arr.ind = T)] <- TRUE}

#---------------------------------------------------------------------- identification of periods with SST Lower than -1*SD
{
  AMM.bool[["Neg"]] <- data.frame(matrix(FALSE, ncol=4, nrow= length(Years))); colnames(AMM.bool[["Neg"]]) <- seasons; rownames(AMM.bool[["Neg"]]) <- Years
  
  AMM.bool[["Neg"]][which(AMM.seasons <=-1, arr.ind = T)] <- TRUE}

a <- t(sapply(AMM.bool,function(x){apply(x,2,sum)}))
write.csv(a,"AMM_events_1980-2020.csv")


Test.comp.f <- function(Phase, neutral){
  
  Result <- rep(NA , nrow(Phase))
  # retriving cells without value and with constant series in which case it cannot be evaluated
  Idx <- which(apply(Phase,2,function(x) sum(!is.na(x)) > 1))
  a <- match(which(apply(Phase,2,function(x) sd(x)==0)), Idx); a <- a[!is.na(a)]; if(length(a!=0)>0) Idx <- Idx[- a]
  b <- match(which(apply(neutral,2,function(x) sd(x)==0)), Idx); b <- b[!is.na(b)]; if(length(b!=0)) Idx <- Idx[-b]
  # needed to use a loop because those are two datasets that cannot change with an apply (unless we built and array and that's more complicated)
  for(j in Idx) Result[j] <- t.test( Phase[,j], neutral[,j ])$p.value
  Result[a] <- -999; if(length(b!=0)) Result[b] <- -999
  return(Result)
}

#------------------------------------------------------------------------ identification MT anomalies in those periods
for ( i in seasons){
  print(paste("Anomalies - ",i))
  E.Pos <- Mt.E.seasons[[i]][AMM.bool[["Pos"]][,i],]; N.Pos <- Mt.N.seasons[[i]][AMM.bool[["Pos"]][,i],]
  Mt.E.anom[["Pos"]][,i] <- apply( E.Pos  ,2,mean,na.rm=T) - Mean.MtE.seasons[,i]
  Mt.N.anom[["Pos"]][,i] <- apply( N.Pos  ,2,mean,na.rm=T) - Mean.MtN.seasons[,i]
  
  Mt.Angle.anom[["Pos"]][,i] <- atan2(dlat(Mt.N.anom[["Pos"]][,i]), dlon(Mt.E.anom[["Pos"]][,i], Coord[,"y"]))*180/pi
  Mt.Mag.anom[["Pos"]][,i] <- Mag(Mt.N.anom[["Pos"]][,i],Mt.E.anom[["Pos"]][,i])
  
  E.Neg <- Mt.E.seasons[[i]][AMM.bool[["Neg"]][,i],]; N.Neg <- Mt.N.seasons[[i]][AMM.bool[["Neg"]][,i],]
  Mt.E.anom[["Neg"]][,i] <- apply(  E.Neg ,2,mean) - Mean.MtE.seasons[,i]
  Mt.N.anom[["Neg"]][,i] <- apply( N.Neg  ,2,mean) - Mean.MtN.seasons[,i]
  
  Mt.Angle.anom[["Neg"]][,i] <- atan2(dlat(Mt.N.anom[["Neg"]][,i]), dlon(Mt.E.anom[["Neg"]][,i], Coord[,"y"]))*180/pi
  Mt.Mag.anom[["Neg"]][,i] <- Mag(Mt.N.anom[["Neg"]][,i],Mt.E.anom[["Neg"]][,i])
  
  print(paste("Statistical significance -", i))
  neutral <- !(AMM.bool[["Pos"]][,i] | AMM.bool[["Neg"]][,i])
  E.neutral <- subset(Mt.E.seasons[[i]],
                      neutral); N.neutral <- subset(Mt.N.seasons[[i]], neutral)
  
  Mt.E.test[["Pos"]][,i] <- Test.comp.f(E.Pos, E.neutral); Mt.N.test[["Pos"]][,i] <- Test.comp.f(N.Pos, N.neutral)
  Mt.E.test[["Neg"]][,i] <- Test.comp.f(E.Neg, E.neutral); Mt.N.test[["Neg"]][,i] <- Test.comp.f(N.Neg, N.neutral)
}
colnames(Coord) <- c("lon", "lat")
Mt.E.anom <- lapply(Mt.E.anom, cbind, Coord)
Mt.N.anom <- lapply(Mt.N.anom, cbind, Coord)
Mt.E.test <- lapply(Mt.E.test, cbind, Coord)
Mt.N.test <- lapply(Mt.N.test, cbind, Coord)

data.anom.uv <- list(Mt.E.anom= Mt.E.anom, Mt.N.anom= Mt.N.anom)
save(data.anom.uv,"data.anom.uv",file="01_AnomComp_AMM_VIMF_uv.RData")
save(Mt.E.test, Mt.N.test, list=c("Mt.E.test","Mt.N.test"), file="02_VIMF_AMM_Comp_Ttest_uv.RData")
##### Data transformation for plotting ------
#------------------------------------------------------------------------ Seasonal mean transformation
{
  Mean.Angle.seasons <- cbind(Coord, Mean.Angle.seasons)
  Mean.Mag.seasons <- cbind(Coord, Mean.Mag.seasons)
  
  data.mean <- melt(Mean.Angle.seasons, id=c("lon","lat"))
  data.mean <- cbind(data.mean, melt(Mean.Mag.seasons, id=c("lon","lat"))[4])
  colnames(data.mean) <- c("lon", "lat","Season","Angle","Magnitud")
}

#------------------------------------------------------------------------ Anomalies in positive AMM
data.anom <- list()
{
  Mt.Angle.anom[["Pos"]] <- cbind(Coord,Mt.Angle.anom[["Pos"]])
  Mt.Mag.anom[["Pos"]] <- cbind(Coord, Mt.Mag.anom[["Pos"]])

  data.anom[["Pos"]] <- melt(Mt.Angle.anom[["Pos"]], id=c("lon","lat"))
  data.anom[["Pos"]] <- cbind(data.anom[["Pos"]], melt(Mt.Mag.anom[["Pos"]], id=c("lon","lat"))[4])
  colnames(data.anom[["Pos"]]) <- c("lon", "lat","Season","Angle","Anomaly")
  }

#------------------------------------------------------------------------ Anomalies in Negative AMM
{
  Mt.Angle.anom[["Neg"]] <- cbind(Coord,Mt.Angle.anom[["Neg"]])
  Mt.Mag.anom[["Neg"]] <- cbind(Coord, Mt.Mag.anom[["Neg"]])
  
  data.anom[["Neg"]] <- melt(Mt.Angle.anom[["Neg"]], id=c("lon","lat"))
  data.anom[["Neg"]] <- cbind(data.anom[["Neg"]], melt(Mt.Mag.anom[["Neg"]], id=c("lon","lat"))[4])
  colnames(data.anom[["Neg"]]) <- c("lon", "lat","Season","Angle","Anomaly")
}

save(data.anom,"data.anom",file="01_AnomComp_AMM_VIMF.RData")

#### plotting ----
paleta <- c(brewer.pal(9,"YlGnBu"),"#81007F"); paleta[1:3] <- c("#9e0142","#f46d43","#fad366")
# ---------------------------------------------------------------------------------------------------------- filtering arrows for better visualization
{
  sele <- list()
  # sele[["lon"]] <- seq(min(data$lon),max(data$lon),by=0.25)
  # sele[["lat"]] <- seq(min(data$lat),max(data$lat),by=0.25)
  sele[["lon"]] <- seq(-83,0,by=0.5)
  sele[["lat"]] <- seq(-21,21,by=0.5)
  data.mean <- data.mean[ which(data.mean$lon %in% sele$lon) ,]; data.mean <- data.mean[ which(data.mean$lat %in% sele$lat) ,]
  data.anom[["Neg"]] <- data.anom[["Neg"]][ which(data.anom[["Neg"]]$lon %in% sele$lon) ,]; data.anom[["Neg"]] <- data.anom[["Neg"]][ which(data.anom[["Neg"]]$lat %in% sele$lat) ,]
  data.anom[["Pos"]] <- data.anom[["Pos"]][ which(data.anom[["Pos"]]$lon %in% sele$lon) ,]; data.anom[["Pos"]] <- data.anom[["Pos"]][ which(data.anom[["Pos"]]$lat %in% sele$lat) ,]
}

# ---------------------------------------------------------------------------------------------------------- plotting both fields (mean & anomalies)
plotting.anom <- function (data.mean, data.anom, Spe.season, title.p, xlim=c(-82,-40), ylim=c(-20,15)){
  # data.anom must be the dataFrame, not the whole list
  
  data.mean.p <- data.mean[data.mean$Season==Spe.season,]
  data.anom.p <- data.anom[data.anom$Season==Spe.season,]
  p1 <- ggplot(data.mean.p, aes(lon,lat))+#facet_wrap(~Season,ncol=2)+
    # scale_y_continuous(position="right")+    geom_contour(aes(z = gh.z)) +
    geom_vector(aes(angle=Angle, mag=Magnitud),
                col="blue", pivot=0.5, skip=1, alpha=0.4)+
    geom_arrow(data= data.anom.p, aes(x=lon,y=lat, angle=Angle, mag=Anomaly),
                col="red", pivot=0.5, skip=1, show.legend = TRUE)+
    # scale_color_gradientn(colours = paleta)+ #brewer.pal(11,"Spectral")
    geom_polygon(data=basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
    coord_fixed(xlim=xlim,ylim=ylim)+
    labs(x="Long",y="Lat",title= title.p)+
    scale_mag()+
    theme_bw()
  
  print(p1)
}

{
#   plotting.anom(data.mean, data.anom = data.anom[["Pos"]], Spe.season = "DJF", title.p = "Anomaly Composite VIMF - DJF (1980-2020, AMM >= 0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Pos"]], Spe.season = "MAM", title.p = "Anomaly Composite VIMF - MAM (1980-2020, AMM >= 0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Pos"]], Spe.season = "JJA", title.p = "Anomaly Composite VIMF - JJA (1980-2020, AMM >= 0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Pos"]], Spe.season = "SON", title.p = "Anomaly Composite VIMF - SON (1980-2020, AMM >= 0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Neg"]], Spe.season = "DJF", title.p = "Anomaly Composite VIMF - DJF (1980-2020, AMM <= -0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Neg"]], Spe.season = "MAM", title.p = "Anomaly Composite VIMF - MAM (1980-2020, AMM <= -0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Neg"]], Spe.season = "JJA", title.p = "Anomaly Composite VIMF - JJA (1980-2020, AMM <= -0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Neg"]], Spe.season = "SON", title.p = "Anomaly Composite VIMF - SON (1980-2020, AMM <= -0.5)")
}

# ---------------------------------------------------------------------------------------------------------- plotting just anomalies
plotting.anom2 <- function (data.anom, Spe.season, title.p, xlim=c(-82,-40), ylim=c(-20,15)){
  # data.anom must be the dataFrame, not the whole list
  
  data.anom.p <- data.anom[data.anom$Season==Spe.season,]
  p2 <- ggplot(data.anom.p, aes(lon,lat))+#facet_wrap(~Season,ncol=2)+
    # scale_y_continuous(position="right")+    geom_contour(aes(z = gh.z)) +
    geom_vector(aes(angle=Angle, mag=Anomaly),
                col="red",pivot=0.5, skip=1)+
    # scale_color_gradientn(colours = paleta)+ #brewer.pal(11,"Spectral")
    geom_polygon(data=basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
    coord_fixed(xlim=xlim,ylim=ylim)+
    labs(x="Long",y="Lat",title= title.p)+
    theme_bw()
  
  print(p2)
}

{plotting.anom2(data.anom = data.anom[["Pos"]], Spe.season = "DJF", title.p = "Anomaly Composite VIMF - DJF (1980-2020, AMM >= 1*SD)")

plotting.anom2(data.anom = data.anom[["Pos"]], Spe.season = "MAM", title.p = "Anomaly Composite VIMF - MAM (1980-2020, AMM >= 1*SD)")

plotting.anom2(data.anom = data.anom[["Pos"]], Spe.season = "JJA", title.p = "Anomaly Composite VIMF - JJA (1980-2020, AMM >= 1*SD)")

plotting.anom2(data.anom = data.anom[["Pos"]], Spe.season = "SON", title.p = "Anomaly Composite VIMF - SON (1980-2020, AMM >= 1*SD)", xlim=c(-80,-20), ylim=c(-20,15))

plotting.anom2(data.anom = data.anom[["Neg"]], Spe.season = "DJF", title.p = "Anomaly Composite VIMF - DJF (1980-2020, AMM <= -1*SD)")

plotting.anom2(data.anom = data.anom[["Neg"]], Spe.season = "MAM", title.p = "Anomaly Composite VIMF - MAM (1980-2020, AMM <= -1*SD)")

plotting.anom2(data.anom = data.anom[["Neg"]], Spe.season = "JJA", title.p = "Anomaly Composite VIMF - JJA (1980-2020, AMM <= -1*SD)")

plotting.anom2(data.anom = data.anom[["Neg"]], Spe.season = "SON", title.p = "Anomaly Composite VIMF - SON (1980-2020, AMM <= -1*SD)")}
