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
dates.m <- seq(as.Date("1980-01-01"),as.Date("2020-12-31"),by="month")
month.Multi <- format(dates.m,format="%m")
years <- factor(format(dates.m,format="%Y")); Years <- unique(years)
Season.y <- paste0(years[-1],"-",time2season(dates.m[-length(dates.m)]))
Year.s <- unique(Season.y)

seasons <- c("DJF","MAM","JJA","SON")

D.e5 <-seq(as.Date("1950-01-01"),as.Date("2020-12-31"),by="month")

Dm.sst <- seq(as.Date("1885-01-01"),as.Date("2020-12-30"),by="month")
## loading files ####
dir.mf <- "../../../01_DataSets/06_Wind/02_ERA5/"
W.e <- brick(paste0(dir.mf,"ERA5_u_v_950hPa-Atl_1950-2020.nc"),varname="u"); W.e <- W.e[[match(dates.m,D.e5)]]
W.n <- brick(paste0(dir.mf,"ERA5_u_v_950hPa-Atl_1950-2020.nc"),varname="v"); W.n <- W.n[[match(dates.m,D.e5)]]

basins <- shapefile("../../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")

Atl3 <- read.zoo(read.csv("../../../01_DataSets/01_SST/04_Cal_Indices/00A_AtlEN_SST.csv")[,-2], index.column = 1); Atl3 <- Atl3[dates.m]

## seasonal preprocessing ----
#------------------------------------------------------------------ Winds accummulation to seasons and separation
{W.e.s <- stackApply(W.e[[-nlayers(W.e)]], Season.y, fun=mean)
W.n.s <- stackApply(W.n[[-nlayers(W.n)]], Season.y, fun=mean)
W.e.s.t <- rasterToPoints(W.e.s)[,-c(1,2)]
W.n.s.t <- rasterToPoints(W.n.s)[,-c(1,2)]

Coord <- as.data.frame(rasterToPoints(W.e.s[[1]])[,c(1,2)])

W.E.seasons <- list()
as.data.frame(matrix(NA, ncol=nrow(W.e.s.t), nrow=length(Years))) -> W.E.seasons[["DJF"]] -> W.E.seasons[["MAM"]] -> W.E.seasons[["JJA"]] -> W.E.seasons[["SON"]]
W.E.seasons[["DJF"]][2:length(Years),] <- t(W.e.s.t[, which(substr(Year.s,6,8) == "DJF")[-1]])
W.E.seasons[["MAM"]][1:length(Years),] <- t(W.e.s.t[, which(substr(Year.s,6,8) == "MAM")])
W.E.seasons[["JJA"]][1:length(Years),] <- t(W.e.s.t[, which(substr(Year.s,6,8) == "JJA")])
W.E.seasons[["SON"]][1:length(Years),] <- t(W.e.s.t[, which(substr(Year.s,6,8) == "SON")])

W.N.seasons <- list()
as.data.frame(matrix(NA, ncol=nrow(W.n.s.t), nrow=length(Years))) -> W.N.seasons[["DJF"]] -> W.N.seasons[["MAM"]] -> W.N.seasons[["JJA"]] -> W.N.seasons[["SON"]]
W.N.seasons[["DJF"]][2:length(Years),] <- t(W.n.s.t[, which(substr(Year.s,6,8) == "DJF")[-1]])
W.N.seasons[["MAM"]][1:length(Years),] <- t(W.n.s.t[, which(substr(Year.s,6,8) == "MAM")])
W.N.seasons[["JJA"]][1:length(Years),] <- t(W.n.s.t[, which(substr(Year.s,6,8) == "JJA")])
W.N.seasons[["SON"]][1:length(Years),] <- t(W.n.s.t[, which(substr(Year.s,6,8) == "SON")])
}#------------------ Winds separation

#------------------------------------------------------------------ seasonal mean
{
  Mean.MtE.seasons <- list()
  for ( i in seasons){Mean.MtE.seasons[[i ]] <- apply(W.E.seasons[[i ]],2, mean, na.rm=T)}
  Mean.MtE.seasons <- simplify2array(Mean.MtE.seasons)
  
  Mean.MtN.seasons <- list()
  for ( i in seasons){Mean.MtN.seasons[[i ]] <- apply(W.N.seasons[[i ]],2, mean, na.rm=T)}
  Mean.MtN.seasons <- simplify2array(Mean.MtN.seasons)
  
  Mean.Angle.seasons <- Mean.MtE.seasons
  Mean.Mag.seasons <- Mean.MtE.seasons
  
  for ( i in seasons){
    Mean.Angle.seasons[,i] <- atan2(dlat(Mean.MtN.seasons[,i]), dlon(Mean.MtE.seasons[,i],Coord[,"y"]))*180/pi
    Mean.Mag.seasons[,i] <- Mag(Mean.MtN.seasons[,i], Mean.MtE.seasons[,i])
  }
}#----------- seasonal mean

#------------------------------------------------------------------ SST accummulation to seasons and separation
{
  # Atl3.s <- aggregate(Atl3[-length(Atl3)], by= list(Season.y), FUN= mean)
  # Atl3.seasons <- data.frame(matrix(NA, ncol=4, nrow= length(Years))); colnames(Atl3.seasons) <- seasons; rownames(Atl3.seasons) <- Years
  # Atl3.seasons[2:length(Years),"DJF"] <- Atl3.s[which(substr(index(Atl3.s),6,8) == "DJF") [-1]]
  # Atl3.seasons[   ,"MAM"] <- Atl3.s[which(substr(index(Atl3.s),6,8) == "MAM")]
  # Atl3.seasons[   ,"JJA"] <- Atl3.s[which(substr(index(Atl3.s),6,8) == "JJA")]
  # Atl3.seasons[   ,"SON"] <- Atl3.s[which(substr(index(Atl3.s),6,8) == "SON")]
  
  Atl3.seasons <- read.csv("../Atl3_std_1980-2020.csv")[,-1]
  
}#--------- SST accummulation to seasons and separation

#### identification of periods SST higher and lower than 1*SD --------
#--------------------------------------------------------------------------- initialization of datasets
{W.E.anom <- list(); W.N.anom <- list(); W.Mag.anom <- list(); W.Angle.anom <- list()
matrix(NA,nrow=nrow(Coord),ncol=4) -> W.Mag.anom[["Pos"]] -> W.Mag.anom[["Neg"]] -> W.Angle.anom [["Pos"]] -> W.Angle.anom [["Neg"]] ->
  W.E.anom[["Pos"]] -> W.E.anom[["Neg"]] -> W.N.anom[["Pos"]] -> W.N.anom[["Neg"]]
W.E.anom <- lapply(W.E.anom, 'colnames<-',seasons); W.N.anom <- lapply(W.N.anom, 'colnames<-',seasons)
W.Mag.anom <- lapply(W.Mag.anom, 'colnames<-',seasons); W.Angle.anom <- lapply(W.Angle.anom, 'colnames<-',seasons)}

Atl3.bool <- list()

#---------------------------------------------------------------------- identification of periods with SST higher than 0.5
{Atl3.bool[["Pos"]] <- data.frame(matrix(NA, ncol=4, nrow= length(Years))); colnames(Atl3.bool[["Pos"]]) <- seasons; rownames(Atl3.bool[["Pos"]]) <- Years
Atl3.bool[["Pos"]][1:length(Years),1:4] <- FALSE
Atl3.bool[["Pos"]][which(Atl3.seasons >=1, arr.ind = T)] <- TRUE}

#---------------------------------------------------------------------- identification of periods with SST Lower than -0.5
{Atl3.bool[["Neg"]] <- data.frame(matrix(NA, ncol=4, nrow= length(Years))); colnames(Atl3.bool[["Neg"]]) <- seasons; rownames(Atl3.bool[["Neg"]]) <- Years
  Atl3.bool[["Neg"]][1:length(Years),1:4] <- FALSE
  Atl3.bool[["Neg"]][which(Atl3.seasons <=-1, arr.ind = T)] <- TRUE}

a <- t(sapply(Atl3.bool,function(x){apply(x,2,sum)}))
write.csv(a,"Atl3_events_1980-2020.csv")
#------------------------------------------------------------------------ identification MT anomalies in those periods
for ( i in seasons){
  W.E.anom[["Pos"]][,i] <- apply( W.E.seasons[[i]][Atl3.bool[["Pos"]][,i],]  ,2,mean) - Mean.MtE.seasons[,i]
  W.N.anom[["Pos"]][,i] <- apply( W.N.seasons[[i]][Atl3.bool[["Pos"]][,i],]  ,2,mean) - Mean.MtN.seasons[,i]
  
  W.Angle.anom[["Pos"]][,i] <- atan2(dlat(W.N.anom[["Pos"]][,i]), dlon(W.E.anom[["Pos"]][,i], Coord[,"y"]))*180/pi
  W.Mag.anom[["Pos"]][,i] <- Mag(W.N.anom[["Pos"]][,i],W.E.anom[["Pos"]][,i])
  
  W.E.anom[["Neg"]][,i] <- apply( W.E.seasons[[i]][Atl3.bool[["Neg"]][,i],]  ,2,mean) - Mean.MtE.seasons[,i]
  W.N.anom[["Neg"]][,i] <- apply( W.N.seasons[[i]][Atl3.bool[["Neg"]][,i],]  ,2,mean) - Mean.MtN.seasons[,i]
  
  W.Angle.anom[["Neg"]][,i] <- atan2(dlat(W.N.anom[["Neg"]][,i]), dlon(W.E.anom[["Neg"]][,i], Coord[,"y"]))*180/pi
  W.Mag.anom[["Neg"]][,i] <- Mag(W.N.anom[["Neg"]][,i],W.E.anom[["Neg"]][,i])
}
colnames(Coord) <- c("lon", "lat")
W.E.anom[["Pos"]] <- cbind(Coord, W.E.anom[["Pos"]])
W.N.anom[["Pos"]] <- cbind(Coord, W.E.anom[["Pos"]])
W.E.anom[["Neg"]] <- cbind(Coord, W.E.anom[["Neg"]])
W.N.anom[["Neg"]] <- cbind(Coord, W.E.anom[["Neg"]])

data.anom.uv <- list(W.E.anom= W.E.anom, W.N.anom= W.N.anom)
save(data.anom.uv,"data.anom.uv",file="01_AnomComp_Atl3_Winds_uv.RData")
##### Data transformation for plotting ------
#------------------------------------------------------------------------ Seasonal mean transformation
{
  Mean.Angle.seasons <- cbind(Coord, Mean.Angle.seasons)
  Mean.Mag.seasons <- cbind(Coord, Mean.Mag.seasons)
  
  data.mean <- melt(Mean.Angle.seasons, id=c("lon","lat"))
  data.mean <- cbind(data.mean, melt(Mean.Mag.seasons, id=c("lon","lat"))[4])
  colnames(data.mean) <- c("lon", "lat","Season","Angle","Magnitud")
}

#------------------------------------------------------------------------ Anomalies in positive Atl3
data.anom <- list()
{
  W.Angle.anom[["Pos"]] <- cbind(Coord,W.Angle.anom[["Pos"]])
  W.Mag.anom[["Pos"]] <- cbind(Coord, W.Mag.anom[["Pos"]])

  data.anom[["Pos"]] <- melt(W.Angle.anom[["Pos"]], id=c("lon","lat"))
  data.anom[["Pos"]] <- cbind(data.anom[["Pos"]], melt(W.Mag.anom[["Pos"]], id=c("lon","lat"))[4])
  colnames(data.anom[["Pos"]]) <- c("lon", "lat","Season","Angle","Anomaly")
  }

#------------------------------------------------------------------------ Anomalies in Negative Atl3
{
  W.Angle.anom[["Neg"]] <- cbind(Coord,W.Angle.anom[["Neg"]])
  W.Mag.anom[["Neg"]] <- cbind(Coord, W.Mag.anom[["Neg"]])
  
  data.anom[["Neg"]] <- melt(W.Angle.anom[["Neg"]], id=c("lon","lat"))
  data.anom[["Neg"]] <- cbind(data.anom[["Neg"]], melt(W.Mag.anom[["Neg"]], id=c("lon","lat"))[4])
  colnames(data.anom[["Neg"]]) <- c("lon", "lat","Season","Angle","Anomaly")
}

save(data.anom,"data.anom",file="01_AnomComp_Atl3_Winds.RData")
load("01_AnomComp_Atl3_Winds.RData")
#### plotting ----
paleta <- c(brewer.pal(9,"YlGnBu"),"#81007F"); paleta[1:3] <- c("#9e0142","#f46d43","#fad366")
# ---------------------------------------------------------------------------------------------------------- filtering arrows for better visualization
{
  sele <- list()
  # sele[["lon"]] <- seq(min(data$lon),max(data$lon),by=0.25)
  # sele[["lat"]] <- seq(min(data$lat),max(data$lat),by=0.25)
  sele[["lon"]] <- seq(-83,10,by=0.5)
  sele[["lat"]] <- seq(-21,20,by=0.5)
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
#   plotting.anom(data.mean, data.anom = data.anom[["Pos"]], Spe.season = "DJF", title.p = "Anomaly Composite Winds - DJF (1980-2020, Atl3 >= 0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Pos"]], Spe.season = "MAM", title.p = "Anomaly Composite Winds - MAM (1980-2020, Atl3 >= 0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Pos"]], Spe.season = "JJA", title.p = "Anomaly Composite Winds - JJA (1980-2020, Atl3 >= 0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Pos"]], Spe.season = "SON", title.p = "Anomaly Composite Winds - SON (1980-2020, Atl3 >= 0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Neg"]], Spe.season = "DJF", title.p = "Anomaly Composite Winds - DJF (1980-2020, Atl3 <= -0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Neg"]], Spe.season = "MAM", title.p = "Anomaly Composite Winds - MAM (1980-2020, Atl3 <= -0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Neg"]], Spe.season = "JJA", title.p = "Anomaly Composite Winds - JJA (1980-2020, Atl3 <= -0.5)")
# 
# plotting.anom(data.mean, data.anom = data.anom[["Neg"]], Spe.season = "SON", title.p = "Anomaly Composite Winds - SON (1980-2020, Atl3 <= -0.5)")
}

# ---------------------------------------------------------------------------------------------------------- plotting just anomalies
plotting.anom2 <- function (data.anom, Spe.season, title.p, xlim=c(-82,10), ylim=c(-25,20)){
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
    theme_bw()+theme(legend.position=c(0.85,0.85),legend.background = element_rect(linetype="solid", colour ="black"))
  
  print(p2)
}

{plotting.anom2(data.anom = data.anom[["Pos"]], Spe.season = "DJF", title.p = "Anomaly Composite Winds - DJF (1980-2020, Atl3 >= 1*SD)" )

plotting.anom2(data.anom = data.anom[["Pos"]], Spe.season = "MAM", title.p = "Anomaly Composite Winds - MAM (1980-2020, Atl3 >= 1*SD)" )

plotting.anom2(data.anom = data.anom[["Pos"]], Spe.season = "JJA", title.p = "Anomaly Composite Winds - JJA (1980-2020, Atl3 >= 1*SD)" )

plotting.anom2(data.anom = data.anom[["Pos"]], Spe.season = "SON", title.p = "Anomaly Composite Winds - SON (1980-2020, Atl3 >= 1*SD)" )

plotting.anom2(data.anom = data.anom[["Neg"]], Spe.season = "DJF", title.p = "Anomaly Composite Winds - DJF (1980-2020, Atl3 <= -1*SD)" )

plotting.anom2(data.anom = data.anom[["Neg"]], Spe.season = "MAM", title.p = "Anomaly Composite Winds - MAM (1980-2020, Atl3 <= -1*SD)" )

plotting.anom2(data.anom = data.anom[["Neg"]], Spe.season = "JJA", title.p = "Anomaly Composite Winds - JJA (1980-2020, Atl3 <= -1*SD)" )

plotting.anom2(data.anom = data.anom[["Neg"]], Spe.season = "SON", title.p = "Anomaly Composite Winds - SON (1980-2020, Atl3 <= -1*SD)" )}
