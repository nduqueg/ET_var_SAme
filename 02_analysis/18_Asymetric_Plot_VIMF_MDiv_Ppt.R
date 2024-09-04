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

## loading ----
VIMF.anom <- list(); MDiv.anom <- list(); Ppt.anom <- list()
load("./01_VIMF/01_AnomComp_Atl3_VIMF_uv.RData"); VIMF.anom[["Atl3"]] <- data.anom.uv; rm(data.anom.uv)
load("./02_VIMFDiv/02_Atl3_divVIM.RData"); MDiv.anom[["Atl3"]] <- data.anom2; rm(data.anom2)
load("./03_Ppt/03_Atl3_Comp_Ppt.RData"); Ppt.anom[["Atl3"]] <- data.anom2; rm(data.anom2)

load("./01_VIMF/01_AnomComp_AMM_VIMF_uv.RData"); VIMF.anom[["AMM"]] <- data.anom.uv; rm(data.anom.uv)
load("./02_VIMFDiv/02_AMM_divVIM.RData"); MDiv.anom[["AMM"]] <- data.anom2; rm(data.anom2)
load("./03_Ppt/03_AMM_Comp_Ppt.RData"); Ppt.anom[["AMM"]] <- data.anom2; rm(data.anom2)

### selection VIMF for ease visualization----

sele <- list()
# sele[["lon"]] <- seq(min(data$lon),max(data$lon),by=0.25)
# sele[["lat"]] <- seq(min(data$lat),max(data$lat),by=0.25)
sele[["lon"]] <- seq(-83,-29,by=0.75)
sele[["lat"]] <- seq(-21,16.5,by=0.75)

sele.fun <- function(VIMF.anom, sele){
  VIMF.anom[["Neg"]] <- VIMF.anom[["Neg"]][ which(VIMF.anom[["Neg"]]$lon %in% sele$lon) ,]; VIMF.anom[["Neg"]] <- VIMF.anom[["Neg"]][ which(VIMF.anom[["Neg"]]$lat %in% sele$lat) ,]
  VIMF.anom[["Pos"]] <- VIMF.anom[["Pos"]][ which(VIMF.anom[["Pos"]]$lon %in% sele$lon) ,]; VIMF.anom[["Pos"]] <- VIMF.anom[["Pos"]][ which(VIMF.anom[["Pos"]]$lat %in% sele$lat) ,]
  return(VIMF.anom)
}

VIMF.anom[["Atl3"]] <- lapply(VIMF.anom[["Atl3"]], sele.fun, sele)
VIMF.anom[["AMM"]] <- lapply(VIMF.anom[["AMM"]], sele.fun, sele)

### Asymetry ----
# Asymetry for each component (u and v)
cal.asym.VIMF <- function(x){
  
  cord <- x$Mt.E.anom$Pos[,5:6]
  VIMF.asym <- lapply(x, function(y,cord.y){
    asym <- y$Pos[,1:4] + y$Neg[,1:4]
    asym <- cbind.data.frame(cord.y,asym)
    return(asym)
  }, cord)
  
  return(VIMF.asym)
}
VIMF.Asym <- lapply(VIMF.anom, cal.asym.VIMF)

# Calculate the Angle and Magnitud of the asymetry
cal.asym.AngMag <- function(VIMF.asym){
  cal.asym.Ang <- list()
  Coord <- VIMF.asym[["Mt.N.anom"]][,1:2]
  cal.asym.Ang[["Angle"]] <- VIMF.asym[["Mt.N.anom"]][,1:2]
  cal.asym.Ang[["Mag"]] <- VIMF.asym[["Mt.N.anom"]][,1:2]
  
  for (i in 1:4){
    cal.asym.Ang[["Angle"]][,i+2] <- atan2(dlat(VIMF.asym[["Mt.N.anom"]][,i+2]), dlon(VIMF.asym[["Mt.E.anom"]][,i+2], Coord[,"lat"]))*180/pi
    
    cal.asym.Ang[["Mag"]][,i+2] <- Mag(VIMF.asym[["Mt.N.anom"]][,i+2],VIMF.asym[["Mt.E.anom"]][,i+2])
  }
  colnames(cal.asym.Ang[["Angle"]])[3:6] <- seasons; colnames(cal.asym.Ang[["Mag"]])[3:6] <- seasons
  
  cal.asym.Ang[["Angle"]] <- melt(cal.asym.Ang[["Angle"]], id=c("lon","lat"))
  cal.asym.Ang[["Mag"]] <- melt(cal.asym.Ang[["Mag"]], id=c("lon","lat"))
  
  cal.asym.return <- cbind(cal.asym.Ang[["Angle"]], cal.asym.Ang[["Mag"]][,"value"])
  colnames(cal.asym.return)[3:5] <- c("Season","Angle","Mag")
  return(cal.asym.return) 
}
VIMF.asym.p <- lapply(VIMF.Asym, cal.asym.AngMag)

# calculate asymetry of Ppt and MDiv

cal.asym.ppt <- function(Ppt.Anom){
  Ppt.asym.p <- Ppt.Anom[Ppt.Anom$Dir=="Pos",-5]
  Ppt.asym.p[,"Anomaly"] <- Ppt.Anom[Ppt.Anom$Dir=="Pos","Anomaly"] + Ppt.Anom[Ppt.Anom$Dir=="Neg","Anomaly"]
  return(Ppt.asym.p)
}
Ppt.Asym.p <- lapply(Ppt.anom, cal.asym.ppt)

cal.asym.MDiv <- function(MDiv.Anom){
  MDiv.asym.p <- MDiv.Anom[ MDiv.Anom$Dir=="Pos",-5]
  MDiv.asym.p[,"Anomaly"] <- MDiv.Anom[MDiv.Anom$Dir=="Pos","Anomaly"] + MDiv.Anom[ MDiv.Anom$Dir=="Neg","Anomaly"]
  return(MDiv.asym.p)
}
MDiv.asym.p <- lapply(MDiv.anom, cal.asym.MDiv)

#### plotting -----
basins <- shapefile("../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")

plotting.asym <- function ( VIMF.asym, MDiv.asym, Ppt.asym, Spe.season, title.p, xlim=c(-82,-40), ylim=c(-20,13), l.pos="bottom", y.axis=TRUE, x.axis=TRUE){
  # data.asym must be the dataFrame, not the whole list
  at.ppt <- c(-250,-150,-100,-75,-50,-25,25,50,75,100,150,250)
  at.ppt.v <- (at.ppt[-length(at.ppt)] - at.ppt[-1])/2 + at.ppt[-1]
  col2alpha <- function(someColor, alpha=1){ newColor <- col2rgb(someColor); if(alpha <=1) alpha <- alpha * 255 else warning("Alpha should be between 0 and 1"); apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}
  paleta <- brewer.pal(11,"BrBG"); paleta <- col2alpha(paleta, alpha = 0.5)
  
  VIMF.asym.p <- VIMF.asym[sapply(Spe.season,function(s,d) which(d == s),VIMF.asym$Season),]
  MDiv.asym.p <- MDiv.asym[sapply(Spe.season,function(s,d) which(d == s),MDiv.asym$Season),]
  
  data.class <- MDiv.asym.p$Anomaly>0; MDiv.asym.p <- cbind(MDiv.asym.p, data.class)
  
  Ppt.asym.p <- Ppt.asym[sapply(Spe.season,function(s,d) which(d == s),Ppt.asym$Season),]
  p <- ggplot(VIMF.asym.p, aes(lon,lat))+facet_wrap(.~Season, ncol=4)+
    # scale_y_continuous(position="right")+    geom_contour(aes(z = gh.z)) +
    geom_raster(data=Ppt.asym.p, aes(lon, lat, fill=Anomaly))+ scale_fill_stepsn(colours=paleta, breaks=at.ppt,
                                                                                   values=scales::rescale(at.ppt.v,from=range(at.ppt)),limits=c(min(at.ppt),max(at.ppt)),
                                                                                   guide=guide_colorsteps(even.steps = F, barheight=unit(7,"cm")))+ # 
    
    geom_polygon(data=basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.2)+
    
    geom_contour(data=dplyr::filter(MDiv.asym.p,data.class==T), aes(lon,lat, z=Anomaly),color="red",binwidth =3,linewidth=0.25)+ #  linetype=1,
    geom_contour(data=dplyr::filter(MDiv.asym.p,data.class==F), aes(lon,lat, z=Anomaly),color="blue",binwidth =3,linewidth=0.25)+
    
    geom_vector(aes(angle=Angle, mag=Mag),
                col="purple4",pivot=0.5, skip=1)+
    scale_mag(max = 50,name="VIMF")+
    # scale_color_gradientn(colours = paleta)+ #brewer.pal(11,"Spectral")
    coord_fixed(xlim=xlim,ylim=ylim)+
    labs(x="Long",y="Lat",title= title.p)+
    theme_bw()
  if (!y.axis){
    p <- p + theme(legend.position = l.pos,
                   axis.text.y= element_blank(), axis.title.y = element_blank()) # c(0.9,0.2)
  }else{
    p <- p + theme(legend.position = l.pos, legend.background = element_rect(linetype="solid", colour ="black"))
  }
  
}

P.Atl3 <- plotting.asym(VIMF.asym = VIMF.asym.p[["Atl3"]], MDiv.asym = MDiv.asym.p[["Atl3"]], Ppt.asym = Ppt.Asym.p[["Atl3"]],
                        Spe.season = c("JJA"), title.p = "Composite Asymetry (1980-2020, |Atl3| >= 1*SD)",xlim=c(-82,-31),l.pos="right", y.axis = FALSE)

# plotting.asym(VIMF.asym = VIMF.asym.p, MDiv.asym = MDiv.asym[["Atl3"]], Ppt.asym = Ppt.asym[["Atl3"]],
P.AMM <- plotting.asym(VIMF.asym = VIMF.asym.p[["AMM"]], MDiv.asym = MDiv.asym.p[["AMM"]], Ppt.asym = Ppt.Asym.p[["AMM"]],
                       Spe.season = c("MAM","JJA","SON"),title.p = "Composite Asymetry (1980-2020, |AMM| >= 1*SD)",xlim=c(-82,-31), l.pos="none")

library(gridExtra)
grid.arrange(P.AMM,P.Atl3, ncol=2, widths= c(2.1,1)) # 1700x500

#### pre-process 2 modes ----
VIMF.anom.p <- list()

VIMF.melt <- function (VIMF.anom){
  VIMF.anom.p <- melt(VIMF.anom, id=c("lon","lat","Season","Angle","Anomaly"))
  colnames(VIMF.anom.p)[6] <- "Dir"
  VIMF.anom.p$Dir <- factor(VIMF.anom.p$Dir)
  return(VIMF.anom.p)
}

VIMF.anom.p <- lapply(VIMF.anom, VIMF.melt)
VIMF.anom.p <- melt(VIMF.anom.p, id= c("lon","lat","Season","Angle","Anomaly","Dir"))
VIMF.anom.p$L1 <- factor(VIMF.anom.p$L1); VIMF.anom.p$Dir <- factor(VIMF.anom.p$Dir, levels = c("Pos","Neg"))

MDiv.anom.p <- melt(MDiv.anom, id= c("lon","lat","Season","Anomaly","Dir")); MDiv.anom.p$L1 <- factor(MDiv.anom.p$L1); MDiv.anom.p$Dir <- factor(MDiv.anom.p$Dir, levels = c("Pos","Neg"))
Ppt.anom.p <- melt(Ppt.anom, id= c("lon","lat","Season","Anomaly","Dir")); Ppt.anom.p$L1 <- factor(Ppt.anom.p$L1); Ppt.anom.p$Dir <- factor(Ppt.anom.p$Dir, levels = c("Pos","Neg"))



#### plotting 2 modes ----
plotting.anom2 <- function (VIMF.anom, MDiv.anom, Ppt.anom, Spe.season, title.p, xlim=c(-82,-40), ylim=c(-20,15)){
  # data.anom must be the dataFrame, not the whole list
  at.ppt <- c(-200,-150,-100,-75,-50,-25,25,50,75,100,150,200)
  
  VIMF.anom.p <- VIMF.anom[VIMF.anom$Season==Spe.season,]
  MDiv.anom.p <- MDiv.anom[MDiv.anom$Season==Spe.season,]
  
  data.class <- MDiv.anom.p$Anomaly>0; MDiv.anom.p <- cbind(MDiv.anom.p, data.class)
  
  Ppt.anom.p <- Ppt.anom[Ppt.anom$Season==Spe.season,]
  p2 <- ggplot(VIMF.anom.p, aes(lon,lat))+facet_grid(Dir~L1, switch = "y")+
    # scale_y_continuous(position="right")+    geom_contour(aes(z = gh.z)) +
    geom_raster(data=Ppt.anom.p, aes(lon, lat, fill=Anomaly),alpha=0.5)+ scale_fill_stepsn(colours=brewer.pal(11,"BrBG"), breaks=at.ppt,
                                                                                           values=c(0:12)/12,limits=c(-200,200),
                                                                                           guide=guide_colorsteps(even.steps = F, barwidth=unit(10,"cm")))+ # 
    
    geom_polygon(data=basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.2)+
    
    geom_contour(data=dplyr::filter(MDiv.anom.p,data.class==T), aes(lon,lat, z=Anomaly),color="red",binwidth =3,linewidth=0.25)+ #  linetype=1,
    geom_contour(data=dplyr::filter(MDiv.anom.p,data.class==F), aes(lon,lat, z=Anomaly),color="blue",binwidth =3,linewidth=0.25)+
    
    geom_vector(aes(angle=Angle, mag=Anomaly),
                col="purple4",pivot=0.5, skip=1)+
    # scale_color_gradientn(colours = paleta)+ #brewer.pal(11,"Spectral")
    coord_fixed(xlim=xlim,ylim=ylim)+
    labs(x="Long",y="Lat",title= title.p)+
    theme_bw()+theme(legend.position="bottom",legend.background = element_rect(linetype="solid", colour ="black")) # c(0.9,0.2)
  
  print(p2)
}

plotting.anom2(VIMF.anom = VIMF.anom.p, MDiv.anom = MDiv.anom.p, Ppt.anom = Ppt.anom.p,
               Spe.season = "DJF", title.p = "Anomaly Composite VIMF - DJF (1980-2020, |modes| >= 1*SD)",xlim=c(-82,-31))

plotting.anom2(VIMF.anom = VIMF.anom.p, MDiv.anom = MDiv.anom.p, Ppt.anom = Ppt.anom.p,
              Spe.season = "MAM", title.p = "Anomaly Composite VIMF - MAM (1980-2020, |modes| >= 1*SD)",xlim=c(-82,-31))

plotting.anom2(VIMF.anom = VIMF.anom.p, MDiv.anom = MDiv.anom.p, Ppt.anom = Ppt.anom.p,
               Spe.season = "JJA", title.p = "Anomaly Composite VIMF - JJA (1980-2020, |modes| >= 1*SD)",xlim=c(-82,-31))

plotting.anom2(VIMF.anom = VIMF.anom.p, MDiv.anom = MDiv.anom.p, Ppt.anom = Ppt.anom.p,
               Spe.season = "JJA", title.p = "Anomaly Composite VIMF - SON (1980-2020, |modes| >= 1*SD)",xlim=c(-82,-31))
