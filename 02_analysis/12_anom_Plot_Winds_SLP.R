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
basins <- shapefile("./01_data/hybas_sa_lev03_v1c.shp")
countries <- shapefile("./01_data/South_America.shp")

WIND.anom <- list(); MDiv.anom <- list(); SLP.anom <- list()
load("./01_data/06_Wind/01_AnomComp_AMM_Winds.RData"); WIND.anom[["AMM"]] <- data.anom; rm(data.anom)
load("./01_data/07_SLP/03_AMM_Comp_SLP.RData"); SLP.anom[["AMM"]] <- data.anom2; rm(data.anom2)
SLP.anom[["AMM"]]$lon[SLP.anom[["AMM"]]$lon>=180] <- SLP.anom[["AMM"]]$lon[SLP.anom[["AMM"]]$lon>=180] - 360

load("./01_data/06_Wind/01_AnomComp_Atl3_Winds.RData"); WIND.anom[["Atl3"]] <- data.anom; rm(data.anom)
load("./01_data/07_SLP/03_Atl3_Comp_SLP.RData"); SLP.anom[["Atl3"]] <- data.anom2; rm(data.anom2)
SLP.anom[["Atl3"]]$lon[SLP.anom[["Atl3"]]$lon>=180] <- SLP.anom[["Atl3"]]$lon[SLP.anom[["Atl3"]]$lon>=180] - 360

### selection WIND for ease visualization----

sele <- list()
# sele[["lon"]] <- seq(min(data$lon),max(data$lon),by=0.25)
# sele[["lat"]] <- seq(min(data$lat),max(data$lat),by=0.25)
sele[["lon"]] <- seq(-83,10,by=0.75)
sele[["lat"]] <- seq(-25,16.5,by=0.75)

sele.fun <- function(WIND.anom, sele){
  WIND.anom[["Neg"]] <- WIND.anom[["Neg"]][ which(WIND.anom[["Neg"]]$lon %in% sele$lon) ,]; WIND.anom[["Neg"]] <- WIND.anom[["Neg"]][ which(WIND.anom[["Neg"]]$lat %in% sele$lat) ,]
  WIND.anom[["Pos"]] <- WIND.anom[["Pos"]][ which(WIND.anom[["Pos"]]$lon %in% sele$lon) ,]; WIND.anom[["Pos"]] <- WIND.anom[["Pos"]][ which(WIND.anom[["Pos"]]$lat %in% sele$lat) ,]
  return(WIND.anom)
}

WIND.anom <- lapply(WIND.anom, sele.fun, sele)

#### plotting -----
WIND.anom.p <- melt(WIND.anom[["Atl3"]], id=c("lon","lat","Season","Angle","Anomaly")); colnames(WIND.anom.p)[6] <- "Dir"; WIND.anom.p$Dir <- factor(WIND.anom.p$Dir, levels = c("Pos", "Neg"))

plotting.anom <- function (WIND.anom, MDiv.anom, SLP.anom, Spe.season, title.p, xlim=c(-82,-30), ylim=c(-20,15), l.pos="bottom", y.axis=TRUE, x.axis=TRUE){
  # data.anom must be the dataFrame, not the whole list
  at.slp <- c(-300,-250,-200,-150,-100,-50,50,100,150,200,250,300)
  
  WIND.anom.p <- WIND.anom[sapply(Spe.season, function(s,d) which(d == s) ,WIND.anom$Season),]

  SLP.anom.p <- SLP.anom[sapply(Spe.season, function(s,d) which(d == s) ,SLP.anom$Season),]
  
  Phase.labs <- c("Positive", "Negative"); names(Phase.labs) <- c("Pos","Neg")
  Season.labs <- c("MAM (Austral Autumn)", "JJA (Austral Winter)", "SON (Austral Spring)"); names(Season.labs) <- c("MAM","JJA","SON")
  
  p2 <- ggplot(WIND.anom.p, aes(lon,lat))+facet_grid(Season ~ Dir ,switch = "y", labeller= labeller(Dir = Phase.labs, Season = Season.labs))+
    # scale_y_continuous(position="right")+    geom_contour(aes(z = gh.z)) +
    geom_raster(data=SLP.anom.p, aes(lon, lat, fill=Anomaly),alpha=0.5)+ scale_fill_stepsn(colours=brewer.pal(11,"PuOr"), breaks=at.slp,
                                                                                           values=c(0:12)/12,limits=c(min(at.slp),max(at.slp)),
                                                                                           guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")))+ # 
    scale_y_continuous(position="right")+
    geom_polygon(data=basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.5)+
    geom_polygon(data=countries,aes(x=long,y=lat, group=group),linetype="dashed", colour="black",fill="NA",size=0.05)+
    
    geom_vector(aes(angle=Angle, mag=Anomaly), col="purple4",pivot=0.5, skip=1)+
    scale_mag(max = 2,name="Wind Speed")+
    # scale_color_gradientn(colours = paleta)+ #brewer.pal(11,"Spectral")
    coord_fixed(xlim=xlim,ylim=ylim)+
    labs(x="Long",y="Lat",title= title.p)+
    theme_bw()
  
  if(y.axis){
    p2 <- p2 +theme(legend.position=l.pos,legend.background = element_rect(linetype="solid", colour ="black")) # c(0.9,0.2)
  }else{
    p2 <- p2 + theme(legend.position=l.pos,legend.background = element_rect(linetype="solid", colour ="black"), # , legend.key.size = unit(1,"cm")
                   axis.text.y= element_blank(), axis.title.y = element_blank()) # c(0.9,0.2)
  }
  
  if(x.axis){
    p2 <- p2 +theme(legend.position=l.pos,legend.background = element_rect(linetype="solid", colour ="black")) # c(0.9,0.2)
  }else{
    p2 <- p2 + theme(legend.position=l.pos,legend.background = element_rect(linetype="solid", colour ="black"), # , legend.key.size = unit(1,"cm")
                   axis.text.x= element_blank(), axis.title.x = element_blank()) # c(0.9,0.2)
  }
  
  p2 <- p2 + theme(strip.text = element_text( size=14))
  # print(p2)
  return(p2)
}

p.Atl3 <- plotting.anom(WIND.anom = WIND.anom.p, SLP.anom = SLP.anom[["Atl3"]],
              Spe.season = "JJA", title.p = "B) |Atl3| >= 1*SD", l.pos= "right")

WIND.anom.p <- melt(WIND.anom[["AMM"]], id=c("lon","lat","Season","Angle","Anomaly")); colnames(WIND.anom.p)[6] <- "Dir"; WIND.anom.p$Dir <- factor(WIND.anom.p$Dir, levels = c("Pos", "Neg"))

p.AMM <- plotting.anom(WIND.anom = WIND.anom.p, SLP.anom = SLP.anom[["AMM"]],
              Spe.season = c("MAM","JJA","SON"), title.p = "A) Anomaly Composite (1980-2020) |AMM| >= 1*SD", l.pos= "right", y.axis = T )

library(gridExtra)
grid.arrange(p.AMM,p.Atl3, ncol=2, widths= c(2.1,1)) # 1400x1400


#### pre-process 2 modes ----
WIND.anom.p <- list()

WIND.melt <- function (WIND.anom){
  WIND.anom.p <- melt(WIND.anom, id=c("lon","lat","Season","Angle","Anomaly"))
  colnames(WIND.anom.p)[6] <- "Dir"
  WIND.anom.p$Dir <- factor(WIND.anom.p$Dir)
  return(WIND.anom.p)
}

WIND.anom.p <- lapply(WIND.anom, WIND.melt)
WIND.anom.p <- melt(WIND.anom.p, id= c("lon","lat","Season","Angle","Anomaly","Dir"))
WIND.anom.p$L1 <- factor(WIND.anom.p$L1); WIND.anom.p$Dir <- factor(WIND.anom.p$Dir, levels = c("Pos","Neg"))

SLP.anom.p <- melt(SLP.anom, id= c("lon","lat","Season","Anomaly","Dir")); SLP.anom.p$L1 <- factor(SLP.anom.p$L1); SLP.anom.p$Dir <- factor(SLP.anom.p$Dir, levels = c("Pos","Neg"))



#### plotting 2 modes ----
plotting.anom2 <- function (WIND.anom,  SLP.anom, Spe.season, title.p, xlim=c(-82,0), ylim=c(-25,15)){
  # data.anom must be the dataFrame, not the whole list
  at.ppt <- c(-300,-250,-200,-150,-100,-50,50,100,150,200,250,300)
  
  WIND.anom.p <- WIND.anom[WIND.anom$Season==Spe.season,]
 
  SLP.anom.p <- SLP.anom[SLP.anom$Season==Spe.season,]
  p2 <- ggplot(WIND.anom.p, aes(lon,lat))+facet_grid(Dir~L1, switch = "y")+
    # scale_y_continuous(position="right")+    geom_contour(aes(z = gh.z)) +
    geom_raster(data=SLP.anom.p, aes(lon, lat, fill=Anomaly),alpha=0.5)+ scale_fill_stepsn(colours=brewer.pal(11,"PuOr"), breaks=at.ppt,
                                                                                           values=c(0:12)/12,limits=c(-300,300),
                                                                                           guide=guide_colorsteps(even.steps = F, barwidth=unit(10,"cm")))+ # 
    
    geom_polygon(data=basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.2)+

    geom_vector(aes(angle=Angle, mag=Anomaly),
                col="purple4",pivot=0.5, skip=1)+
    # scale_color_gradientn(colours = paleta)+ #brewer.pal(11,"Spectral")
    coord_fixed(xlim=xlim,ylim=ylim)+
    labs(x="Long",y="Lat",title= title.p)+
    theme_bw()+theme(legend.position="bottom",legend.background = element_rect(linetype="solid", colour ="black")) # c(0.9,0.2)
  
  print(p2)
}

plotting.anom2(WIND.anom = WIND.anom.p, SLP.anom = SLP.anom.p,
               Spe.season = "DJF", title.p = "Anomaly Composite Winds 950hPa & SLP - DJF (1980-2020, |modes| >= 1*SD)" )

plotting.anom2(WIND.anom = WIND.anom.p, SLP.anom = SLP.anom.p,
              Spe.season = "MAM", title.p = "Anomaly Composite Winds 950hPa & SLP - MAM (1980-2020, |modes| >= 1*SD)" )

plotting.anom2(WIND.anom = WIND.anom.p, SLP.anom = SLP.anom.p,
               Spe.season = "JJA", title.p = "Anomaly Composite Winds 950hPa & SLP - JJA (1980-2020, |modes| >= 1*SD)" )

plotting.anom2(WIND.anom = WIND.anom.p, SLP.anom = SLP.anom.p,
               Spe.season = "SON", title.p = "Anomaly Composite Winds 950hPa & SLP - SON (1980-2020, |modes| >= 1*SD)" )
