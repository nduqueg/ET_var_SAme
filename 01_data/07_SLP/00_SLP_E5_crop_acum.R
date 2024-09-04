rm(list=ls())
cat("\014")

library(reshape)
library(raster)
library(ncdf4)
library(RColorBrewer)
library(hydroTSM)

# DAtes ----
dates.m <- seq(as.Date("1980-01-01"),as.Date("2020-12-31"),by="month")
month.Multi <- format(dates.m,format="%m")
years <- factor(format(dates.m,format="%Y")); Years <- unique(years)
Season.y <- paste0(years[-1],"-",time2season(dates.m[-length(dates.m)]))
Year.s <- unique(Season.y)

D.e5 <-seq(as.Date("1959-01-01"),as.Date("2020-12-31"),by="month")

#### loading, croping & stacking the files ####
P <- brick("ERA5_SLP_1959-2020.nc")
P <- crop(P,extent(c(360-90,360,-40,30)))
writeRaster(P, filename = "ERA5_SLP_1959-2020_crop.nc",format="CDF",overwrite=T)

Colom <- shapefile("../../South_America/South_America.shp")
arrow <- list("SpatialPolygonsRescale", layout.north.arrow(), 
              offset = c(-45,10), scale = 2, first=F)
countries.layer <- list("sp.lines", Colom, col = "Black", first=F)

#### accummulating monthly #####
P.y <- stackApply(P, years, fun=mean, na.rm = T)#; P.y <- P.y *mask_SA
writeRaster(P.y,filename = "ERA5_SLP_yearly_crop.nc",format="CDF",overwrite=T)
P.yM <- stackApply(P.y, rep(1,nlayers(P.y)),fun=mean)


P <- P[[match(dates.m, D.e5)]]

P.mM <- stackApply(P, month.Multi, fun=mean)
writeRaster(P.mM,filename = "ERA5_SLP_MonMulti_crop.nc",format="CDF",overwrite=T)



###### plotting monthly ####
at.m <- seq(min(cellStats(P.mM,stat="min")), max(cellStats(P.mM,stat = "max")), length.out = 12)
#at.m <- round(at.m,2); at.m[1] <- floor(at.m[1]); at.m[12] <- ceiling(at.m[12])
#at.m <- c(0,50,100,200,300,400,600,800,1000,1200,1500,2000)
spplot(P.mM,scales=list(draw=T),col.regions=brewer.pal(11,"PRGn"),sp.layout=list(arrow,countries.layer),
       at= at.m,
       names.attr= c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
       colorkey=list(height=1,
                     labels= list( labels =round(at.m,2), at=at.m)),
       main="SLP monthly long-term mean [ERA5]")



#at.m <- c(10,100,200,300,400,600,800,1000,1500,2000,3000,10000)
spplot(P.yM,scales=list(draw=T),col.regions=brewer.pal(11,"Spectral"),sp.layout=list(arrow,countries.layer),
       at= at.m,
       colorkey=list(height=1,
                     labels= list( labels =round(at.m,2), at=at.m)),
       main="SLP Annual long-term mean [ERA5]")

#### seasonal preprocessing ----
#mask.P <- P[[1]]*0+1
P.s <- stackApply(P[[-nlayers(P)]], Season.y, fun=mean)#; P.s <- P.s*mask.P
P.s.t <- rasterToPoints(P.s)[,-c(1,2)]/100
Coord <- as.data.frame(rasterToPoints(P.s[[1]])[,c(1,2)])

SLP.seasons <- list()
as.data.frame(matrix(NA, ncol=nrow(P.s.t), nrow=length(Years))) -> SLP.seasons[["DJF"]] -> SLP.seasons[["MAM"]] -> SLP.seasons[["JJA"]] -> SLP.seasons[["SON"]]
SLP.seasons[["DJF"]][2:length(Years),] <- t(P.s.t[, which(substr(Year.s,6,8) == "DJF")[-1]])
SLP.seasons[["MAM"]][1:length(Years),] <- t(P.s.t[, which(substr(Year.s,6,8) == "MAM")])
SLP.seasons[["JJA"]][1:length(Years),] <- t(P.s.t[, which(substr(Year.s,6,8) == "JJA")])
SLP.seasons[["SON"]][1:length(Years),] <- t(P.s.t[, which(substr(Year.s,6,8) == "SON")])

Season <- c("DJF","MAM","JJA","SON")
Mean.seasons <- list()
Mean.seasons[["DJF"]] <- cbind(Coord,apply(SLP.seasons[["DJF"]],2, mean, na.rm=T),Season[1])
Mean.seasons[["MAM"]] <- cbind(Coord,apply(SLP.seasons[["MAM"]],2, mean, na.rm=T),Season[2])
Mean.seasons[["JJA"]] <- cbind(Coord,apply(SLP.seasons[["JJA"]],2, mean, na.rm=T),Season[3])
Mean.seasons[["SON"]] <- cbind(Coord,apply(SLP.seasons[["SON"]],2, mean, na.rm=T),Season[4])
Mean.seasons <- lapply(Mean.seasons,`colnames<-`,c("x","y","P","Season"))

Sd.seasons <- list()
Sd.seasons[["DJF"]] <- apply(SLP.seasons[["DJF"]],2, function(x){Sd <- sd(x,na.rm = T); return(Sd)})
Sd.seasons[["MAM"]] <- apply(SLP.seasons[["MAM"]],2, function(x){Sd <- sd(x,na.rm = T); return(Sd)})
Sd.seasons[["JJA"]] <- apply(SLP.seasons[["JJA"]],2, function(x){Sd <- sd(x,na.rm = T); return(Sd)})
Sd.seasons[["SON"]] <- apply(SLP.seasons[["SON"]],2, function(x){Sd <- sd(x,na.rm = T); return(Sd)})
for(i in 1:4){Sd.seasons[[i]] <- cbind(Coord,Sd.seasons[[i]] , Season[i])}
Sd.seasons <- lapply(Sd.seasons,`colnames<-`,c("x","y","P","Season"))

Mean.s.t <- do.call(rbind,Mean.seasons)  # as we don't want to destroy the dataframes that are inside the list, we have to paste them rather than MELT THEM
Var <- "Mean"; Mean.s.t <- cbind(Mean.s.t, Var)

Sd.s.t <- do.call(rbind,Sd.seasons)
Var <- "Coef. Var."; Sd.s.t <- cbind(Sd.s.t, Var)

#data <- rbind(Mean.s.t, Sd.s.t)
Mean.s.t$x <- Mean.s.t$x - 360
Sd.s.t$x <- Sd.s.t$x - 360

#### Seasonal plotting ----

library(maps);library(mapdata); library(ggplot2)
world <- map_data("world")
at.m <-seq(min(Mean.s.t$P),max(Mean.s.t$P),length.out = 12)
p1<-ggplot(Mean.s.t)+
  geom_raster(aes(x=x,y=y,fill=P))+facet_wrap(.~Season,ncol=4)+
  scale_fill_stepsn(colours=brewer.pal(11,"PRGn"), breaks=at.m,
                    values=(at.m-min(at.m))/(max(at.m)-min(at.m)),#limits=c(0,3000),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")))+
  #geom_sf(data=world,fill="transparent")+
  geom_polygon(data=world,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
  coord_fixed(xlim=c(-90,0),ylim=c(-40,30))+labs(x="Long",y="Lat",title="Mean seasonal SLP ERA5 (1980-2020)")+
  theme_bw()


at.m <- round(seq(0,max(Sd.s.t$P),length.out = 12),2)
#Sd.s.t$P <- log(Sd.s.t$P)
p2<-
  ggplot(Sd.s.t)+
  geom_raster(aes(x=x,y=y,fill=P))+facet_wrap(.~Season,ncol=4)+
  scale_fill_stepsn(colours=rev(brewer.pal(11,"Spectral")), breaks=at.m,
                    values=at.m/max(at.m),limits=c(0,max(at.m)),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")))+
  #geom_sf(data=world,fill="transparent")+
  geom_polygon(data=world,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
  coord_fixed(xlim=c(-90,0),ylim=c(-40,30))+labs(x="Long",y="Lat",title="Seasonal Stan.Dev. SLP ERA5 (1980-2020)")+
  theme_bw()

library(gridExtra)
grid.arrange(p1,p2,ncol=1)
