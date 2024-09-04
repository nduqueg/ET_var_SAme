rm(list=ls())
cat("\014")

library(raster)
library(ncdf4)
library(RColorBrewer)
library(hydroTSM)

#### loading, croping & stacking the files ####
dates.m <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"),by="month")
d.e5 <-  seq(as.Date("1950-01-01"),as.Date("2020-12-30"),by="month")

month.Multi <- format(dates.m,format="%m")
years <- format(dates.m,format="%Y"); Years <- unique(years); Years <- Years[-1]
Season.y <- paste0(years[-1],"-",time2season(dates.m)[-length(dates.m)]); Season.y <- c(Season.y, Season.y[length(Season.y)])
Year.s <- unique(Season.y)

E <- brick("ERA5L_TEvap_1950_2020.nc")
E <- E[[match(dates.m, d.e5)]]
mask_SA <- E[[1]]*0+1

SA <- shapefile("../../South_America/South_America.shp")


#### seasonal preprocessing ----

E.s <- stackApply(E, Season.y, fun=sum); E.s <- E.s*mask_SA
writeRaster(E.s,filename = "ERA5L_seasonal_80-DJF_20-SON.nc",format="CDF",overwrite=T)

E.s.t <- rasterToPoints(E.s)[,-c(1,2)]
Coord <- as.data.frame(rasterToPoints(E.s[[1]])[,c(1,2)])

MS.seasons <- list()
as.data.frame(matrix(NA, ncol=nrow(E.s.t), nrow=length(Years))) -> MS.seasons[["DJF"]] -> MS.seasons[["MAM"]] -> MS.seasons[["JJA"]] -> MS.seasons[["SON"]]
MS.seasons[["DJF"]][1:length(Years),] <- t(E.s.t[, which(substr(Year.s,6,8) == "DJF")])
MS.seasons[["MAM"]][1:length(Years),] <- t(E.s.t[, which(substr(Year.s,6,8) == "MAM")])
MS.seasons[["JJA"]][1:length(Years),] <- t(E.s.t[, which(substr(Year.s,6,8) == "JJA")])
MS.seasons[["SON"]][1:length(Years),] <- t(E.s.t[, which(substr(Year.s,6,8) == "SON")])

save(MS.seasons,"MS.seasons",file="ERA5L_seasonal_79-MAM_20-SON.Rdata")

Season <- c("DJF","MAM","JJA","SON")
Mean.seasons <- list()
Mean.seasons[["DJF"]] <- cbind(Coord,apply(MS.seasons[["DJF"]],2, mean, na.rm=T),Season[1])
Mean.seasons[["MAM"]] <- cbind(Coord,apply(MS.seasons[["MAM"]],2, mean, na.rm=T),Season[2])
Mean.seasons[["JJA"]] <- cbind(Coord,apply(MS.seasons[["JJA"]],2, mean, na.rm=T),Season[3])
Mean.seasons[["SON"]] <- cbind(Coord,apply(MS.seasons[["SON"]],2, mean, na.rm=T),Season[4])
Mean.seasons <- lapply(Mean.seasons,`colnames<-`,c("x","y","P","Season"))

CV.seasons <- list()
CV.seasons[["DJF"]] <- apply(MS.seasons[["DJF"]],2, function(x){cv <- sd(x,na.rm = T)/mean(x,na.rm=T); return(cv)})
CV.seasons[["MAM"]] <- apply(MS.seasons[["MAM"]],2, function(x){cv <- sd(x,na.rm = T)/mean(x,na.rm=T); return(cv)})
CV.seasons[["JJA"]] <- apply(MS.seasons[["JJA"]],2, function(x){cv <- sd(x,na.rm = T)/mean(x,na.rm=T); return(cv)})
CV.seasons[["SON"]] <- apply(MS.seasons[["SON"]],2, function(x){cv <- sd(x,na.rm = T)/mean(x,na.rm=T); return(cv)})
for(i in 1:4){CV.seasons[[i]] <- cbind(Coord,CV.seasons[[i]] , Season[i])}
CV.seasons <- lapply(CV.seasons,`colnames<-`,c("x","y","P","Season"))

Mean.s.t <- do.call(rbind,Mean.seasons)  # as we don't want to destroy the dataframes that are inside the list, we have to paste them rather than MELT THEM
Var <- "Mean"; Mean.s.t <- cbind(Mean.s.t, Var)

CV.s.t <- do.call(rbind,CV.seasons)
Var <- "Coef. Var."; CV.s.t <- cbind(CV.s.t, Var)

data <- rbind(Mean.s.t, CV.s.t)

library(maps);library(mapdata); library(ggplot2)
Basins <- shapefile("../../South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")
at.m <- c(0,25,50,75,100,150,200,300,400,500,725,1000)
p1<-ggplot(Mean.s.t)+
  geom_raster(aes(x=x,y=y,fill=P))+facet_wrap(.~Season,ncol=4)+
  scale_fill_stepsn(colours=brewer.pal(11,"RdYlGn"), breaks=at.m,
                    values=at.m/max(at.m),limits=c(0,max(at.m)),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")))+
  #geom_sf(data=world,fill="transparent")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype = "dashed", colour="black",fill="NA",size=0.05)+
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="Mean seasonal ET ERA5L (1980-2020)")+
  theme_bw()


at.m <- c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.75,1,2,3,4)
#CV.s.t$P <- log(CV.s.t$P)
p2<-ggplot(CV.s.t)+
  geom_raster(aes(x=x,y=y,fill=P))+facet_wrap(.~Season,ncol=4)+
  scale_fill_stepsn(colours=rev(brewer.pal(11,"Spectral")), breaks=at.m,
                    values=at.m/max(at.m),limits=c(0,4),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")))+
  #geom_sf(data=world,fill="transparent")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype = "dashed", colour="black",fill="NA",size=0.05)+
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="Seasonal Coef. Var. ET ERA5L (1980-2020)")+
  theme_bw()

library(gridExtra)
grid.arrange(p1,p2,ncol=1)
