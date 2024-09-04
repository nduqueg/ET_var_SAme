rm(list=ls())
cat("\014")

library(raster)
library(ncdf4)
library(RColorBrewer)
library(hydroTSM)

#### loading, croping & stacking the files ####
dates.m <- seq(as.Date("1980-01-01"),as.Date("2020-12-31"),by="month")

month.Multi <- format(dates.m,format="%m")
years <- format(dates.m,format="%Y"); Years <- unique(years)
Season.y <- paste0(years[-1],"-",time2season(dates.m)[-length(dates.m)])
Year.s <- unique(Season.y)

# ext <- extent(-82,-35,-20,15) # Rectangule between -20 to 15°N and -82 to -35°W
cdo.cmd <- "cdo -setgrid,regrid E_1980-2022_GLEAM_v3.8a_MO.nc E_1980-2022_GLEAM_v3.8a_MO_.nc"; system(cdo.cmd)

cdo.cmd <- paste0("cdo -sellonlatbox,-82,-34,-20,15 E_1980-2022_GLEAM_v3.8a_MO_.nc GLEAM_tropSAm.nc"); system(cdo.cmd)

cdo.cmd <- paste0("cdo -seldate,1980-01-01,2020-12-31 GLEAM_tropSAm.nc GLEAM_E_1980-2020.nc"); system(cdo.cmd)
file.remove("GLEAM_tropSAm.nc")

cdo.cmd <- paste0("cdo -seassum GLEAM_E_1980-2020.nc GLEAM_E_1980-2020_seasonal.nc"); system(cdo.cmd)

# E.aux <- brick("E_1980-2021_GLEAM_v3.6a_MO.nc")[[1:length(dates.m)]]
# E <- crop(E.aux, ext)

# names(E@z) <- "Date"
# E@z <- list(Date=as.character(dates.m))
# writeRaster(E,filename="GLEAM_E_1980_2020.nc",format="CDF",overwrite=T)

E <- brick("GLEAM_E_1980_2020.nc")
mask_SA <- E[[1]]*0+1
#### accummulating #####
#system("cdo ymonmean '20CR_TotEvap_1836_2015.nc' '20CR_TotEvap_MonMulti.nc'")
# Problems with croping the file with R, it erases the time dimenssion of the raster file
# TotEvap.mM <- brick("20CR_TotEvap_MonMulti.nc")
SA <- shapefile("../../South_America/South_America.shp")
arrow <- list("SpatialPolygonsRescale", layout.north.arrow(), 
              offset = c(-50,5), scale = 5,which=1, first=F)
countries.layer <- list("sp.lines", Colom, col = "Black", lwd=0.5, first=F)

month.Multi <- factor(format(dates.m,format="%m"))
years <- factor(format(dates.m,format="%Y"))

E.mM <- stackApply(E, month.Multi, fun=mean)
writeRaster(E.mM,filename = "GLEAM_E_MonMulti.nc",format="CDF",overwrite=T)


E.y <- stackApply(E, years, fun=sum, na.rm = T); E.y <- E.y *mask_SA
writeRaster(E.y,filename = "GLEAM_E_yearly.nc",format="CDF",overwrite=T)
E.yM <- stackApply(E.y, rep(1,nlayers(E.y)),fun=mean)

###### plotting ####
#at.m <- seq(min(cellStats(P.mM,stat="min")), max(cellStats(P.mM,stat = "max")), length.out = 12)
#at.m <- round(at.m,2); at.m[1] <- floor(at.m[1]); at.m[12] <- ceiling(at.m[12])
at.m <- c(0,10,20,30,50,75,100,125,150,175,200,250)
spplot(E.mM,scales=list(draw=T),col.regions=brewer.pal(11,"Spectral"),sp.layout=list(arrow,countries.layer),
       at= at.m,names.attr= c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
       colorkey=list(height=1,
                     labels= list( labels =at.m, at=at.m)),
       main="Evaporation monthly long-term mean [GLEAM]")

at.m <- c(0,100,200,300,400,500,725,1000,1500,2000,2500,3000)
spplot(E.yM,scales=list(draw=T),col.regions=brewer.pal(11,"Spectral"),sp.layout=list(arrow,countries.layer),
       at= at.m,names.attr= c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
       colorkey=list(height=1,
                     labels= list( labels =at.m, at=at.m)),
       main="Evaporation annual long-term mean [GLEAM]")


#### seasonal preprocessing ----
E <- E[[-c(1,2,length(dates.m))]]; Season.y <- Season.y[-c(1,2)]; Year.s <- Year.s[-1]
#mask.E <- E[[1]]*0+1
E.s <- stackApply(E, Season.y, fun=sum); E.s <- E.s*mask_SA
writeRaster(E.s,filename = "GLEAM_seasonal_80-MAM_20-SON.nc",format="CDF",overwrite=T)

E.s.t <- rasterToPoints(E.s)[,-c(1,2)]
Coord <- as.data.frame(rasterToPoints(E.s[[1]])[,c(1,2)])

MS.seasons <- list()
as.data.frame(matrix(NA, ncol=nrow(E.s.t), nrow=length(Years))) -> MS.seasons[["DJF"]] -> MS.seasons[["MAM"]] -> MS.seasons[["JJA"]] -> MS.seasons[["SON"]]
MS.seasons[["DJF"]][2:length(Years),] <- t(E.s.t[, which(substr(Year.s,6,8) == "DJF")])
MS.seasons[["MAM"]][1:length(Years),] <- t(E.s.t[, which(substr(Year.s,6,8) == "MAM")])
MS.seasons[["JJA"]][1:length(Years),] <- t(E.s.t[, which(substr(Year.s,6,8) == "JJA")])
MS.seasons[["SON"]][1:length(Years),] <- t(E.s.t[, which(substr(Year.s,6,8) == "SON")])

save(MS.seasons,"MS.seasons",file="GLEAM_seasonal_79-MAM_20-SON.Rdata")

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
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="Mean seasonal ET GLEAM (1980-2020)")+
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
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="Seasonal Coef. Var. ET GLEAM (1980-2020)")+
  theme_bw()

library(gridExtra)
grid.arrange(p1,p2,ncol=1)
