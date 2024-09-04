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

D.e5 <-seq(as.Date("1950-01-01"),as.Date("2020-12-31"),by="month")

# cdo seasmean ERA5_Mdiv_1950-2020.nc ERA5_Mdiv_seasonal.nc


#### loading, croping & stacking the files ####
P <- brick("ERA5_Mdiv_1950-2020.nc")
Colom <- shapefile("../South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")
arrow <- list("SpatialPolygonsRescale", layout.north.arrow(), 
              offset = c(-45,10), scale = 2, first=F)
countries.layer <- list("sp.lines", Colom, col = "Black", first=F)

#### accummulating monthly #####
P <- P[[match(dates.m, D.e5)]]

P.y <- stackApply(P, years, fun=mean, na.rm = T)#; P.y <- P.y *mask_SA
writeRaster(P.y,filename = "ERA5_Fdiv_yearly.nc",format="CDF",overwrite=T)
P.yM <- stackApply(P.y, rep(1,nlayers(P.y)),fun=mean)


P.mM <- stackApply(P, month.Multi, fun=mean)
writeRaster(P.mM,filename = "ERA5_Fdiv_MonMulti.nc",format="CDF",overwrite=T)



###### plotting monthly ####
at.m <- seq(min(cellStats(P.mM,stat="min")), max(cellStats(P.mM,stat = "max")), length.out = 12)
#at.m <- round(at.m,2); at.m[1] <- floor(at.m[1]); at.m[12] <- ceiling(at.m[12])
#at.m <- c(0,50,100,200,300,400,600,800,1000,1200,1500,2000)
spplot(P.mM,scales=list(draw=T),col.regions=brewer.pal(11,"PuOr"),sp.layout=list(arrow,countries.layer),
       at= at.m,
       names.attr= c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
       colorkey=list(height=1,
                     labels= list( labels =round(at.m,2), at=at.m)),
       main="Vert. Integ. Moisture Flux Div. monthly long-term mean [ERA5]")



#at.m <- c(10,100,200,300,400,600,800,1000,1500,2000,3000,10000)
spplot(P.yM,scales=list(draw=T),col.regions=brewer.pal(11,"Spectral"),sp.layout=list(arrow,countries.layer),
       at= at.m,
       colorkey=list(height=1,
                     labels= list( labels =round(at.m,2), at=at.m)),
       main="Vert. Integ. Moisture Flux Div.  Annual long-term mean [ERA5]")

#### seasonal preprocessing ----
#mask.P <- P[[1]]*0+1
P.s <- stackApply(P[[-nlayers(P)]], Season.y, fun=mean)#; P.s <- P.s*mask.P
P.s.t <- rasterToPoints(P.s)[,-c(1,2)]
Coord <- as.data.frame(rasterToPoints(P.s[[1]])[,c(1,2)])

SM.seasons <- list()
as.data.frame(matrix(NA, ncol=nrow(P.s.t), nrow=length(Years))) -> SM.seasons[["DJF"]] -> SM.seasons[["MAM"]] -> SM.seasons[["JJA"]] -> SM.seasons[["SON"]]
SM.seasons[["DJF"]][2:length(Years),] <- t(P.s.t[, which(substr(Year.s,6,8) == "DJF")[-1]])
SM.seasons[["MAM"]][1:length(Years),] <- t(P.s.t[, which(substr(Year.s,6,8) == "MAM")])
SM.seasons[["JJA"]][1:length(Years),] <- t(P.s.t[, which(substr(Year.s,6,8) == "JJA")])
SM.seasons[["SON"]][1:length(Years),] <- t(P.s.t[, which(substr(Year.s,6,8) == "SON")])

Season <- c("DJF","MAM","JJA","SON")
Mean.seasons <- list()
Mean.seasons[["DJF"]] <- cbind(Coord,apply(SM.seasons[["DJF"]],2, mean, na.rm=T),Season[1])
Mean.seasons[["MAM"]] <- cbind(Coord,apply(SM.seasons[["MAM"]],2, mean, na.rm=T),Season[2])
Mean.seasons[["JJA"]] <- cbind(Coord,apply(SM.seasons[["JJA"]],2, mean, na.rm=T),Season[3])
Mean.seasons[["SON"]] <- cbind(Coord,apply(SM.seasons[["SON"]],2, mean, na.rm=T),Season[4])
Mean.seasons <- lapply(Mean.seasons,`colnames<-`,c("x","y","P","Season"))

CV.seasons <- list()
CV.seasons[["DJF"]] <- apply(SM.seasons[["DJF"]],2, function(x){cv <- sd(x,na.rm = T)/mean(x,na.rm=T); return(cv)})
CV.seasons[["MAM"]] <- apply(SM.seasons[["MAM"]],2, function(x){cv <- sd(x,na.rm = T)/mean(x,na.rm=T); return(cv)})
CV.seasons[["JJA"]] <- apply(SM.seasons[["JJA"]],2, function(x){cv <- sd(x,na.rm = T)/mean(x,na.rm=T); return(cv)})
CV.seasons[["SON"]] <- apply(SM.seasons[["SON"]],2, function(x){cv <- sd(x,na.rm = T)/mean(x,na.rm=T); return(cv)})
for(i in 1:4){CV.seasons[[i]] <- cbind(Coord,CV.seasons[[i]] , Season[i])}
CV.seasons <- lapply(CV.seasons,`colnames<-`,c("x","y","P","Season"))

Mean.s.t <- do.call(rbind,Mean.seasons)  # as we don't want to destroy the dataframes that are inside the list, we have to paste them rather than MELT THEM
Var <- "Mean"; Mean.s.t <- cbind(Mean.s.t, Var)

CV.s.t <- do.call(rbind,CV.seasons)
Var <- "Coef. Var."; CV.s.t <- cbind(CV.s.t, Var)

#### seasonal ploting ----

data <- rbind(Mean.s.t, CV.s.t)

library(maps);library(mapdata); library(ggplot2)
world <- map_data("world")
# at.m <- c(0,100,200,500,750,1000,1500,2000,2500,3000)
range <- max(c(abs(min(Mean.s.t$P)),max(Mean.s.t$P)))
at.m <- seq(-range,range,length.out = 12)
at.m[c(6,7)] <- c(-0.01,0.01)

p1<-ggplot(Mean.s.t)+
  geom_raster(aes(x=x,y=y,fill=P))+facet_wrap(.~Season,ncol=4)+
  scale_fill_stepsn(colours=rev(brewer.pal(11,"BrBG")), breaks=at.m,
                    values=c(0:11)/11,limits=c(-range,range),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")))+
  #geom_sf(data=world,fill="transparent")+
  geom_polygon(data=Colom,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="Mean seasonal Mdiv ERA5 (1980-2020)")+
  theme_bw()


at.m <- c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.75,1,1.5,2,3)
#CV.s.t$P <- log(CV.s.t$P)
p2<-
  ggplot(CV.s.t)+
  geom_raster(aes(x=x,y=y,fill=P))+facet_wrap(.~Season,ncol=4)+
  scale_fill_stepsn(colours=rev(brewer.pal(11,"Spectral")), breaks=at.m,
                    values=at.m/max(at.m),limits=c(0,3),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")))+
  #geom_sf(data=world,fill="transparent")+
  geom_polygon(data=Colom,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="Seasonal Coef. Var. Mdiv ERA5-Land (1980-2020)")+
  theme_bw()

library(gridExtra)
grid.arrange(p1,p2,ncol=1)
