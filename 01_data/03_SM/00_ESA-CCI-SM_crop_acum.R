rm(list=ls())
cat("\014")

library(reshape)
library(raster)
library(ncdf4)
library(RColorBrewer)
library(hydroTSM)

meanThresh <- function(x, na.rm =0.5) { na_n <- sum(is.na(x)); if((na_n > 0) & (na_n / length(x) > na.rm)) { return(NA)} else { return( mean(x, na.rm = TRUE) ) }}
#### loading, croping & stacking the files ####
Years <- seq(1979,2020)
Dates.d <- seq(as.Date("1979-01-01"),as.Date("2020-12-31"),by="day")
Dates <- seq(as.Date("1979-01-01"),as.Date("2020-12-31"),by="month")
ext <- extent(-85,-30,-25,15) # Rectangule between -20 to 15°N and -82 to -40°W
# cord <- read.csv("../../00_grids/ESA_CCI_SM_cord.csv")

# SM.m <- list()
# # reading, croping & accumulate montly
# for(j in Years){
#   files <- list.files(path=paste0("./",j,"/"),pattern=".nc")
#   SM <- list(length=length(files))
#   
#   SM.aux <- brick(paste0("./",j,"/",files[1]), varname="sm")
#   SM[[1]] <- crop(SM.aux, ext)
#   
#   pb <- txtProgressBar(min = 0,max=length(files),style = 3)
#   for ( i in 2:length(files)){
#   
#     setTxtProgressBar(pb,i)
#     SM.aux <- brick(paste0("./",j,"/",files[i]),varname="sm")
#     SM[[i]] <- crop(SM.aux, ext)
#   
#   }
#   close(pb)
#   SM <- stack(SM)
#   
#   writeRaster(SM,filename=paste0("ESA_CCI_SM_",j,"_daily.nc"),format="CDF",overwrite=T)
#   SM.t <- rasterToPoints(SM); cord <- SM.t[,1:2]; SM.t <- SM.t[,-c(1,2)]
#   SM.t <- zoo(t(SM.t), order.by = Dates.d[ format(Dates.d,format="%Y")== j])
#   SM.t.m <- daily2monthly(SM.t, FUN= meanThresh, na.rm=0.87) # less than 4 values in a month 
#   SM.m[[j-1979+1]] <- rasterFromXYZ(cbind(cord, t(SM.t.m)))
#   writeRaster(SM.m[[j-1979+1]],filename=paste0("ESA_CCI_SM_monthly_",j,".nc"),format="CDF",overwrite=T)
# }
# SM.m <- list()
# ext <- extent(-85,-35,-25,15)
# for( i in Years){
#   SM.m [[i-1979+1]] <- crop(brick(paste0("ESA_CCI_SM_monthly_",i,".nc")), ext)
#   
# }
# SM.m <- stack(SM.m)
# SM.m.t <- rasterToPoints(SM.m)
# bool <- apply(SM.m.t[,-c(1,2)],1,FUN=function(x){Per <- sum(is.na(x))/length(x); return(Per)})
# SM.m.t <- SM.m.t[bool<=0.4,]
# SM.m <- rasterFromXYZ(SM.m.t)
# writeRaster(SM.m, filename = "ESA_CCI_SM_1979_2020_filtered.nc", format="CDF", overwrite=T)


SM.m <- brick("ESA_CCI_SM_1979_2020_filtered.nc")

#### seasonal preprocessing ----
SM <- SM.m[[-c(1:11,504)]]

Dm.ana <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"),by="month")
DmM.ana <- unique(format(Dm.ana, format="%m"))

years.m <- format(Dm.ana, format="%Y"); Years <- unique(format(Dm.ana, format="%Y"))
Season.y <- paste0(years.m[-1],"-",time2season(Dm.ana)[-length(Dm.ana)]); Season.y <- c(Season.y,Season.y[length(Season.y)])
Year.s <- unique(Season.y)

#mask.P <- P[[1]]*0+1
SM.s <- stackApply(SM, Season.y, fun=meanThresh, na.rm=0.5)#; P.s <- P.s*mask.P
writeRaster(SM.s,filename = "ESACCI_SM_seasonal_80-DJF_20-SON.nc",format="CDF",overwrite=T)

SM.s.t <- rasterToPoints(SM.s)
Coord <- as.data.frame(SM.s.t[,1:2]); SM.s.t <- SM.s.t[,-c(1,2)]

MS.seasons <- list()
as.data.frame(matrix(NA, ncol=nrow(SM.s.t), nrow=length(Year.s)/4)) -> MS.seasons[["DJF"]] -> MS.seasons[["MAM"]] -> MS.seasons[["JJA"]] -> MS.seasons[["SON"]]
MS.seasons[["DJF"]] <- t(SM.s.t[, which(substr(Year.s,6,8) == "DJF")])
MS.seasons[["MAM"]] <- t(SM.s.t[, which(substr(Year.s,6,8) == "MAM")])
MS.seasons[["JJA"]] <- t(SM.s.t[, which(substr(Year.s,6,8) == "JJA")])
MS.seasons[["SON"]] <- t(SM.s.t[, which(substr(Year.s,6,8) == "SON")])
MS.seasons <- lapply(MS.seasons, as.data.frame)

save(MS.seasons,"MS.seasons",file="SM_seasonal_80-DJF_20-SON.Rdata")

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

Mean.s.t <- as.data.frame(do.call(rbind,Mean.seasons))  # as we don't want to destroy the dataframes that are inside the list, we have to paste them rather than MELT THEM
Var <- "Mean"; Mean.s.t <- cbind(Mean.s.t, Var)
Mean.s.t$P <- round(Mean.s.t$P,2)

CV.s.t <- as.data.frame(do.call(rbind,CV.seasons))
Var <- "Coef. Var."; CV.s.t <- cbind(CV.s.t, Var)

data <- rbind(Mean.s.t, CV.s.t)

library(maps);library(mapdata); library(ggplot2)
world <- map_data("world")
Colom <- shapefile("../../South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")
at.m <- seq(0,max(Mean.s.t$P),length.out = 11)
p1<-ggplot(Mean.s.t)+
  geom_raster(aes(x=x,y=y,fill=P))+facet_wrap(.~Season,ncol=4)+
  scale_fill_stepsn(colours=brewer.pal(9,"YlGnBu"), breaks=at.m,
                    values=c(0:12)/12,limits=c(0,max(Mean.s.t$P)),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")))+
  #geom_sf(data=world,fill="transparent")+
  geom_polygon(data=Colom,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
  coord_fixed(xlim=c(-82,-30),ylim=c(-25,20))+labs(x="Long",y="Lat",title="Mean seasonal Soil Moisture - ESA CCI (1980-2020)")+
  theme_bw()


at.m <- c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.75,1,2,3,4)
#CV.s.t$P <- log(CV.s.t$P)
p2<-ggplot(CV.s.t)+
  geom_raster(aes(x=x,y=y,fill=P))+facet_wrap(.~Season,ncol=4)+
  scale_fill_stepsn(colours=rev(brewer.pal(11,"Spectral")), breaks=at.m,
                    values=at.m/max(at.m),limits=c(0,4),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")))+
  #geom_sf(data=world,fill="transparent")+
  geom_polygon(data=Colom,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
  coord_fixed(xlim=c(-82,-30),ylim=c(-25,20))+labs(x="Long",y="Lat",title="Seasonal Coef. Var. Soil Moisture - ESA CCI (1980-2020)")+
  theme_bw()

library(gridExtra)
grid.arrange(p1,p2,ncol=1)
