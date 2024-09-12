rm(list=ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

Dates.ms <- seq(as.Date("1979-01-01"),as.Date("2020-12-31"), by="month")
Dates.a <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month")

years.m <- format(Dates.a, format="%Y"); Years <- unique(format(Dates.a, format="%Y"))[-1]
Season.y <- paste0(years.m[-1] ,"-",time2season(Dates.a)[-length(Dates.a)]); Season.y <- c(Season.y,Season.y[length(Season.y)])
Year.s <- unique(Season.y)
Seasons <- c("DJF","MAM","JJA","SON")

## load data -----
Col.m <- brick("./01_data/05_Rad/CLARA_CloudFrac_1979_2020.nc")
Col.m.t <- rasterToPoints(Col.m); cord <- Col.m.t[,c(1,2)]; Col.m.t <- Col.m.t[,match(Dates.a, Dates.ms)+2]
cord <- as.data.frame(cord); colnames(cord) <- c("lon","lat")
Mask.c <- Col.m[[1]]*0+1
Col.m <- Col.m[[match(Dates.a, Dates.ms)]]

# ---------------------------------SST indices
ELI.m <- read.zoo(read.csv("./01_data/05_Rad/00A_ELI_SST.csv"), 
                  index.column = 1); ELI.m <- ELI.m[Dates.a]

AMM <- read.zoo(read.csv("./01_data/05_Rad/00A_AMM_deltaSST.csv")[,c(1,3)], 
                index.column = 1); AMM <- AMM[Dates.a]

Atl3 <- read.zoo(read.csv("./01_data/05_Rad/00A_AtlEN_SST.csv")[,c(1,3)], 
                 index.column = 1); Atl3 <- Atl3[Dates.a]



Col.m.z <- zoo(t(Col.m.t), order.by = Dates.a)

#### calculate seasonal Cloud Cover ####

Col.s.r <- stackApply(Col.m, Season.y, fun=mean); Col.s.r <- Col.s.r *Mask.c
Col.s.t <- rasterToPoints(Col.s.r)[,-c(1,2)]

MS.seasons <- list()
as.data.frame(matrix(NA, ncol=nrow(Col.s.t), nrow=length(Years))) -> MS.seasons[["DJF"]] -> MS.seasons[["MAM"]] -> MS.seasons[["JJA"]] -> MS.seasons[["SON"]]

for ( i in Seasons) MS.seasons[[i]][1:length(Years),] <- t(Col.s.t[, which(substr(Year.s,6,8) == i)])

# ------------------------------------- Atlantic indices to seasonal time scale
ELI.s <- aggregate(ELI.m, by= list(Season.y), FUN= mean)
AMM.s <- aggregate(AMM, by= list(Season.y), FUN= mean)
Atl3.s <- aggregate(Atl3, by= list(Season.y), FUN= mean)


#### correlation by season ####
Cor.ELI <- list(); Cor.Atl3 <- list(); Cor.AMM <- list(); Cor.TNA <- list(); Cor.TSA <- list()

for ( i in Seasons){
  print(i)
  Cor.ELI[[i]] <- apply(MS.seasons[[i]] ,2, cor, ELI.s[substr(index(ELI.s),6,8) ==i ,])
  Cor.Atl3[[i]] <- apply(MS.seasons[[i]] ,2, cor, Atl3.s[substr(index(Atl3.s),6,8) ==i ,])
  Cor.AMM[[i]] <- apply(MS.seasons[[i]] ,2, cor, AMM.s[substr(index(AMM.s),6,8) ==i ,])

}
Cor.ELI <- melt(lapply(Cor.ELI, function(x,cord) cbind(cord,as.numeric(x)), cord), id=c("lon", "lat")); index <- "ELI"; Cor.ELI <- cbind(Cor.ELI, index)
Cor.Atl3 <- melt(lapply(Cor.Atl3, function(x,cord) cbind(cord,as.numeric(x)), cord), id=c("lon", "lat")); index <- "Atl3"; Cor.Atl3 <- cbind(Cor.Atl3, index)
Cor.AMM <- melt(lapply(Cor.AMM, function(x,cord) cbind(cord,as.numeric(x)), cord), id=c("lon", "lat")); index <- "AMM"; Cor.AMM <- cbind(Cor.AMM, index)


data <- rbind(Cor.ELI, Cor.Atl3, Cor.AMM); data <- data[,-3]



#### plotting the map ####
Basins <- shapefile("./01_data/hybas_sa_lev03_v1c.shp")
SA <- shapefile("./01_data/South_America.shp")

data$index <- factor(data$index, levels = c("ELI","Atl3","AMM"))
data$L1 <- factor(data$L1, levels = c( "DJF", "MAM", "JJA", "SON"))


at.m <- c(-1,-0.9,-0.75,-0.6,-0.45,-0.3,0.3,0.45,0.6,0.75,0.9,1)

ggplot(data)+
  geom_raster(aes(x=lon,y=lat,fill=value))+facet_grid(index ~ L1, switch = "y")+
  scale_fill_stepsn(colours=brewer.pal(11,"RdGy"), breaks=at.m,
                    values=c(0:12)/12,limits=c(-1,1),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")),
                    name="Pearson R")+
  #geom_sf(data=world,fill="transparent")+
  scale_y_continuous(position="right")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",size=0.05)+
  coord_fixed(xlim=c(-82,-35),ylim=c(-20,15))+labs(x="Long",y="Lat",title="Correlation Coef. Indices vs CLARA - Cloud Fraction (1980-2020)")+
  theme_bw()



