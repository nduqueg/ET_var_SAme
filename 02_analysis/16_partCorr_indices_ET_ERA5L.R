rm(list=ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
"%>%"=magrittr::'%>%'

Dates.e5 <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month")
Dates.a <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month")

seasons <- c("DJF","MAM","JJA","SON")
years.m <- format(Dates.a, format="%Y"); Years <- unique(format(Dates.a, format="%Y"))[-1]
Season.y <- paste0(years.m[-1] ,"-",time2season(Dates.a)[-length(Dates.a)]); Season.y <- c(Season.y,Season.y[length(Season.y)])
Year.s <- unique(Season.y)

## load data -----
# cdo -sellonlatbox,-85,-30,-25,20 197902.nc ../MSWEP_TropSAm/197902.nc
# cdo mergetime 197902.nc 197903.nc ....   202010.nc 202011.nc ../MSWEP_Ppt_1979_2020_cdo.nc
Col.m <- brick("../../03_SpatioTemp/01_TempCharac/03_ET/03_detrended_ERA5L_TEvap_1980_2020.nc")
Col.m.t <- rasterToPoints(Col.m)
cord <- Col.m.t[,c(1,2)] %>%  magrittr::set_colnames(.,c("lon","lat"))
Mask.c<- Col.m[[1]]*0+1
Col.m <- Col.m[[match(Dates.a, Dates.e5)]]

# ---------------------------------SST indices
ELI.m <- read.zoo(read.csv("../../01_DataSets/01_SST/04_Cal_Indices/00A_ELI_SST.csv"), 
                  index.column = 1) %>% .[Dates.a]

AMM <- read.zoo(read.csv("../../01_DataSets/01_SST/04_Cal_Indices/00A_AMM_deltaSST.csv")[,c(1,3)], 
                index.column = 1) %>% .[Dates.a]

Atl3 <- read.zoo(read.csv("../../01_DataSets/01_SST/04_Cal_Indices/00A_AtlEN_SST.csv")[,c(1,3)], 
                 index.column = 1) %>% .[Dates.a]

TNA <- read.zoo(read.csv("../../01_DataSets/01_SST/04_Cal_Indices/00A_TNA_SST.csv")[,c(1,3)], 
                index.column = 1) %>% .[Dates.a]

TSA <- read.zoo(read.csv("../../01_DataSets/01_SST/04_Cal_Indices/00A_TSA_SST.csv")[,c(1,3)], 
                index.column = 1) %>% .[Dates.a]

# ---------------------------------------------------------------------------------------------------------------------------- SLP index
SLP.iquitos <- read.zoo(read.csv("../../01_DataSets/07A_SLP/03_indexes/01_SLP_iquitos.csv")[,c(1,2)], 
                        index.column = 1) %>% .[Dates.a]
SLP.PSpain <- read.zoo(read.csv("../../01_DataSets/07A_SLP/03_indexes/01_SLP_PSpain.csv")[,c(1,2)], 
                       index.column = 1) %>% .[Dates.a]

SLP.i <- SLP.PSpain - SLP.iquitos

#### calculate seasonal Ppt ####
Col.s.r <- stackApply(Col.m, Season.y, fun=sum) %>% magrittr::multiply_by(Mask.c)
Col.s.t <- rasterToPoints(Col.s.r)[,-c(1,2)]

MS.seasons <- list()
for ( i in seasons) MS.seasons[[i]] <- Col.s.t[, substr(Year.s,6,8) == i ] %>% t() %>% as.data.frame()

# ------------------------------------- Atlantic indices to seasonal time scale
ELI.s <- aggregate(ELI.m, by= list(Season.y), FUN= mean)
AMM.s <- aggregate(AMM, by= list(Season.y), FUN= mean)
Atl3.s <- aggregate(Atl3, by= list(Season.y), FUN= mean)

TNA.s <- aggregate(TNA, by= list(Season.y), FUN= mean)
TSA.s <- aggregate(TSA, by= list(Season.y), FUN= mean)

SLPb.s <- aggregate(SLP.i, by= list(Season.y), FUN= mean)

#### partial correlation by season ####
library(ppcor)
Cor.ocean <- list(); Cor.SLP <- list()

Part.Cor <- function(P, ELI.c, AMM.c , Atl3.c){
  Pred.ocean <- data.frame(P = P, ELI = ELI.c, AMM = AMM.c, Atl3 = Atl3.c)
  p.Cor <- pcor(Pred.ocean)$estimate[-1,1]
  
  return(p.Cor)
}
  
for (i in seasons){
  print(i)
  Cor.ocean[[i]] <- apply(MS.seasons[[ i ]] ,2, Part.Cor, 
                          ELI.s[  index(ELI.s) %>% substr(.,6,8) == i ],
                          AMM.s[ index(AMM.s) %>% substr(.,6,8) == i ],
                          Atl3.s[ substr(index(Atl3.s),6,8) == i ])
  print("- SLP")
  Cor.SLP[[i]] <- apply(MS.seasons[[ i ]] ,2, Part.Cor, 
                        SLPb.s[ index(SLPb.s) %>% substr(.,6,8) == i ],
                        ELI.s[  index(ELI.s) %>% substr(.,6,8) == i ],
                        AMM.s[ index(AMM.s) %>% substr(.,6,8) == i ])
}

paste.cord <- function(x,cord){ x %>% t() %>% cbind.data.frame( cord,. )}
Cor.ocean <- lapply(Cor.ocean, paste.cord, cord)
Cor.SLP <- lapply(Cor.SLP, paste.cord, cord)
Cor.SLP <- lapply(Cor.SLP, `colnames<-`,c("lon","lat","SLP","ELI","AMM"))

save(Cor.ocean, Cor.SLP, list=c("Cor.ocean","Cor.SLP"),file="02_partCorr_indices_ET_ERA5L.RData")

#### plotting the map ####
load("02_partCorr_indices_ET_ERA5L.RData")
Basins <- shapefile("../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")
SA <- shapefile("../../01_DataSets/South_America/South_America.shp")

data.ocean <- melt(Cor.ocean, id=c("lon","lat")) %>% magrittr::set_colnames(., c("lon","lat", "Index", "value", "Season"))
data.ocean <- within(data.ocean,{
  Index <- factor(Index, levels = c("ELI","Atl3","AMM"))
  Season <-factor(Season, levels = c("DJF","MAM","JJA","SON"))
})

# color scale
at.m <- c(-1,-0.9,-0.75,-0.6,-0.45,-0.3,0.3,0.45,0.6,0.75,0.9,1)
at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]

ggplot(data.ocean)+ facet_grid(Index ~ Season, switch = "y")+
  geom_raster(aes(x=lon,y=lat,fill=value))+
  scale_fill_stepsn(colours=brewer.pal(11,"PiYG"), breaks=at.m,limits=c(-1,1),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")),
                    name="Pearson R")+
  #geom_sf(data=world,fill="transparent")+
  scale_y_continuous(position="right")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), colour="black", linetype="dashed",fill="NA",linewidth=0.15)+
  coord_fixed(xlim=c(-82,-35),ylim=c(-20,15))+labs(x="Long",y="Lat",title="Partial Correlation Coef. Indices vs ET ERA5-Land (1980-2020)", caption="ERA5-Land")+
  theme_bw()


data.SLP <- melt(Cor.SLP, id=c("lon","lat")); colnames(data.SLP) <- c("lon","lat", "Index", "value", "Season")
data.SLP$Index <- factor(data.SLP$Index, levels = c("SLP","ELI","AMM"))
data.SLP$Season <-factor(data.SLP$Season, levels = c("DJF","MAM","JJA","SON"))

ggplot(data.SLP)+ facet_grid(Index ~ Season, switch = "y")+
  geom_raster(aes(x=lon,y=lat,fill=value))+
  scale_fill_stepsn(colours=brewer.pal(11,"PiYG"), breaks=at.m,limits=c(-1,1),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")),
                    name="Pearson R")+
  #geom_sf(data=world,fill="transparent")+
  scale_y_continuous(position="right")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), colour="black", linetype="dashed",fill="NA",linewidth=0.15)+
  coord_fixed(xlim=c(-82,-35),ylim=c(-20,15))+labs(x="Long",y="Lat",title="Partial Correlation Coef. Indices vs ET ERA5-Land (1980-2020)", caption="ERA5-Land")+
  theme_bw()


#------------------------------ individual modes


data.ocean <- within(data.ocean, 
                     Season <- factor(Season, levels = c("MAM","JJA","SON", "DJF"))
                     )

subset(data.ocean, Index=="ELI") %>% 
  ggplot(.)+ facet_grid(Index ~ Season, switch = "y")+
  geom_raster(aes(x=lon,y=lat,fill=value))+
  scale_fill_stepsn(colours=brewer.pal(11,"PiYG"),breaks=at.m,limits=c(-1,1))+
  #geom_sf(data=world,fill="transparent")+
  scale_y_continuous(position="right")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), colour="black", linetype="dashed",fill="NA",linewidth=0.15)+
  coord_fixed(xlim=c(-82,-35),ylim=c(-20,13))+
  labs(x="Longitude [°]",y="Latitude [°]")+
  theme_bw()+theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(),
                   strip.text = element_text(size=20))

subset(data.ocean, Index=="AMM" & Season !="DJF") %>% 
  ggplot(.)+ facet_grid(Index ~ Season, switch = "y")+
  geom_raster(aes(x=lon,y=lat,fill=value))+
  scale_fill_stepsn(colours=brewer.pal(11,"PiYG"),breaks=at.m,limits=c(-1,1))+
  #geom_sf(data=world,fill="transparent")+
  scale_y_continuous(position="right")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), colour="black", linetype="dashed",fill="NA",linewidth=0.15)+
  coord_fixed(xlim=c(-82,-35),ylim=c(-20,13))+
  labs(x="Longitude [°]",y="Latitude [°]")+
  theme_bw()+theme(legend.position = "none", strip.text.x = element_blank(),
                   strip.text = element_text(size=20))

subset(data.ocean, Index=="Atl3" & Season =="JJA") %>% 
  ggplot(.)+ facet_grid(Index ~ Season, switch = "y")+
  geom_raster(aes(x=lon,y=lat,fill=value))+
  scale_fill_stepsn(colours=brewer.pal(11,"PiYG"),breaks=at.m,limits=c(-1,1),
                    guide=guide_colorsteps(even.steps = T, barwidth=unit(10,"cm")),
                    name="Pearson R")+
  #geom_sf(data=world,fill="transparent")+
  scale_y_continuous(position="right")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), colour="black", linetype="dashed",fill="NA",linewidth=0.15)+
  coord_fixed(xlim=c(-82,-35),ylim=c(-20,13))+
  labs(x="Longitude [°]",y="Latitude [°]")+
  theme_bw()+theme(legend.direction = "horizontal", legend.title.position = "top",
                   strip.text = element_text(size=20), legend.text = element_text(size=10), legend.title = element_text(size=20))
