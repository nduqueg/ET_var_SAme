rm(list=ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)


# Dates ----
Dates.ms <- seq(as.Date("1959-01-01"),as.Date("2020-12-30"), by="month")
Dates.a <- seq(as.Date("1980-01-01"),as.Date("2020-12-31"), by="month")

seasons <- c("DJF","MAM","JJA","SON")

years.m <- format(Dates.a, format="%Y"); Years <- unique(format(Dates.a, format="%Y"))
Season.y <- paste0(years.m[-1],"-",time2season(Dates.a)[-length(Dates.a)])
Year.s <- unique(Season.y)
  
## load data -----
P <- brick("./01_data/07_SLP/ERA5_SLP_1959-2020_crop.nc")
Mask.c<- P[[1]]*0+1
cord <- rasterToPoints(P[[1]])[,c(1,2)] 
sele <- match(Dates.a, Dates.ms)[!is.na(match(Dates.a, Dates.ms))]
P <- P[[sele]]

# ---------------------------------SST indices
# ELI.m <- read.zoo(read.csv("../../01_DataSets/01_SST/02_Indices/ELI_m_std_1854_2019.csv")[,-1], index.column = 1)
# Dates.ELI <- index(ELI.m)
# ELI.s <- read.csv("../../01_DataSets/01_SST/02_Indices/ELI_s_1854_2019.csv")
# ELI.s <- ELI.s[match(Year.s,ELI.s$X),]

# AMM <- read.zoo(read.csv("../../../01_DataSets/01_SST/04_Cal_Indices/00A_AMM_deltaSST.csv")[,c(1,3)], 
#                 index.column = 1); AMM <- AMM[Dates.a]

Atl3 <- read.zoo(read.csv("./01_data/01_SSTindices/00A_AtlEN_SST.csv")[,c(1,3)],
                 index.column = 1); Atl3 <- Atl3[Dates.a]
# 
# TNA <- read.zoo(read.csv("../../01_DataSets/01_SST/04_Cal_Indices/00A_TNA_SST.csv")[,c(1,3)], 
#                 index.column = 1); TNA <- TNA[Dates.a]
# 
# TSA <- read.zoo(read.csv("../../01_DataSets/01_SST/04_Cal_Indices/00A_TSA_SST.csv")[,c(1,3)], 
#                 index.column = 1); TSA <- TSA[Dates.a]
# 
# # ---------------------------------------------------------------------------------------------------------------------------- SLP index
# SLP.iquitos <- read.zoo(read.csv("../../01_DataSets/07A_SLP/03_indexes/01_SLP_iquitos.csv")[,c(1,2)], 
#                          index.column = 1); SLP.iquitos <- SLP.iquitos[Dates.a]
# SLP.PSpain <- read.zoo(read.csv("../../01_DataSets/07A_SLP/03_indexes/01_SLP_PSpain.csv")[,c(1,2)], 
#                       index.column = 1); SLP.PSpain <- SLP.PSpain[Dates.a]
# 
# SLP.i <- SLP.PSpain - SLP.iquitos

#### calculate seasonal DIV ####

Col.s.r <- stackApply(P[[-length(Dates.a)]], Season.y, fun=sum); Col.s.r <- Col.s.r *Mask.c
Col.s.t <- rasterToPoints(Col.s.r)[,-c(1,2)]

Pa.seasons <- list()
as.data.frame(matrix(NA, ncol=nrow(Col.s.t), nrow=length(Years))) -> Pa.seasons[["DJF"]] -> Pa.seasons[["MAM"]] -> Pa.seasons[["JJA"]] -> Pa.seasons[["SON"]]

Pa.seasons[["DJF"]][2:length(Years),] <- apply(t(Col.s.t[, which(substr(Year.s,6,8) == "DJF")[-1]]), 2, scale,scale=F)
Pa.seasons[["MAM"]][1:length(Years),] <- apply(t(Col.s.t[, which(substr(Year.s,6,8) == "MAM")]), 2, scale,scale=F)
Pa.seasons[["JJA"]][1:length(Years),] <- apply(t(Col.s.t[, which(substr(Year.s,6,8) == "JJA")]), 2, scale,scale=F)
Pa.seasons[["SON"]][1:length(Years),] <- apply(t(Col.s.t[, which(substr(Year.s,6,8) == "SON")]), 2, scale,scale=F)

# -------------------------------------------------------------------------------------------------------- Atlantic indices to seasonal time scale
Atl3.s <- aggregate(Atl3[-length(Atl3)], by= list(Season.y), FUN= mean)
# Atl3.s <- aggregate(Atl3[-length(Atl3)], by= list(Season.y), FUN= mean)
# 
# TNA.s <- aggregate(TNA[-length(TNA)], by= list(Season.y), FUN= mean)
# TSA.s <- aggregate(TSA[-length(TSA)], by= list(Season.y), FUN= mean)
# 
# SLPb.s <- aggregate(SLP.i[-length(SLP.i)], by= list(Season.y), FUN= mean)

Atl3.seasons <- data.frame(matrix(NA, ncol=4, nrow= length(Years))); colnames(Atl3.seasons) <- seasons; rownames(Atl3.seasons) <- Years
Atl3.seasons[2:length(Years),"DJF"] <- Atl3.s[which(substr(index(Atl3.s),6,8) == "DJF") [-1]]
Atl3.seasons[   ,"MAM"] <- Atl3.s[which(substr(index(Atl3.s),6,8) == "MAM")]
Atl3.seasons[   ,"JJA"] <- Atl3.s[which(substr(index(Atl3.s),6,8) == "JJA")]
Atl3.seasons[   ,"SON"] <- Atl3.s[which(substr(index(Atl3.s),6,8) == "SON")]

#### identifying periods of high and low Atl3 --------
Atl3.bool <- list()
Atl3.seasons <- read.csv("./01_data/01_SSTindices//Atl3_std_1980-2020.csv"); row.names(Atl3.seasons) <- Atl3.seasons[,1]; Atl3.seasons <- Atl3.seasons[,-1]
#---------------------------------------------------------------------- identification of periods with SST higher than 0.5 & lower than -0.5
{
  Atl3.bool[["Pos"]] <- data.frame(matrix(FALSE, ncol=4, nrow= length(Years))); colnames(Atl3.bool[["Pos"]]) <- seasons; rownames(Atl3.bool[["Pos"]]) <- Years
  
  Atl3.bool[["Pos"]][which(Atl3.seasons >=1, arr.ind = T)] <- TRUE
  Atl3.bool[["Neg"]] <- data.frame(matrix(FALSE, ncol=4, nrow= length(Years))); colnames(Atl3.bool[["Neg"]]) <- seasons; rownames(Atl3.bool[["Neg"]]) <- Years
  
  Atl3.bool[["Neg"]][which(Atl3.seasons <=-1, arr.ind = T)] <- TRUE

  a <- t(sapply(Atl3.bool,function(x){apply(x,2,sum)}))
  write.csv(a,"./01_data/07_SLP/Atl3_events_1980-2020.csv")
}

#--------------------------------------------------------------------------- initialization of datasets
{
  P.comp <- list()
  P.comp[["Pos"]] <- P.comp[["Neg"]] <- matrix(NA,nrow=nrow(cord),ncol=4)
  P.comp <- lapply(P.comp, 'colnames<-',seasons)
}

#------------------------------------------------------------------------------------------------------- identification DIV anomalies in those periods
for ( i in seasons){
  print(i)
  P.comp[["Pos"]][,i] <- apply(Pa.seasons[[i]][ Atl3.bool[["Pos"]][,i]  ,], 2,mean, na.rm=T)
  P.comp[["Neg"]][,i] <- apply(Pa.seasons[[i]][ Atl3.bool[["Neg"]][,i]  ,], 2,mean, na.rm=T)
}
#### Data transformation for plotting -----------
data.anom <- list()
{
  P.comp[["Pos"]] <- cbind(cord, P.comp[["Pos"]])
  P.comp[["Neg"]] <- cbind(cord, P.comp[["Neg"]])
  
  data.anom[["Pos"]] <- melt(as.data.frame(P.comp[["Pos"]]), id=c("x","y"))
  colnames(data.anom[["Pos"]]) <- c("lon", "lat","Season","Anomaly")
  Dir <- "Pos"; data.anom[["Pos"]] <- cbind(data.anom[["Pos"]], Dir)
  
  data.anom[["Neg"]] <- melt(as.data.frame(P.comp[["Neg"]]), id=c("x","y"))
  colnames(data.anom[["Neg"]]) <- c("lon", "lat","Season","Anomaly")
  Dir <- "Neg"; data.anom[["Neg"]] <- cbind(data.anom[["Neg"]], Dir)
}

data.anom2 <- rbind(data.anom[["Pos"]], data.anom[["Neg"]])
save(data.anom2,"data.anom2",file="./01_data/07_SLP/03_Atl3_Comp_Ppt.RData")

#### plotting the map ####
library(maps);library(mapdata)
world <- map_data("world")

range <- max(-min(data.anom2$Anomaly),max(data.anom2$Anomaly))

at.m <- round(seq(-range,range,length.out = 12),1)
data.anom2$lon[data.anom2$lon>=180] <- data.anom2$lon[data.anom2$lon>=180] - 360

ggplot(data.anom2)+facet_grid(Dir~Season, switch = "y")+
  geom_raster(aes(x=lon,y=lat,fill=Anomaly))+
  scale_fill_stepsn(colours=brewer.pal(11,"PuOr"), breaks=at.m,
                    values=c(0:12)/12,limits=c(-range,range),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")))+
  # scale_color_stepsn(colours=brewer.pal(11,"BrBG"),breaks=at.m,values=c(0:12)/12,limits=c(min(at.m),max(at.m)))+ 
  scale_y_continuous(position="right")+
  geom_polygon(data=world,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.15)+
  coord_fixed(xlim=c(-82,0),ylim=c(-30,30))+labs(x="Long",y="Lat",title="SLP Composites |Atl3| >=1*SD (1980-2020)")+
  theme_bw()



