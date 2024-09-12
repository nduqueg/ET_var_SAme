rm(list=ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)


# Dates ----
Dates.ms <- seq(as.Date("1979-02-01"),as.Date("2020-11-30"), by="month")
Dates.a <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month")

seasons <- c("DJF","MAM","JJA","SON")

years.m <- format(Dates.a, format="%Y"); Years <- unique(format(Dates.a, format="%Y")); Years <- Years[-1]
Season.y <- paste0(years.m[-1],"-",time2season(Dates.a)[-length(Dates.a)]); Season.y <- c(Season.y,Season.y[length(Season.y)])
Year.s <- unique(Season.y)
  
## load data -----
P <- brick("./01_data/02_Ppt/MSWEP_Ppt_1979_2020.nc")
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

Col.s.r <- stackApply(P, Season.y, fun=sum); Col.s.r <- Col.s.r *Mask.c
Col.s.t <- rasterToPoints(Col.s.r)[,-c(1,2)]

Pa.seasons <- list()
for (i in seasons) Pa.seasons[[i]] <- as.data.frame(matrix(NA, ncol=nrow(Col.s.t), nrow=length(Years)))
for ( i in seasons) Pa.seasons[[i]] <- apply(t(Col.s.t[, which(substr(Year.s,6,8) == i)]), 2, scale,scale=F)

#### identifying periods of high and low Atl3 --------
Atl3.bool <- list()
Atl3.seasons <- read.csv("./01_data/01_SSTindices/Atl3_std_1980-2020.csv"); row.names(Atl3.seasons) <- Atl3.seasons[,1]; Atl3.seasons <- Atl3.seasons[,-1]
#---------------------------------------------------------------------- identification of periods with SST higher than 0.5 & lower than -0.5
{
  Atl3.bool[["Pos"]] <- data.frame(matrix(FALSE, ncol=4, nrow= length(Years))); colnames(Atl3.bool[["Pos"]]) <- seasons; rownames(Atl3.bool[["Pos"]]) <- Years
  
  Atl3.bool[["Pos"]][which(Atl3.seasons >=1, arr.ind = T)] <- TRUE
  Atl3.bool[["Neg"]] <- data.frame(matrix(FALSE, ncol=4, nrow= length(Years))); colnames(Atl3.bool[["Neg"]]) <- seasons; rownames(Atl3.bool[["Neg"]]) <- Years
  
  Atl3.bool[["Neg"]][which(Atl3.seasons <=-1, arr.ind = T)] <- TRUE

  a <- t(sapply(Atl3.bool,function(x){apply(x,2,sum)}))
  write.csv(a,"./01_data/02_Ppt/Atl3_events_1980-2020.csv")
}

#--------------------------------------------------------------------------- initialization of datasets
{
  P.comp <- list()
  P.comp[["Pos"]] <- P.comp[["Neg"]] <- matrix(NA,nrow=nrow(cord),ncol=4)
  P.comp <- lapply(P.comp, 'colnames<-',seasons)
  Test.comp.MS <- list()
  Test.comp.MS <- lapply(Test.comp.MS, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )
}

Test.comp.f <- function(Phase, neutral){
  
  Result <- rep(NA , nrow(Phase))
  # retriving cells without value and with constant series in which case it cannot be evaluated
  Idx <- which(apply(Phase,2,function(x) sum(!is.na(x)) > 1))
  a <- match(which(apply(Phase,2,function(x) sd(x)==0)), Idx); a <- a[!is.na(a)]; if(length(a!=0)>0) Idx <- Idx[- a]
  b <- match(which(apply(neutral,2,function(x) sd(x)==0)), Idx); b <- b[!is.na(b)]; if(length(b!=0)) Idx <- Idx[-b]
  # needed to use a loop because those are two datasets that cannot change with an apply (unless we built and array and that's more complicated)
  for(j in Idx) Result[j] <- wilcox.test( Phase[,j], neutral[,j ])$p.value
  Result[a] <- -999; if(length(b!=0)) Result[b] <- -999
  return(Result)
}
#------------------------------------------------------------------------------------------------------- identification DIV anomalies in those periods

for ( i in seasons){
  print(i)
  Pos <- Pa.seasons[[i]][ Atl3.bool[["Pos"]][,i]  ,]; Neg <- Pa.seasons[[i]][ Atl3.bool[["Neg"]][,i]  ,]
  Neutral <- Pa.seasons[[i]][ !(Atl3.bool[["Pos"]][,i] | Atl3.bool[["Neg"]][,i]) ,]
  
  P.comp[["Pos"]][,i] <- apply(Pos , 2,mean, na.rm=T); P.comp[["Neg"]][,i] <- apply(Neg, 2,mean, na.rm=T)
  
  Test.comp.MS[["Pos"]][[i]] <- Test.comp.f(Pos, Neutral); Test.comp.MS[["Neg"]][[i]] <- Test.comp.f(Neg, Neutral)
}
Test.comp.MS <- lapply(Test.comp.MS, function(x) do.call(cbind,x))

#### Data transformation for plotting -----------
{
  P.comp <- lapply ( P.comp, cbind, cord)
  Test.comp.MS <- lapply ( Test.comp.MS, cbind, cord)
  
  data.anom <- lapply(P.comp, function(x){
    x <- melt(as.data.frame(x), id=c("x","y"));
    colnames(x) <- c("lon", "lat","Season","Anomaly");
    return(x)})
  
  for(Dir in c("Pos","Neg")) data.anom[[Dir]] <- cbind(data.anom[[Dir]], Dir)
  
}

data.anom2 <- do.call(rbind, data.anom); data.anom2$Dir <- factor(data.anom2$Dir, levels = c("Pos","Neg"))

save(data.anom2,"data.anom2",file="./01_data/02_Ppt/03_Atl3_Comp_Ppt.RData")
save(Test.comp.MS,"Test.comp.MS",file="./01_data/02_Ppt/04_Ppt_Atl3_Comp_Ttest.RData")

#### plotting the map ####
SA <- shapefile("./01_data/hybas_sa_lev03_v1c.shp")

at.m <- c(-200,-150,-100,-75,-50,-25,25,50,75,100,150,200)

ggplot(data.anom2)+facet_grid(Dir~Season, switch = "y")+
  geom_raster(aes(x=lon,y=lat,fill=Anomaly))+
  scale_fill_stepsn(colours=brewer.pal(11,"BrBG"), breaks=at.m,
                    values=c(0:12)/12,limits=c(-200,200),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")))+
  # scale_color_stepsn(colours=brewer.pal(11,"BrBG"),breaks=at.m,values=c(0:12)/12,limits=c(min(at.m),max(at.m)))+ 
  scale_y_continuous(position="right")+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.15)+
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="Ppt Composites |Atl3| >=1*SD (1980-2019)")+
  theme_bw()



