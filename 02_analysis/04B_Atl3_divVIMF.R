rm(list=ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)


# Dates ----
Dates.e5 <- seq(as.Date("1950-01-01"),as.Date("2020-12-31"), by="month")
Dates.a <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month")

seasons <- c("DJF","MAM","JJA","SON")

years.m <- format(Dates.a, format="%Y"); Years <- unique(format(Dates.a, format="%Y")); Years <- Years[-1]
Season.y <- paste0(years.m[-1],"-",time2season(Dates.a)[-length(Dates.a)]); Season.y <- c(Season.y,Season.y[length(Season.y)])
Year.s <- unique(Season.y)
  
## load data -----
MDiv <- brick("./01_data/09_MDiv/ERA5_Mdiv_1950-2020.nc")
MDiv <- aggregate(MDiv,2)
MDiv.t <- rasterToPoints(MDiv); cord <- MDiv.t[,c(1,2)]
Mask.c<- MDiv[[1]]*0+1
MDiv <- MDiv[[match(Dates.a, Dates.ms)]]

# ---------------------------------SST indices
# ELI.m <- read.zoo(read.csv("../../01_DataSets/01_SST/02_Indices/ELI_m_std_1854_2019.csv")[,-1], index.column = 1)
# Dates.ELI <- index(ELI.m)
# ELI.s <- read.csv("../../01_DataSets/01_SST/02_Indices/ELI_s_1854_2019.csv")
# ELI.s <- ELI.s[match(Year.s,ELI.s$X),]

# AMM <- read.zoo(read.csv("../../../01_DataSets/01_SST/04_Cal_Indices/00A_AMM_deltaSST.csv")[,c(1,3)], 
#                 index.column = 1); AMM <- AMM[Dates.a]

Atl3 <- read.zoo(read.csv("./01_data/01_SSTindices/00A_AtlEN_SST.csv")[,c(1,2)], 
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

Col.s.r <- stackApply(MDiv, Season.y, fun=sum); Col.s.r <- Col.s.r *Mask.c
Col.s.t <- rasterToPoints(Col.s.r)[,-c(1,2)]

DIVa.seasons <- list()
as.data.frame(matrix(NA, ncol=nrow(Col.s.t), nrow=length(Years))) -> DIVa.seasons[["DJF"]] -> DIVa.seasons[["MAM"]] -> DIVa.seasons[["JJA"]] -> DIVa.seasons[["SON"]]

for ( i in seasons) DIVa.seasons[[i]] <- apply(t(Col.s.t[, which(substr(Year.s,6,8) == i)]), 2, scale,scale=F)

#### identifying periods of high and low Atl3 --------
Atl3.bool <- list()
Atl3.seasons <- read.csv("./01_data/01_SSTindices/Atl3_std_1980-2020.csv"); row.names(Atl3.seasons) <- Atl3.seasons[,1]; Atl3.seasons <- Atl3.seasons[,-1]
#---------------------------------------------------------------------- identification of periods with SST higher than 0.4 & lower than -0.4
{
  Atl3.bool[["Pos"]] <- data.frame(matrix(NA, ncol=4, nrow= length(Years))); colnames(Atl3.bool[["Pos"]]) <- seasons; rownames(Atl3.bool[["Pos"]]) <- Years
  Atl3.bool[["Pos"]][1:length(Years),1:4] <- FALSE
  Atl3.bool[["Pos"]][which(Atl3.seasons >=1, arr.ind = T)] <- TRUE
  Atl3.bool[["Neg"]] <- data.frame(matrix(NA, ncol=4, nrow= length(Years))); colnames(Atl3.bool[["Neg"]]) <- seasons; rownames(Atl3.bool[["Neg"]]) <- Years
  Atl3.bool[["Neg"]][1:length(Years),1:4] <- FALSE
  Atl3.bool[["Neg"]][which(Atl3.seasons <=-1, arr.ind = T)] <- TRUE

  a <- t(sapply(Atl3.bool,function(x){apply(x,2,sum)}))
  write.csv(a,"./01_data/09_MDiv/Atl3_events_1980-2020.csv")
}

#--------------------------------------------------------------------------- initialization of datasets
{
  DIV.comp <- list()
  DIV.comp[["Pos"]] <- DIV.comp[["Neg"]] <- matrix(NA,nrow=nrow(cord),ncol=4)
  DIV.comp <- lapply(DIV.comp, 'colnames<-',seasons)
  Test.comp.E5 <- list()
  Test.comp.E5 <- lapply(Test.comp.E5, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )
}

Test.comp.f <- function(Phase, neutral){
  
  Result <- rep(NA , nrow(Phase))
  # retriving cells without value and with constant series in which case it cannot be evaluated
  Idx <- which(apply(Phase,2,function(x) sum(!is.na(x)) > 1))
  a <- match(which(apply(Phase,2,function(x) sd(x)==0)), Idx); a <- a[!is.na(a)]; if(length(a!=0)>0) Idx <- Idx[- a]
  b <- match(which(apply(neutral,2,function(x) sd(x)==0)), Idx); b <- b[!is.na(b)]; if(length(b!=0)) Idx <- Idx[-b]
  # needed to use a loop because those are two datasets that cannot change with an apply (unless we built and array and that's more complicated)
  for(j in Idx) Result[j] <- t.test( Phase[,j], neutral[,j ])$p.value
  Result[a] <- -999; if(length(b!=0)) Result[b] <- -999
  return(Result)
}
#------------------------------------------------------------------------------------------------------- identification DIV anomalies in those periods
for ( i in seasons){
  print(i)
  Pos <- DIVa.seasons[[i]][ Atl3.bool[["Pos"]][,i]  ,]; Neg <- DIVa.seasons[[i]][ Atl3.bool[["Neg"]][,i]  ,]
  Neutral <- DIVa.seasons[[i]][ !(Atl3.bool[["Pos"]][,i] | Atl3.bool[["Neg"]][,i]) ,]
  
  DIV.comp[["Pos"]][,i] <- apply(Pos , 2,mean, na.rm=T); DIV.comp[["Neg"]][,i] <- apply(Neg, 2,mean, na.rm=T)
  
  Test.comp.E5[["Pos"]][[i]] <- Test.comp.f(Pos, Neutral); Test.comp.E5[["Neg"]][[i]] <- Test.comp.f(Neg, Neutral)
}

Test.comp.E5 <- lapply(Test.comp.E5, function(x) do.call(cbind,x))
#### Data transformation for plotting -----------
{
  DIV.comp <- lapply ( DIV.comp, cbind, cord)
  Test.comp.E5 <- lapply ( Test.comp.E5, cbind, cord)
  
  data.anom <- lapply(DIV.comp, function(x){
    x <- melt(as.data.frame(x), id=c("x","y"));
    colnames(x) <- c("lon", "lat","Season","Anomaly");
    return(x)})
  
  for(Dir in c("Pos","Neg")) data.anom[[Dir]] <- cbind(data.anom[[Dir]], Dir)
  
}

data.anom2 <- do.call(rbind, data.anom); data.anom2$Dir <- factor(data.anom2$Dir, levels = c("Pos","Neg"))

save(data.anom2,"data.anom2",file="./01_data/09_MDiv/02_Atl3_divVIM.RData")
save(Test.comp.E5,"Test.comp.E5",file="./01_data/09_MDiv/03_divVIM_Atl3_Comp_Ttest.RData")

#### plotting the map ####
SA <- shapefile("./01_data/hybas_sa_lev03_v1c.shp")


# at.m <- c(-4,-3,-2,-1.5,-1,-0.5,0.5,1,1.5,2,3,4)*0.1
# at.m <- c(-10,-5,-2.5,-1,-0.5,-0.25,0.25,0.5,1,2.5,5,10)*0.01
at.m <- round(seq(-10.5,10.5,length.out = 12),3)

ggplot(data.anom2)+
  geom_contour(aes(x=lon,y=lat,z=Anomaly,colour=after_stat(level)))+ facet_grid(Dir~Season, switch = "y")+
  scale_color_stepsn(colours=rev(brewer.pal(11,"BrBG")),breaks=at.m,values=c(0:12)/12,limits=c(min(at.m),max(at.m)))+ 
  scale_y_continuous(position="right")+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.05)+
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="VIMDiv Composites |Atl3| >=1*SD (1980-2020)")+
  theme_bw()



