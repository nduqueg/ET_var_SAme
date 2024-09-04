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
# MDiv <- brick("../../../01_DataSets/11_VIMFDiv/ERA5_VIMFdiv_1950-2020.nc")
MDiv <- brick("../../../01_DataSets/11_VIMFDiv/ERA5_Mdiv_1950-2020.nc")
MDiv <- aggregate(MDiv,2)
MDiv.t <- rasterToPoints(MDiv); cord <- MDiv.t[,c(1,2)]
Mask.c<- MDiv[[1]]*0+1
MDiv <- MDiv[[match(Dates.a, Dates.e5)]]

# ---------------------------------SST indices
# ELI.m <- read.zoo(read.csv("../../01_DataSets/01_SST/02_Indices/ELI_m_std_1854_2019.csv")[,-1], index.column = 1)
# Dates.ELI <- index(ELI.m)
# ELI.s <- read.csv("../../01_DataSets/01_SST/02_Indices/ELI_s_1854_2019.csv")
# ELI.s <- ELI.s[match(Year.s,ELI.s$X),]

AMM <- read.zoo(read.csv("../../../01_DataSets/01_SST/04_Cal_Indices/00A_AMM_deltaSST.csv")[,c(1,3)], 
                index.column = 1); AMM <- AMM[Dates.a]

# Atl3 <- read.zoo(read.csv("../../01_DataSets/01_SST/04_Cal_Indices/00A_AtlEN_SST.csv")[,c(1,3)], 
#                  index.column = 1); Atl3 <- Atl3[Dates.a]
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

MDiv.z <- zoo(t(MDiv.t), order.by = Dates.a)

#### calculate seasonal DIV ####

Col.s.r <- stackApply(MDiv, Season.y, fun=sum); Col.s.r <- Col.s.r *Mask.c
Col.s.t <- rasterToPoints(Col.s.r)[,-c(1,2)]

DIVa.seasons <- list()
as.data.frame(matrix(NA, ncol=nrow(Col.s.t), nrow=length(Years))) -> DIVa.seasons[["DJF"]] -> DIVa.seasons[["MAM"]] -> DIVa.seasons[["JJA"]] -> DIVa.seasons[["SON"]]

for ( i in seasons) DIVa.seasons[[i]] <- apply(t(Col.s.t[, which(substr(Year.s,6,8) == i)]), 2, scale,scale=F)

#### identifying periods of high and low AMM --------
AMM.bool <- list()
AMM.seasons <- read.csv("../AMM_std_1980-2020.csv"); row.names(AMM.seasons) <- AMM.seasons[,1]; AMM.seasons <- AMM.seasons[,-1]
#---------------------------------------------------------------------- identification of periods with SST higher than 1*SD & lower than -1*SD
{
  AMM.bool[["Pos"]] <- matrix(FALSE, nrow = length(Years),ncol=4); colnames(AMM.bool[["Pos"]]) <- seasons
  AMM.bool[["Pos"]][which(AMM.seasons >=1.0, arr.ind = T)] <- TRUE
  # AMM.bool[["Pos"]] <- data.frame(matrix(NA, ncol=4, nrow= length(Years))); colnames(AMM.bool[["Pos"]]) <- seasons; rownames(AMM.bool[["Pos"]]) <- Years
  # AMM.bool[["Pos"]][1:length(Years),1:4] <- FALSE
  # AMM.bool[["Pos"]][which(AMM.seasons >=0.3, arr.ind = T)] <- TRUE
  # AMM.bool[["Neg"]] <- data.frame(matrix(NA, ncol=4, nrow= length(Years))); colnames(AMM.bool[["Neg"]]) <- seasons; rownames(AMM.bool[["Neg"]]) <- Years
  # AMM.bool[["Neg"]][1:length(Years),1:4] <- FALSE
  # AMM.bool[["Neg"]][which(AMM.seasons <=-0.3, arr.ind = T)] <- TRUE
  AMM.bool[["Neg"]] <- matrix(FALSE, nrow = length(Years),ncol=4); colnames(AMM.bool[["Neg"]]) <- seasons
  AMM.bool[["Neg"]][which(AMM.seasons <=-1.0, arr.ind = T)] <- TRUE
  
  a <- t(sapply(AMM.bool,function(x){apply(x,2,sum)}))
  write.csv(a,"AMM_events_1980-2020.csv")
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
  Pos <- DIVa.seasons[[i]][ AMM.bool[["Pos"]][,i]  ,]; Neg <- DIVa.seasons[[i]][ AMM.bool[["Neg"]][,i]  ,]
  Neutral <- DIVa.seasons[[i]][ !(AMM.bool[["Pos"]][,i] | AMM.bool[["Neg"]][,i]) ,]
  
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

save(data.anom2,"data.anom2",file="02_AMM_divVIM.RData")
save(Test.comp.E5,"Test.comp.E5",file="03_divVIM_AMM_Comp_Ttest.RData")


#### plotting the map ####
SA <- shapefile("../../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")



# data.anom2$Anomaly <- data.anom2$Anomaly*1000

at.m <- c(-15,-12,-9,-6,-3,3,6,9,12,15)
# at.m <- round(seq(-10.5,10.5,length.out = 12),3)

ggplot(data.anom2)+
  geom_contour(aes(x=lon,y=lat,z=Anomaly,colour=after_stat(level)),binwidth =3)+facet_grid(Dir~Season, switch = "y")+
  scale_color_stepsn(colours=rev(brewer.pal(9,"BrBG")),breaks=at.m,values=c(0:12)/12,limits=c(min(at.m),max(at.m)))+ 
  scale_y_continuous(position="right")+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.05)+
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="VIMFDiv Composites AMM >=1*SD & AMM <=-1*SD (1980-2020)", caption = "0.5°")+
  theme_bw()




##### Positive plus Negative phase ----
data.asymetry <- DIV.comp[["Pos"]]
data.asymetry[,3:6] <- data.asymetry[,3:6] + DIV.comp[["Neg"]][,3:6]

data.asymetry2 <- melt(as.data.frame(data.asymetry), id=c("x","y"))

data.asymetry2$value <- data.asymetry2$value*1000
  
ggplot(data.asymetry2)+
  geom_contour(aes(x=x,y=y,z=value,colour=after_stat(level)))+facet_wrap(.~variable, ncol=4)+
  scale_color_stepsn(colours=rev(brewer.pal(11,"BrBG")),breaks=at.m,values=c(0:12)/12,limits=c(min(at.m),max(at.m)),
                     guide=guide_colorsteps(barwidth=unit(20,"cm")))+ 
  scale_y_continuous(position="right")+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.05)+
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="VIMFDiv asymetry in composites AMM")+
  theme_bw()+theme(legend.position = "bottom")

##### Positive less Negative phase ----
data.enhanced <- DIV.comp[["Pos"]]
data.enhanced[,3:6] <- data.enhanced[,3:6] - DIV.comp[["Neg"]][,3:6]

data.enhanced2 <- melt(as.data.frame(data.enhanced), id=c("x","y"))

data.enhanced2$value <- data.enhanced2$value*1000
data.class <- data.enhanced2$value>=0; data.enhanced2 <- cbind(data.enhanced2,data.class)

ggplot()+facet_wrap(.~variable, ncol=4)+
  geom_contour(data=dplyr::filter(data.enhanced2,data.class==T),aes(x=x,y=y,z=value),linetype=2,color="red",binwidth =0.1)+
  geom_contour(data=dplyr::filter(data.enhanced2,data.class==F),aes(x=x,y=y,z=value),linetype=1,color="blue",binwidth =0.1)+
  # geom_contour(aes(x=x,y=y,z=value,colour=after_stat(level)))+facet_wrap(.~variable, ncol=4)+
  # scale_color_stepsn(colours=rev(brewer.pal(11,"BrBG")),breaks=at.m,values=c(0:12)/12,limits=c(min(at.m),max(at.m)),
  #                    guide=guide_colorsteps(barwidth=unit(20,"cm")))+ 
  scale_y_continuous(position="right")+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.05)+
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="VIMFDiv enhanced in composites AMM")+
  theme_bw()+theme(legend.position = "bottom")
