rm(list=ls())
cat("\014")

library(raster)
library(ncdf4)
library(hydroTSM)
library(ggplot2)
library(ggrepel)
library(reshape)
"%>%"=magrittr::`%>%`

# Dates ----
Dm.h <- seq(as.Date("1850-01-01"),as.Date("2022-11-30"),by="month")
Dm.e <- seq(as.Date("1854-01-01"),as.Date("2020-12-31"),by="month")
Dm.ana <- seq(as.Date("1885-01-01"),as.Date("2020-12-30"),by="month")

# load data ####
E.sst <- brick("../01_ERSSTv5/Anom_ersst_v5_1854_2020.nc")
E.sst.t <- rasterToPoints(E.sst); E.sst.t[E.sst.t[,1]>=180 ,1] <- E.sst.t[E.sst.t[,1]>=180,1] -360
E.sst <- rasterFromXYZ(E.sst.t)
E.sst.Gmean <- zoo(cellStats(E.sst ,stat="mean",na.rm=T), order.by = Dm.e)

{
  aux <- read.csv("../co2_mm_mlo.csv",skip=56); CO2.maunaloa <- zoo(aux[,c("average","deseasonalized")], order.by = as.Date(paste0(aux$year,"-",aux$month,"-01")))
  aux <- read.table("../Meure et al - CO2 concentration",skip=56); aux2 <- aux[match(unique(aux[,2]),aux[,2]),]; for(i in 1:nrow(aux2)){aux2[i,3] <- mean(aux[aux[,2]==aux2[i,2],3])}
  
  CO2.ice <- zoo(aux2[,3], order.by = as.Date(paste0(floor(aux2[,2]),"-0",
                                                    round((aux2[,2]-floor(aux2[,2]))*10+1,0),
                                                    "-01")))
  
  CO2 <- cbind(CO2.maunaloa$deseasonalized,CO2.ice) #plot(CO2,screens=c(1,1),col=c("red","blue"),main="CO2 concentrations merged datasets"); legend("topleft",legend=c("NOAA GML Mauna Loa","Meure et al (2006)"),col=c("red","blue"),lty=c(1,1))
  CO2[is.na(CO2[,1]),1] <- CO2[is.na(CO2[,1]),2]; CO2 <- CO2[,1]
   } #### load and preprocess the CO2


## preprocess -  detrend ####

E.sst <- E.sst[[ which(Dm.e >= Dm.ana[1] & Dm.e <= Dm.ana[length(Dm.ana)] )]]

{
  Pairs.ESST.CO2 <- data.frame(SST=numeric(),CO2=numeric())
  aux <- data.frame(Cor=numeric(),dataset=character(),stringsAsFactors = F)
  
  Pairs.ESST.CO2[1:length(E.sst.Gmean[index(CO2)]),"SST"] <- E.sst.Gmean[index(CO2)]; Pairs.ESST.CO2[1:length(E.sst.Gmean[index(CO2)]),"CO2"] <- CO2[index(E.sst.Gmean)]
  dataset <- "ERSSTv5"; Pairs.ESST.CO2 <- cbind(Pairs.ESST.CO2,dataset)
  aux[1,"Cor"] <- round(cor(Pairs.ESST.CO2[,1],Pairs.ESST.CO2[,2]),3); aux[1,"dataset"] <- dataset
  

  data <- Pairs.ESST.CO2
  ggplot(data , aes(CO2,SST,col=dataset))+geom_jitter()+
    geom_text_repel(data=aux,aes(x=300,y=0.5,col=dataset,label=Cor))+
    scale_color_brewer(type = "qual",palette="Set1")+
    ylab("SST Anom. Global mean [°C]")+xlab("CO2 [ppm]") +labs(title="Global Mean SST vs CO2 concentration")+
    theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
} # create pairs based on Dates
{lm.Esst <- lm(SST ~ CO2,Pairs.ESST.CO2)

  } # linear model between SST & CO2

{
  D.SST.Esst <- zoo(NA, order.by = Dm.e)
 
  for ( i in 1:length(Dm.e)){
    f.x <- Dm.e[i]
    f.y <- as.numeric(CO2[which.min(abs(index(CO2) - f.x))])
    D.SST.Esst[i] <- lm.Esst$coefficients %*% c(1,f.y)
  } # create trend for ERSST v5
  
  data.Gsst.CO2 <- fortify(D.SST.Esst, melt=T)
  levels(data.Gsst.CO2$Series)[levels(data.Gsst.CO2$Series)=="D.SST.Esst"] <-"ERSST v5"
  ggplot(data.Gsst.CO2,aes(Index,Value,col=Series))+
    geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
    ylab("SST Anom. Global mean[°C]")+xlab("year") +labs(title="SST Global Mean with CO2 relationship")+
    theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
  
} # ZOO of CO2 values for the dates in the SST datasets & prediction of SST

{
  E.sst.a <- E.sst - as.numeric(D.SST.Esst[Dm.e >= Dm.ana[1] & Dm.e <= Dm.ana[length(Dm.ana)]])
} # Detrending with CO2 relationship

# writeRaster(H.sst.a,file="../03_HadISST/HadSSTv4_CO2detrend_1885-.nc",format="CDF")
# writeRaster(E.sst.a,file="../01_ERSSTv5/ERSSTv5_CO2detrend_1885-.nc",format="CDF")

lat_Weight.Mean <- function(R, Dm.ana=Dm.ana){
  aux <- rasterToPoints(R)
  R.weight <- apply(
    aux[,3:length(Dm.ana)]*sqrt(cos(3.14159*aux[,2]/180)) 
    ,2, mean, na.rm=T )
  return( zoo(R.weight, order.by = Dm.ana))
} # Special function for calculated the mean with latitude weighting, also for the two datasets


### Trop Atlantic SSTs ----
Dm.ana2 <- seq(as.Date("1950-01-01"),as.Date("2020-12-30"),by="month")
seasons <- c("DJF","MAM","JJA","SON")
years.m <- format(Dm.ana2, format="%Y"); Years <- unique(format(Dm.ana2, format="%Y"))[-1]
Season.y <- paste0(years.m[-1] ,"-",time2season(Dm.ana2)[-length(Dm.ana2)]); Season.y <- c(Season.y,Season.y[length(Season.y)])
Year.s <- unique(Season.y)


E.sst.a <- E.sst.a[[match(Dm.ana2, Dm.ana)]]

E.sst.Atl <- crop(E.sst.a,extent(c(-75,20,-35,30)))

AMM <- rasterToPoints(E.sst.Atl[[which(Dm.ana2=="1981-04-01")]])

Atl3 <- rasterToPoints(E.sst.Atl[[which(Dm.ana2=="1974-06-01")]])

d.SST <- merge.data.frame(AMM,Atl3); colnames(d.SST)[3:4] <- c("AMM","Atl3"); d.SST <- melt(d.SST,id=c("x","y"))

world <- shapefile("../../South_America/World_Continents.shp")

at.m <- round(seq(-1.5,1.5,length.out=12),1)
rec.df <- data.frame(xmi=c(-70,-40,-20),
                     xma=c(-15,0,0),
                     ymi=c(5,-25,-3),
                     yma=c(25,-5,3),
                     variable=c("AMM","AMM","Atl3"))
ggplot(d.SST,aes(x,y)) + facet_wrap(.~variable)+
  geom_raster(aes(fill=value))+
  scale_fill_stepsn(colours=rev(RColorBrewer::brewer.pal(11,"RdBu")), breaks=at.m,
                    limits=c(min(at.m),max(at.m)),
                    guide=guide_colorsteps(even.steps = F, barwidth=unit(10,"cm")),
                    name="anomaly [°C]")+
  geom_polygon(data=world,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.15)+
  geom_rect(data=rec.df,inherit.aes = F,aes(xmin=xmi, xmax=xma, ymin=ymi, ymax=yma),color="black",linetype="longdash", fill=NA)+
  coord_fixed(xlim=c(-62,20),ylim=c(-30,30))+labs(x="Long",y="Lat",title="Sea Surface Temperature Modes")+
  theme_bw()+theme(legend.position = "bottom")

### PCA ----
SST.t <- rasterToPoints(E.sst.Atl)

cells <- SST.t[,c(1,2)]; SST.t <- t(SST.t[,-c(1,2)])

PC <- prcomp(SST.t, retx=T, scale.=F)
Loadings <- cbind.data.frame(cells, PC$rotation)

# Loadings.r <- stack(rasterFromXYZ(Loadings[,c(1,2,3)]))
# Loadings.r[[2]] <- rasterFromXYZ(Loadings[,c(1,2,4)])

data <- Loadings[,1:4] %>% magrittr::set_colnames(c("lon","lat","1st","2nd")) %>% 
  melt(.,id=c("lon","lat"))

rec.df <- data.frame(xmi=c(-70,-40,-20),
                     xma=c(-15,0,0),
                     ymi=c(5,-25,-3),
                     yma=c(25,-5,3),
                     variable=c("2nd","2nd","1st"))

range <- with(data,{  max( min(value) %>% abs(), max(value))})
at.m <- seq(-1*range,range,length.out=12) %>% round(.,3)
at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]


ggplot(data,aes(lon,lat)) + facet_wrap(.~variable)+
  geom_raster(aes(fill=value))+
  scale_fill_stepsn(colours=rev(RColorBrewer::brewer.pal(11,"RdBu")), breaks=at.m, #values=scales::rescale(at.m.v,from=range(at.m)),
                    limits=c(min(at.m),max(at.m)),
                    guide=guide_colorsteps(even.steps = F, barwidth=unit(10,"cm")),
                    name="anomaly [°C]")+
  geom_polygon(data=world,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.15)+
  geom_rect(data=rec.df,inherit.aes = F,aes(xmin=xmi, xmax=xma, ymin=ymi, ymax=yma),color="black",linetype="longdash", fill=NA)+
  coord_fixed(xlim=c(-70,20),ylim=c(-30,30))+labs(x="Long",y="Lat",title="Sea Surface Temperature Modes")+
  theme_bw()+theme(legend.position = "bottom")

#### seasonal PCA ----
SST.s <- stackApply(E.sst.Atl[[-nlayers(E.sst.Atl)]], Season.y[-length(Season.y)], fun=mean)
SST.s.t <- rasterToPoints(SST.s)[,-c(1,2)]

SST.seasons <- list()
for ( i in seasons) SST.seasons[[i]] <- SST.s.t[, substr(Year.s,6,8) == i ] %>% t() %>% as.data.frame()

PC <- Loadings <- list()
for (i in seasons){
  PC[[i]] <- prcomp(SST.seasons[[i]], retx=T, scale.=F)
  Loadings[[i]] <- cbind.data.frame(cells, PC[[i]]$rotation)
}

data <- lapply(Loadings, 
               function(x) x[,1:4] %>% magrittr::set_colnames(c("lon","lat","1st","2nd")) %>% 
                 melt(.,id=c("lon","lat"))
               ) %>% 
  melt(., id=c("lon","lat","variable","value")) %>% 
  within(., {L1 <- factor(L1, levels = seasons)})

rec.df <- data.frame(xmi=c(-70,-40,-20),
                     xma=c(-15,0,0),
                     ymi=c(5,-25,-3),
                     yma=c(25,-5,3))

ggplot(data,aes(lon,lat)) + facet_grid(variable~ L1, switch="y")+
  geom_raster(aes(fill=value))+
  scale_fill_stepsn(colours=rev(RColorBrewer::brewer.pal(11,"RdBu")), breaks=at.m, #values=scales::rescale(at.m.v,from=range(at.m)),
                    limits=c(min(at.m),max(at.m)),
                    guide=guide_colorsteps(even.steps = F, barwidth=unit(10,"cm")),
                    name="Loadings")+
  geom_polygon(data=world,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.15)+
  geom_rect(data=rec.df,inherit.aes = F,aes(xmin=xmi, xmax=xma, ymin=ymi, ymax=yma),color="black",linetype="longdash", fill=NA)+
  coord_fixed(xlim=c(-70,20),ylim=c(-30,30))+labs(x="Long",y="Lat",title="PCA Sea Surface Temperature Modes")+
  theme_bw()+theme(legend.position = "bottom")
