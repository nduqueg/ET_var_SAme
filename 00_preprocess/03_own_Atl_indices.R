rm(list=ls())
cat("\014")

library(raster)
library(ncdf4)
library(hydroTSM)
library(ggplot2)
library(ggrepel)

# Dates ----
Dm.h <- seq(as.Date("1850-01-01"),as.Date("2022-11-30"),by="month")
Dm.e <- seq(as.Date("1854-01-01"),as.Date("2020-12-31"),by="month")
Dm.ana <- seq(as.Date("1885-01-01"),as.Date("2020-12-30"),by="month")

# load data ####
ELI <- read.zoo(read.csv("./01_data/01_SSTindices/00A_ELI_SST.csv")[,-1], index.column = 1)

H.sst <- brick("./00_preprocess//HadSST.4.0.1.0_median.nc")
{
  H.sst.Gmean <- zoo(cellStats(H.sst ,stat="mean",na.rm=T), order.by = Dm.h)
  png(filename = "./01_data/01_SSTindices/01_globalAvg_1850-2022.png", width = 800, height = 600)
  plot(H.sst.Gmean,main="HadSSTv4 global average - based on 1961-1990 clim.",ylab="Temp. Anom. [°C]")
  dev.off()
} #### HadSST - Calculate global mean & plotting

E.sst <- brick("../01_ERSSTv5/Anom_ersst_v5_1854_2020.nc")
{
  E.sst.Gmean <- zoo(cellStats(E.sst ,stat="mean",na.rm=T), order.by = Dm.e)
  png(filename = "./01_data/01_SSTindices/01_globalAvg_1854-2022.png", width = 800, height = 600)
  plot(E.sst.Gmean,main="ERSSTv5 global average - based on 1971-2000 clim.",ylab="Temp. Anom. [°C]")
  dev.off()
} #### ERSST - Calculate global mean & plotting

# join plot
data <- cbind(H.sst.Gmean,E.sst.Gmean)
data.g <- fortify(data,melt=T); levels(data.g$Series) <- c("HadSST v4.0", "ERSST v5")
ggplot(data.g,aes(Index,Value,col=Series))+#facet_wrap(.~Series,ncol=1,scales = "free_y")+
  geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
  ylab("SST Anom. [°C]")+xlab("year") +labs(title="Global Mean")+
  theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

{
  aux <- read.csv("./00_preprocess/co2_mm_mlo.csv",skip=56); CO2.maunaloa <- zoo(aux[,c("average","deseasonalized")], order.by = as.Date(paste0(aux$year,"-",aux$month,"-01")))
  aux <- read.table("./00_preprocess/Meure et al - CO2 concentration",skip=56); aux2 <- aux[match(unique(aux[,2]),aux[,2]),]; for(i in 1:nrow(aux2)){aux2[i,3] <- mean(aux[aux[,2]==aux2[i,2],3])}
  
  CO2.ice <- zoo(aux2[,3], order.by = as.Date(paste0(floor(aux2[,2]),"-0",
                                                    round((aux2[,2]-floor(aux2[,2]))*10+1,0),
                                                    "-01")))
  
  CO2 <- cbind(CO2.maunaloa$deseasonalized,CO2.ice) #plot(CO2,screens=c(1,1),col=c("red","blue"),main="CO2 concentrations merged datasets"); legend("topleft",legend=c("NOAA GML Mauna Loa","Meure et al (2006)"),col=c("red","blue"),lty=c(1,1))
  CO2[is.na(CO2[,1]),1] <- CO2[is.na(CO2[,1]),2]; CO2 <- CO2[,1]
   } #### load and preprocess the CO2


## preprocess -  detrend ####
H.sst <- H.sst[[ which(Dm.h >= Dm.ana[1] & Dm.h <= Dm.ana[length(Dm.ana)] )]]
E.sst <- E.sst[[ which(Dm.e >= Dm.ana[1] & Dm.e <= Dm.ana[length(Dm.ana)] )]]

{# H.sst.a <- H.sst - as.numeric(H.sst.Gmean[Dm.h >= Dm.ana[1] & Dm.h <= Dm.ana[length(Dm.ana)]])
  # E.sst.a <- E.sst - as.numeric(E.sst.Gmean[Dm.e >= Dm.ana[1] & Dm.e <= Dm.ana[length(Dm.ana)]])
  } # could be with linear or could be with the Global average value - COMMENTED OUT

{
  Pairs.ESST.CO2 <- data.frame(SST=numeric(),CO2=numeric())
  Pairs.HSST.CO2 <- data.frame(SST=numeric(),CO2=numeric())
  aux <- data.frame(Cor=numeric(),dataset=character(),stringsAsFactors = F)
  
  Pairs.ESST.CO2[1:length(E.sst.Gmean[index(CO2)]),"SST"] <- E.sst.Gmean[index(CO2)]; Pairs.ESST.CO2[1:length(E.sst.Gmean[index(CO2)]),"CO2"] <- CO2[index(E.sst.Gmean)]
  dataset <- "ERSSTv5"; Pairs.ESST.CO2 <- cbind(Pairs.ESST.CO2,dataset)
  aux[1,"Cor"] <- round(cor(Pairs.ESST.CO2[,1],Pairs.ESST.CO2[,2]),3); aux[1,"dataset"] <- dataset
  
  Pairs.HSST.CO2[1:length(H.sst.Gmean[index(CO2)]),"SST"] <- H.sst.Gmean[index(CO2)]; Pairs.HSST.CO2[1:length(H.sst.Gmean[index(CO2)]),"CO2"] <- CO2[index(H.sst.Gmean)]
  dataset <- "HadSSTv4.0"; Pairs.HSST.CO2 <- cbind(Pairs.HSST.CO2,dataset)
  aux[2,"Cor"] <- round(cor(Pairs.HSST.CO2[,1],Pairs.HSST.CO2[,2]),3); aux[2,"dataset"] <- dataset
  
  
  data <- rbind(Pairs.HSST.CO2, Pairs.ESST.CO2)
  ggplot(data , aes(CO2,SST,col=dataset))+geom_jitter()+
    geom_text_repel(data=aux,aes(x=300,y=0.5,col=dataset,label=Cor))+
    scale_color_brewer(type = "qual",palette="Set1")+
    ylab("SST Anom. Global mean [°C]")+xlab("CO2 [ppm]") +labs(title="Global Mean SST vs CO2 concentration")+
    theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
} # create pairs based on Dates
{lm.Esst <- lm(SST ~ CO2,Pairs.ESST.CO2)
  lm.Hsst <-lm(SST ~ CO2,Pairs.HSST.CO2)
  } # linear model between SST & CO2

{
  D.SST.Hsst <- zoo(NA, order.by = Dm.h)
  D.SST.Esst <- zoo(NA, order.by = Dm.e)
  for ( i in 1:length(Dm.h)){
    f.x <- Dm.h[i]
    f.y <- as.numeric(CO2[which.min(abs(index(CO2) - f.x))])
    D.SST.Hsst[i] <- lm.Hsst$coefficients %*% c(1,f.y)
    } # create trend for HadSST v 4.0
  
  for ( i in 1:length(Dm.e)){
    f.x <- Dm.e[i]
    f.y <- as.numeric(CO2[which.min(abs(index(CO2) - f.x))])
    D.SST.Esst[i] <- lm.Esst$coefficients %*% c(1,f.y)
  } # create trend for ERSST v5
  
  data.Gsst.CO2 <- fortify(cbind(D.SST.Esst,D.SST.Hsst), melt=T)
  levels(data.Gsst.CO2$Series)[levels(data.Gsst.CO2$Series)=="D.SST.Hsst"] <-"HadSST v4.0"
  levels(data.Gsst.CO2$Series)[levels(data.Gsst.CO2$Series)=="D.SST.Esst"] <-"ERSST v5"
  ggplot(data.Gsst.CO2,aes(Index,Value,col=Series))+
    geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
    ylab("SST Anom. Global mean[°C]")+xlab("year") +labs(title="SST Global Mean with CO2 relationship")+
    theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
  
} # ZOO of CO2 values for the dates in the SST datasets & prediction of SST

{
  H.sst.a <- H.sst - as.numeric(D.SST.Hsst[Dm.h >= Dm.ana[1] & Dm.h <= Dm.ana[length(Dm.ana)]])
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
### AMO - not filter ----

H.sst.NAtl <- lat_Weight.Mean( crop(H.sst.a,extent(c(-80,0,0,70))), Dm.ana=Dm.ana)
E.sst.NAtl <- lat_Weight.Mean(crop(E.sst.a,extent(c(360-80,360,0,70))), Dm.ana=Dm.ana)
SST.NAtl <- cbind(H.sst.NAtl,E.sst.NAtl)

write.csv(as.data.frame(SST.NAtl),"./01_data/01_SSTindices/00A_NAtl_SST.csv")

# AMO - filter K=121 - plotting ####
AMO.Atl <- rollmean(SST.NAtl,k=121)
write.csv(as.data.frame(AMO.Atl), "./01_data/01_SSTindices/01A_AMO_1890-2015.csv")

library(ggplot2)

{
  color <- NA
  SST.NAtl.g <- fortify(SST.NAtl,melt=T); Variable <- "SST anom"; SST.NAtl.g <- cbind(SST.NAtl.g, Variable, color)
  SST.NAtl.g$color[which(SST.NAtl.g$Value>0)] <- "Red"; SST.NAtl.g$color[which(SST.NAtl.g$Value<0)] <- "Blue"
  levels(SST.NAtl.g$Series)[levels(SST.NAtl.g$Series)=="H.sst.NAtl"] <- "HadSST v4.0"; levels(SST.NAtl.g$Series)[levels(SST.NAtl.g$Series)=="E.sst.NAtl"] <- "ERSST v5"
  
  AMO.Atl.g <- fortify(AMO.Atl, melt=T);Variable <- "AMO"; AMO.Atl.g <- cbind(AMO.Atl.g, Variable, color)
  AMO.Atl.g$color[which(AMO.Atl.g$Value>0)] <- "Red"; AMO.Atl.g$color[which(AMO.Atl.g$Value<0)] <- "Blue"
  levels(AMO.Atl.g$Series)[levels(AMO.Atl.g$Series)=="H.sst.NAtl"] <- "HadSST v4.0"; levels(AMO.Atl.g$Series)[levels(AMO.Atl.g$Series)=="E.sst.NAtl"] <- "ERSST v5"
  
  data.g <- rbind(SST.NAtl.g, AMO.Atl.g); #data.g$color[which(data.g$Value>0)] <- "Red"; data.g$color[which(data.g$Value<0)] <- "Blue"

}# organize data for plotting

# ggplot()+
#   geom_bar(SST.Atl.g,aes(Index,Value),stat="identity")+
#   geom_bar(AMO.Atl.g,aes(Index,Value,col=color), stat="identity")+
ggplot(data.g,aes(Index,Value,alpha=Variable, col=Series))+#facet_wrap(.~Series,ncol=1,scales = "free_y")+
  geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
  scale_alpha_manual(values=c(0.3,1),guide = 'none')+
  ylab("SST Anom. [°C]")+xlab("year") +labs(title="Atlantic Multidecadal Oscilation",caption = "detrended with CO2 relationship")+
  theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

## Nino3.4 ----
H.sst.N34 <- lat_Weight.Mean( crop(H.sst.a,extent(c(-170,-120,-5,5))), Dm.ana=Dm.ana)
E.sst.N34 <- lat_Weight.Mean(crop(E.sst.a,extent(c(360-170,360-120,-5,5))), Dm.ana=Dm.ana)
SST.N34 <- cbind(H.sst.N34,E.sst.N34)
write.csv(as.data.frame(SST.N34),"./01_data/01_SSTindices/00A_N34_SST.csv")

ONI <- rollmean(SST.N34,k=3)
write.csv(as.data.frame(ONI), "./01_data/01_SSTindices/01A_ONI_1890-2015.csv")

date.comp <- seq(as.Date("1885-01-01"),as.Date("2018-12-31"),by="month")
aux.H.N34 <- H.sst.N34[date.comp]; aux.E.N34 <- E.sst.N34[date.comp]

### PDO ----


# Filter PDO ####


### atlantic EN ----

# SST calculations
H.sst.AEN <- lat_Weight.Mean( crop(H.sst.a,extent(c(-20,0,-3,3))), Dm.ana=Dm.ana)
E.sst.AEN <- lat_Weight.Mean(crop(E.sst.a,extent(c(360-20,360,-3,3))), Dm.ana=Dm.ana)
SST.AEN <- cbind(H.sst.AEN,E.sst.AEN)
write.csv(as.data.frame(SST.AEN),"./01_data/01_SSTindices/00A_AtlEN_SST.csv")

AEN <- rollmean(SST.AEN,k=3)
write.csv(as.data.frame(AEN), "./01_data/01_SSTindices/01A_AEN_1890-2015.csv")

date.comp <- seq(as.Date("1885-01-01"),as.Date("2018-12-31"),by="month")
{cor.ELI.AEN <- numeric(length = 2); names(cor.ELI.AEN) <- c("HadSST v4.0","ERSSTv5")
cor.ELI.AEN[1] <- cor.test(H.sst.AEN[date.comp], ELI$data.m[date.comp])$estimate
cor.ELI.AEN[2] <- cor.test(E.sst.AEN[date.comp], ELI$data.m[date.comp])$estimate

ccf(H.sst.AEN[date.comp], ELI$data.m[date.comp],na.action = na.pass)
ccf(E.sst.AEN[date.comp], ELI$data.m[date.comp],na.action = na.pass)} # all time series correlation

{
  aux.H <- H.sst.AEN[date.comp]; aux.E <- E.sst.AEN[date.comp]; aux.ELI <- ELI$data.m[date.comp]
  cor.ELI.AEN.m <- list();cor.ELI.AEN.m[["HadSST"]] <- numeric();cor.ELI.AEN.m[["ERSST"]] <- numeric()
  date.comp.m <- format(date.comp, format="%m")
  for ( i in unique(date.comp.m)){
    cor.ELI.AEN.m[["HadSST"]][i] <- cor.test(aux.H[date.comp.m==i], aux.ELI[date.comp.m==i])$estimate
    cor.ELI.AEN.m[["ERSST"]][i] <- cor.test(aux.E[date.comp.m==i], aux.ELI[date.comp.m==i])$estimate
  }
  library(reshape2)
  cor.ELI.AEN.m <- data.frame(month=factor(month.abb, levels = month.abb), melt(cor.ELI.AEN.m))
  ggplot(cor.ELI.AEN.m, aes(x=month, y=value,fill=L1, group=L1))+ geom_bar(stat = "identity",position=position_dodge())+
    scale_fill_brewer(type = "qual",palette="Set1",direction = -1)+
    ylab("Correlation")+xlab("Month") +labs(title="Correlation between ELI & Atl3 (1885-2018)")+
    theme_bw()+theme(legend.position=c(0.8,0.2),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
} # correlation by month


#### atlantic EN - plotting ----

{
  color <- NA
  SST.AEN.g <- fortify(SST.AEN, melt =T); Variable <- "SST anom"; SST.AEN.g <- cbind(SST.AEN.g, Variable, color)
  SST.AEN.g$color[which(SST.AEN.g$Value>0)] <- "Red"; SST.AEN.g$color[which(SST.AEN.g$Value<0)] <- "Blue"
  levels(SST.AEN.g$Series)[levels(SST.AEN.g$Series)=="H.sst.AEN"] <- "HadSST v4.0"; levels(SST.AEN.g$Series)[levels(SST.AEN.g$Series)=="E.sst.AEN"] <- "ERSST v5"
  
  AEN.g <- fortify(AEN, melt=T);Variable <- "AEN"; AEN.g <- cbind(AEN.g, Variable, color)
  AEN.g$color[which(AEN.g$Value>0)] <- "Red"; AEN.g$color[which(AEN.g$Value<0)] <- "Blue"
  levels(AEN.g$Series)[levels(AEN.g$Series)=="H.sst.AEN"] <- "HadSST v4.0"; levels(AEN.g$Series)[levels(AEN.g$Series)=="E.sst.AEN"] <- "ERSST v5"
  
  data.g <- rbind(SST.AEN.g, AEN.g);
} # Organize data for plotting

ggplot(data.g,aes(Index,Value,alpha=Variable, col=Series))+#facet_wrap(.~Series,ncol=1,scales = "free_y")+
  geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
  scale_alpha_manual(values=c(0.3,1),guide = 'none')+
  ylab("SST Anom. [°C]")+xlab("year") +labs(title="Atlantic El Niño 3 index",
                                            subtitle=paste(paste("Cor. with ELI =",round(cor.ELI.AEN,2),c("HadSST v4.0","ERSSTv5")),sep="; ",collapse=""),
                                            caption = "detrended with CO2 relationship")+
  theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

{
  data.g <- data.g[which(data.g$Index>=as.Date("1980-01-01") & data.g$Index<=as.Date("2010-12-31")),]
  ggplot(data.g,aes(Index,Value,alpha=Variable, col=Series))+#facet_wrap(.~Series,ncol=1,scales = "free_y")+
    geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
    scale_alpha_manual(values=c(0.3,1),guide = 'none')+
    ylab("SST Anom. [°C]")+xlab("year") +labs(title="Atlantic El Niño 3 index",
                                              subtitle=paste(paste("Cor. with ELI =",round(cor.ELI.AEN,2),c("HadSST v4.0","ERSSTv5")),sep="; ",collapse=""),
                                              caption = "detrended with CO2 relationship")+
    theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
} # plot for the las 30 years


### atlantic meridional mode ----

H.sst.AMM <- lat_Weight.Mean( crop(H.sst.a,extent(c(-70,-15,5,25))), Dm.ana=Dm.ana) - 
  lat_Weight.Mean( crop(H.sst.a,extent(c(-40,0,-25,-5))), Dm.ana=Dm.ana)

E.sst.AMM <- lat_Weight.Mean( crop(E.sst.a,extent(c(360-70,360-15,5,25))), Dm.ana=Dm.ana) - 
  lat_Weight.Mean( crop(E.sst.a,extent(c(360-40,360,-25,-5))), Dm.ana=Dm.ana)

SST.AMM <- cbind(H.sst.AMM,E.sst.AMM)
write.csv(as.data.frame(SST.AMM),"./01_data/01_SSTindices/00A_AMM_deltaSST.csv")

AMM <- rollmean(SST.AMM,k=121, na.rm=T)
write.csv(as.data.frame(AMM), "./01_data/01_SSTindices/01A_AMM_1890-2015.csv")

cor.ELI.AMM<- numeric(length = 2); names(cor.ELI.AMM) <- c("HadSST v4.0","ERSSTv5")
cor.ELI.AMM[1] <- cor.test(H.sst.AMM[date.comp], ELI$data.m[date.comp])$estimate
cor.ELI.AMM[2] <- cor.test(E.sst.AMM[date.comp], ELI$data.m[date.comp])$estimate

ccf(H.sst.AMM[date.comp], ELI$data.m[date.comp],na.action = na.pass)
ccf(E.sst.AMM[date.comp], ELI$data.m[date.comp],na.action = na.pass)


{
  aux.H.AMM <- H.sst.AMM[date.comp]; aux.E.AMM <- E.sst.AMM[date.comp]
  cor.AMM.AEN.m <- list();cor.AMM.AEN.m[["HadSST"]] <- numeric();cor.AMM.AEN.m[["ERSST"]] <- numeric()
  for ( i in unique(date.comp.m)){
    cor.AMM.AEN.m[["HadSST"]][i] <- cor.test(aux.H[date.comp.m==i], aux.H.AMM[date.comp.m==i])$estimate
    cor.AMM.AEN.m[["ERSST"]][i] <- cor.test(aux.E[date.comp.m==i], aux.E.AMM[date.comp.m==i])$estimate
  }
  library(reshape2)
  cor.AMM.AEN.m <- data.frame(month=factor(month.abb, levels = month.abb), melt(cor.AMM.AEN.m))
  ggplot(cor.AMM.AEN.m, aes(x=month, y=value,fill=L1, group=L1))+ geom_bar(stat="identity",position=position_dodge())+
    scale_fill_brewer(type = "qual",palette="Set1",direction = -1)+
    ylab("Correlation")+xlab("Month") +labs(title="Correlation between AMM & Atl3 (1885-2018)")+
    theme_bw()+theme(legend.position=c(0.8,0.2),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
} # correlation by month Atl3


{
  aux.H.AMM <- H.sst.AMM[date.comp]; aux.E.AMM <- E.sst.AMM[date.comp]
  cor.AMM.ELI.m <- list();cor.AMM.ELI.m[["HadSST"]] <- numeric();cor.AMM.ELI.m[["ERSST"]] <- numeric()
  for ( i in unique(date.comp.m)){
    cor.AMM.ELI.m[["HadSST"]][i] <- cor.test(aux.ELI[date.comp.m==i], aux.H.AMM[date.comp.m==i])$estimate
    cor.AMM.ELI.m[["ERSST"]][i] <- cor.test(aux.ELI[date.comp.m==i], aux.E.AMM[date.comp.m==i])$estimate
  }
  library(reshape2)
  cor.AMM.ELI.m <- data.frame(month=factor(month.abb, levels = month.abb), melt(cor.AMM.ELI.m))
  ggplot(cor.AMM.ELI.m, aes(x=month, y=value,fill=L1, group=L1))+ geom_bar(stat="identity",position=position_dodge())+
    scale_fill_brewer(type = "qual",palette="Set1",direction = -1)+
    ylab("Correlation")+xlab("Month") +labs(title="Correlation between AMM & ELI (1885-2018)")+
    theme_bw()+theme(legend.position=c(0.8,0.2),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
} # correlation by month ELI

{
  #aux.H.AMM <- H.sst.AMM[date.comp]; aux.E.AMM <- E.sst.AMM[date.comp]
  cor.AMM.N34.m <- list();cor.AMM.N34.m[["HadSST"]] <- numeric();cor.AMM.N34.m[["ERSST"]] <- numeric()
  for ( i in unique(date.comp.m)){
    cor.AMM.N34.m[["HadSST"]][i] <- cor.test(aux.H.N34[date.comp.m==i], aux.H.AMM[date.comp.m==i])$estimate
    cor.AMM.N34.m[["ERSST"]][i] <- cor.test(aux.E.N34[date.comp.m==i], aux.E.AMM[date.comp.m==i])$estimate
  }
  library(reshape2)
  cor.AMM.N34.m <- data.frame(month=factor(month.abb, levels = month.abb), melt(cor.AMM.N34.m))
  ggplot(cor.AMM.N34.m, aes(x=month, y=value,fill=L1, group=L1))+ geom_bar(stat="identity",position=position_dodge())+
    scale_fill_brewer(type = "qual",palette="Set1",direction = -1)+
    ylab("Correlation")+xlab("Month") +labs(title="Correlation between AMM & N34 (1885-2018)")+
    theme_bw()+theme(legend.position=c(0.8,0.2),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
} # correlation by month N34

### atlantic meridional mode - plotting ----
{
  color <- NA
  SST.AMM.g <- fortify(SST.AMM, melt =T); Variable <- "SST anom"; SST.AMM.g <- cbind(SST.AMM.g, Variable, color)
  SST.AMM.g$color[which(SST.AMM.g$Value>0)] <- "Red"; SST.AMM.g$color[which(SST.AMM.g$Value<0)] <- "Blue"
  levels(SST.AMM.g$Series)[levels(SST.AMM.g$Series)=="H.sst.AMM"] <- "HadSST v4.0"; levels(SST.AMM.g$Series)[levels(SST.AMM.g$Series)=="E.sst.AMM"] <- "ERSST v5"
  
  AMM.g <- fortify(AMM, melt=T);Variable <- "AMM"; AMM.g <- cbind(AMM.g, Variable, color)
  AMM.g$color[which(AMM.g$Value>0)] <- "Red"; AMM.g$color[which(AMM.g$Value<0)] <- "Blue"
  levels(AMM.g$Series)[levels(AMM.g$Series)=="H.sst.AMM"] <- "HadSST v4.0"; levels(AMM.g$Series)[levels(AMM.g$Series)=="E.sst.AMM"] <- "ERSST v5"
  
  data.g <- rbind(SST.AMM.g, AMM.g);
} # Organize data for plotting

ggplot(data.g,aes(Index,Value,alpha=Variable, col=Series))+#facet_wrap(.~Series,ncol=1,scales = "free_y")+
  geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
  scale_alpha_manual(values=c(0.3,1),guide = 'none')+
  ylab("SST Anom. [°C]")+xlab("year") +labs(title="Atlantic Meridional Mode index",caption = "detrended with CO2 relationship")+
  theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

data.g <- data.g[data.g$Index>=as.Date("1975-01-01"),]
ggplot(data.g,aes(Index,Value,alpha=Variable, col=Series))+#facet_wrap(.~Series,ncol=1,scales = "free_y")+
  geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
  scale_alpha_manual(values=c(0.3,1),guide = 'none')+
  ylab("SST Anom. [°C]")+xlab("year") +labs(title="Atlantic Meridional Mode index",caption = "detrended with CO2 relationship")+
  theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

### tropical North atlantic mode ----

H.sst.TNA <- lat_Weight.Mean( crop(H.sst.a,extent(c(-70,-15,5,25))), Dm.ana=Dm.ana)

E.sst.TNA <- lat_Weight.Mean( crop(E.sst.a,extent(c(360-70,360-15,5,25))), Dm.ana=Dm.ana)

SST.TNA <- cbind(H.sst.TNA,E.sst.TNA)
write.csv(as.data.frame(SST.TNA),"./01_data/01_SSTindices/00A_TNA_SST.csv")

D.TNA <- rollmean(SST.TNA,k=121, na.rm=T)
write.csv(as.data.frame(D.TNA), "./01_data/01_SSTindices/01A_decadalTNA_1890-2015.csv")

cor.ELI.TNA<- numeric(length = 2); names(cor.ELI.TNA) <- c("HadSST v4.0","ERSSTv5")
cor.ELI.TNA[1] <- cor.test(H.sst.TNA[date.comp], ELI$data.m[date.comp])$estimate
cor.ELI.TNA[2] <- cor.test(E.sst.TNA[date.comp], ELI$data.m[date.comp])$estimate

ccf(H.sst.TNA[date.comp], ELI$data.m[date.comp],na.action = na.pass)
ccf(E.sst.TNA[date.comp], ELI$data.m[date.comp],na.action = na.pass)

{
  aux.H.TNA <- H.sst.TNA[date.comp]; aux.E.TNA <- E.sst.TNA[date.comp]
  cor.TNA.ELI.m <- list();cor.TNA.ELI.m[["HadSST"]] <- numeric();cor.TNA.ELI.m[["ERSST"]] <- numeric()
  for ( i in unique(date.comp.m)){
    cor.TNA.ELI.m[["HadSST"]][i] <- cor.test(aux.ELI[date.comp.m==i], aux.H.TNA[date.comp.m==i])$estimate
    cor.TNA.ELI.m[["ERSST"]][i] <- cor.test(aux.ELI[date.comp.m==i], aux.E.TNA[date.comp.m==i])$estimate
  }
  library(reshape2)
  cor.TNA.ELI.m <- data.frame(month=factor(month.abb, levels = month.abb), melt(cor.TNA.ELI.m))
  ggplot(cor.TNA.ELI.m, aes(x=month, y=value,fill=L1, group=L1))+ geom_bar(stat="identity",position=position_dodge())+
    scale_fill_brewer(type = "qual",palette="Set1",direction = -1)+
    ylab("Correlation")+xlab("Month") +labs(title="Correlation between TNA & ELI (1885-2018)")+
    theme_bw()+theme(legend.position=c(0.8,0.2),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
} # correlation by month ELI

{
  #aux.H.AMM <- H.sst.AMM[date.comp]; aux.E.AMM <- E.sst.AMM[date.comp]
  cor.TNA.N34.m <- list();cor.TNA.N34.m[["HadSST"]] <- numeric();cor.TNA.N34.m[["ERSST"]] <- numeric()
  for ( i in unique(date.comp.m)){
    cor.TNA.N34.m[["HadSST"]][i] <- cor.test(aux.H.N34[date.comp.m==i], aux.H.TNA[date.comp.m==i])$estimate
    cor.TNA.N34.m[["ERSST"]][i] <- cor.test(aux.E.N34[date.comp.m==i], aux.E.TNA[date.comp.m==i])$estimate
  }
  library(reshape2)
  cor.TNA.N34.m <- data.frame(month=factor(month.abb, levels = month.abb), melt(cor.TNA.N34.m))
  ggplot(cor.TNA.N34.m, aes(x=month, y=value,fill=L1, group=L1))+ geom_bar(stat="identity",position=position_dodge())+
    scale_fill_brewer(type = "qual",palette="Set1",direction = -1)+
    ylab("Correlation")+xlab("Month") +labs(title="Correlation between TNA & N34 (1885-2018)")+
    theme_bw()+theme(legend.position=c(0.8,0.8),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
} # correlation by month N34

### tropical North atlantic  mode - plotting ----
{
  color <- NA
  SST.TNA.g <- fortify(SST.TNA, melt =T); Variable <- "SST anom"; SST.TNA.g <- cbind(SST.TNA.g, Variable, color)
  SST.TNA.g$color[which(SST.TNA.g$Value>0)] <- "Red"; SST.TNA.g$color[which(SST.TNA.g$Value<0)] <- "Blue"
  levels(SST.TNA.g$Series)[levels(SST.TNA.g$Series)=="H.sst.TNA"] <- "HadSST v4.0"; levels(SST.TNA.g$Series)[levels(SST.TNA.g$Series)=="E.sst.TNA"] <- "ERSST v5"
  
  D.TNA.g <- fortify(D.TNA, melt=T);Variable <- "D.TNA"; D.TNA.g <- cbind(D.TNA.g, Variable, color)
  D.TNA.g$color[which(D.TNA.g$Value>0)] <- "Red"; D.TNA.g$color[which(D.TNA.g$Value<0)] <- "Blue"
  levels(D.TNA.g$Series)[levels(D.TNA.g$Series)=="H.sst.TNA"] <- "HadSST v4.0"; levels(D.TNA.g$Series)[levels(D.TNA.g$Series)=="E.sst.TNA"] <- "ERSST v5"
  
  data.g <- rbind(SST.TNA.g, D.TNA.g);
} # Organize data for plotting

ggplot(data.g,aes(Index,Value,alpha=Variable, col=Series))+#facet_wrap(.~Series,ncol=1,scales = "free_y")+
  geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
  scale_alpha_manual(values=c(0.3,1),guide = 'none')+
  ylab("SST Anom. [°C]")+xlab("year") +labs(title="Tropical North Atlantic index",caption = "detrended with CO2 relationship")+
  theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

### tropical South atlantic mode ----

H.sst.TSA <- lat_Weight.Mean( crop(H.sst.a,extent(c(-40,0,-25,-5))), Dm.ana=Dm.ana)

E.sst.TSA <- lat_Weight.Mean( crop(E.sst.a,extent(c(360-40,360,-25,-5))), Dm.ana=Dm.ana)

SST.TSA <- cbind(H.sst.TSA,E.sst.TSA)
write.csv(as.data.frame(SST.TSA),"./01_data/01_SSTindices/00A_TSA_SST.csv")

D.TSA <- rollmean(SST.TSA,k=121, na.rm=T)
write.csv(as.data.frame(D.TSA), "./01_data/01_SSTindices/01A_decadalTSA_1890-2015.csv")

cor.ELI.TSA <- numeric(length = 2); names(cor.ELI.TSA) <- c("HadSST v4.0","ERSSTv5")
cor.ELI.TSA[1] <- cor.test(H.sst.TSA[date.comp], ELI$data.m[date.comp])$estimate
cor.ELI.TSA[2] <- cor.test(E.sst.TSA[date.comp], ELI$data.m[date.comp])$estimate

ccf(H.sst.TSA[date.comp], ELI$data.m[date.comp],na.action = na.pass)
ccf(E.sst.TSA[date.comp], ELI$data.m[date.comp],na.action = na.pass)

### tropical South atlantic  mode - plotting ----
{
  color <- NA
  SST.TSA.g <- fortify(SST.TSA, melt =T); Variable <- "SST anom"; SST.TSA.g <- cbind(SST.TSA.g, Variable, color)
  SST.TSA.g$color[which(SST.TSA.g$Value>0)] <- "Red"; SST.TSA.g$color[which(SST.TSA.g$Value<0)] <- "Blue"
  levels(SST.TSA.g$Series)[levels(SST.TSA.g$Series)=="H.sst.TSA"] <- "HadSST v4.0"; levels(SST.TSA.g$Series)[levels(SST.TSA.g$Series)=="E.sst.TSA"] <- "ERSST v5"
  
  D.TSA.g <- fortify(D.TSA, melt=T);Variable <- "D.TSA"; D.TSA.g <- cbind(D.TSA.g, Variable, color)
  D.TSA.g$color[which(D.TSA.g$Value>0)] <- "Red"; D.TSA.g$color[which(D.TSA.g$Value<0)] <- "Blue"
  levels(D.TSA.g$Series)[levels(D.TSA.g$Series)=="H.sst.TSA"] <- "HadSST v4.0"; levels(D.TSA.g$Series)[levels(D.TSA.g$Series)=="E.sst.TSA"] <- "ERSST v5"
  
  data.g <- rbind(SST.TSA.g, D.TSA.g);
} # Organize data for plotting

ggplot(data.g,aes(Index,Value,alpha=Variable, col=Series))+#facet_wrap(.~Series,ncol=1,scales = "free_y")+
  geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
  scale_alpha_manual(values=c(0.3,1),guide = 'none')+
  ylab("SST Anom. [°C]")+xlab("year") +labs(title="Tropical South Atlantic index",caption = "detrended with CO2 relationship")+
  theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))


### Both Atlantic subtropical modes - plotting ----
aux1 <- SST.TNA.g[SST.TNA.g$Index>=as.Date("1975-01-01"),]; aux2 <- SST.TSA.g[SST.TSA.g$Index>=as.Date("1975-01-01"),]
aux1$Variable <- "TNA";aux2$Variable <- "TSA";
data.g <- rbind(aux1,aux2)

data.g <- data.g[data.g$Series=="ERSST v5",]
ggplot(data.g, aes(Index,Value,col=Variable))+ facet_wrap(.~Series,ncol=1)+
  geom_line()+ scale_color_brewer(type = "qual",palette="Dark2",direction = -1)+
  ylab("SST Anom. [°C]")+xlab("year") +labs(title="Sub-Tropical Atlantic indices",caption = "detrended with CO2 relationship")+
  theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))


### Carebbean mode ----
H.sst.Carb <- lat_Weight.Mean( crop(H.sst.a,extent(c(-85,-65,10,25))), Dm.ana=Dm.ana)

E.sst.Carb <- lat_Weight.Mean( crop(E.sst.a,extent(c(360-85,360-65,10,25))), Dm.ana=Dm.ana)

SST.Carb <- cbind(H.sst.Carb,E.sst.Carb)
write.csv(as.data.frame(SST.Carb),"./01_data/01_SSTindices/00A_Caribbean_SST.csv")

D.Carb <- rollmean(SST.Carb,k=121, na.rm=T)
write.csv(as.data.frame(D.Carb), "./01_data/01_SSTindices/01A_decadalCaribbean_1890-2015.csv")

cor.ELI.Carb <- numeric(length = 2); names(cor.ELI.Carb) <- c("HadSST v4.0","ERSSTv5")
cor.ELI.Carb[1] <- cor.test(H.sst.Carb[date.comp], ELI$data.m[date.comp])$estimate
cor.ELI.Carb[2] <- cor.test(E.sst.Carb[date.comp], ELI$data.m[date.comp])$estimate

ccf(H.sst.Carb[date.comp], ELI$data.m[date.comp],na.action = na.pass)
ccf(E.sst.Carb[date.comp], ELI$data.m[date.comp],na.action = na.pass)

### Carebbean  mode - plotting ----
{
  color <- NA
  SST.Carb.g <- fortify(SST.Carb, melt =T); Variable <- "SST anom"; SST.Carb.g <- cbind(SST.Carb.g, Variable, color)
  SST.Carb.g$color[which(SST.Carb.g$Value>0)] <- "Red"; SST.Carb.g$color[which(SST.Carb.g$Value<0)] <- "Blue"
  levels(SST.Carb.g$Series)[levels(SST.Carb.g$Series)=="H.sst.Carb"] <- "HadSST v4.0"; levels(SST.Carb.g$Series)[levels(SST.Carb.g$Series)=="E.sst.Carb"] <- "ERSST v5"
  
  D.Carb.g <- fortify(D.Carb, melt=T);Variable <- "D.Carb"; D.Carb.g <- cbind(D.Carb.g, Variable, color)
  D.Carb.g$color[which(D.Carb.g$Value>0)] <- "Red"; D.Carb.g$color[which(D.Carb.g$Value<0)] <- "Blue"
  levels(D.Carb.g$Series)[levels(D.Carb.g$Series)=="H.sst.Carb"] <- "HadSST v4.0"; levels(D.Carb.g$Series)[levels(D.Carb.g$Series)=="E.sst.Carb"] <- "ERSST v5"
  
  data.g <- rbind(SST.Carb.g, D.Carb.g);
} # Organize data for plotting

ggplot(data.g,aes(Index,Value,alpha=Variable, col=Series))+#facet_wrap(.~Series,ncol=1,scales = "free_y")+
  geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
  scale_alpha_manual(values=c(0.3,1),guide = 'none')+
  ylab("SST Anom. [°C]")+xlab("year") +labs(title="Caribbean index",caption = "detrended with CO2 relationship")+
  theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

#######################-
## Caribbean and TNA plotting -
#######################-

data.g$ocean <- "Caribbean"; data.TNA.g$ocean <- "TNA"
data.g$Variable <- factor(data.g$Variable, levels=c("SST anom","Decadal"))
data.TNA.g$Variable <- factor(data.TNA.g$Variable, levels=c("SST anom","Decadal"))
data2 <- rbind(data.g,data.TNA.g)
data2 <- subset(data2, Index>="1970-01-01")
ggplot(data2,aes(Index,Value,alpha=Variable, col=Series,linetype=ocean))+#facet_wrap(.~ocean,ncol=1,scales = "free_y")+
  geom_line()+ scale_color_brewer(type = "qual",palette="Set1")+
  scale_alpha_manual(values=c(0.3,1),guide = 'none')+
  ylab("SST Anom. [°C]")+xlab("year") +labs(title="Caribbean and TNA index",caption = "detrended with CO2 relationship")+
  theme_bw()+theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))
