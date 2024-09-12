rm(list=ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)


# Dates ----
Dm.gl <- seq(as.Date("1980-01-01"),as.Date("2020-11-30"), by="month")
Dm.e5 <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month"); Dy.e5 <- format(Dm.e5,format="%Y")
Dates.a <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month")

seasons <- c("DJF","MAM","JJA","SON")

years.m <- format(Dates.a, format="%Y"); Years <- unique(format(Dates.a, format="%Y")); Years <- Years[-1]
Season.y <- paste0(years.m[-1],"-",time2season(Dates.a)[-length(Dates.a)]); Season.y <- c(Season.y, Season.y[length(Season.y)])
Year.s <- unique(Season.y)
  
years.m.e5 <- format(Dm.e5, format="%Y"); Years.e5 <- unique(format(Dm.e5, format="%Y")); Years.e5 <- Years.e5[-1]
Season.y.e5 <- paste0(years.m.e5[-1], "-",time2season(Dm.e5)[-length(Dm.e5)]); Season.y.e5 <- c(Season.y.e5, Season.y.e5[length(Season.y.e5)])
Year.s.e5 <- unique(Season.y.e5)

years.m.gl <- format(Dm.gl, format="%Y"); Years.gl <- unique(format(Dm.gl, format="%Y")); Years.gl <- Years.gl[-1]
Season.y.gl <- paste0(years.m.gl[-1], "-",time2season(Dm.gl)[-length(Dm.gl)])
Year.s.gl <- unique(Season.y.gl)

## load data -----
# ----------------------------------------------------------------------------------------------------------------------- detrended
E.gl <- brick("./01_data/04_Evap/03_detrended_ET_GLEAM_seasonal_1980_2020.nc")
E.gl.t <- rasterToPoints(E.gl); cord.gl <- E.gl.t[,1:2]; E.gl.t <- E.gl.t[,-c(1:2)]

# cdo expr,'e=e*(-1000);' -muldpm [ ERA5L_original.nc ] ERA5L_TEvap_1950_2020.nc
# cdo seassum ERA5L_TEvap_1950_2020.nc ERA5L_TEvap_seasonal.nc
E.e5 <-  brick("./01_data/04_Evap/03_detrended_ET_ERA5L_seasonal_1980_2020.nc")
E.e5.t <- rasterToPoints(E.e5); cord.e5 <- E.e5.t[,1:2]; E.e5.t <- E.e5.t[,-c(1:2)]
E.e5.t <- E.e5.t[,match(Year.s, Year.s.e5)]

# ---------------------------------SST indices
# ELI.m <- read.zoo(read.csv("../../01_DataSets/01_SST/02_Indices/ELI_m_std_1854_2019.csv")[,-1], index.column = 1)
# Dates.ELI <- index(ELI.m)
# ELI.s <- read.csv("../../01_DataSets/01_SST/02_Indices/ELI_s_1854_2019.csv")
# ELI.s <- ELI.s[match(Year.s,ELI.s$X),]

AMM <- read.zoo(read.csv("./01_data/01_SSTindices/AMM_std_1980-2020.csv"),index.column = 1)
Atl3 <- read.zoo(read.csv("./01_data/01_SSTindices/Atl3_std_1980-2020.csv"),index.column = 1)

#### separate seasons ####
E.e5.Mean <- list(); E.e5.s <- list()
E.gl.Mean <- list(); E.gl.s <- list()

for( i in seasons){
  print(paste("seasonal mean - ", i))
  E.e5.s[[i]] <- t(E.e5.t[, substr(Year.s,6,8) == i])
  E.e5.Mean[[i]] <- apply( E.e5.s[[i]], 2, mean)
  
  E.gl.s[[i]] <- t(E.gl.t[, substr(Year.s.gl,6,8) == i])
  E.gl.Mean[[i]] <- apply( E.gl.s[[i]], 2, mean)
}
E.gl.s[["DJF"]] <- rbind(rep(NA,ncol(E.gl.s[["DJF"]])),E.gl.s[["DJF"]])

#### composite of ET in the positive and negative phase of each phenomenon ----
Comp.E5 <- list(); Comp.E5[["AMM"]] <- list(); Comp.E5[["Atl3"]] <- list(); Comp.E5 <- lapply(Comp.E5, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )
Comp.gl <- list(); Comp.gl[["AMM"]] <- list(); Comp.gl[["Atl3"]] <- list(); Comp.gl <- lapply(Comp.gl, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )

Test.comp.E5 <- list(); Test.comp.E5[["AMM"]] <- list(); Test.comp.E5[["Atl3"]] <- list(); Test.comp.E5 <- lapply(Test.comp.E5, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )
Test.comp.gl <- list(); Test.comp.gl[["AMM"]] <- list(); Test.comp.gl[["Atl3"]] <- list(); Test.comp.gl <- lapply(Test.comp.gl, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )

Test.comp.f <- function(Phase, neutral){

  Result <- rep(NA , ncol(Phase))
  # retriving cells without value and with constant series in which case it cannot be evaluated
  Idx <- which(apply(Phase,2,function(x) sum(!is.na(x)) > 1))
  a <- match(which(apply(Phase,2,function(x) sd(x)==0)), Idx); a <- a[!is.na(a)]; if(length(a!=0)>0) Idx <- Idx[- a]
  b <- match(which(apply(neutral,2,function(x) sd(x)==0)), Idx); b <- b[!is.na(b)]; if(length(b!=0)>0) Idx <- Idx[-b]
  # needed to use a loop because those are two datasets that cannot change with an apply (unless we built and array and that's more complicated)
  for(j in Idx) Result[j] <- t.test( Phase[,j], neutral[, j ])$p.value
  if(length(a!=0)) Result[a] <- -999; if(length(b!=0)) Result[b] <- -999
  return(Result)
}

for( i in seasons){
  print(paste(i, "- AMM - ERA5L"))
  Pos <-  E.e5.s [[i]][ AMM[,i]>= 1, ]; Neg <-  E.e5.s [[i]][ AMM[,i]<= -1, ]; Neutral <-  E.e5.s [[i]][ AMM[,i]<= 1 & AMM[,i]>= -1, ]
  Comp.E5[["AMM"]][["Pos"]][[i]] <- apply(Pos ,2, mean, na.rm=T) - E.e5.Mean[[i]]
  Comp.E5[["AMM"]][["Neg"]][[i]] <- apply(Neg ,2, mean, na.rm=T) - E.e5.Mean[[i]]
  Test.comp.E5[["AMM"]][["Pos"]][[i]] <- Test.comp.f(Pos, Neutral); Test.comp.E5[["AMM"]][["Neg"]][[i]] <- Test.comp.f(Neg, Neutral)
  
  print("- Atl3 - ERA5L")
  Pos <-  E.e5.s [[i]][ Atl3[,i]>= 1,]; Neg <-  E.e5.s [[i]][ Atl3[,i]<= -1, ]; Neutral <-  E.e5.s [[i]][ Atl3[,i]<= 1 & Atl3[,i]>= -1, ]
  Comp.E5[["Atl3"]][["Pos"]][[i]] <- apply(Pos ,2, mean, na.rm=T) - E.e5.Mean[[i]]
  Comp.E5[["Atl3"]][["Neg"]][[i]] <- apply(Neg ,2, mean, na.rm=T) - E.e5.Mean[[i]]
  Test.comp.E5[["Atl3"]][["Pos"]][[i]] <- Test.comp.f(Pos, Neutral); Test.comp.E5[["Atl3"]][["Neg"]][[i]] <- Test.comp.f(Neg, Neutral)
  
  print(paste(i, "- AMM - GLEAM"))
  Pos <-  E.gl.s [[i]][ AMM[,i]>= 1, ]; Neg <-  E.gl.s [[i]][ AMM[,i]<= -1, ]; Neutral <-  E.gl.s [[i]][ AMM[,i]<= 1 & AMM[,i]>= -1, ]
  Comp.gl[["AMM"]][["Pos"]][[i]] <- apply(Pos ,2, mean, na.rm=T) - E.gl.Mean[[i]]
  Comp.gl[["AMM"]][["Neg"]][[i]] <- apply(Neg ,2, mean, na.rm=T) - E.gl.Mean[[i]]
  Test.comp.gl[["AMM"]][["Pos"]][[i]] <- Test.comp.f(Pos, Neutral); Test.comp.gl[["AMM"]][["Neg"]][[i]] <- Test.comp.f(Neg, Neutral)
  
  print("- Atl3 - GLEAM")
  Pos <-  E.gl.s [[i]][ Atl3[,i]>= 1, ]; Neg <-  E.gl.s [[i]][ Atl3[,i]<= -1, ]; Neutral <-  E.gl.s [[i]][ Atl3[,i]<= 1 & Atl3[,i]>= -1, ]
  Comp.gl[["Atl3"]][["Pos"]][[i]] <- apply(Pos ,2, mean, na.rm=T) - E.gl.Mean[[i]]
  Comp.gl[["Atl3"]][["Neg"]][[i]] <- apply(Neg ,2, mean, na.rm=T) - E.gl.Mean[[i]]
  Test.comp.gl[["Atl3"]][["Pos"]][[i]] <- Test.comp.f(Pos, Neutral); Test.comp.gl[["Atl3"]][["Neg"]][[i]] <- Test.comp.f(Neg, Neutral)
  
}

Comp.org <- function(x, cord){ x$Pos <- as.data.frame(x$Pos); x$Pos <- cbind(cord,x$Pos); x$Neg <- as.data.frame(x$Neg); x$Neg <- cbind(cord,x$Neg); return(x)}

Comp.E5 <- lapply(Comp.E5, Comp.org, cord.e5)
Comp.gl <- lapply(Comp.gl, Comp.org, cord.gl)

Test.comp.E5 <- lapply(Test.comp.E5, Comp.org, cord.e5)
Test.comp.gl <- lapply(Test.comp.gl, Comp.org, cord.gl)

save(Comp.E5, Comp.gl, list=c("Comp.E5","Comp.gl"), file="./01_data/04_Evap/ET_Composites.RData")
save(Test.comp.E5, Test.comp.gl, list=c("Test.comp.E5", "Test.comp.gl"), file="./01_data/04_Evap/ET_Composites_Ttest.RData")

#### seasonal plotting ####
library(reshape)
library(ggplot2)
library(raster)
library(RColorBrewer)

load("./01_data/04_Evap/ET_Composites.RData")
load("./01_data/04_Evap/ET_Composites_Ttest.RData")
seasons <- c("DJF","MAM","JJA","SON")


# organizing the statistical significance in tables for creating the shapefiles
{
  # to_shape <-function (data.sig, Phase, Mode, dataset){
  #   aux.sig <- list()
  #   for (i in seasons){
  #     aux.t <- data.sig[,c("x","y",i)]; aux.t <- subset(aux.t, !is.na(i))
  #     aux.shp <-  rasterFromXYZ(aux.t , digits=3) <=0.05  # creating the shapefile as a one line dataframe
  #     aux.shp <- sf::st_as_sf(aggregate(rasterToPolygons( aggregate(aux.shp, fact=5, fun="modal"), fun = function(x){x==1})))
  #     aux.shp$Season <- i # attaching the fields to the shapefile
  #     aux.shp$Mode <- Mode; aux.shp$Phase <- Phase; aux.shp$dataset <- dataset;
  #     
  #     aux.sig [[i]] <- aux.shp
  #   }
  #   return(do.call(rbind,aux.sig))
  # }
  # 
  # data.sig <- list()
  # n.modes <- c("AMM","Atl3")
  # dataset <- "ERA5L"
  # data.sig[["ERA5L"]] <- lapply(Test.comp.E5, 
  #                               function(Mode, n.modes){ 
  #                                 n.phase <- names(Mode); mode <- lapply(Mode, function(Phase, n.phase){ 
  #                                   phase <- to_shape(Phase, n.phase, n.mode, dataset); return(phase)}, n.phase); return(do.call(rbind,mode))
  # }, n.modes)
}
data.sig <- list()
data.sig[["ERA5L"]] <- melt(Test.comp.E5, id=c("x","y"))
data.sig[["GLEAM"]] <- melt(Test.comp.gl, id=c("x","y"))

data.sig <- lapply(data.sig, function(x){ y <- x[!is.na(x$value),]; z <- y[y$value!=-999,]; return(z)})
data.sig <- lapply(data.sig,`colnames<-`,c("lon", "lat", "Season", "value", "Phase", "Mode" ))

shape.sig <- list()  # the nested "for"s are necesary because it's not possible to create the shapefiles all together at the same time
for(dataset in c("ERA5L", "GLEAM")){    shape.sig[[dataset]] <- list()
  for(phase in c("Pos","Neg")){   shape.sig[[dataset]][[phase]] <- list()

    for(mode in c("AMM","Atl3")){  shape.sig[[dataset]][[phase]][[mode]] <- list()
      print(paste(dataset,"-",phase,"-",mode))

      for(i in seasons){
        aux.t <- data.sig[[dataset]] 
        aux.t <- subset(aux.t,Phase==phase & Season == i & Mode== mode)# aux.t <- aux.t[aux.t$Phase==Phase & aux.t$Season == i & aux.t$Mode== Mode,]
        aux.shp <- rasterFromXYZ(aux.t[,c("lon","lat","value")], digits=3) <=0.05  # creating the shapefile as a one line dataframe
        if (dataset == "GLEAM") fact.data <- 2 else fact.data <- 5
        aux.shp <- sf::st_as_sf(aggregate(rasterToPolygons( aggregate(aux.shp, fact=fact.data, fun="modal"), fun = function(x){x==1})))
  
        aux.shp$Season <- i # attaching the fields to the shapefile
        aux.shp$Mode <- mode; aux.shp$Phase <- phase; aux.shp$dataset <- dataset;
        
        shape.sig[[dataset]] [[phase]][[mode]] [[i]] <- aux.shp
      }

    shape.sig[[dataset]][[phase]][[mode]] <- do.call(rbind, shape.sig[[dataset]][[phase]][[mode]])  # joining all the lines in the dataframe shapefile
  }

  shape.sig[[dataset]][[phase]] <- do.call(rbind, shape.sig[[dataset]][[phase]])  # joining all the lines in the dataframe shapefile
}

shape.sig[[dataset]] <- do.call(rbind, shape.sig[[dataset]])  # joining all the lines in the dataframe shapefile
}
shape.sig <- do.call(rbind, shape.sig)  # joining all the lines in the dataframe shapefile
shape.sig$Season <- factor(shape.sig$Season, levels = seasons); shape.sig$Phase <- factor(shape.sig$Phase, levels = c("Pos","Neg"))
save(shape.sig,file="./01_data/04_Evap/ET_Composites_Ttest_shape.RData")

load("./01_data/04_Evap/ET_Composites_Ttest_shape.RData")
library(ggpattern)
library(sf)
Basins <- shapefile("./01_data/hybas_sa_lev03_v1c.shp")
SA <- shapefile("./01_data/South_America.shp")

# organizing the composites
data <- melt(Comp.E5, id=c("x","y")); data$dataset <- "ERA5L"; data <- data[!is.nan(data$value),]
aux <- melt(Comp.gl, id=c("x","y")); aux$dataset <- "GLEAM"; aux <- aux[!is.nan(aux$value),]
data <- rbind(data,aux); colnames(data) <- c("lon", "lat", "Season", "value", "Phase", "Mode", "dataset" )
data$Phase <- factor(data$Phase, levels = c("Pos","Neg"))


# plotting
at.m <- c(-80,-50,-25,-15,-10,-5,5,10,15,25,50,80); at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]

plot.sig <- function(data.spec, shape.spec, dataset, lg.pos = "right"){
  mode <- as.character(unique(data.spec$Mode))
  p1 <- ggplot() + facet_grid(Season ~ Phase, switch = "y")+
    geom_raster(data= data.spec, aes(lon, lat, fill=value),na.rm=T)+
    
    scale_fill_stepsn(colours=brewer.pal(11,"PiYG"), breaks=at.m,
                      values=scales::rescale(at.m.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                      guide=guide_colorsteps( barwidth=unit(10,"cm")), name="ET anom. [mm]")+
    
    geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
    geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
    
    geom_sf_pattern(data = shape.spec, pattern= "circle", pattern_density=0.1, pattern_spacing= 0.02, fill="00", colour="00", pattern_colour="blue")+
    
    coord_sf(xlim=c(-82,-35),ylim=c(-20,15))+
    # scale_y_continuous(position="right")+
    
    labs(x="Long",y="Lat",title=paste("Evaporation Composite", mode ,"(1980-2020)"), caption = dataset)+
    theme_bw()
  if(lg.pos != "right") p1 <- p1 + theme(legend.position = lg.pos) else p1 <- p1 + theme(legend.position = "none") 
  
  # print(p1)
  return(p1)
}

AMM.E5 <- plot.sig(data.spec= subset(data, dataset == "ERA5L" & Mode =="AMM" & Season !="DJF"),
                   shape.spec= subset(shape.sig, dataset=="ERA5L" &  Mode=="AMM" & Season !="DJF"),
                   dataset ="ERA5-Land", lg.pos = "none")

AMM.gl <- plot.sig(data.spec= subset(data, dataset == "GLEAM" & Mode =="AMM" & Season !="DJF"),
                    shape.spec= subset(shape.sig, dataset=="GLEAM" &  Mode=="AMM" & Season !="DJF"),
                    dataset ="GLEAM", lg.pos = "none")

Atl3.E5 <- plot.sig(data.spec= subset(data, dataset == "ERA5L" & Mode =="Atl3" & Season =="JJA"),
                    shape.spec= subset(shape.sig, dataset=="ERA5L" &  Mode=="Atl3" & Season =="JJA"),
                    dataset ="ERA5-Land")

Atl3.gl <- plot.sig(data.spec= subset(data, dataset == "GLEAM" & Mode =="Atl3" & Season =="JJA"),
                     shape.spec= subset(shape.sig, dataset=="GLEAM" &  Mode=="Atl3" & Season =="JJA"),
                     dataset ="GLEAM", lg.pos="bottom")

library(gridExtra)
grid.arrange(AMM.E5, Atl3.E5,ncol=2, widths=c(2.05,1))

grid.arrange(AMM.gl, Atl3.gl,ncol=2, widths=c(2.05,1))

grid.arrange(AMM.gl, Atl3.gl,nrow=2, heights=c(2.05,1))

#### asymmetry ----

asym <- function(x){
  cord <- x$Pos[,1:2]
  asy <- x$Pos[,3:6] + x$Neg[,3:6]
  asy <- cbind.data.frame(cord,asy)
  return(asy)
}

Asym.E5 <- lapply(Comp.E5, asym)
data.s <- melt(Asym.E5, id=c("x","y"));  colnames(data.s) <- c("lon", "lat", "Season", "value", "Mode" )


p.1 <- ggplot() + facet_wrap(.~Season, ncol=3)+
  geom_raster(data= subset(data.s, Mode=="AMM" & Season!="DJF"), aes(lon, lat, fill=value),na.rm=T)+
  
  scale_fill_stepsn(colours=brewer.pal(11,"PiYG"), breaks=at.m,
                    values=scales::rescale(at.m.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                    guide=guide_colorsteps( barheight=unit(7,"cm")), name="ET anom. [mm]")+
  
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  
  coord_sf(xlim=c(-82,-35),ylim=c(-20,15))+
  # scale_y_continuous(position="right")+
  
  labs(x="Long",y="Lat",title=paste("Evaporation asymmetry AMM (1980-2020)"))+
  theme_bw()+theme(legend.position = "none")

p.2 <- ggplot() + facet_wrap(.~Season, ncol=3)+
  geom_raster(data= subset(data.s, Mode=="Atl3" & Season=="JJA"), aes(lon, lat, fill=value),na.rm=T)+
  
  scale_fill_stepsn(colours=brewer.pal(11,"PiYG"), breaks=at.m,
                    values=scales::rescale(at.m.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                    guide=guide_colorsteps( barheight=unit(7,"cm")), name="ET anom. [mm]")+
  
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  
  coord_sf(xlim=c(-82,-35),ylim=c(-20,15))+
  # scale_y_continuous(position="right")+
  
  labs(x="Long",y="Lat",title=paste("Atl3"))+
  theme_bw()

grid.arrange(p.1, p.2,ncol=2, widths=c(2.05,1))
