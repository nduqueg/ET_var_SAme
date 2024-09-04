rm(list=ls())
cat("\014")

library(raster)
source("../../../read_brickToPoints.R")
library(hydroTSM)
source("../../../time2season.R")
library(zoo)
library(ncdf4)
library(reshape)

#### Dates ####
Dates <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"),by="month")
Dm.e5 <- seq(as.Date("1950-01-01"),as.Date("2020-12-31"),by="month")
Dm.cci <- seq(as.Date("1979-01-01"),as.Date("2020-12-31"),by="month")

years.m <- format(Dates, format="%Y"); Years <- unique(format(Dates, format="%Y")); Years <- Years[-1]
years.m.e5 <- format(Dm.e5, format="%Y"); Years.e5 <- unique(format(Dm.e5, format="%Y")); Years.e5 <- Years.e5[-1]
years.m.cci <- format(Dm.cci, format="%Y"); Years.cci <- unique(format(Dm.cci, format="%Y")); Years.cci <- Years.cci[-1]

seasons <- c("DJF","MAM","JJA","SON")

Season.y <- paste0(years.m[-1], # eliminate the 1st value of year so as to the Dec of the 1st year is on the next seasonal year
                   "-",time2season(Dates)[-length(Dates)]); Season.y <- c(Season.y,Season.y[length(Season.y)])
Year.s <- unique(Season.y)

Season.y.e5 <- paste0(years.m.e5[-1], "-",time2season(Dm.e5)[-length(Dm.e5)])
Year.s.e5 <- unique(Season.y.e5)
Season.y.cci <- paste0(years.m.cci[-1], "-",time2season(Dm.cci)[-length(Dm.cci)])
Year.s.cci <- unique(Season.y.cci)

##### load data - monthly data #####
cord <- list()

# read SM
# f1 <- nc_open("../../../01_DataSets/03_SM/03_ESA_CCI_SM/ESA_CCI_SM_1979_2020_filtered.nc")
# CCI.m <- ncvar_get(f1, varid = "variable")
CCI.m <- rasterToPoints(brick("../../../01_DataSets/03_SM/03_ESA_CCI_SM/ESA_CCI_SM_1979_2020_filtered.nc"))
# lon.cci <- f1$dim[[which(names(f1$dim)=="longitude" | names(f1$dim)=="lon") ]]$vals; lat.cci <- f1$dim[[which(names(f1$dim)=="latitude" | names(f1$dim)=="lat" )]]$vals
# cord[["CCI"]] <- expand.grid(lon=lon.cci, lat=lat.cci)
# 
# f2 <- nc_open("../../../01_DataSets/03_SM/01_ERA5L/ERA5L_SM_1L_1950-2020.nc")
# # E5.m <- ncvar_get(f2, varid = "swvl1")
# lon <- f2$dim[[which(names(f2$dim)=="longitude" | names(f2$dim)=="lon") ]]$vals; lat <- f2$dim[[which(names(f2$dim)=="latitude" | names(f2$dim)=="lat" )]]$vals
# cord[["E5"]] <- expand.grid(lon=lon, lat=lat)
# 
# # Soils
# fsoil <- nc_open("../../../01_DataSets/01_ERA5L_soil_type.nc")
# E5.soil <- ncvar_get(fsoil, varid="slt")
# lon <- fsoil$dim[[which(names(fsoil$dim)=="longitude" | names(fsoil$dim)=="lon") ]]$vals; lon[lon>=180] <- lon[lon>=180] - 360
# lat <- fsoil$dim[[which(names(fsoil$dim)=="latitude" | names(fsoil$dim)=="lat" )]]$vals
# cord[["E5_soil"]] <- expand.grid(lon=lon, lat=lat)

# trim soils to tropical south america
# E5.soil <- E5.soil[lon >= min(cord[["E5"]]$lon) & lon <= max(cord[["E5"]]$lon),
#                    lat >= min(cord[["E5"]]$lat) & lat <= max(cord[["E5"]]$lat)]; E5.soil <- round(E5.soil,2)
# cord[["E5_soil"]] <- subset(cord[["E5_soil"]], lon >= min(cord[["E5"]]$lon) & lon <= max(cord[["E5"]]$lon) & 
#                               lat >= min(cord[["E5"]]$lat) & lat <= max(cord[["E5"]]$lat)); cord[["E5_soil"]] <- round(cord[["E5_soil"]],2)


# SSTs
AMM <- read.zoo(read.csv("../../03_Composites/AMM_std_1980-2020.csv"),index.column = 1)
Atl3 <- read.zoo(read.csv("../../03_Composites/Atl3_std_1980-2020.csv"),index.column = 1)

#### crop 1980 - onwards####
# CCI.m <- CCI.m[,,match(Dates,Dm.cci)]
# E5.m <- E5.m[,,match(Dates,Dm.e5)]


#### load seasonal preprocessing ####

# cdo seasmean ERA5L_SM_1L_1950-2020.nc ERA5L_SM_1L_seasonal.nc
E5.s.r <- brick("../../../01_DataSets/03_SM/01_ERA5L/ERA5L_SM_1L_seasonal.nc")
soil.r <- shift(brick("../../../01_DataSets/01_ERA5L_soil_type.nc"),dx=-360)
aux <- mask(resample(soil.r, E5.s.r[[1]], method="ngb"),E5.s.r[[1]])
E5.s.r <- mask(resample(E5.s.r,aux),aux)
soil.r <- rasterToPoints(aux); soil.r[,3] <- round(soil.r[,3],0)


E5.s.r <- rasterToPoints(E5.s.r[[match(Year.s,Year.s.e5)]])
cord[["E5"]] <- E5.s.r[,1:2]; E5.s.r <- E5.s.r[,-c(1,2)]

CCI.s.r <- rasterToPoints(brick("../../../01_DataSets/03_SM/03_ESA_CCI_SM/ESACCI_SM_seasonal_80-DJF_20-SON.nc"))
cord[["CCI"]] <- CCI.s.r[,1:2]; CCI.s.r <- CCI.s.r[,-c(1,2)]



# CCI.m <- array(CCI.m,dim=c(length(lon.cci),length(lat.cci),12,length(Dates)/12)); CCI.m <- aperm(CCI.m, c(1,2,4,3))
CCI.s <- list()
# E5.m <- array(E5.m,dim=c(length(lon),length(lat),12,length(Dates)/12)); E5.m <- aperm(E5.m, c(1,2,4,3))
E5.s <- list()
# 
# m.season <- data.frame(DJF=1:3,MAM=4:6,JJA=7:9, SON=10:12) # start with one because the trim of dates was since December
for( i in seasons){
  print(i)
  # CCI.s[[i]] <- apply(CCI.m[,,,m.season[[i]]],c(1,2,3),mean, na.rm=T)
  CCI.s[[i]] <- CCI.s.r[,substr(Year.s, 6,8)==i]

  print("ERA5L")
  # aux <- apply(E5.m[,,,m.season[[i]]],c(1,2,3),mean, na.rm=T)
  E5.s[[i]] <- E5.s.r[,substr(Year.s, 6,8)==i]
}

#### identification of saturation level ####

# CCI.sat <- apply(CCI.m,c(1,2),max, na.rm=T)
# CCI.sat <- cbind(cord[["CCI"]], as.numeric(CCI.sat))
CCI.sat <- apply(CCI.m[,-c(1,2)], 1,max, na.rm=T)
CCI.sat <- cbind(cord[["CCI"]], as.numeric(CCI.sat))
# E5.sat <- apply(E5.m,c(1,2),max, na.rm=T)
# E5.sat <- cbind(cord[["E5"]], as.numeric(E5.sat))

# saturation values from ECMWF documentation (ECMWF, 2023)
soil.features <- data.frame(type=seq(1,7), sat= c(0.403, 0.439, 0.43, 0.52, 0.614, 0.766, 0.472))
E5.sat <- soil.r
for(i in 1:7){
  E5.sat[soil.r[,3]== i,] <- soil.features$sat[i]
}
E5.sat[E5.sat[,3] == 0, ]  <- NA


#### determine Saturation ----
CCI.s.sat <- list(); CCI.s.sat.anom <- list()
E5.s.sat <- list(); E5.s.sat.anom <- list()

for( i in seasons){
  CCI.s.sat[[i]] <- CCI.s[[i]]/ CCI.sat[,3] *100
  E5.s.sat[[i]] <- E5.s[[i]]/ E5.sat[,3] *100
  
}
E5.s.sat.anom <- lapply(E5.s.sat, function(x) t(apply(x, 1, scale, scale=F)))
CCI.s.sat.anom <- lapply(CCI.s.sat, function(x) t(apply(x, 1, scale, scale=F)))

#### composite of SM  in the positive and negative phase of each phenomenon ----
Comp.E5 <- list(); Comp.E5[["AMM"]] <- list(); Comp.E5[["Atl3"]] <- list(); Comp.E5 <- lapply(Comp.E5, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )
Comp.CCI <- list(); Comp.CCI[["AMM"]] <- list(); Comp.CCI[["Atl3"]] <- list(); Comp.CCI <- lapply(Comp.CCI, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )
Comp.E5.anom <- Comp.E5
Comp.CCI.anom <- Comp.CCI

Test.comp.E5 <- list(); Test.comp.E5[["AMM"]] <- list(); Test.comp.E5[["Atl3"]] <- list()
Test.comp.E5 <- lapply(Test.comp.E5, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )

Test.comp.CCI <- list(); Test.comp.CCI[["AMM"]] <- list(); Test.comp.CCI[["Atl3"]] <- list()
Test.comp.CCI <- lapply(Test.comp.CCI, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )

Test.comp.f <- function(Phase, neutral){
  
  Result <- rep(NA , nrow(Phase)) 
  # retriving cells without value and with constant series in which case it cannot be evaluated
  Idx <- which(apply(Phase,1,function(x) sum(!is.na(x)) > 1))
  a <- match(which(apply(Phase,1,function(x) sd(x)==0)), Idx);a <- a[!is.na(a)]; if(length(a!=0)>0) Idx <- Idx[- a]
  b <- match(which(apply(neutral,1,function(x) sd(x)==0)), Idx); b <- b[!is.na(b)]; if(length(b!=0)>0) Idx <- Idx[-b]
  # needed to use a loop because those are two datasets that cannot change with an apply (unless we built and array and that's more complicated)
  for(j in Idx) Result[j] <- t.test( Phase[j,], neutral[j, ])$p.value
  Result[a] <- -999; if(length(b!=0)) Result[b] <- -999
  return(Result)
}
  
for( i in seasons){
  for(j in c("AMM","Atl3")){
    print(paste(i, "-",j,"- ERA5L"))
    Pos <- E5.s.sat[[i]][, get(j)[,i]>= 1]; Neg <- E5.s.sat[[i]][, get(j)[,i]<= -1]; Neutral <- E5.s.sat[[i]][, get(j)[,i]<= 1 & get(j)[,i]>= -1]
    Comp.E5[[j]][["Pos"]][[i]] <- apply(Pos ,1, mean, na.rm=T); Comp.E5[[j]][["Neg"]][[i]] <- apply(Neg,1, mean, na.rm=T)
    Test.comp.E5[[j]][["Pos"]][[i]] <- Test.comp.f(Pos, Neutral); Test.comp.E5[[j]][["Neg"]][[i]] <- Test.comp.f(Neg, Neutral)
    
    Pos <- E5.s.sat.anom[[i]][, get(j)[,i]>= 1]; Neg <- E5.s.sat.anom[[i]][, get(j)[,i]<= -1]; Neutral <- E5.s.sat.anom[[i]][, get(j)[,i]<= 1 & get(j)[,i]>= -1]
    Comp.E5.anom[[j]][["Pos"]][[i]] <- apply(Pos ,1, mean, na.rm=T); Comp.E5.anom[[j]][["Neg"]][[i]] <- apply(Neg,1, mean, na.rm=T)
    
    
    print("- CCI")
    Pos <- CCI.s.sat[[i]][, get(j)[,i]>= 1]; Neg <- CCI.s.sat[[i]][, get(j)[,i]<= -1]; Neutral <- CCI.s.sat[[i]][, get(j)[,i]<= 1 & get(j)[,i]>= -1]
    Comp.CCI[[j]][["Pos"]][[i]] <- apply(Pos ,1, mean, na.rm=T); Comp.CCI[[j]][["Neg"]][[i]] <- apply(Neg,1, mean, na.rm=T)
    Test.comp.CCI[[j]][["Pos"]][[i]] <- Test.comp.f(Pos, Neutral); Test.comp.CCI[[j]][["Neg"]][[i]] <- Test.comp.f(Neg, Neutral)
    
    Pos <- CCI.s.sat.anom[[i]][, get(j)[,i]>= 1]; Neg <- CCI.s.sat.anom[[i]][, get(j)[,i]<= -1]; Neutral <- CCI.s.sat.anom[[i]][, get(j)[,i]<= 1 & get(j)[,i]>= -1]
    Comp.CCI.anom[[j]][["Pos"]][[i]] <- apply(Pos ,1, mean, na.rm=T); Comp.CCI.anom[[j]][["Neg"]][[i]] <- apply(Neg,1, mean, na.rm=T)
  }
  
}

SM.Comp.org <- function(x, cord){ x$Pos <- as.data.frame(x$Pos); x$Pos <- cbind(cord,x$Pos); x$Neg <- as.data.frame(x$Neg); x$Neg <- cbind(cord,x$Neg); return(x)}

Test.comp.E5_2 <- Test.comp.E5
Test.comp.CCI_2 <- Test.comp.CCI

Comp.E5 <- lapply(Comp.E5, SM.Comp.org, cord[["E5"]])
Comp.E5.anom <- lapply(Comp.E5.anom, SM.Comp.org, cord[["E5"]])
Comp.CCI <- lapply(Comp.CCI, SM.Comp.org, cord[["CCI"]])
Comp.CCI.anom <- lapply(Comp.CCI.anom, SM.Comp.org, cord[["CCI"]])
Test.comp.E5 <- lapply(Test.comp.E5, SM.Comp.org, cord[["E5"]])
Test.comp.CCI <- lapply(Test.comp.CCI, SM.Comp.org, cord[["CCI"]])

save(Comp.E5, Comp.CCI , list=c("Comp.E5","Comp.CCI"), file="SatPer_Composites.RData")
save(Comp.E5.anom, Comp.CCI.anom,list=c("Comp.E5.anom","Comp.CCI.anom"), file="SatPer_anom_Composites.RData")
save(Test.comp.E5, Test.comp.CCI , list=c("Test.comp.E5","Test.comp.CCI"), file="SatPer_Composites_Ttest.RData")

#### shapefile composites statistical significance ----
library(reshape)
library(raster)

load("SatPer_Composites.RData")
load("SatPer_anom_Composites.RData")
load("SatPer_Composites_Ttest.RData")

# organizing the statistical significance in tables for creating the shapefiles
data.sig <- list()
data.sig[["ERA5L"]] <- melt(Test.comp.E5, id=c("x","y"))
data.sig[["CCI"]] <- melt(Test.comp.CCI, id=c("x","y"))

data.sig <- lapply(data.sig, function(x){ y <- x[!is.na(x$value),]; z <- y[y$value!=-999,]; return(z)})
data.sig <- lapply(data.sig,`colnames<-`,c("lon", "lat", "Season", "value", "Phase", "Mode" ))

shape.sig <- list()  # the nested "for"s are necesary because it's not possible to create the shapefiles all together at the same time
for(dataset in c("ERA5L", "CCI")){    shape.sig[[dataset]] <- list()
  
  for(phase in c("Pos","Neg")){   shape.sig[[dataset]][[phase]] <- list()
  
    for(mode in c("AMM","Atl3")){  shape.sig[[dataset]][[phase]][[mode]] <- list()
    print(paste(dataset,"-",phase,"-",mode))
    
      for(i in seasons){
        aux.t <- data.sig[[dataset]] 
        aux.t <- subset(aux.t,Phase==phase & Season == i & Mode== mode)  # aux.t <- aux.t[aux.t$Phase==Phase & aux.t$Season == i & aux.t$Mode== Mode,]
        aux.shp <- rasterFromXYZ(aux.t[,c("lon","lat","value")], digits=3) <=0.05  # creating the shapefile as a one line dataframe
        aux.shp <- sf::st_as_sf(aggregate(rasterToPolygons( aggregate(aux.shp, fact=5, fun="modal"), fun = function(x){x==1})))
        
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
save(shape.sig, "shape.sig", file="SatPer_Composites_Ttest_shape.RData")



#### seasonal plotting ----
load("SatPer_Composites.RData")
load("SatPer_anom_Composites.RData")
load("SatPer_Composites_Ttest.RData")
load("SatPer_Composites_Ttest_shape.RData")

library(ggplot2)
library(RColorBrewer)
library(ggpattern)
library(sf)

seasons <- c("DJF","MAM","JJA","SON")
Basins <- shapefile("../../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")
SA <- shapefile("../../../01_DataSets/South_America/South_America.shp")


# organizing the composites
data <- melt(Comp.E5, id=c("x","y")); data$dataset <- "ERA5L"; data <- data[!is.nan(data$value),]
aux <- melt(Comp.CCI, id=c("x","y")); aux$dataset <- "CCI"; aux <- aux[!is.nan(aux$value),]
data <- rbind(data,aux); colnames(data) <- c("lon", "lat", "Season", "value", "Phase", "Mode", "dataset" )
data$Phase <- factor(data$Phase, levels = c("Pos","Neg"))

data.anom <- melt(Comp.E5.anom, c("x","y")); data.anom$dataset <- "ERA5L"; data.anom <- data.anom[!is.nan(data.anom$value),]
aux <- melt(Comp.CCI.anom, id=c("x","y")); aux$dataset <- "CCI"; aux <- aux[!is.nan(aux$value),]
data.anom <- rbind(data.anom,aux); colnames(data.anom) <- c("lon", "lat", "Season", "Anom", "Phase", "Mode", "dataset" )
data.anom$Phase <- factor(data.anom$Phase, levels = c("Pos","Neg"))

# plotting
at.m <- seq(-10,10,length.out=11); at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]
col2alpha <- function(someColor, alpha=100){ newColor <- col2rgb(someColor); apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}
paleta <- brewer.pal(10,"RdYlBu") 
# paleta <- col2alpha(paleta, alpha = 0.4*255)


plot.sig <- function( Data.anom, Shape.spec, mode, df.rect, pat.col="black", lg.pos="vertical"){
  
  data.anom <- within(data.anom, data.class <- Anom >0)

  Phase.labs <- c("Positive", "Negative"); names(Phase.labs) <- c("Pos","Neg")
  p1 <- ggplot() + facet_wrap(. ~ Phase,  labeller= labeller(Phase = Phase.labs))+
    geom_raster(data= Data.anom, aes(lon, lat, fill=Anom),na.rm=T)+
    
    scale_fill_stepsn(colours=paleta, breaks=at.m,
                      values=scales::rescale(at.m.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                      guide=guide_colorsteps(even.steps = F, barwidth=unit(10,"cm")), name="SM\nSat %")+
    
    geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
    geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
    
    geom_rect(data=df.rect,inherit.aes = F,aes(xmin=xmi, xmax=xma, ymin=ymi, ymax=yma),color="black",linetype="dashed",linewidth=1, fill=NA, show.legend = FALSE)+
    
    geom_sf_pattern(data = Shape.spec, pattern= "circle", pattern_density=0.1, pattern_spacing= 0.02, fill="00", colour="0", pattern_colour=pat.col)+
    
    # geom_contour(data=subset(data.anom, data.class), aes(lon,lat, z=Anom, group= data.class), color= "#c51b7d",binwidth = 5, linewidth=1)+
    # geom_contour(data=subset(data.anom, !data.class), aes(lon,lat, z=Anom, group=data.class), color= "#b35806",binwidth = 5, linewidth=1)+

    coord_sf(xlim=c(-80,-35),ylim=c(-20,13))+
    # scale_y_continuous(position="right")+
    
    labs(x="Longitude [°]",y="Latitude [°]")+
    theme_bw() +
    theme(strip.background = element_blank(),strip.text.x = element_blank(),axis.text.x= element_blank(), axis.title.x = element_blank(),
          legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=15),
          legend.position = "bottom")
  
  # print(p1)
  return(p1)
}

rect.Atl3 <- data.frame(xmi=-65, xma=-57, ymi=0, yma=6,Season=as.factor("JJA"), Phase=as.factor("Pos"))

Atl3.E5 <- plot.sig(Data.anom= subset(data.anom, dataset == "ERA5L" & Mode =="Atl3" & Season == "JJA"),
                    Shape.spec= subset(shape.sig, dataset=="ERA5L" & Mode=="Atl3" & Season == "JJA"),
                    df.rect= rect.Atl3,
                    mode="Atl3")
Atl3.E5

rect.AMM <- data.frame(xmi=c(-70,-70,-73), 
                       xma=c(-58,-60,-67), 
                       ymi=c(3,-5,2), 
                       yma=c(10,3,7),
                       Season= factor(c("MAM","JJA","SON"), levels = c("MAM","JJA","SON")),
                       Phase=as.factor(rep("Pos",3)))
for ( i in c("MAM","JJA","SON")){
  p <- plot.sig(Data.anom= subset(data.anom, dataset == "ERA5L" & Mode =="AMM" & Season == i),
                Shape.spec= subset(shape.sig, dataset=="ERA5L" & Mode=="AMM" & Season == i),
                df.rect= subset(rect.AMM, Season == i),
                mode = "AMM")
  print(p)
}


AMM.cci <- plot.sig(data.spec= subset(data, dataset == "CCI" & Mode =="AMM" & Season != "DJF"),
                    data.anom = subset(data.anom, dataset == "CCI" & Mode =="AMM" & Season != "DJF"),
                    shape.spec= subset(shape.sig, dataset=="CCI" & Mode=="AMM" & Season != "DJF"))



Atl3.cci <- plot.sig(data.spec= subset(data, dataset == "CCI" & Mode =="Atl3" & Season == "JJA"),
                     data.anom = subset(data.anom, dataset == "CCI" & Mode =="Atl3" & Season == "JJA"),
                    shape.spec= subset(shape.sig, dataset=="CCI" & Mode=="Atl3" & Season == "JJA"))

library(gridExtra)
grid.arrange(AMM.E5, Atl3.E5, ncol=2, widths=c(2.05,1))

grid.arrange(AMM.cci, Atl3.cci, ncol=2, widths=c(2.05,1))


#### plotting positive less negative ----
load("SatPer_Composites.RData")

seasons <- c("DJF","MAM","JJA","SON")
data.Pos_Neg.E5 <- lapply(Comp.E5, function(x){ y <- x$Pos[,seasons] - x$Neg[,seasons]; y <- cbind(x$Pos[,c("x","y")], y); return(y)})
data.Pos_Neg.CCI <- lapply(Comp.CCI, function(x){ y <- x$Pos[,seasons] - x$Neg[,seasons]; y <- cbind(x$Pos[,c("x","y")], y); return(y)})

# organizing the composites
data.E5 <- melt(data.Pos_Neg.E5, id=c("x","y")); data.E5 <- data.E5[!is.nan(data.E5$value),]; colnames(data.E5) <- c("lon","lat","season", "value", "Mode")
data.CCI <- melt(data.Pos_Neg.CCI, id=c("x","y")); data.CCI <- data.CCI[!is.nan(data.CCI$value),]; colnames(data.CCI) <- c("lon","lat","season", "value", "Mode")

at.m <- seq(-10,10,length.out=11); at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]
col2alpha <- function(someColor, alpha=100){ newColor <- col2rgb(someColor); apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}
paleta <- brewer.pal(10,"RdYlBu") 

p1 <- ggplot(subset(data.E5, Mode=="AMM" & season!="DJF")) + facet_wrap(.~ season)+
  geom_raster( aes(lon, lat, fill=value),na.rm=T)+
  
  scale_fill_stepsn(colours=paleta, breaks=at.m,
                    values=scales::rescale(at.m.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")), name="SM Sat %")+
  
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  
  coord_sf(xlim=c(-82,-35),ylim=c(-20,15))+

  labs(x="Long",y="Lat",title=paste("SM Percentage of saturation Composite (Pos - Neg) AMM (1980-2020)"), caption = "ERA5-Land")+
  theme_bw()+ theme(legend.position = "none")

p2 <- ggplot(subset(data.E5, Mode=="Atl3" & season=="JJA")) + facet_wrap(.~ season)+
  geom_raster(aes(lon, lat, fill=value),na.rm=T)+
  
  scale_fill_stepsn(colours=paleta, breaks=at.m,
                    values=scales::rescale(at.m.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                    guide=guide_colorsteps(even.steps = F, barheight=unit(10,"cm")), name="SM Sat %")+
  
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  
  coord_sf(xlim=c(-82,-35),ylim=c(-20,15))+
  
  labs(x="Long",y="Lat",title=paste("Atl3"), caption = "ERA5-Land")+
  theme_bw()

grid.arrange(p1, p2, ncol=2, widths=c(2.05,1))
