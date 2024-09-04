rm(list=ls())
cat("\014")

library(raster)
source("../../../read_brickToPoints.R")
library(hydroTSM)
source("../../../time2season.R")
library(zoo)
library(ncdf4)

#### Dates ####
Dates <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"),by="month")
Dm.e5 <- seq(as.Date("1950-01-01"),as.Date("2020-12-31"),by="month")

years.m <- format(Dates, format="%Y"); Years <- unique(format(Dates, format="%Y")); Years <- Years[-1]
years.m.e5 <- format(Dm.e5, format="%Y"); Years.e5 <- unique(format(Dm.e5, format="%Y")); Years.e5 <- Years.e5[-1]

seasons <- c("DJF","MAM","JJA","SON")

Season.y <- paste0(years.m[-1], # eliminate the 1st value of year so as to the Dec of the 1st year is on the next seasonal year
                   "-",time2season(Dates)[-length(Dates)]); Season.y <- c(Season.y,Season.y[length(Season.y)])
Year.s <- unique(Season.y)
Season.y.e5 <- paste0(years.m.e5[-1], # eliminate the 1st value of year so as to the Dec of the 1st year is on the next seasonal year
                   "-",time2season(Dm.e5)[-length(Dm.e5)])
Year.s.e5 <- unique(Season.y.e5)

##### load data #####
cord <- list()

# cdo expr,NetRad=(str+ssr)/86400; -merge [ ERA5L_NetThermalRad_1950_2020.nc ERA5L_NetSolarRad_1950_2020.nc ] ERA5L_NetRad_1950-2020.nc
# cdo seasmean ERA5L_NetRad_1950-2020.nc ERA5L_NetRad_seasonal.nc
f2 <- nc_open("../../../01_DataSets/05_SRad/ERA5L_NetRad_seasonal.nc") # already transformed to W/m2

E5.s <- ncvar_get(f2, varid = "NetRad") # To convert to watts per square metre (W m-2), the accumulated values should be divided by the accumulation period expressed in seconds.
# Monthly means of mean daily fluxes, The accumulations in monthly means of daily means have been scaled to have units that include "per day"
# https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790

lon <- f2$dim[[which(names(f2$dim)=="longitude" | names(f2$dim)=="lon") ]]$vals; lat <- f2$dim[[which(names(f2$dim)=="latitude" | names(f2$dim)=="lat" )]]$vals
cord[["E5"]] <- expand.grid(lon=lon, lat=lat)


# SSTs
AMM <- read.zoo(read.csv("../../03_Composites/AMM_std_1980-2020.csv"),index.column = 1)
Atl3 <- read.zoo(read.csv("../../03_Composites/Atl3_std_1980-2020.csv"),index.column = 1)

#### crop 1980 - onwards & seasonal mean ####
E5.s <- E5.s[,,match(Year.s,Year.s.e5)]
Mean.E5 <- list()
E5.s.t <- list()

for( i in seasons){
  print(paste("seasonal mean - ", i))
  Mean.E5[[i]] <- apply( E5.s[,,substr(Year.s,6,8) == i], c(1,2), mean)
  Mean.E5[[i]] <- matrix(Mean.E5[[i]], nrow= prod( dim(Mean.E5[[i]])), ncol=1)
  E5.s.t[[i]] <- matrix( E5.s[,,substr(Year.s,6,8) == i], nrow=prod( dim(E5.s)[1:2]), ncol= length(Years))
}

#### composite of Rn  in the positive and negative phase of each phenomenon ----
Comp.E5 <- list(); Comp.E5[["AMM"]] <- list(); Comp.E5[["Atl3"]] <- list(); Comp.E5 <- lapply(Comp.E5, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )

Test.comp.E5 <- list(); Test.comp.E5[["AMM"]] <- list(); Test.comp.E5[["Atl3"]] <- list()
Test.comp.E5 <- lapply(Test.comp.E5, function(x){ x[["Pos"]]<- list(); x[["Neg"]]<- list(); return(x)} )

Test.comp.f <- function(Phase, neutral){

  Result <- rep(NA , nrow(Phase))
  # retriving cells without value and with constant series in which case it cannot be evaluated
  Idx <- which(apply(Phase,1,function(x) sum(!is.na(x)) > 1))
  a <- match(which(apply(Phase,1,function(x) sd(x)==0)), Idx); a <- a[!is.na(a)]; if(length(a!=0)>0) Idx <- Idx[- a]
  b <- match(which(apply(neutral,1,function(x) sd(x)==0)), Idx); b <- b[!is.na(b)]; if(length(b!=0)) Idx <- Idx[-b]
  # needed to use a loop because those are two datasets that cannot change with an apply (unless we built and array and that's more complicated)
  for(j in Idx) Result[j] <- t.test( Phase[j,], neutral[j, ])$p.value
  Result[a] <- -999; if(length(b!=0)) Result[b] <- -999
  return(Result)
}

for( i in seasons){
  print(paste(i, "- AMM - ERA5L"))
  Pos <-  E5.s.t [[i]][, AMM[,i]>= 1]; Neg <-  E5.s.t [[i]][, AMM[,i]<= -1]; Neutral <-  E5.s.t [[i]][, AMM[,i]<= 1 & AMM[,i]>= -1]
  Comp.E5[["AMM"]][["Pos"]][[i]] <- apply(Pos ,1, mean, na.rm=T) - Mean.E5[[i]]
  Comp.E5[["AMM"]][["Neg"]][[i]] <- apply(Neg ,1, mean, na.rm=T) - Mean.E5[[i]]
  Test.comp.E5[["AMM"]][["Pos"]][[i]] <- Test.comp.f(Pos, Neutral); Test.comp.E5[["AMM"]][["Neg"]][[i]] <- Test.comp.f(Neg, Neutral)


  print("- Atl3 - ERA5L")
  Pos <-  E5.s.t [[i]][, Atl3[,i]>= 1]; Neg <-  E5.s.t [[i]][, Atl3[,i]<= -1]; Neutral <-  E5.s.t [[i]][, Atl3[,i]<= 1 & Atl3[,i]>= -1]
  Comp.E5[["Atl3"]][["Pos"]][[i]] <- apply(Pos ,1, mean, na.rm=T) - Mean.E5[[i]]
  Comp.E5[["Atl3"]][["Neg"]][[i]] <- apply(Neg ,1, mean, na.rm=T) - Mean.E5[[i]]
  Test.comp.E5[["Atl3"]][["Pos"]][[i]] <- Test.comp.f(Pos, Neutral); Test.comp.E5[["Atl3"]][["Neg"]][[i]] <- Test.comp.f(Neg, Neutral)

}

SM.Comp.org <- function(x, cord){ x$Pos <- as.data.frame(x$Pos); x$Pos <- cbind(cord,x$Pos); x$Neg <- as.data.frame(x$Neg); x$Neg <- cbind(cord,x$Neg); return(x)}

Comp.E5 <- lapply(Comp.E5, SM.Comp.org, cord[["E5"]])

Test.comp.E5 <- lapply(Test.comp.E5, SM.Comp.org, cord[["E5"]])


save(Comp.E5, list=c("Comp.E5"), file="Radiation_Composites.RData")
save(Test.comp.E5, list=c("Test.comp.E5"), file="Radiation_Composites_Ttest.RData")

#### seasonal plotting ####
library(reshape)
library(ggplot2)
library(raster)
library(RColorBrewer)

load("Radiation_Composites.RData")
load("Radiation_Composites_Ttest.RData")
seasons <- c("DJF","MAM","JJA","SON")


# organizing the statistical significance in tables for creating the shapefiles

to_polygon <-function (data.sig){ # create polygons dataframes for each phase
  aux.sig <- list()
  for (i in seasons){
    aux.t <- data.sig[,c("lon","lat",i)]; aux.t <- subset(aux.t, !is.na(i))
    aux.shp <-  rasterFromXYZ(aux.t , digits=3) <=0.05  # creating the shapefile as a one line dataframe
    aux.shp <- sf::st_as_sf(aggregate(      rasterToPolygons( aggregate(aux.shp, fact=5, fun="modal"), fun = function(x){x==1})     ))
    aux.shp$Season <- i # attaching the fields to the shapefile

    aux.sig [[i]] <- aux.shp
  }
  return(do.call(rbind,aux.sig))
}

data.sig <- list()
n.modes <- c("AMM","Atl3")
dataset <- "ERA5L"
shape.sig <- lapply(Test.comp.E5,
                              function(Mode){
                                n.phase <- names(Mode); 
                                mode.r <- lapply(Mode, function(Phase){ phase.r <- to_polygon(Phase); # For both phases of the specific mode, create the polygons
                                  return(phase.r)}); 
                                for (phase in n.phase){ mode.r[[phase]]$Phase <- phase}; # attach the identifier of each phase
                                return(do.call(rbind,mode.r))
                                } # return and keep with the other mode
                              )
for (mode in n.modes){shape.sig[[mode]]$Mode <- mode}
shape.sig <- do.call(rbind, shape.sig)
shape.sig$Season <- factor(shape.sig$Season, levels = seasons); shape.sig$Phase <- factor(shape.sig$Phase, levels = c("Pos","Neg"))

library(ggpattern)
library(sf)
Basins <- shapefile("../../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")
SA <- shapefile("../../../01_DataSets/South_America/South_America.shp")

# organizing the composites
data <- melt(Comp.E5, id=c("lon","lat")); data$dataset <- "ERA5L"; data <- data[!is.nan(data$value),]
colnames(data) <- c("lon", "lat", "Season", "value", "Phase", "Mode", "dataset" )
data$Phase <- factor(data$Phase, levels = c("Pos","Neg")); data$Mode <- factor(data$Mode, levels = n.modes)


# plotting
range <- max(abs(min(data$value)),max(data$value ))
at.m <- round(seq(-range, range, length.out = 12),2); at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]

plot.sig <- function(data.spec, shape.spec, dataset, lg.pos = "right"){
  mode <- as.character(unique(data.spec$Mode))
  p1 <- ggplot() + facet_grid(Phase ~ Season, switch = "y")+
    geom_raster(data= data.spec, aes(lon, lat, fill=value),na.rm=T)+
    
    scale_fill_stepsn(colours=rev(brewer.pal(11,"RdBu")), breaks=at.m,
                      limits=c(min(at.m),max(at.m)), # scales::rescale(at.m.v,from=range(at.m)) not necesary because the breaks are evenly positioned
                      guide=guide_colorsteps( barheight=unit(10,"cm")), name="Rn anom. [W/m2]")+
    
    geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
    geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
    
    geom_sf_pattern(data = shape.spec, pattern= "circle", pattern_density=0.1, pattern_spacing=0.02, fill="00", colour="00", pattern_colour="purple")+
    
    coord_sf(xlim=c(-82,-35),ylim=c(-20,15))+
    # scale_y_continuous(position="right")+
    
    labs(x="Long",y="Lat",title=paste("Rn Composite", mode ,"(1980-2020)"), caption = dataset)+
    theme_bw()
  if(lg.pos!="right") p1 <- p1 + theme(legend.position = lg.pos)
  # print(p1)
  return(p1)
}

AMM.E5 <- plot.sig(data.spec= subset(data, dataset == "ERA5L" & Mode =="AMM" & Season !="DJF" ),
                   shape.spec= subset( shape.sig, Mode=="AMM" & Season !="DJF"),
                   dataset ="ERA5-Land", lg.pos="none")

Atl3.E5 <- plot.sig(data.spec= subset(data, dataset == "ERA5L" & Mode =="Atl3" & Season=="JJA"),
                    shape.spec= subset( shape.sig, Mode=="Atl3" & Season=="JJA"),
                    dataset ="ERA5-Land")


library(gridExtra)
grid.arrange(AMM.E5, Atl3.E5,ncol=2, widths=c(2.05,1))

