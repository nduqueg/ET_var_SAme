
rm(list = ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ncdf4)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
"%>%" = magrittr::`%>%`

# Dates ----
Dates.e5 <- seq(as.Date("1980-01-01"),as.Date("2020-12-30"), by="month")
Dates.e5r <- seq(as.Date("1950-01-01"),as.Date("2020-12-30"), by="month")
Dates.gl <- seq(as.Date("1980-01-01"),as.Date("2023-08-31"), by="month")
Dates.a <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month")

seasons <- c("DJF","MAM","JJA","SON")

Year_s <-function(Dates){
  years.m <- format(Dates, format="%Y"); Years <- format(Dates, format="%Y") %>% unique() %>% .[-1]
  Season.y <- paste0(years.m[-1],"-",time2season(Dates)[-length(Dates)]) %>% c(.,.[length(.)])
  Year.s <- unique(Season.y)
  return(list.r <- list(Year.s = Year.s, Season.y = Season.y, Years =Years))
}

Year.s <- Year_s(Dates.a)$Year.s; Years <- Year_s(Dates.a)$Years
Year.s.gl <- Year_s(Dates.gl)$Year.s
Year.s.e5 <- Year_s(Dates.e5)$Year.s
Year.s.e5r <- Year_s(Dates.e5r)$Year.s


## load data ----
# cdo seassum ERA5L_Ppt_1979_2020_cdo.nc ERA5L_Ppt_1979_2020_seasonal.nc
E5.s <- brick("03_detrended_ET_ERA5L_seasonal_1980_2020.nc")
GLEAM.s <- brick("03_detrended_ET_GLEAM_seasonal_1980_2020.nc"); GLEAM.s[[1]] <- GLEAM.s[[1]]*3/2 # add the Dec-1979 as the average of Jan & Feb-1980
Rn.s <- brick("../../../01_DataSets/05_SRad/ERA5L_NetRad_seasonal.nc")
SM.s <- brick("../../../01_DataSets/03_SM/01_ERA5L/ERA5L_SM_1L_seasonal.nc")

# sst standardize indices
AMM.seasons <- read.csv("../../03_Composites/AMM_std_1980-2020.csv"); colnames(AMM.seasons)[1] <- "Years"
Atl3.seasons <- read.csv("../../03_Composites/Atl3_std_1980-2020.csv"); colnames(Atl3.seasons)[1] <- "Years"
ELI.s <- read.csv("../../03_Composites/ELI_std_1980-2020_ERSST.csv"); colnames(ELI.s)[1] <- "Years"

## crop 1980 - onwards ----
E5.s <- E5.s[[match(Year.s, Year.s.e5)]]
GLEAM.s <- GLEAM.s[[ match(Year.s, Year.s.gl) %>% .[!is.na( .)] ]]
Year.s.gl <- Year.s.gl[match(Year.s, Year.s.gl)] %>% .[!is.na( .)]
Rn.s <- Rn.s[[match(Year.s, Year.s.e5r)]]
SM.s <- SM.s[[match(Year.s, Year.s.e5r)]]

#### identification of SM saturation ----
cord <- list()
f2 <- nc_open("../../../01_DataSets/03_SM/01_ERA5L/ERA5L_SM_1L_1950-2020.nc")
lon <- f2$dim[[which(names(f2$dim)=="longitude" | names(f2$dim)=="lon") ]]$vals; lat <- f2$dim[[which(names(f2$dim)=="latitude" | names(f2$dim)=="lat" )]]$vals
cord[["E5"]] <- expand.grid(lon=lon, lat=lat)

# Soils
fsoil <- nc_open("../../../01_DataSets/01_ERA5L_soil_type.nc")
E5.soil <- ncvar_get(fsoil, varid="slt")
lon <- fsoil$dim[[which(names(fsoil$dim)=="longitude" | names(fsoil$dim)=="lon") ]]$vals; lon[lon>=180] <- lon[lon>=180] - 360
lat <- fsoil$dim[[which(names(fsoil$dim)=="latitude" | names(fsoil$dim)=="lat" )]]$vals
cord[["E5_soil"]] <- expand.grid(lon=lon, lat=lat)

# trim soils to tropical south america
E5.soil <- E5.soil[lon >= min(cord[["E5"]]$lon) & lon <= max(cord[["E5"]]$lon),
                   lat >= min(cord[["E5"]]$lat) & lat <= max(cord[["E5"]]$lat)]; E5.soil <- round(E5.soil,2)
cord[["E5_soil"]] <- subset(cord[["E5_soil"]], lon >= min(cord[["E5"]]$lon) & lon <= max(cord[["E5"]]$lon) & 
                              lat >= min(cord[["E5"]]$lat) & lat <= max(cord[["E5"]]$lat)); cord[["E5_soil"]] <- round(cord[["E5_soil"]],2)

# saturation values from ECMWF documentation (ECMWF, 2023)
soil.features <- data.frame(type=seq(1,7), sat= c(0.403, 0.439, 0.43, 0.52, 0.614, 0.766, 0.472))
E5.sat <- E5.soil
for(i in 1:7){
  E5.sat[E5.soil== i] <- soil.features$sat[i]
}
E5.sat <- cbind(cord[["E5_soil"]], as.numeric(E5.sat)); E5.sat[E5.sat[,3] == 0, 3] <- NA
E5.sat <- rasterFromXYZ(E5.sat)

SM.s.sat <- SM.s/E5.sat*100

### crop to boxes ----
GLEAM.box <- E5.box <- Rn.box <- SM.box <- list()
rect <- data.frame(xmi=c(-70,-70,-73,-65),
                   xma=c(-58,-60,-67,-57),
                   ymi=c(3,-5,2,0),
                   yma=c(10,3,7,6),
                   Season= factor(c("MAM","JJA","SON","JJA"), levels = c("MAM","JJA","SON")),
                   mode=c("AMM","AMM","AMM","Atl3"),
                   Phase = "Pos")

for(i in 1:nrow(rect)){
  E5.box[[i]] <- crop(E5.s, extent( rect[i, 1:4] %>% as.numeric())) %>% 
    cellStats(., stat="mean")
  GLEAM.box[[i]] <- crop(GLEAM.s, extent( rect[i, 1:4] %>% as.numeric()) ) %>% 
    cellStats(., stat="mean")
  Rn.box[[i]] <- crop(Rn.s, extent( rect[i, 1:4] %>% as.numeric())) %>% 
    cellStats(., stat="mean")
  SM.box[[i]] <- crop(SM.s.sat, extent( rect[i, 1:4] %>% as.numeric())) %>% 
    cellStats(., stat="mean")
}

### seasonal separation and anomalies ----
E5.box.s <- GLEAM.box.s <- Rn.box.s <- SM.box.s <-list()
for( i in 1:nrow(rect)){
  E5.box.s[[i]] <- matrix(NA,nrow=length(Years), ncol=4) %>% as.data.frame(); colnames(E5.box.s[[i]]) <- seasons
  GLEAM.box.s[[i]] <- matrix(NA,nrow=length(Years), ncol=4) %>% as.data.frame(); colnames(GLEAM.box.s[[i]]) <- seasons
  Rn.box.s[[i]] <- matrix(NA,nrow=length(Years), ncol=4) %>% as.data.frame(); colnames(Rn.box.s[[i]]) <- seasons
  SM.box.s[[i]] <- matrix(NA,nrow=length(Years), ncol=4) %>% as.data.frame(); colnames(SM.box.s[[i]]) <- seasons
  
  for(j in seasons){
    E5.box.s[[i]][1:length(Years) ,j] <- E5.box[[i]][substr(Year.s,6,8) == j]
    GLEAM.box.s[[i]][1:length(Years) ,j] <- GLEAM.box[[i]][substr(Year.s,6,8) == j]
    Rn.box.s[[i]][1:length(Years) ,j] <- Rn.box[[i]][substr(Year.s,6,8) == j]
    SM.box.s[[i]][1:length(Years) ,j] <- SM.box[[i]][substr(Year.s,6,8) == j]
  }
}
# could have found the anomalies with "scale" function
E5.box.mean <- sapply(E5.box.s, function(x) colMeans(x, na.rm = T)) %>% t()
GLEAM.box.mean <- sapply(GLEAM.box.s, function(x) colMeans(x, na.rm = T)) %>% t()
Rn.box.mean <- sapply(Rn.box.s, function(x) colMeans(x, na.rm = T)) %>% t()
SM.box.mean <- sapply(SM.box.s, function(x) colMeans(x, na.rm = T)) %>% t()
Box.anom <- list(); Box.anom[["ERA5L"]] <- Box.anom[["GLEAM"]] <- Box.anom[["Rn"]] <- list()

for( i in 1:nrow(rect)){
  Box.anom[["ERA5L"]][[i]] <- sweep(E5.box.s[[i]], 2, E5.box.mean[i, ], FUN="-") # substrictaion column-wise
  Box.anom[["GLEAM"]][[i]] <- sweep(GLEAM.box.s[[i]], 2, GLEAM.box.mean[i, ], FUN="-")
  Box.anom[["Rn"]][[i]] <- sweep(Rn.box.s[[i]], 2, Rn.box.mean[i, ], FUN="-")
}
Box.anom[["Rn"]] <- lapply(Box.anom[["Rn"]],function(x) apply(x,2,scale))

### correlations by season ----
Corr <- data.frame(ERA5L=numeric(), GLEAM=numeric(), Rn=numeric())
Sig <- data.frame(ERA5L=numeric(), GLEAM=numeric(), Rn=numeric())

for(j in c("ERA5L","GLEAM","Rn")){
  for(i in 1:nrow(rect)){
    a <- Box.anom[[ j]][[ i]][, as.character(rect$Season[i]) ]
    b <- get(paste0(rect$mode[i],".seasons"))[, as.character(rect$Season[i])]
    aux <- cor.test( a, b)
    Corr[i, j] <- round(aux$estimate,2)
    Sig[i, j] <- aux$p.value
  }
}

# text for the graph
mess <- apply(Corr, 2, function(x)paste("Corr =",x,"**"))
mess[which(Sig >=0.05)] <- substr(mess[which(Sig >=0.05)], 1, nchar(mess[which(Sig >=0.05)]) -2 )
mess <- cbind(mess, subset(rect,select=c("Season","mode"))) %>%  melt(., id=c("Season","mode")) %>% 
  within(.,{y.m=rep(c(38,-20,-20,20),3)})

# organize data for plotting ----
Box.anom <- lapply(Box.anom, function(x, Years) lapply(x, cbind.data.frame, Years), as.numeric(Years)) %>% melt(., id="Years")
colnames(Box.anom) <- c("Years", "Season","value","Box","Dataset")

SM.box.s2 <- lapply(SM.box.s, function(x, Years) cbind.data.frame(x, Years), as.numeric(Years)) %>% melt(., id="Years")
colnames(SM.box.s2) <- c("Years", "Season","value","Box")

# filter the time series for the specific panels to be plotted
ET.series <- SM.series  <- data.frame(Years=numeric(), Season=character(),value=numeric(),Box=numeric())
ET.series$Dataset <- character()

for(i in 1:nrow(rect)){
  ET.series <- rbind(ET.series, subset(Box.anom, Season == as.character(rect[i, "Season"]) & Box == i))
  SM.series <- rbind(SM.series, subset(SM.box.s2, Season == as.character(rect[i, "Season"]) & Box == i))
}

Rn.series <- subset(ET.series, Dataset=="Rn")
ET.series <- subset(ET.series, Dataset!="Rn")
SM.series.g <- SM.series %>%
  dplyr::mutate(Years.end = Years + 0.5, Years.ini = Years - 0.5,
                y.b = dplyr::case_when(Box==1 ~ -45, Box==2 ~ -30, Box==3 ~ -30, Box==4 ~ -20)) %>% 
  dplyr::mutate(y.l=y.b + dplyr::case_when(Box==1 ~ -1.7, Box!=1 ~ -1),
                y.u=y.b + dplyr::case_when(Box==1 ~ 1.7, Box!=1 ~ 1))

# SSTs
AMM.series <- AMM.seasons %>% subset(., select=-DJF) %>% melt(., id=c("Years")) ; colnames(AMM.series)[2] <- "Season"
Atl3.series <- Atl3.seasons %>% subset(., select=c("Years","JJA")) %>% melt(., id=c("Years")); colnames(Atl3.series)[2] <- "Season"

ELI.series <- ELI.s %>% dplyr::mutate(Years.end = Years + 0.5, Years = Years - 0.5) %>% 
  melt(.,id=c("Years","Years.end")) %>% 
  dplyr::mutate(.,value = dplyr::case_when(value >=0.75 ~ "Pos", value <=-0.74 ~ "Neg")) %>%
  dplyr::mutate(.,y.graph = dplyr::case_when(variable =="MAM" ~ 45, variable =="JJA" ~ 30, variable =="SON" ~ 30)) %>%#, .default = "Neutral"
  subset(., variable != "DJF" & !is.na(value))
colnames(ELI.series)[3] <- "Season"


## plotting with ENSO on top and SM at bottom with saturation as color----
library(ggrepel)
ELI.series$Years <- ELI.series$Years + 0.5

Season.labs <- c("MAM (Austral Autumn)", "JJA (Austral Winter)", "SON (Austral Spring)"); names(Season.labs) <- c("MAM","JJA","SON")
col2alpha <- function(someColor, alpha=1){ newColor <- col2rgb(someColor); if(alpha <=1) alpha <- alpha * 255 else warning("Alpha should be between 0 and 1"); apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}
at.m <- seq(0,100,length.out=11); at.m.v <- seq(5,95,length.out=10)
paleta.m <- brewer.pal(10,"RdYlBu");  paleta.m <- col2alpha(paleta.m, alpha = 0.8)

f.mult <-  sd(ET.series$value,na.rm=T)/ sd(AMM.series$value)
f.mult.rn  <- sd(ET.series$value,na.rm=T)/ sd(Rn.series$value)


plot.ts <- function(ET.Series, SM.Series.g, Rn.Series, ELI.Series, Mode.Series, Mode="AMM", Mess){
  p <- ggplot()+
    geom_point(data=subset(ELI.Series, value=="Pos"), inherit.aes = F, aes(x=Years, y=y.graph), col="red")+
    geom_point(data=subset(ELI.Series, value=="Neg"), inherit.aes = F, aes(x=Years, y=y.graph), col="Blue")+
    
    geom_rect(data=SM.Series.g, inherit.aes = F, aes(xmin=Years.ini, xmax=Years.end, ymin=y.l, ymax=y.u, fill=value), col="00")+
    scale_fill_stepsn(colours=paleta.m, breaks=at.m,
                      values=scales::rescale(at.m.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                      guide=guide_colorsteps(even.steps = F, barwidth=unit(10,"cm")),
                      name="SM Sat [%]")+
    
    geom_line(data= ET.Series, aes( x=Years,y=value,col= Dataset))+
    scale_x_continuous(breaks = seq(1980,2020,by=4),minor_breaks = seq(1980,2020,by=1))+
    scale_color_manual(values=c("#006837","#78c679","#ff7f00"))+
    
    geom_line(data=Mode.Series,aes(Years, value*f.mult),col="black",linetype="dashed")+
    geom_line(data=Rn.Series, aes(Years, value*f.mult.rn,col= Dataset))+
    scale_y_continuous(name="ET anomaly [mm]",sec.axis = sec_axis(trans=~./f.mult,name=paste(Mode, "& Rn [SD]"), breaks = seq(-4,3,by=1)))+
    
    geom_text_repel(data = Mess,inherit.aes = F, aes(x=2020, y=y.m, label=value, col=variable),segment.color = 'transparent', size=5)+
    
    labs(xlab="Year")+
    theme_bw()+ theme(legend.position = "bottom",#legend.background = element_rect(linewidth = 0.5, linetype="solid", colour ="black"),
                      legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=15),
                      legend.title.position="top",
                      axis.ticks.length = unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1, size=16),axis.title.x = element_blank(),
                      axis.text.y = element_text(margin=margin(0,7,5,0,"pt"), size=20),axis.text.y.left = element_text(color = "#006837"),
                      axis.title.y = element_text(size=20),axis.title.y.left = element_text(color = "#006837",face="bold"),
                      panel.grid = element_line(linetype="dashed",color="lightgrey")
                      )
  return(p)
}

rect <- within(rect, Season <- as.character(Season))
for ( i in 1:4){
  plot.ts(ET.Series= subset(ET.series, Box==i),
          SM.Series.g= subset(SM.series.g, Box==i),
          Rn.Series= subset(Rn.series, Box==i ), 
          ELI.Series= subset(ELI.series, Season== rect$Season[i]),
          Mode.Series= subset(get(paste0(rect$mode[i],".series")), Season ==rect$Season[i]),
          Mode=rect$mode[i], 
          Mess= subset(mess, mode== rect$mode[i] & Season ==rect$Season[i])) %>% print()
}

