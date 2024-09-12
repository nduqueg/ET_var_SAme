
rm(list = ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
"%>%" = magrittr::`%>%`

# Dates ----
Dates.e5 <- seq(as.Date("1950-01-01"),as.Date("2020-12-30"), by="month")
Dates.ms <- seq(as.Date("1979-02-01"),as.Date("2020-11-30"), by="month")
Dates.a <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month")

seasons <- c("DJF","MAM","JJA","SON")

years.m <- format(Dates.a, format="%Y"); Years <- format(Dates.a, format="%Y") %>% unique() %>% .[-1] %>% as.numeric()
Season.y <- paste0(years.m[-1],"-",time2season(Dates.a)[-length(Dates.a)]) %>% c(.,.[length(.)])
Year.s <- unique(Season.y)

years.m.ms <- format(Dates.ms, format="%Y"); Years.ms <- format(Dates.ms, format="%Y") %>% unique() %>% .[-1]
Season.y.ms <- paste0(years.m.ms[-1],"-",time2season(Dates.ms)[-length(Dates.ms)]) %>% c(. , .[length( .)])
Year.s.ms <- unique(Season.y.ms)

years.m.e5 <- format(Dates.e5, format="%Y"); Years.e5 <- format(Dates.e5, format="%Y") %>% unique() %>% .[-1]
Season.y.e5 <- paste0(years.m.e5[-1],"-",time2season(Dates.e5)[-length(Dates.e5)]) %>% c(. , .[length( .)])
Year.s.e5 <- unique(Season.y.e5)

## load data ----
# cdo seassum ERA5L_Ppt_1979_2020_cdo.nc ERA5L_Ppt_1979_2020_seasonal.nc
E5.s <- brick("./01_data/03_SM/ERA5L_SM_1L_seasonal.nc")
ESA.s <- brick("./01_data/03_SM/ESACCI_SM_seasonal_80-DJF_20-SON.nc")
P.s <- brick("./01_data/02_Ppt/MSWEP_Ppt_1979_2020_seasonal.nc")

# sst standardize indices
AMM.seasons <- read.csv("./01_data/01_SSTindices/AMM_std_1980-2020.csv"); colnames(AMM.seasons)[1] <- "Years"
Atl3.seasons <- read.csv("./01_data/01_SSTindices/Atl3_std_1980-2020.csv"); colnames(Atl3.seasons)[1] <- "Years"
ELI.s <- read.csv("./01_data/01_SSTindices/ELI_std_1980-2020_ERSST.csv"); colnames(ELI.s)[1] <- "Years"

## crop 1980 - onwards ----
P.s <- P.s[[match(Year.s, Year.s.ms)]]
E5.s <- E5.s[[match(Year.s, Year.s.e5)]]
# ESA CCI SM is already in the time period specified

### crop to boxes ----
P.box <- E5.box <- ESA.box <- list()
rect <- data.frame(xmin=c(-70,-70,-73,-65), 
                   xmax=c(-58,-60,-67,-57), 
                   ymin=c(3,-5,2,0), 
                   ymax=c(10,3,7,6),
                   Season= factor(c("MAM","JJA","SON","JJA"), levels = c("MAM","JJA","SON")),
                   mode=c("AMM","AMM","AMM","Atl3"), 
                   codes=c("NEOri","C.Ama","SW.Ori","Guaianas"),stringsAsFactors = F)

mess.box <- apply(rect[,1:4],1,function(x){paste0("[",x[1]*(-1),"°W,",x[2]*(-1),"°W]x[",x[3],"°N,",x[4],"°N]")}) %>% 
  cbind( subset(rect, select=c(Season,mode,codes)), .) %>% 
  cbind(., data.frame(x=1984,y=-270))

  
for(i in 1:nrow(rect)){
  P.box[[i]] <- crop(P.s, extent( rect[i, 1:4] %>% as.numeric())) %>% 
    cellStats(., stat="mean")
  E5.box[[i]] <- crop(E5.s, extent( rect[i, 1:4] %>% as.numeric())) %>%
    cellStats(., stat="mean")
  ESA.box[[i]] <- crop(ESA.s, extent( rect[i, 1:4] %>% as.numeric())) %>%
    cellStats(., stat="mean")
}


### seasonal separation and anomalies ----
P.box.s <- E5.box.s <- ESA.box.s <- list()
for( i in 1:nrow(rect)){
  E5.box.s[[i]] <- matrix(NA,nrow=length(Years), ncol=4); colnames(E5.box.s[[i]]) <- seasons
  ESA.box.s[[i]] <- matrix(NA,nrow=length(Years), ncol=4); colnames(ESA.box.s[[i]]) <- seasons
  P.box.s[[i]] <- matrix(NA,nrow=length(Years), ncol=4); colnames(P.box.s[[i]]) <- seasons
  
  for(j in seasons){
    E5.box.s[[i]][1:length(Years) ,j] <- E5.box[[i]][substr(Year.s,6,8) == j]
    ESA.box.s[[i]][1:length(Years) ,j] <- ESA.box[[i]][substr(Year.s,6,8) == j]
    P.box.s[[i]][1:length(Years) ,j] <- P.box[[i]][substr(Year.s,6,8) == j]
  }
}

P.anom <- lapply(P.box.s, function(x) apply(x,2,scale, scale=F)) %>% lapply(., cbind.data.frame, Years)
E5.anom <- lapply(E5.box.s, function(x) apply(x,2,scale,scale=T)) %>% lapply(., cbind.data.frame, Years)
ESA.anom <- lapply(ESA.box.s, function(x) apply(x,2,scale,scale=T)) %>% lapply(., cbind.data.frame, Years)

### correlations by season ----
Corr <- list(); for(i in c("MSWEP","E5","ESA")) Corr[[i]] <- matrix(NA,nrow=nrow(rect), ncol=4) %>% as.data.frame()
Sig <-  list(); for(i in c("MSWEP","E5","ESA")) Sig[[i]] <- matrix(NA,nrow=nrow(rect), ncol=4) %>% as.data.frame()
Corr <- lapply(Corr,`colnames<-`,seasons)
Sig <- lapply(Sig,`colnames<-`,seasons)

for(j in c("MSWEP","E5","ESA")){
  
  for(i in 1:nrow(rect)){
    if(j=="MSWEP"){ supply <- P.box.s[[i]] 
    }else if(j=="E5"){ supply <- E5.box.s[[i]]
    }else if(j=="ESA"){ supply <- ESA.box.s[[i]]
    }
    
      a <- supply[, as.character(rect$Season[i]) ]
      b <- paste0(rect$mode[i],".seasons") %>% get() %>% dplyr::pull( rect$Season[i])
      aux <- cor.test( a, b)
      Corr[[j]][i, as.character(rect$Season[i])] <- round(aux$estimate,2)
      Sig[[j]][i, as.character(rect$Season[i])] <- aux$p.value
      
      Corr[[j]]$Code <- rect$codes
      Sig[[j]]$Code <- rect$codes
      
      Corr[[j]]$mode <- rect$mode
      Sig[[j]]$mode <- rect$mode
  }
}

# organize data for plotting ----

Corr <- Corr %>% melt(id=c("Code","mode")) %>% .[!is.na(.$value),]
Sig <- Sig %>% melt(id=c("Code","mode"))%>% .[!is.na(.$value),]
data <- merge.data.frame(Corr,Sig, by=c("variable","L1","Code","mode")) %>% 
  magrittr::set_colnames(., c("Season","Input","Code","mode","Corr","Sig")) %>% 
  within(., {
    Input[Input=="E5"] <- "SM-ERA5-Land"
    Input[Input=="ESA"] <- "ESA-CCI-SM"
  })

mess <- data %>% dplyr::mutate(., p.value = dplyr::case_when(Sig <=0.05 ~ "**", Sig <=0.1 ~ "*",.default = " "), .keep="unused") %>%
  dplyr::mutate(., y.m=dplyr::case_when(Season =="MAM" ~ -300, Season =="JJA" ~ -300, Season =="SON" ~ -250), .keep="all") %>% 
  dplyr::mutate(., Corr= paste0("Corr.",Input,"= ",Corr,p.value))

P.anom <- melt(P.anom, id="Years"); colnames(P.anom) <- c("Years","Season","value","Box")
E5.anom <- melt(E5.anom, id="Years"); colnames(E5.anom) <- c("Years","Season","value","Box")
ESA.anom <- melt(ESA.anom, id="Years"); colnames(ESA.anom) <- c("Years","Season","value","Box")

# filter the time series for the specific panels to be plotted
P.series <- E5.series <- ESA.series <- data.frame(Years=numeric(), Season=character(),value=numeric(),Box=numeric())
for(i in 1:nrow(rect)){
  
  P.series <- rbind(P.series, subset(P.anom, Season == as.character(rect[i, "Season"]) & Box == i))
  E5.series <- rbind(E5.series, subset(E5.anom, Season == as.character(rect[i, "Season"]) & Box == i))
  ESA.series <- rbind(ESA.series, subset(ESA.anom, Season == as.character(rect[i, "Season"]) & Box == i))
  
}

# sst
AMM.series <- AMM.seasons %>% subset(., select= -DJF) %>% melt(., id=c("Years")) ; colnames(AMM.series)[2] <- "Season"
Atl3.series <- Atl3.seasons %>% subset(., select=c("Years","JJA")) %>% melt(., id=c("Years")); colnames(Atl3.series)[2] <- "Season"

ELI.series <- ELI.s  %>% 
  melt(.,id=c("Years")) %>% 
  dplyr::mutate(.,value = dplyr::case_when(value >=0.75 ~ "Pos", value <=-0.74 ~ "Neg")) %>%
  dplyr::mutate(.,y.graph = dplyr::case_when(variable =="MAM" ~ 250, variable =="JJA" ~ 200, variable =="SON" ~ 175)) %>%#, .default = "Neutral"
  subset(., variable != "DJF" & !is.na(value))
colnames(ELI.series)[2] <- "Season"

## plotting ----
library(ggrepel)
Season.labs <- c("MAM (Austral Autumn)", "JJA (Austral Winter)", "SON (Austral Spring)"); names(Season.labs) <- c("MAM","JJA","SON")

size.t <- 5

f.mult <-  sd(P.series$value,na.rm=T)/ sd(AMM.series$value)
# f.mult.q  <- sd(P.series$value,na.rm=T)/ sd(E5.series$value,na.rm=T)

plot.ts <- function(P.Series, E5.Series, ESA.Series, ELI.Series, Mode.Series, Mode="AMM", Mess, Mess.box){
  p <- ggplot()+
    
    geom_bar(data=P.Series, aes(x=Years,y=value),fill="#a6cee3",col="blue",stat="identity", linewidth=1)+
    
    geom_point(data=subset(ELI.Series, value=="Pos"), inherit.aes = F, aes(x=Years, y=y.graph), col="red")+
    geom_point(data=subset(ELI.Series, value=="Neg"), inherit.aes = F, aes(x=Years, y=y.graph), col="blue")+
    
    scale_x_continuous(breaks = seq(1980,2020,by=4),minor_breaks = seq(1980,2020,by=1))+
    
    geom_line(data=ESA.Series, aes(Years, value*f.mult), col= "#de77ae",linewidth=1)+
    geom_line(data=Mode.Series, aes(Years, value*f.mult),col="black",linetype="dashed",linewidth=1)+
    geom_line(data=E5.Series, aes(Years, value*f.mult), col= "#542788",linewidth=1)+
    
    scale_y_continuous(name="Ppt anomaly [mm]",sec.axis = sec_axis(trans=~./f.mult,name=paste(Mode, "& SM [SD]"), breaks = seq(-4,3,by=1)))+
    
    geom_text_repel(data = Mess, inherit.aes = F, aes(x=2006, y=y.m, label=Corr, col=Input),segment.color = 'transparent', size=size.t)+
    geom_text(data= Mess.box, aes(x,y,label=.), size=size.t)+
    
    scale_color_manual(values = c("#de77ae","blue","#542788"), name="Dataset")+
    
    theme_bw()+ theme(legend.position="bottom", axis.title.x = element_blank(), axis.text.x=element_text(size = 16),
                      axis.ticks.length=unit(-4, "pt"), axis.text.y = element_text(margin=margin(0,7,5,0,"pt"), size=16),axis.text.y.left = element_text(color = "blue",angle=90, hjust = 0.5),
                      axis.title.y = element_text(size=20), axis.title.y.left = element_text(color = "blue",face="bold"),
                      axis.title.y.right = element_text(color = "#542788",face="bold"),
                      panel.grid = element_line(linetype="dashed",color="lightgrey"),
                      legend.text = element_text(size = 16), legend.title = element_text(size=20), legend.background = element_rect(linewidth = 0.5, linetype="solid", colour ="black")
                      )
  
  return(p)
}
rect <- within(rect, Season <- as.character(Season))

for (i in 1:4){
  plot.ts(P.Series= subset(P.series, Box==i ), 
          E5.Series= subset(E5.series, Box==i ), 
          ESA.Series= subset(ESA.series, Box==i ), 
          ELI.Series= subset(ELI.series, Season== rect$Season[i]), 
          Mode.Series= subset(get(paste0(rect$mode[i],".series")), Season== rect$Season[i]), 
          Mode=rect$mode[i], 
          Mess= subset(mess,mode==rect$mode[i] & Season== rect$Season[i]), 
          Mess.box= subset(mess.box,mode==rect$mode[i] & Season== rect$Season[i])) %>% print()
}

  
