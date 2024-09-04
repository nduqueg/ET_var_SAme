rm(list=ls())
cat("\014")

library(reshape)
library(raster)
library(ncdf4)
library(RColorBrewer)
library(hydroTSM)
library(ggplot2)
library(metR)
"%>%" = magrittr::`%>%`

# DAtes ----
dates.m <- seq(as.Date("1980-01-01"),as.Date("2020-12-31"),by="month")
month.Multi <- format(dates.m,format="%m")
years <- factor(format(dates.m,format="%Y")); Years <- unique(years)
Season.y <- paste0(years[-1],"-",time2season(dates.m[-length(dates.m)]))
Year.s <- unique(Season.y)

seasons <- c("DJF","MAM","JJA","SON")

D.e5 <-seq(as.Date("1950-01-01"),as.Date("2020-12-31"),by="month")

Dm.sst <- seq(as.Date("1885-01-01"),as.Date("2020-12-30"),by="month")

## loading ----
basins <- shapefile("../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")

VIMF.anom <- list(); MDiv.anom <- list(); Ppt.anom <- list()
load("./01_VIMF/01_AnomComp_Atl3_VIMF.RData"); VIMF.anom[["Atl3"]] <- data.anom; rm(data.anom)
load("./02_VIMFDiv/02_Atl3_divVIM.RData"); MDiv.anom[["Atl3"]] <- data.anom2; rm(data.anom2)
load("./03_Ppt/03_Atl3_Comp_Ppt.RData"); Ppt.anom[["Atl3"]] <- data.anom2; rm(data.anom2)

load("./01_VIMF/01_AnomComp_AMM_VIMF.RData"); VIMF.anom[["AMM"]] <- data.anom; rm(data.anom)
load("./02_VIMFDiv/02_AMM_divVIM.RData"); MDiv.anom[["AMM"]] <- data.anom2; rm(data.anom2)
load("./03_Ppt/03_AMM_Comp_Ppt.RData"); Ppt.anom[["AMM"]] <- data.anom2; rm(data.anom2)

Test.comp <- list()
VIMF.mag.test <- function(Mt.E.test, Mt.N.test, lev.sig = 0.05){
  Mt.Mag.test <- list()
  Mt.E.test <- lapply(Mt.E.test, round, 3); Mt.N.test <- lapply(Mt.N.test, round, 3)
  Mt.Mag.test[["Pos"]] <- Mt.Mag.test[["Neg"]] <- matrix(NA, nrow=nrow( Mt.E.test[["Pos"]]), ncol= 4)
  Mt.Mag.test <- lapply(Mt.Mag.test, `colnames<-`, seasons)
  for (i in seasons){
    Mt.Mag.test[["Pos"]][,i] <- Mt.E.test[["Pos"]][,i] <lev.sig | Mt.N.test[["Pos"]][,i] <lev.sig
    Mt.Mag.test[["Neg"]][,i] <- Mt.E.test[["Neg"]][,i] <lev.sig | Mt.N.test[["Neg"]][,i] <lev.sig
  }
  Mt.Mag.test <- lapply( Mt.Mag.test, cbind, Mt.E.test$Pos[,c("lon","lat")])
  Mt.Mag.test <- lapply( Mt.Mag.test, melt, id=c("lon","lat"))
  Mt.Mag.test <- lapply( Mt.Mag.test,`colnames<-`, c("lon", "lat","Season","p.value"))
  return( Mt.Mag.test)
}
load("./01_VIMF/02_VIMF_AMM_Comp_Ttest_uv.RData")
Test.comp[["AMM"]] <- VIMF.mag.test(Mt.E.test, Mt.N.test, lev.sig = 0.1)

cord <- Mt.E.test$Pos[,c("lon","lat")]
load("./01_VIMF/02_VIMF_Atl3_Comp_Ttest_uv.RData")
Test.comp[["Atl3"]] <- VIMF.mag.test(Mt.E.test, Mt.N.test, lev.sig = 0.1)

### selection VIMF for ease visualization----

sele <- list()
# sele[["lon"]] <- seq(min(data$lon),max(data$lon),by=0.25)
# sele[["lat"]] <- seq(min(data$lat),max(data$lat),by=0.25)
sele[["lon"]] <- seq(-83.25,-29.25,by=0.75)
sele[["lat"]] <- seq(-21.25,16.25,by=0.75)

sele.fun <- function(VIMF.anom, sele){
  VIMF.anom[["Neg"]] <- VIMF.anom[["Neg"]][ which(VIMF.anom[["Neg"]]$lon %in% sele$lon) ,]; VIMF.anom[["Neg"]] <- VIMF.anom[["Neg"]][ which(VIMF.anom[["Neg"]]$lat %in% sele$lat) ,]
  VIMF.anom[["Pos"]] <- VIMF.anom[["Pos"]][ which(VIMF.anom[["Pos"]]$lon %in% sele$lon) ,]; VIMF.anom[["Pos"]] <- VIMF.anom[["Pos"]][ which(VIMF.anom[["Pos"]]$lat %in% sele$lat) ,]
  return(VIMF.anom)
}

VIMF.anom <- lapply(VIMF.anom, sele.fun, sele)
Test.comp <- lapply(Test.comp, sele.fun, sele)

Ppt.anom <- lapply(Ppt.anom, function(x) {x$Dir <- factor(x$Dir, levels = c("Pos","Neg")); return(x)})
MDiv.anom <- lapply(MDiv.anom, function(x) {x$Dir <- factor(x$Dir, levels = c("Pos","Neg")); return(x)})
#### plotting -----
plotting.anom <- function (VIMF.anom, MDiv.anom, Ppt.anom, df.rect, Spe.season, title.p = NA, xlim=c(-82,-40), ylim=c(-20,13), l.pos="bottom", y.axis=TRUE, x.axis=TRUE){
  # data.anom must be the dataFrame, not the whole list
  at.ppt <- c(-250,-150,-100,-75,-50,-25,25,50,75,100,150,250)
  at.ppt.v <- (at.ppt[-length(at.ppt)] - at.ppt[-1])/2 + at.ppt[-1]
  col2alpha <- function(someColor, alpha=1){ newColor <- col2rgb(someColor); if(alpha <=1) alpha <- alpha * 255 else warning("Alpha should be between 0 and 1"); apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}
  paleta <- brewer.pal(11,"BrBG"); paleta <- col2alpha(paleta, alpha = 0.5)
  
  VIMF.anom.p <- VIMF.anom[sapply(Spe.season,function(s,d) which(d == s),VIMF.anom$Season),]
  MDiv.anom.p <- MDiv.anom[sapply(Spe.season,function(s,d) which(d == s),MDiv.anom$Season),]
  
  MDiv.anom.p$data.class <- MDiv.anom.p$Anomaly>0
  
  Ppt.anom.p <- Ppt.anom[sapply(Spe.season,function(s,d) which(d == s),Ppt.anom$Season),]
  
  Phase.labs <- c("Positive", "Negative"); names(Phase.labs) <- c("Pos","Neg")
  Season.labs <- c("MAM (Austral Autumn)", "JJA (Austral Winter)", "SON (Austral Spring)"); names(Season.labs) <- c("MAM","JJA","SON")
  
  p <- ggplot()+facet_wrap(. ~ Dir, labeller= labeller(Dir = Phase.labs, Season = Season.labs))+
    # scale_y_continuous(position="right")+    geom_contour(aes(z = gh.z)) +
    geom_raster(data=Ppt.anom.p, aes(lon, lat, fill=Anomaly))+ scale_fill_stepsn(colours=paleta, breaks=at.ppt,
                                                                                           values=scales::rescale(at.ppt.v,from=range(at.ppt)),limits=c(min(at.ppt),max(at.ppt)),
                                                                                           guide=guide_colorsteps(even.steps = T,barwidth=unit(12,"cm")),
                                                                                 name = "Anom.\n[mm]")+ # 
    # scale_y_continuous(position="right")+
    geom_polygon(data=basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.5, show.legend = FALSE)+
    
    geom_rect(data=df.rect,inherit.aes = F,aes(xmin=xmi, xmax=xma, ymin=ymi, ymax=yma),color="black",linetype="dashed",linewidth=1, fill=NA, show.legend = FALSE)+
    
    geom_contour(data=subset(MDiv.anom.p,data.class==T), aes(lon,lat, z=Anomaly),color="red",binwidth =3,linewidth=0.35)+ #  linetype=1,
    geom_contour(data=subset(MDiv.anom.p,data.class==F), aes(lon,lat, z=Anomaly),color="blue",binwidth =3,linewidth=0.35)+
    
    geom_vector(data= subset(VIMF.anom.p, !p.value),aes(lon,lat,angle=Angle, mag=Anomaly), col="gray30",alpha=0.5,pivot=0.5, skip=1, show.legend = FALSE)+    
    geom_vector(data= subset(VIMF.anom.p, p.value),aes(lon,lat,angle=Angle, mag=Anomaly), col="purple4",pivot=0.5, skip=1, show.legend = FALSE)+

    scale_mag(max = 50,name="VIMF")+
    
    
    # scale_color_gradientn(colours = paleta)+ #brewer.pal(11,"Spectral")
    coord_fixed(xlim=xlim,ylim=ylim)+
    theme_bw()
  if( !is.na(title.p)) p <- p + labs(x="Longitude [°]",y="Latitude [°]",title= title.p) else p <- p+ labs(x="Longitude [°]",y="Latitude [°]")
  if (y.axis){
    p <- p + theme(legend.position=l.pos#, legend.background = element_rect(linetype="solid", colour ="black")) # c(0.9,0.2), legend.key.size = unit(1,"cm")
    )
  }else{
    p <- p + theme(legend.position=l.pos,#legend.background = element_rect(linetype="solid", colour ="black"), # , legend.key.size = unit(1,"cm")
                   axis.text.y= element_blank(), axis.title.y = element_blank()) # c(0.9,0.2)
  }
  
  if (x.axis){
    p <- p + theme(legend.position=l.pos#,legend.background = element_rect(linetype="solid", colour ="black")) # c(0.9,0.2), legend.key.size = unit(1,"cm")
    )
  }else{
    p <- p + theme(legend.position=l.pos,#legend.background = element_rect(linetype="solid", colour ="black"), # , legend.key.size = unit(1,"cm")
                   axis.text.x= element_blank(), axis.title.x = element_blank()) # c(0.9,0.2)
  }
  
  p <- p + theme(legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=15),
            strip.text = element_text( size=20),
            plot.title = element_text(size=20))
  # print(p)
  return(p)
}

VIMF.anom.p <- melt(VIMF.anom[["Atl3"]], id=c("lon","lat","Season","Angle","Anomaly")); colnames(VIMF.anom.p)[6] <- "Dir"; VIMF.anom.p$Dir <- factor(VIMF.anom.p$Dir, levels = c("Pos", "Neg"))
Test.comp.p <- melt(Test.comp[["Atl3"]], id=c("lon","lat","Season","p.value")); colnames(Test.comp.p)[5] <- "Dir"; Test.comp.p$Dir <- factor(Test.comp.p$Dir, levels = c("Pos", "Neg"))
VIMF.anom.p <- merge.data.frame(VIMF.anom.p, Test.comp.p)

rect.Atl3 <- data.frame(xmi=-65, xma=-57, ymi=0, yma=6,Season=as.factor("JJA"),
                        Dir=as.factor("Pos"))
P.Atl3 <- plotting.anom(VIMF.anom = VIMF.anom.p, MDiv.anom = MDiv.anom[["Atl3"]], Ppt.anom = Ppt.anom[["Atl3"]], df.rect =rect.Atl3,
                        Spe.season = c("JJA"), xlim=c(-82,-33), x.axis = F) #title.p = "Anomaly Composite (1980-2020), |Atl3| >= 1*SD",

P.Atl3 # 1000x480


VIMF.anom.p <- melt(VIMF.anom[["AMM"]], id=c("lon","lat","Season","Angle","Anomaly")); colnames(VIMF.anom.p)[6] <- "Dir"; VIMF.anom.p$Dir <- factor(VIMF.anom.p$Dir, levels = c("Pos", "Neg"))
Test.comp.p <- melt(Test.comp[["AMM"]], id=c("lon","lat","Season","p.value")); colnames(Test.comp.p)[5] <- "Dir"; Test.comp.p$Dir <- factor(Test.comp.p$Dir, levels = c("Pos", "Neg"))
VIMF.anom.p <- merge.data.frame(VIMF.anom.p, Test.comp.p)

rect.AMM <- data.frame(xmi=c(-70,-70,-73), 
                   xma=c(-58,-60,-67), 
                   ymi=c(3,-5,2), 
                   yma=c(10,3,7),
                   Season= factor(c("MAM","JJA","SON"), levels = c("MAM","JJA","SON")),
                   Dir=as.factor(rep("Pos",3)))

for ( i in c("MAM","JJA","SON")){
  P.AMM <- plotting.anom(VIMF.anom = VIMF.anom.p, MDiv.anom = MDiv.anom[["AMM"]], Ppt.anom = Ppt.anom[["AMM"]], df.rect =subset(rect.AMM, Season == i),
                         Spe.season = i,xlim=c(-82,-33), x.axis = F) # , title.p = "Anomaly Composite (1980-2020), |AMM| >= 1*SD"
  print(P.AMM)# 1200x430
}


library(gridExtra)
grid.arrange(P.AMM,P.Atl3, nrow=2, heights= c(2.0,1)) # 1000x1500
# grid.arrange(P.AMM,P.Atl3, nrow=1, widths= c(2.1,1))
# grid.arrange(P.AMM,P.Atl3, nrow=2)


## Pos - Neg ----
VIMF.anom <- list()
Pos.Neg <- function(x){
  y <- x$Pos[,c("lon","lat")]
  y <- cbind.data.frame(y, x$Pos[,seasons] - x$Neg[,seasons])
  return(y)
}

load("./01_VIMF/01_AnomComp_Atl3_VIMF_uv.RData");VIMF.anom[["Atl3"]] <- data.anom.uv %>% 
  
  lapply(., Pos.Neg) # calculate the Composite of Positive minus Negative Phase
rm(data.anom.uv)

load("./01_VIMF/01_AnomComp_AMM_VIMF_uv.RData"); VIMF.anom[["AMM"]] <- data.anom.uv %>% 
  
  lapply(., Pos.Neg) # calculate the Composite of Positive minus Negative Phase
rm(data.anom.uv)

sele.fun <- function(VIMF.anom, sele){
  VIMF.anom[["Mt.E.anom"]] <- VIMF.anom[["Mt.E.anom"]][ which(VIMF.anom[["Mt.E.anom"]]$lon %in% sele$lon) ,]; VIMF.anom[["Mt.E.anom"]] <- VIMF.anom[["Mt.E.anom"]][ which(VIMF.anom[["Mt.E.anom"]]$lat %in% sele$lat) ,]
  VIMF.anom[["Mt.N.anom"]] <- VIMF.anom[["Mt.N.anom"]][ which(VIMF.anom[["Mt.N.anom"]]$lon %in% sele$lon) ,]; VIMF.anom[["Mt.N.anom"]] <- VIMF.anom[["Mt.N.anom"]][ which(VIMF.anom[["Mt.N.anom"]]$lat %in% sele$lat) ,]
  return(VIMF.anom)
}
VIMF.anom <- lapply(VIMF.anom, sele.fun, sele)

# calculate the Angle and magnitude of the vector from the orthogonal components
VIMF.anom.Dif <- lapply(VIMF.anom, function(x){  
  
  Cord <- x$Mt.E.anom[,c("lon","lat")]
  
  Angle <- matrix(NA, ncol=4, nrow= nrow(x$Mt.E.anom)) %>% magrittr::set_colnames(., seasons)
  Magnitud <- matrix(NA, ncol=4, nrow= nrow(x$Mt.E.anom)) %>% magrittr::set_colnames(., seasons)
  
  for ( i in seasons){
    Angle[,i] <- atan2(dlat(x$Mt.N.anom[,i]), dlon(x$Mt.E.anom[,i], Cord[,"lat"]))*180/pi
    Magnitud[,i] <- Mag(x$Mt.N.anom[,i],x$Mt.E.anom[,i])
  }
  
  data.r <- list(Angle= Angle, Magnitude= Magnitud) %>% lapply(., cbind.data.frame, Cord)
  return(data.r)
}) %>% 
  lapply(., function(x){
    y <- lapply(x, melt, id=c("lon","lat"))
    y <- lapply(y, within, variable <- as.character(variable))
    return(y)
  } # and organize it for plotting
  ) %>% 
  lapply(., function(x){
    y <- merge.data.frame(x$Angle, x$Magnitude, by=c("lon","lat","variable")) %>% magrittr::set_colnames(., c("lon","lat","Season","Angle","Anomaly")) %>% 
      within(., Season <- factor(Season, levels= seasons))
    return(y)
    }
         )

# Postive minus Negative phase for Ppt
MDiv.anom.Dif <- lapply(MDiv.anom, function(x){
  y <- dplyr::group_by(x, lon,lat,Season) %>% 
    dplyr::summarise(., Anomaly[Dir=="Pos"] - Anomaly[Dir=="Neg"]) %>% 
    as.data.frame() %>% 
    magrittr::set_colnames(., c("lon","lat","Season","Anomaly")) %>% 
    within(., Season <- factor(Season, levels= seasons))
  return(y)
  })

# Postive minus Negative phase for MDiv
Ppt.anom.Dif <- lapply(Ppt.anom, function(x){
  y <- dplyr::group_by(x, lon,lat,Season) %>% 
    dplyr::summarise(., Anomaly[Dir=="Pos"] - Anomaly[Dir=="Neg"]) %>% 
    as.data.frame() %>% 
    magrittr::set_colnames(., c("lon","lat","Season","Anomaly")) %>% 
    within(., Season <- factor(Season, levels= seasons))
  return(y)
})

#### plotting Pos - Neg -----
plotting.anom2 <- function (VIMF.anom, MDiv.anom, Ppt.anom, df.rect, Spe.season, title.p, xlim=c(-82,-40), ylim=c(-20,13), l.pos="bottom", y.axis=TRUE, x.axis=TRUE, mag.max = 50){
  # data.anom must be the dataFrame, not the whole list
  at.ppt <- c(-250,-150,-100,-75,-50,-25,25,50,75,100,150,250)
  at.ppt.v <- (at.ppt[-length(at.ppt)] - at.ppt[-1])/2 + at.ppt[-1]
  col2alpha <- function(someColor, alpha=1){ newColor <- col2rgb(someColor); if(alpha <=1) alpha <- alpha * 255 else warning("Alpha should be between 0 and 1"); apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}
  paleta <- brewer.pal(11,"BrBG"); paleta <- col2alpha(paleta, alpha = 0.5)
  
  VIMF.anom.p <- VIMF.anom[sapply(Spe.season,function(s,d) which(d == s),VIMF.anom$Season),]
  MDiv.anom.p <- MDiv.anom[sapply(Spe.season,function(s,d) which(d == s),MDiv.anom$Season),]
  
  MDiv.anom.p$data.class <- MDiv.anom.p$Anomaly>0
  
  Ppt.anom.p <- Ppt.anom[sapply(Spe.season,function(s,d) which(d == s),Ppt.anom$Season),]
  
  Season.labs <- c("MAM (Austral Autumn)", "JJA (Austral Winter)", "SON (Austral Spring)"); names(Season.labs) <- c("MAM","JJA","SON")
  
  p <- ggplot()+facet_wrap(.~Season , strip.position = "left", ncol=1,labeller= labeller(Season = Season.labs))+
    # scale_y_continuous(position="right")+    geom_contour(aes(z = gh.z)) +
    geom_raster(data=Ppt.anom.p, aes(lon, lat, fill=Anomaly))+ scale_fill_stepsn(colours=paleta, breaks=at.ppt,
                                                                                 values=scales::rescale(at.ppt.v,from=range(at.ppt)),limits=c(min(at.ppt),max(at.ppt)),
                                                                                 guide=guide_colorsteps(even.steps = T,barwidth=unit(15,"cm")))+ # 
    scale_y_continuous(position="right")+
    geom_polygon(data=basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.5, show.legend = FALSE)+
    
    geom_rect(data=df.rect,inherit.aes = F,aes(xmin=xmi, xmax=xma, ymin=ymi, ymax=yma),color="black",linetype="dashed",linewidth=1, fill=NA, show.legend = FALSE)+
    
    geom_contour(data=subset(MDiv.anom.p,data.class==T), aes(lon,lat, z=Anomaly),color="red",binwidth =3,linewidth=0.35)+ #  linetype=1,
    geom_contour(data=subset(MDiv.anom.p,data.class==F), aes(lon,lat, z=Anomaly),color="blue",binwidth =3,linewidth=0.35)+
    
    # geom_vector(data= subset(VIMF.anom.p, !p.value),aes(lon,lat,angle=Angle, mag=Anomaly), col="gray30",alpha=0.5,pivot=0.5, skip=1, show.legend = FALSE)+    
    geom_vector(data= VIMF.anom.p,aes(lon,lat,angle=Angle, mag=Anomaly), col="purple4",pivot=0.5, skip=1, show.legend = FALSE)+
    
    scale_mag(max = mag.max,name="VIMF")+
    
    
    # scale_color_gradientn(colours = paleta)+ #brewer.pal(11,"Spectral")
    coord_fixed(xlim=xlim,ylim=ylim)+
    labs(x="Longitude",y="Latitude",title= title.p)+
    theme_bw()
  
  if (y.axis){
    p <- p + theme(legend.position=l.pos#, legend.background = element_rect(linetype="solid", colour ="black")) # c(0.9,0.2), legend.key.size = unit(1,"cm")
    )
  }else{
    p <- p + theme(legend.position=l.pos,#legend.background = element_rect(linetype="solid", colour ="black"), # , legend.key.size = unit(1,"cm")
                   axis.text.y= element_blank(), axis.title.y = element_blank()) # c(0.9,0.2)
  }
  
  if (x.axis){
    p <- p + theme(legend.position=l.pos#,legend.background = element_rect(linetype="solid", colour ="black")) # c(0.9,0.2), legend.key.size = unit(1,"cm")
    )
  }else{
    p <- p + theme(legend.position=l.pos,#legend.background = element_rect(linetype="solid", colour ="black"), # , legend.key.size = unit(1,"cm")
                   axis.text.x= element_blank(), axis.title.x = element_blank()) # c(0.9,0.2)
  }
  
  p <- p + theme(legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=15),
                 strip.text = element_text( size=20),
                 plot.title = element_text(size=20))
  # print(p)
  return(p)
}

rect.Atl3 <- data.frame(xmi=-65, xma=-57, ymi=0, yma=6,Season=as.factor("JJA"),
                        Dir=as.factor("Pos"))
P.Atl3 <- plotting.anom2(VIMF.anom = VIMF.anom.Dif[["Atl3"]], 
                        MDiv.anom = MDiv.anom.Dif[["Atl3"]],
                        Ppt.anom = Ppt.anom.Dif[["Atl3"]], df.rect = rect.Atl3,
                        Spe.season = c("JJA"), title.p = "Anomaly Composite (1980-2020), |Atl3| >= 1*SD",xlim=c(-82,-33), l.pos= "bottom")
P.Atl3


rect.AMM <- data.frame(xmi=c(-70,-70,-73), 
                       xma=c(-58,-60,-67), 
                       ymi=c(3,-5,2), 
                       yma=c(10,3,7),
                       Season= factor(c("MAM","JJA","SON"), levels = c("MAM","JJA","SON")),
                       Dir=as.factor(rep("Pos",3)))
P.AMM <- plotting.anom2(VIMF.anom = VIMF.anom.Dif[["AMM"]], 
                       MDiv.anom = MDiv.anom.Dif[["AMM"]],
                       Ppt.anom = Ppt.anom.Dif[["AMM"]], df.rect =rect.AMM,
                       Spe.season = c("MAM","JJA","SON"), title.p = "Anomaly Composite (1980-2020), |AMM| >= 1*SD",xlim=c(-82,-33), l.pos= "none", x.axis = FALSE, mag.max=75)
P.AMM

# grid.arrange(P.AMM,P.Atl3, nrow=2, heights= c(2.0,1)) # 1000x1500
