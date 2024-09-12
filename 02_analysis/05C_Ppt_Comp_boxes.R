rm(list = ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
"%>%" = magrittr::`%>%`

# Dates ----
Dates.ms <- seq(as.Date("1979-02-01"),as.Date("2020-11-30"), by="month")
Dates.chi <- seq(as.Date("1981-01-01"),as.Date("2023-08-31"), by="month")
Dates.a <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month")

seasons <- c("DJF","MAM","JJA","SON")

years.m <- format(Dates.a, format="%Y"); Years <- format(Dates.a, format="%Y") %>% unique() %>% .[-1]
Season.y <- paste0(years.m[-1],"-",time2season(Dates.a)[-length(Dates.a)]) %>% c(.,.[length(.)])
Year.s <- unique(Season.y)

years.m.chi <- format(Dates.chi, format="%Y"); Years.chi <- format(Dates.chi, format="%Y") %>% unique() %>% .[-1]
Season.y.chi <- paste0(years.m.chi[-1],"-",time2season(Dates.chi)[-length(Dates.chi)]) %>% c(. , .[length( .)])
Year.s.chi <- unique(Season.y.chi)

years.m.ms <- format(Dates.ms, format="%Y"); Years.ms <- format(Dates.ms, format="%Y") %>% unique() %>% .[-1]
Season.y.ms <- paste0(years.m.ms[-1],"-",time2season(Dates.ms)[-length(Dates.ms)]) %>% c(. , .[length( .)])
Year.s.ms <- unique(Season.y.ms)

## load data ----
# cdo seassum MSWEP_Ppt_1979_2020_cdo.nc MSWEP_Ppt_1979_2020_seasonal.nc
MS.s <- brick("./01_data/02_Ppt/MSWEP_Ppt_1979_2020_seasonal.nc")
CHIRPS.s <- brick("./01_data/02_Ppt/CHIRPS_TropSAm_seasonal_2023.nc")

# sst standardize indices
AMM.seasons <- read.csv("./01_data/01_SSTindices/AMM_std_1980-2020.csv"); row.names(AMM.seasons) <- AMM.seasons[,1]; AMM.seasons <- AMM.seasons[,-1]
Atl3.seasons <- read.csv("./01_data/01_SSTindices/Atl3_std_1980-2020.csv"); row.names(Atl3.seasons) <- Atl3.seasons[,1]; Atl3.seasons <- Atl3.seasons[,-1]
ELI.seasons <- read.csv("./01_data/01_SSTindices/ELI_std_1980-2020.csv"); row.names(ELI.seasons) <- ELI.seasons[,1]; ELI.seasons <- ELI.seasons[,-1]

## crop 1980 - onwards ----
MS.s <- MS.s[[match(Year.s, Year.s.ms)]]
CHIRPS.s <- CHIRPS.s[[ match(Year.s, Year.s.chi) %>% .[!is.na( .)]]]
Year.s.chi <- Year.s.chi[match(Year.s, Year.s.chi)] %>% .[!is.na( .)]

### crop to boxes ----
CHIRPS.box <- MS.box <- list()

rect <- data.frame(xmin=c(-70,-70,-73,-65), 
                   xmax=c(-58,-60,-67,-57), 
                   ymin=c(3,-5,2,0), 
                   ymax=c(10,3,7,6),
                   Season= factor(c("MAM","JJA","SON","JJA"), levels = c("MAM","JJA","SON")),
                   mode=c("AMM","AMM","AMM","Atl3"))

for(i in 1:4){
  
  MS.box[[i]] <- crop(MS.s, extent( rect[i, 1:4] %>% as.numeric())) %>% 
    cellStats(., stat="mean")
  CHIRPS.box[[i]] <- crop(CHIRPS.s, extent( rect[i, 1:4] %>% as.numeric()) ) %>% 
    cellStats(., stat="mean")
}
CHIRPS.box <- lapply(CHIRPS.box, function(x){ y <- c(rep(NA,4), x); return(y)}) # adding NA for 1980 in CHIRPS

### seasonal separation and anomalies ----
MS.box.s <- CHIRPS.box.s <- list()
for( i in 1:4){
  MS.box.s[[i]] <- matrix(NA,nrow=length(Years), ncol=4); colnames(MS.box.s[[i]]) <- seasons
  CHIRPS.box.s[[i]] <- matrix(NA,nrow=length(Years), ncol=4); colnames(CHIRPS.box.s[[i]]) <- seasons
  
  for(j in seasons){
    MS.box.s[[i]][1:length(Years) ,j] <- MS.box[[i]][substr(Year.s,6,8) == j]
    CHIRPS.box.s[[i]][1:length(Years) ,j] <- CHIRPS.box[[i]][substr(Year.s,6,8) == j]
  }
}
MS.box.mean <- sapply(MS.box.s, function(x) colMeans(x, na.rm = T)) %>% t()
CHIRPS.box.mean <- sapply(CHIRPS.box.s, function(x) colMeans(x, na.rm = T)) %>% t()
Box.anom <- list(); Box.anom[["MSWEP"]] <- Box.anom[["CHIRPS"]] <- list()

for( i in 1:4){
  Box.anom[["MSWEP"]][[i]] <- sweep(MS.box.s[[i]], 2, MS.box.mean[i, ], FUN="-")
  Box.anom[["CHIRPS"]][[i]] <- sweep(CHIRPS.box.s[[i]], 2, CHIRPS.box.mean[i, ], FUN="-")
}

## identifying positive and negative periods ----
AMM.bool <- Atl3.bool <- list()
#---------------------------------------------------------------------- identification of periods with SST higher than or lower than 1*SD
{
  AMM.bool[["Pos"]] <- data.frame(matrix(FALSE, ncol=4, nrow= length(Years))); colnames(AMM.bool[["Pos"]]) <- seasons; rownames(AMM.bool[["Pos"]]) <- Years
  AMM.bool[["Pos"]][which(AMM.seasons >=1, arr.ind = T)] <- TRUE
  AMM.bool[["Neg"]] <- data.frame(matrix(FALSE, ncol=4, nrow= length(Years))); colnames(AMM.bool[["Neg"]]) <- seasons; rownames(AMM.bool[["Neg"]]) <- Years
  AMM.bool[["Neg"]][which(AMM.seasons <=-1, arr.ind = T)] <- TRUE
  
  Atl3.bool[["Pos"]] <- data.frame(matrix(FALSE, ncol=4, nrow= length(Years))); colnames(Atl3.bool[["Pos"]]) <- seasons; rownames(Atl3.bool[["Pos"]]) <- Years
  Atl3.bool[["Pos"]][which(Atl3.seasons >=1, arr.ind = T)] <- TRUE
  Atl3.bool[["Neg"]] <- data.frame(matrix(FALSE, ncol=4, nrow= length(Years))); colnames(Atl3.bool[["Neg"]]) <- seasons; rownames(Atl3.bool[["Neg"]]) <- Years
  Atl3.bool[["Neg"]][which(Atl3.seasons <=-1, arr.ind = T)] <- TRUE
}
#--------------------------------------------------------------------------- initialization of datasets
{
  P.comp <- list(); P.comp[["MSWEP"]] <- P.comp[["CHIRPS"]] <- list()
  P.comp <- lapply(P.comp, function(x){ x[["Pos"]] <- x[["Neg"]] <- x[["Neutral"]] <- matrix(NA, nrow=4, ncol=4); x <- lapply(x, 'colnames<-',seasons); return(x)})
  
  P.comp.se <- list(); P.comp.se[["MSWEP"]] <- P.comp.se[["CHIRPS"]] <- list()
  P.comp.se <- lapply(P.comp.se, function(x){ x[["Pos"]] <- x[["Neg"]] <- x[["Neutral"]] <- matrix(NA, nrow=4, ncol=4); x <- lapply(x, 'colnames<-',seasons); return(x)})
  
  Test.comp <- list(); Test.comp[["MSWEP"]] <- Test.comp[["CHIRPS"]] <- list()
  Test.comp <- lapply(Test.comp, function(x){ x[["Pos"]] <- x[["Neg"]] <- matrix(NA, nrow=4, ncol=4); x <- lapply(x, 'colnames<-',seasons); return(x)})
}

# composites and testing
for ( i in seasons){
  for (j in c("MSWEP","CHIRPS")){
    for(k in 1:4){
      # preselection SST phase and the actual Ppt fot all the modes
      sst.pos <- get( paste0(rect[k,"mode"],".bool") )[["Pos"]][,i]; sst.neg <- get( paste0(rect[k,"mode"],".bool") )[["Neg"]][,i]
      Pos <- Box.anom[[ j ]][[ k ]][ sst.pos ,i]; Neg <- Box.anom[[ j ]][[ k ]][ sst.neg ,i]
      Neutral <- Box.anom[[ j ]][[ k ]][ !(sst.pos | sst.neg) ,i]
      
      P.comp[[j]][["Pos"]][k ,i] <- mean(Pos, na.rm=T); P.comp[[j]][["Neg"]][k ,i] <- mean(Neg, na.rm=T)
      P.comp[[j]][["Neutral"]][k ,i] <- mean(Neutral, na.rm=T)
      P.comp.se[[j]][["Pos"]][k ,i] <- sd(Pos, na.rm=T)/sqrt(length(Pos)) ; P.comp.se[[j]][["Neg"]][k ,i] <- sd(Neg, na.rm=T)/sqrt(length(Neg))
      P.comp.se[[j]][["Neutral"]][k ,i] <- sd(Neutral, na.rm=T)/sqrt(length(Neutral))
      
      Test.comp[[j]][["Pos"]][[k,i]] <- wilcox.test(Pos, Neutral)$p.value; Test.comp[[j]][["Neg"]][[k,i]] <- wilcox.test(Neg, Neutral)$p.value
    }
  }
}

# organize data for plotting ----

rect.comp <- list(); rect.Sig <- list(); rect.SE <- list()

rect.comp[["MSWEP"]] <- rect.comp[["CHIRPS"]] <- rect[,c("Season","mode")] # adding columns of positive and negative phase 
rect.comp <- lapply( rect.comp, function(x){ x$Pos <- NA; x$Neg <- NA; x$Neutral <- NA;  return(x)})
rect.SE[["MSWEP"]] <- rect.SE[["CHIRPS"]] <- rect[,c("Season","mode")]
rect.SE <- lapply( rect.SE, function(x){ x$Pos <- NA; x$Neg <- NA;  x$Neutral; return(x)})
rect.Sig[["MSWEP"]] <- rect.Sig[["CHIRPS"]] <- rect[,c("Season","mode")]
rect.Sig <- lapply( rect.Sig, function(x){ x$Pos <- NA; x$Neg <- NA;  return(x)})

rect$Season <- as.character(rect$Season)

for(i in 1:4){
  for(j in c("MSWEP","CHIRPS")){
    rect.comp[[j]][i,"Pos"] <- P.comp[[j]][["Pos"]][i, rect[i,"Season"]]; rect.comp[[j]][i,"Neg"] <- P.comp[[j]][["Neg"]][i, rect[i,"Season"]]
    rect.comp[[j]][i,"Neutral"] <- P.comp[[j]][["Neutral"]][i, rect[i,"Season"] ]
    
    rect.Sig[[j]][i,"Pos"] <- Test.comp[[j]][["Pos"]][i, rect[i,"Season"]]; rect.Sig[[j]][i,"Neg"] <- Test.comp[[j]][["Neg"]][i, rect[i,"Season"]]
    rect.SE[[j]][i,"Pos"] <- P.comp.se[[j]][["Pos"]][i, rect[i,"Season"]]; rect.SE[[j]][i,"Neg"] <- P.comp.se[[j]][["Neg"]][i, rect[i,"Season"]]
    rect.SE[[j]][i,"Neutral"] <- P.comp.se[[j]][["Neutral"]][i, rect[i,"Season"]]
  }
}
rect.comp <- melt(rect.comp, id=c("Season","mode"))

# reclasify data
rect.Sig2 <- melt(rect.Sig, id=c("Season","mode")) %>%
  dplyr::mutate(value= round(value,3)) %>% 
  dplyr::mutate(., Sig = dplyr::case_when(value <=0.05 ~ "p<0.05", value <=0.1 ~ "p<0.1",.default = NA), .keep="unused")

rect.comp <- merge.data.frame(rect.comp, rect.Sig2, by=c("Season","mode","variable","L1"), all.x = TRUE) %>% 
  within(.,{
  variable <- factor(variable, levels=c("Pos","Neutral", "Neg"))
  L1 <- factor(L1, levels = c("MSWEP","CHIRPS"))
})

rect.SE <- melt(rect.SE, id=c("Season","mode")) %>% 
  within(., {
    variable <- factor(variable, levels=c("Pos","Neutral", "Neg"))
    L1 <- factor(L1, levels = c("MSWEP","CHIRPS"))
  }) %>% 
  merge.data.frame(rect.comp, ., by=c("Season","mode","variable","L1")) %>% 
  dplyr::mutate(., ymi= value.x - value.y, yma= value.x + value.y) %>% 
  within(., { Sig <- factor(Sig)# Sig[is.na(Sig)] <- "None"
  })

mess <- apply(rect[,1:4],1,function(x){paste0("[",x[1]*(-1),"°W,",x[2]*(-1),"°W]x\n[",x[3],"°N,",x[4],"°N]")})
mess <- subset(rect.SE, L1=="MSWEP" & variable=="Neg") %>% dplyr::select(., value.x, variable, mode, Season) %>% dplyr::mutate(variable="Pos") %>% #dplyr::mutate(value.x=value.x*(-1)) %>% 
  .[c(3,1,4,2),] %>% cbind(.,mess)

# actual plotting ----
Season.labs <- c("MAM (Austral Autumn)", "JJA (Austral Winter)", "SON (Austral Spring)"); names(Season.labs) <- c("MAM","JJA","SON")

plot.b <- function(Rect.P, Mess){
  p <- ggplot(data=Rect.P)+
    geom_bar(aes(x=L1, y=value.x, fill=variable),
             stat="identity", color="black", position=position_dodge())+
    scale_fill_manual(values = c("#FF0000","grey","#0000FF"), name="Phase")+
    geom_point(aes(x=L1,y=value.x/2, shape=Sig, group=variable),size=3, position=position_dodge(.9), col="gold",stroke=2)+
    scale_shape_manual(values = c(8,1),na.translate = F,name="Significance") +
    geom_errorbar(aes(x=L1,ymin=value.x - value.y, ymax=value.x + value.y, fill=variable), width=.2, position=position_dodge(.9))+
    geom_text(data = Mess, aes(x="MSWEP", y=value.x, fill=variable, label=mess),hjust=1)+
    # ylim(c(-150,150))+
    labs(x="Dataset",y="Ppt anomaly [mm]")+
    theme_bw()+theme(legend.position ="bottom",legend.spacing = unit(3,"points"), legend.box.margin = margin(l=-58,unit="points"), legend.background =element_rect(linewidth = 0.5, linetype="solid", colour ="black") ,
                     axis.text.x = element_text(size=20),axis.title.x = element_blank(),
                     axis.title.y = element_text(size=20), axis.text.y = element_text(size=19, angle=90, hjust = 0.5),
                     strip.text = element_text(size=20),
                     legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=15), 
                     plot.title = element_text(size=20))
    
  return(p)
}


for ( i in c("MAM","JJA","SON")){
  plot.b(Rect.P= subset(rect.SE, mode=="AMM" & Season==i), Mess= subset(mess, mode=="AMM" & Season==i )) %>% print()
}
plot.b(Rect.P= subset(rect.SE,mode=="Atl3"), Mess= subset(mess,mode=="Atl3"))
