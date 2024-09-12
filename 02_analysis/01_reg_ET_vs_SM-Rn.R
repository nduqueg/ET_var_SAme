rm(list=ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
"%>%"=magrittr::`%>%`

# Dates ----
Dm.sm <- seq(as.Date("1950-01-01"),as.Date("2020-12-31"), by="month")
Dates <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"), by="month")
Dm.e5 <- seq(as.Date("1950-01-01"),as.Date("2020-12-31"),by="month")

years.m <- format(Dates, format="%Y"); Years <- unique(format(Dates, format="%Y")); Years <- Years[-1]
Season.y <- paste0(years.m[-1],"-",time2season(Dates)[-length(Dates)]); Season.y <- c(Season.y, Season.y[length(Season.y)])
Year.s <- unique(Season.y)
seasons <- c("DJF","MAM","JJA","SON")

## load data ----
ET <- brick("./01_data/04_Evap/ERA5L_TEvap_1950_2020.nc")
Rn <- brick("./01_data/05_Rad/ERA5L_NetRad_1950-2020.nc")
SM <- brick("./01_data/03_SM/ERA5L_SM_1L_1950-2020.nc")

masks <- list()
masks[["E5"]] <- ET[[1]]*0+1
masks[["SM"]] <- SM[[1]]*0+1

## crop 1980 -onwards, PPt-resampling  & seasonal pre-processing ----
ET <- ET[[match(Dates,Dm.e5)]];       ET.s <- stackApply(ET, Season.y, fun=mean); ET.s <- ET.s * masks[["E5"]]
Rn <- Rn[[match(Dates,Dm.e5)]];       Rn.s <- stackApply(Rn, Season.y, fun=mean); Rn.s <- Rn.s * masks[["E5"]]
SM <- SM[[match(Dates,Dm.sm)]];         SM.s <- stackApply(SM, Season.y, fun=sum)

# locate both datasets at the same location
SM.s <- mask(resample(aggregate(SM.s, fact=2, fun="mean"), ET.s[[1]], method="ngb"),ET.s[[1]])
SM.s.t <- rasterToPoints(SM.s)
ET.s.t <- rasterToPoints(ET.s)
Rn.s.t <- rasterToPoints(Rn.s)

cells.ET <- ET.s.t[,c(1,2)]; ET.s.t <- ET.s.t[,-c(1,2)]
cells.sm <- SM.s.t[,c(1,2)]; SM.s.t <- SM.s.t[,-c(1,2)]
Rn.s.t <- Rn.s.t[,-c(1,2)]

ET.seasons <- list(); SM.seasons <- list(); Rn.seasons <- list()
for ( i in seasons){
  ET.seasons[[i]] <- as.data.frame(matrix(NA, ncol=nrow(ET.s.t), nrow=length(Years)))
  SM.seasons[[i]] <- as.data.frame(matrix(NA, ncol=nrow(SM.s.t), nrow=length(Years)))
  Rn.seasons[[i]] <- as.data.frame(matrix(NA, ncol=nrow(Rn.s.t), nrow=length(Years)))
}

for( i in seasons){
  ET.seasons[[ i ]] <- t(ET.s.t[, which(substr(Year.s,6,8) == i)])
  Rn.seasons[[ i ]] <- t(Rn.s.t[, which(substr(Year.s,6,8) == i)])
  SM.seasons[[ i ]] <- t(SM.s.t[, which(substr(Year.s,6,8) == i)])
}

# standarization of preditors (Rn & SM)
ET.seasons <- lapply(ET.seasons, function(x) apply(x, 2, scale, scale=T))
Rn.seasons <- lapply(Rn.seasons, function(x) apply(x, 2, scale, scale=T))
SM.seasons <- lapply(SM.seasons, function(x) apply(x, 2, scale, scale=T))
SM.seasons <- lapply(SM.seasons, function(a){ a[, is.na(a[1,])] <- 0; return(a)})

save(Rn.seasons,SM.seasons,ET.seasons,file="./02_analysis/02_standardize_predictors.RData")
load("./02_analysis/02_standardize_predictors.RData")
#### regressions by season ####
# lists
Reg.R2 <- list(); Reg.R2.adj <- list();Reg.Coef <- list();Reg.Coef.sig <- list()

library(parallel)
library(doParallel)
library(tcltk)

for (i in seasons){
  print(i)
  
  numCores <- 20
  
  # cl <- makeSOCKcluster(numCores)
  doParallel::registerDoParallel(cores= numCores)
  
  # pb <- txtProgressBar(min=1, max= 40, style=3)
  # progress <- function(n) setTxtProgressBar(pb, n)
  # opts <- list(progress=progress)
  
  
  # results <- foreach ( j = 1:ncol(ET.seasons[[i]]), .options.snow=opts) %dopar% {
  results <- foreach ( j = 1:ncol(ET.seasons[[i]])) %dopar% {
    # setTxtProgressBar(pb,j)
    predictors <- cbind(SM.seasons[[i]][,j], Rn.seasons[[i]][,j])
    aux <- lm(ET.seasons[[i]][,j] ~ predictors)
    return(aux)
  }
  stopImplicitCluster()
  
  Reg.summary <- lapply(results,function(x) summary(x))
  
  print("extracting R2")
  # Reg.R2[[i]] <- unlist(lapply(Reg.s,function(x) summary(x)$r.squared))
  Reg.R2[[i]] <- unlist(lapply(Reg.summary,function(x) x$r.squared))
  Reg.R2[[i]] <- cbind.data.frame(cells.ET, Reg.R2[[i]])
  
  print("extracting adjusted R2")
  Reg.R2.adj[[i]] <- unlist(lapply(Reg.summary,function(x) x$adj.r.squared))
  Reg.R2.adj[[i]] <- cbind(cells.ET, as.data.frame(Reg.R2.adj[[i]]))
  
  print("extracting coefficients")
  Reg.Coef[[i]] <- as.data.frame(lapply(results,function(x) x$coefficients))
  aux2 <- sapply(Reg.summary,function(x) as.numeric(x$coefficients[,4])) # in the cells that there's no SM variability (i.e. always saturated) 
  ind <- which(lapply(aux2, function(x) length(x)) <= 2) # we need to filter the p.value because it is NA, and when accessing it,
  for(j in ind) aux2[[j]] <- c(aux2[[j]][1], NA, aux2[[j]][2]) # it returns just 2 values
  Reg.Coef.sig[[i]] <- as.data.frame(aux2)
  
  rm(results, Reg.summary)
} 

Reg.R2 <- lapply(Reg.R2, `colnames<-`, c("lon","lat","R2"))
Reg.R2.adj <- lapply(Reg.R2.adj, `colnames<-`, c("lon","lat","R2"))

save(Reg.R2,"Reg.R2",file="./02_analysis/02_regR2_ET_SM_Rn.RData")
save(Reg.R2.adj,"Reg.R2.adj",file="./02_analysis/02_regR2Adj_ET_SM_Rn.RData")
save(Reg.Coef, "Reg.Coef", file="./02_analysis/02_Coef_ET_SM_Rn.RData")
save(Reg.Coef.sig, "Reg.Coef.sig", file="./02_analysis/02_Coef_Significance_ET_SM_Rn.RData")



#### plotting the map ####
load("./02_analysis/02_regR2_ET_SM_Rn.RData")
load("./02_analysis/02_regR2Adj_ET_SM_Rn.RData")
load("./02_analysis/02_Coef_ET_SM_Rn.RData")
load("./02_analysis/02_Coef_Significance_ET_SM_Rn.RData")
seasons <- c("DJF","MAM","JJA","SON")

Basins <- shapefile("./01_data/hybas_sa_lev03_v1c.shp")
SA <- shapefile("./01_data/South_America.shp")

Reg.R2.t <- melt(Reg.R2, id=c("lon","lat")) %>% within(., L1 <- factor(L1, levels = c("MAM","JJA","SON","DJF")))
Reg.R2.adj.t <- melt(Reg.R2.adj, id=c("lon","lat")) %>% within(., L1 <- factor(L1, levels = c("MAM","JJA","SON","DJF")))

cells.ET <- Reg.R2$DJF[,1:2]

Reg.Coef.t <- lapply(Reg.Coef, function(x,y){
  z <- cbind(y,as.data.frame(t(x))); colnames(z)[1:2] <- c("lon","lat"); return(z)
  }, cells.ET) %>% 
  melt(.,id=c("lon","lat")) %>% within(., L1 <- factor(L1, levels = c("MAM","JJA","SON","DJF")))

Reg.Coef.class <- list()
for (i in c("MAM","JJA","SON","DJF")){
  
  aux.sig <- t(Reg.Coef.sig[[i]])[,-1]
  No.Sig <- apply(unname(aux.sig),1,function(x) x > 0.05) %>% t() # No.Sig <- which(aux.sig > 0.05, arr.ind = T) ;   ind.Neither.Sig <- No.Sig[duplicated(No.Sig[,1]), 1]
    
  # bot columns need to be No Significant at the same cell so as to be No significant in that cell
  ind.Neither.Sig <- apply(No.Sig,1, function(x) x[1] & x[2]) %>% which()
  
  aux <- t(Reg.Coef[[i]])[,-1] # do the classification
  class <- apply(abs(aux), 1, which.max)
  
  for (j in 1:nrow(No.Sig)){ # testing that haven't identify the maximum driver wrongly
    if(sum(No.Sig[j,], na.rm=T)==1){
      if(which(No.Sig[j,]) == class[j]) print(paste("cell",j,"is max but not significant"))
    }
  }
  
  class[ind.Neither.Sig] <- 3 # assigning another classification for those No Significant
  Reg.Coef.class[[i]] <- class
}


Reg.Coef.class.t <- lapply(Reg.Coef.class, function(x,y){
  z <- cbind.data.frame(y,x); colnames(z)[1:2] <- c("lon","lat"); return(z)
  }, cells.ET) %>% 
  melt(.,id=c("lon","lat")) %>% 
  within(.,{
    L1 <- factor(L1, levels = c("MAM","JJA","SON","DJF"))
    levels(variable)  <- c("Ranking")
    value <- factor(value); levels(value) <- c("Soil Moisture\n(water-limited)", "Net radiation\n(energy-limited)", "Neither Significant")
  })



at.m <- seq(0,1,length.out = 11); at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]

ggplot(Reg.R2.t)+
  geom_raster(aes(x=lon,y=lat,fill=value))+facet_wrap(.~L1, ncol=4)+
  scale_fill_stepsn(colours=brewer.pal(10,"Spectral"), breaks=at.m,
                    limits=c(min(at.m),max(at.m)), # values=c(0:12)/12, equally spaced
                    guide=guide_colorsteps(barwidth=unit(15,"cm")))+ # even.steps = F, 
  #geom_sf(data=world,fill="transparent")+
  scale_y_continuous(position="right")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  coord_fixed(xlim=c(-82,-40),ylim=c(-20,15))+labs(x="Long",y="Lat",title="Coef. determination ET vs SM and Rn (1980-2020)")+
  theme_bw()+theme(legend.position = "bottom")

p2 <- ggplot(Reg.R2.adj.t)+
  geom_raster(aes(x=lon,y=lat,fill=value))+facet_wrap(.~L1, nrow=1)+
  scale_fill_stepsn(colours=brewer.pal(10,"Spectral"), breaks=at.m,
                    limits=c(min(at.m),max(at.m)), # values=c(0:11)/11,
                    guide=guide_colorsteps( barwidth=unit(10,"cm")), name = "R2")+ # even.steps = F,
  #geom_sf(data=world,fill="transparent")+
  scale_y_continuous(position="right")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  coord_fixed(xlim=c(-82,-33),ylim=c(-20,13))+labs(x="Longitude [°]",y="Latitude [°]")+ # ,title="Adjusted Coef. determination - ET vs SM and Rn (1980-2020)"
  theme_bw() +theme(legend.position = "bottom")


range <- with(Reg.Coef.t, max( abs(min(value,na.rm=T)), max(value, na.rm=T)))
at.m <- seq(-range,range,length.out = 11) %>% round(.,2); at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]
levels(Reg.Coef.t$variable) <-  c("(Intercept)","Soil Moisture","Net Radiation")
data.coef <- subset(Reg.Coef.t, variable!="(Intercept)")
p1 <- ggplot(data.coef)+
  facet_grid(variable~L1, switch = "y")+
  geom_raster(aes(x=lon,y=lat,fill=value))+
  scale_fill_stepsn(colours=brewer.pal(10,"PiYG"), breaks=at.m,
                    limits=c(min(at.m),max(at.m)), # values=c(0:11)/11,, equally spaced
                    guide=guide_colorsteps(even.steps = F, barwidth=unit(10,"cm")))+
  #geom_sf(data=world,fill="transparent")+
  scale_y_continuous(position="right")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  coord_fixed(xlim=c(-82,-35),ylim=c(-20,13))+labs(x="Longitude [°]",y="Latitude [°]",title="Coefficients of regresion - ET vs SM and Rn (1980-2020)")+
  theme_bw()+theme(legend.position = "bottom")

p3 <- ggplot(Reg.Coef.class.t)+
  facet_wrap(.~L1, ncol=4)+
  geom_raster(aes(x=lon,y=lat,fill=value))+
  scale_fill_manual(values=c("#377eb8","#e41a1c","#ffff33"), name="Driver")+
  #geom_sf(data=world,fill="transparent")+
  scale_y_continuous(position="right")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  coord_fixed(xlim=c(-82,-35),ylim=c(-20,13))+labs(x="Longitude [°]",y="Latitude [°]")+ # ,title="Main driver of ET - Highest absolute coefficient value regression standardized ET vs SM and Rn"
  theme_bw()+theme(legend.position = "bottom")

ggplot(subset(Reg.Coef.class.t, L1!="DJF"))+
  facet_grid(variable~L1, switch = "y")+
  geom_raster(aes(x=lon,y=lat,fill=value))+
  scale_fill_manual(values=c("#377eb8","#e41a1c","#ffff33"), name="Driver")+
  #geom_sf(data=world,fill="transparent")+
  scale_y_continuous(position="right")+
  geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
  geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed", colour="black",fill="NA",linewidth=0.05)+
  coord_fixed(xlim=c(-82,-35),ylim=c(-20,13))+labs(x="Longitude [°]",y="Latitude [°]",title="Main driver of ET - Highest absolute coefficient value regression standardized ET vs SM and Rn")+
  theme_bw()+theme(legend.position = "bottom",legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=15),
                   strip.text = element_text(size=20), plot.title = element_text(size=20)) # 1600x550


library(gridExtra)
grid.arrange(p1,p2,nrow=2, heights= c(1.8,1)) #1600 x 1200

grid.arrange(p2,p3,nrow=2, heights= c(1,1)) #1000 x 600
