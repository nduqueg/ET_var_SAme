rm(list=ls())
cat("\014")

library(hydroTSM)
library(ggplot2)
library(ggrepel)
library(reshape2)

# Dates ----
Dm.h <- seq(as.Date("1850-01-01"),as.Date("2022-11-30"),by="month")
Dm.e <- seq(as.Date("1854-01-01"),as.Date("2020-12-31"),by="month")
Dm.ana <- seq(as.Date("1979-12-01"),as.Date("2020-11-30"),by="month")
DmM.ana <- unique(format(Dm.ana, format="%m"))

years.m <- format(Dm.ana, format="%Y"); Years <- unique(format(Dm.ana, format="%Y")); Years <- Years[-1]
Season.y <- paste0(years.m[-1],"-",time2season(Dm.ana)[-length(Dm.ana)]); Season.y <- c(Season.y, Season.y[length(Season.y)])
Year.s <- unique(Season.y)

seasons <- c("DJF","MAM","JJA","SON")

# load data ----

AEN <- read.zoo(read.csv("00A_AtlEN_SST.csv"), index.column = 1); AEN <- AEN[Dm.ana]
AMM <- read.zoo(read.csv("00A_AMM_deltaSST.csv"), index.column = 1); AMM <- AMM[Dm.ana]
ELI <- read.zoo(read.csv("00A_ELI_SST.csv"), index.column = 1)#; ELI <- rollmean(ELI, k=3)
ELI <- ELI[Dm.ana]
Nino34 <- read.zoo(read.csv("00A_N34_SST.csv"), index.column = 1); Nino34 <- Nino34[Dm.ana]

TNA <- read.zoo(read.csv("00A_TNA_SST.csv"), index.column = 1); TNA <- TNA[Dm.ana]
TSA <- read.zoo(read.csv("00A_TSA_SST.csv"), index.column = 1); TSA <- TSA[Dm.ana]

### seasonal distribution ----
AEN.s <- aggregate(AEN[,"E.sst.AEN"], by=list(Season.y), FUN=mean)
AMM.s <- aggregate(AMM[,"E.sst.AMM"], by=list(Season.y), FUN=mean)
ELI.s <- aggregate(ELI, by=list(Season.y), FUN=mean)

TNA.s <- aggregate(TNA[,"E.sst.TNA"], by=list(Season.y), FUN=mean)
TSA.s <- aggregate(TSA[,"E.sst.TSA"], by=list(Season.y), FUN=mean)

### arrange for plotting ----
AEN.d <- list(); AMM.d <- list(); ELI.d <- list(); TNA.d <- list(); TSA.d <- list() # for density
AEN.p <- list(); AMM.p <- list(); ELI.p <- list(); TNA.p <- list(); TSA.p <- list() # for points
for (i in seasons){
  AEN.d[[i]] <- density(AEN.s[substr(index(AEN.s),6,8)== i])
  AMM.d[[i]] <- density(AMM.s[substr(index(AMM.s),6,8)== i])
  ELI.d[[i]] <- density(ELI.s[substr(index(ELI.s),6,8)== i])
  TNA.d[[i]] <- density(TNA.s[substr(index(TNA.s),6,8)== i])
  TSA.d[[i]] <- density(TSA.s[substr(index(TSA.s),6,8)== i])
  
  AEN.d[[i]] <- cbind(as.data.frame(AEN.d[[i]][c("x","y")]),i)
  AMM.d[[i]] <- cbind(as.data.frame(AMM.d[[i]][c("x","y")]),i)
  ELI.d[[i]] <- cbind(as.data.frame(ELI.d[[i]][c("x","y")]),i)
  TNA.d[[i]] <- cbind(as.data.frame(TNA.d[[i]][c("x","y")]),i)
  TSA.d[[i]] <- cbind(as.data.frame(TSA.d[[i]][c("x","y")]),i)
  
  AEN.p[[i]] <- cbind(as.data.frame(AEN.s[substr(index(AEN.s),6,8)== i]),0,i)
  AMM.p[[i]] <- cbind(as.data.frame(AMM.s[substr(index(AMM.s),6,8)== i]),0,i)
  ELI.p[[i]] <- cbind(as.data.frame(ELI.s[substr(index(ELI.s),6,8)== i]),0,i)
  TNA.p[[i]] <- cbind(as.data.frame(TNA.s[substr(index(TNA.s),6,8)== i]),0,i)
  TSA.p[[i]] <- cbind(as.data.frame(TSA.s[substr(index(TSA.s),6,8)== i]),0,i)
}
data <- list()
data[["AEN"]] <- do.call("rbind",AEN.d)
data[["AMM"]] <- do.call("rbind",AMM.d)
data[["ELI"]] <- do.call("rbind",ELI.d)
data[["TNA"]] <- do.call("rbind",TNA.d)
data[["TSA"]] <- do.call("rbind",TSA.d)

Type <- "line"
data <- lapply(data, cbind, Type)

# creating points to add to the plotting
data.p <- list()
data.p[["AEN"]] <- do.call("rbind",AEN.p)
data.p[["AMM"]] <- do.call("rbind",AMM.p)
data.p[["ELI"]] <- do.call("rbind",ELI.p)
data.p[["TNA"]] <- do.call("rbind",TNA.p)
data.p[["TSA"]] <- do.call("rbind",TSA.p)
data.p <- lapply(data.p, 'colnames<-',c("x","y","i"))
Type <- "point"
data.p <- lapply(data.p, cbind,Type)

### standarization ----
AMM.std <- matrix(NA, ncol=4, nrow = length(Years)); colnames(AMM.std) <- seasons; rownames(AMM.std) <- Years
AEN.std <- matrix(NA, ncol=4, nrow = length(Years)); colnames(AEN.std) <- seasons; rownames(AEN.std) <- Years
ELI.std <- matrix(NA, ncol=4, nrow = length(Years)); colnames(ELI.std) <- seasons; rownames(ELI.std) <- Years
TNA.std <- matrix(NA, ncol=4, nrow = length(Years)); colnames(TNA.std) <- seasons; rownames(TNA.std) <- Years
TSA.std <- matrix(NA, ncol=4, nrow = length(Years)); colnames(TSA.std) <- seasons; rownames(TSA.std) <- Years

sha.test <- list(); indices <- c("AMM","AEN","ELI", "TNA","TSA")
for (i in seasons){
  sha.test[[i]] <- list()
  for (j in indices){
    sha.test[[i]][[j]] <- shapiro.test(get(paste0(j,".p"))[[i]][,1])$p.value
  }
  AMM.std[,i] <- scale(AMM.p[[i]][,1])
  AEN.std[,i] <- scale(AEN.p[[i]][,1])
  ELI.std[,i] <- scale(ELI.p[[i]][,1])
  TNA.std[,i] <- scale(TNA.p[[i]][,1])
  TSA.std[,i] <- scale(TSA.p[[i]][,1])
}
sha.test <- do.call(rbind.data.frame,sha.test)
write.csv(sha.test,file="./Indices_distribution/Indices_shapiro_test.csv")
write.csv(AMM.std, "../../../03_SpatioTemp/03_Composites/AMM_std_1980-2020_ERSST.csv")
write.csv(ELI.std, "../../../03_SpatioTemp/03_Composites/ELI_std_1980-2020_ERSST.csv")
write.csv(AEN.std, "../../../03_SpatioTemp/03_Composites/Atl3_std_1980-2020_ERSST.csv")
write.csv(TNA.std, "../../../03_SpatioTemp/03_Composites/TNA_std_1980-2020_ERSST.csv")
write.csv(TSA.std, "../../../03_SpatioTemp/03_Composites/TSA_std_1980-2020_ERSST.csv")

#### plotting -----

ggplot(data[["AEN"]],aes(x,y))+facet_wrap(~i,ncol=2)+
  geom_line()+labs(title="Atl3 distribution", caption="1980-2020",x="SST anomaly [°C]",y="density")+
  theme_bw()

data2 <- rbind(data[["AMM"]],data.p[["AMM"]])
ggplot(data2,aes(x=x,y=y))+
  geom_line(data=dplyr::filter(data2,Type=="line"))+
  geom_point(data=dplyr::filter(data2,Type=="point"))+
  facet_wrap(~i,ncol=2)+
  labs(title="AMM distribution", caption="1980-2020",x="SST anomaly [°C]",y="density")+
  theme_bw()+
  theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

data2 <- rbind(data[["AEN"]],data.p[["AEN"]])
ggplot(data2,aes(x=x,y=y))+
  geom_line(data=dplyr::filter(data2,Type=="line"))+
  geom_point(data=dplyr::filter(data2,Type=="point"))+
  facet_wrap(~i,ncol=2)+
  labs(title="AEN distribution", caption="1980-2020",x="SST anomaly [°C]",y="density")+
  theme_bw()+
  theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

data2 <- rbind(data[["ELI"]],data.p[["ELI"]])
ggplot(data2,aes(x=x,y=y))+
  geom_line(data=dplyr::filter(data2,Type=="line"))+
  geom_point(data=dplyr::filter(data2,Type=="point"))+
  facet_wrap(~i,ncol=2)+
  labs(title="ELI distribution", caption="1980-2020",x="SST anomaly [°C]",y="density")+
  theme_bw()+
  theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

data2 <- rbind(data[["TNA"]],data.p[["TNA"]])
ggplot(data2,aes(x=x,y=y))+
  geom_line(data=dplyr::filter(data2,Type=="line"))+
  geom_point(data=dplyr::filter(data2,Type=="point"))+
  facet_wrap(~i,ncol=2)+
  labs(title="TNA distribution", caption="1980-2020",x="SST anomaly [°C]",y="density")+
  theme_bw()+
  theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

data2 <- rbind(data[["TSA"]],data.p[["TSA"]])
ggplot(data2,aes(x=x,y=y))+
  geom_line(data=dplyr::filter(data2,Type=="line"))+
  geom_point(data=dplyr::filter(data2,Type=="point"))+
  facet_wrap(~i,ncol=2)+
  labs(title="TSA distribution", caption="1980-2020",x="SST anomaly [°C]",y="density")+
  theme_bw()+
  theme(legend.position=c(0.2,0.9),legend.direction = "horizontal",legend.background = element_rect(linetype="solid", colour ="black"),axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(5,5,5,5),vjust = -1),axis.text.y = element_text(margin=margin(0,7,5,0,"pt")),panel.grid = element_line(linetype="dashed",color="lightgray"))

##### composite identification ----
AMM.bool <- matrix(NA, nrow = length(Years),ncol=4); colnames(AMM.bool) <- seasons; rownames(AMM.bool) <- Years
AMM.bool[which(AMM.std >= 1 ,arr.ind = T)] <- 1
AMM.bool[which(AMM.std <= -1 ,arr.ind = T)] <- -1
AMM.bool[which(AMM.std > -1 & AMM.std < 1,arr.ind = T)] <- 0
AMM.bool2 <- as.integer(t(AMM.bool))
AMM.bool <- data.frame(time=Years,AMM.bool)

AEN.bool <- matrix(NA, nrow = length(Years),ncol=4); colnames(AEN.bool) <- seasons; rownames(AEN.bool) <- Years
AEN.bool[which(AEN.std >= 1 ,arr.ind = T)] <- 1
AEN.bool[which(AEN.std <= -1 ,arr.ind = T)] <- -1
AEN.bool[which(AEN.std > -1 & AEN.std < 1,arr.ind = T)] <- 0
AEN.bool2 <- as.integer(t(AEN.bool))
AEN.bool <- data.frame(time=Years,AEN.bool)

ELI.std <- round(ELI.std,2)
ELI.bool <- matrix(NA, nrow = length(Years),ncol=4); colnames(ELI.bool) <- seasons; rownames(ELI.bool) <- Years
ELI.bool[which(ELI.std >= 0.75 ,arr.ind = T)] <- 1
ELI.bool[which(ELI.std <= -0.74 ,arr.ind = T)] <- -1
ELI.bool[which(ELI.std > -0.74 & ELI.std < 0.75,arr.ind = T)] <- 0
ELI.bool2 <- as.integer(t(ELI.bool))
ELI.bool <- data.frame(time=Years,ELI.bool)

TNA.bool <- matrix(NA, nrow = length(Years),ncol=4); colnames(TNA.bool) <- seasons; rownames(TNA.bool) <- Years
TNA.bool[which(TNA.std >= 1 ,arr.ind = T)] <- 1
TNA.bool[which(TNA.std <= -1 ,arr.ind = T)] <- -1
TNA.bool[which(TNA.std > -1 & TNA.std < 1,arr.ind = T)] <- 0
TNA.bool2 <- as.integer(t(TNA.bool))
TNA.bool <- data.frame(time=Years,TNA.bool)

TSA.bool <- matrix(NA, nrow = length(Years),ncol=4); colnames(TSA.bool) <- seasons; rownames(TSA.bool) <- Years
TSA.bool[which(TSA.std >= 1 ,arr.ind = T)] <- 1
TSA.bool[which(TSA.std <= -1 ,arr.ind = T)] <- -1
TSA.bool[which(TSA.std > -1 & TSA.std < 1,arr.ind = T)] <- 0
TSA.bool2 <- as.integer(t(TSA.bool))
TSA.bool <- data.frame(time=Years,TSA.bool)

data <- data.frame(time= seq(1,length(AMM.bool2)),AMM=AMM.bool2, AEN= AEN.bool2, ELI=ELI.bool2, TNA=TNA.bool2, TSA=TSA.bool2)
data <- melt(data,id="time")
data$variable <- factor(data$variable,levels = c("TSA","TNA","AMM","AEN","ELI"))

ggplot(data,aes(time,variable,fill=value))+
  geom_tile()+scale_fill_gradientn(colours=c("blue","00","red"))+
  scale_x_continuous(breaks=seq(1,length(ELI.bool2),by=4),labels=Years)+
  labs(title="SST seasons",y="Mode")+
  theme_bw()

data2 <- list(AMM=AMM.bool, AEN=AEN.bool, ELI=ELI.bool)
data2 <- melt(data2)
data2$L1 <- factor(data2$L1, levels = c("AMM","AEN","ELI"))


ggplot(data2,aes(time,L1,fill=value))+facet_wrap(.~variable,ncol=1)+
  geom_tile()+
  scale_fill_gradientn(colours=c("blue","00","red"))+
  labs(title="SST seasons",y="Mode",x="year")+
  theme_bw()


data2 <- list(TSA= TSA.bool, TNA=TNA.bool, AMM=AMM.bool)
data2 <- melt(data2)
data2$L1 <- factor(data2$L1, levels = c("TSA","TNA","AMM"))


ggplot(data2,aes(time,L1,fill=value))+facet_wrap(.~variable,ncol=1)+
  geom_tile(aes(alpha=L1))+scale_alpha_manual(values=c(0.5,0.5,1))+
  scale_fill_gradientn(colours=c("blue","00","red"))+
  labs(title="SST seasons",y="Mode",x="year")+
  theme_bw()


## calculate variability ----
{
  AEN.sd <- matrix(NA, ncol=2, nrow = 12)
  colnames(AEN.sd) <- colnames(AEN); rownames(AEN.sd) <- month.abb
  
  AMM.sd <- matrix(NA, ncol=2, nrow = 12)
  colnames(AMM.sd) <- colnames(AMM); rownames(AMM.sd) <- month.abb
  
  ELI.sd <- as.data.frame(matrix(NA, ncol=2, nrow = 12))
  colnames(ELI.sd) <- c("Var1","value"); ELI.sd[,1] <- month.abb
  
  Nino34.sd <- matrix(NA, ncol=2, nrow = 12)
  colnames(Nino34.sd) <- colnames(Nino34); rownames(Nino34.sd) <- month.abb
}# ----------------------------------------------------------------------------------------Initializing Standard deviation for each month

{
  AEN.m <- list(); AEN.eof <- list()
  AEN.m[["H.sst.AEN"]] <- matrix(NA, ncol=12, nrow = length(Years))
  AEN.m[["E.sst.AEN"]] <- matrix(NA, ncol=12, nrow = length(Years))
  
  AMM.m <- list(); AMM.eof <- list()
  AMM.m[["H.sst.AMM"]] <- matrix(NA, ncol=12, nrow = length(Years))
  AMM.m[["E.sst.AMM"]] <- matrix(NA, ncol=12, nrow = length(Years))
  
  ELI.m <-as.data.frame(matrix(NA, ncol=12, nrow = length(Years)))
  
  Nino34.m <-list(); Nino34.eof <- list()
  Nino34.m[["H.sst.N34"]] <- matrix(NA, ncol=12, nrow = length(Years))
  Nino34.m[["E.sst.N34"]] <- matrix(NA, ncol=12, nrow = length(Years))
  
}# ----------------------------------------------------------------------------------------Initializing EOF in months

# calculation
for(i in 1:12){
  aux.H <- AEN[format(Dm.ana,format="%m") == DmM.ana[i],"H.sst.AEN"]
  aux.E <- AEN[format(Dm.ana,format="%m") == DmM.ana[i],"E.sst.AEN"]
  AEN.sd[i,"H.sst.AEN"] <- sd(aux.H, na.rm=T)
  AEN.sd[i,"E.sst.AEN"] <- sd(aux.E, na.rm=T)
  AEN.m[["H.sst.AEN"]][,i] <- aux.H
  AEN.m[["E.sst.AEN"]][,i] <- aux.E
  
  aux.H <- AMM[format(Dm.ana,format="%m") == DmM.ana[i],"H.sst.AMM"]
  aux.E <- AMM[format(Dm.ana,format="%m") == DmM.ana[i],"E.sst.AMM"]
  AMM.sd[i,"H.sst.AMM"] <- sd(aux.H, na.rm=T)
  AMM.sd[i,"E.sst.AMM"] <- sd(aux.E, na.rm=T)
  AMM.m[["H.sst.AMM"]][,i] <- aux.H
  AMM.m[["E.sst.AMM"]][,i] <- aux.E
  
  aux.H <- Nino34[format(Dm.ana,format="%m") == DmM.ana[i],"H.sst.N34"]
  aux.E <- Nino34[format(Dm.ana,format="%m") == DmM.ana[i],"E.sst.N34"]
  Nino34.sd[i,"H.sst.N34"] <- sd(aux.H, na.rm=T)
  Nino34.sd[i,"E.sst.N34"] <- sd(aux.E, na.rm=T)
  Nino34.m[["H.sst.N34"]][,i] <- aux.H
  Nino34.m[["E.sst.N34"]][,i] <- aux.E
  
  aux.E <- ELI[format(Dm.ana,format="%m") == DmM.ana[i]]
  ELI.sd[i,2] <- sd(aux.E[,2], na.rm=T)
  ELI.m[,i] <- aux.E[,2]
}

### plotting standard deviations -----
data <- melt(AEN.sd)
ggplot(data,aes(Var1,value,fill=Var2))+geom_bar(stat="identity",position=position_dodge())+
  scale_fill_brewer(type = "qual",palette = "Set1")+
  labs(title="Standard Deviation for Atl.EN (1980-2020)",x="Month",y="SD [°C]")+
  theme_bw()

data <- melt(AMM.sd)
ggplot(data,aes(Var1,value,fill=Var2))+geom_bar(stat="identity",position=position_dodge())+
  scale_fill_brewer(type = "qual",palette = "Set1")+
  labs(title="Standard Deviation for AMM (1980-2020)",x="Month",y="SD [°C]")+
  theme_bw()

data <- melt(Nino34.sd)
ggplot(data,aes(Var1,value,fill=Var2))+geom_bar(stat="identity",position=position_dodge())+
  scale_fill_brewer(type = "qual",palette = "Set1")+
  labs(title="Standard Deviation for Nino34 (1980-2020)",x="Month",y="SD [°C]")+
  theme_bw()

ELI.sd$Var1 <- factor(ELI.sd$Var1, levels = month.abb)
ggplot(ELI.sd,aes(Var1,value))+geom_bar(stat="identity",position=position_dodge())+
  scale_fill_brewer(type = "qual",palette = "Set1")+
  labs(title="Standard Deviation for ELI (1980-2020)",x="Month",y="SD [°Lon]")+
  theme_bw()

### Calculate time based EOF -----

AEN.m <- lapply(AEN.m, function(x) {apply(x,2, scale, scale=F)})
AMM.m <- lapply(AMM.m, function(x) {apply(x,2, scale, scale=F)})
Nino34.m <- lapply(Nino34.m, function(x) {apply(x,2, scale, scale=F)})
ELI.m <- apply(ELI.m, 2, scale, scale=F)

AEN.eof <- lapply(AEN.m, prcomp, retx = T, scale.= F)
AMM.eof <- lapply(AMM.m, prcomp, retx = T, scale.= F)
Nino34.eof <- lapply(Nino34.m, prcomp, retx = T, scale.= F)
ELI.eof <- prcomp(ELI.m, retx=T, scale. = F)

AEN.var <- lapply(AEN.eof, function(x){var <- x$sdev/sum(x$sdev)*100 })
AMM.var <- lapply(AMM.eof, function(x){var <- x$sdev/sum(x$sdev)*100 })
Nino34.var <- lapply(Nino34.eof, function(x){var <- x$sdev/sum(x$sdev)*100 })
ELI.var <- ELI.eof$sdev/ sum(ELI.eof$sdev)*100

ext.eof <- function(x){
  aux <- as.data.frame(x$rotation[,1:3])
  colnames(aux) <- paste0("PC",seq(1,3))
  aux <- cbind(month.abb,aux)
  return(aux)
}
AEN.eof <- lapply(AEN.eof, ext.eof)
AMM.eof <- lapply(AMM.eof, ext.eof)
Nino34.eof <- lapply(Nino34.eof, ext.eof)
ELI.eof <- as.data.frame(ELI.eof$rotation[,1:3]); colnames(ELI.eof) <- paste0("PC",seq(1,3)); ELI.eof <- cbind(month.abb,ELI.eof)

### plotting time based EOF ----
data <- melt(AEN.eof); data$month.abb <- factor(data$month.abb, levels = month.abb)
ggplot(data, aes(x=month.abb, y= value, fill=L1))+ facet_wrap(.~variable)+
  geom_bar(stat="identity",position=position_dodge())+
  scale_fill_brewer(type = "qual",palette = "Set1", direction=-1)+
  labs(title="Loadings for Atl.EN (1980-2020)",x="Month",y="Loading")+
  theme_bw()

data <- melt(AMM.eof); data$month.abb <- factor(data$month.abb, levels = month.abb)
ggplot(data, aes(x=month.abb, y= value, fill=L1))+ facet_wrap(.~variable)+
  geom_bar(stat="identity",position=position_dodge())+
  scale_fill_brewer(type = "qual",palette = "Set1", direction=-1)+
  labs(title="Loadings for AMM (1980-2020)",x="Month",y="Loading")+
  theme_bw()

data <- melt(Nino34.eof); data$month.abb <- factor(data$month.abb, levels = month.abb)
ggplot(data, aes(x=month.abb, y= value, fill=L1))+ facet_wrap(.~variable)+
  geom_bar(stat="identity",position=position_dodge())+
  scale_fill_brewer(type = "qual",palette = "Set1", direction=-1)+
  labs(title="Loadings for Nino34 (1980-2020)",x="Month",y="Loading")+
  theme_bw()

data <- melt(ELI.eof, id="month.abb"); data$month.abb <- factor(data$month.abb, levels = month.abb)
ggplot(data, aes(x=month.abb, y= value))+ facet_wrap(.~variable)+
  geom_bar(stat="identity",position=position_dodge())+
  scale_fill_brewer(type = "qual",palette = "Set1")+
  labs(title="Loadings for ELI (1980-2020)",x="Month",y="Loading")+
  theme_bw()
