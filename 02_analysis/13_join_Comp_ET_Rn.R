rm(list=ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# loading precalculated composites ----
load("./01_data/05_Rad/Radiation_Composites.RData"); Comp.Rn <- Comp.E5
load("./01_data/04_Evap/ET_Composites.RData")

Basins <- shapefile("./01_data/hybas_sa_lev03_v1c.shp")
SA <- shapefile("./01_data/South_America.shp")

## datasets unification ----
data.Rn <- melt(Comp.Rn, id=c("lon","lat")); data.Rn <- data.Rn[!is.na(data.Rn$value),]
Rn.class <- data.Rn$value > 0; data.Rn <- cbind(data.Rn,Rn.class)
colnames(data.Rn) <- c("x","y", "Season", "value", "Phase", "indice", "Rn.class")
data.Rn$Phase <- factor(data.Rn$Phase, levels = c("Pos","Neg"))


data.e5 <- melt(Comp.E5, id=c("x","y"))
colnames(data.e5) <- c("x","y", "Season", "value", "Phase", "indice")
data.e5$Phase <- factor(data.e5$Phase, levels = c("Pos","Neg"))
data.gl <- melt(Comp.gl, id=c("x","y"))
colnames(data.gl) <- c("x","y", "Season", "value", "Phase", "indice")
data.gl$Phase <- factor(data.gl$Phase, levels = c("Pos","Neg"))

#### plotting vertical ----
library(gridExtra)
library(ggpattern)

seasons <- c("DJF","MAM","JJA","SON")

rect <- data.frame(xmi=c(-70,-70,-73,-65),
                   xma=c(-58,-60,-67,-57),
                   ymi=c(3,-5,2,0),
                   yma=c(10,3,7,6),
                   Season= factor(c("MAM","JJA","SON","JJA"), levels = c("MAM","JJA","SON")),
                   mode=c("AMM","AMM","AMM","Atl3"),
                   Phase = factor("Pos"))

load("./01_data/04_Evap/ET_Composites_Ttest_shape.RData")
# ERA5L
plot.sig <- function(data.ET, data.Rns, shape.spec, rect.spec, title.p=NA){
  
  at.m <- c(-80,-50,-25,-15,-10,-5,5,10,15,25,50,80)
  at.v <- c(-65,-37.5,-20,-12.5,-7.5,0,7.5,12.5,20,37.5,65)
  Phase.labs <- c("Positive", "Negative"); names(Phase.labs) <- c("Pos","Neg")
  p <- ggplot()+ facet_wrap(. ~ Phase,  labeller= labeller(Phase = Phase.labs))+
    geom_raster(data= data.ET , aes(x,y,fill=value)) +scale_fill_stepsn(colours=brewer.pal(11,"PiYG"), breaks=at.m,
                                                                               values=scales::rescale(at.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                                                                               guide=guide_colorsteps(even.steps = T,barwidth=unit(12,"cm")),
                                                                               name="Anom.\n[mm]")+
    
    geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.4)+
    geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed",colour="black",fill="NA",linewidth=0.05)+
    
    geom_sf_pattern(data = shape.spec, pattern= "circle", pattern_density=0.1, pattern_spacing= 0.02, fill="00", colour="00", pattern_colour="black")+
    
    geom_contour(data= dplyr::filter(data.Rns, Rn.class==T), aes(x,y,z=value),col="red", breaks=c(3,6,9,12), linewidth=0.6)+
    geom_contour(data= dplyr::filter(data.Rns, Rn.class==F), aes(x,y,z=value),col="blue", breaks=c(-12,-9,-6,-3), linewidth=0.6)+
    
    geom_rect(data= rect.spec,inherit.aes = F,aes(xmin=xmi, xmax=xma, ymin=ymi, ymax=yma),color="black",linetype="longdash", fill=NA, linewidth=1)+
    
    scale_y_continuous(position="right")+
    
    coord_sf(xlim=c(-80,-35),ylim=c(-20,13))+
    
    theme_bw()+theme(strip.background = element_blank(),strip.text.x = element_blank(),
                     legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=15),
                     legend.position = "bottom")
  
  if (is.na(title.p)) p <- p+ labs(x="Longitude [°]",y="Latitude [°]") else p <- p+ labs(x="Longitude [°]",y="Latitude [°]", title = title.p)
  return(p)
}

for ( i in c("MAM","JJA","SON")){
  p <- plot.sig(data.ET = subset(data.e5, indice=="AMM" & Season == i),
                data.Rns = subset(data.Rn, indice=="AMM" & Season == i),
                shape.spec = subset(shape.sig, dataset=="ERA5L" &  Mode=="AMM" & Season ==i),
                rect.spec = subset(rect, mode=="AMM" & Season ==i))
  print(p)
}

plot.sig(data.ET = subset(data.e5, indice=="Atl3" & Season == "JJA"),
         data.Rns = subset(data.Rn, indice=="Atl3" & Season == "JJA"),
         shape.spec = subset(shape.sig, dataset=="ERA5L" &  Mode=="Atl3" & Season =="JJA"),
         rect.spec = subset(rect, mode=="Atl3" & Season =="JJA"))



# GLEAM
{
  data.Gleam.AMM <- data.gl[data.gl$indice=="AMM" & data.gl$Season!="DJF",]
  data.Rn.AMM <- data.Rn[data.Rn$indice=="AMM" & data.Rn$Season!="DJF",]
  
  p1 <- ggplot(data.Gleam.AMM , aes(x, y))+ facet_grid(Season ~ Phase, switch = "y")+
    geom_raster(aes(fill=value)) +scale_fill_stepsn(colours=brewer.pal(11,"PiYG"), breaks=at.m,
                                                    values=scales::rescale(at.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                                                    # guide=guide_colorsteps(even.steps = F),
                                                    name="ET anom. [mm]")+
    
    geom_contour(data= dplyr::filter(data.Rn.AMM, Rn.class==T), aes(z=value),col="red", breaks=c(3,6,9,12), linewidth=0.6)+
    geom_contour(data= dplyr::filter(data.Rn.AMM, Rn.class==F), aes(z=value),col="blue", breaks=c(-12,-9,-6,-3), linewidth=0.6)+
    
    scale_y_continuous(position="right")+
    geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
    geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed",colour="black",fill="NA",linewidth=0.05)+
    coord_fixed(xlim=c(-82,-35),ylim=c(-20,15))+labs(x="Long",y="Lat", title="A) AMM")+ # ,title="Correlation Coef. Indices vs Rn-ERA5L & ET-GLEAM (1980-2020)"
    theme_bw()+theme(legend.position = "none", axis.text.x= element_blank(), axis.title.x = element_blank())
  
  data.Gleam.Atl3 <- data.gl[data.gl$indice=="Atl3" & data.gl$Season=="JJA",]
  data.Rn.Atl3 <- data.Rn[data.Rn$indice=="Atl3" & data.Rn$Season=="JJA",]
  
  p2 <- ggplot(data.Gleam.Atl3 , aes(x, y))+ facet_grid(Season ~Phase, switch = "y")+
    geom_raster(aes(fill=value)) +scale_fill_stepsn(colours=brewer.pal(11,"PiYG"), breaks=at.m,
                                                    values=scales::rescale(at.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                                                    guide=guide_colorsteps(barwidth=unit(10,"cm")),
                                                    name="ET anom. [mm]")+
    
    geom_contour(data= dplyr::filter(data.Rn.Atl3, Rn.class==T), aes(z=value),col="red", breaks=c(3,6,9,12), linewidth=0.6)+
    geom_contour(data= dplyr::filter(data.Rn.Atl3, Rn.class==F), aes(z=value),col="blue", breaks=c(-12,-9,-6,-3), linewidth=0.6)+
    
    scale_y_continuous(position="right")+
    geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.15)+
    geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed",colour="black",fill="NA",linewidth=0.05)+
    coord_fixed(xlim=c(-82,-35),ylim=c(-20,15))+labs(x="Long",y="Lat", title="B) Atl3")+ # , caption = "GLEAM"
    theme_bw()+theme(legend.position = "bottom")
  
  grid.arrange(p1,p2,nrow=2, heights=c(2.05,1)) # 900x1500
}
