rm(list=ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
"%>%"=magrittr::'%>%'

# loading precalculated composites ----
load("../03_Composites/07_Rn/Radiation_Composites.RData"); Comp.Rn <- Comp.E5; rm(Comp.E5)
load("../03_Composites/06_SM/SatPer_anom_Composites.RData"); Comp.SM <- Comp.E5.anom; rm(Comp.E5.anom, Comp.CCI.anom)
load("../03_Composites/08_ET/ET_Composites.RData")

Basins <- shapefile("../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")
SA <- shapefile("../../01_DataSets/South_America/South_America.shp")

## datasets unification ----
data.Rn <- melt(Comp.Rn, id=c("lon","lat")) %>% 
  subset(., !is.na(value)) %>% 
  within(., {
    Rn.class <- value > 0
    L2 <- factor(L2, levels = c("Pos","Neg"))
    })
colnames(data.Rn) <- c("x","y", "Season", "value", "Phase", "indice", "Rn.class")


data.e5 <- melt(Comp.E5, id=c("x","y")) %>% 
  within(., L2 <- factor(L2, levels = c("Pos","Neg")))
colnames(data.e5) <- c("x","y", "Season", "value", "Phase", "indice")

data.SM <- melt(Comp.SM, id=c("x","y")) %>% 
  within(., {
    SM.class <- value > 0
    L2 <- factor(L2, levels = c("Pos","Neg"))
    })
colnames(data.SM) <- c("x","y", "Season", "value", "Phase", "indice","SM.class")


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

load("../03_Composites/08_ET/ET_Composites_Ttest_shape.RData")
# ERA5L
plot.sig <- function(data.ET, data.Rns, data.SMs, shape.spec, rect.spec, title.p=NA){
  
  at.m <- c(-80,-50,-25,-15,-10,-5,5,10,15,25,50,80)
  at.v <- c(-65,-37.5,-20,-12.5,-7.5,0,7.5,12.5,20,37.5,65)
  Phase.labs <- c("Positive", "Negative"); names(Phase.labs) <- c("Pos","Neg")
  p <- ggplot()+ facet_wrap(. ~ Phase,  labeller= labeller(Phase = Phase.labs))+
    geom_raster(data= data.ET , aes(x,y,fill=value)) +scale_fill_stepsn(colours=brewer.pal(11,"PiYG"), breaks=at.m,
                                                                               values=scales::rescale(at.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                                                                               guide=guide_colorsteps(even.steps = T,barheight=unit(6,"cm")),
                                                                               name="ET Anom.\n[mm]")+
    
    geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.4)+
    geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed",colour="black",fill="NA",linewidth=0.05)+
    
    geom_sf_pattern(data = shape.spec, pattern= "circle", pattern_density=0.1, pattern_spacing= 0.02, fill="00", colour="00", pattern_colour="black")+
    
    geom_contour(data= subset(data.SMs, SM.class), aes(x,y, z=value),col="blue", breaks=c(5,10,15,20), linewidth=0.3)+
    geom_contour(data= subset(data.SMs, !SM.class), aes(x,y, z=value),col="red", breaks=c(-20,-15,-10,-5), linewidth=0.3)+
    
    geom_contour(data= subset(data.Rns, Rn.class), aes(x,y, z=value),col="#b35806", breaks=c(3,6,9,12), linewidth=0.3)+
    geom_contour(data= subset(data.Rns, !Rn.class), aes(x,y, z=value),col="#01665e", breaks=c(-12,-9,-6,-3), linewidth=0.3)+
    
    geom_rect(data= rect.spec,inherit.aes = F,aes(xmin=xmi, xmax=xma, ymin=ymi, ymax=yma),color="black",linetype="longdash", fill=NA, linewidth=1)+
    
    scale_y_continuous(position="right")+
    
    coord_sf(xlim=c(-80,-35),ylim=c(-20,13))+
    
    theme_bw()+theme(strip.background = element_blank(),strip.text.x = element_blank(),
                     legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=15),
                     legend.position = "right")
  
  if (is.na(title.p)) p <- p+ labs(x="Longitude [°]",y="Latitude [°]") else p <- p+ labs(x="Longitude [°]",y="Latitude [°]", title = title.p)
  return(p)
}

for ( i in c("MAM","JJA","SON")){
  p <- plot.sig(data.ET = subset(data.e5, indice=="AMM" & Season == i),
                data.Rns = subset(data.Rn, indice=="AMM" & Season == i),
                data.SMs = subset(data.SM, indice=="AMM" & Season == i),
                shape.spec = subset(shape.sig, dataset=="ERA5L" &  Mode=="AMM" & Season ==i),
                rect.spec = subset(rect, mode=="AMM" & Season ==i))
  print(p)
}

plot.sig(data.ET = subset(data.e5, indice=="Atl3" & Season == "JJA"),
         data.Rns = subset(data.Rn, indice=="Atl3" & Season == "JJA"),
         data.SMs = subset(data.SM, indice=="Atl3" & Season == "JJA"),
         shape.spec = subset(shape.sig, dataset=="ERA5L" &  Mode=="Atl3" & Season =="JJA"),
         rect.spec = subset(rect, mode=="Atl3" & Season =="JJA")) # 1200 x 420


