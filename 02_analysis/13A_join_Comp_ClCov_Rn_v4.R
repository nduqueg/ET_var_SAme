rm(list=ls())
cat("\014")

library(hydroTSM)
library(terra)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
"%>%"=magrittr::'%>%'

# loading precalculated composites ----
load("../03_Composites/07_Rn/Radiation_Composites.RData"); Comp.Rn <- Comp.E5; rm(Comp.E5)
load("../03_Composites/07_Rn/20_CldCLARA-day_Composites.RData"); Comp.Cl <- Comp.E5; rm(Comp.E5)

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

data.Cl <- melt(Comp.Cl, id=c("lon","lat")) %>% 
  within(., {
    Cl.class <- value > 0
    L2 <- factor(L2, levels = c("Pos","Neg"))
    })
colnames(data.Cl) <- c("x","y", "Season", "value", "Phase", "indice","Cl.class")


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

load("../03_Composites/07_Rn/Radiation_Composites_Ttest_shape.RData")
# ERA5L
plot.sig <- function(data.Rns, data.Cls,  shape.spec, rect.spec, title.p=NA){
  
  at.m <- seq(-15,15,length.out=11)
  at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]
  Phase.labs <- c("Positive", "Negative"); names(Phase.labs) <- c("Pos","Neg")
  p <- ggplot()+ facet_wrap(. ~ Phase,  labeller= labeller(Phase = Phase.labs))+
    geom_raster(data= data.Rns , aes(x,y,fill=value)) +scale_fill_stepsn(colours=brewer.pal(10,"PuOr") %>% rev(), breaks=at.m,
                                                                               values=scales::rescale(at.m.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                                                                               guide=guide_colorsteps(even.steps = T,barheight=unit(6,"cm")),
                                                                               name="Rn Anom.\n[W/m2]")+
    
    geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.4)+
    geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed",colour="black",fill="NA",linewidth=0.05)+
    
    geom_sf_pattern(data = shape.spec, pattern= "circle", pattern_density=0.1, pattern_spacing= 0.02, fill="00", colour="00", pattern_colour="black")+
    
    geom_contour(data= subset(data.Cls, Cl.class), aes(x,y, z=value),col="blue", breaks=c(4,8,12), linewidth=0.4)+
    geom_contour(data= subset(data.Cls, !Cl.class), aes(x,y, z=value),col="red", breaks=c(-12,-8,-4), linewidth=0.4)+
    
    geom_rect(data= rect.spec,inherit.aes = F,aes(xmin=xmi, xmax=xma, ymin=ymi, ymax=yma),color="black",linetype="longdash", fill=NA, linewidth=1)+
    
    scale_y_continuous(position="right")+
    
    coord_sf(xlim=c(-80,-35),ylim=c(-20,13))+
    
    theme_bw()+theme(strip.background = element_blank(),strip.text.x = element_blank(),
                     legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=15),
                     axis.text.x= element_blank(), axis.title.x = element_blank(),
                     legend.position = "right")
  
  if (is.na(title.p)) p <- p+ labs(x="Longitude [°]",y="Latitude [°]") else p <- p+ labs(x="Longitude [°]",y="Latitude [°]", title = title.p)
  return(p)
}

for ( i in c("MAM","JJA","SON")){
  p <- plot.sig(data.Rns = subset(data.Rn, indice=="AMM" & Season == i),
                data.Cls = subset(data.Cl, indice=="AMM" & Season == i),
                shape.spec = subset(shape.sig, Mode=="AMM" & Season ==i),
                rect.spec = subset(rect, mode=="AMM" & Season ==i))
  print(p)
}

plot.sig(data.Rns = subset(data.Rn, indice=="Atl3" & Season == "JJA"),
         data.Cls = subset(data.Cl, indice=="Atl3" & Season == "JJA"),
         shape.spec = subset(shape.sig, Mode=="Atl3" & Season =="JJA"),
         rect.spec = subset(rect, mode=="Atl3" & Season =="JJA")) # 1200 x 380


