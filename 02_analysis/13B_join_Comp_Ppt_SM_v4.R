rm(list=ls())
cat("\014")

library(hydroTSM)
library(raster)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
"%>%"=magrittr::'%>%'

# loading precalculated composites ----
Comp.Ppt <- list()
load("../03_Composites/03_Ppt/03_AMM_Comp_Ppt.RData"); Comp.Ppt[["AMM"]] <- data.anom2; rm(data.anom2)
load("../03_Composites/03_Ppt/03_Atl3_Comp_Ppt.RData"); Comp.Ppt[["Atl3"]] <- data.anom2; rm(data.anom2)
load("../03_Composites/06_SM/SatPer_anom_Composites.RData"); Comp.SM <- Comp.E5.anom; rm(Comp.E5.anom, Comp.CCI.anom)

Basins <- shapefile("../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")
SA <- shapefile("../../01_DataSets/South_America/South_America.shp")

## datasets unification ----
data.Ppt <- melt(Comp.Ppt, id=c("lon","lat","Season","Anomaly","Dir")) %>% 
  within(., {
    Ppt.class <- Anomaly > 0
    })
colnames(data.Ppt) <- c("x","y", "Season", "value", "Phase", "indice", "Ppt.class")


data.SM <- melt(Comp.SM, id=c("x","y")) %>% 
  within(., {
    SM.class <- value > 0
    L2 <- factor(L2, levels = c("Pos","Neg"))
    })
colnames(data.SM) <- c("x","y", "Season", "value", "Phase", "indice","SM.class")


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

load("../03_Composites/06_SM/SatPer_Composites_Ttest_shape.RData")
# ERA5L
plot.sig <- function(data.SMs, data.Ppts, shape.spec, rect.spec, title.p=NA){
  
  at.m <- seq(-10,10,length.out=11)
  at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]
  col2alpha <- function(someColor, alpha=100){ newColor <- col2rgb(someColor); apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}
  paleta <- brewer.pal(10,"RdYlBu") 
  paleta <- col2alpha(paleta, alpha = 0.7*255)
  
  Phase.labs <- c("Positive", "Negative"); names(Phase.labs) <- c("Pos","Neg")
  p <- ggplot()+ facet_wrap(. ~ Phase,  labeller= labeller(Phase = Phase.labs))+
    geom_raster(data= data.SMs , aes(x,y,fill=value)) +scale_fill_stepsn(colours= paleta, breaks=at.m,
                                                                               values=scales::rescale(at.m.v,from=range(at.m)),limits=c(min(at.m),max(at.m)),
                                                                               guide=guide_colorsteps(even.steps = T,barheight=unit(6,"cm")),
                                                                               name="SM Anom.\n[Sat %]")+
    
    geom_polygon(data=Basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",linewidth=0.4)+
    geom_polygon(data=SA,aes(x=long,y=lat, group=group), linetype="dashed",colour="black",fill="NA",linewidth=0.05)+
    
    geom_sf_pattern(data = shape.spec, pattern= "circle", pattern_density=0.1, pattern_spacing= 0.02, fill="00", colour="00", pattern_colour="black")+
    
    geom_contour(data= subset(data.Ppts, Ppt.class), aes(x,y, z=value),col="#01665e", breaks=c(100,200,300,400), linewidth=0.4)+
    geom_contour(data= subset(data.Ppts, !Ppt.class), aes(x,y, z=value),col="#8c510a", breaks=c(-400,-300,-200,-100), linewidth=0.4)+
    
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
  p <- plot.sig(data.SMs = subset(data.SM, indice=="AMM" & Season == i),
                data.Ppts = subset(data.Ppt, indice=="AMM" & Season == i),
                shape.spec = subset(shape.sig, dataset=="ERA5L" &  Mode=="AMM" & Season ==i),
                rect.spec = subset(rect, mode=="AMM" & Season ==i))
  print(p)
}

plot.sig(data.SMs = subset(data.SM, indice=="Atl3" & Season == "JJA"),
         data.Ppts = subset(data.Ppt, indice=="Atl3" & Season == "JJA"),
         shape.spec = subset(shape.sig, dataset=="ERA5L" &  Mode=="Atl3" & Season =="JJA"),
         rect.spec = subset(rect, mode=="Atl3" & Season =="JJA")) # 1200 x 380


