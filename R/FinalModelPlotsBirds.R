# Author: Melinda de Jonge
# Date: 20-09-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: melinda.dejonge@ru.nl

## -------------------------------------------------------------------------------------------------
#' Script to create figure 4 and S10 based on the final model selected for birds. 
##--------------------------------------------------------------------------------------------------

library(metafor)
library(ggplot2)
library(pals)
library(viridis)

load(file.path('Data','Processed','GroupedData_main.RData'))
Data_birds = DataList[['Birds']]
load(file.path('Output','Models','FinalModelBirds.RData'))

## -----  Create preditions ------------------------------------------------------------------------
Birds_logD <- seq(0,max(Data_birds$logD), by = 0.1)
Birds_Habitat = as.factor(Data_birds$Habitat)
Birds_InfraType = as.factor(Data_birds$InfraType)
Birds_Diet <- seq(0,100, by = 10)

Birds_newgrid <- data.frame(expand.grid(logD = Birds_logD, 
                                        Diet = Birds_Diet,
                                        Habitat = levels(Birds_Habitat),
                                        InfraType = levels(Birds_InfraType)))

Birds_newgrid$logD2 = Birds_newgrid$logD^2
Birds_predgrid <- model.matrix(~logD*Diet+logD*Habitat+InfraType,data=Birds_newgrid)[,-1]
Birds_preds <- predict(final_model, newmods=Birds_predgrid, addx = TRUE)

Birds_newgrid$pred <- Birds_preds$pred 
Birds_newgrid$ci.lb <- Birds_preds$ci.lb 
Birds_newgrid$ci.ub <- Birds_preds$ci.ub 

## ----- Plot for paved roads (Fig. 4) -------------------------------------------------------------
plotTheme = theme(panel.grid.major = element_blank(),
                  panel.background = element_blank(),
                  plot.background = element_blank(),
                  text = element_text(size=8), 
                  title = element_text(size=8),
                  rect=element_rect(size=rel(0.8)),
                  legend.title=element_text(size=8),
                  legend.text=element_text(size=6),
                  strip.background=element_rect(color='transparent',fill='transparent',size=rel(1)),
                  strip.text = element_text(size=8),
                  panel.spacing=unit(2, "mm"))

Birds_newgridPaved = subset(Birds_newgrid,InfraType == "Paved road")
subsCIPaved = subset(Birds_newgridPaved,Diet==0 | Diet==100) # Confidence intervalls only for extremes

Birds.plotPaved<- ggplot(data=Birds_newgridPaved, aes(x=logD, y=pred, colour=Diet, group=as.factor(Diet))) + 
  facet_grid(~Habitat)+
  geom_hline(yintercept = 0, linetype="longdash", color = "grey", size = 1)+
  geom_line(size = 0.7)+
  geom_line(data=subsCIPaved,aes(x=logD, y=ci.ub, colour=Diet, group=as.factor(Diet)),linetype = "dotdash",size=0.5)+
  geom_line(data=subsCIPaved,aes(x=logD, y=ci.lb, colour=Diet, group=as.factor(Diet)),linetype = "dotdash",size=0.5)+
  theme_bw()+
  scale_colour_viridis_c(option="inferno",direction=-1,end = 0.85)+
  ylab(expression(LRR^Delta))+
  ylim(c(-3,3))+
  scale_x_continuous(breaks=c(0,1,2,3,4),labels=c('1','10','100','1000','10000'),expand=c(0,0))+
  xlab("Distance to infrastructure (m)")+
  plotTheme + 
  theme(legend.key.height=unit(0.55,'cm'),
        legend.key.width = unit(0.2,'cm'))

cairo_pdf(filename = 'Figures/Fig4.pdf', height=2, width=4.7)
plot(Birds.plotPaved)
dev.off()

## ----- Plot all infrastructure types (Fig. S10) --------------------------------------------------
subsCI = subset(Birds_newgrid,Diet==0 | Diet==100)
Birds.plotAll <- ggplot(data=Birds_newgrid, aes(x=logD, y=pred, colour=Diet, group=as.factor(Diet))) + 
  facet_grid(InfraType~Habitat)+
  geom_hline(yintercept = 0, linetype="longdash", color = "grey", size = 1)+
  geom_line(size = 0.7)+
  geom_line(data=subsCI,aes(x=logD, y=ci.ub, colour=Diet, group=as.factor(Diet)),linetype = "dotdash",size=0.5)+
  geom_line(data=subsCI,aes(x=logD, y=ci.lb, colour=Diet, group=as.factor(Diet)),linetype = "dotdash",size=0.5)+
  theme_bw()+
  scale_colour_viridis_c(option="inferno",direction=-1,end = 0.85)+
  ylab(expression(LRR^Delta))+
  ylim(c(-3,3))+
  scale_x_continuous(breaks=c(0,1,2,3,4),labels=c('1','10','100','1000','10000'),expand=c(0,0))+
  xlab("Distance to infrastructure (m)") +
  plotTheme + 
  theme(legend.key.height=unit(0.75,'cm'),
        legend.key.width = unit(0.4,'cm'))

cairo_pdf(filename = 'Figures/FigS10.pdf', height=5.6, width=4.7)
plot(Birds.plotAll)
dev.off()