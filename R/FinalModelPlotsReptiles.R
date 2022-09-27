# Author: Melinda de Jonge
# Date: 20-09-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: melinda.dejonge@ru.nl

## -------------------------------------------------------------------------------------------------
#' Script to create figure 5 based on the final model selected for reptiles 
##--------------------------------------------------------------------------------------------------
library(metafor)
library(ggplot2)

load(file.path('Data','Processed','GroupedData_main.RData'))
Data_reptiles = DataList[['Reptiles']]
load(file.path('Output','Models','FinalModelReptiles.RData'))

## -----  Create preditions ------------------------------------------------------------------------
Reptiles_logD <- seq(0,max(Data_reptiles$logD), by = 0.1)
Reptiles_Habitat = as.factor(Data_reptiles$Habitat)
Reptiles.newgrid <- data.frame(expand.grid(logD = Reptiles_logD,
                                           Habitat = levels(Reptiles_Habitat)))
Reptiles.newgrid$logD2 = Reptiles.newgrid$logD^2
Reptiles.predgrid <- model.matrix(~logD*Habitat+logD2,data=Reptiles.newgrid)[,-1]
Reptiles.preds <- predict(final_model, newmods=Reptiles.predgrid, addx = TRUE)

#attach predictions to variables for plotting
Reptiles.newgrid$pred <- Reptiles.preds$pred 
Reptiles.newgrid$ci.lb <- Reptiles.preds$ci.lb 
Reptiles.newgrid$ci.ub <- Reptiles.preds$ci.ub 

## -----  Fig. 5 -----------------------------------------------------------------------------------
plotTheme = theme(panel.grid.major = element_blank(),
                  panel.background = element_blank(),
                  plot.background = element_blank(),
                  text = element_text(size=8), 
                  title = element_text(size=8),
                  legend.text=element_text(size=8),
                  legend.background = element_blank(),
                  legend.direction = 'horizontal',
                  legend.position = 'bottom',
                  legend.title = element_blank())

Reptiles.plot <- ggplot(data=Reptiles.newgrid, aes(x=logD, y=pred, group=Habitat)) + 
  geom_hline(yintercept=0, linetype="longdash", color="grey", size=1) +
  geom_ribbon(aes(ymin=ci.lb, ymax=ci.ub),linetype=0, alpha=0.1) +
  scale_size(guide="none") +
  geom_line(size=0.7, aes(linetype=Habitat)) +
  scale_linetype_manual(breaks=c('Open','Closed'), labels=c('Open habitat','Closed habitat'),
                        values=c("dotdash","solid"), drop=FALSE, name='') +
  theme_bw() +
  ylab(expression(LRR^Delta)) +
  ylim(c(-1.5,1.5)) +
  scale_x_continuous(name ="Distance to infrastructure (m)", 
                     breaks= c(0,1,2,3),
                     labels=c("0" = "1","1" = "10","2" = "100","3" = "1000")) +
  plotTheme

cairo_pdf(file='Figures/Fig5.pdf',
          width = 2.4, height = 2)
plot(Reptiles.plot)
dev.off()