# Author: Melinda de Jonge
# Date: 20-09-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: melinda.dejonge@ru.nl

## -------------------------------------------------------------------------------------------------
#' Script to create figure 3 and S11 based on the final model selected for mammals. 
##--------------------------------------------------------------------------------------------------

library(metafor)
library(ggplot2)
library(pals)
library(viridis)

load(file.path('Data','Processed','GroupedData_main.RData'))
Data_mammals = DataList[['Mammals']]
load(file.path('Output','Models','FinalModelMammals.RData'))

## -----  Create preditions ------------------------------------------------------------------------
Mammals_logD <- seq(0,max(Data_mammals$logD), by = 0.1)
Mammals_Diet = c(0,80)
Mammals_Habitat = as.factor(Data_mammals$Habitat)
Mammals_logBM = seq(min(Data_mammals$logBM),max(Data_mammals$logBM),by=0.1)

Mammals_newgrid <- data.frame(expand.grid(logD = Mammals_logD,
                                          Diet = Mammals_Diet,
                                          logBM = Mammals_logBM,
                                          Habitat = levels(Mammals_Habitat)))

Mammals_newgrid$logD2 = Mammals_newgrid$logD^2
Mammals_predgrid <- model.matrix(~logD*Diet+logD*Habitat+logD*logBM+logD2*logBM,data=Mammals_newgrid)[,-1]
Mammals_preds <- predict(final_model, newmods=Mammals_predgrid, addx = TRUE)

#attach predictions to variables for plotting
Mammals_newgrid$pred <- Mammals_preds$pred
Mammals_newgrid$ci.lb <- Mammals_preds$ci.lb
Mammals_newgrid$ci.ub <- Mammals_preds$ci.ub
Mammals_newgrid$DietName = 'Non-carnivorous'
Mammals_newgrid$DietName[Mammals_newgrid$Diet == 80] = 'Carnivorous'

## ----- Check extrapolation areas -----------------------------------------------------------------
#Add column specifying whether the prediction is extrapolated or not.
checkExtrapolated = function(prediction,Data,BMflex,Dflex){
  maxDiet = ifelse(prediction['DietName']=='Carnivorous',100,20)
  minDiet = ifelse(prediction['DietName']=='Carnivorous',60,0)
  HabType = ifelse(prediction['Habitat']=='Closed','Closed','Open')
  
  SubHabDiet = subset(Data,(Habitat==HabType[1] & Diet >= minDiet[1] & Diet <= maxDiet[1]))
  SubHabDiet = subset(Data,(Habitat==HabType[1] & Diet >= minDiet[1] & Diet <= maxDiet[1]))
  outsideMaxBounds = prediction['logD'] > max(SubHabDiet$logD) | prediction['logD'] < min(SubHabDiet$logD) |
    prediction['logBM'] > max(SubHabDiet$logBM) | prediction['logBM'] < min(SubHabDiet$logBM)
  
  BMupbound = log10(10^as.numeric(prediction['logBM'])/BMflex)
  
  BMlowbound = log10(10^as.numeric(prediction['logBM'])*BMflex)
  
  SubCloser = subset(SubHabDiet,logD>=(as.numeric(prediction['logD'])))
  outsideCloserBMBounds = BMlowbound > max(SubCloser$logBM) | BMupbound < min(SubCloser$logBM)
  
  SubFurther = subset(SubHabDiet,logD<=(as.numeric(prediction['logD'])))
  outsideFurtherBMBounds = BMlowbound > max(SubFurther$logBM) | BMupbound < min(SubFurther$logBM)
  
  outsideDBMBounds = outsideFurtherBMBounds | outsideCloserBMBounds
  
  SubBM = subset(SubHabDiet,logBM>=BMlowbound & logBM<=BMupbound)
  Dupbound = log10(10^as.numeric(prediction['logD'])/Dflex)
  Dlowbound = log10(10^as.numeric(prediction['logD'])*Dflex)
  outsideLimD = (Dlowbound > max(SubBM$logD) | Dupbound < min(SubBM$logD))
  
  extrapolated = outsideMaxBounds | (outsideDBMBounds & outsideLimD)
  return(extrapolated)
}
Mammals_newgrid$Extrapolated = apply(Mammals_newgrid,1,FUN=checkExtrapolated,Data=Data_mammals,BMflex=0.5,Dflex=0.5)
Mammals_newgrid$Extrapolated = ifelse(Mammals_newgrid$Extrapolated,'YES','NO')

## ----- Heeatmap LRR --------------------------------------------------------------------------------
theme.Heatmap = theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border= element_blank(),
                      plot.background = element_blank(),
                      text = element_text(size=8),
                      rect=element_rect(size=rel(0)),
                      legend.key.height=unit(0.75,'cm'),
                      legend.key.width = unit(0.4,'cm'),
                      legend.title=element_text(size=8),
                      legend.text=element_text(size=6),
                      legend.key=element_rect(color='transparent',fill='transparent'),
                      legend.background = element_rect(color='transparent',fill='transparent'),
                      strip.background=element_rect(color='transparent',fill='transparent',size=rel(1)),
                      strip.text = element_text(size=8),
                      panel.spacing=unit(1, "mm"))

plotMHeatmap = ggplot(data=Mammals_newgrid,aes(x=logD,y=logBM))+
  facet_grid(DietName~Habitat)+
  xlab("Distance to infrastructure (m)")+
  ylab("Body mass (kg)") +
  geom_raster(aes(fill=pred))+
  scale_fill_gradientn(colours=rev(ocean.balance(255)),name = expression(LRR^Delta), lim=c(-2,2))+
  geom_point(aes(size=Extrapolated),color='grey25',shape='.')+
  scale_size_manual(values=c(YES=1, NO=NA), guide="none") +
  scale_x_continuous(breaks=c(0,1,2,3),labels=c('1','10','100','1000'),expand=c(0, 0))+
  scale_y_continuous(breaks=c(1,3,5),labels=c('0.01','1','100'),expand=c(0, 0))+
  theme.Heatmap

cairo_pdf(filename = 'Figures/Fig3.pdf', height=2.9, width=3.5)
plot(plotMHeatmap)
dev.off()

## ---- Effect zones -------------------------------------------------------------------------------
parsMams = final_model['b'][[1]]
Mammals_Diet = seq(0,100,by=10)
Mammals_newgrid = data.frame(expand.grid(Diet=Mammals_Diet,
                                         logBM = Mammals_logBM,
                                         Habitat = levels(Mammals_habitat)))

Mammals_coefsum = data.frame(c = parsMams[1,] + Mammals_newgrid$Diet*parsMams[3,] + Mammals_newgrid$logBM*parsMams[5,],
                             b = parsMams[2,] + Mammals_newgrid$Diet*parsMams[7,] + Mammals_newgrid$logBM*parsMams[9,],
                             a = parsMams[6,] + Mammals_newgrid$logBM*parsMams[10,])

Mammals_coefsum$c[Mammals_newgrid$Habitat=="Open"] = Mammals_coefsum$c[Mammals_newgrid$Habitat=="Open"] + parsMams[4,]
Mammals_coefsum$b[Mammals_newgrid$Habitat=="Open"] = Mammals_coefsum$b[Mammals_newgrid$Habitat=="Open"] + parsMams[8,]

abc <- function(X){
  D=X$b^2-4*X$c*X$a
  x1 =(-X$b+sqrt(D))/(2*X$a)
  x2 =(-X$b-sqrt(D))/(2*X$a)
  ans = data.frame(first=x1,
                   second=x2)
  return(ans)
}

intersections = abc(Mammals_coefsum)
Mammals_newgrid$a = Mammals_coefsum$a
Mammals_newgrid$b = Mammals_coefsum$b
Mammals_newgrid$c = Mammals_coefsum$c
Mammals_newgrid$int1 = intersections[,1]
Mammals_newgrid$int2 = intersections[,2]
Mammals_newgrid$logRR1m = Mammals_newgrid$a*log10(2)^2+Mammals_newgrid$b*log10(2) + Mammals_newgrid$c
Mammals_newgrid$locx = -Mammals_newgrid$b/(2*Mammals_newgrid$a)
Mammals_newgrid$locy = Mammals_newgrid$a*Mammals_newgrid$locx^2+Mammals_newgrid$b*Mammals_newgrid$locx + Mammals_newgrid$c
Mammals_newgrid$locy2 = Mammals_newgrid$c - Mammals_newgrid$b^2/(4*Mammals_newgrid$a)


Mammals_newgrid$IEZ = NA
Mammals_newgrid$resp = NA
Mammals_newgrid$IEZ[Mammals_newgrid$a>0 & Mammals_newgrid$locx>0 & Mammals_newgrid$logRR1m>0] = Mammals_newgrid$int2[Mammals_newgrid$a>0 & Mammals_newgrid$locx>0 & Mammals_newgrid$logRR1m>0]
Mammals_newgrid$resp[Mammals_newgrid$a>0 & Mammals_newgrid$locx>0 & Mammals_newgrid$logRR1m>0] = "Increase"
Mammals_newgrid$IEZ[Mammals_newgrid$a>0 & Mammals_newgrid$locx>0 & Mammals_newgrid$logRR1m<0] = Mammals_newgrid$int1[Mammals_newgrid$a>0 & Mammals_newgrid$locx>0 & Mammals_newgrid$logRR1m<0]
Mammals_newgrid$resp[Mammals_newgrid$a>0 & Mammals_newgrid$locx>0 & Mammals_newgrid$logRR1m<0] = "Decrease"
Mammals_newgrid$IEZ[Mammals_newgrid$a>0 & Mammals_newgrid$locx<0] = Mammals_newgrid$int1[Mammals_newgrid$a>0 & Mammals_newgrid$locx<0]
Mammals_newgrid$resp[Mammals_newgrid$a>0 & Mammals_newgrid$locx<0] = "Decrease"
Mammals_newgrid$IEZ[Mammals_newgrid$a<0 & Mammals_newgrid$locx>0 & Mammals_newgrid$int1>0] = Mammals_newgrid$int1[Mammals_newgrid$a<0 & Mammals_newgrid$locx>0& Mammals_newgrid$int1>0]
Mammals_newgrid$IEZ[Mammals_newgrid$a<0 & Mammals_newgrid$locx>0 & Mammals_newgrid$int1<0] = Mammals_newgrid$int2[Mammals_newgrid$a<0 & Mammals_newgrid$locx>0& Mammals_newgrid$int1<0]
Mammals_newgrid$resp[Mammals_newgrid$a<0 & Mammals_newgrid$locx>0 & Mammals_newgrid$int1<0] = "Increase"
Mammals_newgrid$resp[Mammals_newgrid$a<0 & Mammals_newgrid$locx>0 & Mammals_newgrid$int1>0] = "Decrease"
Mammals_newgrid$IEZ[Mammals_newgrid$a<0 & Mammals_newgrid$locx<0] = Mammals_newgrid$int2[Mammals_newgrid$a<0 & Mammals_newgrid$locx<0]
Mammals_newgrid$resp[Mammals_newgrid$a<0 & Mammals_newgrid$locx<0] = "Increase"

## -------- IEZ plots -----------------------------------------------------------
DataZones = Mammals_newgrid
DataZones$IEZ[DataZones$resp=='Increase'] = NA

IEZHeatmap = ggplot(DataZones,aes(x=Diet,y=logBM))+
  facet_grid(~Habitat)+
  ggtitle('Mammals: Infrastructure effect zones')+
  geom_raster(aes(fill=IEZ))+
  scale_fill_viridis_c(na.value = 'transparent',direction = -1,
                       name = 'IEZ (m)', lim=c(0,3),
                       breaks = c(3,2,1,0),labels = c('1000','100','10','1'))+
  xlab('% vertebrates in diet') +
  scale_x_continuous(breaks = c(0,20,40,60,80,100)) +
  scale_y_continuous(breaks = 1:6,labels=c('0.01','0.1','1','10','100','1000')) +
  ylab('Mass (kg)') +
  theme_bw() +
  theme.Heatmap

cairo_pdf(filename = 'Figures/FigS11.pdf', height=2.3, width=4.7)
plot(IEZHeatmap)
dev.off()

