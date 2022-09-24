# Author: Melinda de Jonge
# Date: 20-09-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: melinda.dejonge@ru.nl

## -------------------------------------------------------------------------------------------------
#' 
##--------------------------------------------------------------------------------------------------
library(metafor)
library(ggplot2)
source(file.path('R','Functions.R'))

load(file.path('Data','Processed','GroupedData_main.RData'))

models_main = list('Mammals'=NA, 'Birds'=NA, 'Reptiles'=NA, 'Amphibians'=NA)
models_qualityA = list('Mammals'=NA, 'Birds'=NA, 'Reptiles'=NA, 'Amphibians'=NA)
models_qualityB= list('Mammals'=NA, 'Birds'=NA, 'Reptiles'=NA, 'Amphibians'=NA)

##--------------------------------------------------------------------------------------------------

for(group in 1:length(DataList)) {
  load(file.path('Output','Models',paste0(names(DataList)[[group]],'_intercept_main.RData')))
  models_main[[group]] = m
  
  subsetA = DataList[[group]][DataList[[group]]$StudyQualityA == 1, ]
  subsetA$ControlIDsplit = factor(subsetA$ControlID, levels=unique(as.character(subsetA$ControlID)))
  subsetA_V = metafor::bldiag(lapply(split(subsetA,subsetA$ControlIDsplit), calc.vRRDelta))
  mA = rma.mv(logRR, V=subsetA_V, random=m$random, data=subsetA)
  models_qualityA[[group]] = mA
  
  subsetB = DataList[[group]][DataList[[group]]$StudyQualityB == 1, ]
  subsetB$ControlIDsplit = factor(subsetB$ControlID, levels=unique(as.character(subsetB$ControlID)))
  subsetB_V = metafor::bldiag(lapply(split(subsetB,subsetB$ControlIDsplit), calc.vRRDelta))
  mB = rma.mv(logRR, V=subsetB_V, random=m$random, data=subsetB)
  models_qualityB[[group]] = mB
  
}

##--------------------------------------------------------------------------------------------------
plotSensitivity <- function(model1,model3,model4,title){
  science_theme = theme(panel.grid.major = element_line(size = 0,colour="White"),
                        axis.title.y = element_text(vjust=1.25,size=8,family="Arial"),
                        axis.title.x = element_text(vjust=-0.3,size=8,family="Arial"),
                        axis.line = element_line(size = 0.7, color = "black"),
                        axis.text.y = element_text(size=8),
                        axis.text.x = element_text(size=8),
                        legend.title=element_blank(),
                        legend.position = "none",
                        plot.title = element_text(size=8))
  Data = data.frame(Type = c(paste('all (N=',model1$k,')',sep=''),
                             paste('a == 1 (N=',model3$k,')',sep=''),
                             paste('b == 1 (N=',model4$k,')',sep='')),
                    Estimate = c(model1$beta[1],model3$beta[1],model4$beta[1]),
                    CI_upper = c(model1$ci.ub,model3$ci.ub,model4$ci.ub),
                    CI_lower = c(model1$ci.lb,model3$ci.lb,model4$ci.lb))
  Data$Type = factor(Data$Type,levels=c(paste('b == 1 (N=',model4$k,')',sep=''),
                                        paste('a == 1 (N=',model3$k,')',sep=''),
                                        paste('all (N=',model1$k,')',sep='')))
  DataCI = data.frame(CI = c(Data$CI_upper,Data$CI_lower),
                      Type = c(Data$Type, Data$Type))
  plot <-  ggplot(Data,aes(x=Type,y=Estimate)) +
    ggtitle(title) + 
    geom_point(shape = 15, size=2) +
    geom_line(data = DataCI, aes(x=Type,y=CI, group=Type)) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", size=0.8) +    
    ylab(expression(LRR^Delta))+
    xlab("")+
    theme_bw()+
    science_theme
  return(plot)
}
pm=plotSensitivity(models_main[['Mammals']],models_qualityA[['Mammals']],models_qualityB[['Mammals']], 'a) Mammals')
pb=plotSensitivity(models_main[['Birds']],models_qualityA[['Birds']],models_qualityB[['Birds']], 'b) Birds')
pr=plotSensitivity(models_main[['Reptiles']],models_qualityA[['Reptiles']],models_qualityB[['Reptiles']], 'c) Reptiles')
pa=plotSensitivity(models_main[['Amphibians']],models_qualityA[['Amphibians']],models_qualityB[['Amphibians']],'d) Amphibians')

qualityFig = ggarrange(pm, pb, pr, pa,
                      ncol=2, nrow=2, align='v')

cairo_pdf(file='Figures/FigS9.pdf', width =5, height = 4)
plot(qualityFig)
dev.off()
