# Author: Melinda de Jonge
# Date: 20-09-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: melinda.dejonge@ru.nl

## -------------------------------------------------------------------------------------------------
#' Script to test for publication bias by performing Egger tests and creating funnel plots (Fig S8)
#' Output from Egger tests is stored in 'Output/EggerTests.txt'
##--------------------------------------------------------------------------------------------------

library(metafor)
library(ggplot2)
library(ggpubr)

load(file.path('Data','Processed','GroupedData_main.RData'))
load(file.path('Data','Processed','VarCovarMatrices_main.RData'))

intercepts = c('Mammals'=NA, 'Birds'=NA, 'Reptiles'=NA, 'Amphibians'=NA)

##---- Egger tests ---------------------------------------------------------------------------------
sink(file.path('Output','EggerTests.txt'))
for(group in 1:length(DataList)) {
  load(file.path('Output','Models',paste0(names(DataList)[[group]],'_intercept_main.RData')))

  # Add residuals and precision to dataframes
  DataList[[group]]$residnull = residuals(m)
  DataList[[group]]$invse = 1/sqrt(m$vi)

  # Egger test
  m_egger <- rma.mv(residnull, V=VarCovarList[[group]], mod=invse,
                     random=m$random,
                     data= DataList[[group]], intercept=TRUE, method="REML", level=95, digits=4)
  
  cat(paste0(names(DataList)[[group]],': \n'))
  print(summary(m_egger))
  cat('\n\n\n')

  # Extract intercept to use in funnel plots later
  intercepts[group] = m$b[,1]

}
sink()

##---- Funnel plots ---------------------------------------------------------------------------------

funnelGG <- function(data,title,nullest) {
  science_theme = theme(panel.grid.major = element_line(size=0,colour="White"),
                        axis.title.y = element_text(vjust=1.25,size=8,family="Arial"),
                        axis.title.x = element_text(vjust=-0.3,size=8,family="Arial"),
                        axis.line = element_line(size=0.7, color="black"),
                        axis.text.y = element_text(size=8),
                        axis.text.x = element_text(size=8),
                        legend.title=element_blank(),
                        legend.position = 'none',
                        plot.title = element_text(size=8))
  ylimit = max(c(max(data$residnull),abs(min(data$residnull))))
  plot <-  ggplot(data) +
    ggtitle(title) +
    geom_point(aes(x=invse,y=residnull),shape=20,col='slateblue4',alpha=I(.1)) +
    scale_shape_identity() +
    ylab(expression(LRR^Delta)) +
    xlab("Precision (1/SE)") +
    ylim(c(-ylimit,ylimit)) +
    geom_hline(yintercept=nullest, size=0.8, linetype="dashed", color="grey40") +
    geom_hline(yintercept=0, size=0.8, linetype="dashed", color="grey") +
    theme_bw() +
    science_theme
  return(plot)
}
pm = funnelGG(DataList[['Mammals']],title='a) Mammals',nullest=intercepts[1])
pb = funnelGG(DataList[['Birds']],title='b) Birds',nullest=intercepts[2])
pr = funnelGG(DataList[['Reptiles']],title='c) Reptiles',nullest=intercepts[3])
pa = funnelGG(DataList[['Amphibians']],title='d) Amphibians',nullest=intercepts[4])

funnelFig = ggarrange(pm, pb, pr, pa,
                 ncol=2, nrow=2, align='v')

cairo_pdf(file='Figures/FigS8.pdf', width = 4, height = 4)
plot(funnelFig)
dev.off()
