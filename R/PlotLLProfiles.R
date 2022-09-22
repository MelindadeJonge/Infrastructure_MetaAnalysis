# Author: Melinda de Jonge
# Date: 21-09-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: melinda.dejonge@ru.nl

## -------------------------------------------------------------------------------------------------
#'
##--------------------------------------------------------------------------------------------------

library(metafor)
library(ggpubr)
library(grid)
library(ggplotify)
library(gridGraphics)

ncores = 2

#Load fitted models
load(file.path('Output','Models','FinalModelMammals.RData'))
mammal_model = final_model
load(file.path('Output','Models','FinalModelBirds.RData'))
bird_model = final_model
load(file.path('Output','Models','FinalModelReptiles.RData'))
reptile_model = final_model
load(file.path('Output','Models','FinalModelAmphibians.RData'))
amphibian_model = final_model
modelList = list(Mammals=mammal_model,
                 Birds=bird_model,
                 Reptiles=reptile_model,
                 Ammphibians=amphibian_model)

REList = list(Mammals=c('(~1|ID)','(~1|Order)','(~1|Order/Species)','(~1|Source)','(~1|Source/Study)'),
              Birds=c('(~1|ID)','(~1|Order)','(~1|Order/Species)','(~1|Source)','(~1|Source/Study)'),
              Reptiles=c('(~1|ID)','(~1|Family)','(~1|Family/Species)','(~1|Source)','(~1|Source/Study)'),
              Amphibians=c('(~1|ID)','(~1|Family)','(~1|Family/Species)','(~1|Source)','(~1|Source/Study)'))

##--------------------------------------------------------------------------------------------------
profilePlot = function(profile,title,RE) {
  pdf('Rplots.pdf')
  dev.control(displaylist="enable")
  par(mfrow=c(2,3), mar=c(6,4,0.1,1))
  for(i in 1:profile$comps) {
    ymin = profile[[i]]$ylim[1]
    ymax = profile[[i]]$ylim[2]
    rounding = ifelse(ymax - ymin < 5, 1, 0)
    xticks = seq(from=round(ymin - (rounding * 0.1), rounding),
                 to=round(ymax + 0.1*rounding, rounding), 
                 by=round( ((ymax + 0.1*rounding) - (ymin-0.1*rounding)) / 3, rounding))
    
    plot(profile[[i]]$ll ~ profile[[i]]$sigma,
         type='o', pch=16, bty='l', ylim=c(ymin - 0.1,ymax + 0.2),
         xlab=bquote(sigma[.(i)]^2 ~ .(RE[i])),
         ylab='', yaxt='n',cex.lab=1, main='')
    axis(side=2, at=xticks)
    abline(v=profile[[i]]$vc, lty="dotted")
    abline(h=profile[[i]]$maxll, lty="dotted")
  }
  RecP = recordPlot()
  dev.off()
  basePlot = grid.grabExpr({grid.echo(RecP)})
  
  ggProfPlot=ggarrange(NULL, as.ggplot(basePlot),
                       ncol=1, nrow=2, heights = c(0.075,1), align='hv',
                       labels = title,
                       font.label=list(size=8, face='plain'), vjust=2, hjust=-0.25)
  ggProfPlot = annotate_figure(ggProfPlot,
                               left = text_grob('Restricted log-likelihood \t\t\t\t Restricted log-likelihood',
                                                rot=90, size=9, x=1))
  return(ggProfPlot)
}

##--------------------------------------------------------------------------------------------------

for(group in 1:length(names(REList))) {
  print(paste0('Creating profile likelihood plots for ',names(REList)[group]))
  LLprof = profile(modelList[[group]], parallel="snow",ncpus=ncores,plot=FALSE)
  
  profPlot = profilePlot(LLprof,
                         RE=REList[[group]],
                         title=paste0('Profile plots for ',names(REList)[group]))
  # Save plot
  cairo_pdf(file=file.path('Figures',paste0('FigS',group+3,'.pdf')),
            width = 7.5, height = 5)
  plot(profPlot)
  dev.off()
}


