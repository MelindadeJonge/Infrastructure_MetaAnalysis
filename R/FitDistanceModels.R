# Author: Melinda de Jonge
# Date: 20-09-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: melinda.dejonge@ru.nl

## -------------------------------------------------------------------------------------------------
#'
##--------------------------------------------------------------------------------------------------

library(metafor)
library(ggplot2)
library(ggpubr)
source(file.path('R','Functions.R'))

# Load prepared data
load(file.path('Data','Processed','GroupedData_main.RData'))
load(file.path('Data','Processed','VarCovarMatrices_main.RData'))

# Specify random & fixed effect terms
REList = list(Mammals=list(~1|ID, ~1|Order/SpeciesCode, ~1|Source/Study),
              Birds=list(~1|ID, ~1|Order/SpeciesCode, ~1|Source/Study),
              Reptiles=list(~1|ID, ~1|Family/SpeciesCode, ~1|Source/Study),
              Amphibians=list(~1|ID, ~1|Family/SpeciesCode, ~1|Source/Study))

FEList = list(c(~1),
              c(~logD),
              c(~logD + logD2))

# Specify model names for output files
modelNames= c('intercept','distanceLinear','distanceQuadratic')

## ----- Fit models --------------------------------------------------------------------------------
for(group in 1:length(DataList)) {
  for(i in 1:3) {
    m = rma.mv(logRR, V=VarCovarList[[group]],
               mods=FEList[[i]][[1]],
               random=REList[[group]],
               data=DataList[[group]])
    filename = file.path('Output','Models',
                         paste0(names(DataList)[[group]],'_',modelNames[i],'_main.RData'))
    save(m, file=filename)
  }
}

## ----- Plot quadratic models ---------------------------------------------------------------------
# Load the quadratic model for each group
load(file.path('Output','Models','Mammals_distanceQuadratic_main.RData'))
mammal_modelDistance = m
load(file.path('Output','Models','Birds_distanceQuadratic_main.RData'))
bird_modelDistance = m
load(file.path('Output','Models','Reptiles_distanceQuadratic_main.RData'))
reptile_modelDistance = m
load(file.path('Output','Models','Amphibians_distanceQuadratic_main.RData'))
amphibian_modelDistance = m

# Plotting
pm = distancePlots(mammal_modelDistance,
                   data=DataList[[1]], title='a) Mammals')
pb = distancePlots(bird_modelDistance,
                   data=DataList[[2]], title='b) Birds')
pr = distancePlots(reptile_modelDistance,
                   data=DataList[[3]], title='c) Reptiles')
pa = distancePlots(amphibian_modelDistance,
                   data=DataList[[4]], title='d) Amphibians')

fig2 = ggarrange(pm, pb, pr, pa,
                 ncol=2, nrow=2, align='v')

png(filename=file.path('Figures','Fig2.png'), width=4, height=4,
    units='in',res=600)
plot(fig2)
dev.off()

cairo_pdf(file=file.path('Figures','Fig2.pdf'), width=4, height=4)
plot(fig2)
dev.off()

