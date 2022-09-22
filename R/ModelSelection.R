# Author: Melinda de Jonge
# Date: 20-09-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: melinda.dejonge@ru.nl

## -------------------------------------------------------------------------------------------------
#' 
##--------------------------------------------------------------------------------------------------

library(metafor)
library(doParallel)
library(foreach)

source(file.path('R','Functions.R'))
source(file.path('R','FixedEffects.R'))

ncores = 4

# Load prepared data
load(file.path('Data','Processed','GroupedData_main.RData'))
load(file.path('Data','Processed','VarCovarMatrices_main.RData'))

# Specify random & fixed effect terms
REList = list(Mammals=list(~1|ID, ~1|Order/SpeciesCode, ~1|Source/Study),
              Birds=list(~1|ID, ~1|Order/SpeciesCode, ~1|Source/Study),
              Reptiles=list(~1|ID, ~1|Family/SpeciesCode, ~1|Source/Study),
              Amphibians=list(~1|ID, ~1|Family/SpeciesCode, ~1|Source/Study))

FEListAll = list(FEMB, FEMB, FERep, FEAmph)

## ----- Model Selection & best model fit ----------------------------------------------------------
for(group in 1:length(names(DataList))) {
  print(paste0('Starting fixed effects selection for ',names(DataList)[group]))
  FEList = FEListAll[[group]]

  # Create cluster of workers to run models in parallel
  cl = parallel::makeCluster(ncores)
  registerDoParallel(cl)

  # Fit all fixed effect structures to extract fit statistics and fitted parameters
  Results <- foreach(i=1:length(FEList), .combine='rbind') %dopar% {
    library(metafor)
    FixedEffects = Reduce(paste, deparse(FEList[[i]][[1]]))
    stats=rep(NA, 13)
    try({
      model <- rma.mv(logRR, mods=FEList[[i]][[1]],
                      V=VarCovarList[[group]],
                      random=REList[[group]],
                      data=DataList[[group]])
      FitStats = model$fit.stats
      R2 = R2.func(model)
      stats = c(FitStats['AICc','ML'],FitStats['BIC','ML'],
                model$QM,model$QMp,model$QE,model$QEp,
                round(model$sigma2,4),round(R2,4))
    }, silent=FALSE)
    c(FixedEffects,stats)
  }
  stopImplicitCluster()
  parallel::stopCluster(cl)

  Results = data.frame(Results)
  colnames(Results) = c('FixedEffects','AICc','BIC','Qm','Qmp','QE','QEp',
                        'Sigma2_ID','Sigma2_Order_Family','Sigma2_Species',
                        'Sigma2_Source','Sigma2_Study','R2m','R2c')
  Results$AICc = as.numeric(Results$AICc)
  Results = Results[order(Results$AICc), ]
  Results$dAICc = Results$AICc - Results$AICc[1]
  write.csv(Results, file=file.path('Output',
                                    paste0('ModelSelectionAICc_',names(DataList)[group],'.csv')),
            row.names=FALSE)

  # Refit best model and save in 'Output/Models/'
  final_model <- rma.mv(logRR,mods=formula(Results$FixedEffects[1]),
                        V=VarCovarList[[group]],
                        random=REList[[group]],
                        data=DataList[[group]])
  save(final_model,file=file.path('Output','Models',
                                  paste0('FinalModel',names(DataList)[group],'.RData')))
}


