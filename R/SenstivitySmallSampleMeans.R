# Author: Melinda de Jonge
# Date: 22-09-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: melinda.dejonge@ru.nl

## -------------------------------------------------------------------------------------------------
#' 
##--------------------------------------------------------------------------------------------------
library(metafor)
source(file.path('R','Functions.R'))

load(file.path('Data','Processed','GroupedData_main.RData'))

##--------------------------------------------------------------------------------------------------
fitted_parameters = data.frame('intercept'=c(NA,NA),
                               'CI_up' = c(NA,NA),
                               'CI_low' = c(NA,NA),
                               'QE' = c(NA,NA),
                               'QEp' = c(NA,NA))
row.names(fitted_parameters) = c('AllData','SelectionLargeMeans')

sink(file.path('Output','Results_GearySelection.txt'))
for(group in 1:length(DataList)) {
  # Main intercept model with all data
  load(file.path('Output','Models',paste0(names(DataList)[[group]],'_intercept_main.RData')))
  
  # Get subset of data for which Geary's diagnostic > 3
  DataList[[group]]$GearyInfra = DataList[[group]]$InfraA / DataList[[group]]$InfraSD_imputed * ( (4*DataList[[group]]$InfraN^(3/2)) / (1+4*DataList[[group]]$InfraN) )
  DataList[[group]]$GearyControl = DataList[[group]]$ControlA / DataList[[group]]$ControlSD_imputed * ( (4*DataList[[group]]$ControlN^(3/2)) / (1+4*DataList[[group]]$ControlN) )
  DataList[[group]]$GearyMin = apply(data.frame(DataList[[group]]$GearyInfra,DataList[[group]]$GearyControl),1,min)
  smallMeans = DataList[[group]]$GearyMin <= 3 | is.na(DataList[[group]]$GearyMin)
  subsetLargeMeans = DataList[[group]][!smallMeans, ]

  # Create new variance covariance matrix for the selected subset
  subsetLargeMeans$ControlIDsplit = factor(subsetLargeMeans$ControlID, levels=unique(as.character(subsetLargeMeans$ControlID)))
  subsetLargeMeans_V = metafor::bldiag(lapply(split(subsetLargeMeans,subsetLargeMeans$ControlIDsplit), calc.vRRDelta))

  #Fit intercept model and extract parameters for the subset
  mGearySel = rma.mv(logRR, V=subsetLargeMeans_V, random=m$random, data=subsetLargeMeans)

  # Print parameters and Qc to text file
  fitted_parameters[1,1]=m$b
  fitted_parameters[1,2]=m$ci.ub
  fitted_parameters[1,3]=m$ci.lb
  fitted_parameters[1,4]=m$QE
  fitted_parameters[1,5]=m$QEp
  fitted_parameters[2,1]=mGearySel$b
  fitted_parameters[2,2]=mGearySel$ci.ub
  fitted_parameters[2,3]=mGearySel$ci.lb
  fitted_parameters[2,4]=mGearySel$QE
  fitted_parameters[2,5]=mGearySel$QEp

  cat(paste0(names(DataList)[[group]],': \n'))
  print(fitted_parameters)
  cat('\n\n\n')
}
sink()
