# Author: Melinda de Jonge
# Date: 24-09-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: melinda.dejonge@ru.nl

## -------------------------------------------------------------------------------------------------
#' Script to fit and extract model coefficients for different imputations methods for missing 
#' standard deviations. 
#' Methods tested:
#' - Assuming abundances follow poisson distribution
#' - No imputation, i.e. remove all observations with a missing standard deviation
#' - Bracken1992 method (based on ratio of means to standard deviations of complete observations)
#' - HotDeckNN approach (resampling based on similar observations, creates NN imputations to fit 
#'                      NN models, get mean and median model coefficients)
##--------------------------------------------------------------------------------------------------
library(metafor)
library(metagear)
library(foreach)
library(doParallel)
source(file.path('R','Functions.R'))

# Load data
load(file.path('Data','Processed','GroupedData_main.RData'))

ncores = 4 
nn = 100 # Number of HotDeck imputations

# Data frame for easy output
fitted_parameters = data.frame('intercept'= rep(NA,5),
                               'CI_up' = rep(NA,5),
                               'CI_low' = rep(NA,5))
row.names(fitted_parameters) = c('Poisson','None','Bracken','Hotdeck mean','Hotdeck median')

# Function that sets the controlSD to a single value per controlID, needed to create 
# variance-covariance matrices after HotDeck imputations
harmonizeControlSD_HD <- function(X) {
  sd = rep(X$ControlSD_imputed[1],nrow(X))
  return(sd)
}

##---- Model runs ----------------------------------------------------------------------------------
sink(file.path('Output','Results_SDImputations.txt'))

for(group in 1:length(DataList)) {
  # ---- No imputation -------------------------------------------------------------------
  # Load main model & get fitted coefficients
  load(file.path('Output','Models',paste0(names(DataList)[[group]],'_intercept_main.RData')))
  fitted_parameters[1,1]=m$b
  fitted_parameters[1,2]=m$ci.ub
  fitted_parameters[1,3]=m$ci.lb

  # ---- No imputation -------------------------------------------------------------------
  # Remove all observations without known SD
  Data_noImp = DataList[[group]][!is.na(DataList[[group]]$InfraSD) & !is.na(DataList[[group]]$ControlSD), ]
  # Calculate variance covariance matrix for this subset
  Data_noImp$ControlIDsplit = factor(Data_noImp$ControlID, 
                                     levels=unique(as.character(Data_noImp$ControlID)))
  Data_noImp_V = metafor::bldiag(lapply(split(Data_noImp,Data_noImp$ControlIDsplit),
                                        calc.vRRDelta))
  # Run model and extract fitted coefficients
  m_noImp = rma.mv(logRR, V=Data_noImp_V, random=m$random, data=Data_noImp)
  fitted_parameters[2,1] = m_noImp$b
  fitted_parameters[2,2] = m_noImp$ci.ub
  fitted_parameters[2,3] = m_noImp$ci.lb

  # ---- Bracken imputation --------------------------------------------------------------
  Data_Bracken = DataList[[group]]
  # Impute missing SDs for infrastructure and control sites together 
  dataImpute = data.frame(dev = c(Data_Bracken$InfraSD,Data_Bracken$ControlSD), 
                          abundance = c(Data_Bracken$InfraA,Data_Bracken$ControlA))
  SDdatBracken = impute_SD(dataImpute,'dev','abundance')
  Data_Bracken$InfraSD_imputed = SDdatBracken$dev[1:nrow(Data_Bracken)]
  Data_Bracken$ControlSD_imputed = SDdatBracken$dev[(nrow(Data_Bracken)+1):(2*nrow(Data_Bracken))]
  # Re-calculate effect sizes and variance-covariance matrix with the imputed SDs
  Data_Bracken = calcEffectSize(Data_Bracken)
  Data_Bracken$ControlIDsplit = factor(Data_Bracken$ControlID, 
                                       levels=unique(as.character(Data_Bracken$ControlID)))
  Data_Bracken_V = metafor::bldiag(lapply(split(Data_Bracken,Data_Bracken$ControlIDsplit), 
                                          calc.vRRDelta))
  # Run model and extract fitted coefficients
  m_Bracken = rma.mv(logRR, V=Data_Bracken_V, random=m$random, data=Data_Bracken)
  fitted_parameters[3,1] = m_Bracken$b
  fitted_parameters[3,2] = m_Bracken$ci.ub
  fitted_parameters[3,3] = m_Bracken$ci.lb
  
  # ---- HotDeck NN ----------------------------------------------------------------------
  # Create NN imputations of missing SDs 
  Data_HotDeck = DataList[[group]]
  dataImpute = data.frame(dev = c(Data_HotDeck$InfraSD,Data_HotDeck$ControlSD),
                          abundance = c(Data_HotDeck$InfraA,Data_HotDeck$ControlA))
  SDdatHotDeckNN = impute_SD(dataImpute,'dev','abundance','Hotdeck_NN',range=nn)
  # Create cluster of workers to run models in parallel
  cl = parallel::makeCluster(ncores)
  registerDoParallel(cl)
  # Loop over all nn generated imputations & fit corresponding models
  coefsHD <- foreach(n=1:nn, .combine='rbind') %dopar% {
    library(metafor)
    # Extract the imputed SDs corresponding to this run (n)
    Data_HotDeck$InfraSD_imputed = SDdatHotDeckNN[[n]]$dev[1:nrow(Data_HotDeck)]
    Data_HotDeck$ControlSD_imputed = SDdatHotDeckNN[[n]]$dev[(nrow(Data_HotDeck)+1):(2*nrow(Data_HotDeck))]
    Data_HotDeck$ControlIDsplit = factor(Data_HotDeck$ControlID, 
                                         levels=unique(as.character(Data_HotDeck$ControlID)))
    control_sd = lapply(split(Data_HotDeck,Data_HotDeck$ControlIDsplit),harmonizeControlSD_HD)
    Data_HotDeck$ControlSD_imputed = unlist(control_sd)
    # Re-calculate effect sizes and variance-covariance matrix with the imputed SDs
    Data_HotDeck = calcEffectSize(Data_HotDeck)
    Data_HotDeck_V = metafor::bldiag(lapply(split(Data_HotDeck,Data_HotDeck$ControlIDsplit), 
                                            calc.vRRDelta))
    # Run model and extract fitted coefficients
    coefs = c(NA,NA,NA) # Required in case the model for this specific imputation does not converge
    try({
      m_HD = rma.mv(logRR, V=Data_HotDeck_V, random=m$random, data=Data_HotDeck)
      coefs = c(m_HD$b,m_HD$ci.ub,m$ci.lb)
    }, silent=TRUE)
    coefs
  }
  stopImplicitCluster()
  parallel::stopCluster(cl)
  # Get the mean and median of the HotDeck models coefficients
  fitted_parameters[4,1] = mean(coefsHD[,1],na.rm=TRUE)
  fitted_parameters[4,2] = mean(coefsHD[,2],na.rm=TRUE)
  fitted_parameters[4,3] = mean(coefsHD[,3],na.rm=TRUE)
  fitted_parameters[5,1] = median(coefsHD[,1],na.rm=TRUE)
  fitted_parameters[5,2] = median(coefsHD[,2],na.rm=TRUE)
  fitted_parameters[5,3] = median(coefsHD[,3],na.rm=TRUE)
  
  # ---- Print output to text file -------------------------------------------------------
  cat(paste0(names(DataList)[[group]],': \n'))
  print(fitted_parameters)
  cat('\n\n\n')
}
sink()