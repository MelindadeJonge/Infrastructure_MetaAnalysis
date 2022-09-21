# Author: Melinda de Jonge
# Date: 20-09-2022
# Department of Environmental Sciences, Radboud University Nijmegen
# e-mail: melinda.dejonge@ru.nl

## ---- Script description -------------------------------------------------------------------------
#' Script to prepare datasets for main analyses from raw data. 
#' 1. Remove control points (the full dataset contains these as separate observations)
#' 2. Impute missing standard deviations and calculate effect sizes
#' 3. Add potential moderators: distance, body mass, diet, infrastructure type & habitat type
#' 4. Split the dataset into the four species groups for analyses
#' 5. Calculate the variance-covariance matrix to account for repeated controls for each group
#' 
#' The datasets for each group are stored together in a list as 'Data/Processed/GroupedData_main.RData'.
#' Variance-coverance matrices are stored in the same manner as 'Data/Processed/VarCovarMatrices_main_RData. 
#' 
#' This script also creates worldmaps indicating the locations of the observatoins used in the 
#' analyses for species group and split to infrastructure type. 
##--------------------------------------------------------------------------------------------------

library(ggplot2)
library(viridis)
library(ggpubr)
source(file.path('R','Functions.R'))

# Load raw data
X = read.csv(file.path('Data','Raw','Full_dataset.csv'))
speciesList = read.csv(file.path('Data','Raw','SpeciesList.csv'))

# Add unique ID for each observation
X$ID = paste0('S',rownames(X))

## ---- Remove control points ----------------------------------------------------------------------
# Set distances where only min and max distances are available
X$Dist[is.na(X$Dist)] = (X$MaxDist[is.na(X$Dist)] + X$MinDist[is.na(X$Dist)]) / 2
X$CDist[is.na(X$CDist)] = (X$MaxCDist[is.na(X$CDist)] + X$MinCDist[is.na(X$CDist)]) / 2

#Remove controls
indexControls = rowSums(cbind2((X$CDist == X$Dist), cbind2((X$MaxCDist == X$MaxDist),
                                                          (X$MinCDist == X$MinDist))),
                        na.rm=T) != 0
X = X[!indexControls, ]
print(paste('Number of control points removed:',sum(indexControls)))

## ----- Impute missing standard deviations --------------------------------------------------------
# Set SD to NA when SD == 0 
X$InfraSD[X$InfraSD == 0] = NA
X$ControlSD[X$ControlSD == 0] = NA
X$InfraSD_imputed = X$InfraSD
X$ControlSD_imputed = X$ControlSD

X$InfraSD_imputed[is.na(X$InfraSD)] = sqrt(X$InfraA[is.na(X$InfraSD)] * X$InfraD[is.na(X$InfraSD)]) / 
  X$InfraD[is.na(X$InfraSD)]
X$ControlSD_imputed[is.na(X$ControlSD)] = sqrt(X$ControlA[is.na(X$ControlSD)] * X$ControlD[is.na(X$ControlSD)]) /
  X$ControlD[is.na(X$ControlSD)]

## ---- Calculate effect sizes ---------------------------------------------------------------------
X$logRR_uncorrected = log(X$InfraA / X$ControlA)
X$logRR = X$logRR_uncorrected + 0.5 * (X$InfraSD_imputed^2 / (X$InfraN * X$InfraA^2) - 
                                       X$ControlSD_imputed^2 / (X$ControlN * X$ControlA^2))

X$VAR_uncorrected = (X$InfraSD_imputed^2) / (X$InfraN * X$InfraA^2) + 
  (X$ControlSD_imputed^2) / (X$ControlN * X$ControlA^2)
X$VAR = X$VAR_uncorrected + 0.5 * (X$InfraSD_imputed^4 / (X$InfraN^2 * X$InfraA^4) + 
                                   X$ControlSD_imputed^4 / (X$ControlN^2 * X$ControlA^4))

## ---- Prepare moderators -------------------------------------------------------------------------
# Log distance
X$logD = log10(X$Dist)
X$logD2 = X$logD^2

# Harmonize infra types
InfraTypes = read.csv(file.path('Data','Raw','Crosswalk_Infrastructure.csv'))
X$InfraTypeOriginal = X$InfraType
X$InfraType = InfraTypes[match(X$InfraType, InfraTypes$in_database),'category']
X = X[!is.na(X$InfraType), ]

# Harmonize habitat types
Habitat = read.csv(file.path('Data','Raw','Crosswalk_Habitat.csv'))
X$HabitatOriginal = X$Habitat
X$Habitat = Habitat[match(X$Habitat, Habitat$in_database),'category']
X = X[!is.na(X$Habitat), ]

# Get species traits
X$Order = speciesList$Order[match(X$SpeciesCode, speciesList$Species_code)] 
X$Genus = speciesList$Genus[match(X$SpeciesCode, speciesList$Species_code)] 
X$Family = speciesList$Family[match(X$SpeciesCode, speciesList$Species_code)] 
X$Diet = speciesList$Diet[match(X$SpeciesCode, speciesList$Species_code)] 
X$BM = speciesList$Body_mass_g[match(X$SpeciesCode, speciesList$Species_code)] 
X$logBM = log10(X$BM)

# Save full dataset as csv
write.csv(X, file=file.path('Data','Processed','Dataset_main.csv'))

## ---- Split groups -------------------------------------------------------------------------------
Birds = subset(X, Group == 'Birds') 
Mammals = subset(X, Group == 'Mammals') 
Reptiles = subset(X, Group == 'Reptiles') 
Amphibians = subset(X, Group == 'Amphibians') 

# Set infrastructure type for powerlines to non-traffic for mammals, reptiles and amphibians
Mammals$InfraType[Mammals$InfraType == 'Powerline'] = 'Non-traffic'
Reptiles$InfraType[Reptiles$InfraType == 'Powerline'] = 'Non-traffic'
Amphibians$InfraType[Amphibians$InfraType == 'Powerline'] = 'Non-traffic'

## ---- Variance-Covariance matrices ---------------------------------------------------------------
# Make sure all datapoints are in order of their controlID, otherwise the variance covariance matrix
# won't match with the observations in the dataset.
Mammals = Mammals[order(Mammals$ControlID), ]
Mammals$ControlIDsplit = factor(Mammals$ControlID, levels=unique(as.character(Mammals$ControlID)))
Birds = Birds[order(Birds$ControlID), ]
Birds$ControlIDsplit = factor(Birds$ControlID, levels=unique(as.character(Birds$ControlID)))
Reptiles = Reptiles[order(Reptiles$ControlID), ]
Reptiles$ControlIDsplit = factor(Reptiles$ControlID, levels=unique(as.character(Reptiles$ControlID)))
Amphibians = Amphibians[order(Amphibians$ControlID), ]
Amphibians$ControlIDsplit = factor(Amphibians$ControlID, levels=unique(as.character(Amphibians$ControlID)))

# Calculate variance covariance matrices
Mammals_V = metafor::bldiag(lapply(split(Mammals,Mammals$ControlIDsplit), calc.vRRDelta))
Birds_V = metafor::bldiag(lapply(split(Birds,Birds$ControlIDsplit), calc.vRRDelta))
Reptiles_V = metafor::bldiag(lapply(split(Reptiles,Reptiles$ControlIDsplit), calc.vRRDelta))
Amphibians_V = metafor::bldiag(lapply(split(Amphibians,Amphibians$ControlIDsplit), calc.vRRDelta))

## ---- Save data for further analyses ------------------------------------------------------------- 
DataList = list(Mammals=Mammals,
                Birds=Birds,
                Reptiles=Reptiles,
                Amphibians=Amphibians)

VarCovarList = list(Mammals=Mammals_V,
                    Birds=Birds_V,
                    Reptiles=Reptiles_V,
                    Amphibians=Amphibians_V)

save(VarCovarList, file=file.path('Data','Processed','VarCovarMatrices_main.RData'))
save(DataList, file=file.path('Data','Processed','GroupedData_main.RData'))

## ---- Create Fig. 1 (Maps) -----------------------------------------------------------------------
MammalMap = createMap(Mammals, title='a) Mammals')
BirdMap = createMap(Birds, title='b) Birds')
ReptileMap = createMap(Reptiles, title='c) Reptiles')
AmphibianMap = createMap(Amphibians, title='d) Amhibians')

leg = get_legend(MammalMap)

plotAll = ggarrange(MammalMap, BirdMap ,ReptileMap, AmphibianMap,
                    ncol=2, nrow=2,
                    widths=c(1,1), heights=c(1,1), 
                    common.legend=FALSE, legend='bottom', legend.grob=leg)

cairo_pdf(file=file.path('Figures','Fig1.pdf'),
          width=5, height=3)
plot(plotAll)
dev.off()

png(file=file.path('Figures','Fig1.png'),
    width=5, height=3, units='in', res=600)
plot(plotAll)
dev.off()