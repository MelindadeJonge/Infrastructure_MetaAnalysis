# Infrastructure_MetaAnalysis 

This repository contains the data and scripts to run the analyses for the paper:  
de Jonge, M. M. J., Gallego-Zamorano, J., Huijbregts, M. A. J., Schipper, A. M., & Benítez‐López, A. (in press). The impact of linear infrastructure on terrestrial vertebrate populations: A trait-based approach. *Global Change Biology*

When using this data and code, please cite the abovementioned publication as well as the raw data:  
Jonge, M.M.J. de, Gallego-Zamorano, J., Huijbregts, M.A.J., Schipper, A.M. & Benitez Lopez, A. (2022). Data from: “The impact of linear infrastructure on terrestrial vertebrate populations: A trait-based approach”. DANS EASY [Dataset]. doi: 10.17026/dans-xcw-zyvh. 

## Data overview
### Raw data
- The raw data as extracted from peer-reviewed and grey literature can be found under 'Data/Raw/Dataset.csv'. This contains the sample size, mean abundance at infrastructure and control sites, standard deviation of the means, habitat type, infrastructure type, source, location, species, and quality of each observation.
- A list of all included species and their corresponding traits can be found in: 'Data/Raw/SpeciesList.csv'
- Crosswalks to link the habitat types and infrastructure types as reported in the database to their classifications used in the analysis are found in: 'Data/Raw/Crosswalk_Habitat.csv' and 'Data/Raw/Crosswalk_Infrastructure.csv'

### Processed data
- The processed data including the calculated effect sizes, i.e. the small sample size adjusted log response ratio, the corresponding sampling variances, the species traits and standardizes habitat and infrastructure type used for all main analyses is stored in 'Data/Processed/Dataset_main.csv'. 
- The 'Data/Processed/GroupedData_main.RData' contains the same information but ordered and split to species group as a list object.  
- The variance-covariance matrices calculated to account for repeated controls in the main analyses are stored as a list object in 'Data/Processed/VarCovarMatrices_main.RData'.
- All processed data for the main analyses was created with 'R/PrepareData.R'. 

### Results/Output
- Fitted models are stored in 'Output/Models' for use in further applications in environmental impacts assessments. These are stored as rma.mv objects and require the R package 'metafor' for loading and application. 
