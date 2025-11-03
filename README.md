#**--------Environmental suitability for scrub typhus--------**

Repository contains data, scripts, and key outputs from "Mapping the global distribution and environmental suitability for scrub typhus".

System Requirements

This code was developed and tested using:

R Version: 4.3.0

Operating System: Windows 10 x64 (build 19045)

Software Dependencies

Ensure the following R packages are installed: raster, dplyr, sf, data.table, ranger, parallel, doParallel, foreach, readxl, car, corrplot, tidyverse, pheatmap, mgcv, blockCV, pROC, caret, gbm


Non-Standard Hardware
No special hardware requirements beyond a standard computer with sufficient RAM to handle the scale of analysis.

## **Installation Guide**

### Instructions

Install R: Download and install R from CRAN.

Install Required Packages: Run the following R script to install the required packages:

`install.packages(c("raster", "dplyr", "sf", "data.table", "ranger", 
                   "parallel", "doParallel", "foreach", "readxl", 
                   "car", "corrplot", "tidyverse", "pheatmap", 
                   "mgcv", "blockCV", "pROC", "caret", "gbm"))`
The installation of R and the required packages typically takes less than 30 minutes on a standard desktop computer.

## **Instructions for Use**

Repository Structure
The repository is organized into the following main three folders:

### data/
Contains datasets used in the analysis

- **population_ratio/**: Population information used for population calculations, including country-wise agriculatural ratio, elderly ratio, <15yrs ratio obtained from Worldbank.

- **df_balanced.rds**: Main database for scrub typhus environmental suitability, including presence data and equal amount of absence data with the corresponding coordinations, countries, data source, and values of 30 types of covariates.

- **15_df_balanced.rds**: Database used for sensitivity analysis when buffer distance for obtaining covariates value changed from 30 km to 15km.

- **45_df_balanced.rds**: Database used for sensitivity analysis when buffer distance for obtaining covariates value changed from 30 km to 45km. 

Data Preparation Scripts
- **data_preparation.R**: R code to generate df_balanced.rds from raw data.

- **data_for_sensitivity_analysis.R**: R code to generate base data for sensitivity analysis from raw data.

- **data_for_prediction.R**: R code to generate covaraites data for prediction from raw data.

### scripts/
Contains all analysis code and scripts

**Core Modeling**:
- **weight_mite.R**: R code to generate weights for mite species suitability modeling

- **weight_scrub_typhus.R**: R code to generate weights for scrub typhus environmental suitability.

- **model_mite_suitability.R**: R code to model key mite species suitability.

- **RandomForest.R**: R code of main RandomForest model for environmental suitability

- **Maxent.R**: R code of Maxent method implementation.

**Model Comparison & Validation**:
- **model_comparison_GAM_GLM_BRT_RF.R**: R code to compare GLM, GAM, BRT, RF performance using spatial block cross-validation.

- **spatial_block_cv.R**: R code of spatial block cross-validation function.

- **MOP.R**: R code to generate Mobility-Oriented Parity (MOP) metric for extrapolation analysis.

**Uncertainty & Sensitivity Analysis**:
- **Uncertainty_100fold.R**: R code for uncertainty analysis using 100-fold bootstrap method.

- **Uncertainty_10spatial_block_fold.R**: R code for uncertainty analysis using 10 spatial block folds method.

- **sensitivity_analysis.R**: R code to conduct sensitivity analysis.

**Results Examination**:
- **Examination.R**: R code to compare population-weighted environmental suitability with reported incidence and generate Figure 3.

### outputs/
Contains analysis results

- **spatial_block_cv/**: Spatial block cross-validation results for the final RandomForest model.

- **weight/**: Weight results from modeling process.

- **model_comparison/**: Results comparing different models.

- **Maxent/**: Results from Maxent model.

#### Due to file size limitations, all raster files (including covariates and the main outcome used for prediction) are not included in this repository. These files can be made available upon request by contacting qian@tropmedres.ac)

For additional details, refer to the main manuscript and supplementary documentation provided.
