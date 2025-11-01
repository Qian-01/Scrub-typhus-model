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
The repository is organized into the following main folders:

### data/
Contains datasets used in the analysis

- **population_ratio/**: Population ratios for population calculations

- **df_balanced.rds**: Main database for scrub typhus environmental suitability (30km buffer)

- **15_df_balanced.rds**: Database with 15km buffer distance

- **45_df_balanced.rds**: Database with 45km buffer distance

Data Preparation Scripts:
- **data_preparation.R**: Generates df_balanced.rds from raw data

- **data_for_sensitivity_analysis.R**: Prepares base data for sensitivity analysis

- **data_for_prediction.R**: Generates prediction data from raw data

### scripts/
Contains all analysis code and scripts

**Core Modeling**:
- **weight_mite.R**: Generates weights for mite species suitability modeling

- **weight_scrub_typhus.R**: Generates weights for scrub typhus environmental suitability

- **model_mite_suitability.R**: Models key mite species suitability

- **RandomForest.R**: Main RandomForest model for environmental suitability

- **Maxent.R**: Maxent method implementation

**Model Comparison & Validation**:
- **model_comparison_GAM_GLM_BRT_RF.R**: Compares GLM, GAM, BRT, RF performance using spatial block cross-validation

- **spatial_block_cv.R**: Spatial block cross-validation functions

- **MOP.R**: Generates Mobility-Oriented Parity (MOP) metric for extrapolation analysis

**Uncertainty & Sensitivity Analysis**:
- **Uncertainty_100fold.R**: Uncertainty analysis using 100-fold bootstrap

- **Uncertainty_10spatial_block_fold.R**: Uncertainty analysis using 10 spatial block folds

- **sensitivity_analysis.R**: Conducts sensitivity analysis

**Results Examination**:
- **Examination.R**: Compares population-weighted environmental suitability with reported incidence (generates Figure 3)

### outputs/
Contains analysis results

- **spatial_block_cv/**: Spatial block cross-validation results

- **weight/**: Weight results from modeling process

- **model_comparison/**: Results comparing different models

- **Maxent/**: Results from Maxent model


For additional details, refer to the main manuscript and supplementary documentation provided.
