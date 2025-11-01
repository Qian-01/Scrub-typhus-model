
#---------------------15km,45km buffer------------------------
library(ranger)

# Define predictor variables (make sure this matches your previous definition)
predictors_names <- c(
  "Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
  "Deciduous_Needleleaf_Forests", "Deciduous_Broadleaf_Forests", "Mixed_Forests",
  "Closed_Shrublands", "Open_Shrublands", "Woody_Savannas", "Savannas", "Grasslands",
  "Permanent_Wetlands", "Croplands", "Urban_and_Builtup", "Cropland_Natural_Vegetation_Mosaics",
  "Snow_and_Ice", "Barren_or_Sparsely_Vegetated", "ndvi", "evi", "pop", "prec", 
  "tmin", "tmax", "sp", "ws", "rh", "elevation", "urban_acc_30s", "mite", "rodent_richness"
)

# Function to train and save RF model
train_rf_model <- function(data_file, output_file) {
  # Load balanced data
  df_balanced <- readRDS(data_file)
  # Prepare formula and data
  formula_rf <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))
  df_balanced$occurrence <- as.factor(df_balanced$occurrence)
  
  # Train random forest
  set.seed(123)
  rf_weight <- ranger(
    formula = formula_rf,
    data = df_balanced,
    num.trees = 900,
    keep.inbag = TRUE,
    probability = TRUE,
    importance = 'impurity',
    case.weights = df_balanced$weight
  )
  
  # Save model
  saveRDS(rf_weight, output_file)
  cat("Model saved to:", output_file, "\n")
  
  return(rf_weight)
}

# Train models for both 15km and 45km datasets
rf_15 <- train_rf_model("15_df_balanced_new.rds", "15_rf_weight.rds")
rf_45 <- train_rf_model("45_df_balanced_new.rds", "45_rf_weight.rds")

# Load the saved models if needed
rf_15_loaded <- readRDS("15_rf_weight.rds")
rf_45_loaded <- readRDS("45_rf_weight.rds")

#predict


#------------------spatial-thinned data sensitivity analysis-----------------------
# Train models for spatial-thinned data
rf_thinned <- train_rf_model("df_thinned.rds", "rf_thinned.rds")
rf_thinned <- readRDS("rf_thinned.rds")

#predict


#------------------without weight sensitivity analysis-----------------------
df_balanced <- readRDS("df_balanced.rds")
formula_rf <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))
df_balanced$occurrence <- as.factor(df_balanced$occurrence)
  
# Train random forest
set.seed(123)
rf_noweight <- ranger(
  formula = formula_rf,
  data = df_balanced,
  num.trees = 900,
  keep.inbag = TRUE,
  probability = TRUE,
  importance = 'impurity'
)
saveRDS(rf_noweight,"rf_noweight.rds")
#predict


#------------------without bio info sensitivity analysis-----------------------
# update predictors
predictors_names_nobio <- c(
  "Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
  "Deciduous_Needleleaf_Forests", "Deciduous_Broadleaf_Forests", "Mixed_Forests",
  "Closed_Shrublands", "Open_Shrublands", "Woody_Savannas", "Savannas", "Grasslands",
  "Permanent_Wetlands", "Croplands", "Urban_and_Builtup", "Cropland_Natural_Vegetation_Mosaics",
  "Snow_and_Ice", "Barren_or_Sparsely_Vegetated", "ndvi", "evi", "pop", "prec", 
  "tmin", "tmax", "sp", "ws", "rh", "elevation", "urban_acc_30s"
)
formula_rf_nobio <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))

# Train random forest
set.seed(123)
rf_weight_nobio <- ranger(
  formula = formula_rf_nobio,
  data = df_balanced,
  num.trees = 900,
  keep.inbag = TRUE,
  probability = TRUE,
  importance = 'impurity',
  case.weights = df_balanced$weight
)
  
# Save model
saveRDS(rf_weight_nobio, "rf_weight_nobio.rds")

#predict

