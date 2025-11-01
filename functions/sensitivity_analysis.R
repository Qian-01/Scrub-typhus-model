
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
library(raster)
library(dplyr)
library(sf)
library(data.table)
library(ranger)
library(parallel)
library(doParallel)
library(foreach)
#read scaled stack raster dataframe
predict_data<-readRDS("predict_scaled_covariates_2020.rds")

library(parallel)
n_cores <- 10  
cl <- makeCluster(n_cores)
predict_data$.__index__ <- 1:nrow(predict_data) 

chunk_size <- ceiling(nrow(predict_data) / 600)
data_chunks <- list()
for (i in 1:600) {
  start_idx <- (i-1) * chunk_size + 1
  end_idx <- min(i * chunk_size, nrow(predict_data))
  data_chunks[[i]] <- predict_data[start_idx:end_idx, ]
}
# Create progress file
progress_file <- "progress.txt"
if (file.exists(progress_file)) file.remove(progress_file)

predict_chunk <- function(chunk) {
  progress_msg <- paste("Processing chunk with indices from", min(chunk$.__index__), 
                        "to", max(chunk$.__index__), "at", Sys.time())
  cat(progress_msg, "\n", file = "progress.txt", append = TRUE)
  
  na_rows <- apply(chunk, 1, function(x) any(is.na(x)))
  
  predictions <- rep(NA, nrow(chunk))
  se.fit <- rep(NA, nrow(chunk))
  
  if (any(!na_rows)) {
    valid_predictions <- predict(rf.fit, chunk[!na_rows,], type = "response")
    predictions[!na_rows] <- valid_predictions$predictions[, 2]
    #se.fit[!na_rows] <- valid_predictions[[6]][,2]
  }
  # Write completion to file
  completion_msg <- paste("Completed chunk", min(chunk$.__index__), "-", 
                          max(chunk$.__index__), "at", Sys.time())
  cat(completion_msg, "\n", file = "progress.txt", append = TRUE)
  gc() 
  return(data.frame(index = chunk$.__index__, prob = predictions))
}

clusterExport(cl, varlist = c("rf_weight_nobio", "predict_chunk"))
clusterEvalQ(cl, library(ranger)) 

result <- parLapply(cl, data_chunks, predict_chunk)

stopCluster(cl)
gc()

RF_predictions <- unlist(lapply(result, function(x) x$prob))

template_raster <- raster("annual_prec_2020.tif")
RF_raster <- template_raster
values(RF_raster) <- RF_predictions

writeRaster(RF_raster, filename = "RF_pred_nobio_weight.tif", format = "GTiff", overwrite = TRUE)
#writeRaster(se_raster, filename = "RF_pred_nobio_weight_se.tif", format = "GTiff", overwrite = TRUE)

