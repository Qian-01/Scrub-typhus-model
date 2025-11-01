
# Load required libraries
library(raster)
library(dplyr)
library(sf)
library(data.table)
library(ranger)
library(parallel)
library(doParallel)
library(foreach)

df <- read.csv("humanoccurrence.csv")
df <- df %>% 
  rename(
    Water_bodies = land1, 
    Evergreen_Needleleaf_Forests = land2, 
    Evergreen_Broadleaf_Forests = land3,
    Deciduous_Needleleaf_Forests = land4,
    Deciduous_Broadleaf_Forests = land5,
    Mixed_Forests = land6,
    Closed_Shrublands = land7,
    Open_Shrublands = land8,
    Woody_Savannas = land9,
    Savannas = land10,
    Grasslands = land11,
    Permanent_Wetlands = land12,
    Croplands = land13,
    Urban_and_Builtup = land14,
    Cropland_Natural_Vegetation_Mosaics = land15,
    Snow_and_Ice = land16,
    Barren_or_Sparsely_Vegetated = land17
  )
df$mite<-pmax(df$ld_mean,df$ls_mean,df$lp_mean, na.rm = T)

#re-scale first
predictors_names<- c("Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
                     "Deciduous_Needleleaf_Forests","Deciduous_Broadleaf_Forests","Mixed_Forests",
                     "Closed_Shrublands","Open_Shrublands","Woody_Savannas",
                     "Savannas","Grasslands","Permanent_Wetlands",
                     "Croplands","Urban_and_Builtup","Cropland_Natural_Vegetation_Mosaics",
                     "Snow_and_Ice","Barren_or_Sparsely_Vegetated","ndvi","evi", "pop", "prec", 
                     "tmin", "tmax", "sp","ws","rh", "elevation", "urban_acc_30s",
                     "mite","rodent_richness")
numeric_cols <- predictors_names[sapply(df[, predictors_names], is.numeric)]
scaling_params <- lapply(df[, numeric_cols], function(x) {
  list(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
})

df[, numeric_cols] <- lapply(df[, numeric_cols], scale, center = TRUE, scale = TRUE)

df <- df[complete.cases(df[, predictors_names]), ]

# Read scaled stack raster dataframe
predict_data <- readRDS("predict_scaled_covariates_2020.rds")

# Separate occurrence data
df_0 <- df %>% filter(occurrence == 0)
df_1 <- df %>% filter(occurrence == 1)

# Set up parallel processing
n_cores <- 10
cl <- makeCluster(n_cores)

# Prepare prediction data chunks (same as before)
predict_data$.__index__ <- 1:nrow(predict_data)
chunk_size <- ceiling(nrow(predict_data) / 600)
data_chunks <- list()
for (i in 1:600) {
  start_idx <- (i-1) * chunk_size + 1
  end_idx <- min(i * chunk_size, nrow(predict_data))
  data_chunks[[i]] <- predict_data[start_idx:end_idx, ]
}

# Function to predict on chunks
predict_chunk <- function(chunk, rf_model) {
  # Write progress to file
  progress_msg <- paste("Processing chunk with indices from", min(chunk$.__index__), 
                        "to", max(chunk$.__index__), "at", Sys.time())
  cat(progress_msg, "\n", file = "progress.txt", append = TRUE)
  
  na_rows <- apply(chunk, 1, function(x) any(is.na(x)))
  
  predictions <- rep(NA, nrow(chunk))
  #se.fit <- rep(NA, nrow(chunk))
  
  if (any(!na_rows)) {
    valid_predictions <- predict(rf_model, chunk[!na_rows,], type = "response")
    predictions[!na_rows] <- valid_predictions$predictions[,2]
  }
  
  # Write completion to file
  completion_msg <- paste("Completed chunk", min(chunk$.__index__), "-", 
                          max(chunk$.__index__), "at", Sys.time())
  cat(completion_msg, "\n", file = "progress.txt", append = TRUE)
  
  gc()
  return(data.frame(index = chunk$.__index__, prob = predictions))
}

# Create template raster
template_raster <- raster("annual_prec_2020.tif")

#---------- Bootstrap loop - 100 iterations -----------------------
set.seed(123)  # For reproducibility
bootstrap_seeds <- sample(1:10000, 100)  # Generate 100 different seeds

for (boot_iter in 1:100) {
  cat("Starting bootstrap iteration", boot_iter, "of 50 at", Sys.time(), "\n")
  
  # Set seed for this iteration
  set.seed(bootstrap_seeds[boot_iter])
  
  # Create balanced dataset by sampling from absence data
  df_0_sampled <- df_0 %>% sample_n(nrow(df_1))
  df_balanced <- bind_rows(df_0_sampled, df_1) %>% sample_frac(1)
  
  # Convert to factor
  df_balanced$occurrence <- as.factor(df_balanced$occurrence)
  
  # Create formula
  formula_rf <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))
  
  # Train random forest model
  rf.fit <- ranger(formula_rf,
                   data = df_balanced,
                   num.trees = 900,
                   keep.inbag = TRUE,
                   probability = TRUE,
                   importance = 'impurity',
                   case.weights = df_balanced$weight)
  
  # Create progress file for this iteration
  progress_file <- paste0("progress_boot_", boot_iter, ".txt")
  if (file.exists(progress_file)) file.remove(progress_file)
  
  # Parallel prediction
  clusterExport(cl, varlist = c("rf.fit", "predict_chunk"), envir = environment())
  clusterEvalQ(cl, library(ranger))
  
  result <- parLapply(cl, data_chunks, function(chunk) predict_chunk(chunk, rf.fit))
  
  # Extract predictions
  RF_predictions <- unlist(lapply(result, function(x) x$prob))
  #RF_predictions_se <- unlist(lapply(result, function(x) x$se))
  
  # Create rasters
  RF_raster <- template_raster
  values(RF_raster) <- RF_predictions
  
  #se_raster <- template_raster
  #values(se_raster) <- RF_predictions_se
  
  # Save rasters with bootstrap iteration number
  writeRaster(RF_raster, 
              filename = paste0("RF_pred_bootstrap_", sprintf("%03d", boot_iter), ".tif"), 
              format = "GTiff", overwrite = TRUE)
  #writeRaster(se_raster, 
   #           filename = paste0("RF_pred_bootstrap_se_", sprintf("%03d", boot_iter), ".tif"), 
    #          format = "GTiff", overwrite = TRUE)
  
  # Clean up memory
  rm(rf.fit, RF_raster, RF_predictions,  result)
  gc()
  
  cat("Completed bootstrap iteration", boot_iter, "at", Sys.time(), "\n")
}

# Stop cluster
stopCluster(cl)

cat("All 100 bootstrap iterations completed!\n")
cat("Generated files:\n")
cat("- RF_pred_bootstrap_001.tif to RF_pred_bootstrap_100.tif (predictions)\n")

cat("Creating summary rasters...\n")

# Read all bootstrap prediction rasters
bootstrap_stack <- stack(paste0("RF_pred_bootstrap_", sprintf("%03d", 1:100), ".tif"))

# Calculate mean and standard deviation
mean_prediction <- calc(bootstrap_stack, mean, na.rm = TRUE)
#sd_prediction <- calc(bootstrap_stack, sd, na.rm = TRUE)

# Save summary rasters
writeRaster(mean_prediction, "RF_pred_bootstrap_mean.tif", format = "GTiff", overwrite = TRUE)
#writeRaster(sd_prediction, "RF_pred_bootstrap_sd.tif", format = "GTiff", overwrite = TRUE)

cat("Summary rasters created:\n")
cat("- RF_pred_bootstrap_mean.tif (mean across 100 bootstrap iterations)\n")
