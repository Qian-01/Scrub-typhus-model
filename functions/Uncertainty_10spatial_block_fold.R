
#load library

library(raster)
library(dplyr)
library(sf)
library(data.table)
library(ranger)
library(parallel)
library(doParallel)
library(foreach)

predictors_names<- c("Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
                     "Deciduous_Needleleaf_Forests","Deciduous_Broadleaf_Forests","Mixed_Forests",
                     "Closed_Shrublands","Open_Shrublands","Woody_Savannas",
                     "Savannas","Grasslands","Permanent_Wetlands",
                     "Croplands","Urban_and_Builtup","Cropland_Natural_Vegetation_Mosaics",
                     "Snow_and_Ice","Barren_or_Sparsely_Vegetated","ndvi","evi", "pop", "prec", 
                     "tmin", "tmax", "sp","ws","rh", "elevation", "urban_acc_30s",
                     "mite","rodent_richness")

df_balanced<-readRDS("df_balanced.rds")
# Convert to factor
df_balanced$occurrence <- as.factor(df_balanced$occurrence)
sb<-readRDS("10sb.rds")

# Create formula
formula_rf <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))

library(ranger)
library(parallel)
library(raster)

# Parameters
n_blocks <- 10  
block_size <- 500000

predictors_names <- c("Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
                      "Deciduous_Needleleaf_Forests","Deciduous_Broadleaf_Forests","Mixed_Forests",
                      "Closed_Shrublands","Open_Shrublands","Woody_Savannas",
                      "Savannas","Grasslands","Permanent_Wetlands",
                      "Croplands","Urban_and_Builtup","Cropland_Natural_Vegetation_Mosaics",
                      "Snow_and_Ice","Barren_or_Sparsely_Vegetated","ndvi","evi", "pop", "prec", 
                      "tmin", "tmax", "sp","ws","rh", "elevation", "urban_acc_30s"
                      ,"mite","rodent_richness")

formula_rf <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))

# Create directories for saving models
if (!dir.exists("trained_models")) {
  dir.create("trained_models")
}

cat("=== RF MODEL TRAINING (10 FOLDS) ===\n")
cat("Training weighted and non-weighted RF models for spatial prediction\n")
cat("Models will be saved to: trained_models/ directory\n\n")

# Lists to store trained models
rf_models_weight <- list()

valid_folds <- 0
start_time <- Sys.time()

for (i in 1:n_blocks) {
  fold_start_time <- Sys.time()
  cat(sprintf("[%s] Processing Fold %d/%d", format(Sys.time(), "%H:%M:%S"), i, n_blocks))
  
  tryCatch({
    test_idx <- sb$folds_list[[i]][[2]]
    train_idx <- sb$folds_list[[i]][[1]]
    
    # Check fold size and class balance
    if (length(test_idx) < 5 || length(train_idx) < 5) {
      cat(sprintf(" - SKIPPED (insufficient data: test:%d, train:%d)\n", length(test_idx), length(train_idx)))
      next
    }
    
    train_set <- df_balanced[train_idx, ]
    test_set <- df_balanced[test_idx, ]
    
    # Check if both classes exist in training and test sets
    train_classes <- unique(train_set$occurrence)
    test_classes <- unique(test_set$occurrence)
    
    if (length(train_classes) < 2) {
      cat(" - SKIPPED (training set has only one class)\n")
      next
    }
    
    if (length(test_classes) < 2) {
      cat(" - SKIPPED (test set has only one class)\n")
      next
    }
    
    valid_folds <- valid_folds + 1
    cat(sprintf(" (Valid fold #%d)\n", valid_folds))
    
    train_set$occurrence <- as.factor(train_set$occurrence)
    
    ## ----------- Train Weighted RF -----------
    tryCatch({
      rf_weight <- ranger(formula_rf,
                          data = train_set,
                          num.trees = 900,
                          probability = TRUE,
                          importance = 'impurity',
                          sample.fraction = 1,
                          case.weights = train_set$weight)
      
      # Save model
      model_filename <- sprintf("trained_models/rf_weight_fold%02d.rds", i)
      saveRDS(rf_weight, model_filename)
      rf_models_weight[[valid_folds]] <- rf_weight
      
      cat(sprintf("  ✓ Weighted RF saved: %s\n", model_filename))
      
    }, error = function(e) {
      cat(sprintf("  ✗ Weighted RF failed: %s\n", e$message))
    })
    
    # Fold completion summary
    fold_duration <- as.numeric(difftime(Sys.time(), fold_start_time, units = "mins"))
    cat(sprintf("  Fold %d completed in %.2f minutes\n", i, fold_duration))
    
  }, error = function(e) {
    cat(sprintf("  - FOLD ERROR: %s\n", e$message))
  })
}

training_duration <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("\n=== TRAINING COMPLETED ===\n")
cat(sprintf("Training time: %.2f minutes\n", training_duration))
cat(sprintf("Valid folds: %d/%d\n", valid_folds, n_blocks))
cat(sprintf("Non-weighted models: %d\n", length(rf_models_noweight)))
cat(sprintf("Weighted models: %d\n", length(rf_models_weight)))

## ========================================
## SPATIAL PREDICTION SECTION
## ========================================

cat("\n=== STARTING SPATIAL PREDICTION ===\n")

# Load prediction data
predict_data <- readRDS("predict_scaled_covariates_2020.rds")
template_raster <- raster("annual_prec_2020.tif")

# Prepare data chunks
predict_data$.__index__ <- 1:nrow(predict_data)
chunk_size <- ceiling(nrow(predict_data) / 600)
data_chunks <- list()
for (i in 1:600) {
  start_idx <- (i-1) * chunk_size + 1
  end_idx <- min(i * chunk_size, nrow(predict_data))
  data_chunks[[i]] <- predict_data[start_idx:end_idx, ]
}

cat(sprintf("Prediction data: %d rows\n", nrow(predict_data)))
cat(sprintf("Data chunks: %d chunks of ~%d rows each\n", length(data_chunks), chunk_size))

# Prediction function for ensemble
predict_ensemble <- function(chunk, models_list, model_type) {
  progress_msg <- paste("Processing", model_type, "chunk with indices from", min(chunk$.__index__), 
                        "to", max(chunk$.__index__), "at", Sys.time())
  cat(progress_msg, "\n", file = paste0("progress_", model_type, ".txt"), append = TRUE)
  
  na_rows <- apply(chunk, 1, function(x) any(is.na(x)))
  
  n_models <- length(models_list)
  predictions_matrix <- matrix(NA, nrow = nrow(chunk), ncol = n_models)
  
  # Get predictions from all models for valid rows
  if (any(!na_rows)) {
    for (j in 1:n_models) {
      tryCatch({
        pred_result <- predict(models_list[[j]], chunk[!na_rows,], type = "response")
        predictions_matrix[!na_rows, j] <- pred_result$predictions[,2]  # probability of class 1
      }, error = function(e) {
        cat(sprintf("Error in model %d: %s\n", j, e$message))
      })
    }
  }
  
  # Calculate ensemble mean and uncertainty
  ensemble_prob <- rowMeans(predictions_matrix, na.rm = TRUE)
  
  # Calculate ensemble standard deviation (model uncertainty)
  ensemble_sd <- apply(predictions_matrix, 1, sd, na.rm = TRUE)
  
  completion_msg <- paste("Completed", model_type, "chunk", min(chunk$.__index__), "-", 
                          max(chunk$.__index__), "at", Sys.time())
  cat(completion_msg, "\n", file = paste0("progress_", model_type, ".txt"), append = TRUE)
  
  gc()
  return(data.frame(
    index = chunk$.__index__, 
    prob = ensemble_prob, 
    sd = ensemble_sd
  ))
}

## ----------- Predict with Weighted ensemble -----------
cat("\n--- Weighted RF Ensemble Prediction ---\n")

# Clear progress file
progress_file <- "progress_weight.txt"
if (file.exists(progress_file)) file.remove(progress_file)

cl <- makeCluster(n_cores)

# Export models and function to cluster
clusterExport(cl, varlist = c("rf_models_weight", "predict_ensemble"))
clusterEvalQ(cl, library(ranger))

# Run prediction
result_weight <- parLapply(cl, data_chunks, function(chunk) {
  predict_ensemble(chunk, rf_models_weight, "weight")
})

stopCluster(cl)
gc()

# Combine results
RF_predictions_weight <- unlist(lapply(result_weight, function(x) x$prob))
RF_sd_weight <- unlist(lapply(result_weight, function(x) x$sd))

# Create rasters for weighted ensemble
RF_raster_weight <- template_raster
values(RF_raster_weight) <- RF_predictions_weight

sd_raster_weight <- template_raster
values(sd_raster_weight) <- RF_sd_weight

# Save weighted results
writeRaster(RF_raster_weight, filename = "RF_pred_weight_ensemble.tif", format = "GTiff", overwrite = TRUE)
writeRaster(sd_raster_weight, filename = "RF_pred_weight_ensemble_sd.tif", format = "GTiff", overwrite = TRUE)

cat("✓ Weighted ensemble predictions saved\n")

## ----------- Summary -----------
end_time <- Sys.time()
total_duration <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("\n=== PROCESS COMPLETED ===\n")
cat(sprintf("Total runtime: %.2f minutes (%.2f hours)\n", total_duration, total_duration/60))
cat(sprintf("Models trained: %d non-weighted + %d weighted\n", length(rf_models_noweight), length(rf_models_weight)))

cat("\nFiles created:\n")
cat("Training models:\n")
cat("- trained_models/rf_weight_fold*.rds\n")
cat("\nPrediction rasters:\n")
cat("- RF_pred_weight_ensemble.tif (ensemble mean probability)\n")
cat("- RF_pred_weight_ensemble_se.tif (ensemble mean standard error)\n")
cat("- RF_pred_weight_ensemble_sd.tif (ensemble standard deviation)\n")
cat("\nProgress logs:\n")
cat("- progress_weight.txt\n")

# Save model summary
model_summary <- data.frame(
  Model_Type = c(rep("Weighted", length(rf_models_weight))),
  Fold = c( 1:length(rf_models_weight)),
  Trees = c(rep(900,  length(rf_models_weight))),
  Training_Complete = TRUE
)

write.csv(model_summary, "model_training_summary.csv", row.names = FALSE)
cat("- model_training_summary.csv\n")

