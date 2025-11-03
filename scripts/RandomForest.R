#load library

library(raster)
library(dplyr)
library(sf)
library(data.table)
library(ranger)
library(parallel)
library(doParallel)
library(foreach)
library(raster)
library(dplyr)
library(sf)
library(data.table)
library(ranger)
library(parallel)
library(doParallel)
library(foreach)

df_balanced <- readrds("df_balanced.rds")
df_balanced$occurrence <- as.factor(df_balanced$occurrence)

# Create formula
predictors_names<- c("Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
                     "Deciduous_Needleleaf_Forests","Deciduous_Broadleaf_Forests","Mixed_Forests",
                     "Closed_Shrublands","Open_Shrublands","Woody_Savannas",
                     "Savannas","Grasslands","Permanent_Wetlands",
                     "Croplands","Urban_and_Builtup","Cropland_Natural_Vegetation_Mosaics",
                     "Snow_and_Ice","Barren_or_Sparsely_Vegetated","ndvi","evi", "pop", "prec", 
                     "tmin", "tmax", "sp","ws","rh", "elevation", "urban_acc_30s",
                     "mite","rodent_richness")
formula_rf <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))

# Train random forest model
rf.fit <- ranger(formula_rf,
                 data = df_balanced,
                 num.trees = 900,
                 keep.inbag = TRUE,
                 probability = TRUE,
                 importance = 'impurity',
                 case.weights = df_balanced$weight)


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
  #library(ranger)
  #rf.fit <- readRDS("rf_weight.rds") 
  # Write progress to file
  progress_msg <- paste("Processing chunk with indices from", min(chunk$.__index__), 
                        "to", max(chunk$.__index__), "at", Sys.time())
  cat(progress_msg, "\n", file = "progress.txt", append = TRUE)
  
  na_rows <- apply(chunk, 1, function(x) any(is.na(x)))
  
  predictions <- rep(NA, nrow(chunk))
  se.fit <- rep(NA, nrow(chunk))
  
  if (any(!na_rows)) {
    valid_predictions <- predict(rf.fit, chunk[!na_rows,], type = "se")
    predictions[!na_rows] <- valid_predictions[[1]][,2]
    se.fit[!na_rows] <- valid_predictions[[6]][,2]
  }
  # Write completion to file
  completion_msg <- paste("Completed chunk", min(chunk$.__index__), "-", 
                          max(chunk$.__index__), "at", Sys.time())
  cat(completion_msg, "\n", file = "progress.txt", append = TRUE)
  gc()  
  return(data.frame(index = chunk$.__index__, prob = predictions, se = se.fit))
}

clusterExport(cl, varlist = c("rf.fit", "predict_chunk"))
clusterEvalQ(cl, library(ranger)) 
result <- parLapply(cl, data_chunks, predict_chunk)

stopCluster(cl)
gc()

RF_predictions <- unlist(lapply(result, function(x) x$prob))

template_raster <- raster("annual_prec_2020.tif")

RF_raster <- template_raster
values(RF_raster) <- RF_predictions

writeRaster(RF_raster, filename = "RF_pred.tif", format = "GTiff", overwrite = TRUE)
