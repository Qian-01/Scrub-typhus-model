###spaital block cross validation and cutoff----

library(blockCV)
library(ranger)
library(pROC)
library(sf)
library(caret)

df_balanced<-readRDS("df_balanced.rds")
df_sf <- st_as_sf(df_balanced, coords = c("coordinates_long", "coordinates_lat"), crs = 4326)
trainSet_df_sf <- st_transform(df_sf , crs = 3857)

sb <- cv_spatial(x = trainSet_df_sf,
                   column  = "occurrence",
                   k = 10,
                   size = 500000,
                   selection = "random",
                   iteration = 100)
plot(sb)
saveRDS(sb,"10sb.rds")

sb<-readRDS("10sb.rds")

n_blocks <- 10  
block_size <- 500000
metrics_no <- data.frame(
  AUC = numeric(n_blocks),
  Sensitivity = numeric(n_blocks),
  Specificity = numeric(n_blocks),
  PPV = numeric(n_blocks),  # Positive Predictive Value
  NPV = numeric(n_blocks),  # Negative Predictive Value
  Cutoff = numeric(n_blocks),
  Youden = numeric(n_blocks),
  Accuracy = numeric(n_blocks),
  F1_Score = numeric(n_blocks),
  Kappa = numeric(n_blocks)
)

metrics_w <- data.frame(
  AUC = numeric(n_blocks),
  Sensitivity = numeric(n_blocks),
  Specificity = numeric(n_blocks),
  PPV = numeric(n_blocks),
  NPV = numeric(n_blocks),
  Cutoff = numeric(n_blocks),
  Youden = numeric(n_blocks),
  Accuracy = numeric(n_blocks),
  F1_Score = numeric(n_blocks),
  Kappa = numeric(n_blocks)
)
predictors_names<- c("Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
                     "Deciduous_Needleleaf_Forests","Deciduous_Broadleaf_Forests","Mixed_Forests",
                     "Closed_Shrublands","Open_Shrublands","Woody_Savannas",
                     "Savannas","Grasslands","Permanent_Wetlands",
                     "Croplands","Urban_and_Builtup","Cropland_Natural_Vegetation_Mosaics",
                     "Snow_and_Ice","Barren_or_Sparsely_Vegetated","ndvi","evi", "pop", "prec", 
                     "tmin", "tmax", "sp","ws","rh", "elevation", "urban_acc_30s"
                     ,"mite","rodent_richness")
formula_rf <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))

valid_folds <- 0
for (i in 1:n_blocks) {
  tryCatch({
    test_idx <- sb$folds_list[[i]][[2]]
    train_idx <- sb$folds_list[[i]][[1]]
    
    if (length(test_idx) < 5 || length(train_idx) < 5) {
      cat("Skipping fold", i, "- insufficient data (test:", length(test_idx), "train:", length(train_idx), ")\n")
      next
    }
    
    train_set <- df_balanced[train_idx, ]
    test_set <- df_balanced[test_idx, ]
    
    train_classes <- unique(train_set$occurrence)
    test_classes <- unique(test_set$occurrence)
    
    if (length(train_classes) < 2) {
      cat("Skipping fold", i, "- training set has only one class\n")
      next
    }
    
    if (length(test_classes) < 2) {
      cat("Skipping fold", i, "- test set has only one class\n")
      next
    }
    
    valid_folds <- valid_folds + 1
    
    train_set$occurrence<-as.factor(train_set$occurrence)
    
    ## ----------- Non-weighted model
    
    rf_cv <- ranger(formula_rf,
                    data = train_set,
                    num.trees = 900,
                    probability = TRUE,
                    importance = 'impurity')
    
    pred <- predict(rf_cv, data = test_set, type = "response")
    probs <- pred$predictions[, 2]
    
    actual_value<- test_set$occurrence
    
    if (length(probs) != length(actual_value)) {
      cat("Skipping fold", i, "- prediction and actual_value length mismatch\n")
      next
    }
    
    # ROC analysis
    roc_obj <- roc(actual_value, probs, quiet = TRUE)
    metrics_no$AUC[valid_folds] <- auc(roc_obj)
    
    # Find optimal cutoff using Youden's Index
    coords_result <- coords(roc_obj, "best", best.method = "youden")
    cutoff <- coords_result$threshold
    sensitivity <- coords_result$sensitivity
    specificity <- coords_result$specificity
    
    metrics_no$Cutoff[valid_folds] <- cutoff
    metrics_no$Sensitivity[valid_folds] <- sensitivity
    metrics_no$Specificity[valid_folds] <- specificity
    metrics_no$Youden[valid_folds] <- sensitivity + specificity - 1
    
    # Binary predictions using optimal cutoff
    binary_pred <- ifelse(probs >= cutoff, 1, 0)
    
    
    #if (length(binary_pred) != length(actual_value)) {
     # cat("Warning: binary_pred and actual_value length mismatch in fold", i, "\n")
      #binary_pred <- binary_pred[1:length(actual_value)]
    #}
    
    # Confusion matrix with error handling
    cm <- confusionMatrix(factor(binary_pred, levels = c(0, 1)), 
                          factor(actual_value, levels = c(0, 1)), 
                          positive = "1")
    
    # Extract additional metrics with NA handling
    metrics_no$PPV[valid_folds] <- ifelse(is.na(cm$byClass["Pos Pred Value"]), 0, cm$byClass["Pos Pred Value"])
    metrics_no$NPV[valid_folds] <- ifelse(is.na(cm$byClass["Neg Pred Value"]), 0, cm$byClass["Neg Pred Value"])
    metrics_no$Accuracy[valid_folds] <- cm$overall["Accuracy"]
    metrics_no$F1_Score[valid_folds] <- ifelse(is.na(cm$byClass["F1"]), 0, cm$byClass["F1"])
    metrics_no$Kappa[valid_folds] <- cm$overall["Kappa"]
    
    ## ----------- Weighted model
    rfw_cv <- ranger(formula_rf,
                     data = train_set,
                     num.trees = 900,
                     probability = TRUE,
                     importance = 'impurity',
                     sample.fraction = 1,
                     case.weights = train_set$weight)
    
    predw <- predict(rfw_cv, data = test_set, type = "response")
    probsw <- predw$predictions[, 2]

    # ROC analysis for weighted model
    roc_objw <- roc(actual_value, probsw, quiet = TRUE)
    metrics_w$AUC[valid_folds] <- auc(roc_objw)
    
    # Find optimal cutoff using Youden's Index
    coords_resultw <- coords(roc_objw, "best", best.method = "youden")
    cutoffw <- coords_resultw$threshold
    sensitivityw <- coords_resultw$sensitivity
    specificityw <- coords_resultw$specificity
    
    metrics_w$Cutoff[valid_folds] <- cutoffw
    metrics_w$Sensitivity[valid_folds] <- sensitivityw
    metrics_w$Specificity[valid_folds] <- specificityw
    metrics_w$Youden[valid_folds] <- sensitivityw + specificityw - 1
    
    # Binary predictions using optimal cutoff
    binary_predw <- ifelse(probsw >= cutoffw, 1, 0)

    #if (length(binary_predw) != length(actual_value)) {
    #  binary_predw <- binary_predw[1:length(actual_value)]
    #}
    
    # Confusion matrix for weighted model
    cmw <- confusionMatrix(factor(binary_predw, levels = c(0, 1)), 
                           factor(actual_value, levels = c(0, 1)), 
                           positive = "1")
    
    # Extract additional metrics with NA handling
    metrics_w$PPV[valid_folds] <- ifelse(is.na(cmw$byClass["Pos Pred Value"]), 0, cmw$byClass["Pos Pred Value"])
    metrics_w$NPV[valid_folds] <- ifelse(is.na(cmw$byClass["Neg Pred Value"]), 0, cmw$byClass["Neg Pred Value"])
    metrics_w$Accuracy[valid_folds] <- cmw$overall["Accuracy"]
    metrics_w$F1_Score[valid_folds] <- ifelse(is.na(cmw$byClass["F1"]), 0, cmw$byClass["F1"])
    metrics_w$Kappa[valid_folds] <- cmw$overall["Kappa"]
    
    if (valid_folds %% 10 == 0) {
      cat("Completed", valid_folds, "valid folds out of", i, "total folds\n")
    }
    
  }, error = function(e) {
    cat("Error in fold", i, ":", e$message, "\n")
  })
}

metrics_no <- metrics_no[complete.cases(metrics_no), ]
str(metrics_w)

list_columns <- c("Cutoff", "Sensitivity", "Specificity", "Youden")

for (col in list_columns) {
  metrics_w[[col]] <- unlist(metrics_w[[col]])
}

metrics_w <- metrics_w[complete.cases(metrics_w), ]

calculate_summary <- function(metrics, model_name) {
  summary_stats <- data.frame(
    Metric = names(metrics),
    Mean = sapply(metrics, mean, na.rm = TRUE),
    SD = sapply(metrics, sd, na.rm = TRUE),
    Median = sapply(metrics, median, na.rm = TRUE),
    Q25 = sapply(metrics, quantile, probs = 0.25, na.rm = TRUE),
    Q75 = sapply(metrics, quantile, probs = 0.75, na.rm = TRUE),
    CI_Lower = sapply(metrics, quantile, probs = 0.025, na.rm = TRUE),
    CI_Upper = sapply(metrics, quantile, probs = 0.975, na.rm = TRUE)
  )
  summary_stats$Model <- model_name
  return(summary_stats)
}

summary_no <- calculate_summary(metrics_no, "Non-weighted")
summary_w <- calculate_summary(metrics_w, "Weighted")
combined_summary <- rbind(summary_no, summary_w)

cat("=== MODEL PERFORMANCE SUMMARY ===\n")
cat("Number of valid folds:", nrow(metrics_no), "\n")
cat("Block size:", block_size/1000, "km\n\n")

print(combined_summary)

write.csv(combined_summary, "10spatial_cv_comprehensive_summary-weight1.csv", row.names = FALSE)

detailed_results <- data.frame(
  Fold = rep(1:nrow(metrics_no), 2),
  Model = rep(c("Non-weighted", "Weighted"), each = nrow(metrics_no)),
  rbind(metrics_no, metrics_w)
)

write.csv(detailed_results, "10spatial_cv_detailed_results-weight1.csv", row.names = FALSE)

