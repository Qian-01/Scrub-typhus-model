
#load package----
library(readxl)
library(data.table)
library(dplyr)
df_balanced<-readRDS("df_balanced.rds")

# -------------------- Collinearity Analysis --------------------
cat("Performing collinearity analysis...\n")
library(car)
library(corrplot)
library(tidyverse)

# Extract feature data
predictors_names<- c("Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
                                     "Deciduous_Needleleaf_Forests","Deciduous_Broadleaf_Forests","Mixed_Forests",
                                     "Closed_Shrublands","Open_Shrublands","Woody_Savannas",
                                     "Savannas","Grasslands","Permanent_Wetlands",
                                     "Croplands","Urban_and_Builtup","Cropland_Natural_Vegetation_Mosaics",
                                     "Snow_and_Ice","Barren_or_Sparsely_Vegetated","ndvi","evi", "pop", "prec", 
                                     "tmin", "tmax", "sp","ws","rh", "elevation", "urban_acc_30s"
                                     ,"mite","rodent_richness")
feature_data<-df_balanced[,predictors_names]

# parson correlation coefficient
library(pheatmap)
cor_matrix <- cor(feature_data, use="pairwise.complete.obs") 
pheatmap(cor_matrix, 
         display_numbers = TRUE, 
         angle_col = 45,    
         angle_row = 45,   
         fontsize = 5,
         filename="correlation_heatmap.png")      

# Set correlation threshold
cor_threshold <- 0.8

# Identify highly correlated pairs
high_cor_pairs <- which(abs(cor_matrix) > cor_threshold & upper.tri(cor_matrix), arr.ind = TRUE)
if (nrow(high_cor_pairs) > 0) {
  cat("Highly correlated variable pairs found:\n")
  high_cor_vars <- data.frame(
    Var1 = rownames(cor_matrix)[high_cor_pairs[, 1]],
    Var2 = colnames(cor_matrix)[high_cor_pairs[, 2]],
    Correlation = cor_matrix[high_cor_pairs]
  )
  print(high_cor_vars)
  
  var_variance <- apply(feature_data, 2, var, na.rm = TRUE)
  vars_to_remove <- c()
  
  for (i in 1:nrow(high_cor_pairs)) {
    var1 <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
    var2 <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
    
    if (var_variance[var1] > var_variance[var2]) {
      vars_to_remove <- c(vars_to_remove, var2)
    } else {
      vars_to_remove <- c(vars_to_remove, var1)
    }
  }
  
  vars_to_remove <- unique(vars_to_remove)
  feature_data <- feature_data[, !(names(feature_data) %in% vars_to_remove)]
  
  cat("Variables removed due to high collinearity:", paste(vars_to_remove, collapse = ", "), "\n")
} else {
  cat("No highly correlated variable pairs found, retaining all variables.\n")
}

cat("Final retained feature variables:", paste(colnames(feature_data), collapse = ", "), "\n")
#Final retained feature variables: Water_bodies, Evergreen_Needleleaf_Forests, Evergreen_Broadleaf_Forests, 
#Deciduous_Needleleaf_Forests, Deciduous_Broadleaf_Forests, Mixed_Forests, Closed_Shrublands, Open_Shrublands, 
#Woody_Savannas, Grasslands, Permanent_Wetlands, Croplands, Urban_and_Builtup, Cropland_Natural_Vegetation_Mosaics, 
#Snow_and_Ice, Barren_or_Sparsely_Vegetated, evi, pop, prec, tmax, sp, ws, rh, urban_acc_30s, mite, rodent_richness

# Save collinearity results
write.csv(cor_matrix, "Correlation_Matrix.csv", row.names = TRUE)
cat("Collinearity analysis completed, results saved ")

#concurvity test for GAM---
library(mgcv)
fit_model <- function(var_to_remove, formula, data) {
  reduced_vars <- setdiff(predictors_names, var_to_remove)
  reduced_formula <- as.formula(paste("occurrence ~", paste("s(", reduced_vars, ", k = 10)", collapse = " + ")))
  model <- gam(reduced_formula, family = binomial(link = logit), data = data)
  return(AIC(model))
}
aic_scores <- list()
high_concurvity_pairs <- list(c("ndvi", "evi"), c("tmax", "tmin"), c("Barren_or_Sparsely_Vegetated", "Savannas"),c("elevation","sp"))
for(pair in high_concurvity_pairs) {
  for(var in pair) {
    aic_scores[[paste("Remove", var)]] <- fit_model(var, original_formula, df_balanced)
  }
}
aic_scores
library(mgcv)
library(dplyr)
gam_formula0 <- as.formula(paste("occurrence ~", paste("s(", predictors_names, ", k = 10)", collapse = " + ")))
gam.fit0<- gam(gam_formula0,family= binomial(link = logit), data=df_balanced)
summary(gam.fit0)
print(paste(" AIC :", gam.fit0$aic, "| deviance :", gam.fit0$deviance))
conc0 <-concurvity(gam.fit0,full=FALSE)
saveRDS(gam.fit0,"gam-full.rds")

###visualization of con-curvity of covariates
conc_matrix0 <- conc0$worst
conc_long0 <- as.data.frame(as.table(conc_matrix0))
names(conc_long0) <- c("Term1", "Term2", "Concurvity")
GAMcon0<-ggplot(conc_long0, aes(x=Term1, y=Term2, fill=Concurvity)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="red") + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  labs(fill="Concurvity", 
       title="Concurvity Matrix Visualization",
       x="", 
       y="") 
conc_long0[conc_long0$Concurvity > 0.8 & conc_long0$Term1 != conc_long0$Term2, ]
ggsave("gam_concurvity.png",GAMcon0, height=20, width = 22,units="cm")


###spaital block cross validation for GLM,GAM,BRT,RF----

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

# Initialize metrics data frames for all 8 models
model_names <- c("GLM_no", "GLM_w", "GAM_no", "GAM_w", "BRT_no", "BRT_w", "RF_no", "RF_w")

# Create a list to store metrics for all models
metrics_list <- list()
for (model in model_names) {
  metrics_list[[model]] <- data.frame(
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
}

predictors_names <- c("Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
                      "Deciduous_Needleleaf_Forests","Deciduous_Broadleaf_Forests","Mixed_Forests",
                      "Closed_Shrublands","Open_Shrublands","Woody_Savannas",
                      "Grasslands","Permanent_Wetlands",
                      "Croplands","Urban_and_Builtup","Cropland_Natural_Vegetation_Mosaics",
                      "Snow_and_Ice","Barren_or_Sparsely_Vegetated","evi", "pop", "prec", 
                       "tmax", "sp","ws","rh", "urban_acc_30s"
                      ,"mite","rodent_richness")
# Create formulas for different models
formula_glm <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))
formula_brt <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))
formula_rf <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))

predictors_names <- c("Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
                      "Deciduous_Needleleaf_Forests","Deciduous_Broadleaf_Forests","Mixed_Forests",
                      "Closed_Shrublands","Open_Shrublands","Woody_Savannas",
                      "Grasslands","Permanent_Wetlands",
                      "Croplands","Cropland_Natural_Vegetation_Mosaics",
                      "Snow_and_Ice","Barren_or_Sparsely_Vegetated","evi", "pop", "prec", 
                      "tmin", "sp","ws","rh",  "urban_acc_30s"
                      ,"mite","rodent_richness")
formula_gam <- as.formula(paste("occurrence ~", paste(paste0("s(", predictors_names, ")"), collapse = " + ")))

# Load required libraries
library(mgcv)      # for GAM
library(gbm)       # for BRT
library(ranger)    # for RF
library(pROC)      # for ROC analysis
library(caret)     # for confusionMatrix

valid_folds <- 0

for (i in 1:n_blocks) {
  tryCatch({
    test_idx <- sb$folds_list[[i]][[2]]
    train_idx <- sb$folds_list[[i]][[1]]
    
    # Check fold size and class balance
    if (length(test_idx) < 5 || length(train_idx) < 5) {
      cat("Skipping fold", i, "- insufficient data (test:", length(test_idx), "train:", length(train_idx), ")\n")
      next
    }
    
    train_set <- df_balanced[train_idx, ]
    test_set <- df_balanced[test_idx, ]
    
    # Check if both classes exist in training and test sets
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
    
    actual_value <- test_set$occurrence
    
    # Function to calculate metrics for a model
    calculate_metrics <- function(probs, model_name, fold_num) {
      # ROC analysis
      roc_obj <- roc(actual_value, probs, quiet = TRUE)
      auc_val <- auc(roc_obj)
      
      # Find optimal cutoff using Youden's Index
      coords_result <- coords(roc_obj, "best", best.method = "youden")
      cutoff <- coords_result$threshold
      sensitivity <- coords_result$sensitivity
      specificity <- coords_result$specificity
      
      # Binary predictions using optimal cutoff
      binary_pred <- ifelse(probs >= cutoff, 1, 0)
      
      # Confusion matrix
      cm <- confusionMatrix(factor(binary_pred, levels = c(0, 1)), 
                            factor(actual_value, levels = c(0, 1)), 
                            positive = "1")
      
      # Store metrics
      metrics_list[[model_name]]$AUC[fold_num] <<- auc_val
      metrics_list[[model_name]]$Cutoff[fold_num] <<- cutoff
      metrics_list[[model_name]]$Sensitivity[fold_num] <<- sensitivity
      metrics_list[[model_name]]$Specificity[fold_num] <<- specificity
      metrics_list[[model_name]]$Youden[fold_num] <<- sensitivity + specificity - 1
      metrics_list[[model_name]]$PPV[fold_num] <<- ifelse(is.na(cm$byClass["Pos Pred Value"]), 0, cm$byClass["Pos Pred Value"])
      metrics_list[[model_name]]$NPV[fold_num] <<- ifelse(is.na(cm$byClass["Neg Pred Value"]), 0, cm$byClass["Neg Pred Value"])
      metrics_list[[model_name]]$Accuracy[fold_num] <<- cm$overall["Accuracy"]
      metrics_list[[model_name]]$F1_Score[fold_num] <<- ifelse(is.na(cm$byClass["F1"]), 0, cm$byClass["F1"])
      metrics_list[[model_name]]$Kappa[fold_num] <<- cm$overall["Kappa"]
    }
    
    ## ----------- Model 1: GLM (Non-weighted) -----------
    tryCatch({
      glm_model <- glm(formula_glm, data = train_set, family = binomial)
      glm_pred <- predict(glm_model, newdata = test_set, type = "response")
      calculate_metrics(glm_pred, "GLM_no", valid_folds)
    }, error = function(e) {
      cat("Error in GLM non-weighted for fold", i, ":", e$message, "\n")
    })
    
    ## ----------- Model 2: GLM (Weighted) -----------
    tryCatch({
      glm_model_w <- glm(formula_glm, data = train_set, family = binomial, weights = train_set$weight)
      glm_pred_w <- predict(glm_model_w, newdata = test_set, type = "response")
      calculate_metrics(glm_pred_w, "GLM_w", valid_folds)
    }, error = function(e) {
      cat("Error in GLM weighted for fold", i, ":", e$message, "\n")
    })
    
    ## ----------- Model 3: GAM (Non-weighted) -----------
    tryCatch({
      gam_model <- gam(formula_gam, data = train_set, family = binomial)
      gam_pred <- predict(gam_model, newdata = test_set, type = "response")
      calculate_metrics(gam_pred, "GAM_no", valid_folds)
    }, error = function(e) {
      cat("Error in GAM non-weighted for fold", i, ":", e$message, "\n")
    })
    
    ## ----------- Model 4: GAM (Weighted) -----------
    tryCatch({
      gam_model_w <- gam(formula_gam, data = train_set, family = binomial, weights = train_set$weight)
      gam_pred_w <- predict(gam_model_w, newdata = test_set, type = "response")
      calculate_metrics(gam_pred_w, "GAM_w", valid_folds)
    }, error = function(e) {
      cat("Error in GAM weighted for fold", i, ":", e$message, "\n")
    })
    
    ## ----------- Model 5: BRT (Non-weighted) -----------
    tryCatch({
      # Convert factor to numeric for gbm
      train_set_brt <- train_set
      train_set_brt$occurrence <- as.numeric(as.character(train_set_brt$occurrence))
      
      brt_model <- gbm(formula_brt, data = train_set_brt, 
                       distribution = "bernoulli",
                       n.trees = 10000,
                       shrinkage = 0.005,
                       interaction.depth = 4,
                       bag.fraction = 0.75,
                       n.minobsinnode = 10,
                       cv.folds = 0,
                       verbose = FALSE)
      
      # Find optimal number of trees
      best_iter <- gbm.perf(brt_model, method = "OOB", plot.it = FALSE)
      brt_pred <- predict(brt_model, newdata = test_set, n.trees = best_iter, type = "response")
      calculate_metrics(brt_pred, "BRT_no", valid_folds)
    }, error = function(e) {
      cat("Error in BRT non-weighted for fold", i, ":", e$message, "\n")
    })
    
    ## ----------- Model 6: BRT (Weighted) -----------
    tryCatch({
      # Convert factor to numeric for gbm
      train_set_brt_w <- train_set
      train_set_brt_w$occurrence <- as.numeric(as.character(train_set_brt_w$occurrence))
      
      brt_model_w <- gbm(formula_brt, data = train_set_brt_w, 
                         distribution = "bernoulli",
                         n.trees = 10000,
                         shrinkage = 0.005,
                         interaction.depth = 4,
                         bag.fraction = 0.75,
                         n.minobsinnode = 10,
                         weights = train_set_brt_w$weight,
                         cv.folds = 0,
                         verbose = FALSE)
      
      # Find optimal number of trees
      best_iter_w <- gbm.perf(brt_model_w, method = "OOB", plot.it = FALSE)
      brt_pred_w <- predict(brt_model_w, newdata = test_set, n.trees = best_iter_w, type = "response")
      calculate_metrics(brt_pred_w, "BRT_w", valid_folds)
    }, error = function(e) {
      cat("Error in BRT weighted for fold", i, ":", e$message, "\n")
    })
    
    ## ----------- Model 7: RF (Non-weighted) -----------
    tryCatch({
      train_set$occurrence <- as.factor(train_set$occurrence)
      rf_model <- ranger(formula_rf,
                         data = train_set,
                         num.trees = 900,
                         probability = TRUE,
                         importance = 'impurity')
      
      rf_pred <- predict(rf_model, data = test_set, type = "response")
      rf_probs <- rf_pred$predictions[, 2]
      calculate_metrics(rf_probs, "RF_no", valid_folds)
    }, error = function(e) {
      cat("Error in RF non-weighted for fold", i, ":", e$message, "\n")
    })
    
    ## ----------- Model 8: RF (Weighted) -----------
    tryCatch({
      train_set$occurrence <- as.factor(train_set$occurrence)
      rf_model_w <- ranger(formula_rf,
                           data = train_set,
                           num.trees = 900,
                           probability = TRUE,
                           importance = 'impurity',
                           sample.fraction = 1,
                           case.weights = train_set$weight)
      
      rf_pred_w <- predict(rf_model_w, data = test_set, type = "response")
      rf_probs_w <- rf_pred_w$predictions[, 2]
      calculate_metrics(rf_probs_w, "RF_w", valid_folds)
    }, error = function(e) {
      cat("Error in RF weighted for fold", i, ":", e$message, "\n")
    })
    
    if (valid_folds %% 5 == 0) {
      cat("Completed", valid_folds, "valid folds out of", i, "total folds\n")
    }
    
  }, error = function(e) {
    cat("Error in fold", i, ":", e$message, "\n")
  })
}

# Clean up metrics by removing incomplete cases
for (model in model_names) {
  # Check for list columns and convert to numeric if needed
  for (col in names(metrics_list[[model]])) {
    if (is.list(metrics_list[[model]][[col]])) {
      metrics_list[[model]][[col]] <- unlist(metrics_list[[model]][[col]])
    }
  }
  
  # Remove incomplete cases
  metrics_list[[model]] <- metrics_list[[model]][complete.cases(metrics_list[[model]]), ]
}

# Calculate summary statistics
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

# Generate comprehensive summary
all_summaries <- list()
for (model in model_names) {
  all_summaries[[model]] <- calculate_summary(metrics_list[[model]], model)
}

combined_summary <- do.call(rbind, all_summaries)

# Print results
cat("=== MODEL PERFORMANCE SUMMARY ===\n")
cat("Number of valid folds:", valid_folds, "\n")
cat("Block size:", block_size/1000, "km\n\n")

print(combined_summary)

# Save comprehensive summary
write.csv(combined_summary, "8models_spatial_cv_comprehensive_summary.csv", row.names = FALSE)

# Create detailed results data frame
detailed_results <- data.frame()
for (model in model_names) {
  model_data <- data.frame(
    Fold = 1:nrow(metrics_list[[model]]),
    Model = model,
    metrics_list[[model]]
  )
  detailed_results <- rbind(detailed_results, model_data)
}

write.csv(detailed_results, "8models_spatial_cv_detailed_results.csv", row.names = FALSE)

# Create a comparison table showing mean AUC for each model
auc_comparison <- data.frame(
  Model = model_names,
  Mean_AUC = sapply(model_names, function(x) mean(metrics_list[[x]]$AUC, na.rm = TRUE)),
  SD_AUC = sapply(model_names, function(x) sd(metrics_list[[x]]$AUC, na.rm = TRUE)),
  Mean_Accuracy = sapply(model_names, function(x) mean(metrics_list[[x]]$Accuracy, na.rm = TRUE)),
  Mean_F1 = sapply(model_names, function(x) mean(metrics_list[[x]]$F1_Score, na.rm = TRUE))
)

cat("\n=== MODEL COMPARISON (AUC, Accuracy, F1-Score) ===\n")
print(auc_comparison)

write.csv(auc_comparison, "8models_comparison_summary.csv", row.names = FALSE)

cat("\nAnalysis complete! Files saved:\n")
cat("- 8models_spatial_cv_comprehensive_summary.csv\n")
cat("- 8models_spatial_cv_detailed_results.csv\n")
cat("- 8models_comparison_summary.csv\n")

# visualiztion
library(ggplot2)
library(reshape2)

plot_data <- melt(detailed_results, 
                  id.vars = c("Fold", "Model"),
                  variable.name = "Metric",
                  value.name = "Value")

main_metrics <- c("AUC", "Sensitivity", "Specificity", "PPV", "NPV", "F1_Score")
plot_data_main <- plot_data[plot_data$Metric %in% main_metrics, ]

p <- ggplot(plot_data_main, aes(x = Model, y = Value, fill = Model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  theme_minimal() +
  labs(title = paste("Spatial Block CV Performance Metrics (", nrow(metrics_no), "folds)"),
       subtitle = paste("Block size:", block_size/1000, "km")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("spatial_cv_metrics_comparison-weight1.png", p, width = 8, height = 5, dpi = 300)

cat("\n=== STATISTICAL COMPARISON ===\n")
for (metric in main_metrics) {
  test_result <- wilcox.test(metrics_no[[metric]], metrics_w[[metric]], paired = TRUE)
  cat(metric, "- p-value:", round(test_result$p.value, 4), "\n")
}

cat("\n=== RECOMMENDATIONS FOR GLOBAL DATA ===\n")
cat("For global species distribution modeling, consider:\n")
cat("1. Block number: 100-200 blocks (current:", n_blocks, ")\n")
cat("2. Block size: 25-100 km depending on species mobility\n")
cat("3. Consider environmental stratification for block selection\n")
cat("4. Account for sampling bias and geographic clustering\n")
cat("5. Use buffer zones around protected areas if relevant\n")



