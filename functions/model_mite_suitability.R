# -------------------- 0. Load Packages --------------------
library(sf)
library(raster)
library(dplyr)
library(blockCV)
library(gbm)
library(dismo)
library(pROC)
library(caret)
library(gridExtra)
library(ggplot2)
library(writexl)
library(corrplot)
library(car)

# Clear memory
cat("Clearing memory...\n")
rm(list = ls())
gc()
cat("Memory cleared.\n")

# Set Times New Roman font
windowsFonts(Times = windowsFont("Times New Roman"))

# -------------------- 1. Paths & Parameters --------------------
csv_path <- "./data/env_data/mite_env.csv"
raster_path <- "./data/raster_data/world.tif"
weight_path <- "./data/weight_results/mite_weight_all.csv"
output_dir <- "./outputs/model_results/"

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "PDP_Plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "ROC_Plots"), showWarnings = FALSE)

num_iterations <- 5
num_folds <- 10
seeds <- 1:num_iterations
train_prop <- 0.7  # Proportion for training set

# -------------------- 2. Data Reading and Preparation --------------------
df <- read.csv(csv_path)
df$cases <- as.numeric(as.character(df$cases))
df_weight <- read.csv(weight_path) %>% dplyr::select(XX, YY, weight)

# Merge weights
df <- left_join(df, df_weight, by = c("XX", "YY"))
df$weight[is.na(df$weight)] <- 1

feature_cols <- setdiff(names(df), c("cases", "XX", "YY", "weight"))
shape_r <- raster(raster_path)

# -------------------- 2.5. Collinearity Analysis --------------------
cat("Performing collinearity analysis...\n")

# Extract feature data
feature_data <- df[, feature_cols]

# Compute Pearson correlation matrix
cor_matrix <- cor(feature_data, method = "pearson", use = "complete.obs")

# Visualize correlation matrix
png(file.path(output_dir, "Correlation_Matrix.png"), width = 800, height = 800)
corrplot::corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", 
                   tl.col = "black", tl.srt = 45, addCoef.col = "black", 
                   number.cex = 0.7, diag = FALSE)
dev.off()

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
  feature_cols <- setdiff(feature_cols, vars_to_remove)
  
  cat("Variables removed due to high collinearity:", paste(vars_to_remove, collapse = ", "), "\n")
} else {
  cat("No highly correlated variable pairs found, retaining all variables.\n")
}

cat("Final retained feature variables:", paste(feature_cols, collapse = ", "), "\n")

# Save collinearity results
write.csv(cor_matrix, file.path(output_dir, "Correlation_Matrix.csv"), row.names = TRUE)

cat("Collinearity analysis completed, results saved to", output_dir, "\n")

# -------------------- 3. Initialization --------------------
results_list <- list()
roc_train_list <- list()
roc_test_list <- list()
pdp_avg_list <- lapply(seq_along(feature_cols), function(j) list(Value = NULL, Prediction = NULL, SD = NULL))
variables <- feature_cols

# -------------------- 4. Iterative Modeling (1x Negative Sampling) --------------------
cat("Processing negative sampling ratio: 1x\n")

for (i in seq_len(num_iterations)) {
  set.seed(seeds[i])
  df_pos <- df %>% filter(cases == 1)
  df_neg <- df %>% filter(cases == 0)
  if (nrow(df_neg) < nrow(df_pos)) {
    stop("Insufficient negative points for 1x sampling!")
  }
  df_neg_sample <- sample_n(df_neg, size = nrow(df_pos))  # 1x negative sampling
  df_balanced <- bind_rows(df_pos, df_neg_sample)
  
  set.seed(seeds[i])
  train_idx <- sample(1:nrow(df_balanced), size = round(train_prop * nrow(df_balanced)))
  train_full <- df_balanced[train_idx, ]
  test_full <- df_balanced[-train_idx, ]
  
  spdf <- st_as_sf(train_full, coords = c("XX", "YY"), crs = 4326)
  blocks <- cv_spatial(x = spdf, column = "cases", r = shape_r, k = num_folds,
                       size = 500000, hexagon = FALSE, selection = "random",
                       iteration = 100, biomod2 = TRUE, progress = FALSE)
  train_full$fold <- blocks$folds_ids
  train_full$cases <- as.numeric(as.character(train_full$cases))
  test_full$cases <- as.numeric(as.character(test_full$cases))
  
  fold_results <- list()
  for (fold_id in 1:num_folds) {
    train_data <- train_full %>% filter(fold != fold_id)
    valid_data <- train_full %>% filter(fold == fold_id)
    train_data_selected <- train_data[, c("cases", feature_cols)]
    i.var <- which(colnames(train_data_selected) %in% feature_cols)
    i.y <- which(colnames(train_data_selected) == "cases")
    fold_vector <- match(train_data$fold, sort(unique(train_data$fold)))
    
    cat("Iteration", i, "Fold", fold_id, "Training data samples:", nrow(train_data), "\n")
    cat("Training data cases distribution:", table(train_data$cases), "\n")
    
    brt_model <- tryCatch({
      set.seed(seeds[i])
      gbm.step(
        data = train_data_selected,
        gbm.x = i.var,
        gbm.y = i.y,
        family = "bernoulli",
        tree.complexity = 2,
        learning.rate = 0.005,
        bag.fraction = 0.75,
        max.trees = 10000,
        n.folds = length(unique(fold_vector)),
        fold.vector = fold_vector,
        site.weights = train_data$weight,
        step.size = 50,
        silent = TRUE,
        plot.main = FALSE
      )
    }, error = function(e) {
      cat("Iteration", i, "Fold", fold_id, "Model training failed:", e$message, "\n")
      NULL
    })
    
    if (is.null(brt_model)) {
      cat("Iteration", i, "Fold", fold_id, "Skipped, model is NULL\n")
      next
    }
    best_trees <- brt_model$n.trees
    
    saveRDS(brt_model, file = file.path(output_dir, paste0("BRT_model_iter_", i, "_fold_", fold_id, ".rds")))
    
    pred_prob_train <- predict.gbm(brt_model, train_data[, feature_cols], n.trees = best_trees, type = "response")
    pred_prob_test <- predict.gbm(brt_model, test_full[, feature_cols], n.trees = best_trees, type = "response")
    
    cat("Iteration", i, "Fold", fold_id, "Prediction probability range:", range(pred_prob_train, na.rm = TRUE), "\n")
    
    valid_idx <- !is.na(train_data$cases) & !is.na(pred_prob_train)
    train_data <- train_data[valid_idx, ]
    pred_prob_train <- pred_prob_train[valid_idx]
    
    cm_train <- tryCatch({
      confusionMatrix(as.factor(ifelse(pred_prob_train > 0.5, 1, 0)), as.factor(train_data$cases))
    }, error = function(e) {
      cat("Iteration", i, "Fold", fold_id, "Confusion matrix calculation failed:", e$message, "\n")
      NULL
    })
    
    cm_test <- tryCatch({
      confusionMatrix(as.factor(ifelse(pred_prob_test > 0.5, 1, 0)), as.factor(test_full$cases))
    }, error = function(e) {
      cat("Iteration", i, "Fold", fold_id, "Test confusion matrix calculation failed:", e$message, "\n")
      NULL
    })
    
    roc_train <- tryCatch({
      roc(train_data$cases, pred_prob_train, quiet = TRUE)
    }, error = function(e) {
      cat("Iteration", i, "Fold", fold_id, "ROC calculation failed:", e$message, "\n")
      NULL
    })
    
    roc_test <- tryCatch({
      roc(test_full$cases, pred_prob_test, quiet = TRUE)
    }, error = function(e) {
      cat("Iteration", i, "Fold", fold_id, "Test ROC calculation failed:", e$message, "\n")
      NULL
    })
    
    if (is.null(roc_train)) {
      cat("Iteration", i, "Fold", fold_id, "Skipped, ROC is NULL\n")
      next
    }
    
    roc_train_list[[length(roc_train_list) + 1]] <- data.frame(
      FPR = 1 - roc_train$specificities,
      TPR = roc_train$sensitivities,
      Iteration = i,
      Fold = fold_id
    )
    if (!is.null(roc_test)) {
      roc_test_list[[length(roc_test_list) + 1]] <- data.frame(
        FPR = 1 - roc_test$specificities,
        TPR = roc_test$sensitivities,
        Iteration = i,
        Fold = fold_id
      )
    }
    
    var_imp <- summary(brt_model, n.trees = best_trees, plotit = FALSE) %>%
      rename(Variable = var, RelativeInfluence = rel.inf)
    var_imp$Iteration <- i
    var_imp$Fold <- fold_id
    
    cutoff_value <- tryCatch({
      coords <- pROC::coords(roc_train, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
      cat("Iteration", i, "Fold", fold_id, "Cutoff:", coords$threshold[1], 
          "Sensitivity:", coords$sensitivity[1], "Specificity:", coords$specificity[1], "\n")
      as.numeric(coords$threshold[1])
    }, error = function(e) {
      cat("Iteration", i, "Fold", fold_id, "Cutoff calculation failed:", e$message, "Using default 0.5\n")
      0.5
    })
    
    if (length(cutoff_value) > 1 || is.na(cutoff_value) || !is.finite(cutoff_value) || cutoff_value < 0 || cutoff_value > 1) {
      cat("Iteration", i, "Fold", fold_id, "Invalid cutoff:", cutoff_value, "Using default 0.5\n")
      cutoff_value <- 0.5
    }
    
    for (j in seq_along(variables)) {
      pdp_obj <- tryCatch({ 
        plot.gbm(brt_model, i.var = variables[j], n.trees = best_trees, return.grid = TRUE) 
      }, error = function(e) NULL)
      if (!is.null(pdp_obj)) {
        pdp_prob <- 1 / (1 + exp(-pdp_obj$y))
        if (i == 1 && fold_id == 1) {
          pdp_avg_list[[j]]$Value <- pdp_obj[[variables[j]]]
          pdp_avg_list[[j]]$Prediction <- pdp_prob
          pdp_avg_list[[j]]$SD <- pdp_prob^2
        } else {
          interp_pred <- approx(pdp_obj[[variables[j]]], pdp_prob, xout = pdp_avg_list[[j]]$Value, rule = 2)$y
          pdp_avg_list[[j]]$Prediction <- pdp_avg_list[[j]]$Prediction + interp_pred
          pdp_avg_list[[j]]$SD <- pdp_avg_list[[j]]$SD + interp_pred^2
        }
      }
    }
    
    fold_results[[fold_id]] <- list(
      Model = brt_model,
      Accuracy_Train = if (!is.null(cm_train)) cm_train$overall["Accuracy"] else NA,
      Accuracy_Test = if (!is.null(cm_test)) cm_test$overall["Accuracy"] else NA,
      AUC_Train = if (!is.null(roc_train)) as.numeric(roc_train$auc) else NA,
      AUC_Test = if (!is.null(roc_test)) as.numeric(roc_test$auc) else NA,
      F1 = if (!is.null(cm_test)) cm_test$byClass["F1"] else NA,
      Sensitivity = if (!is.null(cm_test)) cm_test$byClass["Sensitivity"] else NA,
      Specificity = if (!is.null(cm_test)) cm_test$byClass["Specificity"] else NA,
      VarImportance = var_imp,
      Cutoff = cutoff_value
    )
    
    cat("Iteration", i, "Fold", fold_id, "Stored Cutoff:", fold_results[[fold_id]]$Cutoff, "\n")
  }
  
  results_list[[i]] <- list(Seed = seeds[i], Fold_Results = fold_results)
  cat("Iteration", i, "Stored folds in results_list:", length(results_list[[i]]$Fold_Results), "\n")
}

# -------------------- 5. PDP Average and Plotting --------------------
for (j in seq_along(variables)) {
  pdp_avg_list[[j]]$Prediction <- pdp_avg_list[[j]]$Prediction / (num_iterations * num_folds)
  pdp_avg_list[[j]]$SD <- sqrt(pdp_avg_list[[j]]$SD / (num_iterations * num_folds) - pdp_avg_list[[j]]$Prediction^2)
  pdp_avg_list[[j]] <- data.frame(Value = pdp_avg_list[[j]]$Value, 
                                  Prediction = pdp_avg_list[[j]]$Prediction, 
                                  SD = pdp_avg_list[[j]]$SD)
  
  pdp_avg_list[[j]]$lower <- pdp_avg_list[[j]]$Prediction - 1.96 * pdp_avg_list[[j]]$SD
  pdp_avg_list[[j]]$upper <- pdp_avg_list[[j]]$Prediction + 1.96 * pdp_avg_list[[j]]$SD
  
  smoothed_lower <- predict(loess(pdp_avg_list[[j]]$lower ~ pdp_avg_list[[j]]$Value, span = 0.2))
  smoothed_upper <- predict(loess(pdp_avg_list[[j]]$upper ~ pdp_avg_list[[j]]$Value, span = 0.2))
  
  hist_data <- df[[variables[j]]]
  
  pdp_range <- range(pdp_avg_list[[j]]$Value, na.rm = TRUE)
  pdp_breaks <- pretty(pdp_avg_list[[j]]$Value, n = 10)
  
  pdp_plot <- ggplot(pdp_avg_list[[j]], aes(x = Value)) +
    geom_smooth(aes(y = Prediction), method = "loess", span = 0.2, color = "blue", se = FALSE) +
    geom_ribbon(aes(ymin = smoothed_lower, ymax = smoothed_upper), fill = "lightblue", alpha = 0.3) +
    annotate("text", x = Inf, y = Inf, label = "95% CI", hjust = 1.1, vjust = 1.1, size = 5, family = "Times", color = "black") +
    xlab(variables[j]) +
    ylab("Predicted Probability") +
    scale_x_continuous(limits = pdp_range, breaks = pdp_breaks) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      text = element_text(family = "Times", size = 25),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    )
  
  hist_plot <- ggplot(data.frame(x = hist_data), aes(x = x)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "grey", color = "black", alpha = 0.5) +
    xlab(NULL) +
    ylab("Density") +
    scale_x_continuous(limits = pdp_range, breaks = pdp_breaks) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      text = element_text(family = "Times", size = 20),
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
    )
  
  combined_plot <- grid.arrange(pdp_plot, hist_plot, nrow = 2, heights = c(3, 1))
  
  ggsave(
    filename = file.path(output_dir, "PDP_Plots", paste0("PDP_", variables[j], "_with_histogram.png")),
    plot = combined_plot,
    width = 10,
    height = 10,
    dpi = 300
  )
}
save(pdp_avg_list, file = file.path(output_dir, "pdp_avg_list.RData"))

# -------------------- 6. ROC Curves --------------------
roc_train_all <- do.call(rbind, roc_train_list)
roc_test_all <- do.call(rbind, roc_test_list)

fpr_grid <- seq(0, 1, by = 0.01)
interpolate_roc <- function(roc_data, fpr_grid) {
  tpr_interpolated <- matrix(NA, nrow = num_iterations * num_folds, ncol = length(fpr_grid))
  idx <- 1
  for (i in 1:num_iterations) {
    for (f in 1:num_folds) {
      roc_k <- roc_data %>% filter(Iteration == i, Fold = f)
      if (nrow(roc_k) > 0) {
        tpr_interpolated[idx, ] <- approx(roc_k$FPR, roc_k$TPR, xout = fpr_grid, rule = 2)$y
        idx <- idx + 1
      }
    }
  }
  return(tpr_interpolated[1:(idx-1), ])
}

tpr_train_interpolated <- interpolate_roc(roc_train_all, fpr_grid)
tpr_test_interpolated <- interpolate_roc(roc_test_all, fpr_grid)

roc_train_summary <- data.frame(
  FPR = fpr_grid,
  Mean_TPR = colMeans(tpr_train_interpolated, na.rm = TRUE),
  SD_TPR = apply(tpr_train_interpolated, 2, sd, na.rm = TRUE)
) %>%
  mutate(
    Lower_CI = Mean_TPR - 1.96 * SD_TPR / sqrt(num_iterations * num_folds),
    Upper_CI = Mean_TPR + 1.96 * SD_TPR / sqrt(num_iterations * num_folds)
  )

roc_test_summary <- data.frame(
  FPR = fpr_grid,
  Mean_TPR = colMeans(tpr_test_interpolated, na.rm = TRUE),
  SD_TPR = apply(tpr_test_interpolated, 2, sd, na.rm = TRUE)
) %>%
  mutate(
    Lower_CI = Mean_TPR - 1.96 * SD_TPR / sqrt(num_iterations * num_folds),
    Upper_CI = Mean_TPR + 1.96 * SD_TPR / sqrt(num_iterations * num_folds)
  )

# -------------------- 7. Save Models --------------------
saveRDS(results_list, file = file.path(output_dir, "results_list.rds"))
cat("results_list saved to", file.path(output_dir, "results_list.rds"), "\n")

# -------------------- 8. Output Variable Importance, Accuracy, and Cutoff Tables --------------------
var_imp_all <- do.call(rbind, lapply(seq_len(num_iterations), function(i) {
  fold_results <- results_list[[i]]$Fold_Results
  if (is.null(fold_results) || length(fold_results) == 0) return(data.frame())
  do.call(rbind, lapply(seq_len(num_folds), function(fold_id) {
    if (fold_id <= length(fold_results) && !is.null(fold_results[[fold_id]]$VarImportance)) {
      fold_results[[fold_id]]$VarImportance
    } else {
      data.frame()
    }
  }))
}))

eval_results <- data.frame(
  Iteration = integer(),
  Fold = integer(),
  AUC_Train = numeric(),
  AUC_Test = numeric(),
  Accuracy_Train = numeric(),
  Accuracy_Test = numeric(),
  F1 = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  Cutoff = numeric()
)

for (i in seq_len(num_iterations)) {
  fold_results <- results_list[[i]]$Fold_Results
  if (is.null(fold_results) || length(fold_results) == 0) next
  for (fold_id in 1:num_folds) {
    if (fold_id <= length(fold_results) && !is.null(fold_results[[fold_id]])) {
      fr <- fold_results[[fold_id]]
      tmp_eval <- data.frame(
        Iteration = i,
        Fold = fold_id,
        AUC_Train = as.numeric(fr$AUC_Train %||% NA),
        AUC_Test = as.numeric(fr$AUC_Test %||% NA),
        Accuracy_Train = as.numeric(fr$Accuracy_Train %||% NA),
        Accuracy_Test = as.numeric(fr$Accuracy_Test %||% NA),
        F1 = as.numeric(fr$F1 %||% NA),
        Sensitivity = as.numeric(fr$Sensitivity %||% NA),
        Specificity = as.numeric(fr$Specificity %||% NA),
        Cutoff = if (!is.null(fr$Cutoff) && !is.na(fr$Cutoff) && is.finite(fr$Cutoff)) as.numeric(fr$Cutoff) else NA
      )
    } else {
      tmp_eval <- data.frame(
        Iteration = i,
        Fold = fold_id,
        AUC_Train = NA,
        AUC_Test = NA,
        Accuracy_Train = NA,
        Accuracy_Test = NA,
        F1 = NA,
        Sensitivity = NA,
        Specificity = NA,
        Cutoff = NA
      )
    }
    eval_results <- rbind(eval_results, tmp_eval)
  }
}

var_imp_summary <- var_imp_all %>%
  group_by(Variable) %>%
  summarise(
    Mean_RelativeInfluence = mean(RelativeInfluence, na.rm = TRUE),
    SD_RelativeInfluence = sd(RelativeInfluence, na.rm = TRUE)
  ) %>%
  arrange(desc(Mean_RelativeInfluence))

var_imp_per_model <- var_imp_all %>%
  select(Iteration, Fold, Variable, RelativeInfluence) %>%
  arrange(Iteration, Fold, desc(RelativeInfluence))

cutoff_summary <- eval_results %>%
  filter(!is.na(Cutoff) & is.finite(Cutoff)) %>%
  summarise(
    Mean_Cutoff = mean(Cutoff, na.rm = TRUE),
    SD_Cutoff = sd(Cutoff, na.rm = TRUE),
    Min_Cutoff = min(Cutoff, na.rm = TRUE),
    Max_Cutoff = max(Cutoff, na.rm = TRUE)
  )

cutoff_per_model <- eval_results %>%
  select(Iteration, Fold, Cutoff) %>%
  arrange(Iteration, Fold) %>%
  filter(!is.na(Cutoff) & is.finite(Cutoff))

accuracy_summary <- eval_results %>%
  summarise(
    Mean_Accuracy_Train = mean(Accuracy_Train, na.rm = TRUE),
    SD_Accuracy_Train = sd(Accuracy_Train, na.rm = TRUE),
    Mean_Accuracy_Test = mean(Accuracy_Test, na.rm = TRUE),
    SD_Accuracy_Test = sd(Accuracy_Test, na.rm = TRUE),
    Mean_AUC_Train = mean(AUC_Train, na.rm = TRUE),
    SD_AUC_Train = sd(AUC_Train, na.rm = TRUE),
    Mean_AUC_Test = mean(AUC_Test, na.rm = TRUE),
    SD_AUC_Test = sd(AUC_Test, na.rm = TRUE),
    Mean_F1 = mean(F1, na.rm = TRUE),
    SD_F1 = sd(F1, na.rm = TRUE),
    Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE),
    SD_Sensitivity = sd(Sensitivity, na.rm = TRUE),
    Mean_Specificity = mean(Specificity, na.rm = TRUE),
    SD_Specificity = sd(Specificity, na.rm = TRUE)
  )

output_list <- list(
  Variable_Importance_Per_Model = var_imp_per_model,
  Variable_Importance_Summary = var_imp_summary,
  Cutoff_Per_Model = cutoff_per_model,
  Cutoff_Summary = cutoff_summary,
  Accuracy_Summary = accuracy_summary
)

write_xlsx(
  output_list,
  path = file.path(output_dir, "Variable_Importance_Accuracy_and_Cutoff.xlsx")
)

# Visualize cutoff distribution
if (nrow(cutoff_per_model) > 0) {
  p_cutoff <- ggplot(cutoff_per_model, aes(x = Cutoff)) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.5, color = "black") +
    xlab("Cutoff Value") +
    ylab("Frequency") +
    theme_minimal() +
    theme(
      text = element_text(family = "Times", size = 20),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  ggsave(
    filename = file.path(output_dir, "Cutoff_Distribution.png"),
    plot = p_cutoff,
    width = 8,
    height = 6,
    dpi = 300
  )
  cat("Cutoff distribution plot saved to", file.path(output_dir, "Cutoff_Distribution.png"), "\n")
}

cat("\nAll PDP plots, ROC curves, model RDS files, variable importance tables, and accuracy tables saved to", output_dir, "\n")
cat("\n1x negative sampling analysis completed.\n")
