# -------------------- 0. Load Packages --------------------
library(tidyverse)
library(pROC)
library(broom)
library(arm)
library(gridExtra)

# Clear memory
cat("Clearing memory...\n")
rm(list = ls())
gc()
cat("Memory cleared.\n")

# -------------------- 1. Setup and Parameters --------------------
output_dir <- "./outputs/weight_results/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- 2. Data Preparation --------------------
data <- read_csv("./data/soc_data/mite_soc.csv") %>%
  mutate(log_gdp = log(gdp + 1),
         log_pop = log(pop + 1)) %>%
  drop_na()

var_name <- c("log_gdp", "log_pop", "urban_accessibility", "cases", "XX", "YY")
data_model <- data[, c("code", var_name)]
positive <- data_model %>% filter(cases == 1)
negative <- data_model %>% filter(cases == 0)

# -------------------- 3. Iterative Modeling --------------------
n_iter <- 50
auc_list <- numeric(n_iter)
acc_list <- numeric(n_iter)
weight_list_positive <- list()
weight_list_all <- list()
coef_all <- list()
roc_list <- list()
model_list <- character(n_iter)

for (i in 1:n_iter) {
  set.seed(123 + i)
  cat("Running iteration", i, "of", n_iter, "...\n")
  
  # Sample negative points (1:1 ratio)
  neg_sample <- sample_n(negative, size = nrow(positive))
  survey <- bind_rows(positive, neg_sample) %>% mutate(if_survey = cases)
  
  # Fit Bayesian logistic regression
  glm_formula <- as.formula("if_survey ~ log_gdp + log_pop + urban_accessibility")
  glm <- tryCatch({
    bayesglm(glm_formula, data = survey, family = binomial(link = "logit"),
             maxit = 100, prior.scale = 2.5, prior.df = 1)
  }, error = function(e) {
    cat("Iteration", i, "model fitting failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(glm)) next
  
  # Save model
  model_file <- file.path(output_dir, paste0("bayesglm_model_iter_", i, ".rds"))
  saveRDS(glm, model_file)
  model_list[i] <- model_file
  
  # Predict probabilities
  prob <- tryCatch({
    predict(glm, newdata = survey, type = "response")
  }, error = function(e) {
    cat("Iteration", i, "prediction failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(prob)) next
  
  # Calculate ROC and AUC
  roc_obj <- tryCatch({
    roc(survey$if_survey, prob, quiet = TRUE)
  }, error = function(e) {
    cat("Iteration", i, "ROC calculation failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (!is.null(roc_obj)) {
    auc_list[i] <- auc(roc_obj)
    roc_list[[i]] <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities,
      iter = i
    )
  } else {
    auc_list[i] <- NA
  }
  
  # Calculate accuracy
  acc_list[i] <- mean(ifelse(prob > 0.5, 1, 0) == survey$if_survey, na.rm = TRUE)
  
  # Calculate weights
  prob_pos <- predict(glm, newdata = positive, type = "response") %>% pmin(0.9) %>% pmax(0.1)
  weight_list_positive[[i]] <- (1 / prob_pos) / mean(1 / prob_pos, na.rm = TRUE)
  
  prob_all <- predict(glm, newdata = data_model, type = "response") %>% pmin(0.9) %>% pmax(0.1)
  weight_all <- ifelse(data_model$cases == 1, 1 / prob_all, 1 / (1 - prob_all))
  weight_list_all[[i]] <- weight_all / mean(weight_all, na.rm = TRUE)
  
  # Extract standardized coefficients
  survey_scaled <- survey
  survey_scaled[, c("log_gdp", "log_pop", "urban_accessibility")] <-
    scale(survey_scaled[, c("log_gdp", "log_pop", "urban_accessibility")])
  glm_scaled <- tryCatch({
    glm(glm_formula, data = survey_scaled, family = binomial(link = "logit"))
  }, error = function(e) {
    cat("Iteration", i, "standardized model fitting failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (!is.null(glm_scaled)) {
    coef_df <- tidy(glm_scaled) %>%
      dplyr::select(term, estimate) %>%
      filter(term != "(Intercept)") %>%
      mutate(iter = i)
    coef_all[[i]] <- coef_df
  }
}

# -------------------- 4. Save Weights --------------------
positive$weight <- rowMeans(do.call(cbind, weight_list_positive), na.rm = TRUE)
data_model$weight <- rowMeans(do.call(cbind, weight_list_all), na.rm = TRUE)
write_csv(positive, file.path(output_dir, "mite_weight_positive.csv"))
write_csv(data_model, file.path(output_dir, "mite_weight_all.csv"))

# Save model list
model_df <- data.frame(Iteration = 1:n_iter, Model_File = model_list)
write_csv(model_df, file.path(output_dir, "Model_File_List.csv"))

# -------------------- 5. ROC Curve Plot --------------------
roc_all <- bind_rows(roc_list)
roc_summary <- roc_all %>%
  group_by(fpr) %>%
  summarise(tpr_mean = mean(tpr, na.rm = TRUE), tpr_sd = sd(tpr, na.rm = TRUE), .groups = "drop")

p_roc <- ggplot(roc_summary, aes(x = fpr, y = tpr_mean)) +
  geom_line(color = "blue", size = 1.2) +
  geom_ribbon(aes(ymin = pmax(tpr_mean - tpr_sd, 0), ymax = pmin(tpr_mean + tpr_sd, 1)), alpha = 0.2, fill = "blue") +
  labs(title = "Mean ROC Curve (50 Iterations)", x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal(base_family = "Times") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    text = element_text(size = 12)
  )
ggsave(file.path(output_dir, "AUC_ROC_Combined.png"), width = 6, height = 6, dpi = 300)

# -------------------- 6. Accuracy Plot --------------------
accuracy_df <- data.frame(Iteration = 1:n_iter, Accuracy = acc_list)
p_acc <- ggplot(accuracy_df, aes(x = "", y = Accuracy)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Accuracy across 50 Iterations", y = "Accuracy", x = "") +
  theme_minimal(base_family = "Times") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    text = element_text(size = 12)
  )
ggsave(file.path(output_dir, "Accuracy_Boxplot.png"), width = 4, height = 6, dpi = 300)

# -------------------- 7. Variable Importance Plot --------------------
coef_df_all <- bind_rows(coef_all)
coef_summary <- coef_df_all %>%
  group_by(term) %>%
  summarise(mean_coef = mean(estimate, na.rm = TRUE), sd_coef = sd(estimate, na.rm = TRUE), .groups = "drop") %>%
  mutate(term = reorder(term, mean_coef))

p_coef <- ggplot(coef_summary, aes(x = term, y = mean_coef)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_coef - sd_coef, ymax = mean_coef + sd_coef), width = 0.2) +
  coord_flip() +
  labs(title = "Variable Importance (Standardized Coefficients)", x = "", y = "Mean Coefficient ± SD") +
  theme_minimal(base_family = "Times") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    text = element_text(size = 12)
  )
ggsave(file.path(output_dir, "Variable_Importance.png"), width = 8, height = 6, dpi = 300)

# -------------------- 8. Save Metrics --------------------
write_csv_safe <- function(data, path) {
  tryCatch({
    write_csv(data, path)
  }, error = function(e) {
    cat("Failed to write CSV file", path, ":", conditionMessage(e), "\n")
  })
}

auc_df <- data.frame(Iteration = 1:n_iter, AUC = round(auc_list, 4))
accuracy_df <- data.frame(Iteration = 1:n_iter, Accuracy = round(acc_list, 4))
mean_auc_df <- data.frame(Metric = "Mean_AUC", Value = round(mean(auc_list, na.rm = TRUE), 4))

write_csv_safe(auc_df, file.path(output_dir, "AUC_values.csv"))
write_csv_safe(accuracy_df, file.path(output_dir, "Accuracy_values.csv"))
write_csv_safe(mean_auc_df, file.path(output_dir, "Mean_AUC.csv"))

cat("Mean AUC:", round(mean(auc_list, na.rm = TRUE), 4), "\n")
cat("Analysis completed, results saved to", output_dir, "\n")
