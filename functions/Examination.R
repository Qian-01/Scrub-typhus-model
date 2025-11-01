
library(ggplot2)
library(dplyr)
library(readxl)
library(broom)

df <- read_excel("calculation of result.xlsx")
df <- read_excel("weighted es.xls", sheet = 1)
df_plot <- df %>%
  rename(
    country = Country
  ) %>%
  filter(country %in% c("Mainland China","Japan"))%>%
  mutate(
    incidence_plot = ifelse(incidence <= 0, 0.01, incidence),
    ES_Group = ifelse(esw < 0.421527266, "< cutoff (Blue)", ">= cutoff (Red)")
  )


library(ggpubr)

summary_table <- df_plot %>%
  group_by(country, ES_Group) %>%
  summarise(mean_incidence = mean(incidence_plot), .groups = "drop")

p_values <- df_plot %>%
  group_by(country) %>%
  summarise(
    p_value = tryCatch({
      wilcox.test(incidence_plot ~ ES_Group)$p.value
    }, error = function(e) NA_real_)
  )

annotations <- left_join(summary_table, p_values, by = "country") %>%
  group_by(country) %>%
  mutate(
    label = paste0(ES_Group, ": ", format(mean_incidence, digits = 2)),
    p_text = ifelse(row_number() == 1, paste0("p = ", format(p_value[1], digits = 2)), NA)
  )
    
# ----------  A ----------
p1 <- ggplot(df_plot, aes(x = incidence_plot, y = country, fill = ES_Group)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "black") +
  scale_x_log10() +
  scale_fill_manual(values = c("< cutoff (Blue)" = "blue", ">= cutoff (Red)" = "red")) +
  labs(
    x = "Mean annual incidence (per 100,000 population)",
    y = NULL,
    fill = "Environmental Suitability\n(cutoff = 0.42)",
    title = "A"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.45, 0.80),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = NA)
  )
p1

annotations <- left_join(summary_table, p_values, by = "country") %>%
  group_by(country) %>%
  mutate(
    label = paste0(ES_Group, ": ", format(mean_incidence, digits = 2)),
    p_text = case_when(
      row_number() == 1 & !is.na(p_value) ~ ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", format(p_value, digits = 2))),
      TRUE ~ NA_character_
    )
  )



# ----------  B ----------
# Poisson
# Spearman correlation
cor_test <- cor.test(df_plot$incidence_plot, df_plot$esw, method = "spearman")

cor_test$estimate   # correlation coefficient (rho)
cor_test$p.value    # p-value

glm_model <- glm(incidence_plot~ esw, data = df_plot, family = "poisson")
summary(glm_model)
newdata <- data.frame(esw = seq(0, 1, length.out = 200))
pred <- predict(glm_model, newdata, se.fit = TRUE, type = "link")

newdata <- newdata %>%
  mutate(
    fit = exp(pred$fit),
    lower = exp(pred$fit - 1.96 * pred$se.fit),
    upper = exp(pred$fit + 1.96 * pred$se.fit)
  )


library(ggplot2)

p2 <- ggplot() +
  geom_point(data = df_plot, aes(x = esw, y = incidence_plot, color = "Data points"), alpha = 0.6) +
  geom_line(data = newdata, aes(x = esw, y = fit, color = "GLM Regression line"), size = 1.1) +
  geom_ribbon(data = newdata, aes(x = esw, ymin = lower, ymax = upper, fill = "95% CI"), alpha = 0.3) +
  scale_color_manual(name = NULL, values = c(
    "Data points" = "blue",
    "GLM Regression line" = "red"
  )) +
  scale_fill_manual(name = NULL, values = c("95% CI" = "red")) +
  labs(
    x = "Population Weighted Environmental Suitability",
    y = "Reported Incidence (per 100,000)",
    title = "B") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.25, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_blank()
  )

p2


# ---------- combine----------
library(patchwork)
fig2 <- p1 + p2
ggsave("Figure3.png", fig2, width = 14.65, height = 5.27, units = "in", dpi = 300)

