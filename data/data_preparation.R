library(readxl)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(glmnet)
library(glmnetUtils)
library(caret)
library(pROC)
library(DescTools)

##read occurrence data
final_dt1 <- read_excel("occurrence.xlsx")

final_dt1 <- final_dt1 %>%
  mutate(occurrence = 1)%>%
  mutate(year = as.integer(year))
num_rows <- nrow(final_dt1)

##read absence data
final_dt0 <- read.csv("absence.csv")

#data balance
final_dt0 <- final_dt0 %>%
  sample_n(size = num_rows)

final_dt0 <- final_dt0 %>%
  mutate(
    occurrence = 0,  # 为occurrence赋值为0
    ID = as.character(ID),  
    temp_id = row_number(),
    ID = if_else(is.na(ID), paste0("0_", temp_id), ID) 
  ) %>%
  select(-temp_id)

##combine occurrence and absence
final_dt <- bind_rows(final_dt0, final_dt1)
saveRDS(final_dt,"final_dtused.rds")

final_dt<-as.data.frame(final_dt)
final_dt <- final_dt %>% 
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

#re-scale
predictors_names<- c("Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
               "Deciduous_Needleleaf_Forests","Deciduous_Broadleaf_Forests","Mixed_Forests",
               "Closed_Shrublands","Open_Shrublands","Woody_Savannas",
               "Savannas","Grasslands","Permanent_Wetlands",
               "Croplands","Urban_and_Builtup","Cropland_Natural_Vegetation_Mosaics",
               "Snow_and_Ice","Barren_or_Sparsely_Vegetated","ndvi","evi", "pop", "prec", 
                     "tmin", "tmax", "sp","ws","rh", "elevation", "urban_acc_30s")
numeric_cols <- predictors_names[sapply(final_dt[, predictors_names], is.numeric)]
final_dt[, numeric_cols] <- lapply(final_dt[, numeric_cols], scale, center = TRUE, scale = TRUE)

#set train data (70%) and test data (30%)
set.seed(123)
trainIndex <- createDataPartition(final_dt$occurrence, p = .7, 
                                  list = FALSE, 
                                  times = 1)
saveRDS(final_dt,"final_dtused.rds")


trainSet <- final_dt[trainIndex, ]
testSet <- final_dt[-trainIndex, ]
saveRDS(testSet,"testset30%.rds")
saveRDS(trainSet,"trainset70%.rds")
