#---------------- 15km/45km data preparation ----------------

# Load libraries
library(data.table)
library(dplyr)
library(ranger)

# Read data
final_dt_15 <- readRDS("final_dtused-15km.rds")
final_dt_45 <- readRDS("final_dtused-45km.rds")

# Convert to data.frame and calculate mite variable
process_data <- function(dt) {
  df <- as.data.frame(dt)
  df$mite <- pmax(df$ld_mean, df$ls_mean, df$lp_mean, na.rm = TRUE)
  return(df)
}

df_15 <- process_data(final_dt_15)
df_45 <- process_data(final_dt_45)

# Define land cover names
land_cover_names <- c(
  "Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
  "Deciduous_Needleleaf_Forests", "Deciduous_Broadleaf_Forests", "Mixed_Forests",
  "Closed_Shrublands", "Open_Shrublands", "Woody_Savannas", "Savannas", "Grasslands",
  "Permanent_Wetlands", "Croplands", "Urban_and_Builtup", "Cropland_Natural_Vegetation_Mosaics",
  "Snow_and_Ice", "Barren_or_Sparsely_Vegetated"
)

# Rename land cover columns
names(df_15)[1:17] <- land_cover_names
names(df_45)[1:17] <- land_cover_names

# Define predictor variables
predictors_names <- c(
  land_cover_names, "ndvi", "evi", "pop", "prec", "tmin", "tmax", 
  "sp", "ws", "rh", "elevation", "urban_acc_30s", "mite", "rodent_richness"
)

# Scale numeric predictors and save parameters
scale_and_save <- function(df, file_suffix) {
  numeric_cols <- predictors_names[sapply(df[, predictors_names], is.numeric)]
  
  # Save scaling parameters
  scaling_params <- lapply(df[, numeric_cols], function(x) {
    list(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
  })
  saveRDS(scaling_params, paste0(file_suffix, "_scaling_params.rds"))
  
  # Scale the data
  df[, numeric_cols] <- scale(df[, numeric_cols], center = TRUE, scale = TRUE)
  return(df)
}

df_15 <- scale_and_save(df_15, "15")
df_45 <- scale_and_save(df_45, "45")

# Balance datasets and save
balance_and_save <- function(df, file_name) {
  df_clean <- df[complete.cases(df[, predictors_names]), ]
  df_0 <- df_clean %>% filter(occurrence == 0)
  df_1 <- df_clean %>% filter(occurrence == 1)  
  df_0_sampled <- df_0 %>% sample_n(nrow(df_1))
  df_balanced <- bind_rows(df_0_sampled, df_1) %>% sample_frac(1) 
  saveRDS(df_balanced, file_name)
}
balance_and_save(df_15, "15_df_balanced_new.rds")
balance_and_save(df_45, "45_df_balanced_new.rds")

# data preparation for prediction

library(raster)
library(dplyr)

#new data with three mite species info and rodent richness
all_rawcovariates_stack<-readRDS("predictrawdata_raster2020.rds")

#need to stack ld_repro.tif, lp_repro.tif, ls_repro.tif, richness_repro.tif
ld_mean<-raster("ld_repro.tif")
ls_mean<-raster("ls_repro.tif")
lp_mean<-raster("lp_repro.tif")
rodent_richness<-raster("richness_repro.tif")

library(raster)
ref_raster <- all_rawcovariates_stack[[1]]  

ld_mean_resampled <- resample(ld_mean, ref_raster, method = "bilinear")
lp_mean_resampled <- resample(lp_mean, ref_raster, method = "bilinear")
ls_mean_resampled <- resample(ls_mean, ref_raster, method = "bilinear")
rodent_richness_resampled <- resample(rodent_richness, ref_raster, method = "bilinear")

names(ld_mean_resampled) <- "ld_mean"
names(lp_mean_resampled) <- "lp_mean"
names(ls_mean_resampled) <- "ls_mean"
names(rodent_richness_resampled) <- "rodent_richness"

mite_max_prob <- overlay(ld_mean_resampled, lp_mean_resampled, ls_mean_resampled, fun = function(x, y, z) {
  pmax(x, y, z, na.rm = TRUE)
})

names(mite_max_prob) <- "mite"

full_stack <- stack(all_rawcovariates_stack,
                    mite_max_prob,
                    rodent_richness_resampled)

library(data.table)
full_values <- as.data.table(raster::as.data.frame(full_stack, xy = TRUE, na.rm = FALSE))
setnames(full_values, 
         old = c("A2020001_subdataset3_band1", "A2020001_subdataset3_band2", "A2020001_subdataset3_band3", "A2020001_subdataset3_band4", "A2020001_subdataset3_band5", 
                 "A2020001_subdataset3_band6", "A2020001_subdataset3_band7", "A2020001_subdataset3_band8", "A2020001_subdataset3_band9", 
                 "A2020001_subdataset3_band10", "A2020001_subdataset3_band11", "A2020001_subdataset3_band12", "A2020001_subdataset3_band13", 
                 "A2020001_subdataset3_band14", "A2020001_subdataset3_band15", "A2020001_subdataset3_band16", "A2020001_subdataset3_band17",
                 "X2020_subdataset1_average_5km_clip" ,        
                 "X2020_subdataset2_average_5km_clip", "pop_2020_5km",                          
                 "annual_prec_2020",                       
                 "annual_average_tmax_2020",            
                 "annual_average_tmin_2020",               
                 "wc2.1_2.5m_elev_5km" ,                  
                 "ws_annual_mean_2020_5km_clip", 
                 "sp_annual_mean_2020_5km_clip", 
                 "annual_mean_rh_2020_5km_clip",
                 "urban_accessibility_2015_5km_clip"),
         new = c("Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests", 
                 "Deciduous_Needleleaf_Forests", "Deciduous_Broadleaf_Forests", "Mixed_Forests",
                 "Closed_Shrublands", "Open_Shrublands", "Woody_Savannas",
                 "Savannas", "Grasslands", "Permanent_Wetlands",
                 "Croplands", "Urban_and_Builtup", "Cropland_Natural_Vegetation_Mosaics",
                 "Snow_and_Ice", "Barren_or_Sparsely_Vegetated","ndvi","evi","pop",
                 "prec","tmax","tmin","elevation","ws","sp","rh","urban_acc_30s"))

head(colnames(full_values))

#use scaling parameters to rescale it
scaling_params<-readRDS("45_scaling_params.rds")
vars_to_scale <- intersect(names(scaling_params), colnames(full_values))
for (col in vars_to_scale) {
  mu <- scaling_params[[col]]$mean
  sigma <- scaling_params[[col]]$sd
  full_values[[col]] <- (full_values[[col]] - mu) / sigma
}
saveRDS(full_values, "45_predict_scaled_covariates_2020.rds")

scaling_params<-readRDS("15_scaling_params.rds")
vars_to_scale <- intersect(names(scaling_params), colnames(full_values))
for (col in vars_to_scale) {
  mu <- scaling_params[[col]]$mean
  sigma <- scaling_params[[col]]$sd
  full_values[[col]] <- (full_values[[col]] - mu) / sigma
}
saveRDS(full_values, "15_predict_scaled_covariates_2020.rds")



#---------------- thinning ----------------

# Load libraries
library(spatstat)
library(sf)
library(sf)
library(dplyr)
library(FNN)

#load data
df_balanced<-readRDS("df_balanced.rds")
df_sf <- st_as_sf(df_balanced, coords = c("coordinates_long", "coordinates_lat"), crs = 4326)
df_sf <- st_transform(df_sf, crs = 3857)

# mini distance -5km
min_distance <- 5000  
coords <- st_coordinates(df_sf)

distances <- dist(coords)
keep_indices <- c()
remaining_indices <- 1:nrow(df_sf)

simple_grid_thin <- function(df_sf, cell_size = 5000) {
  coords <- st_coordinates(df_sf)
  
  x_grid <- floor(coords[,1] / (cell_size/111000))
  y_grid <- floor(coords[,2] / (cell_size/111000))
  grid_id <- paste(x_grid, y_grid, sep = "_")

  set.seed(123)
  keep_indices <- c()
  for(gid in unique(grid_id)) {
    grid_points <- which(grid_id == gid)
    keep_indices <- c(keep_indices, sample(grid_points, 1))
  }
  
  return(keep_indices)
}

keep_indices <- simple_grid_thin(df_sf, cell_size = 5000)
df_thinned <- df_balanced[keep_indices, ]
saveRDS(df_thinned,"df_thinned.rds")


