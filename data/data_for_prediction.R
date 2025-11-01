
#stack all the global covariates 
library(mgcv)
library(raster)
library(dplyr)
library(sf)

#load 27 covariates
Water_bodies <- raster("A2020001_subdataset3_band1.tif")
Evergreen_Needleleaf_Forests <- raster("A2020001_subdataset3_band2.tif")
Evergreen_Broadleaf_Forests <- raster("A2020001_subdataset3_band3.tif")
Deciduous_Needleleaf_Forests <- raster("A2020001_subdataset3_band4.tif")
Deciduous_Broadleaf_Forests <- raster("A2020001_subdataset3_band5.tif")
Mixed_Forests <- raster("A2020001_subdataset3_band6.tif")
Closed_Shrublands <- raster("A2020001_subdataset3_band7.tif")
Open_Shrublands <- raster("A2020001_subdataset3_band8.tif")
Woody_Savannas <- raster("A2020001_subdataset3_band9.tif")
Savannas <- raster("A2020001_subdataset3_band10.tif")
Grasslands <- raster("A2020001_subdataset3_band11.tif")
Permanent_Wetlands <- raster("A2020001_subdataset3_band12.tif")
Croplands <- raster("A2020001_subdataset3_band13.tif")
Urban_and_Builtup <- raster("A2020001_subdataset3_band14.tif")
Cropland_Natural_Vegetation_Mosaics <- raster("A2020001_subdataset3_band15.tif")
Snow_and_Ice <- raster("A2020001_subdataset3_band16.tif")
Barren_or_Sparsely_Vegetated <- raster("A2020001_subdataset3_band17.tif")

evi <- raster("2020_subdataset2_average_5km_clip.tif")  
ndvi <- raster("2020_subdataset1_average_5km_clip.tif")  
pop <- raster("pop_2020_5km.tif") 
prec <- raster("annual_prec_2020.tif")
tmax <- raster("annual_average_tmax_2020.tif")
tmin <- raster("annual_average_tmin_2020.tif")
elevation <- raster("wc2.1_2.5m_elev_5km.tif")

ws <- raster("ws_annual_mean_2020_5km_clip.tif")
sp <- raster("sp_annual_mean_2020_5km_clip.tif")
rh <- raster("annual_mean_rh_2020_5km_clip.tif")

urban_acc_30s <- raster("urban_accessibility_2015_5km_clip.tif")


# Combine all raster into a single stack
all_rawcovariates_stack <- stack(Water_bodies, Evergreen_Needleleaf_Forests, Evergreen_Broadleaf_Forests,
                              Deciduous_Needleleaf_Forests, Deciduous_Broadleaf_Forests, Mixed_Forests,
                              Closed_Shrublands, Open_Shrublands, Woody_Savannas, Savannas, Grasslands,
                              Permanent_Wetlands, Croplands, Urban_and_Builtup, Cropland_Natural_Vegetation_Mosaics,
                              Snow_and_Ice, Barren_or_Sparsely_Vegetated,
                              evi, ndvi, pop, prec, tmax, tmin,ws,sp,rh, elevation,urban_acc_30s)
#save RDS of stacked raster
saveRDS(all_rawcovariates_stack,"predictrawdata_raster2020.rds")

library(mgcv)
library(raster)
library(dplyr)
library(sf)

all_rawcovariates_stack<-readRDS("predictrawdata_raster2020.rds")

library(data.table)
predict_rawdata <- as.data.table(getValues(all_rawcovariates_stack))
#rename
setnames(predict_rawdata, 
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
saveRDS(predict_rawdata ,"cor_rawdata_2020.rds")


#scale based on final_dt1
#load final_dt1
library(data.table)
final_dt_p<-readRDS("final_dtused-0527.rds")

final_dt_p <- final_dt_p %>% 
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

predictors_names<- c("Water_bodies", "Evergreen_Needleleaf_Forests", "Evergreen_Broadleaf_Forests",
                     "Deciduous_Needleleaf_Forests","Deciduous_Broadleaf_Forests","Mixed_Forests",
                     "Closed_Shrublands","Open_Shrublands","Woody_Savannas",
                     "Savannas","Grasslands","Permanent_Wetlands",
                     "Croplands","Urban_and_Builtup","Cropland_Natural_Vegetation_Mosaics",
                     "Snow_and_Ice","Barren_or_Sparsely_Vegetated","ndvi","evi", "pop", "prec", 
                     "tmin", "tmax", "sp","ws","rh", "elevation", "urban_acc_30s")
numeric_cols <- predictors_names[sapply(final_dt_p[, predictors_names], is.numeric)]

# final_dt1 is data frame with covariates used for scale for the model
means <- sapply(final_dt_p[, numeric_cols], mean, na.rm = TRUE)
sds <- sapply(final_dt_p[, numeric_cols], sd, na.rm = TRUE)

#predict_rawdata<-readRDS("cor_rawdata_2020.rds")
for (col in numeric_cols) {
  if (col %in% names(predict_rawdata)) {  # Check if the column exists in predict_rawdata
    predict_rawdata[[col]] <- (predict_rawdata[[col]] - means[[col]]) / sds[[col]]
  }
}
saveRDS(predict_rawdata ,"cor_data.rds")
