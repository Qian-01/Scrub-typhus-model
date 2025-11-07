
rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(raster)
  library(dplyr)
  library(sf)
  library(data.table)
  library(ranger)
  library(purrr)
})

# ---------------- load data ----------------
DF_RDS   <- "df_with_clusters.rds"
PRED_RDS <- "predict_scaled_covariates_2020.rds"
TPL_TIF  <- "annual_prec_2020.tif"

rasterOptions(overwrite = TRUE)  

df <- readRDS(DF_RDS)
predict_data <- readRDS(PRED_RDS)
template_raster <- raster(TPL_TIF)

predictors_names <- c(
  "Water_bodies","Evergreen_Needleleaf_Forests","Evergreen_Broadleaf_Forests",
  "Deciduous_Needleleaf_Forests","Deciduous_Broadleaf_Forests","Mixed_Forests",
  "Closed_Shrublands","Open_Shrublands","Woody_Savannas",
  "Savannas","Grasslands","Permanent_Wetlands",
  "Croplands","Urban_and_Builtup","Cropland_Natural_Vegetation_Mosaics",
  "Snow_and_Ice","Barren_or_Sparsely_Vegetated","ndvi","evi","pop","prec",
  "tmin","tmax","sp","ws","rh","elevation","urban_acc_30s","mite","rodent_richness"
)

prediction_data_full <- as.data.frame(predict_data)[, predictors_names, drop = FALSE]

# occurrence=0/1 
df_0 <- dplyr::filter(df, occurrence %in% c(0, "0"))
df_1 <- dplyr::filter(df, occurrence %in% c(1, "1"))
cat("df_0 rows:", nrow(df_0), "| df_1 rows:", nrow(df_1), "\n")

# ---------------- ----------------
TOTAL_ITERS <- 100
bootstrap_seeds <- 10000 + seq_len(TOTAL_ITERS)  

existing <- list.files(pattern = "^RF_pred_bootstrap_\\d{3}\\.tif$")
existing_idx <- as.integer(sub("^RF_pred_bootstrap_(\\d{3})\\.tif$", "\\1", existing))

todo_idx <- setdiff(seq_len(TOTAL_ITERS), existing_idx)
if (length(todo_idx) == 0) {
  message("test all predictions are done")
  quit(save = "no")
}
message("running following：", paste(sprintf("%03d", todo_idx), collapse = ", "))

# ---------------- parrallel ----------------
NUM_THREADS <- max(1, parallel::detectCores() - 2)  #
message("ranger cores used：", NUM_THREADS)

GTiff_OPTS <- c("TILED=YES", "BIGTIFF=YES", "COMPRESS=NONE")

.safe_unlink <- function(fn) {
  aux <- c(fn,
           paste0(fn, ".aux.xml"),
           paste0(fn, ".ovr"),
           sub("\\.tif$", ".tfw", fn))
  invisible(try(unlink(aux[file.exists(aux)], force = TRUE), silent = TRUE))
}

.move_or_copy <- function(src, dst, max_try = 3, sleep_sec = 0.5) {
  for (k in 1:max_try) {
    .safe_unlink(dst)
    ok <- try(file.rename(src, dst), silent = TRUE)
    if (isTRUE(ok)) return(invisible(TRUE))
    ok2 <- try(file.copy(src, dst, overwrite = TRUE), silent = TRUE)
    if (isTRUE(ok2)) {
      unlink(src, force = TRUE)
      return(invisible(TRUE))
    }
    Sys.sleep(sleep_sec)
  }
  stop(sprintf("cannot move object to location：%s -> %s", src, dst))
}

make_blocks <- function(r, target_blocks = 12) {
  nr <- nrow(r)
  n  <- max(1, min(target_blocks, nr))
  base <- floor(nr / n)
  extra <- nr %% n
  nrows <- rep(base, n)
  if (extra > 0) nrows[1:extra] <- nrows[1:extra] + 1
  row <- cumsum(c(1, head(nrows, -1)))
  data.frame(row = row, nrows = nrows, n = n)
}
bs <- make_blocks(template_raster, target_blocks = 12)
ncols_tr <- ncol(template_raster)

# ---------------- main fucntion, run, prediction, write----------------
for (boot_iter in todo_idx) {
  cat("Starting CLUSTER bootstrap iteration", boot_iter, "of", TOTAL_ITERS, "at", Sys.time(), "\n")
  set.seed(bootstrap_seeds[boot_iter])
  
  unique_clusters_0 <- unique(df_0$cluster)
  boot_clusters_0 <- sample(unique_clusters_0, size = length(unique_clusters_0), replace = TRUE)
  df_0_sampled <- map_dfr(boot_clusters_0, function(clust) dplyr::filter(df_0, cluster == clust))
  df_1_sampled <- df_1
  
  df_balanced <- bind_rows(df_1_sampled, df_0_sampled)
  df_balanced$occurrence <- as.factor(df_balanced$occurrence)
  
  formula_rf <- as.formula(paste("occurrence ~", paste(predictors_names, collapse = " + ")))
  
  wts <- if ("weight" %in% names(df_balanced)) df_balanced$weight else rep(1, nrow(df_balanced))
  rf.fit <- ranger(
    formula_rf,
    data         = df_balanced,
    num.trees    = 900,
    probability  = TRUE,
    importance   = "impurity",
    case.weights = wts,
    num.threads  = NUM_THREADS,
    keep.inbag   = FALSE
  )
  
  # ====== output ======
  out_fn  <- sprintf("RF_pred_bootstrap_%03d.tif", boot_iter)
  .safe_unlink(out_fn)  
  
  out_tmp <- tempfile(
    pattern = sprintf("RF_pred_bootstrap_%03d_", boot_iter),
    tmpdir  = tempdir(),
    fileext = ".tif"
  )
  .safe_unlink(out_tmp)
  try(raster::removeTmpFiles(h = 0), silent = TRUE)
  gc()
  
  # ====== write ======
  out_r <- NULL
  tryCatch({
    out_r <<- writeStart(template_raster, filename = out_tmp, format = "GTiff",
                         options = GTiff_OPTS, overwrite = TRUE)
  }, error = function(e) {
    try(raster::removeTmpFiles(h = 0), silent = TRUE)
    .safe_unlink(out_tmp)
    out_r <<- writeStart(template_raster, filename = out_tmp, format = "GTiff",
                         options = GTiff_OPTS, overwrite = TRUE)
  })
  
  cat("Making block-wise predictions for bootstrap", boot_iter, "...\n")
  for (i in seq_len(bs$n)) {
    start_row <- bs$row[i]
    nrows_blk <- bs$nrows[i]
    cell_start <- (start_row - 1) * ncols_tr + 1
    cell_end   <- (start_row + nrows_blk - 1) * ncols_tr
    idx <- cell_start:cell_end
    
    pred_block_df <- prediction_data_full[idx, , drop = FALSE]
    
    pred_block <- predict(
      rf.fit, data = pred_block_df, type = "response", num.threads = NUM_THREADS
    )$predictions[, 2]
    
    out_r <- writeValues(out_r, pred_block, start_row)
  }
  
  out_r <- writeStop(out_r)
  
  ok <- try(.move_or_copy(out_tmp, out_fn), silent = TRUE)
  if (inherits(ok, "try-error")) {
    Sys.sleep(0.5)
    ok <- try(.move_or_copy(out_tmp, out_fn), silent = TRUE)
    if (inherits(ok, "try-error")) {
      stop("fail to move：", out_tmp, " -> ", out_fn)
    }
  }
  cat("Saved", out_fn, "\n")
  
  rm(rf.fit, out_r, pred_block, pred_block_df); gc()
  cat("Completed bootstrap iteration", boot_iter, "at", Sys.time(), "\n\n")
}

message("all the prediction done：")
print(sort(list.files(pattern = "^RF_pred_bootstrap_\\d{3}\\.tif$")))

# ===================== everage and get the sd =====================
library(raster)
library(tidyverse)

# setwd("your/path/here")

pred_files <- list.files(pattern = "^RF_pred_bootstrap_\\d{3}\\.tif$")
cat("Found", length(pred_files), "prediction files\n")

pred_files <- pred_files[order(as.numeric(gsub("RF_pred_bootstrap_(\\d+)\\.tif", "\\1", pred_files)))]

if (length(pred_files) == 0) {
  stop("No prediction files found!")
}

template <- raster(pred_files[1])
n_layers <- length(pred_files)

cat("Loading prediction rasters...\n")
pred_stack <- stack()

for (i in seq_along(pred_files)) {
  cat("Loading:", pred_files[i], "(", i, "of", length(pred_files), ")\n")
  pred_stack <- addLayer(pred_stack, raster(pred_files[i]))
}

cat("Calculating mean prediction...\n")
mean_pred <- calc(pred_stack, fun = mean, na.rm = TRUE)

cat("Calculating standard deviation...\n")
sd_pred <- calc(pred_stack, fun = sd, na.rm = TRUE)

cat("Saving results...\n")

# save mean
writeRaster(mean_pred, 
           filename = "RF_pred_mean.tif",
           format = "GTiff",
           options = c("COMPRESS=LZW", "PREDICTOR=2"),
           overwrite = TRUE)

# save sd
writeRaster(sd_pred,
           filename = "RF_pred_sd.tif", 
           format = "GTiff",
           options = c("COMPRESS=LZW", "PREDICTOR=2"),
           overwrite = TRUE)


rm(pred_stack)
gc()

cat("All done! Created:\n")
cat("- RF_pred_mean.tif (平均预测)\n") 
cat("- RF_pred_sd.tif (预测标准差)\n")

