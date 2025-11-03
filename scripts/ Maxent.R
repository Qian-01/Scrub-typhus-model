
# =========================================================
# MaxEnt modeling pipeline (paths anonymized; English notes)
# =========================================================

rm(list = ls()); gc()

# ---- Java heap (set BEFORE loading rJava) ----
options(java.parameters = "-Xmx16g")  # adjust if needed

# ---- Dependencies ----
if (!requireNamespace("ENMeval", quietly = TRUE)) install.packages("ENMeval")
if (!requireNamespace("rJava",   quietly = TRUE)) install.packages("rJava")
if (!requireNamespace("raster",  quietly = TRUE)) install.packages("raster")
if (!requireNamespace("dismo",   quietly = TRUE)) install.packages("dismo")

library(raster)
library(dismo)
library(ENMeval)
library(rJava)

# ---- Anonymized paths (REPLACE placeholders) ----
presence_csv <- "<PATH_TO_PRESENCE_CSV>"        # CSV with cols: XX (lon), YY (lat) in WGS84
env_dir      <- "<PATH_TO_ENV_TIF_DIR>"         # directory with .tif environmental layers
out_dir      <- "<PATH_TO_OUTPUT_DIR>"          # output directory for model & rasters

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Read presence data ----
presence <- read.csv(presence_csv)
stopifnot(all(c("XX", "YY") %in% names(presence)))
p <- presence[, c("XX", "YY")]

# ---- Load environmental rasters & align to a reference ----
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(env_files) >= 1)

env_rasters <- lapply(env_files, function(f) {
  r <- raster(f)
  crs(r) <- CRS("+proj=longlat +datum=WGS84")  # ensure WGS84
  r
})

ref_raster <- env_rasters[[1]]

# Resample to match reference (extent/resolution)
env_rasters_aligned <- lapply(env_rasters, function(r) {
  same <- try(compareRaster(ref_raster, r, extent = TRUE, rowcol = TRUE, res = TRUE,
                            stopiffalse = FALSE), silent = TRUE)
  if (inherits(same, "try-error") || isFALSE(same)) {
    r <- resample(r, ref_raster, method = "bilinear")
  }
  r
})

# ---- Mask for valid (non-NA across all layers) area ----
env_stack_temp <- stack(env_rasters_aligned)
mask_r <- calc(env_stack_temp, fun = function(x) ifelse(all(!is.na(x)), 1, NA))

# Apply mask
env_rasters_masked <- lapply(env_rasters_aligned, function(r) raster::mask(r, mask_r))

# Optional: downscale resolution to reduce memory
env_rasters_lowres <- lapply(env_rasters_masked, function(r) {
  aggregate(r, fact = 2, fun = mean, na.rm = TRUE)  # adjust 'fact' as needed
})

# Final predictor stack
env_stack <- stack(env_rasters_lowres)
crs(env_stack) <- CRS("+proj=longlat +datum=WGS84")

# ---- Filter presence points with NA predictors removed ----
p_sp   <- SpatialPoints(p, proj4string = CRS("+proj=longlat +datum=WGS84"))
p_vals <- extract(env_stack, p_sp)
keep_p <- rowSums(is.na(p_vals)) == 0
p_sp_valid <- p_sp[keep_p, ]
p <- as.data.frame(coordinates(p_sp_valid))
colnames(p) <- c("XX", "YY")

# ---- Background points within valid mask ----
set.seed(123)
bg_mat <- randomPoints(mask_r, n = 10000)  # background in non-NA area
bg_sp  <- SpatialPoints(bg_mat, proj4string = CRS("+proj=longlat +datum=WGS84"))
bg_vals <- extract(env_stack, bg_sp)
keep_bg <- rowSums(is.na(bg_vals)) == 0
bg_sp_valid <- bg_sp[keep_bg, ]
bg <- as.data.frame(coordinates(bg_sp_valid))
colnames(bg) <- c("lon", "lat")

# ---- Train/test split (70/30) ----
set.seed(123)
train_idx    <- sample(seq_len(nrow(p)), size = floor(0.7 * nrow(p)))
train_points <- p[train_idx, ]
test_points  <- p[-train_idx, ]

# ---- MaxEnt (Java-based dismo::maxent) ----
# NOTE: Ensure maxent is available. dismo::maxent() may require maxent.jar in the path.
# You can place maxent.jar in: system.file("java", package="dismo")
me <- maxent(
  x = env_stack,
  p = train_points,
  a = bg,
  removeDuplicates = TRUE,
  args = c(
    "betamultiplier=1.0",
    "responsecurves",
    "jackknife",
    "pictures=TRUE",
    "writeplotdata=TRUE",
    "randomtestpoints=30",
    "maximumiterations=500"
  ),
  path = file.path(out_dir, "maxent_outputs")
)

# ---- Save model ----
saveRDS(me, file = file.path(out_dir, "model.rds"))

# ---- Predict (cloglog output) ----
pred_path <- file.path(out_dir, "prediction_cloglog.tif")
pred <- predict(me, env_stack,
                args = c("outputformat=cloglog"),
                filename = pred_path,
                overwrite = TRUE,
                options = "COMPRESS=DEFLATE")

# ---- Evaluate (test set) ----
eval_test <- evaluate(p = test_points, a = bg, model = me, x = env_stack)
print(eval_test)

# ---- Evaluate (train set) ----
eval_train <- evaluate(p = train_points, a = bg, model = me, x = env_stack)
print(eval_train@auc)

# ---- Optional: 5-fold CV (simple evaluate-based folds) ----
# For rigorous spatial CV, consider ENMeval or blockCV.
eval_cv <- evaluate(p = p, a = bg, model = me, x = env_stack, nfold = 5)
print(eval_cv)

# ---- Quick visualization ----
plot(raster(pred_path), main = "MaxEnt prediction (cloglog)")
points(p, pch = 16, col = "red")  # presence
