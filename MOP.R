
# =========================================
# MOP via terra (tiled; standardized by M)
# Paths anonymized with placeholders
# =========================================

rm(list = ls()); gc()

# ---- Dependencies ----
pkgs <- c("terra","sf","units","fields")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(terra); library(sf); library(units); library(fields)

# ---- Paths (ANONYMIZED) ----
# Replace the placeholders with your actual paths before running.
var_dir      <- "<PATH_TO_VARS_DIR>"            # e.g., "D:/project/mop/vars"
presence_csv <- "<PATH_TO_PRESENCE_CSV>"                # CSV with cols: XX, YY (WGS84)
out_dir      <- "<PATH_TO_OUTPUT_DIR>"                  # e.g., "D:/project/mop/outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Parameters ----
BUFFER_KM   <- 30
M_MAX       <- 30000
TILES_X     <- 8
TILES_Y     <- 8
G_BATCH     <- 2000
PERC_VEC    <- c(1,10,50,100)
MASK_STRICT <- FALSE

# ---- Helper: compute MOP for a batch of G (z-scored) ----
mop_batch <- function(M_vars_z, G_vars_z, perc_vec, g_batch = 2000) {
  nG <- nrow(G_vars_z)
  out <- matrix(NA_real_, nrow = nG, ncol = length(perc_vec))
  colnames(out) <- paste0("p", sprintf("%02d", perc_vec))
  i <- 1
  while (i <= nG) {
    j <- min(i + g_batch - 1, nG)
    D <- rdist(M_vars_z, G_vars_z[i:j, , drop = FALSE])
    for (cc in seq_len(ncol(D))) {
      di  <- D[, cc]
      ord <- sort(di, na.last = NA)
      nM  <- length(ord)
      for (k in seq_along(perc_vec)) {
        kk <- max(1L, floor(perc_vec[k] / 100 * nM))
        out[i + cc - 1, k] <- mean(ord[1:kk], na.rm = TRUE)
      }
    }
    i <- j + 1
  }
  as.data.frame(out)
}

# ---- Helper: write tile results back into rasters (4 layers) ----
write_tile <- function(out_rast, tile_ext, xy_df, vals_df) {
  tile_template <- crop(out_rast[[1]], tile_ext)
  dat <- cbind(x = xy_df$x, y = xy_df$y, vals_df)
  pts <- vect(dat, geom = c("x","y"), crs = crs(out_rast))
  nmz <- names(vals_df)  # "p01","p10","p50","p100"
  for (k in seq_along(nmz)) {
    r_tile <- rasterize(
      pts, tile_template, field = nmz[k],
      fun = "mean", background = NA, touches = FALSE
    )
    out_rast[[k]] <- mosaic(out_rast[[k]], r_tile, fun = "mean")
  }
  out_rast
}

# ---- Read variables; build M (buffer) and G (full) ----
var_files <- list.files(var_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(var_files) >= 3)
r <- rast(var_files)
crs_r <- crs(r)

pts <- read.csv(presence_csv)
stopifnot(all(c("XX","YY") %in% names(pts)))
pts_sf   <- st_as_sf(pts, coords = c("XX","YY"), crs = 4326)
pts_3857 <- st_transform(pts_sf, 3857)
buf      <- st_union(st_buffer(pts_3857, set_units(BUFFER_KM, "km")))
buf_rcrs <- st_transform(buf, crs = st_crs(crs_r))
buf_v    <- vect(buf_rcrs)

r_M <- mask(crop(r, buf_v), buf_v)  # reference region (M)
r_G <- r                            # projection region (G)

# ---- Stats on M (for standardization and strict extrapolation) ----
mu_df  <- global(r_M, mean, na.rm = TRUE)
sd_df  <- global(r_M, sd,   na.rm = TRUE)
min_df <- global(r_M, min,  na.rm = TRUE)
max_df <- global(r_M, max,  na.rm = TRUE)

mu_vec  <- as.numeric(mu_df$mean)
sd_vec  <- as.numeric(sd_df$sd)
min_vec <- as.numeric(min_df$min)
max_vec <- as.numeric(max_df$max)

sd_vec[!is.finite(sd_vec) | sd_vec == 0] <- 1e-9

stopifnot(
  length(mu_vec)  == nlyr(r),
  length(sd_vec)  == nlyr(r),
  length(min_vec) == nlyr(r),
  length(max_vec) == nlyr(r)
)

cat(sprintf("Layers: %d\n", nlyr(r)))
cat("M-based mean/sd/min/max computed.\n")

# ---- Reference cloud M: sampling + standardization ----
set.seed(1)
M_samp <- spatSample(
  r_M, size = M_MAX, method = "random",
  as.points = TRUE, values = TRUE, na.rm = TRUE
)
M_df     <- as.data.frame(M_samp, geom = "XY")
M_vars   <- as.matrix(M_df[, names(r_M), drop = FALSE])
M_vars_z <- sweep(sweep(M_vars, 2, mu_vec, "-"), 2, sd_vec, "/")

cat(sprintf("Reference cloud size (M): %d\n", nrow(M_vars_z)))

# ---- Prepare outputs: raw distances (4 layers) & strict mask (1 layer) ----
out_raw <- rast(r_G[[1]])
out_raw <- c(out_raw, out_raw, out_raw, out_raw)
names(out_raw) <- c("MOPraw_1pct","MOPraw_10pct","MOPraw_50pct","MOPraw_100pct")
values(out_raw) <- NA

out_strict <- rast(r_G[[1]])
names(out_strict) <- "MOP_strict_extrap"
values(out_strict) <- NA

# ---- Tiled processing over G ----
extG <- ext(r_G)
xs <- seq(extG$xmin, extG$xmax, length.out = TILES_X + 1)
ys <- seq(extG$ymin, extG$ymax, length.out = TILES_Y + 1)

cat(sprintf("Tiling: %dx%d; G batch=%d; M_MAX=%d\n", TILES_X, TILES_Y, G_BATCH, M_MAX))

for (ix in 1:TILES_X) {
  for (iy in 1:TILES_Y) {
    tile_ext <- ext(xs[ix], xs[ix+1], ys[iy], ys[iy+1])
    r_tile   <- crop(r_G, tile_ext)
    if (ncell(r_tile) == 0) next
    
    pts_tile <- as.points(r_tile, values = TRUE, na.rm = TRUE)
    if (nrow(pts_tile) == 0) next
    
    df_tile <- as.data.frame(pts_tile, geom = "XY")
    G_vars  <- as.matrix(df_tile[, names(r_G), drop = FALSE])
    xy      <- df_tile[, c("x","y")]
    
    below <- t(t(G_vars) < min_vec)
    above <- t(t(G_vars) > max_vec)
    strict_num <- as.integer(rowSums(below | above) > 0)
    
    G_vars_z <- sweep(sweep(G_vars, 2, mu_vec, "-"), 2, sd_vec, "/")
    
    vals_df <- mop_batch(M_vars_z, G_vars_z, PERC_VEC, g_batch = G_BATCH)
    colnames(vals_df) <- c("p01","p10","p50","p100")
    
    out_raw <- write_tile(out_raw, tile_ext, xy, vals_df)
    
    dat_se <- data.frame(x = xy$x, y = xy$y, strict = strict_num)
    pts_se <- vect(dat_se, geom = c("x","y"), crs = crs(out_strict))
    r_se   <- rasterize(
      pts_se, crop(out_strict, tile_ext), field = "strict",
      fun = "mean", background = NA, touches = FALSE
    )
    out_strict <- mosaic(out_strict, r_se, fun = "mean")
    
    cat(sprintf("  - Tile (%d,%d): n=%d; strict=%.2f%%\n",
                ix, iy, nrow(df_tile), 100 * mean(strict_num == 1)))
    gc()
  }
}

# ---- Write raw distances ----
f_raw <- c(
  file.path(out_dir, "MOPraw_1pct.tif"),
  file.path(out_dir, "MOPraw_10pct.tif"),
  file.path(out_dir, "MOPraw_50pct.tif"),
  file.path(out_dir, "MOPraw_100pct.tif")
)
writeRaster(out_raw, filename = f_raw, overwrite = TRUE, datatype = "FLT4S")

# ---- Normalize (0–1) per layer; similarity (=1 - dist01) ----
gmin  <- global(out_raw, min, na.rm = TRUE)$min
gmax  <- global(out_raw, max, na.rm = TRUE)$max
mdmin <- as.numeric(gmin)
mdmax <- as.numeric(gmax)
rng   <- mdmax - mdmin
rng[rng == 0 | !is.finite(rng)] <- 1e-9

out_normDist <- app(out_raw, fun = function(v) { (v - mdmin) / rng })
out_normDist <- clamp(out_normDist, 0, 1)
names(out_normDist) <- c("MOPdist01_1pct","MOPdist01_10pct","MOPdist01_50pct","MOPdist01_100pct")

out_similarity <- 1 - out_normDist
names(out_similarity) <- c("MOPsim01_1pct","MOPsim01_10pct","MOPsim01_50pct","MOPsim01_100pct")

if (MASK_STRICT) {
  out_normDist   <- mask(out_normDist, out_strict, maskvalues = 1)
  out_similarity <- mask(out_similarity, out_strict, maskvalues = 1)
}

# ---- Write normalized distance, similarity, and strict mask ----
f_dist01 <- c(
  file.path(out_dir, "MOPdist01_1pct.tif"),
  file.path(out_dir, "MOPdist01_10pct.tif"),
  file.path(out_dir, "MOPdist01_50pct.tif"),
  file.path(out_dir, "MOPdist01_100pct.tif")
)
f_sim01  <- c(
  file.path(out_dir, "MOPsim01_1pct.tif"),
  file.path(out_dir, "MOPsim01_10pct.tif"),
  file.path(out_dir, "MOPsim01_50pct.tif"),
  file.path(out_dir, "MOPsim01_100pct.tif")
)
f_strict <- file.path(out_dir, "MOP_strict_extrapolation.tif")

writeRaster(out_normDist,   filename = f_dist01, overwrite = TRUE, datatype = "FLT4S")
writeRaster(out_similarity, filename = f_sim01,  overwrite = TRUE, datatype = "FLT4S")
writeRaster(out_strict,     filename = f_strict, overwrite = TRUE, datatype = "INT1U")

cat("Done. Outputs:\n")
for (f in c(f_raw, f_dist01, f_sim01, f_strict)) cat(" - ", f, "\n")
