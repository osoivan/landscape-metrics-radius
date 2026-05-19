# ============================================================
# ONLY MISSING METRICS: ENN_MN + PARA_MN
# Robust PSOCK version
# Reuses the same strategy as your working code
# ============================================================

library(terra)
library(landscapemetrics)
library(data.table)
library(parallel)

# -----------------------------
# USER SETTINGS
# -----------------------------
in_lc   <- "C:/Zoomal/Landcover/test_uni.tif"
out_dir <- "C:/Zoomal/test_uni"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

crs_out          <- "EPSG:3857"
radius_m         <- 500
out_res_m        <- 100
agg_res_m        <- 100
chunk_size       <- 500     # lower because ENN/PARA are heavier
n_workers        <- 3       # safer for RAM
save_chunks      <- TRUE

exclude_clouds <- TRUE
cloud_code     <- 7

aoi_path <- NULL
# aoi_path <- "D:/PK_MALARIA/aoi/study_area.shp"

classes <- c(
  forest                  = 1,
  plantations             = 2,
  grass_shrub_agriculture = 3,
  urban_built_up          = 4,
  bare_open_land          = 5,
  water_wetland           = 6,
  clouds                  = 7
)

metrics_missing <- c(
  "lsm_c_enn_mn",
  "lsm_c_para_mn"
)

gdal_opts <- c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER")

terraOptions(progress = 1, memfrac = 0.85, tempdir = tempdir())

# -----------------------------
# HELPER FUNCTIONS
# -----------------------------
safe_modal <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

write_metric_raster <- function(dt_sub, out_path, template_rast, aoi = NULL, gdal_opts) {
  r <- template_rast
  v <- rep(NA_real_, ncell(r))
  v[dt_sub$cell] <- dt_sub$value
  values(r) <- v
  
  if (!is.null(aoi)) {
    r <- mask(r, aoi)
  }
  
  r[is.na(r)] <- 0
  
  writeRaster(
    r,
    out_path,
    overwrite = TRUE,
    NAflag = -9999,
    gdal = gdal_opts
  )
}

run_chunk_worker <- function(i, chunk_ids, raster_file, pts_df, radius_m, metric_set,
                             out_dir, save_chunks) {
  library(terra)
  library(landscapemetrics)
  library(data.table)
  
  idx <- chunk_ids[[i]]
  sub_pts <- pts_df[idx, , drop = FALSE]
  
  landscape <- rast(raster_file)
  coords <- as.matrix(sub_pts[, c("x", "y")])
  
  smp <- suppressWarnings(
    sample_lsm(
      landscape   = landscape,
      y           = coords,
      plot_id     = sub_pts$cell,
      shape       = "circle",
      size        = radius_m,
      transform   = FALSE,
      all_classes = TRUE,
      progress    = FALSE,
      what        = metric_set,
      level       = "class"
    )
  )
  
  if (is.null(smp) || nrow(smp) == 0) {
    message("Finished chunk ", i, " of ", length(chunk_ids), " (empty)")
    return(NULL)
  }
  
  dt <- as.data.table(smp)[, .(
    cell         = as.integer(plot_id),
    class        = as.integer(class),
    metric_short = as.character(metric),
    value        = as.numeric(value)
  )]
  
  if (save_chunks && nrow(dt) > 0) {
    fwrite(dt, file.path(out_dir, paste0("chunk_", i, ".csv")))
  }
  
  message("Finished chunk ", i, " of ", length(chunk_ids))
  dt
}

run_metric_set_parallel <- function(metric_set, tag, chunk_ids, raster_file, pts_df, radius_m,
                                    out_dir, save_chunks, n_workers,
                                    exclude_clouds = TRUE, cloud_code = 7) {
  cat("\nRunning metric set:", tag, "\n")
  
  tag_dir <- file.path(out_dir, tag)
  dir.create(tag_dir, showWarnings = FALSE)
  
  cl <- makeCluster(n_workers, type = "PSOCK")
  on.exit(stopCluster(cl), add = TRUE)
  
  clusterEvalQ(cl, {
    library(terra)
    library(landscapemetrics)
    library(data.table)
    NULL
  })
  
  clusterExport(
    cl,
    varlist = c(
      "chunk_ids", "raster_file", "pts_df", "radius_m",
      "metric_set", "tag_dir", "save_chunks", "run_chunk_worker"
    ),
    envir = environment()
  )
  
  res_list <- parLapplyLB(
    cl,
    X = seq_along(chunk_ids),
    fun = function(i) {
      run_chunk_worker(
        i = i,
        chunk_ids = chunk_ids,
        raster_file = raster_file,
        pts_df = pts_df,
        radius_m = radius_m,
        metric_set = metric_set,
        out_dir = tag_dir,
        save_chunks = save_chunks
      )
    }
  )
  
  res_dt <- rbindlist(res_list, fill = TRUE)
  
  if (nrow(res_dt) == 0) return(res_dt)
  
  if (exclude_clouds) {
    res_dt <- res_dt[class != cloud_code]
  }
  
  fwrite(res_dt, file.path(out_dir, paste0("all_", tag, "_metrics_long.csv")))
  res_dt
}

# -----------------------------
# 1) LOAD
# -----------------------------
cat("\n[1/8] Loading raster...\n")
lu <- rast(in_lc)

# -----------------------------
# 2) PROJECT
# -----------------------------
cat("[2/8] Projecting raster...\n")
lu <- project(lu, crs_out, method = "near", by_util = TRUE)
lu <- round(lu)

if (exclude_clouds) {
  lu[lu == cloud_code] <- NA
}

# -----------------------------
# 3) AOI
# -----------------------------
cat("[3/8] Preparing AOI...\n")
aoi <- NULL
if (!is.null(aoi_path)) {
  aoi <- vect(aoi_path)
  aoi <- project(aoi, crs_out)
  aoi <- makeValid(aoi)
  lu <- crop(lu, aoi)
  lu <- mask(lu, aoi)
}

# -----------------------------
# 4) AGGREGATE
# -----------------------------
cat("[4/8] Aggregating raster...\n")
res_in <- res(lu)[1]
fact   <- max(1, round(agg_res_m / res_in))

cat("Original resolution:", res_in, "m\n")
cat("Aggregation factor :", fact, "\n")
cat("Approx agg. res.   :", fact * res_in, "m\n")

lu_fast <- try(
  aggregate(lu, fact = fact, fun = "modal"),
  silent = TRUE
)

if (inherits(lu_fast, "try-error")) {
  cat("Falling back to safe custom modal function...\n")
  lu_fast <- aggregate(lu, fact = fact, fun = safe_modal)
}

lu_fast_file <- file.path(out_dir, "lu_fast_agg.tif")
writeRaster(
  lu_fast,
  lu_fast_file,
  overwrite = TRUE,
  gdal = gdal_opts
)

rm(lu_fast)
gc()

# -----------------------------
# 5) TEMPLATE + POINTS
# -----------------------------
cat("[5/8] Creating 1 km template...\n")
lu_fast_tmp <- rast(lu_fast_file)

tmpl <- rast(ext(lu_fast_tmp), resolution = out_res_m, crs = crs_out)
values(tmpl) <- NA_real_

if (!is.null(aoi)) {
  aoi_mask <- rasterize(aoi, tmpl, field = 1, background = NA)
  keep_cells <- which(!is.na(values(aoi_mask)))
  rm(aoi_mask)
} else {
  vals_fast <- !is.na(resample(lu_fast_tmp, tmpl, method = "near"))
  keep_cells <- which(values(vals_fast) == 1)
  rm(vals_fast)
}

xy <- xyFromCell(tmpl, keep_cells)

pts_df <- data.frame(
  x = xy[, 1],
  y = xy[, 2],
  cell = keep_cells
)

rm(xy)
gc()

cat("Number of 1 km cells to process:", nrow(pts_df), "\n")

# -----------------------------
# 6) CHUNKS
# -----------------------------
cat("[6/8] Building chunks...\n")
chunk_ids <- split(
  seq_len(nrow(pts_df)),
  ceiling(seq_len(nrow(pts_df)) / chunk_size)
)
cat("Number of chunks:", length(chunk_ids), "\n")

# -----------------------------
# 7) RUN ONLY MISSING METRICS
# -----------------------------
cat("[7/8] Running missing metrics...\n")

res_missing <- run_metric_set_parallel(
  metric_set      = metrics_missing,
  tag             = "missing",
  chunk_ids       = chunk_ids,
  raster_file     = lu_fast_file,
  pts_df          = pts_df,
  radius_m        = radius_m,
  out_dir         = out_dir,
  save_chunks     = save_chunks,
  n_workers       = n_workers,
  exclude_clouds  = exclude_clouds,
  cloud_code      = cloud_code
)

if (nrow(res_missing) == 0) {
  stop("No missing metrics were computed.")
}

fwrite(res_missing, file.path(out_dir, "missing_metrics_long.csv"))

# -----------------------------
# 8) EXPORT FINAL RASTERS
# -----------------------------
cat("[8/8] Exporting rasters...\n")
metric_names <- sort(unique(res_missing$metric_short))

for (m in metric_names) {
  for (cls_name in names(classes)) {
    cls_code <- classes[[cls_name]]
    if (exclude_clouds && cls_code == cloud_code) next
    
    dt_sub <- res_missing[metric_short == m & class == cls_code, .(cell, value)]
    out_file <- file.path(out_dir, paste0(m, "_", cls_name, ".tif"))
    
    if (nrow(dt_sub) == 0) {
      r0 <- tmpl
      if (!is.null(aoi)) r0 <- mask(r0, aoi)
      r0[is.na(r0)] <- 0
      writeRaster(
        r0,
        out_file,
        overwrite = TRUE,
        NAflag = -9999,
        gdal = gdal_opts
      )
      cat("Wrote (zeros):", out_file, "\n")
    } else {
      write_metric_raster(
        dt_sub        = dt_sub,
        out_path      = out_file,
        template_rast = tmpl,
        aoi           = aoi,
        gdal_opts     = gdal_opts
      )
      cat("Wrote:", out_file, "\n")
    }
  }
}

cat("\nDONE ✅ Missing metrics exported to:\n", out_dir, "\n", sep = "")