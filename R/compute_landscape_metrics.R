#!/usr/bin/env Rscript
# ============================================================
# Landscape Metrics by Moving Circular Window
# Author: Cesar Ivan Alvarez Mendoza
# Project: ZOOMAL / landscape metrics workflow
# Description:
#   Computes class-level landscape metrics from a categorical land-cover raster
#   using circular sampling windows around a regular output grid.
#   Optionally clips/masks the raster using an AOI shapefile.
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(landscapemetrics)
  library(data.table)
  library(parallel)
})

# -----------------------------
# Argument parser: --key value
# -----------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --key value pairs. See README.md.")
  }
  keys <- gsub("^--", "", args[seq(1, length(args), by = 2)])
  vals <- args[seq(2, length(args), by = 2)]
  as.list(stats::setNames(vals, keys))
}

get_arg <- function(args, key, default = NULL, type = "character") {
  value <- args[[key]]
  if (is.null(value) || value == "NULL" || value == "") value <- default
  if (is.null(value)) return(NULL)
  switch(type,
         numeric = as.numeric(value),
         integer = as.integer(value),
         logical = tolower(value) %in% c("true", "t", "1", "yes", "y"),
         character = as.character(value),
         as.character(value))
}

# -----------------------------
# Helper functions
# -----------------------------
safe_modal <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

parse_class_map <- function(class_string) {
  # Example: "forest=1,plantations=2,urban=4,water=6"
  parts <- unlist(strsplit(class_string, ","))
  kv <- strsplit(parts, "=")
  class_values <- as.integer(vapply(kv, `[`, character(1), 2))
  class_names <- vapply(kv, `[`, character(1), 1)
  stats::setNames(class_values, class_names)
}

write_metric_raster <- function(dt_sub, out_path, template_rast, aoi = NULL, gdal_opts) {
  r <- template_rast
  v <- rep(NA_real_, ncell(r))
  v[dt_sub$cell] <- dt_sub$value
  values(r) <- v
  if (!is.null(aoi)) r <- mask(r, aoi)
  r[is.na(r)] <- 0
  writeRaster(r, out_path, overwrite = TRUE, NAflag = -9999, gdal = gdal_opts)
}

run_chunk_worker <- function(i, chunk_ids, raster_file, pts_df, radius_m, metric_set,
                             out_dir, save_chunks) {
  suppressPackageStartupMessages({
    library(terra)
    library(landscapemetrics)
    library(data.table)
  })
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

  if (is.null(smp) || nrow(smp) == 0) return(NULL)

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
                                    exclude_class = NULL) {
  cat("\nRunning metric set:", tag, "\n")
  tag_dir <- file.path(out_dir, tag)
  dir.create(tag_dir, recursive = TRUE, showWarnings = FALSE)

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
    varlist = c("chunk_ids", "raster_file", "pts_df", "radius_m", "metric_set",
                "tag_dir", "save_chunks", "run_chunk_worker"),
    envir = environment()
  )

  res_list <- parLapplyLB(cl, seq_along(chunk_ids), function(i) {
    run_chunk_worker(i, chunk_ids, raster_file, pts_df, radius_m,
                     metric_set, tag_dir, save_chunks)
  })

  res_dt <- rbindlist(res_list, fill = TRUE)
  if (nrow(res_dt) == 0) return(res_dt)
  if (!is.null(exclude_class)) res_dt <- res_dt[!class %in% exclude_class]
  fwrite(res_dt, file.path(out_dir, paste0("all_", tag, "_metrics_long.csv")))
  res_dt
}

# -----------------------------
# Main workflow
# -----------------------------
args <- parse_args()

in_lc      <- get_arg(args, "input", type = "character")
out_dir    <- get_arg(args, "out_dir", default = "outputs", type = "character")
aoi_path   <- get_arg(args, "aoi", default = NULL, type = "character")
crs_out    <- get_arg(args, "crs", default = "EPSG:3857", type = "character")
radius_m   <- get_arg(args, "radius", default = 500, type = "numeric")
out_res_m  <- get_arg(args, "out_res", default = 1000, type = "numeric")
agg_res_m  <- get_arg(args, "agg_res", default = 100, type = "numeric")
chunk_size <- get_arg(args, "chunk_size", default = 500, type = "integer")
n_workers  <- get_arg(args, "workers", default = 3, type = "integer")
save_chunks <- get_arg(args, "save_chunks", default = TRUE, type = "logical")
exclude_codes <- get_arg(args, "exclude_codes", default = NULL, type = "character")
metrics_string <- get_arg(args, "metrics", default = "lsm_c_pland,lsm_c_ed,lsm_c_clumpy,lsm_c_ai,lsm_c_enn_mn,lsm_c_para_mn", type = "character")
classes_string <- get_arg(args, "classes", default = "forest=1,plantations=2,grass_shrub_agriculture=3,urban_built_up=4,bare_open_land=5,water_wetland=6,clouds=7", type = "character")

if (is.null(in_lc)) stop("Missing required argument: --input path/to/landcover.tif")
if (!file.exists(in_lc)) stop("Input raster not found: ", in_lc)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
classes <- parse_class_map(classes_string)
metric_set <- trimws(unlist(strsplit(metrics_string, ",")))
gdal_opts <- c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER")
terraOptions(progress = 1, memfrac = 0.85, tempdir = tempdir())

exclude_class <- NULL
if (!is.null(exclude_codes)) {
  exclude_class <- as.integer(trimws(unlist(strsplit(exclude_codes, ","))))
}

cat("\n[1/8] Loading raster...\n")
lu <- rast(in_lc)

cat("[2/8] Projecting raster...\n")
lu <- project(lu, crs_out, method = "near", by_util = TRUE)
lu <- round(lu)
if (!is.null(exclude_class)) lu[lu %in% exclude_class] <- NA

cat("[3/8] Preparing AOI...\n")
aoi <- NULL
if (!is.null(aoi_path)) {
  if (!file.exists(aoi_path)) stop("AOI file not found: ", aoi_path)
  aoi <- vect(aoi_path)
  aoi <- project(aoi, crs_out)
  aoi <- makeValid(aoi)
  lu <- crop(lu, aoi)
  lu <- mask(lu, aoi)
}

cat("[4/8] Aggregating raster...\n")
res_in <- res(lu)[1]
fact <- max(1, round(agg_res_m / res_in))
cat("Original resolution:", res_in, "m\n")
cat("Aggregation factor:", fact, "\n")
cat("Approx. aggregation resolution:", fact * res_in, "m\n")

lu_fast <- try(aggregate(lu, fact = fact, fun = "modal"), silent = TRUE)
if (inherits(lu_fast, "try-error")) lu_fast <- aggregate(lu, fact = fact, fun = safe_modal)

lu_fast_file <- file.path(out_dir, "landcover_aggregated.tif")
writeRaster(lu_fast, lu_fast_file, overwrite = TRUE, gdal = gdal_opts)
rm(lu_fast); gc()

cat("[5/8] Creating output template and sampling points...\n")
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
pts_df <- data.frame(x = xy[, 1], y = xy[, 2], cell = keep_cells)
rm(xy); gc()
cat("Number of output cells:", nrow(pts_df), "\n")

cat("[6/8] Building chunks...\n")
chunk_ids <- split(seq_len(nrow(pts_df)), ceiling(seq_len(nrow(pts_df)) / chunk_size))
cat("Number of chunks:", length(chunk_ids), "\n")

cat("[7/8] Computing metrics...\n")
res_metrics <- run_metric_set_parallel(
  metric_set = metric_set,
  tag = paste0("r", radius_m, "m"),
  chunk_ids = chunk_ids,
  raster_file = lu_fast_file,
  pts_df = pts_df,
  radius_m = radius_m,
  out_dir = out_dir,
  save_chunks = save_chunks,
  n_workers = n_workers,
  exclude_class = exclude_class
)

if (nrow(res_metrics) == 0) stop("No metrics were computed.")
fwrite(res_metrics, file.path(out_dir, "landscape_metrics_long.csv"))

cat("[8/8] Exporting metric rasters...\n")
metric_names <- sort(unique(res_metrics$metric_short))
for (m in metric_names) {
  for (cls_name in names(classes)) {
    cls_code <- classes[[cls_name]]
    if (!is.null(exclude_class) && cls_code %in% exclude_class) next
    dt_sub <- res_metrics[metric_short == m & class == cls_code, .(cell, value)]
    out_file <- file.path(out_dir, paste0(m, "_", cls_name, ".tif"))
    if (nrow(dt_sub) == 0) {
      r0 <- tmpl
      if (!is.null(aoi)) r0 <- mask(r0, aoi)
      r0[is.na(r0)] <- 0
      writeRaster(r0, out_file, overwrite = TRUE, NAflag = -9999, gdal = gdal_opts)
    } else {
      write_metric_raster(dt_sub, out_file, tmpl, aoi, gdal_opts)
    }
    cat("Wrote:", out_file, "\n")
  }
}

cat("\nDONE. Results exported to: ", out_dir, "\n", sep = "")
