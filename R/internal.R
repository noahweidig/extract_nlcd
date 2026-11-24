NLCD_WCS <- "https://dmsdata.cr.usgs.gov/geoserver/mrlc_Land-Cover-Native_conus_year_data/wcs"
NLCD_COVERAGE <- "Land-Cover-Native_conus_year_data"
NLCD_NODATA <- 250
# Annual NLCD CONUS grid origin (EPSG:5070, upper-left corner)
GRID_X0 <- -2415585
GRID_Y0 <- 3314805
RES <- 30
TILE_PX <- 1024L
TILE_M <- TILE_PX * RES

NLCD_CLASSES <- c(
  "11" = "Open Water", "12" = "Perennial Ice/Snow",
  "21" = "Developed, Open Space", "22" = "Developed, Low Intensity",
  "23" = "Developed, Medium Intensity", "24" = "Developed, High Intensity",
  "31" = "Barren Land", "41" = "Deciduous Forest", "42" = "Evergreen Forest",
  "43" = "Mixed Forest", "52" = "Shrub/Scrub", "71" = "Grassland/Herbaceous",
  "81" = "Pasture/Hay", "82" = "Cultivated Crops", "90" = "Woody Wetlands",
  "95" = "Emergent Herbaceous Wetlands"
)

# Official palette, from the color table embedded in the USGS tiles
NLCD_COLORS <- c(
  "11" = "#466B9F", "12" = "#D1DEF8", "21" = "#DEC5C5", "22" = "#D99282",
  "23" = "#EB0000", "24" = "#AB0000", "31" = "#B3AC9F", "41" = "#68AB5F",
  "42" = "#1C5F2C", "43" = "#B5C58F", "52" = "#CCB879", "71" = "#DFDFC2",
  "81" = "#DCD939", "82" = "#AB6C28", "90" = "#B8D9EB", "95" = "#6C9FB8"
)

read_geom <- function(x) {
  if (is.character(x) && length(x) == 1L) x <- sf::st_read(x, quiet = TRUE)
  if (inherits(x, "SpatVector")) x <- sf::st_as_sf(x)
  if (inherits(x, "sfc")) x <- sf::st_sf(geometry = x)
  if (!inherits(x, "sf")) stop("`geometry` must be sf, sfc, SpatVector, or a file path.", call. = FALSE)
  if (is.na(sf::st_crs(x))) stop("`geometry` has no CRS.", call. = FALSE)
  types <- unique(as.character(sf::st_geometry_type(x)))
  if (!all(types %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("`geometry` must be POLYGON/MULTIPOLYGON.", call. = FALSE)
  }
  x
}

# Per-feature tile index ranges, vectorised (no per-feature R loop).
feature_bboxes <- function(g) {
  b <- vapply(g, function(p) {
    if (length(p) == 0L) return(rep(NA_real_, 4))
    as.numeric(sf::st_bbox(p))
  }, numeric(4))
  col0 <- floor((b[1, ] - GRID_X0) / TILE_M)
  col1 <- floor((b[3, ] - GRID_X0) / TILE_M)
  row0 <- floor((GRID_Y0 - b[4, ]) / TILE_M)
  row1 <- floor((GRID_Y0 - b[2, ]) / TILE_M)
  cbind(col0 = col0, col1 = col1, row0 = row0, row1 = row1)
}

tiles_for <- function(bb) {
  keys <- unlist(lapply(seq_len(nrow(bb)), function(i) {
    cc <- bb[i, "col0"]:bb[i, "col1"]
    rr <- bb[i, "row0"]:bb[i, "row1"]
    paste(rep(cc, times = length(rr)), rep(rr, each = length(cc)), sep = "_")
  }), use.names = FALSE)
  unique(keys)
}

tile_url <- function(key, year) {
  cr <- as.numeric(do.call(rbind, strsplit(key, "_", fixed = TRUE)))
  cr <- matrix(cr, ncol = 2)
  xmin <- GRID_X0 + cr[, 1] * TILE_M
  ymax <- GRID_Y0 - cr[, 2] * TILE_M
  sprintf(paste0(
    "%s?service=WCS&version=1.0.0&request=GetCoverage&coverage=%s",
    "&crs=EPSG:5070&response_crs=EPSG:5070&bbox=%.0f,%.0f,%.0f,%.0f",
    "&width=%d&height=%d&format=GeoTIFF&time=%d-01-01T00:00:00Z"
  ), NLCD_WCS, NLCD_COVERAGE, xmin, ymax - TILE_M, xmin + TILE_M, ymax,
  TILE_PX, TILE_PX, year)
}

is_tiff <- function(f) {
  if (!file.exists(f) || file.size(f) < 8) return(FALSE)
  m <- readBin(f, "raw", 4)
  identical(m, as.raw(c(0x49, 0x49, 0x2a, 0x00))) || identical(m, as.raw(c(0x4d, 0x4d, 0x00, 0x2a)))
}

fetch_tiles <- function(keys, year, dir, batch = 100L, tries = 6L) {
  files <- file.path(dir, paste0(keys, ".tif"))
  for (attempt in seq_len(tries)) {
    need <- which(!vapply(files, is_tiff, logical(1)))
    if (!length(need)) return(files)
    for (b in split(need, ceiling(seq_along(need) / batch))) {
      tmp <- paste0(files[b], ".part")
      res <- tryCatch(curl::multi_download(
        tile_url(keys[b], year), tmp, progress = FALSE,
        multiplex = TRUE, timeout = 300,
        useragent = "fastnlcd R package"
      ), error = function(e) NULL)
      ok <- if (is.null(res)) rep(FALSE, length(b)) else
        res$success %in% TRUE & res$status_code %in% 200L
      good <- ok & vapply(tmp, is_tiff, logical(1))
      file.rename(tmp[good], files[b][good])
      unlink(tmp[!good])
      if (!all(good)) Sys.sleep(min(60, 2^attempt))
    }
  }
  bad <- !vapply(files, is_tiff, logical(1))
  if (any(bad)) {
    stop(sum(bad), " NLCD tile(s) failed to download from USGS WCS after ", tries,
         " tries. Re-run to resume (downloaded tiles are cached).", call. = FALSE)
  }
  files
}

to_long <- function(res, ids, year, cell_area) {
  fcols <- grep("^frac_", names(res), value = TRUE)
  if (!length(fcols)) {
    return(data.frame(id = ids[0], year = integer(), class = integer(),
                      description = character(), nlcd_color = character(),
                      area = numeric(), pct = numeric()))
  }
  cls <- as.integer(sub("^frac_", "", fcols))
  fr <- as.matrix(res[fcols])
  fr[is.na(fr)] <- 0
  keep <- cls %in% as.integer(names(NLCD_CLASSES))
  fr <- fr[, keep, drop = FALSE]
  cls <- cls[keep]
  nz <- which(fr > 0, arr.ind = TRUE)
  valid_area <- res$count * cell_area
  data.frame(
    id = ids[nz[, 1]],
    year = year,
    class = cls[nz[, 2]],
    description = unname(NLCD_CLASSES[as.character(cls[nz[, 2]])]),
    nlcd_color = unname(NLCD_COLORS[as.character(cls[nz[, 2]])]),
    area = fr[nz] * valid_area[nz[, 1]],
    pct = fr[nz] * 100,
    stringsAsFactors = FALSE
  )
}

to_wide <- function(df, ids) {
  cls <- sort(unique(df$class))
  w <- data.frame(id = ids, year = if (nrow(df)) df$year[1] else NA_integer_)
  i <- match(df$id, ids)
  for (cl in cls) {
    s <- df$class == cl
    p <- a <- numeric(length(ids))
    p[i[s]] <- df$pct[s]
    a[i[s]] <- df$area[s]
    w[[paste0("pct_", cl)]] <- p
    w[[paste0("area_", cl)]] <- a
  }
  w
}
