#' Extract Annual NLCD land cover percent and area for polygons
#'
#' Returns the percent and area of each Annual NLCD (Collection 1.2) land
#' cover class inside every feature of `geometry`.
#'
#' Speed comes from four things:
#' * Only native 30 m tiles touching your features are requested from the
#'   USGS WCS. Nothing national is downloaded.
#' * Tiles are fetched in parallel (HTTP/2 multiplexing via `curl`) and cached
#'   on disk, so repeat calls and other years over the same area are instant
#'   or near instant.
#' * Requests are snapped to the NLCD grid, so no resampling happens.
#' * Extraction uses exact pixel coverage fractions (`exactextractr`, C++).
#'
#' Memory stays flat for large layers: features are sorted spatially and run
#' in chunks of `chunk_size`, each chunk reads only its own tiles through a
#' lazy VRT, and per-feature cell reads are capped.
#'
#' @param year Integer, 1985 to 2025.
#' @param geometry An `sf`/`sfc` object, a `terra::SpatVector`, or a path to
#'   any vector file `sf::st_read()` can open. Polygons only. Must fall within
#'   CONUS.
#' @param units Area units: one of `"ha"`, `"km2"`, `"m2"`, `"acres"`, `"mi2"`.
#' @param id Optional column name in `geometry` to use as the feature id.
#'   Defaults to row number.
#' @param wide If `TRUE`, return one row per feature with `pct_<class>` and
#'   `area_<class>` columns (class colors in `attr(x, "nlcd_colors")`). Default `FALSE` returns long (tidy) format.
#' @param cache_dir Directory for cached tiles. Defaults to a persistent user
#'   cache (`tools::R_user_dir("fastnlcd", "cache")`).
#' @param chunk_size Features per chunk. Lower it if memory is tight.
#' @param cores Worker processes for extraction (forked; Unix/macOS only,
#'   ignored on Windows). Tiles are always downloaded in parallel.
#' @param quiet Suppress progress messages.
#'
#' @return A `data.frame`. Long format columns: `id`, `year`, `class`,
#'   `description`, `nlcd_color` (official NLCD hex color), `area` (in `units`), `pct` (0-100, of the feature's
#'   area with valid NLCD data). Classes absent from a feature are omitted.
#' @export
#' @examples
#' \dontrun{
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"))
#' extract_nlcd(year = 2025, geometry = nc[1:5, ], units = "km2", id = "NAME")
#' }
extract_nlcd <- function(year, geometry, units = c("ha", "km2", "m2", "acres", "mi2"),
                         id = NULL, wide = FALSE, cache_dir = NULL,
                         chunk_size = 500L, cores = 1L, quiet = FALSE) {
  units <- match.arg(units)
  year <- as.integer(year)
  if (length(year) != 1L || is.na(year) || year < 1985L || year > 2025L) {
    stop("`year` must be a single year from 1985 to 2025.", call. = FALSE)
  }
  if (is.null(cache_dir)) cache_dir <- tools::R_user_dir("fastnlcd", "cache")
  ydir <- file.path(cache_dir, year)
  dir.create(ydir, recursive = TRUE, showWarnings = FALSE)

  g <- read_geom(geometry)
  ids <- if (is.null(id)) seq_len(nrow(g)) else {
    if (!id %in% names(g)) stop("`id` column '", id, "' not found.", call. = FALSE)
    g[[id]]
  }
  g <- sf::st_transform(sf::st_geometry(g), 5070)
  n <- length(g)
  if (n == 0L) stop("`geometry` has no features.", call. = FALSE)

  bb <- feature_bboxes(g)
  ok <- stats::complete.cases(bb)
  ord <- order(bb[, "row0"], bb[, "col0"])
  ord <- ord[ok[ord]]
  chunks <- split(ord, ceiling(seq_along(ord) / chunk_size))

  cell_m2 <- 900
  conv <- c(m2 = 1, ha = 1e-4, km2 = 1e-6, acres = 1 / 4046.8564224, mi2 = 1 / 2589988.110336)[[units]]
  out <- vector("list", length(chunks))
  t0 <- Sys.time()

  all_tiles <- tiles_for(bb[ord, , drop = FALSE])
  if (!quiet) message("[extract_nlcd] ", length(all_tiles), " tiles (",
                      sum(!file.exists(file.path(ydir, paste0(all_tiles, ".tif")))),
                      " to download)")
  fetch_tiles(all_tiles, year, ydir)

  run_chunk <- function(k) {
    idx <- chunks[[k]]
    files <- file.path(ydir, paste0(tiles_for(bb[idx, , drop = FALSE]), ".tif"))
    vrt <- tempfile(fileext = ".vrt")
    on.exit(unlink(vrt))
    r <- terra::vrt(files, filename = vrt)
    terra::NAflag(r) <- NLCD_NODATA
    res <- exactextractr::exact_extract(
      r, sf::st_sf(geometry = g[idx]), fun = c("frac", "count"),
      max_cells_in_memory = 3e7, progress = FALSE
    )
    to_long(res, ids[idx], year, cell_m2 * conv)
  }

  cores <- max(1L, as.integer(cores))
  if (cores > 1L && .Platform$OS.type == "unix") {
    out <- parallel::mclapply(seq_along(chunks), run_chunk, mc.cores = cores,
                              mc.preschedule = FALSE)
    err <- vapply(out, inherits, logical(1), "try-error")
    if (any(err)) stop(out[[which(err)[1]]], call. = FALSE)
  } else {
    for (k in seq_along(chunks)) {
      out[[k]] <- run_chunk(k)
      if (!quiet) {
        message(sprintf("\r[extract_nlcd] %d/%d features | %.1fs",
                        sum(lengths(chunks[seq_len(k)])), length(ord),
                        as.numeric(Sys.time() - t0, units = "secs")), appendLF = FALSE)
      }
    }
  }
  if (!quiet) message("")
  if (any(!ok)) warning(sum(!ok), " empty geometries skipped.", call. = FALSE)

  df <- do.call(rbind, out)
  df <- df[order(match(df$id, ids), df$class), , drop = FALSE]
  rownames(df) <- NULL
  if (wide) {
    pal <- unique(df[c("class", "nlcd_color")])
    df <- to_wide(df, ids[ok])
    attr(df, "nlcd_colors") <- stats::setNames(pal$nlcd_color, pal$class)
  }
  attr(df, "units") <- units
  df
}
