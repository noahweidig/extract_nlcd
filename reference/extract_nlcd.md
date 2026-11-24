# Extract Annual NLCD land cover percent and area for polygons

Returns the percent and area of each Annual NLCD (Collection 1.2) land
cover class inside every feature of `geometry`.

## Usage

``` r
extract_nlcd(
  year,
  geometry,
  units = c("ha", "km2", "m2", "acres", "mi2"),
  id = NULL,
  wide = FALSE,
  cache_dir = NULL,
  chunk_size = 500L,
  cores = 1L,
  quiet = FALSE
)
```

## Arguments

- year:

  Integer, 1985 to 2025.

- geometry:

  An `sf`/`sfc` object, a
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html),
  or a path to any vector file
  [`sf::st_read()`](https://r-spatial.github.io/sf/reference/st_read.html)
  can open. Polygons only. Must fall within CONUS.

- units:

  Area units: one of `"ha"`, `"km2"`, `"m2"`, `"acres"`, `"mi2"`.

- id:

  Optional column name in `geometry` to use as the feature id. Defaults
  to row number.

- wide:

  If `TRUE`, return one row per feature with `pct_<class>` and
  `area_<class>` columns (class colors in `attr(x, "nlcd_colors")`).
  Default `FALSE` returns long (tidy) format.

- cache_dir:

  Directory for cached tiles. Defaults to a persistent user cache
  (`tools::R_user_dir("extractnlcd", "cache")`).

- chunk_size:

  Features per chunk. Lower it if memory is tight.

- cores:

  Worker processes for extraction (forked; Unix/macOS only, ignored on
  Windows). Tiles are always downloaded in parallel.

- quiet:

  Suppress progress messages.

## Value

A `data.frame`. Long format columns: `id`, `year`, `class`,
`description`, `nlcd_color` (official NLCD hex color), `area` (in
`units`), `pct` (0-100, of the feature's area with valid NLCD data).
Classes absent from a feature are omitted.

## Details

Speed comes from four things:

- Only native 30 m tiles touching your features are requested from the
  USGS WCS. Nothing national is downloaded.

- Tiles are fetched in parallel (HTTP/2 multiplexing via `curl`) and
  cached on disk, so repeat calls and other years over the same area are
  instant or near instant.

- Requests are snapped to the NLCD grid, so no resampling happens.

- Extraction uses exact pixel coverage fractions (`exactextractr`, C++).

Memory stays flat for large layers: features are sorted spatially and
run in chunks of `chunk_size`, each chunk reads only its own tiles
through a lazy VRT, and per-feature cell reads are capped.

## Examples

``` r
if (FALSE) { # \dontrun{
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"))
extract_nlcd(year = 2025, geometry = nc[1:5, ], units = "km2", id = "NAME")
} # }
```
