test_that("input validation", {
  sq <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0),c(1,0),c(1,1),c(0,0)))), crs = 5070)
  expect_error(extract_nlcd(1984, sq), "1985")
  expect_error(extract_nlcd(2026, sq), "1985")
  expect_error(extract_nlcd(2020, sf::st_sfc(sf::st_point(c(0, 0)), crs = 5070)), "POLYGON")
})

test_that("tile keys cover bbox", {
  g <- sf::st_sfc(sf::st_buffer(sf::st_point(c(1e6, 1.8e6)), 40000), crs = 5070)
  bb <- fastnlcd:::feature_bboxes(g)
  expect_true(length(fastnlcd:::tiles_for(bb)) >= 4)
})

test_that("live extraction sums to 100%", {
  skip_on_cran()
  skip_if_offline("dmsdata.cr.usgs.gov")
  g <- sf::st_sfc(
    sf::st_buffer(sf::st_point(c(1015000, 1815000)), 2000),
    sf::st_buffer(sf::st_point(c(1050000, 1850000)), 5000),
    crs = 5070)
  x <- extract_nlcd(2025, g, units = "ha", cache_dir = tempfile(), quiet = TRUE)
  s <- tapply(x$pct, x$id, sum)
  expect_equal(as.vector(s), c(100, 100), tolerance = 1e-6)
  a <- tapply(x$area, x$id, sum)
  expect_equal(as.vector(a), as.numeric(sf::st_area(g)) / 1e4, tolerance = 0.01)
  w <- extract_nlcd(2025, g, units = "km2", wide = TRUE, cache_dir = tempfile(), quiet = TRUE)
  expect_equal(nrow(w), 2)
  expect_true(all(grepl("^#[0-9A-F]{6}$", x$nlcd_color)))
})
