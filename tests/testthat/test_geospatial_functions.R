# Perform unit tests on geospatial functions

library(testthat)

source('../../code/geospatial-functions.R')

# Data to test functions with
data_hull <- data.frame(
  genus = 'Angulithes',
  p_lng_PALEOMAP = c(22.4714),
  p_lat_PALEOMAP = c(13.3926)
)

data_grid <- data.frame(
  genus = 'Epicymatoceras',
  p_lng_PALEOMAP = c(46.1666, 44.9509, 45.9943,  8.1176,  8.3009),
  p_lat_PALEOMAP = c(36.6893, 37.4041, 36.6062, 49.2477, 48.9229)
)

# Here are the tests

test_that(
  'Convex hulls are being drawn correctly', {
  expect_equal(
    round(draw.convex.hull(data_hull, 'p_lat_PALEOMAP', 'p_lng_PALEOMAP',
      buffer = TRUE), digits = 4), 
    312.8689)
   expect_type(draw.convex.hull(data_hull, 'p_lat_PALEOMAP', 'p_lng_PALEOMAP', 
       buffer = TRUE, return.hull = TRUE), 
     'list')
  }
)

test_that(
  'Grid sub-sampling works', {
  expect_setequal(
    unique.within.grid(data_grid, 'p_lat_PALEOMAP', 'p_lng_PALEOMAP',
      'EPSG:4326', 'EPSG:8857', c(25000) )$cell,
    c(8429, 7981, 8577, 303, 452))
  }
)

test_that(
  'Jackknifing works', {
  expect_setequal(
    round(jackknife.area(data_grid, 'p_lat_PALEOMAP', 'p_lng_PALEOMAP', 
      'EPSG:4326'), digits = 2),
    c(83998.39, 85266.79, 56893.29, 30442.10, 28357.46)
  )
  withr::local_seed(5)
  expect_setequal(
    round(jackknife.area(data_grid, 'p_lat_PALEOMAP', 'p_lng_PALEOMAP',
      'EPSG:4326', maxiter = 4), digits = 4),
    c(28513.5528, 156.0908, 27245.1564, 1112.3056)
  )
  }
)

test_that(
  'Boostrapping works', {
  withr::local_seed(5)
  expect_setequal(
    round(bootstrap.area(data_grid, 'p_lat_PALEOMAP', 'p_lng_PALEOMAP',
      'EPSG:4326'), digits = 4),
    c(1112.3056, 156.0908, 2224.6884, 28513.5528, 28217.4122, 57049.3759,
      27245.1564, 29329.7950, 54668.5966, 56753.2353)
  )
  withr::local_seed(5)
  expect_setequal(
      round(bootstrap.area(data_grid, 'p_lat_PALEOMAP', 'p_lng_PALEOMAP',
        'EPSG:4326', maxiter = 4), digits = 4),
      c(156.0908, 54668.5966, 27245.1564, 2224.6884)
  )
  }
)