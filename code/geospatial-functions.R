# Collection of functions for geospatial analysis part, including 
# palaeorotations, fitting convex hulls and an equal area grid, sub-sampling
# via bootstrapping and jackknifing

library(terra)
library(sf)
library(bootstrap)

draw.convex.hull <- function(data, collat, collng, proj = 'EPSG:4326',
  return.hull = FALSE, buffer = FALSE, buffer.size = 10000) {
  # Draw convex hulls around point data in dataframe and return its area (km2) 
  # Variables: 
  #   data        = data frame containing points to draw hull around
  #   collat      = column name containing latitude
  #   collng      = column name containing longitude
  #   proj        = projection data is in. Default in R is EPSG:4326
  #   return_hull = return hull as well as area
  #   buffer      = use a buffer around points
  #   buffer.size = size of buffer as radius of circle around points. Default is 
  #                 10000 map units, which is metres for most map projections
  # Returns:
  #   area of hull if return.hull = FALSE
  #   list with area and hull object if return.hull = TRUE
  
  # Convert to Spatvector
  data_vect <- terra::vect(data, geom = c(collng, collat), crs = crs(proj),
    keepgeom = FALSE)
  
  # Create a buffer around locations
  # width should be in meters. TO DO: Check if true for 'EPSG:4326'
  if (buffer == TRUE) {
    data_vect <- terra::buffer(data_vect, width = buffer.size) 
  }

  # Draw convex hull
  data_hull <- terra::hull(data_vect, type = 'convex')
  # is the hull valid?
  if (!is.valid(data_hull)) {
    stop("The convex hull is not valid")
  }
  
  # Get areas of polygons as estimate of geographic distribution area
  distr <- terra::expanse(data_hull, unit = 'km', transform = TRUE)
  
  # Return only the polygon area, or the area and the hull for plotting purposes
  if (return.hull) {
    return( list(distr, data_hull) )
  } else {
    return(distr)
  }
}

unique.within.grid <- function(data, collat, collng, proj = 'EPSG:4326',
  eq.proj = 'EPSG:8857', resolution = 25000, buffer = TRUE, 
  buffer.size = 10000) {
  # Sub-sample fossil occurrences by making occurrence coordinates within the 
  # same equal area grid cell synonymous. Specifically, retain only one 
  # occurrence per cell.
  # Variables: 
  #   data        = data frame containing points/occurrences to fit grid to
  #   collat      = column name containing latitude
  #   collng      = column name containing longitude
  #   proj        = projection data is in. Default in R is EPSG:4326
  #   eq_proj     = projection to fit equal-area grids. Default is EPSG:8857,
  #                 the equal earth projection
  #   resolution  = the resolution to use, specified as the size of grids in 
  #                 meters
  #   buffer      = use a buffer around singleton points
  #   buffer.size = size of buffer as radius of circle around points. Default is 
  #                 10000 map units, which is metres for most map projections
  # Returns:
  #   Data frame of sub-sampled fossil occurrences.
  
  # Convert to Spatvector
  data_vect <- terra::vect(data, geom = c(collng, collat), keepgeom = TRUE, 
    crs = crs(proj))
  
  # Change to equal area projection
  data_proj <- terra::project(data_vect, crs(eq.proj))
  
  # Add a buffer if taxon is only known from a single location. This stops 
  # problems with the data having zero extent
  # TODO Add a test for this feature
  # Check how many unique locations there are for the genus
  uniq.loc <- length(unique(paste(data[ , collat], data[ , collng])))
  if (buffer == TRUE & uniq.loc == 1) {
    data_proj <- terra::buffer(data_proj, width = buffer.size) 
  }
  
  # Make equal area grid cells covering the area with fossil occurrences
  grid <- terra::rast(extent = 1.1 * ext(data_proj), crs = crs(data_proj),
    res = resolution)
  grid[1:ncell(grid)] <- 1:ncell(grid)
  
  # Extract cell ids from raster
  cellids <- terra::extract(grid, data_proj, bind = TRUE, cells = TRUE)
  
  # Use this dataframe for future processing
  cellids_df <- as.data.frame(cellids)
  
  # Retain only a single observation for each cell
  cellids_df <- cellids_df[!duplicated(cellids_df[ , 'cell']), ]
  
  return(cellids_df)
}

jackknife.area <- function(data, collat, collng, proj, maxiter = 1000){
  # Jackknife occurrence data and estimate resulting area by drawing a convex
  # hull around the sub-sampled data. This function always takes away one 
  # observation per sub-sample. 
  # Variables:
  #   data    = data frame containing occurrence data
  #   collat  = column name containing latitude
  #   collng  = column name containing longitude
  #   proj    = projection data is in. Default in R is EPSG:4326
  #   maxiter = maximum number of sub-sampling iterations to perform
  # Returns:
  #   vector of area estimates based on jackknifed data
  
  # Create convenience function for jackknife
  theta <- function(v, data, collat, collng, proj) {
    area <- draw.convex.hull(data[v, ], collat, collng, proj) 
    return(area)
  }
  
  # Limit number of iterations to maxiter
  if (nrow(data) > maxiter) {
    samples <- sample(1:nrow(data), size = maxiter, replace = FALSE)
  } else {
    samples <- 1:nrow(data)
  }
  
  # Perform jackknifing
  results <- jackknife(samples, theta, data, collat, collng, proj)
  return(results$jack.values)
}

bootstrap.area <- function(data, collat, collng, proj, maxiter = 1000,
  n.sample = 3) {
  # Perform boostrap subsampling on occurrence data and estimate resulting
  # area by fitting a convex hull around the sub-sampled data. The size of the
  # sub-sample is set by n.sample.
  # Variables:
  #   data      = data frame containing occurrence data
  #   collat    = column name containing latitude
  #   collng    = column name containing longitude
  #   proj      = projection data is in. Default in R is EPSG:4326
  #   maxiter   = maximum number of sub-sampling iterations to perform
  #   n.sample  = size of subsample
  # Returns:
  #   vector of area estimates based on bootstrapped data
  
  # Get all unique combinations of locations/rows
  combs <- combn(1:nrow(data), n.sample)
  
  # Limit number of iterations to maxiter
  if (ncol(combs) > maxiter) {
    idx <- sample(1:ncol(combs), size = maxiter, replace = FALSE)
    samples <- combs[ , idx]
  } else {
    samples <- combs
  }
  
  # Calculate the areas on the subsampled occurrences
  areas <- vector(mode = 'numeric', length = ncol(samples)) 
  for (i in 1:ncol(samples)) {
    areas[i] <- draw.convex.hull(data[samples[ , i], ], collat, collng, proj)
  }
 
  return(areas)
}
