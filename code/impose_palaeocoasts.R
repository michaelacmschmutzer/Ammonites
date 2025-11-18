# Impose reconstructed palaeocoastlines onto convex hull geographic range 
# approximations
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(chronosphere)
library(via)
library(sf)
library(terra)
library(ggplot2)

# TO DO: The palaeo-coastlines from Kocsis & Scotese (2021) are reconstructed for 
# the end of the Maastrichtian, but the palaeo-coordinates I use to draw convex
# hulls and plot palaeo-positions are located more in the middle of the 
# Maastrichtian. The difference is small but visible in e.g. the position of 
# India. What to do?

# Switch off sf working with spherical geometry. It causes strange distortions
sf_use_s2(FALSE)

# Command line arguments if run from bash
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0){
  cephas <- 'ammonoids' # 'nautilids' 'ammonoids'
} else {
  cephas <- args[1]
}

# Get the latest version of the Kocsis & Scotese (2021) PaleoMAP coastline data
# Check if already downloaded...
coasts_downloaded <- file.exists(
  '../results/mapping_files/palaeocoastlines/Maastrichtian_coastlines.shp')

if (coasts_downloaded == FALSE) {
  # ... if not, download the coastline files and save...
  coastlines <-
    chronosphere::fetch(src = "paleomap", ser = "paleocoastlines", ver = "7")
  
  maast_coast_shx <- coastlines@stack$'65Ma_CS_v7.shx'
  maast_coast <- sf::st_as_sf(maast_coast_shx)
  
  # Save Maastrichtian coastlines for later use
  st_write(
    maast_coast,
    '../results/mapping_files/palaeocoastlines/Maastrichtian_coastlines.shp',
    quiet = TRUE
  ) 
} else {
  # ... if yes, just read them in. 
  maast_coast <- st_read(
    '../results/mapping_files/palaeocoastlines/Maastrichtian_coastlines.shp',
    quiet = TRUE)
}

# Get all unique genera
cepha_rota <- read.csv('../data/cephalopods_palaeorotated.csv')
# Use the same age regardless of group
age <- unique(cepha_rota$age)

# Get either nautilids or ammonoids
if (cephas == 'nautilids') {
  genera <- unique(cepha_rota[cepha_rota$order=='Nautilida', 'genus'])
} else {
  genera <- unique(cepha_rota[cepha_rota$order!='Nautilida', 'genus'])
}

# Prepare plates for plotting
# plates <- st_read(
#   '../results/mapping_files/palaeocoastlines/PALEOMAP_coastlines.shp',
#   quiet = TRUE)

marine_area <- data.frame(genus = genera, marine.area.km2 = 0)

for (g in genera) {
  
  # Load convex hull for the genus
  genus_area <- st_read(paste('../results/mapping_files/convex_hulls/', 
    'PALEOMAP_', g,'_convex_hull.shp', sep = ''), quiet = TRUE)
  
  # Cut the convex hull so that it only encompasses regions that were marine
  # at the end of the Maastrichtian
  trim_area <- st_difference(genus_area, st_union(st_combine(maast_coast)))
  
  # What is the area of the trimmed-down convex hull? 
  trim_vect <- terra::vect(trim_area)
  distr <- terra::expanse(trim_vect, unit = 'km', transform = TRUE)
  
  # Some low occurrence genera have fossils that fall entirely outside the area
  # that was marine at the the end of the Maastrichtian. The result is that 
  # distr is an empty vector. Keep 'marine.area.km2' zero in that case
  if (length(distr) > 0) {
    marine_area[marine_area$genus == g, 'marine.area.km2'] <- distr  
  }
}

# Save
resfile <- paste('../results/geographic_distributions/', cephas, 
  '_marine_areas.csv', sep = '')
write.csv(marine_area, resfile, row.names = FALSE)

# For plotting palaeo-coastlines and underlying continents
# 
# ggplot() +
#   geom_sf() +
#   geom_sf(data = plates, colour = 'black', fill = 'black') +
#   geom_sf(data = maast_coast, colour = 'gray30', fill = 'gray30') +
#   geom_sf(data = genus_area, colour = 'lightgray', fill = 'lightgray',
#     alpha = 0.8) +
#   geom_sf(data = trim_area, colour = 'red4', fill = 'red4', alpha = 0.6)
