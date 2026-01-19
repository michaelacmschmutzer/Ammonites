# Rotate occurrence locations to their geographic position in the Maastrichtian.
# Generate consensus postions using all five rotation models in the palaeoverse
# palaeorotate() function
setwd("~/Documents/Projects/Ammonites/code/")

#library(palaeoverse)
library(stringr)
library(dplyr)
library(rgplates)
library(sf)
library(ggplot2)


cepha <- read.csv('../data/cephalopods.csv')

# Calculate age midpoints
cepha$age <- median((cepha$max_ma + cepha$min_ma) / 2)

# TO DO: Check all points have the same age ("Maastrichtian")

# Rotate using all models. TO DO: Examine uncertainties
rotamodels <- c('MERDITH2021', 'TorsvikCocks2017',
                'GOLONKA', 'MATTHEWS2016_pmag_ref',
                'PALEOMAP')

# Rotation of locations using palaeoverse
#cepha_rota <- palaeorotate(occdf = cepha,
#                          method = 'point',
#                          model = rotamodels,
#                          uncertainty = TRUE)

# Prepare dataframe to insert data
cepha_rota <- cepha


################### PLOTTING ###################

# Define plotting colours
oceanblue <- '#1A6BB0'

for (m in rotamodels) {
  
  # Set anchor to switch off mantle reference frame from 0 to 1 for Torsvik
  if (m == 'TorsvikCocks2017') { anchor <- 1 } else { anchor <- 0 }
  
  # Fossil location rotation with rgplates
  coords <- select(cepha, c(lng, lat)) 
  coords_rota <- reconstruct(coords, age = unique(cepha$age), model = m, 
    anchor = anchor)
  collng <- paste('p_lng', m, sep = '_')
  collat <- paste('p_lat', m, sep = '_')
  cepha_rota[, collng] <- coords_rota[, 1]
  cepha_rota[, collat] <- coords_rota[, 2]
  
  # Plate rotation with rgplates
  coastlines <- reconstruct('coastlines', age = unique(cepha$age), model = m, 
    anchor = anchor)
  
  # Save the rotated coastlines for later plotting
  st_write(coastlines, 
    paste('../results/mapping_files/palaeocoastlines/', m,
    '_coastlines.shp', sep = ''), 
    append = FALSE
  )

  # Get map edges 
  edge <- mapedge()

  # Transform to Robinson
  epsg <- 'ESRI:54030'
  coasts <- sf::st_transform(coastlines, crs = epsg)
  edges <- sf::st_transform(edge, crs = epsg)

  # Also transform locations
  locations <- cepha_rota %>% 
    filter(!is.na(!!sym(collng)) | !is.na(!!sym(collat))) %>%
    sf::st_as_sf(coords = c(collng, collat), crs = 4326) %>%
    sf::st_transform(crs = epsg)

  p <- ggplot() +
    geom_sf() +
    geom_sf(data = edges, colour = 'gray30', fill = oceanblue) +
    geom_sf(data = coasts, colour = 'gray90') + 
    geom_sf(data = locations, colour = 'darkred', size = 1) +
    theme_minimal()
  filename = paste(c('../results/palaeorotations/', m,'_rotation.pdf'), 
                   collapse = '')
  ggsave(filename, plot = p, width = 12, height = 6, units = 'cm')
}

# Get median palaeoposition from p_lat and p_lng predicted by all five models
latcol <- paste('p_lat', rotamodels, sep = '_')
lngcol <- paste('p_lng', rotamodels, sep = '_')

# Check that there is at least one successful rotation for each location
na_num <- rowSums(is.na(cepha_rota[append(latcol, lngcol)]))
if (max(na_num) >= 2 * length(rotamodels) - 2) {
  warning("Some rows contain no successful rotations. The median p_lng and 
          p_lat will be zero for these locations.")
}

# Calculate the medians using all available data. If rotation failed for some 
# models (na values), ignore and calculate median from rotations that worked. 
cepha_rota <- cepha_rota %>% 
  rowwise() %>%
  mutate(
    p_lat = median(c_across(all_of(latcol)), na.rm = TRUE),
    p_lng = median(c_across((all_of(lngcol))), na.rm = TRUE))

# Save rotated positions
write.csv(cepha_rota, '../data/cephalopods_palaeorotated.csv', 
  row.names = FALSE)

# Plot median/consensus locations with the PALEOMAP rotation
# Transform locations
locations <- cepha_rota %>% 
  filter(!is.na('p_lng') | !is.na('p_lat')) %>%
  sf::st_as_sf(coords = c('p_lng', 'p_lat'), crs = 4326) %>%
  sf::st_transform(crs = epsg)

p <- ggplot() +
  geom_sf() +
  geom_sf(data = edges, colour = 'gray30', fill = oceanblue) +
  geom_sf(data = coasts, colour = 'gray90') + 
  geom_sf(data = locations, colour = 'darkred', size = 1) +
  theme_minimal()
filename = '../results/palaeorotations/Consensus_rotation.pdf'
ggsave(filename, plot = p, width = 12, height = 6, units = 'cm')
