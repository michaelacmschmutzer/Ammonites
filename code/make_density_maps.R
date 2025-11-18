# Make kernel density maps to show the bias in fossil sampling
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rgplates)
library(tmaptools)

cepha <- read.csv('../data/cephalopods.csv')

orders <- unique(cepha$order)

survloc <- read.csv('../data/ammonoids_survivor_locations.csv')

# Define plotting colours
oceanblue <- '#1A6BB0'
colours = c('lightblue', 'aquamarine3', '#12255A', 'red3')
darkblue <- '#12255A'
orange <- '#E38C59'
darkred <- '#CA282C'

world <- ne_countries(scale = 'medium', returnclass = 'sf')

# Get map edges 
edge <- mapedge()

# Transform to Robinson
epsg <- 'ESRI:54030'
coasts <- sf::st_transform(world, crs = epsg)
edges <- sf::st_transform(edge, crs = epsg)

# Also transform locations
locations <- cepha %>% 
  filter(!is.na('lng') | !is.na('lat')) %>%
  sf::st_as_sf(coords = c('lng', 'lat'), crs = 4326) %>%
  sf::st_transform(crs = epsg)
survpos <- survloc %>%
  sf::st_as_sf(coords = c('lng', 'lat'), crs = 4326) %>%
  sf::st_transform(crs = epsg)

# Plot map
map <- ggplot() +
  geom_sf() +
  geom_sf(data = edges, colour = 'white', fill = oceanblue) +
  geom_sf(data = coasts, colour = 'gray90') + 
  stat_density_2d(
    data = locations,
    mapping = ggplot2::aes(
      x = purrr::map_dbl(geometry, ~.[1]),
      y = purrr::map_dbl(geometry, ~.[2])),
    geom = 'polygon',
    adjust = 0.8,
    contour = TRUE,
    alpha = 0.5) +
  geom_sf(data = locations, colour = 'black', size = 0.5) +
  geom_sf(data = survpos, colour = darkred, shape = 4, size = 1.5) +
  ylab('') +
  xlab('') +
  theme_void()
map_ratio <- get_asp_ratio(edges)
ggsave('../results/diagnostic_plots/Fossil_occurrence_density.pdf', plot = map,
       width = 12, height = 12/map_ratio, units = 'cm')  

