# Estimate geographic range by drawing a convex hull around all fossil
# occurrances of a given taxon. Use consensus palaeorotated locations.
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(terra)
#library(palaeoverse)
library(rgplates)
library(sf)
library(ggplot2)
library(ggtext)
library(reshape2)
library(wesanderson)
library(forcats)

source('geospatial-functions.R')

# Command line arguments if run from bash
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0){
  cephas <- 'ammonoids' # 'nautilids' 'ammonoids'
} else {
  cephas <- args[1]
}

# Import palaeorotated data
cepha_rota <- read.csv('../data/cephalopods_palaeorotated.csv')
# Use the same age regardless of group
age <- unique(cepha_rota$age)

# Get either nautilids or ammonoids
if (cephas == 'nautilids') {
  cornu_rota <- cepha_rota[cepha_rota$order=='Nautilida', ]
} else {
  cornu_rota <- cepha_rota[cepha_rota$order!='Nautilida', ]
}

# Use all rotation models.
rotamodels <- c('MERDITH2021', 'TorsvikCocks2017',
                'GOLONKA', 'MATTHEWS2016_pmag_ref',
                'PALEOMAP')

# Count number of occurrences for each genus
genera_nocc <- cornu_rota %>% count(genus)

# Data frame for results
geo_distr <- data.frame(genus = genera_nocc$genus, n.occ = genera_nocc$n)

# Define plotting colours 
oceanblue <- '#1A6BB0'
darjeeling_palette <- wes_palette('Darjeeling2', n = length(rotamodels))

for (m in rotamodels) {
  
  # Set anchor to switch off mantle reference frame from 0 to 1 for Torsvik
  if (m == 'TorsvikCocks2017') { anchor <- 1 } else { anchor <- 0 }
  
  # Columns with rotated geographic position of fossil occurrence for model m
  collng <- paste('p_lng', m, sep = '_')
  collat <- paste('p_lat', m, sep = '_')

  # Remove rows with NAs in geographic position for model m
  cornu_rota_geo <- cornu_rota %>% 
    filter(!is.na(!!sym(collng)) & !is.na(!!sym(collat)))
  
  # Keep only those genera with three unique locations for each reconstruction
  tally_loc_genus <- cornu_rota_geo %>%
    group_by(genus) %>%
    summarise(count = n_distinct(!!sym(collng), !!sym(collat)))
  tally_loc_genus <- as.data.frame(tally_loc_genus)
  genera <- tally_loc_genus %>% 
    filter(count > 2) %>%
    pull(genus)
  genera_buffer <- tally_loc_genus %>%
    filter(count < 3) %>%
    pull(genus)

  # Make column to store area in
  geo_distr[ , paste(m, 'area.km2', sep = '.') ] <- NA
  
  # Prepare plotting
  
  # Prepare plate rotation for plotting
  coastlines <- st_read(paste('../results/mapping_files/palaeocoastlines/', 
    m, '_coastlines.shp', sep = ''), quiet = TRUE)
  
  # Get map edges 
  edge <- mapedge()
  
  # Transform to Robinson
  epsg <- 'ESRI:54030'
  coasts <- sf::st_transform(coastlines, crs = epsg)
  edges <- sf::st_transform(edge, crs = epsg)
  
  # For each genus and each reconstruction, draw a convex hull
  for (g in tally_loc_genus$genus) {
    
    # Get data for genus
    genus <- cornu_rota_geo %>%
      filter(genus == g)
    
    # Use buffer for genera with fewer than three unique locations
    if (tally_loc_genus[tally_loc_genus$genus == g, 'count'] < 3) {
      use.buffer <- TRUE
    } else {
      use.buffer <- FALSE
    }
    
    # draw hull
    hull_list <- draw.convex.hull(genus, collat, collng, 
      return.hull = TRUE, buffer = use.buffer)
    distr <- hull_list[[1]]
    genus_hull <- hull_list[[2]]
    
    # Convert convex hull to shape file
    genus_area <- sf::st_as_sf(genus_hull)
    
    # Save convex hull for later plotting and processing
    st_write(genus_area, 
      paste('../results/mapping_files/convex_hulls/', m, '_', g,
      '_convex_hull.shp', sep = ''),
      append = FALSE, # To overwrite in case code is run again
      quiet = TRUE
    )
    
    # Store area
    geo_distr[geo_distr$genus == g, paste(m, 'area.km2', sep = '.')] <- distr
    
    # Get percentage of occurrences used for each area estimate
    genus_perc <- nrow(genus) / geo_distr[geo_distr$genus == g, 'n.occ']
    percol <- paste(m, 'percent.occ', sep = '.')
    geo_distr[geo_distr$genus == g, percol] <- genus_perc * 100
    
    
    ################### PLOTTING ###################
    
    # Set projection on covex hull
    genus_area <- genus_area %>% 
#      st_as_sfc() %>%
#      st_segmentize(1e7) %>%
      st_set_crs(4326)
    
    # Transform area to Robinson
    area <- sf::st_transform(genus_area, crs = epsg)
    
    # Also transform locations
    locations <- genus %>% 
      sf::st_as_sf(coords = c(collng, collat), crs = 4326) %>%
      sf::st_transform(crs = epsg)
    
    # Plot
    p <- ggplot() +
      geom_sf() +
      geom_sf(data = edges, colour = 'gray30', fill = oceanblue) +
      geom_sf(data = coasts, colour = 'gray90') +
      geom_sf(data = area, colour = 'darkred', fill = 'darkred', alpha = 0.3) +
      geom_sf(data = locations, colour = 'gray20', size = 1) +
      labs(x = bquote(italic(.(g))), y = '') +
      theme_minimal() +
      theme(axis.text.x = element_text(size = 7))
    filename = paste('../results/geographic_distributions/', cephas, '/', m, 
                     '_', g, '_distribution.pdf', sep = '')
    ggsave(filename, plot = p, width = 12, height = 7, units = 'cm')
    
  }
}

################### PLOTTING ###################

# Prepare for ggplot

# Filter out unneeded columns and reshape
keepcols <- append('genus', paste(rotamodels, 'area.km2', sep = '.'))
geo_distr_plot <- select(geo_distr, all_of(keepcols))
geo_distr_melt <- melt(geo_distr_plot, id.vars = 'genus', 
                       variable.name = 'model', value.name = 'area')
geo_distr_melt <- filter(geo_distr_melt, !is.na(area))

# Add information about survival
survfile <- paste('../data/', cephas, '_extinction_genus.csv', sep = '')
survival <- read.csv(survfile)
survival <- select(survival, c('genus', 'survival'))
geo_distr_merge <- merge(geo_distr_melt, survival, by = 'genus')

extinct <- survival[survival$survival==FALSE, 'genus']

# Need different genus name rotation for ammonoids/nautiloids
if (cephas == 'ammonoids') {
  rot <- 90
  hadj <- 1
  vadj <- 0.5
} else {
  rot <- 15
  vadj <- 0.9
  hadj <- 0.9
}

# Plot areas on log10 axes
p <- geo_distr_merge %>%
  mutate(genus = fct_reorder(genus, survival)) %>%
  ggplot(na.rem = TRUE, 
  aes(x = genus, y = area, colour = model)) +
  geom_jitter(position = position_jitter(0.4), cex = 3) +
  stat_summary(fun = median,
#               fun.max = function(z) { quantile(z, 0.75) },
#               fun.min = function(z) { quantile(z, 0.25) },
    geom = 'point', shape = 18, size = 3.5, # For interquantile range use 'pointrange'
    color = 'red3') +
  scale_y_log10('Area (km²)',
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x))) +
  labs(x = 'Genus') +
  scale_x_discrete(labels = ~ if_else(
    .x %in% extinct, paste0("<span style='color: red3'>", .x, "</span>"), .x
  )) +
  scale_colour_manual(labels = c('Merdith 2021', 
    'Torsvik & Cocks 2017', 'Golonka', 'Matthews 2016', 'Paleomap'), 
    values = darjeeling_palette) + 
  theme_classic() +
  theme(
    axis.text.x = element_markdown(size = 8, face = 'italic',  colour = 'black',
        angle = rot, hjust = hadj, vjust = vadj),
    axis.text.y = element_text(size = 8, colour = 'black'),
    axis.title = element_text(size = 8),
    legend.title = element_blank(), 
    legend.text = element_text(size= 8),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0., 'mm'),
    legend.position = 'bottom',
    legend.justification = 'left',
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-10, 1, 1, -30))
plotfile <- paste('../results/geographic_distributions/', cephas,
                  '/', cephas, '_geographic_ranges.pdf', sep = '')
ggsave(plotfile, width = 12, height = 8, units = 'cm', plot = p)



################### SUMMARY STATS ###################

# Estimate median geographic range and interquartile range
area_col <- paste(rotamodels, 'area.km2', sep = '.')
geo_distr <- geo_distr %>%
  rowwise() %>%
  mutate('median area (km2)' = median(c_across(all_of(area_col)), na.rm = TRUE),
         'interquartile range' = IQR(c_across(all_of(area_col)), na.rm =TRUE),
         'n.models' = sum(!is.na(c_across(all_of(area_col)))))

# Save geographic areas
resfile <- paste('../results/geographic_distributions/', cephas,
  '_distributions_genus.csv', sep = '')
write.csv(geo_distr, resfile, row.names = FALSE)
