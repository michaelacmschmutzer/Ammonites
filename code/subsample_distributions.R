# Sub-sample geographic locations to check robustness of geographic range 
# reconstructions
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)

source('geospatial-functions.R')

# Command line arguments if run from bash
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0){
  cephas <- 'nautilids'  # 'nautilids' 'ammonoids'
} else {
  cephas <- args[1]
}

cepha_rota <- read.csv('../data/cephalopods_palaeorotated.csv')

# Get either nautilids or ammonoids
if (cephas == 'nautilids') {
  cornu_rota <- cepha_rota[cepha_rota$order=='Nautilida', ]
} else {
  cornu_rota <- cepha_rota[cepha_rota$order!='Nautilida', ]
}

# Rotation models
rotamodels <- c('MERDITH2021', 'TorsvikCocks2017',
                'GOLONKA', 'MATTHEWS2016_pmag_ref',
                'PALEOMAP')

# Use the Equal Earth projection for area estimation
eq_proj <- 'EPSG:8857'
# Use the Robinson projection for plotting
epsg <- 'ESRI:54030'

# set resolution to 10 arc minutes, roughly 25km. Specify in meter
# this sets the length of the side of a grid cell
resolution <- c(25000)

# Create data frames
bootstrap.df <- data.frame(data.frame(matrix(ncol = 3, nrow = 0)))
names(bootstrap.df) <- c('genus', 'area.km2', 'model')
jackknife.df <- data.frame(data.frame(matrix(ncol = 3, nrow = 0)))
names(jackknife.df) <- c('genus', 'area.km2', 'model')

for (m in rotamodels) {
  
  # Columns with rotated geographic position of fossil occurrence for model m
  collng <- paste('p_lng', m, sep = '_')
  collat <- paste('p_lat', m, sep = '_')
  
  # Remove rows with NAs in geographic postion for model m
  cornu_rota_geo <- cornu_rota %>% 
    filter(!is.na(!!sym(collng)) & !is.na(!!sym(collat)))
  
  # Keep only those genera with more than three unique locations for each 
  # reconstruction
  tally_loc_genus <- cornu_rota_geo %>%
    group_by(genus) %>%
    summarise(count = n_distinct(!!sym(collng), !!sym(collat)))
  tally_loc_genus <- as.data.frame(tally_loc_genus)
  genera <- tally_loc_genus %>% 
    filter(count > 3) %>%
    pull(genus)
  
  for (g in genera) {
    
    # Get data for genus
    genus <- cornu_rota_geo %>%
      filter(genus == g)
    
    # Get unique occurrences within the cells of an equal-area grid
    occ_uniq <- unique.within.grid(genus, collat, collng, 'EPSG:4326', eq_proj, 
      resolution)
    
    # if the number of unique locations is more than three, do bootstrapping and
    # jackknifing on the unique occurrences. Perform only maxiter = 1000 
    # iterations on a random set of occurrences (no replacement)
    if (nrow(occ_uniq) > 3) {
      boot.areas <- bootstrap.area(occ_uniq, collat, collng, 'EPSG:4326')
      jack.areas <- jackknife.area(occ_uniq, collat, collng, 'EPSG:4326')
      
      # append to data frames
      bt.df <- data.frame(genus = g, area.km2 = boot.areas, model = m)
      jk.df <- data.frame(genus = g, area.km2 = jack.areas, model = m)
      bootstrap.df <- rbind(bootstrap.df, bt.df)
      jackknife.df <- rbind(jackknife.df, jk.df)
    }
    
  }
  
}

resdir <- paste('../results/subsampling_distributions/', cephas, '/', sep = '')
write.csv(bootstrap.df, paste(resdir, 'bootstrap.csv', sep = ''),
  row.names = FALSE)
write.csv(jackknife.df,paste(resdir, 'jackknife.csv', sep = ''),
  row.names = FALSE)
