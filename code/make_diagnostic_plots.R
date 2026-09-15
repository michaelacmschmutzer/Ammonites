# Make some initial diagnostic plots
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rgplates)
library(DataExplorer)

# TODO: Due to taxonomic mismatches in data from different sources,
# several ammonoid genera are assigned to different orders and appear twice
# TODO: Try Corrplot

###################  Read in data ################### 

# Occurrence data
cepha <- read.csv("../data/cephalopods.csv")
cepha_species <- read.csv("../data/cephalopods_palaeorotated_species.csv")

# Genus-level data
body_sizes <- read.csv(
  '../data/doi_10_5061_dryad_zpc866t5b__v20210211/genus.sizes.ranges.rev.csv')
hatch_nauti_raw <- read.csv('../data/nautilids_embryonic_shell_size.csv')
hatch_ammon_raw <- read.csv('../data/ammonoids_embryonic_shell_size.csv')

# Simplify body and hatching size data, 
# take median of species to get at genus level
hatch_nauti_genus <- hatch_nauti_raw %>%
  group_by(genus) %>%
  summarise(med.hatching.size = median(hatching.size..mm.))
hatch_ammon_genus <- hatch_ammon_raw %>%
  group_by(genus) %>%
  summarise(med.hatching.size = median(hatching.size..mm.))
body_sizes <- dplyr::select(body_sizes, c('genus', 'logvol'))

orders <- unique(cepha$order)
suborders <- unique(cepha[!is.na(cepha$suborder), 'suborder'])

# Define plotting colours
oceanblue <- '#1A6BB0'
colours = c('lightblue', 'aquamarine3', '#12255A', 'red3')
orange <- '#E38C59'
# from https://www.color-hex.com/color-palette/23101
nautilus <- c('#e38c59', '#fcd4ba', '#1789a9', '#71b8c7', '#c2f3fd')

################### Check if taxonomy is unambiguous ################### 

# Get out the unique taxonomy of each genus
genus_taxon <- cepha %>%
  select(c('order', 'suborder', 'superfamily', 'family', 'genus')) %>%
  distinct() 

# Sort by genus
genus_taxon <- genus_taxon[order(genus_taxon$genus), ]

# Tally if each genus is unique. Check if that is the case. Each genus listed
# more than once has ambiguous classification
tally_genus_taxon <- genus_taxon %>%
  count(genus)

if (any(tally_genus_taxon$n > 1)) {
  cont <- readline(
    'Taxonomy is not resolved. Enter "y" if you want to continue anyway')
  if (cont != 'y') {
    stop('Good. Fix the taxonomy.')
  }
}

################### Some preliminary explorations ###################
# DataExplorer

introduce(cepha)
plot_intro(cepha)
missing <- plot_missing(cepha)
ggsave("../results/diagnostic_plots/Missing_rows.pdf", plot = missing,
       width = 30, height = 27, units = 'cm')
plot_bar(cepha)

# Percent of species-level occurrences
nrow(cepha_species) / nrow(cepha)

# Number of ammonoid occurrences genus
sum(cepha$is.nautilid == FALSE)
# Number of nautilid occurrences genus
sum(cepha$is.nautilid == TRUE)
# Number of ammonoid occurrences species
sum(cepha_species$is.nautilid == FALSE)
# Number of nautilid occurrences species
sum(cepha_species$is.nautilid == TRUE)

# Number of ammonoid genera
ammon_genus <- unique(cepha[cepha$is.nautilid == FALSE, 'genus'])
length(ammon_genus)
# Number of nautilid genera
nauti_genus <- unique(cepha[cepha$is.nautilid == TRUE, 'genus'])
length(nauti_genus)
# Number of ammonoid species
ammon_species <- unique(cepha_species[cepha_species$is.nautilid == FALSE, 'species'])
length(ammon_species)
# Number of nautilid species
nauti_species <- unique(cepha_species[cepha_species$is.nautilid == TRUE, 'species'])
length(nauti_species)

# Number of ammonoid genera belonging to each suborder
for (s in suborders) {
  n <- length(unique(cepha[cepha$suborder == s, 'genus']))
  print(paste(n, 'number of genera in', s))
}

# Percentage of species with hatching size data
hatch_ammon_raw <- hatch_ammon_raw %>%
  mutate(binomial.name = paste(genus, species))
sum(ammon_species %in% hatch_ammon_raw$binomial.name) / length(ammon_species) 
hatch_nauti_raw <- hatch_nauti_raw %>%
  mutate(binomial.name = paste(genus, species))
sum(nauti_species %in% hatch_nauti_raw$binomial.name) / length(nauti_species) 
# Percentage of genera with hatching size data
sum(ammon_genus %in% hatch_ammon_genus$genus)/ length(ammon_genus)
sum(nauti_genus %in% hatch_nauti_genus$genus) / length(nauti_genus)
# Percentage of genera with body size data
sum(ammon_genus %in% body_sizes$genus) / length(ammon_genus)
sum(nauti_genus %in% body_sizes$genus) / length(nauti_genus)
# Number of genera with both hatching and body size data
sum(
  (ammon_genus %in% hatch_ammon_genus$genus) & 
  (ammon_genus %in% body_sizes$genus)
)
sum(
  (nauti_genus %in% hatch_nauti_genus$genus) & 
  (nauti_genus %in% body_sizes$genus)
) 


################### Some more explorations ###################

# Plot midpoint time of fossil occurrences
cepha <- cepha %>% mutate(mid_ma = 0.5 * (max_ma + min_ma))
pdf("../results/diagnostic_plots/Time_midpoints.pdf")
par(mfrow=c(2, 2))
for (i in 1:length(orders) ) {
  hist(cepha[cepha$order==orders[i], 'mid_ma'], 
       main=orders[i],
       xlab='Age (midpoint, mya)', 
       ylab='Frequency',
       xlim=rev(range(cepha[cepha$order==orders[i], 'mid_ma'])),
       col=colours[i])
}
dev.off()

# Plot time-ranges of occurrences (histogram)
cepha <- cepha %>% mutate(range_ma = max_ma - min_ma)
pdf("../results/diagnostic_plots/Occurrence_time_ranges.pdf")
par(mfrow=c(2, 2))
for (i in 1:length(orders) ) {
  hist(cepha[cepha$order==orders[i], 'range_ma'], 
       main=orders[i],
       xlab='Occurrence time range (mya)', 
       ylab='Frequency',
       col=colours[i])
}
dev.off()

# Plot present-day distribution of occurrences
colours_named <- colours
names(colours_named) <- orders
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

# Plot map
map <- ggplot() +
  geom_sf() +
  geom_sf(data = edges, colour = 'gray30', fill = oceanblue) +
  geom_sf(data = coasts, colour = 'gray90') + 
  geom_sf(data = locations, colour = 'darkred', size = 1) +
  theme_minimal()
ggsave("../results/diagnostic_plots/World_distribution.pdf", plot = map,
       width = 12, height = 7, units = 'cm')

# Plot number of occurrences per genus
tally_genus <- cepha %>% 
  group_by(order) %>%
  count(genus)
tally_genus <- as.data.frame(tally_genus)
pdf("../results/diagnostic_plots/Num_occurrences_genus.pdf")
par(mfrow=c(2, 2))
for (i in 1:length(orders) ) {
  hist(tally_genus[tally_genus$order==orders[i], 'n'], 
       main=orders[i],
       xlab='Number of occurrences (genus)', 
       ylab='Frequency',
       col=colours[i])
}
dev.off()

# Plot number of occurrences per species
tally_species <- cepha %>% 
  # Remove all entries where species is unknown
  filter(!is.na(species)) %>% 
  group_by(order, genus, species) %>%
  count(genus, species) %>%
  ungroup()
tally_species <- as.data.frame(tally_species)
pdf("../results/diagnostic_plots/Num_occurrences_species.pdf")
par(mfrow=c(2, 2))
for (i in 1:length(orders) ) {
  hist(tally_species[tally_species$order==orders[i], 'n'], 
       main=orders[i],
       xlab='Number of occurrences (species)', 
       ylab='Frequency',
       col=colours[i])
}
dev.off()

# Plot number of unique geographic locations per genus
tally_loc_genus <- cepha %>%
  group_by(order, genus) %>%
  summarise(count = n_distinct(lat, lng)) %>%
  ungroup()
tally_loc_genus <- as.data.frame(tally_loc_genus)
pdf("../results/diagnostic_plots/Num_unique_loc_genus.pdf")
par(mfrow=c(2, 2))
for (i in 1:length(orders) ) {
  hist(tally_loc_genus[tally_loc_genus$order==orders[i], 'count'], 
       main=orders[i],
       xlab='Number of locations (genus)', 
       ylab='Frequency',
       col=colours[i])
}
dev.off()

# Plot number of unique geographic locations per species
tally_loc_species <- cepha %>%
  # Remove all entries where species is unknown
  filter(!is.na(species)) %>% 
  group_by(order, genus, species) %>%
  summarise(count = n_distinct(lat, lng)) %>%
  ungroup()
tally_loc_species <- as.data.frame(tally_loc_species)
pdf("../results/diagnostic_plots/Num_unique_loc_species.pdf")
par(mfrow=c(2, 2))
for (i in 1:length(orders) ) {
  hist(tally_loc_species[tally_loc_species$order==orders[i], 'count'], 
       main=orders[i],
       xlab='Number of locations (species)', 
       ylab='Frequency',
       col=colours[i])
}
dev.off()

# Print some useful summary statistics
suff_loc_genus <- sum(tally_loc_genus$count >= 3)
print(paste('Total number of genera with >= 3 locations:', suff_loc_genus))
nautindex <- tally_loc_genus$order == 'Nautilida'
naut_loc_genus <- sum(tally_loc_genus$count[nautindex] >= 3)
print(paste('... of which are nautiloids:',  naut_loc_genus))
suff_loc_species <- sum(tally_loc_species$count >= 3)
print(paste('Total number of species with >= 3 locations:', suff_loc_species))
nautindex <- tally_loc_species$order == 'Nautilida'
naut_loc_species <- sum(tally_loc_species$count[nautindex] >= 3)
print(paste('... of which are nautiloids:',  naut_loc_species))
