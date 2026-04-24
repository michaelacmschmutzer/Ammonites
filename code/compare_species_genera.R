# Compare abundances between species and genus level analysis. Why do they 
# differ so markedly?
# Also compare geographic areas. How good a proxy for species-level trends
# are genera?

setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(scales)
library(stringr)
library(systemfonts)
library(latex2exp)

source('statistical-functions.R')

# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9


###################  Read in data ###################

cepha <- read.csv('../data/cephalopods.csv')

# Abundances (counts) at genus and species level
nauti_genus_abun <- read.csv(
  '../results/genus/analyse_survivalship/nautilids_abundance_raw.csv'
  )
ammon_genus_abun <- read.csv(
  '../results/genus/analyse_survivalship/ammonoids_abundance_raw.csv'
  )
nauti_species_abun <- read.csv(
  '../results/species/analyse_survivalship/nautilids_abundance_raw_species.csv'
  )
ammon_species_abun <- read.csv(
  '../results/species/analyse_survivalship/ammonoids_abundance_raw_species.csv'
  )

# Geographic ranges at genus and species level
nauti_genus_area <- read.csv(
  '../results/genus/geographic_distributions/nautilids_distributions_genus.csv'
  )
ammon_genus_area <- read.csv(
  '../results/genus/geographic_distributions/ammonoids_distributions_genus.csv'
  )
nauti_species_area <- read.csv(
  '../results/species/geographic_distributions/nautilids_distributions_species.csv'
  ) 
ammon_species_area <- read.csv(
  '../results/species/geographic_distributions/ammonoids_distributions_species.csv'
  ) 

# Bootstrapped geographic ranges at genus and species level
nauti_genus_boot <- read.csv(
  '../results/genus/subsampling_distributions/nautilids/bootstrap.csv')
ammon_genus_boot <- read.csv(
  '../results/genus/subsampling_distributions/ammonoids/bootstrap.csv')
nauti_species_boot <- read.csv(
  '../results/species/subsampling_distributions/nautilids/bootstrap.csv')
ammon_species_boot <- read.csv(
  '../results/species/subsampling_distributions/ammonoids/bootstrap.csv')

# Read in survival data for ammonoids at genus level
surv_ammon <- read.csv('../data/ammonoids_extinction_genus.csv')
# Simplify
surv_ammon <- select(surv_ammon, c('genus', 'survival'))


###################  How species-rich is each genus? ###################

# How many species are there per genus in the all-Maastrichtian data?
count_per_genus <- cepha %>%
  filter(!is.na(species)) %>% # Ignoring unassigned at species
  group_by(genus) %>%
  summarise(num.species = n_distinct(species)) %>%
  arrange(desc(num.species))
count_per_genus <- merge(count_per_genus, surv_ammon)

wilcox.test(num.species ~ survival, data = count_per_genus)

genus_nspecies_abun <- merge(count_per_genus, ammon_genus_abun)

plot(log(genus_nspecies_abun$num.species), log(genus_nspecies_abun$n))

p <- ggplot(data = genus_nspecies_abun, 
       aes(x = num.species, y = n)) +
  geom_point(cex = 3, color = darkblue) +
  #  labs(x = 'Ammonoids', y = 'Area (km²)') +
  scale_x_log10(
    'Number of species per genus', 
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x)),
    limits = c(1, 1.2 * max(genus_nspecies_abun$num.species))) +
  scale_y_log10(
    'Number of occurrences in genus',
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x)),
    limits = c(1, 1.2 * max(genus_nspecies_abun$n))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
ggsave(
  '../results/compare_taxon_level/Ammonoids_genus_speciesvsoccurrences.png',
  width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

cor(genus_nspecies_abun$num.species, genus_nspecies_abun$n, method = 'spearman')

###################  Compare abundance ###################

# Plot a histogram of the number of occurrences per species

p <- ggplot(ammon_species_abun, aes(x = n)) +
  geom_histogram(binwidth = 5, fill = darkblue) +
  xlab('Number of occurrences') +
  ylab('Number of species') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
ggsave(
  '../results/compare_taxon_level/Ammonoid_abundance_species.png',
  width = 6, height = 6, units = 'cm', dpi = 600, plot = p)

# Get the genus of each species to act as key in next step
nauti_species_abun <- nauti_species_abun %>%
  mutate(genus = str_split_i(species, ' ', 1))
ammon_species_abun <- ammon_species_abun %>%
  mutate(genus = str_split_i(species, ' ', 1))

# Combine genus and species into same data frame
nauti_abun <- merge(nauti_species_abun, nauti_genus_abun, by = 'genus', 
  suffixes = c('.species', '.genus'))
ammon_abun <- merge(ammon_species_abun, ammon_genus_abun, by = 'genus', 
  suffixes = c('.species', '.genus'))

# Plotting
pamabun <- ggplot(data = ammon_abun, 
    aes(x = n.genus, y = n.species)) +
  geom_point(cex = 3, color = darkblue) +
  xlab('Genus number of occurrences') +
  ylab('Species number of occurrences') +
  stat_summary(
    fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
  stat_summary(
    fun.max = function(z) { quantile(z, 0.75) },
    fun.min = function(z) { quantile(z, 0.25) },
    geom = 'errorbar', color = errbarcol, width = 0.1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
pnaabun <- ggplot(data = nauti_abun, 
    aes(x = n.genus, y = n.species)) +
  geom_point(cex = 3, color = orange) +
  xlab('Genus number of occurrences') +
  ylab('Species number of occurrences') +
  stat_summary(
    fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
  stat_summary(
    fun.max = function(z) { quantile(z, 0.75) },
    fun.min = function(z) { quantile(z, 0.25) },
    geom = 'errorbar', color = errbarcol, width = 0.1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
p <- grid.arrange(pamabun, pnaabun, ncol = 2)
ggsave(
  '../results/compare_taxon_level/Raw_abundance_comparison.png',
  width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

###################  Compare area ###################

# Get the genus of each species to act as key in next step
nauti_species_area <- nauti_species_area %>%
  mutate(genus = str_split_i(species, ' ', 1))
ammon_species_area <- ammon_species_area %>%
  mutate(genus = str_split_i(species, ' ', 1))

# Combine genus and species into same data frame
nauti_area <- merge(nauti_species_area, nauti_genus_area, by = 'genus', 
  suffixes = c('.species', '.genus'))
ammon_area <- merge(ammon_species_area, ammon_genus_area, by = 'genus', 
  suffixes = c('.species', '.genus'))

# Plot genus area vs species area (PALEOMAP only for now)
pamarea <- ggplot(data = ammon_area, 
    aes(x = PALEOMAP.area.km2.genus, y = PALEOMAP.area.km2.species)) +
  #geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
    geom_point(cex = 3, color = darkblue) +
  #  labs(x = 'Ammonoids', y = 'Area (km²)') +
  scale_x_log10(
    TeX('Genus area (km$^{2}$)'), 
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x)),
    limits = c(1, 1.2 * max(ammon_area$PALEOMAP.area.km2.genus))) +
  scale_y_log10(
    TeX('Species area (km$^{2}$)'),
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x)),
    limits = c(1, 1.2 * max(ammon_area$PALEOMAP.area.km2.genus))) +
  stat_summary(
    fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
  stat_summary(
    fun.max = function(z) { quantile(z, 0.75) },
    fun.min = function(z) { quantile(z, 0.25) },
    geom = 'errorbar', color = errbarcol, width = 0.1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
pnaarea <- ggplot(data = nauti_area, 
         aes(x = PALEOMAP.area.km2.genus, y = PALEOMAP.area.km2.species)) +
    #geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
    geom_point(cex = 3, color = orange) +
    #  labs(x = 'Ammonoids', y = 'Area (km²)') +
    scale_x_log10(
      TeX('Genus area (km$^{2}$)'), 
      breaks = scales::trans_breaks('log10', function(x) 10^x),
      labels = scales::trans_format('log10', scales::math_format(10^.x)),
      limits = c(1, 1.2 * max(ammon_area$PALEOMAP.area.km2.genus))) +
    scale_y_log10(
      TeX('Species area (km$^{2}$)'),
      breaks = scales::trans_breaks('log10', function(x) 10^x),
      labels = scales::trans_format('log10', scales::math_format(10^.x)),
      limits = c(1, 1.2 * max(ammon_area$PALEOMAP.area.km2.genus))) +
    stat_summary(
      fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
    stat_summary(
      fun.max = function(z) { quantile(z, 0.75) },
      fun.min = function(z) { quantile(z, 0.25) },
      geom = 'errorbar', color = errbarcol, width = 0.1) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = fsize, colour = 'black', family = font),
      axis.text.y = element_text(size = fsize, colour = 'black', family = font),
      axis.title.y = element_text(size = fsize, family = font),
      axis.title.x = element_text(size = fsize, family = font),
    )
p <- grid.arrange(pamarea, pnaarea, ncol = 2)
ggsave(
  '../results/compare_taxon_level/Raw_area_comparison.png',
  width = 12, height = 6, units = 'cm', dpi = 600, plot = p)
  
###################  Compare bootstrap ###################

# For both genus and species, summarise bootstrap replicates: Use median
nauti_genus_medboot <- nauti_genus_boot %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(PALEOMAP.area.km2 = median(area.km2))
ammon_genus_medboot <- ammon_genus_boot %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(PALEOMAP.area.km2 = median(area.km2))
nauti_species_medboot <- nauti_species_boot %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(PALEOMAP.area.km2 = median(area.km2))
ammon_species_medboot <- ammon_species_boot %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(PALEOMAP.area.km2 = median(area.km2))

# Get the genus of each species to act as key in next step
nauti_species_medboot <- nauti_species_medboot %>%
  mutate(genus = str_split_i(species, ' ', 1))
ammon_species_medboot <- ammon_species_medboot %>%
  mutate(genus = str_split_i(species, ' ', 1))

# Combine genus and species into same data frame
nauti_medboot <- merge(nauti_species_medboot, nauti_genus_medboot, by = 'genus', 
  suffixes = c('.species', '.genus'))
ammon_medboot <- merge(ammon_species_medboot, ammon_genus_medboot, by = 'genus', 
  suffixes = c('.species', '.genus'))

# Plot genus area vs species area (PALEOMAP only for now)
pamboot <- ggplot(data = ammon_medboot, 
                  aes(x = PALEOMAP.area.km2.genus, y = PALEOMAP.area.km2.species)) +
  #geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  geom_point(cex = 3, color = darkblue) +
  #  labs(x = 'Ammonoids', y = 'Area (km²)') +
  scale_x_log10(
    TeX('Genus area (km$^{2}$)'), 
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x)),
    limits = c(1, 1.2 * max(ammon_medboot$PALEOMAP.area.km2.genus))) +
  scale_y_log10(
    TeX('Species area (km$^{2}$)'),
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x)),
    limits = c(1, 1.2 * max(ammon_medboot$PALEOMAP.area.km2.genus))) +
  stat_summary(
    fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
  stat_summary(
    fun.max = function(z) { quantile(z, 0.75) },
    fun.min = function(z) { quantile(z, 0.25) },
    geom = 'errorbar', color = errbarcol, width = 0.1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
pnaboot <- ggplot(data = nauti_medboot, 
                  aes(x = PALEOMAP.area.km2.genus, y = PALEOMAP.area.km2.species)) +
  #geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  geom_point(cex = 3, color = orange) +
  #  labs(x = 'Ammonoids', y = 'Area (km²)') +
  scale_x_log10(
    TeX('Genus area (km$^{2}$)'), 
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x)),
    limits = c(1, 1.2 * max(ammon_medboot$PALEOMAP.area.km2.genus))) +
  scale_y_log10(
    TeX('Species area (km$^{2}$)'),
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x)),
    limits = c(1, 1.2 * max(ammon_medboot$PALEOMAP.area.km2.genus))) +
  stat_summary(
    fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
  stat_summary(
    fun.max = function(z) { quantile(z, 0.75) },
    fun.min = function(z) { quantile(z, 0.25) },
    geom = 'errorbar', color = errbarcol, width = 0.1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
p <- grid.arrange(pamboot, pnaboot, ncol = 2)
ggsave(
  '../results/compare_taxon_level/Bootstrap_comparison.png',
  width = 12, height = 6, units = 'cm', dpi = 600, plot = p)
