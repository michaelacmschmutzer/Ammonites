# Plot and analyse geographic ranges
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(ggplot2)
library(gridExtra)
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

###################  Read in data  ################### 

# Read in geographic distribution data
geo_distr_nauti <- read.csv(
  '../results/geographic_distributions/nautilids_distributions_genus.csv') 
geo_distr_ammon <- read.csv(
  '../results/geographic_distributions/ammonoids_distributions_genus.csv') 

# Read in subsampled geographic distribution data
boot_distr_nauti <- read.csv(
  '../results/subsampling_distributions/nautilids/bootstrap.csv')
boot_distr_ammon <- read.csv(
  '../results/subsampling_distributions/ammonoids/bootstrap.csv')
jack_distr_nauti <- read.csv(
  '../results/subsampling_distributions/nautilids/jackknife.csv')
jack_distr_ammon <- read.csv(
  '../results/subsampling_distributions/ammonoids/jackknife.csv')

# Read in marine areas geographic distribution data
mari_nauti <- read.csv(
  '../results/geographic_distributions/nautilids_marine_areas.csv')
mari_ammon <- read.csv(
  '../results/geographic_distributions/ammonoids_marine_areas.csv')

# Read in survival data 
surv_nauti <- read.csv('../data/nautilids_extinction_genus.csv') 
surv_ammon <- read.csv('../data/ammonoids_extinction_genus.csv')
# Simplify
surv_nauti <- select(surv_nauti, c('genus', 'survival'))
surv_ammon <- select(surv_ammon, c('genus', 'survival'))

###################  Plotting Area vs Survival  ################### 

# Just to show... The area occupied by survivors is not larger than the area
# occupied by genera that went extinct. It is just that survivors are more
# abundant in the areas they occupied

distr_ammon_surv <- merge(geo_distr_ammon, surv_ammon)
distr_nauti_surv <- merge(geo_distr_nauti, surv_nauti)

# Plotting min and max
minarea <- 1 # min(geo_distr_nauti$PALEOMAP.area.km2, na.rm = TRUE) 
maxarea <- max(geo_distr_ammon$PALEOMAP.area.km2, na.rm = TRUE) + 1e4

pdistram <- ggplot(data = distr_ammon_surv, 
                   aes(x = survival, y = PALEOMAP.area.km2)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  #  labs(x = 'Ammonoids', y = 'Area (km²)') +
  scale_x_discrete('Ammonoids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  # scale_y_log10('Area (km²)',
  #   breaks = scales::trans_breaks('log10', function(x) 10^x),
  #   labels = scales::trans_format('log10', scales::math_format(10^.x)),
  #   limits = c(minarea, maxarea)) +
  scale_y_continuous(TeX('Area (km$^{2}$)'),
                     label = scientific_10,
                     limits = c(minarea, maxarea)) +
  #  ylim(0, max(distr_ammon_surv$PALEOMAP.area.km2) + 10) +
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
pdistrna <- ggplot(data = distr_nauti_surv, 
                   aes(x = survival, y = PALEOMAP.area.km2)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = orange) +
  #  labs(x = 'Nautilids', y = 'Area (km²)') +
  scale_x_discrete('Nautilids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  # scale_y_log10('Area (km²)',
  #   breaks = scales::trans_breaks('log10', function(x) 10^x),
  #   labels = scales::trans_format('log10', scales::math_format(10^.x)),
  #   limits = c(minarea, maxarea)) +
  scale_y_continuous(TeX('Area (km$^{2}$)'),
                     label = scientific_10,
                     limits = c(minarea, maxarea)) +
  #  ylim(0, max(distr_ammon_surv$PALEOMAP.area.km2) + 10) +
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
p <- grid.arrange(pdistram, pdistrna, ncol = 2)
ggsave('../results/abundance_survivalship/Genus_area_survival.png',
       width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

# Marine area and survival

mari_ammon_surv <- merge(mari_ammon, surv_ammon)
mari_nauti_surv <- merge(mari_nauti, surv_nauti)

# Plotting min and max
minarea <- 1 # min(mari_nauti$marine.area.km2, na.rm = TRUE) 
maxarea <- max(mari_ammon$marine.area.km2, na.rm = TRUE) + 1e4

pmariam <- ggplot(data = mari_ammon_surv, 
                   aes(x = survival, y = marine.area.km2)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  #  labs(x = 'Ammonoids', y = 'Area (km²)') +
  scale_x_discrete('Ammonoids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  # scale_y_log10('Area (km²)',
  #   breaks = scales::trans_breaks('log10', function(x) 10^x),
  #   labels = scales::trans_format('log10', scales::math_format(10^.x)),
  #   limits = c(minarea, maxarea)) +
  scale_y_continuous(TeX('Marine area (km$^{2}$)'),
                     label = scientific_10,
                     limits = c(minarea, maxarea)) +
  #  ylim(0, max(distr_ammon_surv$PALEOMAP.area.km2) + 10) +
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
pmarina <- ggplot(data = mari_nauti_surv, 
                   aes(x = survival, y = marine.area.km2)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = orange) +
  #  labs(x = 'Nautilids', y = 'Area (km²)') +
  scale_x_discrete('Nautilids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  # scale_y_log10('Area (km²)',
  #   breaks = scales::trans_breaks('log10', function(x) 10^x),
  #   labels = scales::trans_format('log10', scales::math_format(10^.x)),
  #   limits = c(minarea, maxarea)) +
  scale_y_continuous(TeX('Marine area (km$^{2}$)'),
                     label = scientific_10,
                     limits = c(minarea, maxarea)) +
  #  ylim(0, max(distr_ammon_surv$PALEOMAP.area.km2) + 10) +
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
p <- grid.arrange(pmariam, pmarina, ncol = 2)
ggsave('../results/abundance_survivalship/Genus_marine_area_survival.png',
       width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

################### Pairwise comparison Area Estimates ################### 

# For marine data, it only really makes sense to compare PALEOMAP areas,
# because the coastline data from Kocsis & Scotese (2021) are based on the 
# PALEOMAP plate tectonic model

# Merge dataframes for plotting
raw_mari_ammon <- merge(mari_ammon, geo_distr_ammon, by = 'genus') %>%
  select(c(genus, marine.area.km2, PALEOMAP.area.km2))
raw_mari_nauti <- merge(mari_nauti, geo_distr_nauti, by = 'genus') %>%
  select(c(genus, marine.area.km2, PALEOMAP.area.km2)) 

# Plotting min and max
minarea <- 0 # min(mari_nauti$marine.area.km2, na.rm = TRUE) 
maxarea_x <- max(raw_mari_ammon$PALEOMAP.area.km2, na.rm = TRUE) + 1e4
maxarea_y <- max(raw_mari_ammon$marine.area.km2, na.rm = TRUE) + 1e4

pramam <- ggplot(raw_mari_ammon, 
  aes(x = PALEOMAP.area.km2, y = marine.area.km2)) +
  geom_point(cex = 3, color = darkblue) +
  scale_x_continuous(
    TeX('Convex hull area (km$^{2}$)'),
    label = scientific_10,
    limits = c(minarea, maxarea_x)) +
  scale_y_continuous(
    TeX('Marine area (km$^{2}$)'),
    label = scientific_10,
    limits = c(minarea, maxarea_y)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
    legend.position = 'none'
  )
pramna <- ggplot(raw_mari_nauti, 
  aes(x = PALEOMAP.area.km2, y = marine.area.km2)) +
  geom_point(cex = 3, color = orange) +
  scale_x_continuous(
    TeX('Convex hull area (km$^{2}$)'),
    label = scientific_10,
    limits = c(minarea, maxarea_x)) +
  scale_y_continuous(
    TeX('Marine area (km$^{2}$)'),
    label = scientific_10,
    limits = c(minarea, maxarea_y)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
    legend.position = 'none'
  )
p <- grid.arrange(pramam, pramna, ncol = 2)
# ggsave('../results/abundance_survivalship/Genus_marine_area_survival.png',
#        width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

################### Statistical tests Area vs Survival ################### 

aream.exti <- 
  distr_ammon_surv[distr_ammon_surv$survival == FALSE, 'PALEOMAP.area.km2']
aream.surv <- 
  distr_ammon_surv[distr_ammon_surv$survival ==  TRUE, 'PALEOMAP.area.km2']
arena.exti <- 
  distr_nauti_surv[distr_nauti_surv$survival == FALSE, 'PALEOMAP.area.km2']
arena.surv <- 
  distr_nauti_surv[distr_nauti_surv$survival ==  TRUE, 'PALEOMAP.area.km2']

median(aream.exti, na.rm = TRUE)
median(aream.surv, na.rm = TRUE)
median(arena.exti, na.rm = TRUE)
median(arena.surv, na.rm = TRUE)

wilcox.test(PALEOMAP.area.km2 ~ survival, data = distr_ammon_surv,
            na.action = na.pass)
wilcox.test(PALEOMAP.area.km2 ~ survival, data = distr_nauti_surv,
            na.action = na.pass)

# Same for marine-only areas
wilcox.test(marine.area.km2 ~ survival, data = mari_ammon_surv)
wilcox.test(marine.area.km2 ~ survival, data = mari_nauti_surv)

# For bootstrap values. Use medians to avoid pseudo-replication, and only
# consider PALEOMAP for now
boot_ammon_medi <- boot_distr_ammon %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(median.area.km2 = median(area.km2))
boot_nauti_medi <- boot_distr_nauti %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(median.area.km2 = median(area.km2))
# Add survival info
boot_ammon_surv <- merge(boot_ammon_medi, surv_ammon)
boot_nauti_surv <- merge(boot_nauti_medi, surv_nauti)

wilcox.test(median.area.km2 ~ survival, data = boot_ammon_surv)
wilcox.test(median.area.km2 ~ survival, data = boot_nauti_surv)

# For Jackknifing values. Use medians to avoid pseudoreplication, and only
# consider PALEOMAP for now
jack_ammon_medi <- jack_distr_ammon %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(median.area.km2 = median(area.km2))
jack_nauti_medi <- jack_distr_nauti %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(median.area.km2 = median(area.km2))
# Add survival info
jack_ammon_surv <- merge(jack_ammon_medi, surv_ammon)
jack_nauti_surv <- merge(jack_nauti_medi, surv_nauti)

wilcox.test(median.area.km2 ~ survival, data = jack_ammon_surv)
wilcox.test(median.area.km2 ~ survival, data = jack_nauti_surv)

