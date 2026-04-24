# Plot and analyse geographic ranges at the species level
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
  '../results/species/geographic_distributions/nautilids_distributions_species.csv') 
geo_distr_ammon <- read.csv(
  '../results/species/geographic_distributions/ammonoids_distributions_species.csv') 

# Read in subsampled geographic distribution data
boot_distr_nauti <- read.csv(
  '../results/species/subsampling_distributions/nautilids/bootstrap.csv')
boot_distr_ammon <- read.csv(
  '../results/species/subsampling_distributions/ammonoids/bootstrap.csv')
jack_distr_nauti <- read.csv(
  '../results/species/subsampling_distributions/nautilids/jackknife.csv')
jack_distr_ammon <- read.csv(
  '../results/species/subsampling_distributions/ammonoids/jackknife.csv')

# Read in survival data 
surv_nauti <- read.csv('../data/nautilids_extinction_species.csv') 
surv_ammon <- read.csv('../data/ammonoids_extinction_species.csv')
# Use binomial species name as identifier and simplify
surv_nauti <- surv_nauti %>%
  mutate(species = paste(genus, species)) %>%
  select(c(species, survival))
surv_ammon <- surv_ammon %>%
  mutate(species = paste(genus, species)) %>%
  select(c(species, survival))

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
  scale_y_log10(TeX('Area (km$^{2}$)'),
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x)),
    limits = c(minarea, maxarea)) +
  # scale_y_continuous(TeX('Area (km$^{2}$)'),
  #                    label = scientific_10,
  #                    #limits = c(minarea, maxarea)
  #                    ) +
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
  scale_y_log10(TeX('Area (km$^{2}$)'),
    breaks = scales::trans_breaks('log10', function(x) 10^x),
    labels = scales::trans_format('log10', scales::math_format(10^.x)),
    limits = c(minarea, maxarea)) +
  # scale_y_continuous(TeX('Area (km$^{2}$)'),
  #                    label = scientific_10,
  #                    #limits = c(minarea, maxarea)
  #                    ) +
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
ggsave(
  '../results/species/analyse_survivalship/Species_area_survival.png',
  width = 12, height = 6, units = 'cm', dpi = 600, plot = p)


wilcox.test(PALEOMAP.area.km2 ~ survival, data = distr_ammon_surv)
wilcox.test(PALEOMAP.area.km2 ~ survival, data = distr_nauti_surv)

###################  Plotting Bootstrap vs Survival  ################### 

# Get medians for each species. Use only PALEOMAP for now
boot_ammon_median <- boot_distr_ammon %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(median.area.km2 = median(area.km2))
boot_nauti_median <- boot_distr_nauti %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(median.area.km2 = median(area.km2))

boot_ammon_surv <- merge(boot_ammon_median, surv_ammon)
boot_nauti_surv <- merge(boot_nauti_median, surv_nauti)


# Plotting min and max
minarea <- 1 #
maxarea <- max(boot_ammon_surv$median.area.km2, na.rm = TRUE) + 1e4

pdistram <- ggplot(data = boot_ammon_surv, 
                   aes(x = survival, y = median.area.km2)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  #  labs(x = 'Ammonoids', y = 'Area (km²)') +
  scale_x_discrete('Ammonoids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  scale_y_log10(TeX('Area (km$^{2}$)'),
                breaks = scales::trans_breaks('log10', function(x) 10^x),
                labels = scales::trans_format('log10', scales::math_format(10^.x)),
                limits = c(minarea, maxarea)) +
  # scale_y_continuous(TeX('Area (km$^{2}$)'),
  #                    label = scientific_10,
  #                    #limits = c(minarea, maxarea)
  #                    ) +
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
pdistrna <- ggplot(data = boot_nauti_surv, 
                   aes(x = survival, y = median.area.km2)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = orange) +
  #  labs(x = 'Nautilids', y = 'Area (km²)') +
  scale_x_discrete('Nautilids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  scale_y_log10(TeX('Area (km$^{2}$)'),
                breaks = scales::trans_breaks('log10', function(x) 10^x),
                labels = scales::trans_format('log10', scales::math_format(10^.x)),
                limits = c(minarea, maxarea)) +
  # scale_y_continuous(TeX('Area (km$^{2}$)'),
  #                    label = scientific_10,
  #                    #limits = c(minarea, maxarea)
  #                    ) +
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
ggsave(
  '../results/species/analyse_survivalship/Species_area_survival_boot.png',
  width = 12, height = 6, units = 'cm', dpi = 600, plot = p)


###################  Plotting Jackknife vs Survival  ################### 

# Get medians for each species. Use only PALEOMAP for now
jack_ammon_median <- jack_distr_ammon %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(median.area.km2 = median(area.km2))
jack_nauti_median <- jack_distr_nauti %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(median.area.km2 = median(area.km2))

jack_ammon_surv <- merge(jack_ammon_median, surv_ammon)
jack_nauti_surv <- merge(jack_nauti_median, surv_nauti)

# Plotting min and max
minarea <- 1 # 
maxarea <- max(jack_ammon_surv$median.area.km2, na.rm = TRUE) + 1e4

pdistram <- ggplot(data = jack_ammon_surv, 
                   aes(x = survival, y = median.area.km2)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  #  labs(x = 'Ammonoids', y = 'Area (km²)') +
  scale_x_discrete('Ammonoids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  scale_y_log10(TeX('Area (km$^{2}$)'),
                breaks = scales::trans_breaks('log10', function(x) 10^x),
                labels = scales::trans_format('log10', scales::math_format(10^.x)),
                limits = c(minarea, maxarea)) +
  # scale_y_continuous(TeX('Area (km$^{2}$)'),
  #                    label = scientific_10,
  #                    #limits = c(minarea, maxarea)
  #                    ) +
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
pdistrna <- ggplot(data = jack_nauti_surv, 
                   aes(x = survival, y = median.area.km2)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = orange) +
  #  labs(x = 'Nautilids', y = 'Area (km²)') +
  scale_x_discrete('Nautilids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  scale_y_log10(TeX('Area (km$^{2}$)'),
                breaks = scales::trans_breaks('log10', function(x) 10^x),
                labels = scales::trans_format('log10', scales::math_format(10^.x)),
                limits = c(minarea, maxarea)) +
  # scale_y_continuous(TeX('Area (km$^{2}$)'),
  #                    label = scientific_10,
  #                    #limits = c(minarea, maxarea)
  #                    ) +
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
ggsave(
  '../results/species/analyse_survivalship/Species_area_survival_jack.png',
  width = 12, height = 6, units = 'cm', dpi = 600, plot = p)
