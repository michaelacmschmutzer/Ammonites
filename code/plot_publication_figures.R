# Plot publication figures
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(ggplot2)
library(forcats)
library(ggpubr)
library(ggguides)
library(scales)
library(stringr)
library(systemfonts)
library(latex2exp)
library(ggtext)
library(reshape2)

source('statistical-functions.R')

# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Times'
fsize <- 8
cex_size <- 2
jitter <- 0.2

################### Compile plotting data sets ###################
cepha <- read.csv('../data/cephalopods.csv')

# Genus-level data
body_sizes <- read.csv(
  '../data/doi_10_5061_dryad_zpc866t5b__v20210211/genus.sizes.ranges.rev.csv')
geo_nauti_genus <- read.csv(
  '../results/genus/geographic_distributions/nautilids_distributions_genus.csv')
geo_ammon_genus <- read.csv(
  '../results/genus/geographic_distributions/ammonoids_distributions_genus.csv')
geo_nauti_boot_genus <- read.csv(
  '../results/genus/subsampling_distributions/nautilids/bootstrap.csv')
geo_ammon_boot_genus <- read.csv(
  '../results/genus/subsampling_distributions/ammonoids/bootstrap.csv')
geo_nauti_jack_genus <- read.csv(
  '../results/genus/subsampling_distributions/nautilids/jackknife.csv')
geo_ammon_jack_genus <- read.csv(
  '../results/genus/subsampling_distributions/ammonoids/jackknife.csv')

hatch_nauti <- read.csv('../data/nautilids_embryonic_shell_size.csv')
hatch_ammon <- read.csv('../data/ammonoids_embryonic_shell_size.csv')
surv_nauti_genus <- read.csv('../data/nautilids_extinction_genus.csv')
surv_ammon_genus <- read.csv('../data/ammonoids_extinction_genus.csv')

# Species-level data
geo_nauti_species <- read.csv(
  '../results/species/geographic_distributions/nautilids_distributions_species.csv')
geo_ammon_species <- read.csv(
  '../results/species/geographic_distributions/ammonoids_distributions_species.csv')
geo_nauti_boot_species <- read.csv(
  '../results/species/subsampling_distributions/nautilids/bootstrap.csv')
geo_ammon_boot_species <- read.csv(
  '../results/species/subsampling_distributions/ammonoids/bootstrap.csv')
geo_nauti_jack_species <- read.csv(
  '../results/species/subsampling_distributions/nautilids/jackknife.csv')
geo_ammon_jack_species <- read.csv(
  '../results/species/subsampling_distributions/ammonoids/jackknife.csv')

surv_nauti_species <- read.csv('../data/nautilids_extinction_species.csv')
surv_ammon_species <- read.csv('../data/ammonoids_extinction_species.csv')

# Read in effect size and power analysis results
effect.sizes.ammon <- read.csv(
  '../results/comparing_hypotheses/ammonoids_effect_sizes.csv')
effect.sizes.nauti <- read.csv(
  '../results/comparing_hypotheses/nautilids_effect_sizes.csv')
am.power.df <- read.csv('../results/comparing_hypotheses/ammonoids_power.csv',
  check.names = FALSE)
na.power.df <- read.csv('../results/comparing_hypotheses/nautilids_power.csv',
  check.names = FALSE)

# Ammonoid genera
ammon_genus <- unique(cepha[cepha$is.nautilid == FALSE, 'genus'])
# Nautilid genera
nauti_genus <- unique(cepha[cepha$is.nautilid == TRUE, 'genus'])

# Place genus and species name in single column
surv_nauti_species <- surv_nauti_species %>%
  mutate(species = paste(genus, species))
surv_ammon_species <- surv_ammon_species %>%
  mutate(species = paste(genus, species))

# Simplify data, take median of species to get at genus level
hatch_nauti_genus <- hatch_nauti %>%
  group_by(genus) %>%
  summarise(med.hatching.size = median(hatching.size..mm.))
hatch_ammon_genus <- hatch_ammon %>%
  group_by(genus) %>%
  summarise(med.hatching.size = median(hatching.size..mm.))
body_sizes <- dplyr::select(body_sizes, c('genus', 'logvol'))

surv_nauti_genus <- dplyr::select(surv_nauti_genus, c('genus', 'survival'))
surv_ammon_genus <- dplyr::select(surv_ammon_genus, c('genus', 'survival'))
surv_nauti_species <- dplyr::select(surv_nauti_species, c('species', 'survival'))
surv_ammon_species <- dplyr::select(surv_ammon_species, c('species', 'survival'))

# Simplify sub-sampling results
# Use PALEOMAP for now
boot_nauti_genus <- geo_nauti_boot_genus %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(boot.median.area = median(area.km2))
boot_ammon_genus <- geo_ammon_boot_genus %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(boot.median.area = median(area.km2))
jack_nauti_genus <- geo_nauti_jack_genus %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(jack.median.area = median(area.km2))
jack_ammon_genus <- geo_ammon_jack_genus %>%
  filter(model == 'PALEOMAP') %>%
  group_by(genus) %>%
  summarise(jack.median.area = median(area.km2))

boot_nauti_species <- geo_nauti_boot_species %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(boot.median.area = median(area.km2))
boot_ammon_species <- geo_ammon_boot_species %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(boot.median.area = median(area.km2))
jack_nauti_species <- geo_nauti_jack_species %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(jack.median.area = median(area.km2))
jack_ammon_species <- geo_ammon_jack_species %>%
  filter(model == 'PALEOMAP') %>%
  group_by(species) %>%
  summarise(jack.median.area = median(area.km2))

# Handle body sizes
bodyvol_nauti <- body_sizes %>% filter(genus %in% nauti_genus)
bodyvol_ammon <- body_sizes %>% filter(genus %in% ammon_genus)

# Merge with survival data
hatch_nauti <- na.omit(merge(hatch_nauti_genus, surv_nauti_genus))
hatch_ammon <- na.omit(merge(hatch_ammon_genus, surv_ammon_genus))

bodyvol_nauti <- merge(bodyvol_nauti, surv_nauti_genus)
bodyvol_ammon <- merge(bodyvol_ammon, surv_ammon_genus)

geo_nauti_genus <- merge(geo_nauti_genus, surv_nauti_genus)
geo_ammon_genus <- merge(geo_ammon_genus, surv_ammon_genus)

boot_nauti_genus <- merge(boot_nauti_genus, surv_nauti_genus)
boot_ammon_genus <- merge(boot_ammon_genus, surv_ammon_genus)

jack_nauti_genus <- merge(jack_nauti_genus, surv_nauti_genus)
jack_ammon_genus <- merge(jack_ammon_genus, surv_ammon_genus)

geo_nauti_species <- merge(geo_nauti_species, surv_nauti_species)
geo_ammon_species <- merge(geo_ammon_species, surv_ammon_species)

boot_nauti_species <- merge(boot_nauti_species, surv_nauti_species)
boot_ammon_species <- merge(boot_ammon_species, surv_ammon_species)

jack_nauti_species <- merge(jack_nauti_species, surv_nauti_species)
jack_ammon_species <- merge(jack_ammon_species, surv_ammon_species)

# Reformat power analysis results
# Ammonoids
am.power.df.plt <- melt(am.power.df)
colnames(am.power.df.plt) <- c('variable', 'eff.size', 'power')
am.power.df.plt$variable <- fct_rev(factor(am.power.df.plt$variable, 
  levels = effect.sizes.ammon$variable))
# Reframe effect sizes so that they are all > 0.5 
plot.ef.ammon <- effect.sizes.ammon %>%
  mutate(rel.eff.size = if_else(eff.size < 0.5, 1 - eff.size, eff.size))
# rounded for plotting
plot.ef.ammon$eff.size <- 
  as.factor(round(plot.ef.ammon$rel.eff.size, digits = 1))
# Nautilids
na.power.df.plt <- melt(na.power.df)
colnames(na.power.df.plt) <- c('variable', 'eff.size', 'power')
na.power.df.plt$variable <- fct_rev(factor(na.power.df.plt$variable, 
  levels = effect.sizes.nauti$variable))
# Reframe effect sizes so that they are all > 0.5 
plot.ef.nauti <- effect.sizes.nauti %>%
  mutate(rel.eff.size = if_else(eff.size < 0.5, 1 - eff.size, eff.size))
# rounded for plotting
plot.ef.nauti$eff.size <- 
  as.factor(round(plot.ef.nauti$rel.eff.size, digits = 1))

################### Plot main results figure: DATA ###################

# Plotting min and max
minarea <- 1
maxarea <- max(geo_ammon_genus$PALEOMAP.area.km2, na.rm = TRUE) + 1e4

# HATCHING SIZES

plt_hat_am_gen <- ggplot(data = hatch_ammon, 
         aes(x = survival, y = med.hatching.size)) +
  geom_jitter(position = position_jitter(jitter), cex = cex_size, color = darkblue) +
  labs(x = 'Ammonoids', y = 'Median hatching\nsize (mm)') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
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
plt_hat_na_gen <- ggplot(data = hatch_nauti,
                   aes(x = survival, y = med.hatching.size)) +
  geom_jitter(position = position_jitter(jitter), cex = cex_size, color = orange) +
  labs(x = 'Nautilids', y = 'Median hatching\nsize (mm)') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
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

# BODY VOLUMES

plt_bod_am_gen <- ggplot(data = bodyvol_ammon, aes(x = survival, y = logvol)) +
  geom_jitter(position = position_jitter(jitter), cex = cex_size, color = darkblue) +
  xlab('Ammonoids') +
  ylab('Body volume<br>(log<sub>10</sub> mm<sup>3</sup>)') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  ylim(0, 7) +
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
    axis.title.y = element_markdown(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
plt_bod_na_gen <- ggplot(data = bodyvol_nauti, aes(x = survival, y = logvol)) +
  geom_jitter(position = position_jitter(jitter), cex = cex_size, color = orange) +
  xlab('Nautilids') +
  ylab('Body volume<br>(log<sub>10</sub> mm<sup>3</sup>)') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  ylim(0, 7) +
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
    axis.title.y = element_markdown(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )

# GEOGRAPHIC RANGE SIZES 

plt_geo_am_gen <- ggplot(data = geo_ammon_genus, 
    aes(x = survival, y = PALEOMAP.area.km2)) +
  geom_jitter(position = position_jitter(jitter), cex = cex_size, color = darkblue) +
  scale_x_discrete('Ammonoids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  scale_y_log10(TeX('Genus area (km$^{2}$)'),
                breaks = scales::trans_breaks('log10', function(x) 10^x),
                labels = scales::trans_format('log10', scales::math_format(10^.x)),
                limits = c(minarea, maxarea)) +
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
plt_geo_na_gen <- ggplot(data = geo_nauti_genus, 
    aes(x = survival, y = PALEOMAP.area.km2)) +
  geom_jitter(position = position_jitter(jitter), cex = cex_size, color = orange) +
  scale_x_discrete('Nautilids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  scale_y_log10(TeX('Genus area (km$^{2}$)'),
                breaks = scales::trans_breaks('log10', function(x) 10^x),
                labels = scales::trans_format('log10', scales::math_format(10^.x)),
                limits = c(minarea, maxarea)) +
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


plt_geo_am_spp <- ggplot(data = geo_ammon_species, 
                         aes(x = survival, y = PALEOMAP.area.km2)) +
  geom_jitter(position = position_jitter(jitter), cex = cex_size, color = darkblue) +
  scale_x_discrete('Ammonoids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  scale_y_log10(TeX('Species area (km$^{2}$)'),
                breaks = scales::trans_breaks('log10', function(x) 10^x),
                labels = scales::trans_format('log10', scales::math_format(10^.x)),
                limits = c(minarea, maxarea)) +
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
plt_geo_na_spp <- ggplot(data = geo_nauti_species, 
                         aes(x = survival, y = PALEOMAP.area.km2)) +
  geom_jitter(position = position_jitter(jitter), cex = cex_size, color = orange) +
  scale_x_discrete('Nautilids', 
                   labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  scale_y_log10(TeX('Species area (km$^{2}$)'),
                breaks = scales::trans_breaks('log10', function(x) 10^x),
                labels = scales::trans_format('log10', scales::math_format(10^.x)),
                limits = c(minarea, maxarea)) +
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

# COMBINE ALL PLOTS

p <- ggarrange(
  plt_hat_am_gen,
  plt_hat_na_gen,
  plt_geo_am_gen, 
  plt_geo_na_gen,
  plt_bod_am_gen,
  plt_bod_na_gen,
  plt_geo_am_spp, 
  plt_geo_na_spp, 
  nrow = 2, ncol = 4, 
  labels = c('A)', 'B)', 'E)', 'F)', 'C)', 'D)', 'G)', 'H)'),
  font.label = list(size = fsize + 3, family = font))
ggsave('../writing/figures/Data_figure.png',
  width = 16, height = 8, units = 'cm', dpi = 600, plot = p)


################# Plot main results figure: EFFECT AND POWER #################

# Make sure variables are ordered as specified in data frame row order
effect.sizes.ammon$variable <- fct_rev(factor(effect.sizes.ammon$variable, 
  levels = effect.sizes.ammon$variable))
effect.sizes.nauti$variable <- fct_rev(factor(effect.sizes.nauti$variable,
  levels = effect.sizes.nauti$variable))


# Plotting ... ammonoids
plt_eff_am <- ggplot(effect.sizes.ammon, aes(x = eff.size, y = variable)) +
  geom_hline(yintercept = nrow(effect.sizes.ammon) + 0.6, linewidth = 1) +
  geom_vline(xintercept = 0.5, linetype = 'longdash', linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci.lower, xmax = ci.upper, width = 0.3),
                color = 'black') +
  geom_point(size = cex_size, color = darkblue) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1), 
                     expand = c(0,0)) +
  labs(x = 'Effect size', y = '') +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.x = element_text(size = fsize, colour = 'black', family = font),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 2, r = 10, b = 0, l = 0, unit = 'pt')
  )
# And nautilids
plt_eff_na <- ggplot(effect.sizes.nauti, aes(x = eff.size, y = variable)) +
  geom_hline(yintercept = nrow(effect.sizes.nauti) + 0.6, linewidth = 1) +
  geom_vline(xintercept = 0.5, linetype = 'longdash', linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci.lower, xmax = ci.upper, width = 0.3),
                color = 'black') +
  geom_point(size = cex_size, color = orange) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1), 
                     expand = c(0,0)) +
  labs(x = 'Effect size', y = '') +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.x = element_text(size = fsize, colour = 'black', family = font),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 2, r = 10, b = 0, l = 0, unit = 'pt')
  )

# and power analysis
plt_pwr_am <- ggplot(am.power.df.plt, aes(x = eff.size, y = variable)) +
  geom_tile(aes(fill = power)) +
  scale_fill_distiller(name = 'Power', palette = 'Reds', direction = 1,
    limits = c(0,1)) +
  geom_point(data = plot.ef.ammon, aes(x = eff.size, y = variable)) +
  labs(x = 'Effect size', y = '') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.x = element_text(size = fsize, colour = 'black', family = font),
    legend.title = element_text(size = fsize, colour = 'black', family = font),
    legend.text = element_text(size = fsize, colour = 'black', family = font),
    legend.position = 'none',
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
  )
plt_pwr_na <- ggplot(na.power.df.plt, aes(x = eff.size, y = variable)) +
  geom_tile(aes(fill = power)) +
  scale_fill_distiller(name = 'Power', palette = 'Reds', direction = 1, 
#    breaks = c(0.1, 0.5, 0.9), 
    expand = c(0,0),
    limits = c(0, 1), guide = guide_colourbar(barwidth = 12, barheight = 0.5,
    direction = 'horizontal')) +
#  colorbar_style(width = 3, height = 0.5, frame = TRUE, aesthetic = "fill") +
  geom_point(data = plot.ef.nauti, aes(x = eff.size, y = variable)) +
  labs(x = 'Effect size', y = '') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.x = element_text(size = fsize, colour = 'black', family = font),
    legend.title = element_text(size = fsize, colour = 'black', family = font),
    legend.text = element_text(size = fsize, colour = 'black', family = font),
    legend.position = 'bottom',
    legend.margin = margin(-5.5, 5.5, 5.5, -115),
    plot.margin = margin(5.5, 5.5, -5.5, 5.5),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
  ) 

p <- ggarrange(
  plt_eff_am,
  plt_eff_na,
  plt_pwr_am,
  plt_pwr_na,
  nrow = 2, ncol = 2, 
  labels = c('A)', 'B)', 'C)', 'D)'),
  font.label = list(size = fsize + 3, family = font))
ggsave('../writing/figures/Effect_size_power_figure.png',
       width = 16, height = 8, units = 'cm', dpi = 600, plot = p)
