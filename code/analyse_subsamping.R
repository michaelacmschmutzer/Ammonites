# Analyse sub-sampling results
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(forcats)
library(ggplot2)
library(ggtext)
library(reshape2)
library(wesanderson)

source('statistical-functions.R')

# Command line arguments if run from bash
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0){
  cephas <- 'ammonoids'  # 'nautilids' 'ammonoids'
} else {
  cephas <- args[1]
}

# Get the raw geographic distributions
distrfile <- paste('../results/geographic_distributions/', cephas,
  '_distributions_genus.csv', sep = '')
geo.distr <- read.csv(distrfile)

# Get bootstrapping and jackknifing results
sampldir <- paste('../results/subsampling_distributions/', cephas, '/', sep = '')
bootstrap.df <- read.csv(paste(sampldir, 'bootstrap.csv', sep = ''))
jackknife.df <- read.csv(paste(sampldir, 'jackknife.csv', sep = ''))

# Add information about survival
survfile <- paste('../data/', cephas, '_extinction_genus.csv', sep = '')
survival <- read.csv(survfile)
survival <- select(survival, c('genus', 'survival'))
bootstrap.df.surv <- merge(bootstrap.df, survival, by = 'genus')
jackknife.df.surv <- merge(jackknife.df, survival, by = 'genus')

extinct <- survival[survival$survival==FALSE, 'genus']

# Rotation models
rotamodels <- c('MERDITH2021', 'TorsvikCocks2017',
                'GOLONKA', 'MATTHEWS2016_pmag_ref',
                'PALEOMAP')

# Plotting colours
darjeeling_palette <- wes_palette('Darjeeling2', n = length(rotamodels))

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


for (m in rotamodels) {
  
  ##### BOOTSTRAPPING 
  
  p <- bootstrap.df.surv %>%
    mutate(genus = fct_reorder(genus, survival)) %>%
    filter(model == m) %>%
    ggplot(aes(x = genus, y = area.km2)) +
    geom_jitter(position = position_jitter(0.4), cex = 3) +
    stat_summary(fun = median, geom = 'point', shape = 18, size = 3.5,
                 color = 'red3') +
    stat_summary(fun.max = function(z) { quantile(z, 0.75) },
                 fun.min = function(z) { quantile(z, 0.25) },
                 geom = 'errorbar', color = 'red3', width = 0.1) +
    scale_y_log10('Area (km²)',
                  breaks = scales::trans_breaks('log10', function(x) 10^x),
                  labels = scales::trans_format('log10', scales::math_format(10^.x))) +
    labs(x = 'Genus') +
    scale_x_discrete(labels = ~ if_else(
      .x %in% extinct, paste0("<span style='color: red3'>", .x, "</span>"), .x
    )) +
    theme_classic() +
    theme(
      axis.text.x = element_markdown(size = 8, face = 'italic',
                                     colour = 'black', angle = rot, hjust = hadj, vjust = vadj),
      axis.text.y = element_text(size = 8, colour = 'black'),
      axis.title = element_text(size = 8))
  plotfile <- paste('../results/subsampling_distributions/', cephas,
                    '/', m, '_bootstrap.png', sep = '')
  ggsave(plotfile, width = 12, height = 8, units = 'cm', plot = p)
  
  ##### JACKKNIFING
  
  
  p <- jackknife.df.surv %>%
    mutate(genus = fct_reorder(genus, survival)) %>%
    filter(model == m) %>%
    ggplot(aes(x = genus, y = area.km2)) +
    geom_jitter(position = position_jitter(0.4), cex = 3) +
    stat_summary(fun = median, geom = 'point', shape = 18, size = 3.5,
                 color = 'red3') +
    stat_summary(
      fun.max = function(z) { quantile(z, 0.75) },
      fun.min = function(z) { quantile(z, 0.25) },
      geom = 'errorbar', color = 'red3', width = 0.1) +
    scale_y_log10('Area (km²)',
                  breaks = scales::trans_breaks('log10', function(x) 10^x),
                  labels = scales::trans_format('log10', scales::math_format(10^.x))) +
    labs(x = 'Genus') +
    scale_x_discrete(labels = ~ if_else(
      .x %in% extinct, paste0("<span style='color: red3'>", .x, "</span>"), .x
    )) +
    theme_classic() +
    theme(
      axis.text.x = element_markdown(size = 8, face = 'italic',  
                                     colour = 'black', angle = rot, hjust = hadj, vjust = vadj),
      axis.text.y = element_text(size = 8, colour = 'black'),
      axis.title = element_text(size = 8))
  plotfile <- paste('../results/subsampling_distributions/', cephas,
                    '/', m, '_jackknife.png', sep = '')
  ggsave(plotfile, width = 12, height = 8, units = 'cm', plot = p)
}

# Assess the robustness of the geographic distribution to sub-sampling and 
# indicate the variation caused by different plate tectonic reconstructions

# Correlate the raw area estimates with the median bootstrap and jackknife ones

# Simplify the raw area estimates and reformat
areacol <- paste(rotamodels, '.area.km2', sep = '')
names(areacol) <- rotamodels
geo <- geo.distr %>%
  select(all_of(c('genus', c('genus', areacol)))) 
geo_melt <- melt(geo, id = 'genus')
names(geo_melt) <- c('genus', 'model', 'area')

# Calculate medians for boostrapped samples
median_boot <- bootstrap.df %>%
  group_by(genus, model) %>%
  summarise(boot.area = median(area.km2), .groups = 'keep')

# Calculate medians for jackknifed samples
median_jack <- jackknife.df %>%
  group_by(genus, model) %>%
  summarise(jack.area = median(area.km2), .groups = 'keep')

# Combine the data-frames so I can easily compare the medians
geo_boot <- merge(geo_melt, median_boot, by = c('genus', 'model'))
geo_comb <- merge(geo_boot, median_jack, by = c('genus', 'model'))

# Summary dataframe. Compare the correlations between the raw and sub-sampled
# areas to determine the effect of sampling bias. Compare the variation between
# area (medians) between different plate tectonic models to see the influence
# of plate tectonic model assumptions
geo_summary <- data.frame(genus = unique(geo_comb$genus))

boot.corr <- geo_comb %>%
  group_by(genus) %>%
  summarise(rho.boot = get.spearmans.rho(area, boot.area), .groups = 'keep')
jack.corr <- geo_comb %>%
  group_by(genus) %>%
  summarise(rho.jack = get.spearmans.rho(area, jack.area), .groups = 'keep')

geo_summary <- merge(geo_summary, boot.corr, by = 'genus')
geo_summary <- merge(geo_summary, jack.corr, by = 'genus')

geo_summary$med.dev.area = 
  with(geo_comb, tapply(area, genus, scaled.median.deviation))
geo_summary$med.dev.boot = 
  with(geo_comb, tapply(boot.area, genus, scaled.median.deviation))
geo_summary$med.dev.jack = 
  with(geo_comb, tapply(jack.area, genus, scaled.median.deviation))

resdir <- paste('../results/subsampling_distributions/', cephas, '/', sep = '')
write.csv(geo_summary, paste(resdir, 'geo_sub_summary.csv', sep = ''),
row.names = FALSE)

