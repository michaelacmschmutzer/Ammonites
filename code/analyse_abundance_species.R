# Analysis of abundance datasets at the species level
setwd("~/Documents/Projects/Ammonites/code/")

library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(scales)
library(stringr)
library(systemfonts)

source('statistical-functions.R')
source('geospatial-functions.R')


# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9

###################  Read in data ################### 

cepha_rota <- read.csv('../data/cephalopods_palaeorotated_species.csv')

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

################### Analyse raw abundance data ###################

# Count raw abundance in dataset vs survival
tally_species_nauti <- cepha_rota %>% 
  filter(order == 'Nautilida') %>%
  count(species)
tally_species_ammon <- cepha_rota %>% 
  filter(!order == 'Nautilida') %>%
  count(species)

# Add survivalship info
tally_nauti_surv <- merge(tally_species_nauti, surv_nauti)
tally_ammon_surv <- merge(tally_species_ammon, surv_ammon)


set.seed(36)
pammon <- ggplot(data = tally_ammon_surv, aes(x = survival, y = n)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
#  geom_text_repel(data = subset(tally_ammon_surv, survival == TRUE),
#                  aes(label = species),
#                  nudge_x = - 0.1,
#                  family = font,
#                  size = 2.5,
#                  fontface = 'italic',
#                  box.padding   = unit(0.5, 'lines'),
#                  point.padding = unit(0.6, 'lines')) +
  stat_summary(
    fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
  stat_summary(
    fun.max = function(z) { quantile(z, 0.75) },
    fun.min = function(z) { quantile(z, 0.25) },
    geom = 'errorbar', color = errbarcol, width = 0.1) +
  labs(x = 'Ammonoids', y = 'Number of occurrences per species') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )

pnauti <- ggplot(data = tally_nauti_surv, aes(x = survival, y = n)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, colour = orange) +
  #  geom_text_repel(data = subset(tally_ammon_surv, survival == TRUE),
  #                  aes(label = species),
  #                  nudge_x = - 0.1,
  #                  family = font,
  #                  size = 2.5,
  #                  fontface = 'italic',
  #                  box.padding   = unit(0.5, 'lines'),
  #                  point.padding = unit(0.6, 'lines')) +
  stat_summary(
    fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
  stat_summary(
    fun.max = function(z) { quantile(z, 0.75) },
    fun.min = function(z) { quantile(z, 0.25) },
    geom = 'errorbar', color = errbarcol, width = 0.1) +
  labs(x = 'Nautilids', y = 'Number of occurrences per species') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )
p <- grid.arrange(pammon, pnauti, ncol = 2)

wilcox.test(n ~ survival, data = tally_ammon_surv)
wilcox.test(n ~ survival, data = tally_nauti_surv)

################### Raw abundance Europe & N. America ###################

# Consider only European occurrences
tally_species_ammon_eur <- cepha_rota %>% 
  filter(!order == 'Nautilida') %>%
  filter(lng >= -10 & lng <= 60) %>%
  filter(lat <= 70 & lat >= 30) %>%
  count(species)

# Add survivalship info
tally_ammon_surv_eur <- merge(tally_species_ammon_eur, surv_ammon)

ggplot(data = tally_ammon_surv_eur, aes(x = survival, y = n)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  #  geom_text_repel(data = subset(tally_ammon_surv, survival == TRUE),
  #                  aes(label = species),
  #                  nudge_x = - 0.1,
  #                  family = font,
  #                  size = 2.5,
  #                  fontface = 'italic',
  #                  box.padding   = unit(0.5, 'lines'),
  #                  point.padding = unit(0.6, 'lines')) +
  stat_summary(
    fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
  stat_summary(
    fun.max = function(z) { quantile(z, 0.75) },
    fun.min = function(z) { quantile(z, 0.25) },
    geom = 'errorbar', color = errbarcol, width = 0.1) +
  labs(x = 'Ammonoids', y = 'Number of occurrences per species') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )

wilcox.test(n ~ survival, data = tally_ammon_surv_eur)

# Consider only North American occurrences
tally_species_ammon_nam <- cepha_rota %>% 
  filter(!order == 'Nautilida') %>%
  filter(lng >= -135 & lng <= -60) %>%
  filter(lat <= 70 & lat >= 15) %>%
  count(species)

# Add survivalship info
tally_ammon_surv_nam <- merge(tally_species_ammon_nam, surv_ammon)

ggplot(data = tally_ammon_surv_nam, aes(x = survival, y = n)) +
  geom_jitter(position = position_jitter(0.1), cex = 3, color = darkblue) +
  #  geom_text_repel(data = subset(tally_ammon_surv, survival == TRUE),
  #                  aes(label = species),
  #                  nudge_x = - 0.1,
  #                  family = font,
  #                  size = 2.5,
  #                  fontface = 'italic',
  #                  box.padding   = unit(0.5, 'lines'),
  #                  point.padding = unit(0.6, 'lines')) +
  stat_summary(
    fun = median, geom = 'point', shape = 18, size = 3.5, color = errbarcol) +
  stat_summary(
    fun.max = function(z) { quantile(z, 0.75) },
    fun.min = function(z) { quantile(z, 0.25) },
    geom = 'errorbar', color = errbarcol, width = 0.1) +
  labs(x = 'Ammonoids', y = 'Number of occurrences per species') +
  scale_x_discrete(labels = c('FALSE' = 'Extinct', 'TRUE' = 'Survived')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.y = element_text(size = fsize, family = font),
    axis.title.x = element_text(size = fsize, family = font),
  )

wilcox.test(n ~ survival, data = tally_ammon_surv_nam)
